// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The turn policy (copilot/turnPolicy.js) at its interface: rows of events →
// steps, with a fake clock and a stub judge — every ending, the last-request
// note once, 'first' against 'final', a reminder and the last-request note in
// one request, the warmer retry, the waits off the turn's time, when a
// send-back is allowed, the model that refuses tools, the one done shape, the
// errors, stoppedByFailure and STOPPED_BY. Then the policy inside the loop
// (agentLoop.js): the three ending paths share one key set, and a reminder
// request that becomes the last carries both notes in order. The dock and the
// harness reading the endings alike: src/__tests__/turnEndings.test.jsx.

import { describe, expect, it, vi } from 'vitest';
import { runCopilotTurn } from '../copilot/agentLoop';
import { ModelTimeoutError } from '../copilot/robustness';
import { applyTurnEvent, newAssistantTurn } from '../copilot/turnState';
import { STOPPED_BY, createTurnPolicy, stoppedByFailure } from '../copilot/turnPolicy';
import { createVerification } from '../copilot/verification';
import { toolBudgetNote } from '../prompts/copilot';

// The done event's keys, in order: one shape on every path.
const DONE_KEYS = ['type', 'text', 'rawText', 'outcome', 'rounds', 'stoppedBy', 'cut', 'reason', 'audit', 'references', 'gates'];

const LIMITS = Object.freeze({ maxRounds: 6, maxTurnMs: 600000, maxNoProgressRounds: 2 });

// A clock the test sets; it counts its reads.
const fakeClock = (time = 0) => {
    const clock = { time, reads: 0, now: () => { clock.reads += 1; return clock.time; } };
    return clock;
};

// A judge that records each call and answers `answer(text, options)` (by
// default the verification's done for that text, references unfilled).
const stubJudge = (answer = null) => vi.fn((text, options) => (answer ? answer(text, options) : {
    done: {
        type: 'done', text, outcome: { status: 'achieved', note: 'n.', body: text, claimed: null, computed: true, reasons: [] },
        rounds: options.rounds, stoppedBy: options.stoppedBy, cut: options.cut, audit: [], references: { written: 0, filled: 0, unknown: [] }, gates: []
    }
}));

const policyOf = ({ limits = {}, clock = fakeClock(), tools = true, judge = stubJudge() } = {}) => ({
    policy: createTurnPolicy({ limits: { ...LIMITS, ...limits }, clock, tools, judge }),
    clock,
    judge
});

const START = Object.freeze({ kind: 'start', waitedMs: 0 });
const reply = (text = '', { calls = [], looped = false, cut = false } = {}) => ({ kind: 'reply', text, calls, looped, cut });
const CALLS = Object.freeze([{ id: 'c0', name: 'echo', arguments: '{}' }]);
const round = (progressed = true, waitedMs = 0) => ({ kind: 'round', progressed, waitedMs });

// The steps for a row of events.
const drive = (policy, events) => events.map((event) => policy.next(event));
const asks = (steps) => steps.filter((step) => step.next === 'ask');

const EMPTY_FINAL = 'The analysis budget for this question was used before the model wrote an answer. '
    + 'The analyses it ran are listed above.';
const EMPTY_STALLED = 'The model kept repeating analyses that added nothing new, and wrote no answer. '
    + 'The analyses it ran are listed above.';

describe('the turn policy: events → steps', () => {
    it('an answer: the first request offers the tools, the reply is judged, and the turn ends on its done', () => {
        const { policy, judge } = policyOf();
        const [first, end] = drive(policy, [START, reply('The value is 4.')]);
        expect(first).toEqual({ next: 'ask', round: 0, why: 'first', toolChoice: 'auto', retry: null, notes: [], sendBack: null });
        expect(judge).toHaveBeenCalledWith('The value is 4.', { canSendBack: true, rounds: 0, stoppedBy: 'answer', cut: false });
        expect(end).toEqual({
            next: 'end',
            say: null,
            done: {
                type: 'done', text: 'The value is 4.', rawText: 'The value is 4.',
                outcome: { status: 'achieved', note: 'n.', body: 'The value is 4.', claimed: null, computed: true, reasons: [] },
                rounds: 0, stoppedBy: 'answer', cut: false, reason: null, audit: [], references: { written: 0, filled: 0, unknown: [] }, gates: []
            }
        });
    });

    it('a reply with calls runs them; their round makes the next request, with the results', () => {
        const { policy } = policyOf();
        const steps = drive(policy, [START, reply('Checking.', { calls: CALLS }), round(), reply('', { calls: CALLS }), round(), reply('Done.')]);
        expect(steps.map((step) => step.next)).toEqual(['ask', 'run', 'ask', 'run', 'ask', 'end']);
        expect(steps[1]).toEqual({ next: 'run' });
        expect(asks(steps).map(({ round: index, why, toolChoice }) => [index, why, toolChoice])).toEqual([[0, 'first', 'auto'], [1, 'tools', 'auto'], [2, 'tools', 'auto']]);
        expect(steps.at(-1).done).toMatchObject({ rounds: 2, stoppedBy: 'answer' });
    });

    it('the round budget: the next request asks for the answer — why "final", tool_choice none, the note once', () => {
        const { policy, judge } = policyOf({ limits: { maxRounds: 1 } });
        const steps = drive(policy, [START, reply('', { calls: CALLS }), round(), reply('Best effort.')]);
        expect(steps[2]).toEqual({ next: 'ask', round: 1, why: 'final', toolChoice: 'none', retry: null, notes: [toolBudgetNote(1)], sendBack: null });
        // No send-back on the last request: the answer is judged as it stands.
        expect(judge).toHaveBeenLastCalledWith('Best effort.', { canSendBack: false, rounds: 1, stoppedBy: 'max_rounds', cut: false });
        expect(steps.at(-1).done).toMatchObject({ stoppedBy: 'max_rounds', rounds: 1 });
    });

    it('the first request stays "first" when it is already the last (maxRounds 0), with the note', () => {
        const { policy } = policyOf({ limits: { maxRounds: 0 } });
        expect(policy.next(START)).toEqual({ next: 'ask', round: 0, why: 'first', toolChoice: 'none', retry: null, notes: [toolBudgetNote(0)], sendBack: null });
        expect(policy.next(reply('Answer.')).done.stoppedBy).toBe('max_rounds');
    });

    it('rounds that add nothing: two in a row end the analyses (no_progress); a round that adds resets the count', () => {
        const stalled = policyOf().policy;
        const steps = drive(stalled, [START, reply('', { calls: CALLS }), round(false), reply('', { calls: CALLS }), round(false), reply('Answer.')]);
        expect(steps[4]).toMatchObject({ why: 'final', toolChoice: 'none' });
        expect(steps[4].notes).toEqual([{
            role: 'user',
            content: 'Your last 2 steps repeated earlier analyses or failed, so they added nothing new. Do not call them again. '
                + 'Answer now from the results above, say what is still missing, or ask the user with ask_user if a choice blocks you. '
                + 'End with the Outcome line.'
        }]);
        expect(steps.at(-1).done).toMatchObject({ stoppedBy: 'no_progress', rounds: 2 });
        const reset = policyOf().policy;
        const later = drive(reset, [START, reply('', { calls: CALLS }), round(false), reply('', { calls: CALLS }), round(true), reply('', { calls: CALLS }), round(false)]);
        expect(asks(later).map((step) => step.toolChoice)).toEqual(['auto', 'auto', 'auto', 'auto']);
    });

    it('an empty reply gets one reminder, then the turn ends empty with why; the judge never reads blank text', () => {
        const cases = [
            [{}, 'had no answer and no analysis call. ', 'the model wrote no answer, even after a reminder'],
            [{ cut: true }, 'had no answer and no analysis call (it reached the length limit while reasoning). ', 'the model spent its reply on reasoning and reached the length limit'],
            [{ looped: true }, 'kept repeating the same passage, so it was stopped. ', 'the model kept repeating the same passage, even after a reminder']
        ];
        for (const [flags, said, reason] of cases) {
            const { policy, judge } = policyOf();
            const [, reminder, end] = drive(policy, [START, reply('  ', flags), reply('', flags)]);
            expect(reminder).toMatchObject({ next: 'ask', round: 0, why: 'reminder', toolChoice: 'auto', retry: flags.looped ? 'warmer' : null, sendBack: null });
            expect(reminder.notes).toEqual([{
                role: 'user',
                content: `Your last reply ${said}Answer the question now from the context and the results above, briefly, or call one analysis if `
                    + 'one is still missing. End with the Outcome line.'
            }]);
            expect(end).toEqual({
                next: 'end', say: null,
                done: { type: 'done', text: '', rawText: '', outcome: null, rounds: 0, stoppedBy: 'empty', cut: Boolean(flags.cut), reason, audit: [], references: null, gates: [] }
            });
            expect(judge).not.toHaveBeenCalled();
        }
    });

    it('the last request came back empty, with nothing standing: a fixed answer that says so (max_rounds, no_progress)', () => {
        const budget = policyOf({ limits: { maxRounds: 1 } });
        const steps = drive(budget.policy, [START, reply('', { calls: CALLS }), round(), reply(''), reply('')]);
        // The reminder on the last request carries no second budget note.
        expect(steps[3]).toEqual({
            next: 'ask', round: 1, why: 'reminder', toolChoice: 'none', retry: null, notes: [expect.objectContaining({ role: 'user' })], sendBack: null
        });
        expect(steps[4]).toEqual({
            next: 'end',
            say: EMPTY_FINAL,
            done: {
                type: 'done', text: EMPTY_FINAL, rawText: EMPTY_FINAL,
                outcome: { status: 'not achieved', note: 'no answer was written within the analysis budget', body: EMPTY_FINAL },
                rounds: 1, stoppedBy: 'max_rounds', cut: false, reason: null, audit: [], references: null, gates: []
            }
        });
        expect(budget.judge).not.toHaveBeenCalled();
        const stalled = policyOf();
        const end = drive(stalled.policy, [START, reply('', { calls: CALLS }), round(false), reply('', { calls: CALLS }), round(false), reply(''), reply('')]).at(-1);
        expect(end.say).toBe(EMPTY_STALLED);
        expect(end.done).toMatchObject({
            text: EMPTY_STALLED, rawText: EMPTY_STALLED, stoppedBy: 'no_progress', reason: null, references: null, gates: [],
            outcome: { status: 'not achieved', note: 'the analyses stopped adding anything new and no answer was written', body: EMPTY_STALLED }
        });
    });

    it('a reminder request that becomes the last carries the reminder, then the last-request note, and why "final"', () => {
        const { policy, clock } = policyOf({ limits: { maxTurnMs: 100 } });
        policy.next(START);
        policy.next(reply('', { calls: CALLS }));
        clock.time = 50;
        expect(policy.next(round())).toMatchObject({ why: 'tools', toolChoice: 'auto', notes: [] });
        clock.time = 150;
        const reminder = policy.next(reply(''));
        expect(reminder).toMatchObject({ why: 'final', toolChoice: 'none' });
        expect(reminder.notes.map((note) => note.content.slice(0, 40))).toEqual([
            'Your last reply had no answer and no ana',
            'You have used the analysis budget for th'
        ]);
        expect(reminder.notes[1]).toEqual(toolBudgetNote(1));
    });

    it('one reminder a turn: a later empty reply ends the turn', () => {
        const { policy } = policyOf();
        const steps = drive(policy, [START, reply(''), reply('', { calls: CALLS }), round(), reply('')]);
        expect(asks(steps).map((step) => step.why)).toEqual(['first', 'reminder', 'tools']);
        expect(steps.at(-1).done).toMatchObject({ stoppedBy: 'empty', reason: 'the model wrote no answer, even after a reminder', rounds: 1 });
    });

    it('after a looping reply every later request is asked warmer, the retry without tools included', () => {
        const { policy } = policyOf();
        const steps = drive(policy, [START, reply('', { looped: true }), reply('', { calls: CALLS }), round(), reply('Answer.')]);
        expect(asks(steps).map((step) => step.retry)).toEqual([null, 'warmer', 'warmer']);
        const again = policyOf().policy;
        expect(drive(again, [START, reply('', { looped: true }), { kind: 'no_tools' }]).map((step) => step.retry)).toEqual([null, 'warmer', 'warmer']);
        // A loop after the reminder ends the turn: there is no request left to warm.
        const late = policyOf().policy;
        expect(drive(late, [START, reply(''), reply('', { looped: true })]).at(-1).done.reason)
            .toBe('the model kept repeating the same passage, even after a reminder');
    });

    it('the user\'s waits are off the turn\'s time, which is read at most once an event and never in round 0', () => {
        const clock = fakeClock(0);
        const { policy } = policyOf({ clock, limits: { maxTurnMs: 100 } });
        expect(clock.reads).toBe(1);
        // A Workbench question took 1,000 ms before the first request.
        clock.time = 1000;
        policy.next({ kind: 'start', waitedMs: 1000 });
        policy.next(reply('', { calls: CALLS }));
        expect(clock.reads).toBe(1);
        // An ask_user call in the round waited 500 ms: 1,550 − 1,000 − 500 = 50 ms of turn time.
        clock.time = 1550;
        expect(policy.next(round(true, 500))).toMatchObject({ toolChoice: 'auto' });
        expect(clock.reads).toBe(2);
        policy.next(reply('', { calls: CALLS }));
        clock.time = 1700;
        expect(policy.next(round(true, 0))).toMatchObject({ why: 'final', toolChoice: 'none' });
        expect(clock.reads).toBe(3);
        // The last request: no send-back to decide, so no reading.
        policy.next(reply('Late.'));
        expect(clock.reads).toBe(3);
    });

    it('a send-back is allowed only with tools on, before the last request and with turn time left', () => {
        const options = (events, setup = {}) => {
            const { policy, judge, clock } = policyOf(setup);
            events.forEach((event) => { if (typeof event === 'function') event(clock); else policy.next(event); });
            return judge.mock.calls.at(-1)[1];
        };
        expect(options([START, reply('A.')])).toMatchObject({ canSendBack: true, stoppedBy: 'answer' });
        expect(options([START, reply('A.')], { tools: false })).toMatchObject({ canSendBack: false, stoppedBy: 'answer' });
        expect(options([START, reply('', { calls: CALLS }), round(), reply('A.')], { limits: { maxRounds: 1 } }))
            .toMatchObject({ canSendBack: false, stoppedBy: 'max_rounds' });
        // Asked with time left (50 ms), answered after it was spent (150 ms): judged as an answer, never sent back.
        expect(options([
            START, reply('', { calls: CALLS }), (clock) => { clock.time = 50; }, round(), (clock) => { clock.time = 150; }, reply('A.')
        ], { limits: { maxTurnMs: 100 } })).toMatchObject({ canSendBack: false, stoppedBy: 'answer' });
    });

    it('at the verification\'s own cap the answer stands: the policy allows, the verification declines', () => {
        const verification = createVerification({ maxSendBacks: 1 });
        const judge = vi.fn(verification.judge);
        const { policy } = policyOf({ judge });
        const text = 'Ti forms a superstructure. Site 2 has a preferred displacement direction.';
        const steps = drive(policy, [START, reply(text), reply(text)]);
        expect(steps[1]).toMatchObject({ next: 'ask', why: 'gate', notes: [], sendBack: { id: 'order-between-cells' } });
        expect(judge.mock.calls.map(([, { canSendBack }]) => canSendBack)).toEqual([true, true]);
        expect(steps[2].done).toMatchObject({ stoppedBy: 'answer', gates: [{ round: 0, id: 'order-between-cells', missing: ['analyze_modulations'] }] });
        expect(steps[2].done.audit.map((finding) => finding.id)).toEqual(['order-between-cells', 'direction']);
    });

    it('the answer the gate withdrew stands when the model then writes nothing, judged with canSendBack false', () => {
        const sendBack = { id: 'clusters', what: 'clusters', missing: ['analyze_clusters'], notice: 'n', nudge: 'Run analyze_clusters.' };
        const judge = stubJudge((text, options) => (options.canSendBack ? { sendBack } : { done: { text: `[${text}]`, outcome: null } }));
        const { policy } = policyOf({ judge });
        const steps = drive(policy, [START, reply('The clusters breathe.'), reply(''), reply('')]);
        expect(steps[1]).toEqual({ next: 'ask', round: 0, why: 'gate', toolChoice: 'auto', retry: null, notes: [], sendBack });
        // The reminder first, then the withdrawn answer stands.
        expect(steps[2]).toMatchObject({ why: 'reminder' });
        expect(judge).toHaveBeenLastCalledWith('The clusters breathe.', { canSendBack: false, rounds: 0, stoppedBy: 'answer', cut: false });
        expect(steps[3].done).toEqual({
            type: 'done', text: '[The clusters breathe.]', rawText: '[The clusters breathe.]', outcome: null,
            rounds: 0, stoppedBy: 'answer', cut: false, reason: null, audit: [], references: null, gates: []
        });
        // With the reminder spent before the send-back, it stands at once.
        const early = policyOf({ judge: stubJudge((text, options) => (options.canSendBack ? { sendBack } : { done: { text } })) });
        expect(drive(early.policy, [START, reply(''), reply('The clusters breathe.'), reply('')]).at(-1).done.text).toBe('The clusters breathe.');
    });

    it('a model that refuses tools: no_tools makes the retry without them; every later request has none', () => {
        const { policy, judge } = policyOf();
        const steps = drive(policy, [START, { kind: 'no_tools' }, reply(''), reply('Plain answer.')]);
        expect(steps[1]).toEqual({ next: 'ask', round: 0, why: 'no_tools', toolChoice: undefined, retry: null, notes: [], sendBack: null });
        expect(steps[2]).toMatchObject({ why: 'reminder', toolChoice: undefined });
        expect(judge).toHaveBeenCalledWith('Plain answer.', { canSendBack: false, rounds: 0, stoppedBy: 'answer', cut: false });
        // A request with no tools has no last request: never 'final', never a budget note.
        const off = policyOf({ tools: false, limits: { maxRounds: 0 } }).policy;
        expect(off.next(START)).toEqual({ next: 'ask', round: 0, why: 'first', toolChoice: undefined, retry: null, notes: [], sendBack: null });
    });

    it('judge returning null gives the empty ending; its own errors propagate', () => {
        const before = policyOf({ judge: vi.fn(() => null) }).policy;
        expect(drive(before, [START, reply('Text.')]).at(-1).done).toMatchObject({ stoppedBy: 'empty', reason: 'the model wrote no answer, even after a reminder', text: '' });
        const last = policyOf({ judge: vi.fn(() => undefined), limits: { maxRounds: 0 } }).policy;
        expect(drive(last, [START, reply('Text.')]).at(-1)).toMatchObject({ say: EMPTY_FINAL, done: { stoppedBy: 'max_rounds' } });
        const failure = new Error('verification.judge: this turn\'s answer is judged; one verification serves one turn');
        const throwing = policyOf({ judge: () => { throw failure; } }).policy;
        throwing.next(START);
        expect(() => throwing.next(reply('Text.'))).toThrow(failure);
    });

    it('the one done shape: the 11 keys, in order, on every end', () => {
        const ends = [];
        const collect = (setup, events) => ends.push(drive(policyOf(setup).policy, events).at(-1));
        collect({}, [START, reply('An answer.')]);
        collect({ judge: stubJudge((text) => ({ done: { text: 'filled', rawText: text, outcome: null, audit: [], references: { written: 1, filled: 1, unknown: [] }, gates: [] } })) }, [START, reply('{{r1.x}}')]);
        collect({ limits: { maxRounds: 1 } }, [START, reply('', { calls: CALLS }), round(), reply('Late answer.')]);
        collect({}, [START, reply('', { calls: CALLS }), round(false), reply('', { calls: CALLS }), round(false), reply('Answer.')]);
        collect({ limits: { maxRounds: 1 } }, [START, reply('', { calls: CALLS }), round(), reply(''), reply('')]);
        collect({}, [START, reply('', { calls: CALLS }), round(false), reply('', { calls: CALLS }), round(false), reply(''), reply('')]);
        collect({}, [START, reply(''), reply('')]);
        collect({ tools: false }, [START, reply(''), reply('')]);
        expect(ends.map((end) => end.next)).toEqual(Array(8).fill('end'));
        ends.forEach((end) => expect(Object.keys(end.done)).toEqual(DONE_KEYS));
        expect(ends.map((end) => end.done.stoppedBy)).toEqual(['answer', 'answer', 'max_rounds', 'no_progress', 'max_rounds', 'no_progress', 'empty', 'empty']);
        expect(ends[1].done).toMatchObject({ text: 'filled', rawText: '{{r1.x}}' });
        // Every ending the policy gives, and every one a thrown turn gets, is a row of STOPPED_BY.
        const endings = new Set([...ends.map((end) => end.done.stoppedBy), stoppedByFailure(new DOMException('Aborted', 'AbortError')), stoppedByFailure(new Error('x'))]);
        expect([...endings].sort()).toEqual(Object.keys(STOPPED_BY).sort());
    });

    it('every ask step: tool_choice undefined exactly when tools are off, at most two notes', () => {
        const rows = [
            [{}, [START, reply(''), reply('', { calls: CALLS }), round(), reply('', { looped: true })]],
            [{ limits: { maxRounds: 1 } }, [START, reply('', { calls: CALLS }), round(), reply('')]],
            [{ tools: false }, [START, reply('')]],
            [{}, [START, { kind: 'no_tools' }, reply('')]]
        ];
        for (const [setup, events] of rows) {
            const steps = asks(drive(policyOf(setup).policy, events));
            const toolsOff = (step) => step.why === 'no_tools' || setup.tools === false || steps.slice(0, steps.indexOf(step)).some((earlier) => earlier.why === 'no_tools');
            steps.forEach((step) => {
                expect(step.toolChoice === undefined).toBe(toolsOff(step));
                expect(step.notes.length).toBeLessThanOrEqual(2);
            });
        }
    });
});

describe('the turn policy: errors', () => {
    it('bad configuration is a TypeError; Infinity limits are allowed', () => {
        const ok = { limits: LIMITS, clock: fakeClock(), tools: true, judge: stubJudge() };
        expect(() => createTurnPolicy(ok)).not.toThrow();
        expect(() => createTurnPolicy({ ...ok, limits: { maxRounds: Infinity, maxTurnMs: Infinity, maxNoProgressRounds: Infinity } })).not.toThrow();
        for (const bad of [
            undefined,
            { ...ok, limits: null },
            { ...ok, limits: { ...LIMITS, maxRounds: undefined } },
            { ...ok, limits: { ...LIMITS, maxTurnMs: Number.NaN } },
            { ...ok, limits: { ...LIMITS, maxNoProgressRounds: '2' } },
            { ...ok, clock: {} },
            { ...ok, clock: null },
            { ...ok, tools: 'yes' },
            { ...ok, judge: null }
        ]) {
            expect(() => createTurnPolicy(bad)).toThrow(TypeError);
        }
    });

    it('an event out of order is a TypeError, and a turn ends exactly once', () => {
        const fresh = () => policyOf().policy;
        expect(() => fresh().next(reply('x'))).toThrow(/a 'reply' event is out of order: expected start/);
        expect(() => fresh().next(round())).toThrow(TypeError);
        const started = fresh();
        started.next(START);
        expect(() => started.next(START)).toThrow(/expected reply or no_tools/);
        expect(() => started.next(round())).toThrow(TypeError);
        const running = fresh();
        drive(running, [START, reply('', { calls: CALLS })]);
        expect(() => running.next(reply('x'))).toThrow(/expected round/);
        const ended = fresh();
        drive(ended, [START, reply('Answer.')]);
        for (const event of [START, reply('Again.'), round(), { kind: 'no_tools' }]) {
            expect(() => ended.next(event)).toThrow(/the turn has ended/);
        }
        expect(() => fresh().next({ kind: 'later' })).toThrow(TypeError);
        expect(() => fresh().next(null)).toThrow(TypeError);
    });

    it('a malformed event, calls a request did not allow, or no_tools after tools ran is a TypeError', () => {
        const at = (setup, events) => {
            const { policy } = policyOf(setup);
            events.slice(0, -1).forEach((event) => policy.next(event));
            return () => policy.next(events.at(-1));
        };
        expect(at({}, [{ kind: 'start' }])).toThrow(/waitedMs/);
        expect(at({}, [{ kind: 'start', waitedMs: Number.NaN }])).toThrow(/waitedMs/);
        expect(at({}, [START, { kind: 'reply', text: 'x' }])).toThrow(/a reply must be/);
        expect(at({}, [START, { ...reply('x'), text: 4 }])).toThrow(/a reply must be/);
        expect(at({}, [START, { ...reply('x'), looped: 'no' }])).toThrow(/a reply must be/);
        expect(at({}, [START, reply('', { calls: CALLS }), { kind: 'round', waitedMs: 0 }])).toThrow(/a round must be/);
        expect(at({}, [START, reply('', { calls: CALLS }), { kind: 'round', progressed: true }])).toThrow(/waitedMs/);
        expect(at({ tools: false }, [START, reply('', { calls: CALLS })])).toThrow(/allowed none/);
        expect(at({ limits: { maxRounds: 0 } }, [START, reply('', { calls: CALLS })])).toThrow(/allowed none/);
        expect(at({ tools: false }, [START, { kind: 'no_tools' }])).toThrow(/no_tools/);
        expect(at({}, [START, reply('', { calls: CALLS }), round(), { kind: 'no_tools' }])).toThrow(/no_tools/);
    });

    it('a judge that sends back what the policy did not allow, or answers otherwise, is a TypeError', () => {
        const sendBack = { id: 'x', what: 'x', missing: [], notice: 'n', nudge: 'm' };
        const last = policyOf({ judge: vi.fn(() => ({ sendBack })), limits: { maxRounds: 0 } }).policy;
        last.next(START);
        expect(() => last.next(reply('A.'))).toThrow(/no send-back was allowed/);
        for (const answer of ['done', {}, { done: 'text' }]) {
            const odd = policyOf({ judge: vi.fn(() => answer) }).policy;
            odd.next(START);
            expect(() => odd.next(reply('A.'))).toThrow(/judge must return/);
        }
    });
});

describe('the endings: STOPPED_BY and stoppedByFailure', () => {
    it('STOPPED_BY holds every ending with its note — the dock\'s notes, word for word — and cutShort only for a budget', () => {
        expect(STOPPED_BY).toEqual({
            answer: { cutShort: false, note: null },
            max_rounds: { cutShort: true, note: 'Answered after using the analysis budget for one question.' },
            no_progress: { cutShort: true, note: 'Answered when repeated analyses stopped adding anything new.' },
            empty: { cutShort: false, note: 'The model returned an empty reply.' },
            stopped: { cutShort: false, note: 'Stopped.' },
            error: { cutShort: false, note: 'Stopped by an error.' }
        });
        expect(Object.isFrozen(STOPPED_BY)).toBe(true);
        expect(Object.values(STOPPED_BY).every((row) => Object.isFrozen(row))).toBe(true);
    });

    it('stoppedByFailure: Stop (an AbortError) is stopped; a stall or any other error, error', () => {
        expect(stoppedByFailure(new DOMException('Aborted', 'AbortError'))).toBe('stopped');
        expect(stoppedByFailure(Object.assign(new Error('x'), { name: 'AbortError' }))).toBe('stopped');
        expect(stoppedByFailure(new ModelTimeoutError('stall', 180000))).toBe('error');
        expect(stoppedByFailure(new Error('HTTP 500'))).toBe('error');
        expect(stoppedByFailure(null)).toBe('error');
        expect(stoppedByFailure(undefined)).toBe('error');
    });
});

// The policy inside the loop (agentLoop.js), against a scripted model port.
describe('the policy in the loop', () => {
    const scripted = (replies, { onRequest = null } = {}) => {
        const requests = [];
        const model = (intent) => {
            requests.push({ ...intent, messages: intent.messages.map((message) => ({ ...message })) });
            onRequest?.(requests.length);
            const reply = replies.shift();
            return (async function* stream() {
                if (reply.text) yield { content: reply.text };
                for (const [index, call] of (reply.calls || []).entries()) {
                    yield { toolCall: { index, id: `c${requests.length}_${index}`, name: call.name, arguments: JSON.stringify(call.args ?? {}) } };
                }
                yield { finishReason: reply.calls ? 'tool_calls' : 'stop' };
            }());
        };
        return { model, requests };
    };
    const echo = (run = async ({ value }) => ({ data: { value }, checks: [] })) => ({
        name: 'echo', description: 'Echo.', parameters: { type: 'object', properties: { value: { type: 'number' } }, required: ['value'] }, run
    });
    const MESSAGES = [{ role: 'system', content: 'S' }, { role: 'user', content: 'question' }];
    const call = { calls: [{ name: 'echo', args: { value: 1 } }] };
    const doneOf = async (replies, options = {}) => {
        const { model, requests } = scripted([...replies]);
        const events = [];
        for await (const event of runCopilotTurn({ model, messages: MESSAGES, tools: [echo()], ...options })) events.push(event);
        return { done: events.at(-1), events, requests };
    };

    it('the three ending paths share one key set: an answer, the last request empty, an empty reply', async () => {
        const answered = await doneOf([call, { text: 'The value is {{r1.value}}.' }]);
        const emptyFinal = await doneOf([call, { text: '' }, { text: '' }], { limits: { maxRounds: 1 } });
        const empty = await doneOf([{ text: '' }, { text: '' }]);
        expect([answered, emptyFinal, empty].map(({ done }) => done.stoppedBy)).toEqual(['answer', 'max_rounds', 'empty']);
        [answered, emptyFinal, empty].forEach(({ done }) => expect(Object.keys(done)).toEqual(DONE_KEYS));
        expect(answered.done).toMatchObject({ text: 'The value is 1.', rawText: 'The value is {{r1.value}}.', reason: null, references: { written: 1, filled: 1, unknown: [] }, gates: [] });
        expect(emptyFinal.done).toMatchObject({ rawText: emptyFinal.done.text, reason: null, references: null, gates: [] });
        expect(empty.done).toMatchObject({ text: '', rawText: '', references: null, gates: [] });
        // The fixed answer is said before the done, as the turn's text.
        expect(emptyFinal.events.at(-2)).toEqual({ type: 'text', delta: emptyFinal.done.text });
    });

    it('the dock reads each path as it read the older shapes (turnState): the added keys change nothing on screen', async () => {
        const fold = (events) => events.reduce((turn, event) => applyTurnEvent(turn, event, 0), newAssistantTurn(0));
        const older = (done) => {
            // Before the turn policy: no reason on an answer, no rawText unless it differed;
            // none of rawText, references or gates on the two empty paths.
            const { rawText, reason, references, gates, ...rest } = done;
            if (done.stoppedBy === 'empty') return { ...rest, reason };
            if (done.outcome && !done.outcome.computed) return rest;
            return { ...rest, ...(rawText !== done.text ? { rawText } : {}), references, gates };
        };
        for (const replies of [[call, { text: 'The value is {{r1.value}}.' }], [{ text: 'Plain.' }]]) {
            const { events } = await doneOf(replies);
            expect(fold(events)).toEqual(fold([...events.slice(0, -1), older(events.at(-1))]));
        }
        const { events: emptyFinal } = await doneOf([call, { text: '' }, { text: '' }], { limits: { maxRounds: 1 } });
        expect(fold(emptyFinal)).toEqual(fold([...emptyFinal.slice(0, -1), older(emptyFinal.at(-1))]));
        const { events: empty } = await doneOf([{ text: '' }, { text: '' }]);
        expect(fold(empty)).toEqual(fold([...empty.slice(0, -1), older(empty.at(-1))]));
    });

    it('a reminder request that becomes the last carries both notes, in order, in one message', async () => {
        const clock = fakeClock(0);
        const tool = echo(async ({ value }) => { clock.time = 50; return { data: { value }, checks: [] }; });
        const { model, requests } = scripted([call, { text: '' }, { text: 'Late.' }], {
            // The empty reply to the second request arrives after the turn's time.
            onRequest: (count) => { if (count === 2) clock.time = 150; }
        });
        const events = [];
        for await (const event of runCopilotTurn({ model, messages: MESSAGES, tools: [tool], clock, limits: { maxTurnMs: 100 } })) events.push(event);
        const asked = events.filter((event) => event.type === 'request').map((event) => [event.why, event.sent.toolChoice]);
        expect(asked).toEqual([['first', 'auto'], ['tools', 'auto'], ['final', 'none']]);
        const last = requests[2].messages.at(-1);
        expect(last.role).toBe('user');
        expect(last.content).toMatch(/^Your last reply had no answer and no analysis call\. [^\n]*\n\nYou have used the analysis budget for this question \(1 rounds of tool calls\)\./);
        expect(events.at(-1)).toMatchObject({ type: 'done', stoppedBy: 'max_rounds', text: 'Late.' });
    });
});
