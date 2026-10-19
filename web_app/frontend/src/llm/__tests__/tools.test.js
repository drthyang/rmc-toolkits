// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The tool contract (copilot/tools.js) at its interface: how every call ends
// (the result statuses and who decides each), the invariant table over one
// turn that reaches every status, the constructors' refusals, the tool
// message the model reads and its compaction, and the helpers. The exact
// texts the model reads are pinned here, at the constructors that write
// them; src/__tests__/modelVisible/toolMessages.test.js pins them as the
// loop sends them, byte for byte.

import { describe, expect, it, vi } from 'vitest';
import { runCopilotTurn } from '../copilot/agentLoop';
import { normalizeAskRequest } from '../copilot/askUser';
import {
    budgetJson, CHOICE_SOURCES, compactToolMessage, ERROR_STATUSES, executeToolCall, fitToBudget, questionResult, ran, repeatedCall,
    REPEAT_NOTE, RESULT_STATUSES, severityOf, SKIP_REASONS, skippedCall, toolMessage, toolResult, unknownToolCall, usableTools,
    withNote
} from '../copilot/tools';
import { applyTurnEvent, closeTurn, newAssistantTurn } from '../copilot/turnState';
import { createVerification } from '../copilot/verification';

const echo = (overrides = {}) => ({
    name: 'echo',
    description: 'Echo the value.',
    parameters: { type: 'object', properties: { value: { type: 'number', maximum: 10 } }, required: ['value'] },
    run: vi.fn(async ({ value }) => ({ data: { value, list: [1, 2, 3] }, checks: [{ level: 'ok', text: 'fine' }] })),
    ...overrides
});
const refuser = echo({
    name: 'refuser',
    description: 'Refuses.',
    run: async ({ value }) => ({ data: { value, refused: 'no' }, checks: [{ level: 'fail', text: 'Refused: not this way.' }], facts: { seen: value } })
});
const call = (tool, args, options = {}) => executeToolCall(tool, typeof args === 'string' ? args : JSON.stringify(args), { timeoutMs: 1000, ...options });
const REFUSAL = [{ level: 'fail', text: 'refused' }];

describe('how a call ends: the result statuses', () => {
    it('is one of five, closed and frozen', () => {
        expect(RESULT_STATUSES).toEqual(['ok', 'refused', 'error', 'unknown_tool', 'skipped']);
        expect(ERROR_STATUSES).toEqual(['error', 'unknown_tool', 'skipped']);
        expect(CHOICE_SOURCES).toEqual(['option', 'other', 'skipped']);
        expect(Object.isFrozen(RESULT_STATUSES) && Object.isFrozen(ERROR_STATUSES) && Object.isFrozen(SKIP_REASONS)).toBe(true);
    });

    it('ok: the tool answered — its data, its checks and facts, no error, no note', async () => {
        const result = await call(echo({ run: async () => ({ data: { v: 1 }, checks: [{ level: 'warn', text: 'a caveat' }], facts: { x: 1 } }) }), { value: 4 });
        expect(result).toEqual({
            resultStatus: 'ok', args: { value: 4 }, data: { v: 1 }, checks: [{ level: 'warn', text: 'a caveat' }], error: null,
            facts: { x: 1 }, ms: expect.any(Number), cached: false, note: null
        });
        expect(Object.isFrozen(result)).toBe(true);
        expect(severityOf(result)).toBe('warn');
        expect(ran(result)).toBe(true);
    });

    it('refused: the tool answered with a fail check (no error) — it ran, and keeps its data and facts', async () => {
        const result = await call(refuser, { value: 2 });
        expect(result).toMatchObject({ resultStatus: 'refused', error: null, data: { value: 2, refused: 'no' }, facts: { seen: 2 } });
        expect(severityOf(result)).toBe('failed');
        expect(ran(result)).toBe(true);
    });

    it('error: arguments that do not parse or that the schema refuses, a throw, the time limit — the message the model reads', async () => {
        const unparsed = await call(echo(), '{value:1}');
        expect(unparsed).toMatchObject({ resultStatus: 'error', args: '{value:1}', data: null, checks: [] });
        expect(unparsed.error).toMatch(/^arguments are not valid JSON \(/);
        expect(await call(echo(), { value: 99, colour: 'red' })).toMatchObject({
            resultStatus: 'error', args: { value: 99, colour: 'red' }, error: 'value must be ≤ 10, got 99',
            checks: [{ level: 'note', text: 'ignored unknown argument "colour"' }]
        });
        expect(await call(echo({ run: async () => { throw new Error('engine said no'); } }), { value: 1, extra: 1 })).toMatchObject({
            resultStatus: 'error', args: { value: 1 }, error: 'engine said no', checks: [{ level: 'note', text: 'ignored unknown argument "extra"' }]
        });
        expect(await call(echo({ run: () => new Promise(() => {}) }), { value: 1 }, { timeoutMs: 20 }))
            .toMatchObject({ resultStatus: 'error', error: 'the analysis did not finish within 0 s', data: null });
        // A throw with no text still names an error (the invariant holds).
        expect((await call(echo({ run: async () => { throw ''; } }), { value: 1 })).error).toBe('the tool failed without saying why');
        const failed = await call(echo({ run: async () => { throw new Error('x'); } }), { value: 1 });
        expect(Object.keys(failed)).not.toContain('facts');
        expect(ran(failed)).toBe(true);
    });

    it('an aborted turn throws, never a result', async () => {
        const controller = new AbortController();
        const pending = call(echo({ run: () => new Promise(() => {}) }), { value: 1 }, { signal: controller.signal, timeoutMs: null });
        controller.abort();
        await expect(pending).rejects.toMatchObject({ name: 'AbortError' });
    });

    it('unknown_tool: the names the turn offers, or none', () => {
        expect(unknownToolCall({}, 'nope', ['echo', 'boom'])).toMatchObject({
            resultStatus: 'unknown_tool', args: {}, data: null, error: 'unknown tool "nope"; available: echo, boom', ms: 0, cached: false
        });
        expect(unknownToolCall('{x', 'nope', []).error).toBe('unknown tool "nope"; available: none');
        expect(ran(unknownToolCall({}, 'nope', []))).toBe(false);
    });

    it('skipped: the loop did not run it — a closed set of reasons, each with its exact text', () => {
        expect(Object.keys(SKIP_REASONS)).toEqual(['call_limit', 'one_question', 'question_limit']);
        expect(skippedCall({ value: 5 }, 'call_limit', 4)).toMatchObject({
            resultStatus: 'skipped', args: { value: 5 }, data: null,
            error: 'skipped: at most 4 tool calls per step — call it again in the next step if it is still needed'
        });
        expect(skippedCall({}, 'one_question').error).toBe('one question at a time: ask the next one after this answer, if it is still needed');
        expect(skippedCall({}, 'question_limit', 2).error)
            .toBe('already asked 2 questions for this request: proceed with your default and say in the answer what you assumed');
        expect(skippedCall({}, 'question_limit', 1).error).toMatch(/^already asked 1 question for this request/);
        expect(ran(skippedCall({}, 'one_question'))).toBe(false);
        expect(() => skippedCall({}, 'too_many')).toThrow(TypeError);
        expect(() => skippedCall({}, 'toString')).toThrow(TypeError);
        expect(() => skippedCall({}, 'call_limit')).toThrow(/names its number/);
    });

    it('a repeat: an ok or refused result again, at no time, cached, with the repeat\'s note', async () => {
        const first = await call(echo(), { value: 2 });
        const again = repeatedCall(first);
        expect(again).toEqual({ ...first, ms: 0, cached: true, note: REPEAT_NOTE });
        expect(REPEAT_NOTE).toBe('identical call earlier in this question — the same result is repeated; change the arguments or answer');
        expect(repeatedCall(await call(refuser, { value: 1 }))).toMatchObject({ resultStatus: 'refused', cached: true, facts: { seen: 1 } });
        for (const failed of [skippedCall({}, 'one_question'), unknownToolCall({}, 'x', []), await call(echo(), '{')]) {
            expect(() => repeatedCall(failed)).toThrow(TypeError);
        }
    });

    it('a note: the app\'s goes before the result\'s own; none leaves the result as it is', async () => {
        const first = await call(echo(), { value: 2 });
        expect(withNote(first, null)).toBe(first);
        expect(withNote(first, 'the Workbench applied the user\'s choice: x').note).toBe('the Workbench applied the user\'s choice: x');
        expect(withNote(repeatedCall(first), 'the Workbench applied the user\'s choice: x').note)
            .toBe(`the Workbench applied the user's choice: x; ${REPEAT_NOTE}`);
        // An error keeps its note too (a settled window on a call that failed).
        expect(withNote(await call(echo(), '{'), 'n')).toMatchObject({ resultStatus: 'error', note: 'n' });
    });

    it('a question the user answered: the answer the model reads, and the choice beside it', () => {
        const request = normalizeAskRequest({ question: 'Which shell?', options: [{ label: 'Whole shell', description: '2.21–3.23 Å' }, 'First peak'], default: 'Whole shell' });
        expect(questionResult(request, { answer: 'First peak', source: 'option', index: 1 }, { args: { a: 1 }, ms: 12 })).toEqual({
            resultStatus: 'ok', args: { a: 1 }, data: { question: 'Which shell?', answer: 'First peak' },
            checks: [{ level: 'ok', text: 'The user chose "First peak". Use exactly this choice and say in the answer which one you used.' }],
            error: null, ms: 12, cached: false, note: null, choice: { source: 'option', index: 1, confirmed: null }
        });
        expect(questionResult(request, { source: 'skipped' }, { confirmed: false })).toMatchObject({
            data: { answer: 'Whole shell', skipped: true }, choice: { source: 'skipped', index: null, confirmed: false }
        });
        expect(questionResult(request, { answer: ' 2.3 to 3 ', source: 'other', index: null }, { confirmed: true }))
            .toMatchObject({ data: { answer: '2.3 to 3', typed: true }, choice: { source: 'other', index: null, confirmed: true } });
        expect(() => questionResult(request, { source: 'option', index: 0 }, { confirmed: 'yes' })).toThrow(TypeError);
    });

    it('an interactive tool\'s reply passes through as the choice (confirmed: null); never another tool\'s, a malformed one or a refusal\'s', async () => {
        const asking = (choice, checks = [{ level: 'ok', text: 'answered' }]) => echo({
            name: 'asking', interactive: true, run: async () => ({ data: { answer: 'A' }, checks, choice })
        });
        expect((await call(asking({ source: 'option', index: 0 }), { value: 1 })).choice).toEqual({ source: 'option', index: 0, confirmed: null });
        expect((await call(asking({ source: 'skipped' }), { value: 1 })).choice).toEqual({ source: 'skipped', index: null, confirmed: null });
        expect(Object.keys(await call(asking({ source: 'guess', index: 0 }), { value: 1 }))).not.toContain('choice');
        expect(Object.keys(await call(asking({ source: 'option', index: 0 }, REFUSAL), { value: 1 }))).not.toContain('choice');
        const plain = echo({ run: async () => ({ data: {}, checks: [], choice: { source: 'option', index: 0 } }) });
        expect(Object.keys(await call(plain, { value: 1 }))).not.toContain('choice');
    });
});

describe('the constructors refuse a result that breaks an invariant', () => {
    it('a known status; an error exactly for error, unknown_tool and skipped; a fail check exactly for refused', () => {
        expect(() => toolResult({ resultStatus: 'failed' })).toThrow(/resultStatus must be one of/);
        expect(() => toolResult({ resultStatus: 'error' })).toThrow(/names its error/);
        expect(() => toolResult({ resultStatus: 'skipped', error: '' })).toThrow(/names its error/);
        expect(() => toolResult({ error: 'boom' })).toThrow(/has no error/);
        expect(() => toolResult({ resultStatus: 'refused' })).toThrow(/carries the tool's fail check/);
        expect(() => toolResult({ checks: REFUSAL })).toThrow(/it is 'refused'/);
        expect(toolResult({ resultStatus: 'error', error: 'x', checks: [{ level: 'note', text: 'ignored unknown argument "a"' }] }).checks).toHaveLength(1);
    });

    it('no data or facts with an error; a repeat only of ok or refused; a choice only on ok, well formed; a note that says something', () => {
        expect(() => toolResult({ resultStatus: 'error', error: 'x', data: {} })).toThrow(/has no data/);
        expect(() => toolResult({ resultStatus: 'unknown_tool', error: 'x', facts: {} })).toThrow(/has no facts/);
        expect(() => toolResult({ resultStatus: 'skipped', error: 'x', cached: true })).toThrow(/repeated from the cache/);
        expect(() => toolResult({ resultStatus: 'refused', checks: REFUSAL, choice: { source: 'option', index: 0, confirmed: null } })).toThrow(/has no choice/);
        expect(() => toolResult({ choice: { source: 'typed', index: null, confirmed: null } })).toThrow(/choice.source/);
        expect(() => toolResult({ choice: { source: 'option', index: 1.5, confirmed: null } })).toThrow(/choice.index/);
        expect(() => toolResult({ note: '' })).toThrow(/note/);
        // The defaults: a step that ran and answered.
        expect(toolResult()).toEqual({ resultStatus: 'ok', args: {}, data: null, checks: [], error: null, ms: 0, cached: false, note: null });
        expect(Object.isFrozen(toolResult({ choice: { source: 'other', index: null, confirmed: true } }).choice)).toBe(true);
    });
});

// One turn that reaches every status, a repeat, a note, the model's
// question and the Workbench's (a stub clarifier), read off the loop's
// events, the messages it sent and what the verification recorded.
describe('the invariant table, over one turn', () => {
    const scripted = (replies) => {
        const requests = [];
        const stream = (intent) => {
            requests.push({ ...intent, messages: intent.messages.map((message) => ({ ...message })) });
            const reply = replies[Math.min(requests.length - 1, replies.length - 1)];
            return (async function* fragments() {
                if (reply.text) yield { content: reply.text };
                for (const [index, [name, args]] of (reply.calls || []).entries()) {
                    yield { toolCall: { index, id: `c${requests.length}_${index}`, name, arguments: typeof args === 'string' ? args : JSON.stringify(args) } };
                }
                yield { finishReason: reply.calls ? 'tool_calls' : 'stop' };
            }());
        };
        return { stream, requests };
    };
    const clarifier = {
        // The Workbench asks before the turn (a stub of clarify.js's item).
        beforeTurn: () => ({
            kind: 'neighbours',
            request: { question: 'Which neighbours?', options: [{ label: 'A' }, { label: 'B' }], allowOther: true, defaultIndex: 0 },
            choices: [{ label: 'A', apply: () => 'the user chose A.', isDefault: true, confirmed: null }, { label: 'B', apply: () => 'the user chose B.', confirmed: null }]
        }),
        beforeCall: () => null
    };
    const asked = { name: 'ask_user', interactive: true, description: 'Ask.', parameters: { type: 'object', properties: { question: { type: 'string' } } },
        run: async (args, { askUser }) => { const reply = await askUser({ question: args.question, options: [{ label: 'X' }, { label: 'Y' }], defaultIndex: 0 }); return { data: { answer: reply.answer }, checks: [{ level: 'ok', text: 'answered' }], choice: reply }; } };

    const turn = async () => {
        const model = scripted([
            { calls: [['echo', { value: 1 }], ['refuser', { value: 1 }], ['echo', '{'], ['nope', {}], ['echo', { value: 5 }]] },
            { calls: [['echo', { value: 1 }], ['ask', { question: 'Which?' }]] },
            { text: 'Done.\nOutcome: done.' }
        ]);
        const verification = createVerification();
        const recorded = [];
        const spy = { ...verification, record: (step, show) => { const sent = verification.record(step, show); recorded.push({ step, sent }); return sent; } };
        const events = [];
        for await (const event of runCopilotTurn({
            model: model.stream,
            messages: [{ role: 'system', content: 'SYS' }, { role: 'user', content: 'question' }],
            tools: [echo(), refuser, { ...asked, name: 'ask' }],
            toolContext: { askUser: async () => ({ answer: 'X', source: 'option', index: 0 }) },
            clarifier,
            verification: spy
        })) events.push(event);
        return { events, recorded, requests: model.requests };
    };

    it('reaches every status, a repeat, the model\'s question and the Workbench\'s', async () => {
        const { events } = await turn();
        const ends = events.filter((event) => event.type === 'tool_end');
        expect(ends.map((end) => [end.name, end.resultStatus, end.cached])).toEqual([
            ['ask_user', 'ok', false], ['echo', 'ok', false], ['refuser', 'refused', false], ['echo', 'error', false],
            ['nope', 'unknown_tool', false], ['echo', 'skipped', false], ['echo', 'ok', true], ['ask', 'ok', false]
        ]);
    });

    it('1. every tool_end has exactly one result status; a step has none only while running or once stopped', async () => {
        const { events } = await turn();
        const ends = events.filter((event) => event.type === 'tool_end');
        ends.forEach((end) => expect(RESULT_STATUSES).toContain(end.resultStatus));
        // The dock's fold: null while running, the status at its end, null when Stop closes a running step.
        let folded = newAssistantTurn(0);
        for (const event of events) {
            folded = applyTurnEvent(folded, event, 0);
            if (event.type === 'tool_start') expect(folded.items.at(-1)).toMatchObject({ status: 'running', resultStatus: null });
            if (event.type === 'tool_end') expect(folded.items.find((item) => item.id === event.id).resultStatus).toBe(event.resultStatus);
        }
        const running = applyTurnEvent(newAssistantTurn(0), { type: 'tool_start', id: 't', name: 'echo', label: 'Echo', args: {} }, 0);
        expect(closeTurn(running, 'stopped').items[0]).toMatchObject({ status: 'stopped', resultStatus: null });
    });

    it('2. an error exactly for error, unknown_tool and skipped — and then no ref, no data and its facts dropped', async () => {
        const { events, recorded } = await turn();
        const ends = events.filter((event) => event.type === 'tool_end' && event.id.startsWith('c'));
        ends.forEach((end) => {
            const failed = ERROR_STATUSES.includes(end.resultStatus);
            expect(Boolean(end.error), end.id).toBe(failed);
            expect(end.ref === null, end.id).toBe(failed);
            if (failed) expect(end.data, end.id).toBeNull();
        });
        recorded.forEach(({ step, sent }) => {
            if (ERROR_STATUSES.includes(step.resultStatus)) expect(Object.keys(step)).not.toContain('facts');
            expect(sent.ref === null).toBe(Boolean(step.error));
        });
        // A refusal's facts are kept: the answer checks read them.
        expect(recorded.find(({ step }) => step.resultStatus === 'refused').step.facts).toEqual({ seen: 1 });
    });

    it('3. the step\'s status is its severity (severityOf); 4. it is an analysis exactly when it ran (ran)', async () => {
        const { events, recorded } = await turn();
        events.filter((event) => event.type === 'tool_end').forEach((end) => expect(end.status).toBe(severityOf(end)));
        recorded.forEach(({ step, sent }) => expect(sent.analysis).toBe(ran(step) && !step.interactive && !step.reading && !step.view));
        expect(recorded.map(({ sent }) => sent.analysis)).toEqual([true, true, true, false, false, true, false]);
    });

    it('5. a repeat is ok or refused', async () => {
        const { events } = await turn();
        events.filter((event) => event.type === 'tool_end' && event.cached).forEach((end) => expect(['ok', 'refused']).toContain(end.resultStatus));
    });

    it('6. the message bytes depend only on the ref, the error, the note, the checks and the data, in that order; the status and the choice are never sent', async () => {
        const { events, requests } = await turn();
        const messages = new Map(requests.at(-1).messages.filter((message) => message.role === 'tool').map((message) => [message.tool_call_id, message.content]));
        const ends = events.filter((event) => event.type === 'tool_end' && event.id.startsWith('c'));
        ends.forEach((end) => {
            const expected = JSON.stringify({
                ...(end.ref && !end.error ? { ref: end.ref } : {}),
                ...(end.error ? { error: end.error } : {}),
                ...(end.note ? { note: end.note } : {}),
                ...(end.checks.length ? { checks: end.checks } : {}),
                ...(!end.error ? { data: end.data } : {})
            });
            expect(messages.get(end.id.split('#')[0])).toBe(expected);
        });
        const sent = [...messages.values()].join('\n');
        expect(sent).not.toMatch(/resultStatus|result_status|"choice"|confirmed/);
        // The same message from results that differ in nothing the model reads.
        const base = await call(echo(), { value: 2 });
        const message = toolMessage(base, { budget: 3000, ref: 'r9' }).content;
        expect(toolMessage({ ...base, resultStatus: 'refused', choice: { source: 'option' }, ms: 99, cached: true, facts: { x: 1 }, args: 'other' }, { budget: 3000, ref: 'r9' }).content).toBe(message);
    });

    it('7. a Workbench question takes no ref: the first call after it is r1', async () => {
        const { events, recorded } = await turn();
        const [question, first] = events.filter((event) => event.type === 'tool_end');
        expect(question).toMatchObject({ name: 'ask_user', ref: null, analysis: false, note: '(The Workbench asked: "Which neighbours?" — the user chose A.)', choice: { source: 'option', index: 0, confirmed: null } });
        expect(first.ref).toBe('r1');
        // Not recorded with the verification at all.
        expect(recorded.map(({ step }) => step.name)).not.toContain('ask_user');
    });
});

describe('the tool message and its compaction', () => {
    it('the ref (with no error), the error, the note, the checks, then the data fitted to the budget less the rest', () => {
        const ok = toolResult({ data: { v: 1 }, checks: [{ level: 'ok', text: 'fine' }], note: 'a note' });
        expect(toolMessage(ok, { budget: 3000, ref: 'r1' })).toEqual({
            content: '{"ref":"r1","note":"a note","checks":[{"level":"ok","text":"fine"}],"data":{"v":1}}', data: { v: 1 }
        });
        const failed = toolResult({ resultStatus: 'error', error: 'boom', checks: [{ level: 'note', text: 'n' }], note: 'a note' });
        expect(toolMessage(failed, { budget: 3000, ref: 'r1' })).toEqual({ content: '{"error":"boom","note":"a note","checks":[{"level":"note","text":"n"}]}', data: null });
        // No ref, no checks: the data alone; the data null is still sent.
        expect(toolMessage(toolResult(), { budget: 3000 }).content).toBe('{"data":null}');
        expect(() => toolMessage(ok, {})).toThrow(/budget/);
    });

    it('fits large results to the budget without dropping the checks', () => {
        const big = toolResult({ data: { rows: Array.from({ length: 500 }, (_, index) => ({ index, value: index * 1.5 })) }, checks: [{ level: 'warn', text: 'big' }] });
        const { content, data } = toolMessage(big, { budget: 1200, ref: 'r1' });
        expect(content.length).toBeLessThanOrEqual(1300);
        const parsed = JSON.parse(content);
        expect(parsed.checks).toEqual([{ level: 'warn', text: 'big' }]);
        expect(parsed.data).toEqual(data);
        expect(parsed.data.rows.length + parsed.data.rows_omitted).toBe(500);
    });

    it('compacted: the ref, the checks and the headline numbers kept; an error kept with its checks; never twice', () => {
        const full = toolMessage(toolResult({ data: { value: 1, rows: Array.from({ length: 200 }, (_, index) => index) }, checks: [{ level: 'ok', text: 'round 1' }], note: 'n' }), { budget: 900, ref: 'r1' }).content;
        const compact = compactToolMessage(full);
        expect(JSON.parse(compact)).toEqual({
            ref: 'r1',
            compacted: 'this earlier result was shortened to fit the context; its checks and headline numbers are kept',
            checks: [{ level: 'ok', text: 'round 1' }],
            key: { value: 1 }
        });
        expect(compactToolMessage(compact)).toBe(compact);
        expect(compactToolMessage('{"error":"engine said no","note":"a long note"}'))
            .toBe('{"compacted":"this earlier result was shortened to fit the context; its checks are kept","error":"engine said no"}');
        expect(compactToolMessage('not json')).toBe('{"compacted":"this earlier result was shortened to fit the context"}');
    });
});

describe('tool helpers', () => {
    it('keeps only well-formed, uniquely named tools', () => {
        const ok = echo();
        expect(usableTools([ok, { ...ok }, { ...ok, name: 'Bad Name' }, { name: 'x', description: 'd' }, null])).toEqual([ok]);
    });

    it('the step\'s severity, from the error and the checks: failed > warn > ok', () => {
        expect(severityOf({ error: 'x', checks: [] })).toBe('failed');
        expect(severityOf({ error: null, checks: [{ level: 'fail', text: 'a' }] })).toBe('failed');
        expect(severityOf({ error: null, checks: [{ level: 'warn', text: 'a' }, { level: 'ok', text: 'b' }] })).toBe('warn');
        expect(severityOf({ error: null, checks: [{ level: 'note', text: 'a' }] })).toBe('ok');
    });

    it('ran: not an unknown tool, not a skipped call, and only a result', () => {
        expect(['ok', 'refused', 'error', 'unknown_tool', 'skipped'].map((resultStatus) => ran({ resultStatus }))).toEqual([true, true, true, false, false]);
        expect(ran({})).toBe(false);
        expect(ran(null)).toBe(false);
    });

    it('trims nested arrays and hard-cuts as a last resort, never silently', () => {
        const value = { a: { b: Array.from({ length: 100 }, (_, i) => i) }, c: 'x' };
        const fitted = fitToBudget(value, 120);
        expect(JSON.stringify(fitted).length).toBeLessThanOrEqual(120);
        expect(fitted.a.b_omitted + fitted.a.b.length).toBe(100);
        const cut = JSON.parse(budgetJson({ s: 'y'.repeat(500) }, 100));
        expect(cut.truncated).toBe(true);
        expect(fitToBudget(null, 10)).toBeNull();
    });
});
