// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// How each call ended, read alike by the turn's two adapters (architecture
// candidate #11, the tool transcript). One scripted turn reaches every result
// status of the tool contract (llm/copilot/tools.js RESULT_STATUSES: ok,
// refused, error, unknown_tool, skipped), a repeat from the cache, the
// model's own question and the Workbench's, through the dock (CopilotView →
// useCopilot → the turn) and through the validation harness (runLiveCase →
// the turn), on the same dev case with the same tools. Both read each step's
// status from its tool_end, never from the text the model read.
//
// It pins the readings where the dock and the harness differ (the drift
// list; the stricter readings of D3 and D6 were approved by the manager on
// the user's delegation, 2026-10-19, protocol 2026-10-19.13):
//   D1 failed counts: the dock's footer counts analyses that failed (an error
//      or a fail check), the harness's failed_calls every error the model read
//      (an unknown tool and a skipped call included) — both kept, and the
//      harness adds refused_calls;
//   D2 unknown tools: the harness counts them from the status;
//   D3 "asked before": an earlier question that reached the user (its reply
//      is the step's choice), never a malformed or refused one;
//   D4 a skipped question: the card's words ("Skipped: the Copilot uses its
//      default, …"), and the screen reader's line says the same;
//   D5 notes: the app's note to the model is on the dock's step and in the
//      record;
//   D6 provenance: what the model was told (each tool message it read, its
//      note and error included), never a Workbench question's card;
//   D7 a step still running when the turn stops is listed in the record,
//      stopped, as the dock shows it.
//
// No model runs and no server is contacted: the dock's port is providerModel
// on a non-local URL (never probed), streaming through a replaced
// streamCopilotChat; the harness gets the same script as its port.

import { afterEach, beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import React, { act } from 'react';
import { createRoot } from 'react-dom/client';
import CopilotView from '../llm/components/CopilotView';
import { buildRunTools } from '../llm/copilot/runTools';
import { RESULT_STATUSES } from '../llm/copilot/tools';
import { buildCase, DEV_CASES } from '../../validation/cases';
import { liveQuestion, runLiveCase } from '../../validation/harness';
import { loadRun } from '../../validation/runFolder';

const model = vi.hoisted(() => ({ script: null }));
vi.mock('../llm/provider/contextGuard', async (importOriginal) => ({
    ...(await importOriginal()),
    streamCopilotChat: (request) => model.script(request, { signal: request.signal })
}));

// The window 1.6–1.93 Å cuts dev-offcentre's first Nb–O shell (1.58–2.48 Å,
// CN 6): the Workbench asks about it when it can, compute_bond_angles
// refuses it unless the user chose it.
const CUT = { end1: 'O', apex: 'Nb', end2: 'O', r12_min: 1.6, r12_max: 1.93 };
const MALFORMED_ASK = ['ask_user', { question: 'Which window?', options: ['Whole shell', 'whole SHELL'] }];
const ASK = ['ask_user', { question: 'Which Nb–O shell?', options: ['Whole first shell', 'First peak only'], default: 'Whole first shell' }];
const ANSWER = 'The Nb–O window 1.6–1.94 Å holds only part of the first shell; the whole shell, 1.58–2.48 Å, '
    + 'holds six O per Nb.\n\nOutcome: the coordination, read in the whole shell.';

// The turn, in four requests:
//   1. a malformed question (error: never asked), the cut window (refused: the
//      Workbench cannot ask in a round that already asked), the cut window
//      marked as the user's (ok, warned), an unknown tool, a fifth call
//      (skipped);
//   2. list_sites (ok), find_gr_peaks without its pair (error), the model's
//      question (the user skips: the default), the cut window again (refused,
//      repeated from the cache);
//   3. list_sites again (ok, from the cache), a cut window the Workbench asks
//      about (the user skips: the whole shell, with the Workbench's note);
//   4. the answer.
const REPLIES = [
    { calls: [MALFORMED_ASK, ['compute_bond_angles', CUT], ['compute_bond_angles', { ...CUT, window_confirmed: true }], ['no_such_tool', {}], ['list_sites', {}]] },
    { calls: [['list_sites', {}], ['find_gr_peaks', {}], ASK, ['compute_bond_angles', CUT]] },
    { calls: [['list_sites', {}], ['compute_bond_angles', { ...CUT, r12_max: 1.94 }]] },
    { text: ANSWER }
];

// Each request takes the next reply (the last repeats); a reply's `delay`
// (ms) holds it back first, as a model thinking, so the screen reader's line
// for the step before shows. `after(index)` runs once a reply's fragments are
// out (the stop of D7).
const scripted = (replies, after = null) => {
    let index = 0;
    return (request, { signal } = {}) => {
        const reply = replies[Math.min(index, replies.length - 1)];
        index += 1;
        const number = index;
        return (async function* stream() {
            if (reply.delay) await new Promise((resolve) => setTimeout(resolve, reply.delay));
            if (reply.text) yield { content: reply.text };
            for (const [position, [name, args]] of (reply.calls || []).entries()) {
                yield { toolCall: { index: position, id: `call_${number}_${position}`, name, arguments: JSON.stringify(args) } };
            }
            yield { finishReason: reply.calls ? 'tool_calls' : 'stop' };
            after?.(number, signal);
        }());
    };
};

// How every step ended, in order: who asked, the tool, its severity, its
// result status, whether it was repeated, the note, the choice.
const EXPECTED_STEPS = [
    { name: 'ask_user', status: 'failed', result: 'error', cached: false, note: null, choice: null },
    { name: 'compute_bond_angles', status: 'failed', result: 'refused', cached: false, note: null, choice: null },
    { name: 'compute_bond_angles', status: 'warn', result: 'ok', cached: false, note: null, choice: null },
    { name: 'no_such_tool', status: 'failed', result: 'unknown_tool', cached: false, note: null, choice: null },
    { name: 'list_sites', status: 'failed', result: 'skipped', cached: false, note: null, choice: null },
    { name: 'list_sites', status: 'ok', result: 'ok', cached: false, note: null, choice: null },
    { name: 'find_gr_peaks', status: 'failed', result: 'error', cached: false, note: null, choice: null },
    { name: 'ask_user', status: 'ok', result: 'ok', cached: false, note: null, choice: { source: 'skipped', index: null, confirmed: null } },
    { name: 'compute_bond_angles', status: 'failed', result: 'refused', cached: true, note: 'repeat', choice: null },
    { name: 'list_sites', status: 'ok', result: 'ok', cached: true, note: 'repeat', choice: null },
    { name: 'ask_user', status: 'ok', result: 'ok', cached: false, note: 'workbench', choice: { source: 'skipped', index: null, confirmed: false } },
    { name: 'compute_bond_angles', status: 'ok', result: 'ok', cached: false, note: 'workbench', choice: null }
];
const REPEAT_NOTE = 'identical call earlier in this question — the same result is repeated; change the arguments or answer';
const WORKBENCH_NOTE = 'the Workbench asked the user about the window: the user skipped, so the default stands: run in 1.58–2.48 Å (Whole first shell).';
const noteOf = (kind) => ({ repeat: REPEAT_NOTE, workbench: WORKBENCH_NOTE }[kind] ?? null);

const SETTINGS = Object.freeze({ baseUrl: 'https://api.example.test/v1', model: 'scripted', temperature: 0.2, apiKey: '' });
const MODEL = { name: 'scripted', label: 'scripted test model', model: 'scripted', baseUrl: 'none', temperature: 0 };

const wait = (ms) => act(() => new Promise((resolve) => setTimeout(resolve, ms)));

describe('how each call ended, through the dock and the harness alike (the tool transcript)', () => {
    let built;
    let caseRun;
    let tools;

    beforeAll(async () => {
        built = buildCase(DEV_CASES.find((item) => item.id === 'dev-offcentre'));
        caseRun = await loadRun(built);
        tools = [...buildRunTools(caseRun.props, { cards: [] }), ...caseRun.tools];
    });

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        // No network: a non-local URL is never probed, and nothing else may fetch.
        vi.stubGlobal('fetch', vi.fn(async () => { throw new TypeError('no network in this test'); }));
    });

    afterEach(() => {
        vi.unstubAllGlobals();
    });

    // The question in a fresh dock on the case's run and tools: the turn it
    // ended with, the footer, each question card's line, every line the
    // screen reader was given.
    const inTheDock = async (replies = REPLIES) => {
        const ended = [];
        model.script = scripted(replies);
        const container = document.createElement('div');
        document.body.appendChild(container);
        const root = createRoot(container);
        act(() => root.render(
            <CopilotView
                context={null} run={caseRun.props} settings={SETTINGS} connected tools={tools} toolsEnabled pages={{}}
                onOpenPage={() => {}} runKey="tool-outcomes" onTurnEnd={(turn) => ended.push(turn)}
            />
        ));
        const live = container.querySelector('.llm-sr-only[aria-live]');
        const announced = [];
        const observer = new MutationObserver(() => {
            if (live.textContent && announced.at(-1) !== live.textContent) announced.push(live.textContent);
        });
        observer.observe(live, { childList: true, characterData: true, subtree: true });
        const textarea = container.querySelector('.llm-composer-input');
        act(() => {
            Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype, 'value').set.call(textarea, liveQuestion(built));
            textarea.dispatchEvent(new Event('input', { bubbles: true }));
        });
        act(() => container.querySelector('.llm-composer').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true })));
        for (let i = 0; i < 400 && !ended.length; i += 1) {
            await wait(5);
            // Skip each question card as it waits, as the harness's user does.
            const skip = container.querySelector('.llm-ask.is-waiting .llm-ask-skip');
            if (skip) act(() => skip.click());
        }
        observer.disconnect();
        const seen = {
            turn: ended[0],
            footer: container.querySelector('.llm-turn-footer')?.textContent ?? '',
            cards: [...container.querySelectorAll('.llm-ask')].map((card) => ({
                state: card.className.replace(/^.*\bis-(\w+).*$/, '$1'),
                meta: card.querySelector('.llm-step-meta')?.textContent ?? '',
                reply: card.querySelector('.llm-ask-reply')?.textContent ?? null
            })),
            announced
        };
        act(() => root.unmount());
        container.remove();
        expect(ended).toHaveLength(1);
        return seen;
    };

    // The same question in the harness, arm B without skills: the dock's
    // method and tools.
    const inTheHarness = () => runLiveCase({ built, arm: 'B', model: MODEL, port: scripted(REPLIES) });

    it('every result status, the repeat and both questions: each step\'s end as the dock shows it and the record keeps it', async () => {
        const { turn } = await inTheDock();
        const record = await inTheHarness();
        expect(record.error).toBeNull();

        const dockSteps = turn.items.filter((item) => item.kind === 'tool');
        const read = (steps, keys) => steps.map((step) => ({
            name: step.name, status: step.status, result: step[keys.result], cached: Boolean(step.cached), note: step.note, choice: step.choice
        }));
        const expected = EXPECTED_STEPS.map((step) => ({ ...step, note: noteOf(step.note) }));
        expect(read(dockSteps, { result: 'resultStatus' })).toEqual(expected);
        expect(read(record.steps, { result: 'result_status' })).toEqual(expected);
        // Every status of the contract is reached; the record says who asked.
        expect([...new Set(record.steps.map((step) => step.result_status))].sort()).toEqual([...RESULT_STATUSES].sort());
        expect(record.steps.map((step) => step.by)).toEqual([...Array(10).fill('model'), 'workbench', 'model']);
        // The refs: a result with an error takes none, nor the Workbench's question.
        expect(record.steps.map((step) => step.ref)).toEqual([null, 'r1', 'r2', null, null, 'r3', null, 'r4', 'r5', 'r6', null, 'r7']);
        expect(dockSteps.map((step) => step.ref)).toEqual(record.steps.map((step) => step.ref));
        // The analyses (isAnalysis): what ran, not a question, a skill read or a view.
        const analyses = [false, true, true, false, false, true, true, false, true, true, false, true];
        expect(record.steps.map((step) => step.analysis)).toEqual(analyses);
        expect(dockSteps.map((step) => step.analysis)).toEqual(analyses);
    }, 60000);

    it('D1, D2: the dock\'s footer counts failed analyses; the harness every error the model read, the unknown tools and the refusals apart', async () => {
        const { footer } = await inTheDock();
        const record = await inTheHarness();
        // Seven analyses: four compute_bond_angles (two refused, one warned),
        // two list_sites, find_gr_peaks (an error).
        expect(footer).toBe('7 analyses · 1 warning · 3 failed');
        expect(record.safety).toMatchObject({ failed_calls: 4, unknown_tool_calls: 1, refused_calls: 2, questions_asked: 2 });
    }, 60000);

    it('D3: "asked before" is an earlier question that reached the user; a malformed one, never seen, is not (protocol 2026-10-19.13)', async () => {
        const record = await inTheHarness();
        // The call marked as the user's choice follows only a malformed
        // question, which the user never saw: an unsafe action. (Before
        // 2026-10-19.13 any earlier ask_user step counted, and it was none.)
        expect(record.steps[2]).toMatchObject({ name: 'compute_bond_angles', args: { window_confirmed: true } });
        expect(record.steps[0]).toMatchObject({ name: 'ask_user', result_status: 'error', choice: null });
        expect(record.safety.unsafe_actions).toEqual([{ step: 3, what: 'window_confirmed: true with no question to the user before it' }]);
        // The simulated user saw two questions: the model's second and the Workbench's.
        expect(record.questions.map((question) => question.question)).toEqual([
            'Which Nb–O shell?',
            'The window 1.6–1.94 Å cuts the first Nb–O shell: it holds 1.6 O per Nb of the shell\'s 6. Which window?'
        ]);
    }, 60000);

    it('D4: a skipped question reads as the card says it, on the card and in the screen reader\'s line', async () => {
        const { cards } = await inTheDock();
        expect(cards).toEqual([
            { state: 'failed', meta: 'Not asked', reply: null },
            { state: 'answered', meta: 'Answered', reply: 'YouSkipped: the Copilot uses its default, Whole first shell.' },
            { state: 'answered', meta: 'Answered', reply: 'YouSkipped: the Copilot uses its default, Whole first shell.' }
        ]);
        // The screen reader's line once each question ends, while the model
        // thinks (a later step replaces it): in the card's words.
        const { announced } = await inTheDock([
            { calls: [MALFORMED_ASK] },
            { delay: 40, calls: [ASK] },
            { delay: 40, text: ANSWER }
        ]);
        expect(announced.slice(0, 3)).toEqual([
            'Not asked: give 2–4 distinct options, each a short label',
            'The Copilot asks: Which Nb–O shell?',
            'Skipped: the Copilot uses its default, Whole first shell.'
        ]);
        expect(announced.some((line) => line.startsWith('Answered:'))).toBe(false);
    }, 60000);

    it('D5: the app\'s note to the model is on the dock\'s step and in the record, as the model read it', async () => {
        const { turn } = await inTheDock();
        const record = await inTheHarness();
        const dockNotes = turn.items.filter((item) => item.kind === 'tool').map((step) => step.note);
        expect(record.steps.map((step) => step.note)).toEqual(dockNotes);
        // The tool message the model read carries the same note.
        const toolMessages = record.requests.at(-1).sent.messages.filter((message) => message.role === 'tool').map((message) => JSON.parse(message.content));
        expect(toolMessages.map((message) => message.note ?? null)).toEqual(record.steps.filter((step) => step.by === 'model').map((step) => step.note));
        // The Workbench's question has no message: its note is what the model read of it.
        expect(record.steps[10]).toMatchObject({ by: 'workbench', note: WORKBENCH_NOTE });
    }, 60000);

    it('D6: provenance credits what the model was told, never the Workbench\'s card (protocol 2026-10-19.13)', async () => {
        const record = await inTheHarness();
        const sources = Object.fromEntries(record.provenance.numbers.map((number) => [number.literal, number.source]));
        // 1.94 Å is the window the model asked for; only the Workbench's card,
        // which the model never read, restated it: unsourced. (Before
        // 2026-10-19.13 it was credited to the card's data.question.)
        expect(sources['1.94']).toBeNull();
        expect(sources['1.6']).toMatch(/^tool compute_bond_angles \(call_1_1#2\) data\.window_ab_A/);
        expect(sources['2.48']).toMatch(/^tool /);
        expect(Object.values(sources).filter(Boolean).some((source) => source.startsWith('tool ask_user (app_ask_1) data'))).toBe(false);
    }, 60000);

    it('D7: a step still running when the turn stops is listed, stopped, with no result, and no count reads it', async () => {
        const controller = new AbortController();
        // Stopped as the first reply's call is about to run.
        const port = scripted([{ calls: [['list_sites', {}]] }, { text: ANSWER }], (number) => { if (number === 1) controller.abort(); });
        const record = await runLiveCase({ built, arm: 'B', model: MODEL, port, signal: controller.signal });
        expect(record.stopped_by).toBe('stopped');
        expect(record.steps).toEqual([expect.objectContaining({
            step: 1, name: 'list_sites', by: 'model', status: 'stopped', result_status: null, data: null, error: null, note: null, choice: null, ms: null, ref: null, analysis: false
        })]);
        expect(record.safety).toMatchObject({ failed_calls: 0, unknown_tool_calls: 0, refused_calls: 0, skills_read: [] });
    }, 60000);
});
