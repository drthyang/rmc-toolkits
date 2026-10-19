// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The dock's turn reducer (copilot/turnState.js): the gate's send-backs, the
// answer and the earlier segments with their references filled, a stopped
// turn, and the footer's count of analyses.

import { describe, expect, it } from 'vitest';
import { applyTurnEvent, closeTurn, newAssistantTurn, questionReply, stepSummary } from '../copilot/turnState';

const fold = (events) => events.reduce((turn, event) => applyTurnEvent(turn, event, 0), newAssistantTurn(0));
const tool = (id, end = {}) => [
    { type: 'tool_start', id, name: id, label: id, args: {} },
    ...(end ? [{ type: 'tool_end', id, name: id, status: 'ok', args: {}, data: null, checks: [], error: null, ms: 1, ref: null, ...end }] : [])
];

describe('the turn shows the send-back and the filled answer', () => {
    it('keeps each gate\'s notice, withdraws the sent-back text, and shows the answer with its numbers filled', () => {
        const turn = fold([
            { type: 'round', round: 0 },
            { type: 'text', delta: 'Ti forms a superstructure.' },
            { type: 'gate', id: 'order-between-cells', notice: 'The answer concluded order between cells without analyze_modulations; the app asked for it first.' },
            { type: 'retract' },
            { type: 'round', round: 1 },
            { type: 'text', delta: 'Ordered at {{r1.x}}×.\nOutcome: ordered.' },
            {
                type: 'done', text: 'Ordered at 498×.\nOutcome: ordered.', rawText: 'Ordered at {{r1.x}}×.\nOutcome: ordered.',
                outcome: { status: 'achieved', note: '1 analysis ran.', body: 'Ordered at 498×.', claimed: null, computed: true },
                stoppedBy: 'answer', audit: [], references: { written: 1, filled: 1, unknown: [] }, gates: []
            }
        ]);
        expect(turn.gates).toEqual([{ id: 'order-between-cells', notice: expect.stringMatching(/^The answer concluded order between cells/), text: '' }]);
        expect(turn.items.filter((item) => item.kind === 'text').map((item) => item.text)).toEqual(['', 'Ordered at 498×.']);
        expect(turn.references).toEqual({ written: 1, filled: 1, unknown: [] });
    });

    it('an ended or stopped turn shows every segment with its references filled, and keeps an answer the gate withdrew', () => {
        const turn = fold([
            { type: 'round', round: 0 },
            { type: 'tool_start', id: 't1', name: 'echo', label: 'Echo', args: {} },
            { type: 'tool_end', id: 't1', name: 'echo', status: 'ok', args: {}, data: { value: 4 }, checks: [], error: null, ms: 1, ref: 'r1', analysis: true },
            { type: 'round', round: 1 },
            { type: 'text', delta: 'The value is {{r1.value}}; checking more.' }
        ]);
        expect(closeTurn(turn, 'stopped').items.find((item) => item.kind === 'text' && item.text).text).toBe('The value is 4; checking more.');
        const gated = fold([
            { type: 'round', round: 0 },
            { type: 'text', delta: 'Ti forms a superstructure.' },
            { type: 'gate', id: 'order-between-cells', notice: 'n', text: 'Ti forms a superstructure.' },
            { type: 'retract' }
        ]);
        expect(closeTurn(gated, 'stopped').text).toBe('Ti forms a superstructure.\n\n(Interrupted before the answer was finished.)');
    });

    // Fixed with the verification (candidate #4): the earlier segments were
    // filled without the context, so a {{ctx.…}} stayed raw there.
    it('fills an earlier segment\'s {{ctx.…}} from the context as the model read it (the start event\'s shown)', () => {
        const events = [
            { type: 'start', shown: { run: 'Demo', symmetry: { noise_floor_A: 0.031 } } },
            { type: 'round', round: 0 },
            { type: 'text', delta: 'The noise floor is {{ctx.symmetry.noise_floor_A}} Å; checking the sites.' },
            ...tool('list_sites', { analysis: true }),
            { type: 'round', round: 1 }
        ];
        const stopped = closeTurn(fold(events), 'stopped');
        expect(stopped.items.find((item) => item.kind === 'text' && item.text).text).toBe('The noise floor is 0.031 Å; checking the sites.');
        const ended = fold([...events, { type: 'text', delta: 'Done.' }, { type: 'done', text: 'Done.', outcome: null, stoppedBy: 'answer', audit: [], references: null }]);
        expect(ended.items.filter((item) => item.kind === 'text').map((item) => item.text)).toEqual(['The noise floor is 0.031 Å; checking the sites.', 'Done.']);
        // Without a start event (an older caller), it stays as written.
        expect(closeTurn(fold(events.slice(1)), 'stopped').items.find((item) => item.kind === 'text' && item.text).text)
            .toBe('The noise floor is {{ctx.symmetry.noise_floor_A}} Å; checking the sites.');
    });
});

// Behaviour change (candidate #4): the footer counts the steps that ended as
// analyses (tool_end's `analysis`, the verification's isAnalysis, which the
// computed Outcome counts). Before, it counted every tool step that was not a
// question or a skill read: a view tool, an unknown tool, a call skipped by
// the budget and a step stopped mid-run counted too.
describe('the footer counts analyses as the Outcome does', () => {
    it('counts tool_end.analysis: never a view, an unknown or skipped call, a question, a skill read or a stopped step', () => {
        const turn = fold([
            { type: 'round', round: 0 },
            ...tool('find_gr_peaks', { analysis: true, status: 'warn' }),
            ...tool('compute_bond_angles', { analysis: true, status: 'failed', error: null }),
            ...tool('highlight_cards', { analysis: false }),
            ...tool('analyse_sites', { analysis: false, status: 'failed', error: 'unknown tool "analyse_sites"' }),
            ...tool('list_sites', { analysis: false, status: 'failed', error: 'skipped: at most 4 tool calls per step' }),
            { type: 'tool_start', id: 'q', name: 'ask_user', label: 'Question for you', args: {}, interactive: true },
            { type: 'tool_end', id: 'q', name: 'ask_user', status: 'ok', args: {}, data: null, checks: [], error: null, ms: 1, ref: 'r3', analysis: false },
            { type: 'tool_start', id: 's', name: 'read_skill', label: 'read_skill', args: {}, reading: true },
            { type: 'tool_end', id: 's', name: 'read_skill', status: 'ok', args: {}, data: null, checks: [], error: null, ms: 1, ref: 'r4', analysis: false },
            ...tool('analyze_clusters', null)
        ]);
        expect(stepSummary(turn.items)).toBe('2 analyses · 1 warning · 1 failed');
        expect(stepSummary(closeTurn(turn, 'stopped').items)).toBe('2 analyses · 1 warning · 1 failed');
        expect(stepSummary(fold(tool('highlight_cards', { analysis: false })).items)).toBe('');
        expect(stepSummary(fold(tool('list_sites', { analysis: true })).items)).toBe('1 analysis');
    });
});

// The tool transcript (candidate #11): a step keeps how its call ended, the
// app's note and the user's choice, as tool_end carries them; a question's
// card and the screen reader read its state from them.
describe('a step keeps its result: the status, the note, the choice', () => {
    it('null while it runs and once stopped; tool_end\'s at its end', () => {
        const running = fold(tool('list_sites', null));
        expect(running.items[0]).toMatchObject({ status: 'running', resultStatus: null, note: null, choice: null });
        expect(closeTurn(running, 'stopped').items[0]).toMatchObject({ status: 'stopped', resultStatus: null });
        const ended = fold(tool('list_sites', { resultStatus: 'ok', cached: true, note: 'identical call earlier in this question', analysis: true }));
        expect(ended.items[0]).toMatchObject({ status: 'ok', resultStatus: 'ok', cached: true, note: 'identical call earlier in this question', choice: null });
    });

    it('a question\'s state, from its result and the user\'s choice', () => {
        const question = (end) => fold([
            { type: 'tool_start', id: 'q', name: 'ask_user', label: 'Question for you', args: { question: 'Which?' }, interactive: true },
            ...(end ? [{ type: 'tool_end', id: 'q', name: 'ask_user', status: 'ok', resultStatus: 'ok', args: {}, checks: [], error: null, ms: 1, ref: null, ...end }] : [])
        ]).items[0];
        expect(questionReply(question(null))).toEqual({ state: 'waiting', text: '' });
        expect(questionReply(closeTurn({ ...newAssistantTurn(0), items: [question(null)] }, 'stopped').items[0])).toEqual({ state: 'stopped', text: 'Not answered' });
        expect(questionReply(question({ status: 'failed', resultStatus: 'skipped', error: 'one question at a time', data: null })))
            .toEqual({ state: 'not_asked', text: 'one question at a time' });
        expect(questionReply(question({ status: 'failed', resultStatus: 'error', error: 'the user cannot be asked here', data: null })).state).toBe('not_asked');
        expect(questionReply(question({ data: { answer: 'Whole shell', skipped: true }, choice: { source: 'skipped', index: null, confirmed: false } })))
            .toEqual({ state: 'skipped', text: 'Skipped: the Copilot uses its default, Whole shell.' });
        expect(questionReply(question({ data: { answer: null, skipped: true }, choice: { source: 'skipped', index: null, confirmed: null } })))
            .toEqual({ state: 'skipped', text: 'Skipped.' });
        expect(questionReply(question({ data: { answer: 'First peak' }, choice: { source: 'option', index: 1, confirmed: null } })))
            .toEqual({ state: 'answered', text: 'First peak' });
        // The words the model read play no part: "skipped" in the data alone is no skip.
        expect(questionReply(question({ data: { answer: 'A', skipped: true }, choice: { source: 'option', index: 0, confirmed: null } })).state).toBe('answered');
    });
});
