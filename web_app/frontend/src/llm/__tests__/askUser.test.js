// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { MAX_OPTIONS, answerResult, buildAskUserTool, normalizeAskRequest, replyChoice } from '../copilot/askUser';
import { validateToolArguments } from '../copilot/toolCalls';

describe('normalizeAskRequest', () => {
    it('keeps a question, its options with their meaning, and the default', () => {
        expect(normalizeAskRequest({
            question: '  Which   shell? ',
            options: [{ label: 'Whole shell', description: '2.33–2.96 Å' }, { label: 'First peak', description: '2.33–2.71 Å' }],
            default: 'whole SHELL'
        })).toEqual({
            question: 'Which shell?',
            options: [{ label: 'Whole shell', description: '2.33–2.96 Å' }, { label: 'First peak', description: '2.33–2.71 Å' }],
            allowOther: true,
            defaultIndex: 0
        });
    });

    it('reads plain strings, splits a description written into the label, and drops repeats', () => {
        const request = normalizeAskRequest({
            question: 'Which run?',
            options: ['250 K — the Demo', '5 K', '5 k', ''],
            allow_other: false,
            default: '2'
        });
        expect(request.options).toEqual([{ label: '250 K', description: 'the Demo' }, { label: '5 K' }]);
        expect(request.allowOther).toBe(false);
        expect(request.defaultIndex).toBe(1);
        // A range is not a description: "2.33–2.96 Å" keeps its dash.
        expect(normalizeAskRequest({ question: 'q', options: ['2.33–2.96 Å', '2.33–2.71 Å'] }).options[0]).toEqual({ label: '2.33–2.96 Å' });
    });

    it('never splits a spaced range into a label and a description (review 2026-10-12 #9)', () => {
        const { options } = normalizeAskRequest({ question: 'Which window?', options: ['2.2 – 2.7 Å', 'Whole shell 2.2 - 3.2 Å — CN 6', 'First peak — 2.2 – 2.7 Å'] });
        expect(options).toEqual([
            { label: '2.2 – 2.7 Å' },
            { label: 'Whole shell 2.2 - 3.2 Å', description: 'CN 6' },
            { label: 'First peak', description: '2.2 – 2.7 Å' }
        ]);
    });

    it('ignores a default that names no option', () => {
        expect(normalizeAskRequest({ question: 'q', options: ['a', 'b'], default: 'c' }).defaultIndex).toBeNull();
    });

    it('refuses what the model can fix, in words it can read', () => {
        expect(() => normalizeAskRequest({ question: ' ', options: ['a', 'b'] })).toThrow(/question is empty/);
        expect(() => normalizeAskRequest({ question: 'q', options: ['a', 'A'] })).toThrow(/2–4 distinct options/);
        expect(() => normalizeAskRequest({ question: 'q', options: ['a', 'b', 'c', 'd', 'e'] })).toThrow(`at most ${MAX_OPTIONS}`);
    });
});

describe('answerResult', () => {
    const request = normalizeAskRequest({
        question: 'Which shell?',
        options: [{ label: 'Whole shell', description: '2.33–2.96 Å' }, 'First peak'],
        default: 'Whole shell'
    });

    it('returns the chosen option with its meaning, to be used exactly', () => {
        expect(answerResult(request, { answer: 'Whole shell', source: 'option', index: 0 })).toEqual({
            data: { question: 'Which shell?', answer: 'Whole shell', means: '2.33–2.96 Å' },
            checks: [{ level: 'ok', text: 'The user chose "Whole shell" (2.33–2.96 Å). Use exactly this choice and say in the answer which one you used.' }]
        });
    });

    it('returns a typed answer as typed', () => {
        expect(answerResult(request, { answer: ' 2.3–3.1 Å ', source: 'other', index: null }).data)
            .toEqual({ question: 'Which shell?', answer: '2.3–3.1 Å', typed: true });
    });

    it('turns a skip into the default, said as an assumption', () => {
        const skipped = answerResult(request, { source: 'skipped' });
        expect(skipped.data).toEqual({ question: 'Which shell?', answer: 'Whole shell', skipped: true, means: '2.33–2.96 Å' });
        expect(skipped.checks[0].text).toMatch(/say in the answer that you assumed it/);
        const noDefault = answerResult(normalizeAskRequest({ question: 'q', options: ['a', 'b'] }), { source: 'skipped' });
        expect(noDefault.data.answer).toBeNull();
        expect(noDefault.checks[0].text).toMatch(/most standard option/);
    });
});

describe('the reply as the step\'s choice (the tool transcript)', () => {
    it('keeps where the reply came from and the option\'s index; a reply of no known source reads as an option', () => {
        expect(replyChoice({ answer: 'A', source: 'option', index: 1 })).toEqual({ source: 'option', index: 1 });
        expect(replyChoice({ answer: ' typed ', source: 'other', index: null })).toEqual({ source: 'other', index: null });
        expect(replyChoice({ source: 'skipped' })).toEqual({ source: 'skipped', index: null });
        expect(replyChoice(undefined)).toEqual({ source: 'option', index: null });
    });

    it('ask_user returns it beside the answer the model reads', async () => {
        const reply = { answer: null, source: 'skipped', index: null };
        const result = await buildAskUserTool().run({ question: 'Which?', options: ['a', 'b'], default: 'a' }, { askUser: async () => reply });
        expect(result).toEqual({ ...answerResult(normalizeAskRequest({ question: 'Which?', options: ['a', 'b'], default: 'a' }), reply), choice: { source: 'skipped', index: null } });
    });
});

describe('the ask_user tool', () => {
    it('is interactive, not a view tool, and its schema takes strings or objects as options', () => {
        const tool = buildAskUserTool();
        expect(tool.interactive).toBe(true);
        expect(tool.view).toBeUndefined();
        const strings = validateToolArguments(tool.parameters, { question: 'q', options: ['a', 'b'] });
        expect(strings.ok).toBe(true);
        expect(strings.value.options).toEqual([{ label: 'a' }, { label: 'b' }]);
        expect(validateToolArguments(tool.parameters, { question: 'q', options: ['a'] }).errors)
            .toEqual(['options must hold at least 2 item(s), got 1']);
    });
});
