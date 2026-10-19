// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The method switch's text rules (copilot/method.js, plan item 1.8): a tool's
// description with its advice put back (method 'full') or alone ('off'), a
// check's text with or without its advice, and a tool as the turn offers it.
// The turn's use of them is in turn.test.js; what the dock, arm B and arm C
// are sent is in src/__tests__/modelVisible/.

import { describe, expect, it, vi } from 'vitest';
import { TURN_METHODS, checkText, methodCheck, toolAsSent, toolDescription } from '../copilot/method';

const tool = (fields) => ({ name: 'probe', description: 'Measure the shell (it is whole). Returns a number.', run: async () => ({}), ...fields });

describe('a tool\'s description as sent', () => {
    it('without advice: the description under either method', () => {
        for (const method of TURN_METHODS) expect(toolDescription(tool({}), method)).toBe('Measure the shell (it is whole). Returns a number.');
    });

    it('a string: after the description and a space with method full; left out with off', () => {
        const probe = tool({ method: 'Use it first.' });
        expect(toolDescription(probe, 'full')).toBe('Measure the shell (it is whole). Returns a number. Use it first.');
        expect(toolDescription(probe, 'off')).toBe(probe.description);
    });

    it('parts: anchored text goes back where it was in the description as written, strings after it in order', () => {
        const probe = tool({
            method: [
                'Then compare.',
                { after: 'Returns a number.', text: ' Read it twice.' },
                { after: 'Measure the shell (', text: 'from the g(r): ' },
                'Use it first.'
            ]
        });
        expect(toolDescription(probe, 'full'))
            .toBe('Measure the shell (from the g(r): it is whole). Returns a number. Read it twice. Then compare. Use it first.');
        expect(toolDescription(probe, 'off')).toBe(probe.description);
    });

    it('a malformed method fails under either method: an anchor not found or found twice, an empty or foreign part', () => {
        const bad = [
            [{ after: 'not there', text: ' x' }],
            [{ after: 'e', text: ' x' }],
            [{ after: 'Returns', text: '' }],
            [''],
            [{ text: 'no anchor' }],
            42
        ];
        for (const method of bad) {
            for (const mode of TURN_METHODS) expect(() => toolDescription(tool({ method }), mode)).toThrow(TypeError);
        }
    });
});

describe('a check\'s text as sent', () => {
    it('methodCheck: the text is every part, as the dock reads it; method lists the advice', () => {
        const check = methodCheck('fail', 'Refused: the window cuts the shell.', { method: ' Use the whole shell.' }, ' Ask the user to keep it.');
        expect(check).toEqual({
            level: 'fail',
            text: 'Refused: the window cuts the shell. Use the whole shell. Ask the user to keep it.',
            method: [' Use the whole shell.']
        });
        expect(checkText(check, 'full')).toBe(check.text);
        expect(checkText(check, 'off')).toBe('Refused: the window cuts the shell. Ask the user to keep it.');
        expect(() => methodCheck('note', 'A fact.', { method: '' })).toThrow(TypeError);
        expect(() => methodCheck('note', 'A fact.', 7)).toThrow(TypeError);
    });

    it('removes each piece of advice once, in order, even when the same words appear twice', () => {
        const check = methodCheck('note', 'Widen it', { method: ' (see A)' }, ' or narrow it', { method: ' (see A)' }, '.');
        expect(check.text).toBe('Widen it (see A) or narrow it (see A).');
        expect(checkText(check, 'off')).toBe('Widen it or narrow it.');
    });

    it('a check without advice is sent as written; advice not in the text, in order, is an error that does not quote it', () => {
        expect(checkText({ level: 'ok', text: 'Fine.' }, 'off')).toBe('Fine.');
        expect(checkText({ level: 'ok', text: 'Fine.', method: [] }, 'off')).toBe('Fine.');
        const broken = { level: 'note', text: 'A fact.', method: [' Run X next.'] };
        expect(checkText(broken, 'full')).toBe('A fact.');
        expect(() => checkText(broken, 'off')).toThrow(TypeError);
        try {
            checkText(broken, 'off');
        } catch (error) {
            expect(error.message).not.toContain('Run X');
        }
    });
});

describe('a tool as the turn offers it', () => {
    const result = {
        data: { value: 1 },
        checks: [methodCheck('fail', 'No value.', { method: ' Try the other tool.' }), { level: 'ok', text: 'Read.' }],
        facts: { value: 1 }
    };

    it('full: the description with its advice, no method field, the tool\'s own run()', () => {
        const run = vi.fn(async () => result);
        const probe = tool({ method: 'Use it first.', run, label: () => 'Probe', page: 'dashboard' });
        const sent = toolAsSent(probe, 'full');
        expect(sent).not.toHaveProperty('method');
        expect(sent.description).toBe('Measure the shell (it is whole). Returns a number. Use it first.');
        expect(sent.run).toBe(run);
        expect(sent.label()).toBe('Probe');
        expect(sent.page).toBe('dashboard');
        // The original is untouched.
        expect(probe.method).toBe('Use it first.');
    });

    it('off: the description alone, and run() answers with each check\'s text without its advice; data and facts as returned', async () => {
        const run = vi.fn(async () => result);
        const sent = toolAsSent(tool({ method: 'Use it first.', run }), 'off');
        expect(sent).not.toHaveProperty('method');
        expect(sent.description).toBe('Measure the shell (it is whole). Returns a number.');
        const context = { signal: new AbortController().signal };
        await expect(sent.run({ a: 1 }, context)).resolves.toEqual({
            data: { value: 1 },
            checks: [{ level: 'fail', text: 'No value.' }, { level: 'ok', text: 'Read.' }],
            facts: { value: 1 }
        });
        expect(run).toHaveBeenCalledWith({ a: 1 }, context);
    });
});
