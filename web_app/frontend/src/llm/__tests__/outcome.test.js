// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { parseOutcome } from '../copilot/verification/outcome';

describe('parseOutcome', () => {
    it('splits a plain verdict line from the body', () => {
        expect(parseOutcome('Answer text.\n\nOutcome: achieved — Rwp quoted for all 3 datasets.')).toEqual({
            status: 'achieved',
            note: 'Rwp quoted for all 3 datasets.',
            body: 'Answer text.'
        });
    });

    it('tolerates Markdown decoration, other dashes and casing', () => {
        expect(parseOutcome('x\n**Outcome:** partial – the Ta site lacks a PCA table.')).toMatchObject({
            status: 'partial', note: 'the Ta site lacks a PCA table.'
        });
        expect(parseOutcome('x\n- *Outcome: NOT ACHIEVED* - no partials file.')).toMatchObject({
            status: 'not achieved', note: 'no partials file.'
        });
        expect(parseOutcome('x\nOutcome: not_achieved: missing data')).toMatchObject({
            status: 'not achieved', note: 'missing data'
        });
        expect(parseOutcome('x\n**Outcome: achieved**')).toMatchObject({ status: 'achieved', note: '' });
    });

    it('accepts the line in backticks, as the prompt shows it, and as a heading', () => {
        expect(parseOutcome('Answer.\n\n`Outcome: achieved — the slope is −0.0005 per row.`')).toEqual({
            status: 'achieved',
            note: 'the slope is −0.0005 per row.',
            body: 'Answer.'
        });
        expect(parseOutcome('x\n`Outcome: partial` — no Ga site table.')).toMatchObject({
            status: 'partial', note: 'no Ga site table.'
        });
        expect(parseOutcome('x\n### Outcome: not achieved — no partials file.')).toMatchObject({
            status: 'not achieved', note: 'no partials file.'
        });
    });

    it('takes the last verdict line', () => {
        const parsed = parseOutcome('Outcome: partial — draft\nmore\nOutcome: achieved — final');
        expect(parsed.status).toBe('achieved');
        expect(parsed.body).toBe('Outcome: partial — draft\nmore');
    });

    it('reads a line with no status word as the explanation (the verdict is the app\'s: THY-18)', () => {
        expect(parseOutcome('Answer.\n\nOutcome: the split site is ordered at (0, 0, ½), from analyze_modulations.')).toEqual({
            status: null, note: 'the split site is ordered at (0, 0, ½), from analyze_modulations.', body: 'Answer.'
        });
        expect(parseOutcome('x\n**Outcome:** — no partials file, so no shell.')).toMatchObject({ status: null, note: 'no partials file, so no shell.' });
        // A word that is no status stays in the sentence.
        expect(parseOutcome('Outcome: excellent — hm')).toMatchObject({ status: null, note: 'excellent — hm' });
    });

    it('returns null without an Outcome line', () => {
        expect(parseOutcome('No verdict here.')).toBeNull();
        expect(parseOutcome('The outcome: achieved is in the text')).toBeNull();
        expect(parseOutcome('')).toBeNull();
        expect(parseOutcome(null)).toBeNull();
    });
});
