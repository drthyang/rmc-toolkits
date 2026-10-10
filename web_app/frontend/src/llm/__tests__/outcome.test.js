// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { parseOutcome } from '../copilot/outcome';

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

    it('takes the last verdict line', () => {
        const parsed = parseOutcome('Outcome: partial — draft\nmore\nOutcome: achieved — final');
        expect(parsed.status).toBe('achieved');
        expect(parsed.body).toBe('Outcome: partial — draft\nmore');
    });

    it('returns null without a verdict or for an unknown status', () => {
        expect(parseOutcome('No verdict here.')).toBeNull();
        expect(parseOutcome('Outcome: excellent — hm')).toBeNull();
        expect(parseOutcome('The outcome: achieved is in the text')).toBeNull();
        expect(parseOutcome('')).toBeNull();
        expect(parseOutcome(null)).toBeNull();
    });
});
