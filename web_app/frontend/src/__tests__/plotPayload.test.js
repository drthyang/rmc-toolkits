// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Flask plot payloads carry masked regions as JSON null (the app-wide strict
// JSON provider). The chart must treat those as gaps, and must reject a body
// that did not parse as JSON (axios then returns the raw string) with a
// message instead of throwing while rendering.

import { describe, expect, it } from 'vitest';
import { nearestFiniteIndex, niceDomain, plotPayloadError } from '../plotDomain';

describe('plotPayloadError', () => {
    it('rejects a raw string body (invalid JSON returned by axios as text)', () => {
        expect(plotPayloadError('{"series":[{"y":[NaN]}]}')).toMatch(/not valid JSON/);
    });

    it('rejects a missing or series-less payload', () => {
        expect(plotPayloadError(null)).toBeTruthy();
        expect(plotPayloadError([])).toBeTruthy();
        expect(plotPayloadError({ title: 'x' })).toMatch(/no series/);
    });

    it('accepts a payload whose series carry null gaps', () => {
        const payload = JSON.parse('{"xLabel":"Q","yLabel":"F(Q)","series":[{"label":"a","x":[1,2,3],"y":[0.1,null,0.3]}]}');
        expect(plotPayloadError(payload)).toBeNull();
    });
});

describe('null gaps in series', () => {
    it('are skipped by the axis domain', () => {
        expect(niceDomain([0.1, null, 0.3])).toEqual(niceDomain([0.1, 0.3]));
    });

    it('are never the hover target', () => {
        // x = null would coerce to 0 and win for a target near 0.
        expect(nearestFiniteIndex([null, 1, 2], [5, 6, 7], 0.1)).toBe(1);
        // A finite x with a null y is not a drawable point either.
        expect(nearestFiniteIndex([0, 1, 2], [null, 6, 7], 0)).toBe(1);
        expect(nearestFiniteIndex([0, 1], [null, NaN], 0)).toBe(-1);
    });
});
