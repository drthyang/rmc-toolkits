// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The finder runs synchronously on the main thread (ModelSummary's useMemo), so it must
// refuse a basis too large to analyse interactively — a box with one reference site per
// atom (a glass, or an imported P1 configuration) — and say so, instead of freezing.

import { describe, it, expect } from 'vitest';

import { describeSymmetry, toleranceLadder, MAX_SYMMETRY_SITES } from '../symmetryModel.js';

const glass = (n) => {
    let s = 11;
    const rand = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
    const L = Math.cbrt(n / 0.07);
    return {
        latticeVectors: [[L, 0, 0], [0, L, 0], [0, 0, L]],
        supercell: [1, 1, 1],
        basis: Array.from({ length: n }, (_, i) => ({ el: i % 3 ? 'O' : 'Si', frac: [rand(), rand(), rand()] })),
    };
};

describe('basis-size cap', () => {
    it('is 2000 sites', () => {
        expect(MAX_SYMMETRY_SITES).toBe(2000);
    });

    it('skips a basis above the cap with an explanation and no ladder', () => {
        const structure = glass(MAX_SYMMETRY_SITES + 1);
        const t0 = performance.now();
        const found = describeSymmetry(structure, 0.2);
        const ladder = toleranceLadder(structure, 1.0);
        expect(performance.now() - t0).toBeLessThan(200);
        expect(found).toMatchObject({ skipped: true, spaceGroup: 'not analysed', spaceGroupNumber: null, orbits: [] });
        expect(found.pointGroup).toMatch(/2001 sites/);
        expect(found.reason).toMatch(/2000/);
        expect(ladder).toEqual([]);
    });

    it('still analyses a basis at the cap, quickly', { timeout: 20000 }, () => {
        const structure = glass(MAX_SYMMETRY_SITES);
        const t0 = performance.now();
        const found = describeSymmetry(structure, 0.2);
        expect(found.spaceGroup).toBe('P1');
        expect(toleranceLadder(structure, 1.0).map((b) => b.spaceGroup)).toEqual(['P1']);
        expect(performance.now() - t0).toBeLessThan(5000);
    });
});
