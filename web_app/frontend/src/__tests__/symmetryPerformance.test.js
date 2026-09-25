// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The finder runs synchronously on the main thread (ModelSummary's useMemo), so the page
// freezes for as long as it takes. The basis-size cap (symmetryLimits.test.js) does not
// bound that time on its own: a crystalline box declared as a 1×1×1 supercell has one
// pure translation per repeat unit, and every rotation of the lattice then holds with
// every one of them, so the candidate operations grow with the square of the box: a
// 4×4×4 rocksalt box has 12 288, and the ladder of a noisy 3×3×3 one (5184) took over a
// minute.
// Such a cell is a supercell of the structure's repeat unit — the .rmc6f "Supercell
// dimensions" do not describe it — and the card must say so instead of freezing.

import { describe, it, expect } from 'vitest';

import { describeSymmetry, toleranceLadder, MAX_SYMMETRY_OPS } from '../symmetryModel.js';
import { operationEstimate } from '../symmetry.js';
import { STRUCTURES, redescribe, withNoise } from './fixtures/symmetryStructures.js';

const box = ({ A, basis }) => ({ latticeVectors: A, supercell: [1, 1, 1], basis });
const timed = (f) => { const t0 = performance.now(); const out = f(); return [out, performance.now() - t0]; };

describe('the operation budget', () => {
    it('counts the pure translations and lattice rotations a pass would combine', () => {
        const rocksalt = STRUCTURES.rocksalt();
        expect(operationEstimate(rocksalt.A, rocksalt.basis, 0.2)).toEqual({ translations: 4, rotations: 48, operations: 192, exceeds: false });
        const big = redescribe(rocksalt, [[4, 0, 0], [0, 4, 0], [0, 0, 4]]);
        expect(operationEstimate(big.A, big.basis, 1.0)).toEqual({ translations: 256, rotations: 48, operations: 12288, exceeds: false });
        // With a limit it stops counting as soon as the product exceeds it.
        expect(operationEstimate(big.A, big.basis, 1.0, MAX_SYMMETRY_OPS)).toEqual({ translations: 9, rotations: 48, operations: 432, exceeds: true });
    });

    it('allows every correctly declared cell (at most 48 rotations × 4 centring translations)', () => {
        expect(MAX_SYMMETRY_OPS).toBeGreaterThanOrEqual(192);
    });

    it('declines a crystalline box that is a supercell of its repeat unit, quickly and with a reason', () => {
        const structure = box(redescribe(STRUCTURES.rocksalt(), [[4, 0, 0], [0, 4, 0], [0, 0, 4]]));
        const [found, t1] = timed(() => describeSymmetry(structure, 0.2));
        const [ladder, t2] = timed(() => toleranceLadder(structure, 1.0));
        expect(found).toMatchObject({ skipped: true, spaceGroup: 'not analysed', spaceGroupNumber: null, orbits: [] });
        expect(Number.isNaN(found.maxResidual)).toBe(true);
        expect(found.pointGroup).toMatch(/translations per cell/);
        expect(found.reason).toMatch(/supercell/);
        expect(found.reason).toMatch(String(MAX_SYMMETRY_OPS));
        expect(ladder).toEqual([]);
        expect(t1 + t2).toBeLessThan(1500);
    });

    it('declines the noisy box too, headline and ladder alike', () => {
        const structure = box(withNoise(redescribe(STRUCTURES.rocksalt(), [[3, 0, 0], [0, 3, 0], [0, 0, 3]]), 0.05, 1));
        const [found, t1] = timed(() => describeSymmetry(structure, 0.2));
        const [ladder, t2] = timed(() => toleranceLadder(structure, 1.0));
        expect(found.skipped).toBe(true);
        expect(ladder).toEqual([]);
        expect(t1 + t2).toBeLessThan(1500);
    });

    it('still analyses a doubled cell within the budget', () => {
        const structure = box(withNoise(redescribe(STRUCTURES.perovskite(), [[2, 0, 0], [0, 2, 0], [0, 0, 2]]), 0.03, 1));
        const found = describeSymmetry(structure, 0.2);
        expect(found.skipped).toBeUndefined();
        expect(found).toMatchObject({ spaceGroup: 'Pm-3m', spaceGroupNumber: 221 });
        expect(toleranceLadder(structure, 1.0).at(-1)).toMatchObject({ spaceGroup: 'Pm-3m' });
    });
});
