// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The finder only tests rotations that are integer matrices in the GIVEN cell. On an
// oblique description of a lattice (a 60° or sheared cell of a tetragonal, orthorhombic or
// hexagonal crystal) some lattice rotations need entries of ±2 there, and a search over
// {-1, 0, 1} missed them: the proper subgroup that was left got its own ITA number
// (rutile on a 60° cell read "Cmmm No. 65"). And on a supercell of the crystal's own
// translation lattice, rotations of that lattice that do not map the supercell onto
// itself cannot be tested at all, so the group found is only a lower bound and must never
// carry a number (rocksalt in a 1×1×2 cell read "I4/mmm No. 139").

import { describe, it, expect } from 'vitest';

import { latticePointOps, spaceGroupAtTolerance } from '../symmetry.js';
import { centeringOfOps, pureTranslations } from '../spaceGroupSymbol.js';
import { STRUCTURES, redescribe, withNoise, cubicCell, tetragonalCell, hexagonalCell } from './fixtures/symmetryStructures.js';

const matMul = (X, Y) => X.map((r) => [0, 1, 2].map((j) => r[0] * Y[0][j] + r[1] * Y[1][j] + r[2] * Y[2][j]));

const named = (structure, tol = 0.02) => {
    const found = spaceGroupAtTolerance(structure.A, structure.basis, tol);
    return { spaceGroup: found.spaceGroup, spaceGroupNumber: found.spaceGroupNumber };
};

/** The n₁×n₂×n₃ supercell of a structure (every site repeated along each axis). */
function tiled({ A, basis }, n) {
    const out = [];
    for (const s of basis) {
        for (let i = 0; i < n[0]; i++) for (let j = 0; j < n[1]; j++) for (let k = 0; k < n[2]; k++) {
            out.push({ el: s.el, frac: [(s.frac[0] + i) / n[0], (s.frac[1] + j) / n[1], (s.frac[2] + k) / n[2]] });
        }
    }
    return { A: A.map((row, i) => row.map((v) => v * n[i])), basis: out };
}

// Unimodular basis changes (rows = new basis vectors in the old basis) whose cells are
// far from reduced.
const OBLIQUE = [
    [[1, 0, 1], [1, 1, 2], [1, 0, 2]],
    [[1, 1, -3], [0, 1, -4], [0, 0, 1]],
    [[0, -1, 1], [0, 1, 0], [-1, -1, 1]],
    [[1, 2, 0], [0, 1, -1], [-3, -6, 1]],
];
const SIXTY = [[1, 0, 0], [1, 1, 0], [0, 0, 1]];
const SHEARED = [[1, 0, 0], [1, 1, 0], [1, 0, 1]];

describe('lattice rotations on any basis of the lattice', () => {
    it.each([
        ['primitive cubic', cubicCell(4), 48],
        ['F-cubic (primitive cell)', [[0, 2, 2], [2, 0, 2], [2, 2, 0]], 48],
        ['I-cubic (primitive cell)', [[-2, 2, 2], [2, -2, 2], [2, 2, -2]], 48],
        ['hexagonal', hexagonalCell(4, 6.5), 24],
        ['tetragonal', tetragonalCell(4, 6), 16],
        ['orthorhombic', [[4, 0, 0], [0, 5, 0], [0, 0, 6]], 8],
    ])('finds every rotation of a %s lattice on oblique cells', (_, A, order) => {
        expect(latticePointOps(A, 1e-6)).toHaveLength(order);
        for (const U of OBLIQUE) expect(latticePointOps(matMul(U, A), 1e-6), JSON.stringify(U)).toHaveLength(order);
    });
});

describe('groups on oblique cells are named in full', () => {
    it.each([
        ['rutile on a 60° cell', () => redescribe(STRUCTURES.rutile(), SIXTY), 'P4_2/mnm', 136],
        ['rutile on a sheared cell', () => redescribe(STRUCTURES.rutile(), SHEARED), 'P4_2/mnm', 136],
        ['Pnma perovskite on a 60° cell', () => redescribe(STRUCTURES.pnmaPerovskite(), SIXTY), 'Pnma', 62],
        ['Pnma perovskite on a sheared cell', () => redescribe(STRUCTURES.pnmaPerovskite(), SHEARED), 'Pnma', 62],
        ['hcp on a sheared cell', () => redescribe(STRUCTURES.hcp(), SHEARED), 'P6_3/mmc', 194],
        ['hcp on a γ = 60° cell', () => redescribe(STRUCTURES.hcp(), SIXTY), 'P6_3/mmc', 194],
        ['wurtzite on a γ = 60° cell', () => redescribe(STRUCTURES.wurtzite(), SIXTY), 'P6_3mc', 186],
        ['bismuth on a sheared cell', () => redescribe(STRUCTURES.bismuth(), SHEARED), 'R-3m', 166],
    ])('%s', (_, make, symbol, number) => {
        expect(named(make())).toEqual({ spaceGroup: symbol, spaceGroupNumber: number });
        expect(named(withNoise(make(), 0.01, 3), 0.1)).toEqual({ spaceGroup: symbol, spaceGroupNumber: number });
    });
});

describe('groups on strongly oblique cells are named in their conventional cell', () => {
    // The naming step builds the conventional cell from the shortest lattice vectors along
    // the symmetry axes. Searched in a ±2 box of the given basis, those vectors are out of
    // reach on a strongly oblique cell (coefficients up to 6 here), and the group — found
    // in full — was left as its crystal class.
    it.each([
        ['CsCl', 'cscl', [[1, 1, -3], [0, 1, -4], [0, 0, 1]]],
        ['rocksalt', 'rocksalt', [[1, 0, 0], [3, 1, 0], [2, 2, 1]]],
        ['diamond', 'diamond', [[1, 2, 0], [0, 1, -1], [-3, -6, 1]]],
        ['zincblende', 'zincblende', [[1, 0, 0], [3, 1, 0], [2, 2, 1]]],
        ['bcc iron', 'bccIron', [[1, 1, -3], [0, 1, -4], [0, 0, 1]]],
        ['perovskite', 'perovskite', [[1, 2, 0], [0, 1, -1], [-3, -6, 1]]],
        ['hcp', 'hcp', [[1, 2, 0], [0, 1, -1], [-3, -6, 1]]],
    ])('%s', (_, key, U) => {
        const base = STRUCTURES[key]();
        const cell = redescribe(base, U);
        expect(cell.basis).toHaveLength(base.basis.length);
        expect(named(cell)).toEqual({ spaceGroup: base.symbol, spaceGroupNumber: base.number });
        expect(named(withNoise(cell, 0.01, 3), 0.1)).toEqual({ spaceGroup: base.symbol, spaceGroupNumber: base.number });
    });
});

describe('a supercell of the translation lattice gives a lower bound, never a number', () => {
    it('does not name rocksalt in a 1×1×2 cell I4/mmm', () => {
        // The F lattice of rocksalt has cubic 3-folds that do not map the 5.6 × 5.6 × 11.3 Å
        // cell onto itself, so they are never tried; the verified group is tetragonal and
        // only a lower bound on Fm-3m.
        const cell = redescribe(STRUCTURES.rocksalt(), [[1, 0, 0], [0, 1, 0], [0, 0, 2]]);
        expect(named(cell)).toEqual({ spaceGroup: '≥ I4/mmm', spaceGroupNumber: null });
        expect(named(withNoise(cell, 0.01, 3), 0.1)).toEqual({ spaceGroup: '≥ I4/mmm', spaceGroupNumber: null });
    });

    it('does not name a Pnma crystal on a triple oblique supercell P-1', () => {
        // Only ±1 map this cell's lattice onto itself; the 2-folds and mirrors of the
        // orthorhombic translation lattice cannot be tested there.
        const cell = redescribe(STRUCTURES.pnmaPerovskite(), [[1, 1, 0], [0, 1, 1], [1, 0, 2]]);
        expect(named(cell)).toEqual({ spaceGroup: '≥ P-1', spaceGroupNumber: null });
    });

    it('marks an unnamed class as a lower bound too', () => {
        // Translations of 1/26 are beyond the snapping grid (orders ≤ 25): the group cannot
        // be named, and the cubic 3-folds were not testable either.
        expect(named(tiled(STRUCTURES.cscl(), [26, 1, 1]))).toEqual({ spaceGroup: '≥ 4/mmm class', spaceGroupNumber: null });
    });
});

describe('supercell translations are snapped to their exact fractions', () => {
    it('names CsCl in a 5×5×5 cell Pm-3m', () => {
        // Fifths are not on a 1/24 grid (0.2 is 0.008 from 5/24): snapped there, the
        // translations did not form a lattice and the group went unnamed.
        expect(named(tiled(STRUCTURES.cscl(), [5, 5, 5]))).toEqual({ spaceGroup: 'Pm-3m', spaceGroupNumber: 221 });
    });

    it('names CsCl in a 7×1×1 cell as a lower bound, with sevenths snapped', () => {
        expect(named(tiled(STRUCTURES.cscl(), [7, 1, 1]))).toEqual({ spaceGroup: '≥ P4/mmm', spaceGroupNumber: null });
    });

    it('snaps noisy centring translations whole, from the order of the translation group', () => {
        // Refined translations of noisy fixtures at a loose tolerance (R-3m and I-42m of the
        // 230-group set, 0.08 Å noise, τ = 0.6 Å): up to 0.0075 off their fractions. The group
        // order fixes the grid (thirds for three translations, halves for two), so they snap.
        const I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        const ops = (ts) => ts.map((t) => ({ R: I, t }));
        const R = ops([[0, 0, 0], [0.6624, 0.3287, 0.3333], [0.3346, 0.6682, 0.6659]]);
        expect(pureTranslations(R)).toEqual([[0, 0, 0], [1 / 3, 2 / 3, 2 / 3], [2 / 3, 1 / 3, 1 / 3]]);
        expect(centeringOfOps(R)).toBe('R');
        expect(centeringOfOps(ops([[0, 0, 0], [0.5025, 0.5049, 0.5012]]))).toBe('I');
        // A translation that is on no grid of the group order does not snap.
        expect(pureTranslations(ops([[0, 0, 0], [0.47, 0.5, 0.5]]))).toBeNull();
        // Fifths of a 5-fold supercell, with noise.
        expect(pureTranslations(ops([0, 1, 2, 3, 4].map((k) => [k / 5 + 0.004, 0, 0]))).map((t) => t[0]))
            .toEqual([0, 0.2, 0.4, 0.6, 0.8]);
    });

    it('gives CsCl in a 5×1×1 cell the tetragonal lower bound', () => {
        expect(named(tiled(STRUCTURES.cscl(), [5, 1, 1]))).toEqual({ spaceGroup: '≥ P4/mmm', spaceGroupNumber: null });
    });
});
