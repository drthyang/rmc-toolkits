// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Detected SG card names a group in its ITA standard setting whatever cell the RMC
// model happens to use: another axis order, a centred cell choice, the principal axis
// along a or b, a primitive cell, or a supercell of the true cell. Where it cannot, it
// reports the crystal class with no number — never another group's symbol or number.

import { describe, it, expect } from 'vitest';

import { spaceGroupAtTolerance } from '../symmetry.js';
import {
    STRUCTURES, AXIS_SETTINGS, redescribe, withNoise, lacunarSpinel, at, cubicCell, tetragonalCell, centred,
} from './fixtures/symmetryStructures.js';
import { SPACE_GROUP_FIXTURES, structureFor } from './fixtures/spaceGroups.js';

const named = (structure, tol = 0.02) => {
    const found = spaceGroupAtTolerance(structure.A, structure.basis, tol);
    return { spaceGroup: found.spaceGroup, spaceGroupNumber: found.spaceGroupNumber };
};

describe('every fixture group in every axis order', () => {
    it.each(SPACE_GROUP_FIXTURES)('$symbol (#$number) is named correctly or not at all', (fixture) => {
        const base = structureFor(fixture);
        for (const [name, M] of Object.entries(AXIS_SETTINGS)) {
            const found = named(redescribe(base, M));
            if (found.spaceGroupNumber !== null) expect(found, name).toEqual({ spaceGroup: fixture.symbol, spaceGroupNumber: fixture.number });
            // Every group with its principal axis moved is still recognised.
            expect(found.spaceGroupNumber, `${fixture.symbol} in ${name}`).toBe(fixture.number);
        }
    });
});

describe('centred and non-conventional cells', () => {
    it('names Amm2 BaTiO3 on its pseudo-cubic cell', () => {
        const d = 0.02;
        const basis = [at('Ba', 0, 0, 0), at('Ti', 0.5 + d, 0.5 + d, 0.5), at('O', 0.5, 0.5, 0), at('O', 0.5, 0, 0.5), at('O', 0, 0.5, 0.5)];
        expect(named({ A: cubicCell(4.0), basis }, 0.01)).toEqual({ spaceGroup: 'Amm2', spaceGroupNumber: 38 });
    });

    it('names Cmm2 given on primitive tetragonal axes', () => {
        const basis = [at('X', 0, 0, 0), at('Y', 0.2, 0.2, 0.3), at('Y', -0.2, -0.2, 0.3)];
        expect(named({ A: tetragonalCell(4.1, 5.3), basis }, 0.01)).toEqual({ spaceGroup: 'Cmm2', spaceGroupNumber: 35 });
    });

    it('names Amm2 in every axis order (the centering letter moves with the axes)', () => {
        const amm2 = structureFor(SPACE_GROUP_FIXTURES.find((f) => f.number === 38));
        for (const M of Object.values(AXIS_SETTINGS)) expect(named(redescribe(amm2, M))).toEqual({ spaceGroup: 'Amm2', spaceGroupNumber: 38 });
    });

    it('names rocksalt and bcc iron on their primitive cells', () => {
        const fccPrim = [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]];
        const bccPrim = [[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]];
        expect(named(redescribe(STRUCTURES.rocksalt(), fccPrim))).toEqual({ spaceGroup: 'Fm-3m', spaceGroupNumber: 225 });
        expect(named(redescribe(STRUCTURES.bccIron(), bccPrim))).toEqual({ spaceGroup: 'Im-3m', spaceGroupNumber: 229 });
        expect(named(redescribe(STRUCTURES.diamond(), fccPrim))).toEqual({ spaceGroup: 'Fd-3m', spaceGroupNumber: 227 });
    });

    it('names C2/c given in its I2/a cell choice and P2_1/c given as P2_1/n', () => {
        const c2c = structureFor(SPACE_GROUP_FIXTURES.find((f) => f.number === 15));
        expect(named(redescribe(c2c, [[1, 0, 1], [0, 1, 0], [0, 0, 1]]))).toEqual({ spaceGroup: 'C2/c', spaceGroupNumber: 15 });
        const p21c = structureFor(SPACE_GROUP_FIXTURES.find((f) => f.number === 14));
        expect(named(redescribe(p21c, [[1, 0, 0], [0, 1, 0], [1, 0, 1]]))).toEqual({ spaceGroup: 'P2_1/c', spaceGroupNumber: 14 });
    });

    it('names a mirror on a cell diagonal Cm, not Pm', () => {
        // A mirror ⊥ [1-10] of a primitive cubic lattice: the monoclinic cell with b along
        // [1-10] is C-centred. Accepting the given cell with its unique axis on the diagonal
        // spelled "Pm" (No. 6).
        const basis = [at('X', 0, 0, 0), at('Y', 0.1, 0.1, 0.3)];
        expect(named({ A: cubicCell(4.0), basis }, 0.01)).toEqual({ spaceGroup: 'Cm', spaceGroupNumber: 8 });
    });

    it('names R-3m bismuth given on rhombohedral axes', () => {
        const rhombohedral = [[2 / 3, 1 / 3, 1 / 3], [-1 / 3, 1 / 3, 1 / 3], [-1 / 3, -2 / 3, 1 / 3]];
        expect(named(redescribe(STRUCTURES.bismuth(), rhombohedral))).toEqual({ spaceGroup: 'R-3m', spaceGroupNumber: 166 });
    });

    it('tells I2_12_12_1 from I222 and I2_13 from I23 at any origin', () => {
        // Same element types, different arrangement: only I222 and I23 have a point where
        // the three 2-folds meet. An origin shift must not change the answer.
        for (const [number, symbol] of [[23, 'I222'], [24, 'I2_12_12_1'], [197, 'I23'], [199, 'I2_13']]) {
            const base = structureFor(SPACE_GROUP_FIXTURES.find((f) => f.number === number));
            for (const shift of [[0, 0, 0], [0.13, 0.29, 0.41]]) {
                const moved = redescribe(base, AXIS_SETTINGS.abc, shift);
                expect(named(moved), `${symbol} shifted ${shift}`).toEqual({ spaceGroup: symbol, spaceGroupNumber: number });
            }
        }
    });
});

describe('principal axis not along c', () => {
    it('names rutile with its 4_2 axis along a or b', () => {
        expect(named(redescribe(STRUCTURES.rutile(), AXIS_SETTINGS.cab))).toEqual({ spaceGroup: 'P4_2/mnm', spaceGroupNumber: 136 });
        expect(named(redescribe(STRUCTURES.rutile(), AXIS_SETTINGS.bca))).toEqual({ spaceGroup: 'P4_2/mnm', spaceGroupNumber: 136 });
    });

    it('names hcp and wurtzite with c along a, and wurtzite on a 60° cell', () => {
        expect(named(redescribe(STRUCTURES.hcp(), AXIS_SETTINGS.cab))).toEqual({ spaceGroup: 'P6_3/mmc', spaceGroupNumber: 194 });
        expect(named(redescribe(STRUCTURES.wurtzite(), AXIS_SETTINGS.bca))).toEqual({ spaceGroup: 'P6_3mc', spaceGroupNumber: 186 });
        expect(named(redescribe(STRUCTURES.wurtzite(), [[1, 0, 0], [1, 1, 0], [0, 0, 1]]))).toEqual({ spaceGroup: 'P6_3mc', spaceGroupNumber: 186 });
    });

    it('names the Pnma perovskite in every axis order', () => {
        for (const M of Object.values(AXIS_SETTINGS)) expect(named(redescribe(STRUCTURES.pnmaPerovskite(), M))).toEqual({ spaceGroup: 'Pnma', spaceGroupNumber: 62 });
    });
});

describe('subgroups of a centred cubic parent keep its cell', () => {
    it('names the R3m lacunar spinel in the F-43m cell', () => {
        // Ga moved along [111]: the low-temperature distortion of GaV4S8-type spinels.
        expect(named(lacunarSpinel(0.01), 0.005)).toEqual({ spaceGroup: 'R3m', spaceGroupNumber: 160 });
    });

    it('names a tetragonally strained rocksalt in its F cell as I4/mmm', () => {
        const { basis } = STRUCTURES.rocksalt();
        expect(named({ A: tetragonalCell(5.64, 5.64 * 1.05), basis }, 0.01)).toEqual({ spaceGroup: 'I4/mmm', spaceGroupNumber: 139 });
    });

    it('names a noisy R3m lacunar spinel at a tolerance above the noise', () => {
        const noisy = withNoise(lacunarSpinel(0.01), 0.01, 5);
        expect(named(noisy, 0.06)).toEqual({ spaceGroup: 'R3m', spaceGroupNumber: 160 });
    });
});

describe('supercells of the true cell', () => {
    // The declared cell is a multiple of the true one: the pure translations then form a
    // finer lattice than any Bravais centering, and must not be read as F, I or C.
    it('names CsCl and perovskite in doubled cells as Pm-3m', () => {
        const double = [[2, 0, 0], [0, 2, 0], [0, 0, 2]];
        expect(named(redescribe(STRUCTURES.cscl(), double))).toEqual({ spaceGroup: 'Pm-3m', spaceGroupNumber: 221 });
        expect(named(redescribe(STRUCTURES.perovskite(), double))).toEqual({ spaceGroup: 'Pm-3m', spaceGroupNumber: 221 });
    });

    it('does not read a 2×2×1 cell of a primitive structure as C-centred', () => {
        // The pure translations (½,0,0), (0,½,0), (½,½,0) contain the C vector but are not
        // a C lattice. The tetragonal cell cannot test the cubic 3-folds (they do not map a
        // 7.8 × 7.8 × 3.9 Å lattice onto itself), so the verified P4/mmm is a lower bound.
        const found = named(redescribe(STRUCTURES.perovskite(), [[2, 0, 0], [0, 2, 0], [0, 0, 1]]));
        expect(found).toEqual({ spaceGroup: '≥ P4/mmm', spaceGroupNumber: null });
    });

    it('marks a group found in a rotated supercell as a lower bound', () => {
        // In a √2×√2×2 cell the cubic 3-folds do not map the cell's lattice onto itself,
        // so they are never tried: the verified group, P4/mmm in the 3.9 Å cell, is only
        // a lower bound on the true Pm-3m, and is shown as one, with no number.
        const found = named(redescribe(STRUCTURES.perovskite(), [[1, 1, 0], [-1, 1, 0], [0, 0, 2]]));
        expect(found).toEqual({ spaceGroup: '≥ P4/mmm', spaceGroupNumber: null });
    });

    it('names rocksalt in a doubled F cell as Fm-3m', () => {
        expect(named(redescribe(STRUCTURES.rocksalt(), [[2, 0, 0], [0, 2, 0], [0, 0, 2]]))).toEqual({ spaceGroup: 'Fm-3m', spaceGroupNumber: 225 });
    });

    it('keeps a genuine superstructure: a doubled cell with one site changed is not the small cell', () => {
        const doubled = redescribe(STRUCTURES.cscl(), [[2, 0, 0], [0, 2, 0], [0, 0, 2]]);
        const basis = doubled.basis.map((s, i) => (i === 0 ? { ...s, el: 'Rb' } : s));
        const found = spaceGroupAtTolerance(doubled.A, basis, 0.02);
        expect(found.spaceGroup).toBe('Pm-3m');   // Rb on one corner of the 2×2×2 cell: still cubic
        expect(found.nSpace).toBe(48);            // but with no extra translations
    });
});

describe('noise does not change the name', () => {
    it.each([
        ['rutile, 4_2 along a', () => redescribe(STRUCTURES.rutile(), AXIS_SETTINGS.cab), 'P4_2/mnm', 136],
        ['Pnma perovskite, Pbnm axes', () => redescribe(STRUCTURES.pnmaPerovskite(), AXIS_SETTINGS.bca), 'Pnma', 62],
        ['zincblende', () => STRUCTURES.zincblende(), 'F-43m', 216],
        ['diamond', () => STRUCTURES.diamond(), 'Fd-3m', 227],
        ['wurtzite', () => STRUCTURES.wurtzite(), 'P6_3mc', 186],
        ['bismuth', () => STRUCTURES.bismuth(), 'R-3m', 166],
    ])('%s', (_, make, symbol, number) => {
        expect(named(withNoise(make(), 0.01, 9), 0.1)).toEqual({ spaceGroup: symbol, spaceGroupNumber: number });
    });
});

// silence unused-helper lint for helpers kept for future cases
void centred;
