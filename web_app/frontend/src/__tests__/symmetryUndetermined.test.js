// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// When nothing can be analysed — a lattice with a non-finite entry, a singular lattice, or
// no basis — no operation survives, not even the identity. The finder then reported
// "P1 No. 1" with 0 operations: a space-group number for a structure it never saw. The
// honest answer is "undetermined", with no number and no orbits.

import { describe, it, expect } from 'vitest';

import { findSpaceGroupOps, spaceGroupAtTolerance, symmetryLadder } from '../symmetry.js';
import { describeSymmetry, toleranceLadder } from '../symmetryModel.js';

const basis = [{ el: 'Cs', frac: [0, 0, 0] }, { el: 'Cl', frac: [0.5, 0.5, 0.5] }];

describe('a structure that cannot be analysed is undetermined, never P1', () => {
    it.each([
        ['a NaN lattice entry', [[Number.NaN, 0, 0], [0, 4, 0], [0, 0, 4]]],
        ['an infinite lattice entry', [[Number.POSITIVE_INFINITY, 0, 0], [0, 4, 0], [0, 0, 4]]],
        ['a singular lattice', [[4, 0, 0], [0, 4, 0], [0, 0, 0]]],
    ])('%s', (_, A) => {
        expect(spaceGroupAtTolerance(A, basis, 0.2)).toMatchObject({ spaceGroup: 'undetermined', spaceGroupNumber: null, nSpace: 0 });
        expect(findSpaceGroupOps(A, basis, 0.2)).toMatchObject({ spaceGroup: 'undetermined', spaceGroupNumber: null, nSpace: 0 });
        expect(symmetryLadder(A, basis, 1.0)).toEqual([]);
        const structure = { latticeVectors: A, supercell: [1, 1, 1], basis };
        expect(describeSymmetry(structure, 0.2)).toMatchObject({ spaceGroup: 'undetermined', spaceGroupNumber: null, nSpace: 0, orbits: [] });
        expect(toleranceLadder(structure, 1.0)).toEqual([]);
    });

    it('an empty basis', () => {
        const A = [[4, 0, 0], [0, 4, 0], [0, 0, 4]];
        expect(spaceGroupAtTolerance(A, [], 0.2)).toMatchObject({ spaceGroup: 'undetermined', spaceGroupNumber: null, nSpace: 0 });
        expect(findSpaceGroupOps(A, [], 0.2)).toMatchObject({ spaceGroup: 'undetermined', spaceGroupNumber: null, nSpace: 0 });
    });

    it('still names a valid structure', () => {
        expect(spaceGroupAtTolerance([[4, 0, 0], [0, 4, 0], [0, 0, 4]], basis, 0.2)).toMatchObject({ spaceGroup: 'Pm-3m', spaceGroupNumber: 221 });
    });
});
