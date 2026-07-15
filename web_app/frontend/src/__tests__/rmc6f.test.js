// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, it, expect } from 'vitest';
import { parseAtomLine } from '../rmc6f.js';
import { siteDisplacementsFromRmc6f } from '../workers/pcaKde.js';

describe('parseAtomLine', () => {
    it('parses the current 10-field format (with a bracketed type label)', () => {
        const parts = '1 Ga [1] 0.024375 0.073409 0.074344 4 0 0 0'.trim().split(/\s+/);
        expect(parseAtomLine(parts)).toEqual({
            element: 'Ga',
            coords: [0.024375, 0.073409, 0.074344],
            referenceNumber: 4,
            cellIndices: [0, 0, 0],
        });
    });

    it('parses the older 9-field format (no type label)', () => {
        const parts = '3 Pb 0.998456155893590 0.001266296130703 0.049993349135976 5 0 0 2'
            .trim().split(/\s+/);
        const atom = parseAtomLine(parts);
        expect(atom.element).toBe('Pb');
        expect(atom.referenceNumber).toBe(5);
        expect(atom.cellIndices).toEqual([0, 0, 2]);
        expect(atom.coords[0]).toBeCloseTo(0.998456155893590, 12);
        expect(atom.coords[2]).toBeCloseTo(0.049993349135976, 12);
    });

    it('rejects too-few fields and non-numeric coords / ref / cell', () => {
        expect(parseAtomLine('1 Pb 0.1 0.2 0.3 5 0 0'.split(/\s+/))).toBeNull(); // 8 fields
        expect(parseAtomLine('1 Pb x y z 5 0 0 0'.split(/\s+/))).toBeNull();      // coords NaN
        expect(parseAtomLine('1 Pb 0.1 0.2 0.3 x 0 0 0'.split(/\s+/))).toBeNull(); // ref NaN
    });
});

// End-to-end: the site extractor must load an older, label-less `.rmc6f` file
// (this is the 2018-era format shared by users) the same as the current one.
describe('siteDisplacementsFromRmc6f — older 9-field .rmc6f', () => {
    const OLD_FORMAT = `(Version 6f format configuration file)
Metadata owner:     Maxim
Metadata date:      17-07-2018
Number of atoms:                     4
Supercell dimensions:                2  2  2
Cell (Ang/deg):   10.000000  10.000000  10.000000   90.000000   90.000000   90.000000
Lattice vectors (Ang):
  10.000000   0.000000   0.000000
   0.000000  10.000000   0.000000
   0.000000   0.000000  10.000000
Atoms:
       1   Pb  0.010  0.010  0.010     1   0   0   0
       2   Pb  0.020  0.480  0.020     1   0   1   0
       3   Pb  0.490  0.010  0.480     1   1   0   1
       4   Pb  0.510  0.510  0.510     1   1   1   1
`;

    it('reads supercell, lattice, and the reference site from an old-format file', () => {
        const parsed = siteDisplacementsFromRmc6f(OLD_FORMAT);
        expect(parsed.supercell).toEqual([2, 2, 2]);
        expect(parsed.latticeVectors[0]).toEqual([10, 0, 0]);
        expect(parsed.sites).toHaveLength(1);
        expect(parsed.sites[0].referenceNumber).toBe(1);
        expect(parsed.sites[0].element).toBe('Pb');
        expect(parsed.sites[0].count).toBe(4);
        expect(parsed.sites[0].displacements).toHaveLength(4);
    });
});
