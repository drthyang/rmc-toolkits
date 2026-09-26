// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A Wyckoff label is a multiplicity AND a letter of the same standard cell. When the group
// is named in a cell other than the one the model uses (a primitive or rhombohedral cell,
// an R subgroup kept in its F-cubic parent cell), the letter is read in the naming cell,
// so the multiplicity printed beside it must be that cell's too — '3a' for the Ga of an R3m
// lacunar spinel, never the given F-cubic cell's '4a', which R3m does not have.

import { describe, it, expect } from 'vitest';

import { describeSymmetry, orbitLabel } from '../symmetryModel.js';
import { WYCKOFF_DATA } from '../wyckoffTable.js';
import { STRUCTURES, AXIS_SETTINGS, redescribe, lacunarSpinel } from './fixtures/symmetryStructures.js';

const describe0 = ({ A, basis }, tol = 0.02) => describeSymmetry({ latticeVectors: A, supercell: [1, 1, 1], basis }, tol);
const tabulated = (number) => new Set(WYCKOFF_DATA[number].split(';').map((row) => {
    const [letter, multiplicity] = row.split(':');
    return `${multiplicity}${letter}`;
}));
const labels = (found) => Object.fromEntries(found.orbits.map((o) => [o.element, orbitLabel(o)]));

// Every orbit that has a letter is labelled with a (multiplicity, letter) pair the group's
// table actually lists.
const expectTabulated = (found) => {
    const table = tabulated(found.spaceGroupNumber);
    for (const o of found.orbits) {
        if (!o.wyckoff) continue;
        expect(table.has(orbitLabel(o)), `${found.spaceGroup} ${o.element}: ${orbitLabel(o)}`).toBe(true);
    }
};

describe('Wyckoff labels carry the multiplicity of the cell the letter is read in', () => {
    it('labels the Ga of an R3m lacunar spinel in its F-cubic cell 3a', () => {
        const found = describe0(lacunarSpinel(0.01), 0.005);
        expect(found.spaceGroup).toBe('R3m');
        const ga = found.orbits.find((o) => o.element === 'Ga');
        expect(ga.size).toBe(4);                 // 4 Ga in the given F-cubic cell ...
        expect(ga.wyckoffMultiplicity).toBe(3);  // ... 3 in the hexagonal cell R3m is named in
        expect(orbitLabel(ga)).toBe('3a');
        expectTabulated(found);
    });

    it('labels rocksalt on its primitive cell 4a / 4b, as in the Fm-3m cell', () => {
        const fccPrim = [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]];
        const found = describe0(redescribe(STRUCTURES.rocksalt(), fccPrim));
        expect(found.spaceGroup).toBe('Fm-3m');
        expect(found.orbits.map((o) => o.size)).toEqual([1, 1]);
        const l = labels(found);
        expect(new Set([l.Na, l.Cl])).toEqual(new Set(['4a', '4b']));
        expectTabulated(found);
    });

    it('keeps the given multiplicity when the group is named in the given cell', () => {
        const found = describe0(STRUCTURES.rutile());
        expect(labels(found)).toEqual({ Ti: '2a', O: '4f' });
        for (const o of found.orbits) expect(o.wyckoffMultiplicity).toBe(o.size);
        expect(labels(describe0(redescribe(STRUCTURES.rutile(), AXIS_SETTINGS.cab)))).toEqual({ Ti: '2a', O: '4f' });
    });

    it('labels an orbit with no letter by its given multiplicity and site symmetry, and gives it no naming-cell multiplicity', () => {
        const found = describe0(redescribe(STRUCTURES.perovskite(), [[2, 0, 0], [0, 2, 0], [0, 0, 1]]));
        expect(found.spaceGroupNumber).toBeNull();
        for (const o of found.orbits) {
            expect(o.wyckoffMultiplicity).toBeNull();
            expect(orbitLabel(o)).toBe(`${o.size} (${o.site})`);
        }
    });
});
