// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Wyckoff letters are read in the standard cell the group is NAMED in, not in the cell the
// model happens to use: the coordinates are transformed to that cell first, and letters
// are withheld when there is no such cell.

import { describe, it, expect } from 'vitest';

import { describeSymmetry } from '../symmetryModel.js';
import { WYCKOFF_DATA } from '../wyckoffTable.js';
import { STRUCTURES, AXIS_SETTINGS, redescribe, lacunarSpinel, orbits } from './fixtures/symmetryStructures.js';

const describe0 = ({ A, basis }, tol = 0.02) => describeSymmetry({ latticeVectors: A, supercell: [1, 1, 1], basis }, tol);
const letters = (symmetry) => Object.fromEntries(symmetry.orbits.map((o) => [o.element, `${o.size}${o.wyckoff ?? '?'}`]));

// P2_1/c (#14) with B on 2a (0,0,0), X on 2b (1/2,0,0) and G on the general position 4e
// (without G the special positions alone have the higher symmetry C2/m).
const p21c = () => {
    const ops = [(x, y, z) => [x, y, z], (x, y, z) => [-x, y + 0.5, -z + 0.5], (x, y, z) => [-x, -y, -z], (x, y, z) => [x, -y + 0.5, z + 0.5]];
    const A = [[7.1, 0, 0], [0, 5.3, 0], [9.7 * Math.cos((104 * Math.PI) / 180), 0, 9.7 * Math.sin((104 * Math.PI) / 180)]];
    return { A, basis: orbits(ops, [['B', 0, 0, 0], ['X', 0.5, 0, 0], ['G', 0.137, 0.213, 0.061]]) };
};

describe('letters in the naming cell', () => {
    it('labels P2_1/c in its standard setting', () => {
        const found = describe0(p21c());
        expect(found.spaceGroup).toBe('P2_1/c');
        expect(letters(found)).toEqual({ B: '2a', X: '2b', G: '4e' });
    });

    it('gives a consistent description when the axes are relabelled', () => {
        // The B–X separation is a/2, perpendicular to the glide, so X can be 2b or 2d
        // relative to B on 2a/2c — never 2c, which is what the untransformed coordinates gave.
        const valid = new Set(['2a|2b', '2c|2d', '2b|2a', '2d|2c']);
        for (const [name, M] of Object.entries(AXIS_SETTINGS)) {
            const found = describe0(redescribe(p21c(), M));
            expect(found.spaceGroup, name).toBe('P2_1/c');
            const l = letters(found);
            expect(valid.has(`${l.B}|${l.X}`), `${name}: ${l.B} ${l.X}`).toBe(true);
            expect(l.G).toBe('4e');
        }
    });

    it('labels the R3m lacunar spinel on hexagonal axes with hexagonal multiplicities', () => {
        const found = describe0(lacunarSpinel(0.01), 0.005);
        expect(found.spaceGroup).toBe('R3m');
        // Ga sits on the 3-fold axis: 3a in the hexagonal cell (4 Ga in the F-cubic cell).
        const ga = found.orbits.find((o) => o.element === 'Ga');
        expect(ga.size).toBe(4);
        expect(ga.wyckoff).toBe('a');
        // Every orbit gets a letter whose multiplicity matches the hexagonal cell.
        for (const o of found.orbits) expect(o.wyckoff, `${o.element} ${o.size}`).not.toBeNull();
    });

    it('labels rutile the same whichever axis carries the 4_2', () => {
        const reference = letters(describe0(STRUCTURES.rutile()));
        expect(reference).toEqual({ Ti: '2a', O: '4f' });
        expect(letters(describe0(redescribe(STRUCTURES.rutile(), AXIS_SETTINGS.cab)))).toEqual(reference);
    });

    it('withholds letters for a lower bound or an unnamed class', () => {
        const found = describe0(redescribe(STRUCTURES.perovskite(), [[2, 0, 0], [0, 2, 0], [0, 0, 1]]));
        expect(found.spaceGroupNumber).toBeNull();
        for (const o of found.orbits) expect(o.wyckoff).toBeNull();
    });
});

describe('special positions of P222_1 and I2_12_12_1 at the ITA origin', () => {
    it('tabulates all 1731 ITA positions', () => {
        const rows = Object.values(WYCKOFF_DATA).reduce((n, packed) => n + packed.split(';').length, 0);
        expect(rows).toBe(1731);
        expect(WYCKOFF_DATA[17].split(';').map((r) => r.split(':')[0]).join('')).toBe('abcde');
        expect(WYCKOFF_DATA[24].split(';').map((r) => r.split(':')[0]).join('')).toBe('abcd');
    });

    // ITA #17: (x,y,z) (-x,-y,z+1/2) (-x,y,-z+1/2) (x,-y,-z). 2a x,0,0 · 2c 0,y,1/4 · 4e.
    it('labels P222_1 2a, 2c and 4e', () => {
        const ops = [(x, y, z) => [x, y, z], (x, y, z) => [-x, -y, z + 0.5], (x, y, z) => [-x, y, -z + 0.5], (x, y, z) => [x, -y, -z]];
        const basis = orbits(ops, [['A', 0.21, 0, 0], ['C', 0, 0.31, 0.25], ['G', 0.137, 0.213, 0.061]]);
        const found = describe0({ A: [[8.1, 0, 0], [0, 5.3, 0], [0, 0, 6.7]], basis });
        expect(found.spaceGroup).toBe('P222_1');
        expect(letters(found)).toEqual({ A: '2a', C: '2c', G: '4e' });
    });

    // ITA #24: (x,y,z) (-x+1/2,-y,z+1/2) (-x,y+1/2,-z+1/2) (x+1/2,-y+1/2,-z) + I. 4a x,0,1/4 · 4c 0,1/4,z.
    it('labels I2_12_12_1 4a, 4c and 8d', () => {
        const gens = [(x, y, z) => [x, y, z], (x, y, z) => [-x + 0.5, -y, z + 0.5], (x, y, z) => [-x, y + 0.5, -z + 0.5], (x, y, z) => [x + 0.5, -y + 0.5, -z]];
        const ops = gens.flatMap((g) => [g, (x, y, z) => g(x, y, z).map((v) => v + 0.5)]);
        const basis = orbits(ops, [['A', 0.21, 0, 0.25], ['C', 0, 0.25, 0.37], ['G', 0.137, 0.213, 0.061]]);
        const found = describe0({ A: [[8.1, 0, 0], [0, 5.3, 0], [0, 0, 6.7]], basis });
        expect(found.spaceGroup).toBe('I2_12_12_1');
        expect(letters(found)).toEqual({ A: '4a', C: '4c', G: '8d' });
    });
});
