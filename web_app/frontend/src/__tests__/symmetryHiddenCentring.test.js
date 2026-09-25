// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A centred lattice on a primitive cell whose unique or principal axis is a basis vector.
// The element DIRECTIONS of such a cell look conventional — the 2-fold of C2 on
// ((a+b)/2, b, c) runs along the second basis vector, the 4-fold of I4 on
// ((a+b+c)/2, (−a+b+c)/2, c) along the third — but the other basis vectors lean on the
// axis, and the cell is primitive where the conventional one is centred. Named in it, the
// centring letter came from its (primitive) pure translations: C2 read P2 (No. 3), I4 read
// P4 (No. 75), R3 read P3 (No. 143), each with the P group's Wyckoff letters, and every
// monoclinic rung of a hexagonal ladder was a P group (a 2-fold along a_hex in a
// P-hexagonal lattice has the orthohexagonal C lattice: C2, not P2).

import { describe, it, expect } from 'vitest';

import { spaceGroupAtTolerance, symmetryLadder } from '../symmetry.js';
import { describeSymmetry, orbitLabel } from '../symmetryModel.js';
import { SPACE_GROUP_FIXTURES, closeGroup, orbitBasis, cellForSystem, structureFor } from './fixtures/spaceGroups.js';
import { STRUCTURES, redescribe, orbits, withNoise } from './fixtures/symmetryStructures.js';

const byNumber = new Map(SPACE_GROUP_FIXTURES.map((f) => [f.number, f]));
const third = 1 / 3;
// Primitive cells (rows = new basis vectors in the conventional fractional basis) that keep
// the unique / principal axis as a basis vector.
const PRIMITIVE = {
    C: [[0.5, 0.5, 0], [0, 1, 0], [0, 0, 1]],
    I: [[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [0, 0, 1]],
    R: [[2 * third, third, third], [-third, third, third], [0, 0, 1]],
};

const labels = (st) => {
    const d = describeSymmetry({ latticeVectors: st.A, supercell: [1, 1, 1], basis: st.basis }, 0.05);
    return { symbol: d.spaceGroup, number: d.spaceGroupNumber, sites: Object.fromEntries(d.orbits.map((o) => [o.element, orbitLabel(o)])) };
};

describe('a centred group on a primitive cell is named in its centred cell', () => {
    it.each([
        [5, 'C'], [8, 'C'], [9, 'C'], [12, 'C'], [15, 'C'],
        [79, 'I'], [80, 'I'], [82, 'I'], [87, 'I'],
        [146, 'R'], [148, 'R'],
    ])('No. %i on the primitive %s-lattice cell', (number, lattice) => {
        const f = byNumber.get(number);
        const ops = closeGroup(f);
        const conventional = {
            A: cellForSystem(f.system),
            basis: [...orbitBasis(ops, 'A', [0.137, 0.213, 0.061]), ...orbitBasis(ops, 'B', [0, 0, 0])],
        };
        const primitive = redescribe(conventional, PRIMITIVE[lattice]);
        expect(primitive.basis.length).toBeLessThan(conventional.basis.length);
        const want = labels(conventional);
        expect(want).toMatchObject({ symbol: f.symbol, number });
        // Same group, same number and the same Wyckoff labels (multiplicity + letter of the
        // standard cell) as the conventional description.
        expect(labels(primitive)).toEqual(want);
    });

    it('labels the inversion centres of P2/c on the sheared cell (a+b, b, c) as in the conventional cell', () => {
        // ITA No. 13: 2a 0,0,0 · 2b ½,½,0 · 2c 0,½,0 · 2d ½,0,0 · 4g general. In the sheared
        // cell X at ½,0,0 reads ½,½,0 — 2b.
        const ops = [(x, y, z) => [x, y, z], (x, y, z) => [-x, y, -z + 0.5], (x, y, z) => [-x, -y, -z], (x, y, z) => [x, -y, z + 0.5]];
        const st = { A: cellForSystem('monoclinic'), basis: orbits(ops, [['Ba', 0, 0, 0], ['X', 0.5, 0, 0], ['O', 0.21, 0.33, 0.07]]) };
        const want = { symbol: 'P2/c', number: 13, sites: { Ba: '2a', X: '2d', O: '4g' } };
        expect(labels(st)).toEqual(want);
        expect(labels(redescribe(st, [[1, 1, 0], [0, 1, 0], [0, 0, 1]]))).toEqual(want);
    });

    it('reads P2_1/c letters in the reduced cell choice when the given cell is not one', () => {
        // ITA No. 14: 2a 0,0,0 · 2b ½,0,0 · 2c 0,0,½ · 2d ½,0,½. Cell choices a' = a + k·c
        // all spell P2_1/c and label the inversion centres differently; the finder picks the
        // shortest a and c (with β non-acute) among the cells it builds, which here is the
        // conventional cell (7.1, 5.3, 9.7 Å, β = 104°).
        const ops = [(x, y, z) => [x, y, z], (x, y, z) => [-x, y + 0.5, -z + 0.5], (x, y, z) => [-x, -y, -z], (x, y, z) => [x, -y + 0.5, z + 0.5]];
        const st = { A: cellForSystem('monoclinic'), basis: orbits(ops, [['Ba', 0, 0, 0], ['X', 0.5, 0, 0], ['O', 0.21, 0.33, 0.07]]) };
        const want = { symbol: 'P2_1/c', number: 14, sites: { Ba: '2a', X: '2b', O: '4e' } };
        expect(labels(st)).toEqual(want);
        expect(labels(redescribe(st, [[1, 1, 0], [0, 1, 0], [0, 0, 1]]))).toEqual(want);   // sheared
        expect(labels(redescribe(st, [[1, 0, 0], [0, 1, 0], [1, 0, 1]]))).toEqual(want);   // P2_1/n cell
    });
});

// Independent check of a P-named group: in a primitive lattice the component of every lattice
// vector along the principal rotation axis — (1/n)·Σ Rᵏ·g for the highest-order proper
// rotation (a rotoinversion's rotation part) — is itself a lattice vector. In a centred
// lattice it is not: the axial part of (a+b)/2 along b is b/2.
const mm = (A, B) => A.map((r) => [0, 1, 2].map((j) => r[0] * B[0][j] + r[1] * B[1][j] + r[2] * B[2][j]));
const mv = (R, v) => [0, 1, 2].map((i) => R[i][0] * v[0] + R[i][1] * v[1] + R[i][2] * v[2]);
const isI = (R) => R.every((r, i) => r.every((x, j) => x === (i === j ? 1 : 0)));
const det = (m) => m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
const order = (R) => { let M = R; for (let n = 1; n <= 6; n++) { if (isI(M)) return n; M = mm(M, R); } return 0; };
function hiddenCentring(found) {
    if (!found.spaceGroupNumber || !/^P/.test(found.spaceGroup) || !found.setting?.ops) return null;   // P1, P-1: no axis
    const { ops, translations } = found.setting;
    let axis = null;
    for (const { R } of ops) {
        const proper = det(R) === 1 ? R : R.map((row) => row.map((x) => -x));
        const n = order(proper);
        if (n > 1 && (!axis || n > axis.n)) axis = { R: proper, n };
    }
    if (!axis) return null;
    const inLattice = (v) => translations.some((tau) => v.every((x, i) => Math.abs(x - tau[i] - Math.round(x - tau[i])) < 1e-6));
    for (const g of [[1, 0, 0], [0, 1, 0], [0, 0, 1], ...translations]) {
        const p = [0, 0, 0];
        let M = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        for (let k = 0; k < axis.n; k++) { const v = mv(M, g); for (let i = 0; i < 3; i++) p[i] += v[i] / axis.n; M = mm(M, axis.R); }
        if (!inLattice(p)) return `${found.spaceGroup}: axial part of [${g}] is [${p.map((x) => +x.toFixed(3))}]`;
    }
    return null;
}

const ladderNames = (st) => symmetryLadder(st.A, st.basis, 1.0).map((b) => {
    const found = spaceGroupAtTolerance(st.A, st.basis, (b.from + b.to) / 2);
    return { found, name: `${found.spaceGroup}${found.spaceGroupNumber ? ` #${found.spaceGroupNumber}` : ''}` };
});

describe('monoclinic subgroups of a hexagonal lattice are C-centred', () => {
    // The {E, 2[100]} subgroup of noisy P6_3/mmc (C2, No. 5) is in symmetryRepresentatives.test.js.
    it('names the monoclinic rungs of wurtzite and bismuth ladders with C symbols', () => {
        const names = (key) => ladderNames(withNoise(STRUCTURES[key](), 0.03, 3)).map((r) => r.name);
        expect(names('wurtzite')).toContain('Cm #8');
        expect(names('wurtzite')).not.toContain('Pm #6');
        const bismuth = names('bismuth');
        expect(bismuth).toEqual(expect.arrayContaining(['C2 #5', 'C2/m #12']));
        expect(bismuth.filter((n) => /^P2/.test(n))).toEqual([]);
    });

    it('no P-named rung of a noisy trigonal or hexagonal fixture ladder hides a centring', { timeout: 60000 }, () => {
        const wrong = [];
        let checked = 0;
        for (const f of SPACE_GROUP_FIXTURES.filter((g) => g.number >= 143 && g.number <= 194)) {
            for (const { found } of ladderNames(withNoise(structureFor(f), 0.03, f.number))) {
                if (!found.spaceGroupNumber || !/^P/.test(found.spaceGroup) || found.spaceGroupNumber <= 2) continue;
                checked++;
                const why = hiddenCentring(found);
                if (why) wrong.push(`${f.symbol}: ${why}`);
            }
        }
        expect(checked).toBeGreaterThan(100);
        expect(wrong).toEqual([]);
    });
});
