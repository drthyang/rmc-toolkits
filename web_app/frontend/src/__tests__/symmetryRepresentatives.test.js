// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// An operation {R|t} stands for its whole coset {R | t + ℓ}, ℓ a lattice vector, and which
// representative the finder happens to hold is an accident: a refined translation is
// wrapped into [0, 1), so noise turns an exact 0 into 0.9995. Where the lattice projects
// onto an axis or plane in a FRACTION of its period, representatives differ in kind: the
// [100] 2-fold of a hexagonal cell and the 2_1 half a cell away, a cubic [111] 3-fold
// and a 3_1, a mirror on a tetragonal diagonal and an n-glide. The symbol must not depend
// on which one was held: a noisy hexagonal subgroup {E, 2[100]} was named P2_1 (No. 4) or
// P2 (No. 3) by the representative — it is C2 (No. 5), whose orthohexagonal C cell holds
// the 2 and the 2_1 alike (symmetryHiddenCentring.test.js has the centring side).

import { describe, it, expect } from 'vitest';

import { spaceGroupAtTolerance } from '../symmetry.js';
import { cosetElements, classifyElement, screwIndex } from '../spaceGroupSymbol.js';
import { SPACE_GROUP_FIXTURES, closeGroup, structureFor } from './fixtures/spaceGroups.js';
import { withNoise } from './fixtures/symmetryStructures.js';

const hex2along100 = [[1, -1, 0], [0, -1, 0], [0, 0, -1]];     // x−y, −y, −z
const cubic3along111 = [[0, 0, 1], [1, 0, 0], [0, 1, 0]];      // z, x, y
const diagonalMirror = [[0, 1, 0], [1, 0, 0], [0, 0, 1]];      // y, x, z
const labels = (R, t) => cosetElements(R, t).map((e) => e.label).sort();

describe('every lattice representative of an operation is offered', () => {
    it('a hexagonal [100] 2-fold coset holds the rotation whichever representative is held', () => {
        for (const t of [[0, 0, 0], [0, 1, 0], [1, 1, 0], [0.9995, 0.999, 0.0002]]) {
            expect(labels(hex2along100, t)).toEqual(['2', '2_1']);
        }
    });

    it('a cubic [111] 3-fold coset holds the rotation and both screws', () => {
        expect(labels(cubic3along111, [1, 0, 0])).toEqual(['3', '3_1', '3_2']);
    });

    it('a diagonal mirror coset holds the mirror and the n-glide', () => {
        expect(labels(diagonalMirror, [0.9996, 0.0003, 0])).toEqual(['m', 'n']);
    });

    it('a representative with no fractional projection keeps its kind', () => {
        const fourfold = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        expect(labels(fourfold, [0, 0, 0.25])).toEqual(['4_1']);
        expect(labels(fourfold, [0.9997, 0.0002, 0.2499])).toEqual(['4_1']);
        expect(labels([[1, 0, 0], [0, 1, 0], [0, 0, -1]], [0.5, 0, 0])).toEqual(['a']);
    });

    it('reads the screw from the unwrapped intrinsic translation, whatever the axis signs', () => {
        // A 3_1 along [1 −1 −1] and the same screw after the proper change of basis
        // S = diag(1, −1, −1), which carries the axis to [1 1 1]. Wrapped component by
        // component, the intrinsic translation (1, −1, −1)/3 becomes (1/3, 2/3, 2/3), whose
        // projection onto the axis is −1/3: the screw read as the other hand.
        const R = [[0, -1, 0], [0, 0, 1], [-1, 0, 0]];            // −y, z, −x
        const S = [[1, 0, 0], [0, -1, 0], [0, 0, -1]];
        const Rs = mm(mm(S, R), S);                               // S⁻¹ = S
        const t = [1 / 3, -1 / 3, -1 / 3];
        expect(classifyElement(R, [0, 0, 0]).label).toBe('3');
        expect(classifyElement(R, t).kind).toBe('screw');
        expect(classifyElement(R, t).label).toBe(classifyElement(Rs, mv(S, t)).label);
        expect(screwIndex(R, t, [1, -1, -1], 3)).toBe(screwIndex(Rs, mv(S, t), [1, 1, 1], 3));
    });
});

// Element-type invariant of a group, independent of cell, origin and representative: for
// each coset of the pure translations, its rotation type (det, trace) and whether it
// contains an operation with no intrinsic translation (a rotation or mirror rather than
// only screws or glides), counted per translation.
const mm = (A, B) => A.map((r) => [0, 1, 2].map((j) => r[0] * B[0][j] + r[1] * B[1][j] + r[2] * B[2][j]));
const mv = (A, v) => A.map((r) => r[0] * v[0] + r[1] * v[1] + r[2] * v[2]);
const isI = (R) => R.every((r, i) => r.every((x, j) => x === (i === j ? 1 : 0)));
function invariant(ops) {
    const T = ops.filter((o) => isI(o.R)).map((o) => o.t);
    const lattice = [];
    for (let a = -2; a <= 2; a++) for (let b = -2; b <= 2; b++) for (let c = -2; c <= 2; c++) for (const t of T) lattice.push([a + t[0], b + t[1], c + t[2]]);
    const counts = {};
    for (const { R, t } of ops) {
        let n = 1, P = R;
        while (!isI(P) && n < 7) { P = mm(P, R); n++; }
        const proj = (v) => {
            const s = [0, 0, 0];
            let Rk = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
            for (let k = 0; k < n; k++) { const u = mv(Rk, v); for (let i = 0; i < 3; i++) s[i] += u[i] / n; Rk = mm(Rk, R); }
            return s;
        };
        const w = proj(t);
        const pure = lattice.some((L) => proj(L).every((x, i) => Math.abs(x - w[i]) < 0.03));
        const det = R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1]) - R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0]) + R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]);
        const key = `${det}/${R[0][0] + R[1][1] + R[2][2]}/${pure ? 'p' : 's'}`;
        counts[key] = (counts[key] ?? 0) + 1 / T.length;
    }
    return Object.keys(counts).sort().map((k) => `${k}:${counts[k].toFixed(3)}`).join(' ');
}

describe('noisy subgroups are named after the elements they hold', () => {
    const byNumber = new Map(SPACE_GROUP_FIXTURES.map((f) => [f.number, f]));
    const trigonalToCubic = SPACE_GROUP_FIXTURES.filter((f) => f.number >= 143);

    it('names the {E, 2[100]} subgroup of noisy P6_3/mmc C2 (No. 5), not P2 or P2_1', () => {
        // A 2-fold along a_hex of a P-hexagonal lattice: the conventional monoclinic cell is
        // the orthohexagonal (a, a+2b, c) cell, C-centred, twice the hexagonal volume.
        const s = withNoise(structureFor(byNumber.get(194)), 0.01, 2);
        const found = spaceGroupAtTolerance(s.A, s.basis, 0.03);
        expect(found.nSpace).toBe(2);
        expect(found).toMatchObject({ spaceGroup: 'C2', spaceGroupNumber: 5 });
        expect(found.setting.ratio).toBeCloseTo(2, 9);
    });

    it.each([[0.01, 0.03], [0.01, 0.045], [0.02, 0.08]])('σ = %f Å at τ = %f Å: every named trigonal, hexagonal and cubic subgroup', { timeout: 60000 }, (sigma, tol) => {
        const wrong = [];
        for (const f of trigonalToCubic) {
            const s = withNoise(structureFor(f), sigma, 2);
            const found = spaceGroupAtTolerance(s.A, s.basis, tol);
            if (!found.spaceGroupNumber) continue;
            const named = invariant(closeGroup(byNumber.get(found.spaceGroupNumber)));
            if (invariant(found.ops) !== named) wrong.push(`${f.symbol} → ${found.spaceGroup} (No. ${found.spaceGroupNumber})`);
        }
        expect(wrong).toEqual([]);
    });
});
