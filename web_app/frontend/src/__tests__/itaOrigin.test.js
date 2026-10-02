// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The shift onto ITA's standard origin (itaOrigin.js): exact congruence solving by a
// unimodular diagonal form, every equivalent ITA origin considered, the one nearest the
// .rmc6f origin returned, no shift along a polar axis.

import { describe, expect, it } from 'vitest';

import { diagonalForm, itaOriginShift } from '../itaOrigin.js';
import { itaOperations } from '../itaOperations.js';
import { SPACE_GROUP_FIXTURES, cellForSystem } from './fixtures/spaceGroups.js';

const cyc = (x) => x - Math.round(x);
const mul = (X, Y) => X.map((row) => Y[0].map((_, j) => row.reduce((s, v, k) => s + v * Y[k][j], 0)));
const det = (m) => m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

// The group described about the origin p0 (x_new = x − p0): t_new = t + (R − I)·p0.
const describedAbout = (ops, p0) => ops.map(({ R, t }) => ({
    R, t: t.map((v, i) => v + R[i][0] * p0[0] + R[i][1] * p0[1] + R[i][2] * p0[2] - p0[i]),
}));
const shifted = (ops, p) => ops.map(({ R, t }) => ({
    R, t: t.map((v, i) => v + p[i] - (R[i][0] * p[0] + R[i][1] * p[1] + R[i][2] * p[2])),
}));
const sameGroup = (a, b) => a.length === b.length && a.every((o) => b.some((q) => q.R.flat().join() === o.R.flat().join()
    && o.t.every((v, i) => Math.abs(cyc(v - q.t[i])) < 1e-9)));

describe('diagonalForm', () => {
    it('diagonalizes by unimodular operations', () => {
        const M = [[2, -1, 0], [1, 1, 0], [0, 0, 0], [1, 1, 0], [-1, 2, 0], [0, 0, 2]];
        const { D, U, V, rank } = diagonalForm(M);
        expect(mul(mul(U, M), V)).toEqual(D);
        expect(Math.abs(det(V))).toBe(1);
        expect(rank).toBe(3);
        D.forEach((row, i) => row.forEach((v, j) => { if (i !== j) expect(v).toBe(0); }));
    });
});

describe('itaOriginShift', () => {
    it('brings every one of the 230 groups back to ITA from a random origin', () => {
        let s = 11;
        const rand = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
        for (const f of SPACE_GROUP_FIXTURES) {
            const ita = itaOperations(f.number);
            const ops = describedAbout(ita, [rand(), rand(), rand()]);
            const found = itaOriginShift(ops, f.number, cellForSystem(f.system));
            expect(found, f.symbol).not.toBeNull();
            expect(sameGroup(shifted(ops, found.shift), ita), f.symbol).toBe(true);
        }
    });

    it('returns the ITA origin nearest the given offset, and none along a polar axis', () => {
        const A = cellForSystem('cubic');
        // Pm-3m about (0.49, 0.5, 0.51): its origin is at (0.51, 0.5, 0.49) and the body
        // centre, an equivalent ITA origin, at (0.01, 0, -0.01); moving onto that one is the
        // shift (-0.01, 0, 0.01).
        const pm3m = describedAbout(itaOperations(221), [0.49, 0.5, 0.51]);
        const found = itaOriginShift(pm3m, 221, A);
        found.total.forEach((v, i) => expect(Math.abs(v - [-0.01, 0, 0.01][i])).toBeLessThan(1e-9));
        // P4mm (polar c) about (0.2, 0.3, 0.37): x, y move onto the 4-fold axis, z stays.
        const p4mm = describedAbout(itaOperations(99), [0.2, 0.3, 0.37]);
        const polar = itaOriginShift(p4mm, 99, cellForSystem('tetragonal'));
        expect(polar.shift[2]).toBeCloseTo(0, 12);
        expect(sameGroup(shifted(p4mm, polar.shift), itaOperations(99))).toBe(true);
    });

    it('declines a set whose rotations or centring are not the group\'s', () => {
        const A = cellForSystem('cubic');
        expect(itaOriginShift(itaOperations(221), 225, A)).toBeNull();   // P, not F
        expect(itaOriginShift(itaOperations(195), 221, A)).toBeNull();   // a subgroup's rotations
        expect(itaOriginShift([], 221, A)).toBeNull();
    });
});
