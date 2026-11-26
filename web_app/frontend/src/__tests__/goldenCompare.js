// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Compare a browser engine's result with a Python golden (a fixture written by
// a tests/generate_*_fixture.py) -- the rule tests/golden_compare.py applies to
// the regenerated fixture on the Python side:
//
// - strings, booleans and null exactly;
// - a number that is an integer on both sides exactly (counts, cell indices,
//   frequencies: a point in another cell is a different result, not round-off);
// - any other number to RTOL relative plus ATOL absolute (the engines sum in a
//   different order, and the JS eigensolver and incomplete gamma are not
//   LAPACK's or SciPy's), or to RTOL relative alone for a `relative` key
//   (p-values reach 1e-142, where any absolute floor would check nothing);
// - an eigenvector set (a `bases` key, given the descending eigenvalues its
//   rows belong to) eigenspace by eigenspace, to AXIS_TOL: an axis whose
//   eigenvalue is separated from its neighbours by more than DEGENERATE_GAP of
//   the largest up to its sign, and a (numerically) degenerate group of axes by
//   the projector onto the space it spans -- within that space any orthonormal
//   basis is an eigenbasis, and LAPACK and the JS Jacobi solver pick different
//   ones. A separated axis is fixed to eps / DEGENERATE_GAP ~ 2e-10, inside
//   AXIS_TOL;
// - objects by their whole key set.

import { expect } from 'vitest';

export const RTOL = 1e-9;
export const ATOL = 1e-12;
export const AXIS_TOL = 1e-9;
export const DEGENERATE_GAP = 1e-6;

const isArrayLike = (value) => Array.isArray(value) || ArrayBuffer.isView(value);

const expectNumber = (actual, expected, path, relative) => {
    expect(typeof actual, path).toBe('number');
    if (Number.isInteger(expected) && Number.isInteger(actual)) {
        // + 0 folds -0 into 0: the same integer.
        expect(actual + 0, path).toBe(expected + 0);
        return;
    }
    const bound = RTOL * Math.abs(expected) + (relative ? 1e-300 : ATOL);
    expect(Math.abs(actual - expected), `${path}: ${actual} vs ${expected}`).toBeLessThanOrEqual(bound);
};

/** Groups of row indices whose eigenvalues (descending) are degenerate. */
export const eigenspaces = (spectrum) => {
    const scale = Math.max(Math.abs(spectrum[0]), 1e-300);
    const groups = [[0]];
    for (let i = 1; i < spectrum.length; i += 1) {
        if (spectrum[i - 1] - spectrum[i] > DEGENERATE_GAP * scale) groups.push([i]);
        else groups[groups.length - 1].push(i);
    }
    return groups;
};

const projector = (rows) => [0, 1, 2].map((a) => [0, 1, 2].map(
    (b) => rows.reduce((sum, row) => sum + row[a] * row[b], 0)
));

const expectEigenbasis = (actual, expected, spectrum, path) => {
    expect(isArrayLike(actual), `${path} is an array`).toBe(true);
    expect(actual.length, `${path} length`).toBe(expected.length);
    eigenspaces(spectrum).forEach((group) => {
        let misfit;
        if (group.length === 1) {
            const [i] = group;
            const differ = (sign) => Math.max(...expected[i].map((value, k) => Math.abs(sign * actual[i][k] - value)));
            misfit = Math.min(differ(1), differ(-1));
        } else {
            const want = projector(group.map((i) => expected[i]));
            const got = projector(group.map((i) => actual[i]));
            misfit = Math.max(...want.flatMap((row, a) => row.map((value, b) => Math.abs(got[a][b] - value))));
        }
        expect(misfit, `${path} rows ${group.join(',')}`).toBeLessThanOrEqual(AXIS_TOL);
    });
};

export const expectGolden = (actual, expected, path, options = {}) => {
    const { relative = new Set(), bases = {}, key = null } = options;
    if (expected === null) {
        expect(actual, path).toBeNull();
        return;
    }
    if (typeof expected === 'number') {
        expectNumber(actual, expected, path, relative.has(key));
        return;
    }
    if (typeof expected !== 'object') {
        expect(actual, path).toBe(expected);
        return;
    }
    if (Array.isArray(expected)) {
        if (key !== null && bases[key]) {
            expectEigenbasis(actual, expected, bases[key], path);
            return;
        }
        expect(isArrayLike(actual), `${path} is an array`).toBe(true);
        expect(actual.length, `${path} length`).toBe(expected.length);
        expected.forEach((value, i) => expectGolden(actual[i], value, `${path}[${i}]`, { relative, key }));
        return;
    }
    expect(actual !== null && typeof actual === 'object', `${path} is an object`).toBe(true);
    expect(Object.keys(actual).sort(), `${path} keys`).toEqual(Object.keys(expected).sort());
    Object.entries(expected).forEach(([name, value]) => {
        expectGolden(actual[name], value, `${path}.${name}`, { relative, bases, key: name });
    });
};
