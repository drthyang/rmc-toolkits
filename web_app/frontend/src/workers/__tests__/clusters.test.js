// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Parity of the static-mode cluster engine against the Python reference
// (rmc_toolkits/clusters.py) via the golden fixture written by
// tests/generate_clusters_fixture.py, plus the worker's request checks.

import { describe, expect, it } from 'vitest';
import { clusterDistortion, mix32, shufflePermutation } from '../clusters';
import { siteDisplacementsFromRmc6f } from '../pcaKde';
import { handlePcaMessage } from '../pcaKdeWorker';
import fixture from '../../__tests__/fixtures/clusters_fixture.json';

// Integers, strings and booleans must match exactly; every other float to
// 1e-9 (summation order differs between numpy and the port).

const expectClose = (actual, expected, path) => {
    if (expected === null || typeof expected !== 'object') {
        if (typeof expected === 'number' && !Number.isInteger(expected)) {
            expect(actual, path).toBeCloseTo(expected, 9);
        } else {
            expect(actual, path).toEqual(expected);
        }
        return;
    }
    if (Array.isArray(expected)) {
        expect(actual, path).toHaveLength(expected.length);
        expected.forEach((value, index) => expectClose(actual[index], value, `${path}[${index}]`));
        return;
    }
    expect(Object.keys(actual).sort(), path).toEqual(Object.keys(expected).sort());
    Object.entries(expected).forEach(([key, value]) => expectClose(actual[key], value, `${path}.${key}`));
};

describe('Python parity (clusters_fixture.json)', () => {
    it('hashes and permutes like the Python engine', () => {
        Object.entries(fixture.hashes.mix32).forEach(([value, hashed]) => expect(mix32(Number(value))).toBe(hashed));
        const { n, shuffle, template, member, seed, order } = fixture.hashes.permutation;
        expect(shufflePermutation(n, shuffle, template, member, seed)).toEqual(order);
    });

    for (const testCase of fixture.cases) {
        const parsed = siteDisplacementsFromRmc6f(testCase.text);
        testCase.queries.forEach(({ element, window }, index) => {
            it(`${testCase.name}: ${element} ${window.join('–')} Å`, () => {
                expectClose(clusterDistortion(parsed, element, window), testCase.expected[index], 'result');
            });
        });
    }
});

describe('the worker\'s clusters request', () => {
    const text = fixture.cases[0].text;
    const ask = (params) => handlePcaMessage({ kind: 'clusters', ...params }, async () => text);

    it('answers with the engine payload', async () => {
        const result = await ask({ element: 'nb', rMin: '2.5', rMax: 3.5 });
        expect(result.element).toBe('Nb');
        expect(result.kind).toMatchObject({ members: 4, edges: 6, clusters: 64 });
    });

    it('rejects missing, inverted and oversized windows with the route\'s messages', async () => {
        await expect(ask({ element: 'Nb', rMin: 2.5 })).rejects.toThrow('rMax required');
        await expect(ask({ element: 'Nb', rMin: 2.5, rMax: 20 })).rejects.toThrow('rMax must be <= 8, got 20');
        await expect(ask({ element: 'Nb', rMin: 3.5, rMax: 2.5 })).rejects.toThrow('window must be 0 <= min < max');
        await expect(ask({ element: 'Xx', rMin: 2.5, rMax: 3.5 })).rejects.toThrow('no site holds Xx');
    });
});
