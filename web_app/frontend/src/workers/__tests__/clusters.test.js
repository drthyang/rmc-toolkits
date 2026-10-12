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

describe('a site at fraction 0 (review 2026-10-12 #1)', () => {
    it('pairs each atom with its partner of the same cell when the site mean is a hair below 0', () => {
        // The rigid V2 dimers' first site averages −0.0001 Å: wrapped to
        // 0.9999 it pulled the partner from the next cell, which hid the
        // dimers' correlation and broke edges.
        const parsed = siteDisplacementsFromRmc6f(fixture.cases.find((entry) => entry.name === 'origin-dimer').text);
        expect(parsed.sites[0].siteFractional[0]).toBeGreaterThan(0.99);
        expect(parsed.sites[0].siteCellFractional[0]).toBeCloseTo(-2e-5, 9);
        const result = clusterDistortion(parsed, 'V', [2.0, 2.5]);
        expect(result.templates[0].offsets).toEqual([[0, 0, 0], [0, 0, 0]]);
        expect(result.outliers.brokenEdgeFraction).toBe(0);
        expect(result.z.breathingSd).toBeLessThan(-10);
    });
});

describe('templates with their own edges (review 2026-10-12 #2)', () => {
    it('measures each chain with its own template\'s edges', () => {
        // Chain B (bonded 0-1, 0-2) measured with chain A's edges (0-1, 1-2)
        // counted its 4 Å end-to-end span as an edge: ladderWidth 1.03 Å.
        const parsed = siteDisplacementsFromRmc6f(fixture.cases.find((entry) => entry.name === 'chains').text);
        const result = clusterDistortion(parsed, 'Mo', [1.8, 2.2]);
        expect(result.templates.map((template) => template.edges)).toEqual([[[0, 1], [1, 2]], [[0, 1], [0, 2]]]);
        expect(result.observed.ladderWidth).toBeLessThan(0.2);
        expect(result.outliers.brokenEdgeFraction).toBe(0);
    });
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
