// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Parity of the static-mode modulation engine against the Python reference
// (rmc_toolkits/modulations.py) via the golden fixture written by
// tests/generate_modulations_fixture.py, plus the worker's request.

import { describe, expect, it } from 'vitest';
import { displacementModulations, zoneBoundaryPoints } from '../modulations';
import { siteDisplacementsFromRmc6f } from '../pcaKde';
import { handlePcaMessage } from '../pcaKdeWorker';
import fixture from '../../__tests__/fixtures/modulations_fixture.json';

// Integers, strings and booleans exactly; other floats to 1e-9.
const expectClose = (actual, expected, path) => {
    if (expected === null || typeof expected !== 'object') {
        if (typeof expected === 'number' && !Number.isInteger(expected)) expect(actual, path).toBeCloseTo(expected, 9);
        else expect(actual, path).toEqual(expected);
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

describe('Python parity (modulations_fixture.json)', () => {
    it('halves only even axes', () => {
        expect(zoneBoundaryPoints([4, 4, 4])).toHaveLength(7);
        expect(zoneBoundaryPoints([3, 4, 5])).toEqual([[0, 1, 0]]);
    });

    for (const testCase of fixture.cases) {
        const parsed = siteDisplacementsFromRmc6f(testCase.text);
        testCase.elements.forEach((element, index) => {
            it(`${testCase.name}: ${element ?? 'every element'}`, () => {
                expectClose(displacementModulations(parsed, element), testCase.expected[index], 'result');
            });
        });
    }
});

describe('the worker\'s modulations request', () => {
    const text = fixture.cases[0].text;
    const ask = (params) => handlePcaMessage({ kind: 'modulations', ...params }, async () => text);

    it('answers every element by default and one when asked', async () => {
        expect((await ask({})).elements.map((entry) => entry.element)).toEqual(['Ga', 'Nb', 'Se']);
        expect((await ask({ element: 'all' })).elements).toHaveLength(3);
        const nb = await ask({ element: 'nb' });
        expect(nb.elements.map((entry) => entry.element)).toEqual(['Nb']);
        expect(nb.elements[0].strongest[0].q).toEqual([0, 0, 0.5]);
        await expect(ask({ element: 'Xx' })).rejects.toThrow('no site holds Xx');
    });
});
