// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Axis-domain regression tests. niceDomain once did Math.min(...finite), which
// pushes every point onto the argument stack and blew up on long RMCProfile
// series; these pin both the size it used to fail at and the range it returns.

import { describe, expect, it } from 'vitest';
import { niceDomain } from '../plotDomain';

// Comfortably past V8's ~124k argument-stack limit.
const LARGE = 2e5;

// Extremes planted mid-array so a single pass has to scan the whole series.
const MIN = -37.5;
const MAX = 91.25;

const largeSeries = () => {
    const values = new Array(LARGE);
    for (let i = 0; i < LARGE; i += 1) values[i] = Math.sin(i) * 10;
    values[LARGE / 2] = MIN;
    values[LARGE / 4] = MAX;
    return values;
};

describe('niceDomain', () => {
    it('handles a series long enough to overflow the argument stack', () => {
        const values = largeSeries();
        // Guards the fixture: if this stops throwing, LARGE no longer covers the bug.
        expect(() => Math.min(...values)).toThrow(RangeError);

        const pad = (MAX - MIN) * 0.05;
        expect(niceDomain(values)).toEqual([MIN - pad, MAX + pad]);
    });

    it('pads the data extent by 5% on each side', () => {
        expect(niceDomain([0, 10])).toEqual([-0.5, 10.5]);
    });

    it('ignores NaN and Infinity', () => {
        expect(niceDomain([NaN, 2, Infinity, 4, -Infinity])).toEqual(niceDomain([2, 4]));
    });

    it('falls back to [0, 1] with no finite values', () => {
        expect(niceDomain([])).toEqual([0, 1]);
        expect(niceDomain([NaN, Infinity])).toEqual([0, 1]);
    });

    it('widens a flat series so the axis has span', () => {
        expect(niceDomain([5, 5, 5])).toEqual([4 - 0.1, 6 + 0.1]);
    });
});
