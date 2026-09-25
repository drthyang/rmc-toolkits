// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Regression tests for the 1.0 audit of the displacement-direction engine —
// the JS twin of tests/test_orientation_fixes.py. Engine-level expectations
// are shared verbatim with the Python suite so the two engines cannot drift.

import { describe, it, expect } from 'vitest';

import {
    MIN_FREQUENCY,
    orientationHistogram,
    recommendedFrequency
} from '../orientation.js';

// Deterministic Gaussian sampler (same construction as orientation.test.js).
const makeRng = (seed) => {
    let value = seed >>> 0;
    const next = () => {
        value += 0x6d2b79f5;
        let mixed = value;
        mixed = Math.imul(mixed ^ (mixed >>> 15), mixed | 1);
        mixed ^= mixed + Math.imul(mixed ^ (mixed >>> 7), mixed | 61);
        return ((mixed ^ (mixed >>> 14)) >>> 0) / 4294967296;
    };
    return () => {
        let u = 0;
        let v = 0;
        while (u === 0) u = next();
        while (v === 0) v = next();
        return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
    };
};

const cloud = (n = 200, seed = 0) => {
    const gauss = makeRng(seed);
    const points = [];
    for (let i = 0; i < n; i += 1) points.push([gauss(), gauss(), gauss()]);
    return points;
};

// orientation.numerics.5/.15/.29, orientation.parity.12/.18
describe('non-finite input', () => {
    it('rejects a NaN row with a clear message in both frames', () => {
        for (const frame of ['cartesian', 'pca']) {
            const vectors = cloud();
            vectors[5][0] = NaN;
            expect(() => orientationHistogram(vectors, { frame, frequency: 3 })).toThrow(/non-finite.*row 5/);
        }
    });

    it('rejects an inf row with a clear message in both frames', () => {
        for (const frame of ['cartesian', 'pca']) {
            const vectors = cloud();
            vectors[7][2] = -Infinity;
            expect(() => orientationHistogram(vectors, { frame, frequency: 3 })).toThrow(/1 non-finite.*row 7/);
        }
    });

    it('rejects non-finite or out-of-range options', () => {
        const vectors = cloud();
        for (const options of [
            { smoothing: NaN },
            { smoothing: -1 },
            { targetPerCell: NaN },
            { targetPerCell: 0 },
            { minAmplitude: NaN },
            { frequency: NaN }
        ]) {
            expect(() => orientationHistogram(vectors, options)).toThrow();
        }
    });

    it('rejects invalid recommendedFrequency bounds', () => {
        expect(() => recommendedFrequency(1000, { maxFrequency: 0 })).toThrow();
        expect(() => recommendedFrequency(NaN)).toThrow();
        expect(() => recommendedFrequency(1000, { targetPerCell: Infinity })).toThrow();
        expect(recommendedFrequency(1000, { maxFrequency: MIN_FREQUENCY })).toBe(MIN_FREQUENCY);
    });
});

// orientation.physics.10/.25, orientation.numerics.30 — shared verbatim with
// RECOMMENDED_FREQUENCY_PINS in tests/test_orientation_fixes.py.
const RECOMMENDED_FREQUENCY_PINS = [
    [0, 1], [294, 1], [300, 1], [503, 1], [504, 2], [774, 2], [1000, 2], [1103, 2],
    [1104, 3], [12000, 9], [12023, 9], [12024, 10], [10000000, 24]
];

describe('recommendedFrequency floors to the target occupancy', () => {
    it('never drops below 12 points per cell (except at the dodecahedron floor)', () => {
        for (let n = 1; n < 20000; n += 7) {
            const frequency = recommendedFrequency(n);
            if (frequency > MIN_FREQUENCY) expect(n / (10 * frequency * frequency + 2)).toBeGreaterThanOrEqual(12);
            if (frequency < 24) expect(n / (10 * (frequency + 1) ** 2 + 2)).toBeLessThan(12);
        }
    });

    it('matches the Python pins', () => {
        RECOMMENDED_FREQUENCY_PINS.forEach(([n, expected]) => expect(recommendedFrequency(n)).toBe(expected));
        expect(recommendedFrequency(1000, { targetPerCell: 5 })).toBe(4);
        expect(recommendedFrequency(5000, { maxFrequency: 3 })).toBe(3);
    });
});
