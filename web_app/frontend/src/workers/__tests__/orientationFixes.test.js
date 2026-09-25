// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Regression tests for the 1.0 audit of the displacement-direction engine —
// the JS twin of tests/test_orientation_fixes.py. Engine-level expectations
// are shared verbatim with the Python suite so the two engines cannot drift.

import { describe, it, expect } from 'vitest';

import {
    MIN_FREQUENCY,
    assignCells,
    goldbergTiling,
    orientationHistogram,
    recommendedFrequency,
    siteOrientationHistogram
} from '../orientation.js';
import { siteDisplacementsFromRmc6f } from '../pcaKde.js';
import { handlePcaMessage } from '../pcaKdeWorker.js';

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

// orientation.numerics.28 — exact Voronoi ties must resolve centrosymmetrically.
const axisCloud = (copies = 500, length = 0.1) => {
    const points = [];
    for (const axis of [[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1]]) {
        for (let i = 0; i < copies; i += 1) points.push(axis.map((value) => value * length));
    }
    return points;
};

const bodyDiagonalCloud = (copies = 200, length = 0.1) => {
    const points = [];
    for (const a of [1, -1]) {
        for (const b of [1, -1]) {
            for (const c of [1, -1]) {
                for (let i = 0; i < copies; i += 1) {
                    points.push([a, b, c].map((value) => value * length / Math.sqrt(3)));
                }
            }
        }
    }
    return points;
};

// Shared verbatim with AXIS_CELLS_NU3 in tests/test_orientation_fixes.py.
const AXIS_CELLS_NU3 = [36, 15, 16, 86, 71, 62];

describe('centrosymmetric tie-break', () => {
    it('keeps an exactly centrosymmetric axis cloud at zero asymmetry', () => {
        for (const frequency of [1, 2, 3, 5, 6, 9, 10, 11]) {
            const result = orientationHistogram(axisCloud(), { frequency, geometry: false });
            result.counts.forEach((count, cell) => expect(count).toBe(result.counts[result.antipode[cell]]));
            expect(result.antipodalAsymmetry).toBe(0);
        }
    });

    it('keeps a body-diagonal cloud at zero asymmetry', () => {
        for (const frequency of [4, 5, 7, 11]) {
            expect(orientationHistogram(bodyDiagonalCloud(), { frequency, geometry: false }).antipodalAsymmetry).toBe(0);
        }
    });

    it('is exactly inversion-equivariant, ties included', () => {
        const tiling = goldbergTiling(5);
        const directions = [
            ...cloud(2000, 31),
            [1, 0, 0], [0, 1, 0], [0, 0, 1],
            ...bodyDiagonalCloud(1),
            ...tiling.polygons.map((polygon) => polygon[0]),
            ...tiling.centers
        ];
        const plus = assignCells(tiling, directions);
        const minus = assignCells(tiling, directions.map((u) => u.map((value) => -value)));
        minus.forEach((cell, index) => expect(cell).toBe(tiling.antipode[plus[index]]));
    });

    it('pins the axis cells shared with the Python engine', () => {
        const tiling = goldbergTiling(3);
        expect(assignCells(tiling, [[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1]]))
            .toEqual(AXIS_CELLS_NU3);
    });
});

// orientation.parity.17 — cyclic neighbour / polygon order must not depend on
// the atan2 branch cut. Shared verbatim with tests/test_orientation_fixes.py.
const NEIGHBOR_CHECKSUMS = { 4: 13107810, 6: 68604861, 10: 524199205, 38: 108147532718 };
const POLYGON_STARTS_NU38 = {
    2052: [0.008313449390497447, 0.01625704168659895, 0.9998332836802503],
    4275: [0.9998332836802503, 0.00831344939049745, -0.016257041686598955]
};

const neighborChecksum = (neighbors) => {
    let total = 0;
    neighbors.forEach((row, cell) => {
        row.forEach((neighbor, slot) => { total += (slot + 1) * (neighbor + 1) * ((cell % 97) + 1); });
    });
    return total;
};

describe('canonical cyclic order', () => {
    it('starts every neighbour row at its smallest index', () => {
        for (const frequency of [3, 6, 10]) {
            goldbergTiling(frequency).neighbors.forEach((row) => {
                const valid = row.filter((value) => value >= 0);
                expect(valid[0]).toBe(Math.min(...valid));
            });
        }
    });

    it('matches the Python neighbour checksums', () => {
        Object.entries(NEIGHBOR_CHECKSUMS).forEach(([frequency, expected]) => {
            expect(neighborChecksum(goldbergTiling(Number(frequency)).neighbors)).toBe(expected);
        });
    });

    it('matches the Python polygon start vertices', () => {
        const tiling = goldbergTiling(38);
        Object.entries(POLYGON_STARTS_NU38).forEach(([cell, vertex]) => {
            vertex.forEach((value, axis) => expect(tiling.polygons[Number(cell)][0][axis]).toBeCloseTo(value, 12));
        });
    });
});

// orientation.parity.19 — element "all" / "" is the pooled default in both
// transports and must not be stamped onto the payload (Python omits the key).
const twoSiteRmc6f = () => {
    const gauss = makeRng(41);
    const lines = [
        'Supercell dimensions 4 4 4',
        'Lattice vectors (Ang):',
        '32 0 0',
        '0 32 0',
        '0 0 32',
        'Atoms:'
    ];
    let atom = 0;
    for (const [element, reference, offset] of [['Ga', 1, 0], ['Se', 2, 0.5]]) {
        for (let ix = 0; ix < 4; ix += 1) {
            for (let iy = 0; iy < 4; iy += 1) {
                for (let iz = 0; iz < 4; iz += 1) {
                    atom += 1;
                    const coord = [ix, iy, iz].map((index) => (index + offset) / 4 + gauss() * 0.003);
                    lines.push(`${atom} ${element} [${reference}] ${coord.map((value) => value.toFixed(10)).join(' ')} ${reference} ${ix} ${iy} ${iz}`);
                }
            }
        }
    }
    return lines.join('\n');
};

describe('element "all" normalisation', () => {
    it('does not stamp element on a pooled library result', () => {
        const parsed = siteDisplacementsFromRmc6f(twoSiteRmc6f());
        for (const element of ['all', '', null]) {
            const result = siteOrientationHistogram(parsed, { element, frequency: 2, geometry: false });
            expect(result.totalPoints).toBe(128);
            expect('element' in result).toBe(false);
        }
        expect(siteOrientationHistogram(parsed, { element: 'Se', frequency: 2, geometry: false }).element).toBe('Se');
    });

    it('normalises element "all" in the worker transport like the Flask route', async () => {
        const text = twoSiteRmc6f();
        const result = await handlePcaMessage(
            { kind: 'orientation', element: 'all', frequency: 2, geometry: false },
            async () => text
        );
        expect(result.totalPoints).toBe(128);
        expect('element' in result).toBe(false);
    });
});
