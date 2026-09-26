// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Regression tests for the 1.0 audit of the displacement-direction engine —
// the JS twin of tests/test_orientation_fixes.py. Engine-level expectations
// are shared verbatim with the Python suite so the two engines cannot drift.

import { describe, it, expect } from 'vitest';

import {
    MAP_TEST_MIN_PAIRS,
    MIN_FREQUENCY,
    assignCells,
    goldbergTiling,
    logGamma,
    normalDeviate,
    normalQuantile,
    mapTest,
    orientationHistogram,
    pearsonNullMoments,
    recommendedFrequency,
    regularizedGamma,
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
// +x and +y are exact two-way ties at nu = 3; since review:4 they go to the
// lower index of the pair (35, 14; previously 36, 15 by the walk's path).
const AXIS_CELLS_NU3 = [35, 14, 16, 85, 70, 62];

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

// <111>, <100> and <110> (26 directions), in the order of _special_directions
// in tests/test_orientation_fixes.py.
const specialDirections = () => {
    const signs = [1, -1];
    const directions = [];
    signs.forEach((a) => signs.forEach((b) => signs.forEach((c) => directions.push([a, b, c]))));
    for (let axis = 0; axis < 3; axis += 1) {
        signs.forEach((sign) => {
            const v = [0, 0, 0];
            v[axis] = sign;
            directions.push(v);
        });
    }
    [[0, 1], [0, 2], [1, 2]].forEach(([i, j]) => signs.forEach((a) => signs.forEach((b) => {
        const v = [0, 0, 0];
        v[i] = a;
        v[j] = b;
        directions.push(v);
    })));
    return directions;
};

// Brute-force reference for the tie rule over all cells: fold into the
// canonical hemisphere, take the lowest index within `tolerance` of the best
// dot product, map folded directions back through the antipode table.
const lowestIndexRule = (tiling, directions, tolerance = 1e-12) => directions.map((direction) => {
    const length = Math.hypot(direction[0], direction[1], direction[2]);
    const u = direction.map((value) => value / length);
    const flip = u[0] < 0 || (u[0] === 0 && (u[1] < 0 || (u[1] === 0 && u[2] < 0)));
    const w = flip ? u.map((value) => -value) : u;
    const dots = tiling.centers.map((c) => c[0] * w[0] + c[1] * w[1] + c[2] * w[2]);
    const best = Math.max(...dots);
    const cell = dots.findIndex((value) => value >= best - tolerance);
    return flip ? tiling.antipode[cell] : cell;
});

// Shared verbatim with SPECIAL_DIRECTION_CELLS in tests/test_orientation_fixes.py.
const SPECIAL_DIRECTION_CELLS = {
    1: [0, 1, 2, 6, 5, 9, 11, 10, 2, 9, 1, 11, 0, 10, 1, 8, 3, 11, 2, 6, 5, 9, 0, 7, 4, 10],
    2: [1, 16, 10, 26, 15, 27, 39, 38, 18, 40, 8, 34, 9, 30, 2, 24, 7, 35, 5, 17, 14, 29, 0, 21, 11, 31],
    3: [5, 43, 49, 77, 28, 56, 83, 91, 35, 85, 14, 70, 16, 62, 6, 47, 57, 89, 7, 45, 26, 68, 1, 40, 75, 80],
    4: [6, 71, 81, 130, 48, 98, 145, 159, 59, 148, 23, 120, 29, 104, 8, 79, 96, 155, 12, 76, 43, 115,
        1, 68, 128, 137],
    10: [33, 405, 460, 761, 257, 561, 858, 983, 320, 880, 110, 690, 155, 590, 287, 525, 919, 927, 128,
        635, 223, 656, 77, 377, 736, 792],
    11: [37, 485, 551, 915, 306, 673, 1033, 1193, 381, 1061, 128, 830, 182, 710, 343, 633, 1110, 1119,
        152, 765, 268, 788, 90, 454, 887, 954]
};
// Shared verbatim with WELL_MODEL_PEAKS in tests/test_orientation_fixes.py.
// Before the tie rule Python put the nu=2 peak at [0.5, 0.309, 0.809] and JS
// at [0.309, 0.809, 0.5].
const WELL_MODEL_PEAKS = {
    2: { peakCell: 1, peakDirection: [0.3090169943749474, 0.8090169943749473, 0.5], occupiedCells: 8 },
    4: { peakCell: 6, peakDirection: [0.42532540417601994, 0.5877852522924731, 0.6881909602355868], occupiedCells: 8 }
};

// Atoms on the eight <111> wells: 30 at +[111], 30 at -[111], 20 at each other.
const wellModel = (length = 0.1) => {
    const rows = [];
    specialDirections().slice(0, 8).forEach((direction) => {
        const weight = Math.abs(direction[0] + direction[1] + direction[2]) === 3 ? 30 : 20;
        for (let k = 0; k < weight; k += 1) rows.push(direction.map((value) => value * length / Math.sqrt(3)));
    });
    return rows;
};

// review:4
describe('exact tie resolution (lowest index)', () => {
    it('resolves <111>, <100> and <110> by the lowest-index rule at every frequency', () => {
        const directions = specialDirections();
        for (const frequency of [...Array.from({ length: 30 }, (_, i) => i + 1), 38]) {
            const tiling = goldbergTiling(frequency);
            expect(assignCells(tiling, directions), `nu=${frequency}`).toEqual(lowestIndexRule(tiling, directions));
        }
    });

    it('resolves Voronoi vertices (three-way ties) by the same rule', () => {
        for (const frequency of [3, 5, 10]) {
            const tiling = goldbergTiling(frequency);
            const vertices = tiling.polygons.flatMap((polygon) => polygon);
            expect(assignCells(tiling, vertices), `nu=${frequency}`).toEqual(lowestIndexRule(tiling, vertices));
        }
    });

    it('matches the Python pins for the special directions', () => {
        Object.entries(SPECIAL_DIRECTION_CELLS).forEach(([frequency, expected]) => {
            expect(assignCells(goldbergTiling(Number(frequency)), specialDirections()), `nu=${frequency}`).toEqual(expected);
        });
    });

    it('gives an 8-well <111> model the same peak as the Python engine', () => {
        Object.entries(WELL_MODEL_PEAKS).forEach(([frequency, expected]) => {
            const result = orientationHistogram(wellModel(), { frequency: Number(frequency), geometry: false });
            expect(result.peakCell, `nu=${frequency}`).toBe(expected.peakCell);
            result.peakDirection.forEach((value, i) => expect(Math.abs(value - expected.peakDirection[i])).toBeLessThan(1e-12));
            expect(result.counts.filter((count) => count > 0).length).toBe(expected.occupiedCells);
        });
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

// orientation.numerics.2/.14, orientation.parity.11/.16, orientation.physics.21
describe('tie-tolerant peak', () => {
    const tiedCloud = () => {
        const tiling = goldbergTiling(3);
        const points = [];
        for (const cell of [2, 10]) for (let i = 0; i < 20; i += 1) points.push(tiling.centers[cell]);
        for (const cell of [50, 60, 70]) points.push(tiling.centers[cell]);
        return points;
    };

    it('resolves a tied maximum to the lowest index and counts the ties', () => {
        const result = orientationHistogram(tiedCloud(), { frequency: 3, geometry: false });
        expect(result.peakCell).toBe(2);
        expect(result.peakTieCount).toBe(2);
        expect(result.peakDirection).toEqual(goldbergTiling(3).centers[2]);
    });

    it('reports a single cell for an untied maximum', () => {
        const gauss = makeRng(1);
        const points = [];
        for (let i = 0; i < 20000; i += 1) points.push([gauss() * 0.5, gauss() * 0.1, gauss() * 0.1]);
        const result = orientationHistogram(points, { frequency: 8, geometry: false });
        expect(result.peakTieCount).toBe(1);
        expect(result.enhancement[result.peakCell]).toBe(result.vmax);
    });
});

// --- significance statistics -----------------------------------------------

// scipy.special references: [a, x, gammainc(a, x), gammaincc(a, x)].
const SCIPY_GAMMA = [
    [1, 0.2, 0.18126924692201815, 0.8187307530779818],
    [2, 0.9029912627020413, 0.22861236738599036, 0.7713876326140097],
    [5, 0.15932045206675688, 7.492499892444724e-07, 0.9999992507500107],
    [39, 2.5649892741371745, 3.6267987030691646e-32, 1.0],
    [40, 0.2, 1.1087134276149818e-76, 1.0],
    [2.5, 3.0, 0.6937810815867218, 0.30621891841327825],
    [2.5, 40.0, 0.9999999999999991, 8.391825114831597e-16],
    [500.5, 480.0, 0.18030685630218252, 0.8196931436978174],
    [500.5, 560.0, 0.9949983047395444, 0.005001695260455566],
    [2880.5, 2600.0, 3.330569543243583e-08, 0.9999999666943046],
    [2880.5, 3300.0, 0.9999999999999606, 3.94213016105765e-14],
    [20480.5, 20000.0, 0.00036000517881559747, 0.9996399948211844],
    [20480.5, 22000.0, 1.0, 1.7291900418622563e-25],
    [0.5, 0.001, 0.035670591729679894, 0.9643294082703201],
    // Non-half-integer shapes, as the map test's gamma reference uses them.
    [126.07368146574035, 170.0, 0.999815486123718, 0.00018451387628206668],
    [208.39, 240.5, 0.9839109318641394, 0.016089068135860633],
    [2.4567, 1.66, 0.3607701545219544, 0.6392298454780456],
    [2.4567, 14.0, 0.9999666494786777, 3.335052132222686e-05],
    [1.2, 0.05, 0.024258165295933994, 0.9757418347040661],
    [19.72, 60.0, 0.9999999995402341, 4.5976590602670967e-10],
    [19.72, 8.0, 0.0003310527835421063, 0.9996689472164579]
];
// scipy.special.ndtri references: [p, ndtri(p)].
const SCIPY_NDTRI = [
    [1e-300, -37.0470962993612], [1e-100, -21.273453560965322], [1e-20, -9.262340089798409],
    [1e-05, -4.264890793922825], [0.02, -2.053748910631823], [0.3, -0.5244005127080409],
    [0.5, 0.0], [0.8, 0.8416212335729143], [0.999, 3.090232306167813],
    [0.999999999999, 7.0344869100478356]
];

const relClose = (actual, expected, rtol) => {
    expect(Math.abs(actual - expected)).toBeLessThanOrEqual(rtol * Math.max(1, Math.abs(expected)));
};

describe('special functions match scipy', () => {
    it('regularized incomplete gamma, small tail to 1e-11 relative', () => {
        SCIPY_GAMMA.forEach(([a, x, lower, upper]) => {
            const result = regularizedGamma(a, x);
            const [small, reference] = lower < upper ? [result.lower, lower] : [result.upper, upper];
            expect(Math.abs(small - reference) / reference).toBeLessThan(1e-11);
            relClose(result.lower + result.upper, 1, 1e-15);
        });
    });

    it('normal quantile to 1e-14', () => {
        SCIPY_NDTRI.forEach(([p, z]) => relClose(normalQuantile(p), z, 1e-14));
        expect(Math.abs(normalDeviate(0.5, 0.5))).toBe(0);
        relClose(normalDeviate(1e-5, 1 - 1e-5), 4.264890793922825, 1e-14);
    });

    it('logGamma to 1e-14', () => {
        relClose(logGamma(0.5), 0.5723649429247, 1e-14);
        relClose(logGamma(3.7), 1.428072326665388, 1e-14);
        relClose(logGamma(10), 12.801827480081469, 1e-14);
        relClose(logGamma(20480.5), 182830.05847756125, 1e-14);
    });
});

// Deterministic, RNG-free cloud — identical to golden_cloud() in
// tests/test_orientation_fixes.py.
const goldenCloud = ({ nSphere = 900, nLobe = 60, lobe = [0.3, -0.5, 0.81] } = {}) => {
    const points = [];
    for (let i = 0; i < nSphere; i += 1) {
        const z = 1 - (2 * i + 1) / nSphere;
        const r = Math.sqrt(1 - z * z);
        const phi = i * 2.399963229728653;
        const radius = 0.05 + 0.03 * ((i * 7) % 11) / 11;
        points.push([r * Math.cos(phi) * radius, r * Math.sin(phi) * radius, z * radius]);
    }
    const length = Math.hypot(...lobe);
    const unit = lobe.map((value) => value / length);
    for (let j = 0; j < nLobe; j += 1) {
        const jitter = [Math.cos(j * 1.3), Math.sin(j * 1.7), Math.cos(j * 0.9)].map((value) => 0.02 * value);
        points.push(unit.map((value, axis) => (value + jitter[axis]) * 0.12));
    }
    return points;
};

// Shared verbatim with GOLDEN_PEAK in tests/test_orientation_fixes.py:
// [frequency, smoothing, nLobe, expected fields].
const GOLDEN_PEAK = [
    [6, 1, 60, { peakCell: 181, peakCount: 39, peakExpected: 2.5649892741371745,
        peakLocalPValue: 3.6267987030691646e-32, peakPValue: 1.3129011305110377e-29,
        peakSignificance: 11.238918217145391 }],
    [10, 2, 60, { peakCell: 480, peakCount: 45, peakExpected: 0.8960392805524379,
        peakLocalPValue: 2.4905767221833515e-59, peakPValue: 2.4955578756277183e-56,
        peakSignificance: 15.770175177825124 }],
    [null, 0, 60, { peakCell: 10, peakCount: 72, peakExpected: 23.631937108780438,
        peakLocalPValue: 1.0240227941836545e-15, peakPValue: 4.300895735571259e-14,
        peakSignificance: 7.460770577845533 }],
    [10, 2, 0, { peakCell: 436, peakCount: 2, peakExpected: 0.9029912627020413,
        peakLocalPValue: 0.22861236738599036, peakPValue: 1.0,
        peakSignificance: -22.629329003077444 }],
    [2, 0, 6, { peakCell: 10, peakCount: 27, peakExpected: 22.30264064641154,
        peakLocalPValue: 0.18477705092278982, peakPValue: 0.999812237582693,
        peakSignificance: -3.5567124319438927 }]
];

const assertGolden = (table) => {
    table.forEach(([frequency, smoothing, nLobe, expected]) => {
        const result = orientationHistogram(goldenCloud({ nLobe }), { frequency, smoothing, geometry: false });
        Object.entries(expected).forEach(([key, value]) => {
            if (typeof value === 'boolean' || Number.isInteger(value) || value === null) {
                expect(result[key], `${frequency}/${smoothing}/${nLobe} ${key}`).toBe(value);
            } else {
                // Relative tolerance: p-values reach 1e-142.
                expect(Math.abs(result[key] - value), `${frequency}/${smoothing}/${nLobe} ${key}`)
                    .toBeLessThanOrEqual(1e-9 * Math.abs(value));
            }
        });
    });
};

const isotropicUnits = (gauss, n) => {
    const points = [];
    for (let i = 0; i < n; i += 1) points.push([gauss(), gauss(), gauss()]);
    return points;
};

// orientation.numerics.1/.26, orientation.physics.20
describe('peak significance (look-elsewhere-corrected Poisson tail)', () => {
    it('matches the Python golden values', () => assertGolden(GOLDEN_PEAK));

    it('reads isotropic clouds as not significant at the UI defaults', () => {
        const gauss = makeRng(2024);
        const significance = [];
        const local = [];
        for (let k = 0; k < 80; k += 1) {
            const result = orientationHistogram(isotropicUnits(gauss, 216), { frequency: 10, smoothing: 2, geometry: false });
            significance.push(result.peakSignificance);
            local.push(result.peakZScore);
        }
        // review:3 -- a missing or NaN field makes every `> 2` count below
        // false, so the rate checks alone pass vacuously (they did on main).
        expect(significance.every(Number.isFinite)).toBe(true);
        expect(local.every(Number.isFinite)).toBe(true);
        expect(local.filter((value) => value >= 3).length / local.length).toBeGreaterThan(0.5);
        expect(significance.filter((value) => value > 2).length / significance.length).toBeLessThanOrEqual(0.05);
        expect(significance.filter((value) => value > 3).length).toBeLessThanOrEqual(1);
    });

    it('reads a real lobe as significant (positive control)', () => {
        const result = orientationHistogram(goldenCloud(), { frequency: 6, geometry: false });
        expect(result.peakSignificance).toBeGreaterThan(5);
        expect(result.peakPValue).toBeLessThan(1e-6);
    });
});

// orientation.numerics.3, orientation.physics.7 — shared verbatim with
// GOLDEN_MAP in tests/test_orientation_fixes.py.
const GOLDEN_MAP = [
    [6, 1, 60, { mapChiSquare: 800.1657022469657, mapDegreesOfFreedom: 361,
        mapNullSd: 26.89394588964023, mapNullSkewness: 0.178121379052147,
        mapExpectedPairs: 1288.2313214115104,
        mapPValue: 8.22254884213926e-33, mapSignificance: 11.87251368423015 }],
    [10, 2, 60, { mapChiSquare: 2596.980147082603, mapDegreesOfFreedom: 1001,
        mapNullSd: 44.93009360956353, mapNullSkewness: 0.13854679111724788,
        mapExpectedPairs: 466.6502011243447,
        mapPValue: 5.115170092676409e-113, mapSignificance: 22.56203637660662 }],
    [null, 0, 60, { mapChiSquare: 109.05527451398147, mapDegreesOfFreedom: 41,
        mapNullSd: 9.050975709801982, mapNullSkewness: 0.4503623992125751,
        mapExpectedPairs: 10991.483242707845,
        mapPValue: 4.869112256098684e-08, mapSignificance: 5.331542086891653 }],
    // A Fibonacci sphere is far more uniform than counting noise: X^2 lies
    // below the reference's support, so the deviate sits at the floor.
    [10, 2, 0, { mapChiSquare: 206.08354289942974, mapDegreesOfFreedom: 1001,
        mapNullSd: 44.942491392316214, mapNullSkewness: 0.14177901437232995,
        mapExpectedPairs: 410.1132665642458,
        mapPValue: 1.0, mapSignificance: -37.0470962993612 }],
    [2, 0, 6, { mapChiSquare: 2.8616952082936087, mapDegreesOfFreedom: 41,
        mapNullSd: 9.05071282835656, mapNullSkewness: 0.450877377382182,
        mapExpectedPairs: 9789.110678651203,
        mapPValue: 1.0, mapSignificance: -8.886765298073716 }]
];

// Deterministic unequal cell probabilities, built identically in Python
// (_pinned_probabilities in tests/test_orientation_fixes.py).
const pinnedProbabilities = (cells) => {
    const raw = [];
    for (let i = 1; i <= cells; i += 1) raw.push(1 + 0.3 * Math.sin(i));
    const total = raw.reduce((sum, value) => sum + value, 0);
    return raw.map((value) => value / total);
};

// Shared verbatim with PEARSON_MOMENTS in tests/test_orientation_fixes.py:
// [N, C, exact mean, variance, third central moment of X^2] for
// pinnedProbabilities(C), from full multinomial enumeration.
const PEARSON_MOMENTS = [
    [3, 4, 3.0, 4.218531116561122, 13.943752126285027],
    [6, 5, 4.0, 6.91373687143512, 29.01155630009905],
    [9, 4, 3.0, 5.406177038853709, 20.39545791721937],
    [12, 3, 2.0, 3.6728421787498275, 12.957998592922594],
    [2, 7, 6.0, 7.183146093091783, 43.811413401314184]
];

// Exact moments of X^2 by enumerating every multinomial outcome.
const exactPearsonMoments = (used, p) => {
    const outcomes = [];
    const walk = (cell, left, counts) => {
        if (cell === p.length - 1) {
            outcomes.push([...counts, left]);
            return;
        }
        for (let k = 0; k <= left; k += 1) walk(cell + 1, left - k, [...counts, k]);
    };
    walk(0, used, []);
    const factorial = (n) => { let value = 1; for (let i = 2; i <= n; i += 1) value *= i; return value; };
    let mean = 0;
    const rows = outcomes.map((counts) => {
        let weight = factorial(used);
        let x2 = 0;
        counts.forEach((c, m) => {
            weight *= p[m] ** c / factorial(c);
            x2 += (c - used * p[m]) ** 2 / (used * p[m]);
        });
        mean += weight * x2;
        return [weight, x2];
    });
    let variance = 0;
    let third = 0;
    rows.forEach(([weight, x2]) => {
        variance += weight * (x2 - mean) ** 2;
        third += weight * (x2 - mean) ** 3;
    });
    return [mean, variance, third];
};

describe('map significance (Pearson X^2, exact-moment reference)', () => {
    it('matches the Python golden values', () => assertGolden(GOLDEN_MAP));

    it('has the exact multinomial null moments', () => {
        PEARSON_MOMENTS.forEach(([used, cells, ...expected]) => {
            const p = pinnedProbabilities(cells);
            const moments = pearsonNullMoments(used, p);
            const exact = exactPearsonMoments(used, p);
            expected.forEach((value, i) => {
                relClose(moments[i], value, 1e-12);
                relClose(exact[i], value, 1e-12);
            });
        });
    });

    it('reads isotropic clouds as noise and a one-sided cloud as overwhelming', () => {
        const gauss = makeRng(77);
        const values = [];
        for (let k = 0; k < 80; k += 1) {
            values.push(orientationHistogram(isotropicUnits(gauss, 216), { frequency: 10, smoothing: 2, geometry: false }).mapSignificance);
        }
        expect(values.every(Number.isFinite)).toBe(true);
        expect(values.filter((value) => value > 2).length / values.length).toBeLessThanOrEqual(0.06);
        expect(values.filter((value) => value > 3).length).toBeLessThanOrEqual(1);
        const hemisphere = isotropicUnits(gauss, 1000).map(([x, y, z]) => [Math.abs(x), y, z]);
        const oneSided = orientationHistogram(hemisphere, { frequency: 10, smoothing: 2, geometry: false });
        expect(oneSided.significance).toBeLessThan(1.5);
        expect(oneSided.mapSignificance).toBeGreaterThan(10);
    });

    it('stays calibrated on sparse isotropic maps (27 copies at nu=10)', () => {
        // chi^2_{C-1} put ~2% of these above 3 sigma (15x nominal).
        const gauss = makeRng(2718);
        const values = [];
        for (let k = 0; k < 600; k += 1) {
            values.push(orientationHistogram(isotropicUnits(gauss, 27), { frequency: 10, smoothing: 2, geometry: false }).mapSignificance);
        }
        expect(values.every(Number.isFinite)).toBe(true);
        expect(values.filter((value) => value > 3).length).toBeLessThanOrEqual(4);
        expect(values.filter((value) => value > 2).length / values.length).toBeLessThanOrEqual(0.04);
    });

    it('detects a tight lobe even when the map is sparse', () => {
        const gauss = makeRng(5);
        const cap = [];
        while (cap.length < 27) {
            const [x, y, z] = isotropicUnits(gauss, 1)[0];
            const length = Math.hypot(x, y, z);
            if (z / length > Math.cos(10 * Math.PI / 180)) cap.push([x, y, z]);
        }
        expect(orientationHistogram(cap, { frequency: 10, geometry: false }).mapSignificance).toBeGreaterThan(8);
    });

    it('withholds the test when fewer than MAP_TEST_MIN_PAIRS coincidences are expected', () => {
        const gauss = makeRng(4);
        [1, 2, 8].forEach((n) => {
            const result = orientationHistogram(isotropicUnits(gauss, n), { frequency: 10, geometry: false });
            expect(result.mapExpectedPairs).toBeLessThan(MAP_TEST_MIN_PAIRS);
            expect(result.mapPValue).toBeNull();
            expect(result.mapSignificance).toBeNull();
            expect(Number.isFinite(result.mapChiSquare)).toBe(true);
        });
        const populated = orientationHistogram(isotropicUnits(gauss, 27), { frequency: 10, geometry: false });
        expect(populated.mapExpectedPairs).toBeGreaterThanOrEqual(MAP_TEST_MIN_PAIRS);
        expect(Number.isFinite(populated.mapSignificance)).toBe(true);
        // The helper is the engine's path.
        const p = pinnedProbabilities(12);
        expect(mapTest(30, 12, p).deviate).toBeGreaterThan(0);
    });
});

// orientation.numerics.4/.13, orientation.physics.6 — shared verbatim with
// GOLDEN_ASYMMETRY in tests/test_orientation_fixes.py.
const GOLDEN_ASYMMETRY = [
    [6, 1, 60, { antipodalAsymmetry: 0.2, antipodalAsymmetryNull: 0.3401168600567151,
        antipodalAsymmetryNullSd: 0.019345746381926036,
        antipodalAsymmetryZ: -7.2427735425923245, antipodalAsymmetrySignificant: false }],
    [null, 0, 60, { antipodalAsymmetry: 0.08333333333333333,
        antipodalAsymmetryNull: 0.11738132743663726,
        antipodalAsymmetryNullSd: 0.01946344130131603,
        antipodalAsymmetryZ: -1.7493306335813166,
        antipodalAsymmetrySignificant: false }],
    [10, 2, 0, { antipodalAsymmetry: 0.19111111111111112,
        antipodalAsymmetryNull: 0.5672222222222222,
        antipodalAsymmetryNullSd: 0.020957040126511128,
        antipodalAsymmetryZ: -17.94676675907692, antipodalAsymmetrySignificant: false }],
    [4, 0, 400, { antipodalAsymmetry: 0.3415384615384615,
        antipodalAsymmetryNull: 0.1752318529176018,
        antipodalAsymmetryNullSd: 0.01673666421475673,
        antipodalAsymmetryZ: 9.936663990320547, antipodalAsymmetrySignificant: true }]
];

describe('antipodal asymmetry null (conditional binomial split)', () => {
    it('matches the Python golden values', () => assertGolden(GOLDEN_ASYMMETRY));

    it('stays below the statistic bound and rarely flags a centrosymmetric site', () => {
        const gauss = makeRng(11);
        const small = orientationHistogram(isotropicUnits(gauss, 216), { frequency: 10, smoothing: 2, geometry: false });
        expect(small.antipodalAsymmetryNull).toBeLessThan(0.9);
        let flagged = 0;
        const zValues = [];
        for (let k = 0; k < 60; k += 1) {
            const rod = isotropicUnits(gauss, 1000).map(([x, y, z]) => [x * 0.2, y * 0.03, z * 0.03]);
            const result = orientationHistogram(rod, { frequency: 10, geometry: false });
            expect(typeof result.antipodalAsymmetrySignificant).toBe('boolean');
            if (result.antipodalAsymmetrySignificant) flagged += 1;
            zValues.push(result.antipodalAsymmetryZ);
        }
        expect(zValues.every(Number.isFinite)).toBe(true);
        expect(flagged).toBeLessThanOrEqual(1);
        expect(Math.abs(zValues.reduce((sum, value) => sum + value, 0) / zValues.length)).toBeLessThan(0.4);
    });

    it('flags a one-sided cloud at the UI default resolution', () => {
        const gauss = makeRng(12);
        const hemisphere = isotropicUnits(gauss, 1000).map(([x, y, z]) => [Math.abs(x), y, z]);
        const result = orientationHistogram(hemisphere, { frequency: 10, smoothing: 2, geometry: false });
        expect(result.antipodalAsymmetrySignificant).toBe(true);
        expect(result.antipodalAsymmetryZ).toBeGreaterThan(10);
    });

    it('reports no z when every occupied pair holds one atom', () => {
        const tiling = goldbergTiling(4);
        const cells = tiling.centers.map((center, cell) => cell).filter((cell) => cell < tiling.antipode[cell]).slice(0, 20);
        const result = orientationHistogram(cells.map((cell) => tiling.centers[cell].map((value) => value * 0.1)), { frequency: 4, geometry: false });
        expect(result.antipodalAsymmetryNullSd).toBe(0);
        expect(result.antipodalAsymmetryZ).toBeNull();
        expect(result.antipodalAsymmetrySignificant).toBe(false);
        expect(result.antipodalAsymmetryNull).toBe(1);
    });
});

// orientation.physics.9/.24, orientation.numerics.27 — shared verbatim with
// GOLDEN_ANISOTROPY in tests/test_orientation_fixes.py.
const GOLDEN_ANISOTROPY = [
    [6, 1, 60, { orientationAnisotropy: 0.12495614840544977, orientationEffectivePoints: 960.0,
        orientationAnisotropyNull: 0.05182412242070032,
        orientationBinghamStatistic: 18.73684781045718,
        orientationBinghamPValue: 0.0021515346560243643,
        orientationAnisotropySignificance: 2.855045282725758 }],
    [10, 2, 0, { orientationEffectivePoints: 900.0, orientationAnisotropyNull: 0.05352372348458313,
        orientationBinghamPValue: 0.9999999999999987,
        orientationAnisotropySignificance: -7.9084691847240896 }],
    [4, 0, 400, { orientationAnisotropy: 0.6150392816231949, orientationEffectivePoints: 1300.0,
        orientationAnisotropyNull: 0.04453442987940475,
        orientationBinghamStatistic: 614.6941425286471,
        orientationBinghamPValue: 1.3514132534260834e-130,
        orientationAnisotropySignificance: 24.286800999672977 }],
    [2, 0, 6, { orientationAnisotropy: 0.013271734268148538, orientationEffectivePoints: 906.0,
        orientationAnisotropyNull: 0.05334619820786263,
        orientationBinghamStatistic: 0.19947859577278887,
        orientationBinghamPValue: 0.9991194624727917,
        orientationAnisotropySignificance: -3.127820820781734 }]
];
const GOLDEN_ANISOTROPY_AMPLITUDE2 = {
    orientationEffectivePoints: 725.3165446406116,
    orientationAnisotropyNull: 0.05962162122187657,
    orientationBinghamStatistic: 129.85364394492348,
    orientationBinghamPValue: 2.5564292957448824e-26,
    orientationAnisotropySignificance: 10.549388830489368
};

describe('orientation anisotropy null (Bingham test)', () => {
    it('matches the Python golden values', () => {
        assertGolden(GOLDEN_ANISOTROPY);
        const weighted = orientationHistogram(goldenCloud(), { frequency: 6, smoothing: 1, weight: 'amplitude2', geometry: false });
        Object.entries(GOLDEN_ANISOTROPY_AMPLITUDE2).forEach(([key, value]) => {
            expect(Math.abs(weighted[key] - value), key).toBeLessThanOrEqual(1e-9 * Math.abs(value));
        });
    });

    it('reports the isotropic expectation and stays calibrated on noise', () => {
        const gauss = makeRng(21);
        const anisotropy = [];
        let above2 = 0;
        let expected = 0;
        for (let k = 0; k < 120; k += 1) {
            const result = orientationHistogram(isotropicUnits(gauss, 216), { frequency: 2, geometry: false });
            anisotropy.push(result.orientationAnisotropy);
            expected = result.orientationAnisotropyNull;
            expect(Number.isFinite(result.orientationAnisotropySignificance)).toBe(true);
            if (result.orientationAnisotropySignificance > 2) above2 += 1;
        }
        expect(expected).toBeCloseTo(9 / Math.sqrt(10 * Math.PI * 216), 12);
        const mean = anisotropy.reduce((sum, value) => sum + value, 0) / anisotropy.length;
        expect(Math.abs(mean / expected - 1)).toBeLessThan(0.12);
        expect(above2 / anisotropy.length).toBeLessThanOrEqual(0.06);
    });
});
