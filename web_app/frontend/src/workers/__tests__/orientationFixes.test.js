// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Regression tests for the 0.6.0 audit of the displacement-direction engine —
// the JS twin of the rules in tests/test_orientation_fixes.py. The values the
// two engines must share (recommendedFrequency boundaries, tie-family cells,
// neighbour checksums, polygon starts, the 8-well peaks and every significance
// statistic on the golden cloud) are not copied here: they come from one
// generated golden, orientation_fixture.json, checked in
// src/__tests__/orientationParity.test.js.

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
import fixture from '../../__tests__/fixtures/orientation_fixture.json';

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

// orientation.physics.10/.25, orientation.numerics.30 (the boundary table
// both engines share is in orientation_fixture.json).
describe('recommendedFrequency floors to the target occupancy', () => {
    it('never drops below 12 points per cell (except at the dodecahedron floor)', () => {
        for (let n = 1; n < 20000; n += 7) {
            const frequency = recommendedFrequency(n);
            if (frequency > MIN_FREQUENCY) expect(n / (10 * frequency * frequency + 2)).toBeGreaterThanOrEqual(12);
            if (frequency < 24) expect(n / (10 * (frequency + 1) ** 2 + 2)).toBeLessThan(12);
        }
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
});

// <111>, <100> and <110> (26 directions).
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

// review:4 (the cells both engines give the 26 directions, and the 8-well
// <111> model's peaks, are in orientation_fixture.json).
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
});

// orientation.parity.17 — cyclic neighbour / polygon order must not depend on
// the atan2 branch cut (the neighbour checksums and nu = 38 polygon starts
// both engines give are in orientation_fixture.json).
describe('canonical cyclic order', () => {
    it('starts every neighbour row at its smallest index', () => {
        for (const frequency of [3, 6, 10]) {
            goldbergTiling(frequency).neighbors.forEach((row) => {
                const valid = row.filter((value) => value >= 0);
                expect(valid[0]).toBe(Math.min(...valid));
            });
        }
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
    // The worker transport's 'all': the request contract
    // (fixtures/analysis_requests.json, src/__tests__/analysisRequests.test.js).
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

// The golden cloud (tests/generate_orientation_fixture.py, stored in the
// fixture): an RNG-free Fibonacci sphere plus a tight lobe of 60 points.
const goldenCloud = () => [...fixture.inputs['golden-sphere'], ...fixture.inputs['golden-lobe'].slice(0, 60)];

const isotropicUnits = (gauss, n) => {
    const points = [];
    for (let i = 0; i < n; i += 1) points.push([gauss(), gauss(), gauss()]);
    return points;
};

// orientation.numerics.1/.26, orientation.physics.20
describe('peak significance (look-elsewhere-corrected Poisson tail)', () => {
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

// orientation.numerics.3, orientation.physics.7

// Deterministic unequal cell probabilities.
const unequalProbabilities = (cells) => {
    const raw = [];
    for (let i = 1; i <= cells; i += 1) raw.push(1 + 0.3 * Math.sin(i));
    const total = raw.reduce((sum, value) => sum + value, 0);
    return raw.map((value) => value / total);
};

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
    it('has the exact multinomial null moments', () => {
        // The golden's cases (whose moments Python's engine gives, checked in
        // orientationParity.test.js), against full enumeration here.
        fixture.pearsonMoments.forEach(({ used, probabilities }) => {
            const moments = pearsonNullMoments(used, probabilities);
            exactPearsonMoments(used, probabilities).forEach((value, i) => relClose(moments[i], value, 1e-12));
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
        const p = unequalProbabilities(12);
        expect(mapTest(30, 12, p).deviate).toBeGreaterThan(0);
    });
});

// orientation.numerics.4/.13, orientation.physics.6
describe('antipodal asymmetry null (conditional binomial split)', () => {
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

// orientation.physics.9/.24, orientation.numerics.27
describe('orientation anisotropy null (Bingham test)', () => {
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

// Plan item 2.5: the direction resolution, a measurement limit (the twin of
// PeakCellRadiusTests). peakDirection is the peak cell's centre, so it names a
// direction no more finely than the cell's circumscribed radius,
// peakCellRadiusDeg: every direction assigned to the cell lies within it.
describe('the peak cell radius (peakCellRadiusDeg)', () => {
    const angle = (u, v) => Math.acos(Math.max(-1, Math.min(1, (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / Math.hypot(...u))));

    it('is the dodecahedron\'s closed form at nu = 1', () => {
        const exact = Math.acos(Math.sqrt((5 + 2 * Math.sqrt(5)) / 15));
        goldbergTiling(1).radii.forEach((radius) => expect(Math.abs(radius - exact)).toBeLessThan(1e-12));
        expect(exact * 180 / Math.PI).toBeCloseTo(37.3774, 4);
    });

    it('bounds every direction assigned to a cell', () => {
        const directions = cloud(6000, 31);
        for (const nu of [1, 2, 3, 6]) {
            const tiling = goldbergTiling(nu);
            const cells = assignCells(tiling, directions);
            cells.forEach((cell, index) => {
                expect(angle(directions[index], tiling.centers[cell])).toBeLessThanOrEqual(tiling.radii[cell] + 1e-12);
            });
        }
    });

    it('puts the 15 degree line between nu = 2 and nu = 3', () => {
        const degrees = (nu) => goldbergTiling(nu).radii.map((radius) => radius * 180 / Math.PI);
        expect(Math.min(...degrees(2))).toBeGreaterThan(15);
        expect(Math.max(...degrees(3))).toBeLessThan(15);
    });

    it('is reported for the peak cell', () => {
        for (const frequency of [1, 2, 4, null]) {
            const result = orientationHistogram(goldenCloud(), { frequency, geometry: false });
            const tiling = goldbergTiling(result.frequency);
            expect(result.peakCellRadiusDeg).toBe(tiling.radii[result.peakCell] * (180 / Math.PI));
        }
    });
});
