// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The worker's kernel comes from the slab's source atoms (one row per atom),
// not from the periodic images or the fit subsample, exactly as in kde.py
// (tests/test_kde_bandwidth_source.py is the Python twin).

import { describe, expect, it } from 'vitest';
import { augmentPeriodicImages, computeKde, covariance, imageRank } from '../localKdeWorker';

const makeRandom = (seed) => {
    let value = seed >>> 0;
    return () => {
        value += 0x6D2B79F5;
        let mixed = value;
        mixed = Math.imul(mixed ^ (mixed >>> 15), mixed | 1);
        mixed ^= mixed + Math.imul(mixed ^ (mixed >>> 7), mixed | 61);
        return ((mixed ^ (mixed >>> 14)) >>> 0) / 4294967296;
    };
};
const gaussian = (random) => Math.sqrt(-2 * Math.log(1 - random())) * Math.cos(2 * Math.PI * random());
const fold = (value) => ((value % 1) + 1) % 1;

// Two isotropic sites on the anti-diagonal at z = 0.25 (like Ga in GaNb4Se8).
const twoSiteLayer = () => {
    const random = makeRandom(11);
    const points = [];
    [[0.25, 0.75], [0.75, 0.25]].forEach(([u, v]) => {
        for (let i = 0; i < 400; i += 1) {
            points.push({
                x: fold(u + 0.006 * gaussian(random)),
                y: fold(v + 0.006 * gaussian(random)),
                z: fold(0.25 + 0.006 * gaussian(random))
            });
        }
    });
    return points;
};

const cSlice = (points, overrides = {}) => computeKde({
    points,
    normal: [0, 0, 1],
    uVector: [1, 0, 0],
    vVector: [0, 1, 0],
    range: [0, 1],
    zCenter: 0.25,
    thickness: 0.08,
    bandwidth: 0.03,
    gridSize: 32,
    logScale: false,
    ...overrides
});

const expectSameMatrix = (actual, expected, tolerance = 1e-12) => {
    const scale = Math.max(...expected.flat().map(Math.abs));
    actual.flat().forEach((value, index) => {
        expect(Math.abs(value - expected.flat()[index]) / scale).toBeLessThan(tolerance);
    });
};

const scaled = ({ c00, c01, c11 }, factor) => [
    [c00 * factor * factor, c01 * factor * factor],
    [c01 * factor * factor, c11 * factor * factor]
];

describe('kernel from the source atoms', () => {
    it('ranks the unshifted atom first and every image uniquely', () => {
        const ranks = new Set();
        for (let ox = -1; ox <= 1; ox += 1) {
            for (let oy = -1; oy <= 1; oy += 1) {
                for (let oz = -1; oz <= 1; oz += 1) ranks.add(imageRank(ox, oy, oz));
            }
        }
        expect(ranks.size).toBe(27);
        expect(Math.min(...ranks)).toBe(imageRank(0, 0, 0));
        const { imageRanks } = augmentPeriodicImages([{ x: 0.02, y: 0.5, z: 0.5 }], 0.1);
        expect(imageRanks).toEqual([imageRank(0, 0, 0), imageRank(1, 0, 0)]);
    });

    it('is bw^2 times the covariance of the folded source atoms', async () => {
        const points = twoSiteLayer();
        const result = await cSlice(points);
        expect(result.slabCount).toBe(points.length);
        expectSameMatrix(result.kernel.covariance, scaled(covariance(points.map((p) => [p.x, p.y])), 0.03));
    });

    it('does not change when the thickness moves the image margin', async () => {
        // dz 0.08 -> margin 0.1; dz 0.30 -> margin 0.3 admits the x/y images at 1.25 / -0.25.
        const points = twoSiteLayer();
        const thin = await cSlice(points, { thickness: 0.08 });
        const thick = await cSlice(points, { thickness: 0.3 });
        expect(thick.slabCount).toBe(thin.slabCount);
        expect(thick.fitCount).toBeGreaterThan(thin.fitCount);
        expectSameMatrix(thick.kernel.covariance, thin.kernel.covariance);
        expect(thick.vmax / thin.vmax).toBeCloseTo(1, 6);
    });

    it('scales exactly as bw^2 across the bandwidth margin step', async () => {
        const points = twoSiteLayer();
        const below = await cSlice(points, { bandwidth: 0.12 });
        const above = await cSlice(points, { bandwidth: 0.13 });
        below.kernel.covariance.flat().forEach((value, index) => {
            expect(above.kernel.covariance.flat()[index] / value).toBeCloseTo((0.13 / 0.12) ** 2, 12);
        });
    });

    it('represents depth-wrapped atoms by their nearest image', async () => {
        const random = makeRandom(2);
        const points = Array.from({ length: 60 }, () => ({ x: random(), y: random(), z: 0.97 + 0.02 * random() }));
        const result = await cSlice(points, { zCenter: 0, thickness: 0.1 });
        expect(result.slabCount).toBe(60);
        expectSameMatrix(result.kernel.covariance, scaled(covariance(points.map((p) => [p.x, p.y])), 0.03));
    });

    it('does not depend on the 6000-point fit subsample', async () => {
        const random = makeRandom(3);
        const points = Array.from({ length: 8000 }, () => ({ x: random(), y: random(), z: 0.5 }));
        const result = await cSlice(points, { zCenter: 0.5, gridSize: 16 });
        expect(result.fitCount).toBe(6000);
        expectSameMatrix(result.kernel.covariance, scaled(covariance(points.map((p) => [p.x, p.y])), 0.03));
    });
});
