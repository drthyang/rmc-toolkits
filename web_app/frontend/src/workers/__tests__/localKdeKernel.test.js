// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The browser KDE kernel is scipy's gaussian_kde kernel H = f^2 C evaluated
// exactly (no ridge, no isotropic fallback), and the worker declines exactly
// the slabs kde.py declines. Python-vs-JS numbers live in kdeParity.test.js;
// these checks pin the JS engine to an in-test brute force.

import { describe, expect, it } from 'vitest';
import {
    KDE_MESSAGES,
    computeDensityCpu,
    computeKde,
    covariance,
    hasTwoDimensionalSpread,
    makeKernel
} from '../localKdeWorker';

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

// Brute-force bivariate normal mixture with H = f^2 C, via the explicit inverse.
// Terms below exp(-60) are dropped, as in the worker (and the WGSL shader):
// SciPy keeps them, but they are < 1e-26 of a kernel's peak and matter only
// when every node sits more than 11 sigma from every atom (an empty map).
const bruteForce = (samples, factor, { grid, xMin, yMin, xStep, yStep }) => {
    const { c00, c01, c11 } = covariance(samples);
    const h00 = factor * factor * c00;
    const h01 = factor * factor * c01;
    const h11 = factor * factor * c11;
    const det = h00 * h11 - h01 * h01;
    const out = [];
    for (let y = 0; y < grid; y += 1) {
        const row = [];
        for (let x = 0; x < grid; x += 1) {
            const gx = xMin + x * xStep;
            const gy = yMin + y * yStep;
            let sum = 0;
            samples.forEach(([u, v]) => {
                const dx = gx - u;
                const dy = gy - v;
                const q = (h11 * dx * dx - 2 * h01 * dx * dy + h00 * dy * dy) / det;
                if (-0.5 * q > -60) sum += Math.exp(-0.5 * q);
            });
            row.push(sum / (2 * Math.PI * Math.sqrt(det) * samples.length));
        }
        out.push(row);
    }
    return out;
};

const cluster = (random, centre, sigma, count) => Array.from({ length: count }, () => [
    centre[0] + sigma * gaussian(random),
    centre[1] + sigma * gaussian(random)
]);

describe('makeKernel is exactly H = f^2 C', () => {
    const cases = {
        // One compact site: det(H) ~ 1e-15, where the old kernel swapped in an
        // inflated isotropic fallback.
        'single site': (random) => cluster(random, [0.5, 0.5], 0.006, 400),
        // Two sites on the anti-diagonal: a needle kernel whose minor variance
        // (~3e-9 at f = 0.01) the old fixed 1e-8 ridge swamped.
        'two-site needle': (random) => [
            ...cluster(random, [0.25, 0.75], 0.006, 200),
            ...cluster(random, [0.75, 0.25], 0.006, 200)
        ],
        'cell-filling': (random) => Array.from({ length: 400 }, () => [random(), random()])
    };
    for (const [name, build] of Object.entries(cases)) {
        for (const factor of [0.005, 0.01, 0.03]) {
            it(`${name}, bw=${factor}`, () => {
                const samples = build(makeRandom(7));
                const kernel = makeKernel(samples, factor);
                // A 41 x 41 window of +-6 major kernel sigmas around the first
                // atom, so even the narrowest kernel is resolved by the nodes.
                const half = 6 * kernel.sigmaMajor;
                const window = {
                    grid: 41,
                    xMin: samples[0][0] - half,
                    yMin: samples[0][1] - half,
                    xStep: (2 * half) / 40,
                    yStep: (2 * half) / 40
                };
                const density = computeDensityCpu({ samples, kernel, ...window });
                const reference = bruteForce(samples, factor, window);
                const peak = Math.max(...reference.flat());
                let worst = 0;
                reference.forEach((row, y) => row.forEach((value, x) => {
                    worst = Math.max(worst, Math.abs(density[y][x] - value) / peak);
                }));
                expect(worst).toBeLessThan(1e-9);
                const { c00, c01, c11 } = covariance(samples);
                const f2 = factor * factor;
                expect(kernel.covariance).toEqual([[c00 * f2, c01 * f2], [c01 * f2, c11 * f2]]);
            });
        }
    }
});

describe('degenerate slabs are declined, not regularized', () => {
    it('rank test follows numpy: exact lines are rank 1, tiny spread is rank 2', () => {
        const line = Array.from({ length: 33 }, (_, k) => [(16 + k) / 64, 0.25 + 0.5 * (16 + k) / 64]);
        expect(hasTwoDimensionalSpread(line)).toBe(false);
        const random = makeRandom(3);
        const near = line.map(([u, v]) => [u, v + 1e-9 * gaussian(random)]);
        expect(hasTwoDimensionalSpread(near)).toBe(true);
        expect(hasTwoDimensionalSpread([[0.1, 0.2], [0.1, 0.2], [0.1, 0.2]])).toBe(false);
    });

    const plane = (uv) => uv.map(([x, y]) => ({ x, y, z: 0.5 }));
    const slice = (points, overrides = {}) => computeKde({
        points,
        normal: [0, 0, 1],
        uVector: [1, 0, 0],
        vVector: [0, 1, 0],
        range: [0, 1],
        zCenter: 0.5,
        thickness: 0.1,
        bandwidth: 0.03,
        gridSize: 16,
        logScale: false,
        ...overrides
    });

    it('declines a collinear slab instead of drawing a ridge-regularized kernel', async () => {
        const line = Array.from({ length: 33 }, (_, k) => [(16 + k) / 64, 0.25 + 0.5 * (16 + k) / 64]);
        const result = await slice(plane(line));
        expect(result.slabCount).toBe(33);
        expect(result.fitCount).toBe(0);
        expect(result.vmax).toBe(0);
        expect(result.kernel).toBeNull();
        expect(result.message).toBe(KDE_MESSAGES.collinear);
    });

    it('declines fewer than three distinct positions', async () => {
        const pair = [...Array(10).fill([0.3, 0.4]), ...Array(10).fill([0.6, 0.7])];
        const result = await slice(plane(pair));
        expect(result.fitCount).toBe(0);
        expect(result.message).toBe(KDE_MESSAGES.fewUnique);
    });

    it('uses the bandwidth as given: no substitution for 0 and no 1e-4 floor', async () => {
        const random = makeRandom(9);
        const points = plane(Array.from({ length: 50 }, () => [random(), random()]));
        for (const bandwidth of [0, -0.03, Number.NaN, Number.POSITIVE_INFINITY]) {
            const result = await slice(points, { bandwidth });
            expect(result.fitCount).toBe(0);
            expect(result.bw).toBeNull();
            expect(result.message).toBe(KDE_MESSAGES.bandwidth);
        }
        // Below the old 1e-4 floor the factor still scales H as f^2.
        const tiny = await slice(points, { bandwidth: 5e-5 });
        const base = await slice(points, { bandwidth: 0.03 });
        expect(tiny.kernel.covariance[0][0] / base.kernel.covariance[0][0]).toBeCloseTo((5e-5 / 0.03) ** 2, 15);
    });
});
