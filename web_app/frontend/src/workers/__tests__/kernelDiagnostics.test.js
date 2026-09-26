// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Kernel diagnostics: the worker flags a kernel narrower than half a grid step
// (an aliased map) exactly as kde.py does (tests/test_kde_kernel_diagnostics.py),
// and the page's Angstrom readout maps H through the cell metric.

import { describe, expect, it } from 'vitest';
import { KDE_WARNINGS, computeKde } from '../localKdeWorker';
import { kernelSigmaAngstrom } from '../slabSelection';

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

const cSlice = (points, overrides = {}) => computeKde({
    points,
    normal: [0, 0, 1],
    uVector: [1, 0, 0],
    vVector: [0, 1, 0],
    range: [0, 1],
    zCenter: 0.5,
    thickness: 0.08,
    bandwidth: 0.03,
    gridSize: 120,
    logScale: false,
    ...overrides
});

describe('sub-grid kernel warning', () => {
    it('flags the needle kernel of a two-site layer', async () => {
        const random = makeRandom(5);
        const points = [];
        [[0.25, 0.75], [0.75, 0.25]].forEach(([u, v]) => {
            for (let i = 0; i < 300; i += 1) {
                points.push({ x: u + 0.006 * gaussian(random), y: v + 0.006 * gaussian(random), z: 0.5 });
            }
        });
        const result = await cSlice(points);
        expect(result.kernel.sigmaMinor).toBeLessThan(0.5 / 119);
        expect(result.warnings).toEqual([{ code: 'subgrid', message: KDE_WARNINGS.subgrid }]);
    });

    it('stays quiet for a cell-filling slab at grid 120 but not at grid 16', async () => {
        const random = makeRandom(6);
        const points = Array.from({ length: 3000 }, () => ({ x: random(), y: random(), z: 0.5 }));
        expect((await cSlice(points)).warnings).toEqual([]);
        expect((await cSlice(points, { gridSize: 16 })).warnings.map((warning) => warning.code)).toEqual(['subgrid']);
    });
});

describe('kernel sigma in Angstrom', () => {
    it('scales by the cell edge in a cubic cell', () => {
        const sigma = kernelSigmaAngstrom([[4e-4, 0], [0, 1e-4]], [10, 0, 0], [0, 10, 0]);
        expect(sigma.major).toBeCloseTo(0.2, 12);
        expect(sigma.minor).toBeCloseTo(0.1, 12);
    });

    it('uses the metric of an oblique plane', () => {
        // Hexagonal a = b = 5, gamma = 120: an isotropic fractional H is not
        // isotropic in Angstrom. Compare with M H M^T computed directly.
        const eu = [5, 0, 0];
        const ev = [-2.5, 5 * Math.sqrt(3) / 2, 0];
        const h = [[1e-4, 0], [0, 1e-4]];
        const m = [[eu[0], ev[0]], [eu[1], ev[1]]];
        const c = [0, 1].map((i) => [0, 1].map((j) => (
            m[i][0] * (h[0][0] * m[j][0] + h[0][1] * m[j][1]) + m[i][1] * (h[1][0] * m[j][0] + h[1][1] * m[j][1])
        )));
        const mean = (c[0][0] + c[1][1]) / 2;
        const radius = Math.hypot((c[0][0] - c[1][1]) / 2, c[0][1]);
        const sigma = kernelSigmaAngstrom(h, eu, ev);
        expect(sigma.major).toBeCloseTo(Math.sqrt(mean + radius), 12);
        expect(sigma.minor).toBeCloseTo(Math.sqrt(mean - radius), 12);
        expect(sigma.major / sigma.minor).toBeCloseTo(Math.sqrt(3), 10);
    });
});
