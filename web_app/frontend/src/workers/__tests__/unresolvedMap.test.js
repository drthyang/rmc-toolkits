// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A map whose kernel falls between all grid nodes holds only Gaussian tails
// (1e-14 on the GaNb4Se8 AVERAGE Nb layer). The worker, like kde.py
// (tests/test_kde_unresolved.py, same point sets), flags it `unresolved` when
// the grid-summed linear density is below UNRESOLVED_MASS_LIMIT (about 1 for a
// resolved map) and draws no contours over it, in either scale.

import { existsSync, readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { structureFromRmc6f } from '../../browserData.js';
import { KDE_WARNINGS, KERNEL_MIN_SIGMA, UNRESOLVED_MASS_LIMIT, computeKde } from '../localKdeWorker';

const REPO_ROOT = fileURLToPath(new URL('../../../../../', import.meta.url));
const AVERAGE_5K = `${REPO_ROOT}data/5K_try1/GaNb4Se8_5KAVERAGE.rmc6f`;

// Two sites on y = 0.5, midway between the node rows of a 32-point grid.
const twoSiteNeedle = (spread) => {
    const points = [];
    for (const x0 of [0.2, 0.8]) {
        for (const dx of [-0.01, 0.01]) {
            for (const dy of [-spread, spread]) points.push({ x: x0 + dx, y: 0.5 + dy, z: 0.5, element: 'X' });
        }
    }
    return points;
};

const cSlice = (points, overrides = {}) => computeKde({
    points,
    normal: [0, 0, 1],
    uVector: [1, 0, 0],
    vVector: [0, 1, 0],
    range: [0, 1],
    zCenter: 0.5,
    thickness: 0.08,
    bandwidth: 0.076,
    gridSize: 32,
    logScale: false,
    ...overrides
});

const codes = (result) => result.warnings.map((warning) => warning.code);

const gridMass = (result) => {
    const [xMin, xMax, yMin, yMax] = result.extent;
    const cell = ((xMax - xMin) / (result.grid - 1)) * ((yMax - yMin) / (result.grid - 1));
    let sum = 0;
    for (const row of result.density) {
        for (const value of row) sum += result.log ? 10 ** value - 1e-12 : value;
    }
    return sum * cell;
};

describe('unresolved (between-node) kernels', () => {
    it('flags a kernel that misses every node and draws no contours, in both scales', async () => {
        for (const logScale of [false, true]) {
            const result = await cSlice(twoSiteNeedle(0.01), { logScale });
            expect(result.message).toBeNull();
            expect(result.kernel).not.toBeNull();
            expect(codes(result)).toEqual(['subgrid', 'unresolved']);
            expect(result.warnings[1].message).toBe(KDE_WARNINGS.unresolved);
            expect(result.contours).toEqual([]);
            expect(result.vmax).toBeGreaterThan(result.vmin);
            expect(gridMass(result)).toBeLessThan(UNRESOLVED_MASS_LIMIT);
        }
    });

    it('flags a fully underflowed map too', async () => {
        const result = await cSlice(twoSiteNeedle(0.004), { bandwidth: 0.03 });
        expect(result.vmax).toBe(0);
        expect(codes(result)).toEqual(['subgrid', 'unresolved']);
        expect(result.contours).toEqual([]);
    });

    it('returns a kernel below the evaluable floor as a flagged zero map, never NaN', async () => {
        // kde.py: test_a_kernel_below_the_evaluable_floor_is_a_flagged_zero_map.
        // At bw = 1e-200 the normalizer 1 / (2 pi det L) overflows and every
        // node used to come out Inf * 0 = NaN.
        for (const bandwidth of [1e-30, 1e-200]) {
            const result = await cSlice(twoSiteNeedle(0.01), { bandwidth });
            expect(result.density.flat().every((value) => value === 0)).toBe(true);
            expect(result.message).toBeNull();
            expect(result.kernel).not.toBeNull();
            expect(result.kernel.sigmaMinor).toBeLessThan(KERNEL_MIN_SIGMA);
            expect(result.fitCount).toBe(8);
            expect(codes(result)).toEqual(['subgrid', 'unresolved']);
            expect(result.contours).toEqual([]);
        }
    });

    it('evaluates a narrow kernel above the floor', async () => {
        const points = twoSiteNeedle(0.01);
        points[0] = { ...points[0], x: 0, y: 0 };
        const result = await cSlice(points, { bandwidth: 1e-6 });
        expect(result.kernel.sigmaMinor).toBeGreaterThan(KERNEL_MIN_SIGMA);
        expect(result.vmax).toBeGreaterThan(1e9);
        expect(result.density.flat().every(Number.isFinite)).toBe(true);
    });

    it('leaves the same layer alone when a wider kernel reaches the nodes', async () => {
        const result = await cSlice(twoSiteNeedle(0.01), { bandwidth: 0.5 });
        expect(codes(result)).not.toContain('unresolved');
        expect(gridMass(result)).toBeGreaterThan(0.1);
        expect(result.contours).toHaveLength(8);
    });

    it('attaches no warning to a declined slab', async () => {
        const result = await cSlice(twoSiteNeedle(0.01), { bandwidth: 0 });
        expect(result.message).not.toBeNull();
        expect(result.warnings).toEqual([]);
    });

    (existsSync(AVERAGE_5K) ? it : it.skip)('flags the GaNb4Se8 AVERAGE Nb layer needle', async () => {
        const text = readFileSync(AVERAGE_5K, 'utf8');
        const points = structureFromRmc6f({ text, path: AVERAGE_5K }, 1_000_000).points
            .filter((point) => point.element === 'Nb');
        const result = await cSlice(points, { zCenter: 0.15, bandwidth: 0.03, gridSize: 120, logScale: true });
        expect(result.kernel.sigmaMinor).toBeLessThan(1e-5);
        expect(codes(result)).toEqual(['subgrid', 'unresolved']);
        expect(result.contours).toEqual([]);
    });
});
