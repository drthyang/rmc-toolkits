// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The PCA page's "Shell" paints the KDE density sampled on the p% ellipsoid. Where
// that ellipsoid pokes out of the sampled box (small Box, high Level) the old
// sampler clamped to the box face and painted the nearer, denser face value as
// hot caps along PC1 (pca.numerics.34). The sampler must report "no data" there.

import { describe, it, expect } from 'vitest';

import { sampleFieldTrilinear } from '../marchingCubes.js';

const n = 5;
// A linear field f = i + 2j + 3k is reproduced exactly by trilinear interpolation.
const field = new Float64Array(n * n * n);
for (let i = 0; i < n; i += 1) {
    for (let j = 0; j < n; j += 1) {
        for (let k = 0; k < n; k += 1) field[(i * n + j) * n + k] = i + 2 * j + 3 * k;
    }
}

describe('sampleFieldTrilinear', () => {
    it('interpolates inside the box, including on its faces', () => {
        expect(sampleFieldTrilinear(field, n, n, n, 1.25, 2.5, 3.75)).toBeCloseTo(1.25 + 5 + 11.25, 12);
        expect(sampleFieldTrilinear(field, n, n, n, 4, 4, 4)).toBeCloseTo(24, 12);
        expect(sampleFieldTrilinear(field, n, n, n, 0, 0, 0)).toBe(0);
    });

    it('returns NaN outside the box instead of the clamped face value', () => {
        expect(Number.isNaN(sampleFieldTrilinear(field, n, n, n, 4.5, 2, 2))).toBe(true);
        expect(Number.isNaN(sampleFieldTrilinear(field, n, n, n, 2, -0.2, 2))).toBe(true);
        expect(Number.isNaN(sampleFieldTrilinear(field, n, n, n, 2, 2, 7))).toBe(true);
    });

    it('supports a non-cubic grid', () => {
        const [nx, ny, nz] = [3, 4, 2];
        const f = new Float64Array(nx * ny * nz).map((_, idx) => idx);
        // index = (i*ny + j)*nz + k is itself linear in (i, j, k).
        expect(sampleFieldTrilinear(f, nx, ny, nz, 1.5, 2.5, 0.5)).toBeCloseTo((1.5 * ny + 2.5) * nz + 0.5, 12);
    });
});
