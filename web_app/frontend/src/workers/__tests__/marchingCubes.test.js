// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, it, expect } from 'vitest';

import { marchingCubes } from '../marchingCubes.js';

// Sample a scalar field on an n^3 grid in C order over (x, y, z).
const sampleField = (n, fn) => {
    const field = new Float64Array(n * n * n);
    for (let i = 0; i < n; i += 1) {
        for (let j = 0; j < n; j += 1) {
            for (let k = 0; k < n; k += 1) {
                field[(i * n + j) * n + k] = fn(i, j, k);
            }
        }
    }
    return field;
};

describe('marchingCubes', () => {
    it('contours a sphere into a closed surface of the right radius', () => {
        const n = 32;
        const center = (n - 1) / 2;
        const radius = 10;
        // Field high inside, low outside, so the surface encloses the ball.
        const field = sampleField(n, (i, j, k) => {
            const r = Math.hypot(i - center, j - center, k - center);
            return radius - r;
        });
        const { positions, count } = marchingCubes(field, n, n, n, 0);
        expect(count).toBeGreaterThan(100);

        // Every vertex sits ~radius from the center (grid units).
        let maxErr = 0;
        for (let v = 0; v < count; v += 1) {
            const r = Math.hypot(
                positions[3 * v] - center,
                positions[3 * v + 1] - center,
                positions[3 * v + 2] - center
            );
            maxErr = Math.max(maxErr, Math.abs(r - radius));
        }
        // Linear edge interpolation places vertices within a fraction of a cell.
        expect(maxErr).toBeLessThan(0.6);
    });

    it('normals point outward from a density blob', () => {
        const n = 24;
        const center = (n - 1) / 2;
        // A Gaussian blob: high at the center, decaying outward.
        const field = sampleField(n, (i, j, k) => {
            const r2 = (i - center) ** 2 + (j - center) ** 2 + (k - center) ** 2;
            return Math.exp(-r2 / 40);
        });
        const iso = Math.exp(-(6 ** 2) / 40);
        const { positions, normals, count } = marchingCubes(field, n, n, n, iso);
        expect(count).toBeGreaterThan(100);

        // Outward radial direction should align with the surface normal.
        let aligned = 0;
        for (let v = 0; v < count; v += 1) {
            const rx = positions[3 * v] - center;
            const ry = positions[3 * v + 1] - center;
            const rz = positions[3 * v + 2] - center;
            const rlen = Math.hypot(rx, ry, rz) || 1;
            const dotProduct = (rx * normals[3 * v] + ry * normals[3 * v + 1] + rz * normals[3 * v + 2]) / rlen;
            if (dotProduct > 0) aligned += 1;
        }
        expect(aligned / count).toBeGreaterThan(0.95);
    });

    it('returns an empty surface when the level is never crossed', () => {
        const n = 8;
        const field = sampleField(n, () => 1);
        const { count } = marchingCubes(field, n, n, n, 5);
        expect(count).toBe(0);
    });
});
