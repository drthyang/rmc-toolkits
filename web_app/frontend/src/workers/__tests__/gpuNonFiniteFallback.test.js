// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// GPU KDE must always degrade gracefully: a density map read back from the GPU
// with a non-finite node (a float32 overflow, a misbehaving device) is dropped
// and the float64 CPU loop draws the map instead, so a NaN never reaches the
// page. WebGPU is absent under node, so the GPU entry points are mocked.

import { describe, expect, it, vi } from 'vitest';

vi.mock('../gpuKde', () => ({
    shouldUseGpu: () => true,
    computeDensityGpu: async ({ grid }) => Array.from({ length: grid }, () => new Array(grid).fill(Number.NaN))
}));

const { computeKde } = await import('../localKdeWorker');

const points = [];
for (const x0 of [0.2, 0.5, 0.8]) {
    for (const dx of [-0.03, 0.02]) {
        for (const dy of [-0.02, 0.03]) points.push({ x: x0 + dx, y: 0.4 + x0 * 0.2 + dy, z: 0.5, element: 'X' });
    }
}

describe('GPU density fallback', () => {
    it('replaces a GPU map with a non-finite node by the CPU map', async () => {
        const result = await computeKde({
            points,
            normal: [0, 0, 1],
            uVector: [1, 0, 0],
            vVector: [0, 1, 0],
            range: [0, 1],
            zCenter: 0.5,
            thickness: 0.08,
            bandwidth: 0.3,
            gridSize: 32,
            logScale: false
        });
        expect(result.backend).toBe('cpu');
        expect(result.density.flat().every(Number.isFinite)).toBe(true);
        expect(result.vmax).toBeGreaterThan(0);
        expect(result.warnings.map((warning) => warning.code)).not.toContain('unresolved');
    });
});
