// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Atoms exactly on a slab face are in the slab: the worker, the Slab-In-Cell
// highlight (both through isInSlab) and kde.py share one inclusive test with a
// small tolerance. Python twin: tests/test_kde_slab_faces.py.

import { describe, expect, it } from 'vitest';
import { computeKde } from '../localKdeWorker';
import { SLAB_FACE_TOLERANCE, isInSlab } from '../slabSelection';

// 16 atoms on each of the layers z = k/8 (exact binary fractions).
const idealLayers = () => {
    const points = [];
    for (let k = 0; k < 8; k += 1) {
        for (let i = 0; i < 4; i += 1) {
            for (let j = 0; j < 4; j += 1) points.push({ x: 0.1 + 0.2 * i, y: 0.15 + 0.2 * j, z: k / 8 });
        }
    }
    return points;
};

// Atoms in the slab by exact integer arithmetic in thousandths (with wrap).
const exactCount = (centerMilli, thicknessMilli) => {
    let count = 0;
    for (let k = 0; k < 8; k += 1) {
        const depth = 125 * k;
        if ([-1, 0, 1].some((shift) => 2 * Math.abs(depth + shift * 1000 - centerMilli) <= thicknessMilli)) count += 16;
    }
    return count;
};

const cSlice = (points, zCenter, thickness) => computeKde({
    points,
    normal: [0, 0, 1],
    uVector: [1, 0, 0],
    vVector: [0, 1, 0],
    range: [0, 1],
    zCenter,
    thickness,
    bandwidth: 0.05,
    gridSize: 16,
    logScale: false
});

describe('slab faces', () => {
    it('includes the face atoms of an ideal configuration at every slider position', async () => {
        const points = idealLayers();
        const misses = [];
        for (const thicknessMilli of [80, 100, 150, 250]) {
            for (let centerMilli = 0; centerMilli <= 1000; centerMilli += 5) {
                const result = await cSlice(points, centerMilli / 1000, thicknessMilli / 1000);
                const expected = exactCount(centerMilli, thicknessMilli);
                if (result.slabCount !== expected) misses.push([centerMilli / 1000, thicknessMilli / 1000, result.slabCount, expected]);
            }
        }
        expect(misses).toEqual([]);
    });

    it('the Slab-In-Cell test includes the same face atoms', () => {
        // z = 0.125 against zCenter = 0.165, thickness = 0.08: |0.125 - 0.165| rounds to 0.04000000000000001.
        expect(Math.abs(0.125 - 0.165) <= 0.08 / 2).toBe(false);
        expect(isInSlab(0.125, 0.165, 0.08)).toBe(true);
        expect(isInSlab(0.25, 0.21, 0.08)).toBe(true);
        expect(isInSlab(0.54 + 2 * SLAB_FACE_TOLERANCE, 0.5, 0.08)).toBe(false);
    });

    it('echoes the slider fractions as z/dz', async () => {
        const result = await cSlice(idealLayers(), 0.165, 0.08);
        expect(result.z).toBe(0.165);
        expect(result.dz).toBe(0.08);
        expect(result.slabCount).toBe(16);
    });
});
