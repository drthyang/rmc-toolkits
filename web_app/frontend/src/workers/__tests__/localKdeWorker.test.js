// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { augmentPeriodicImages, computeKde } from '../localKdeWorker';

// Deterministic pseudo-random stream so the fixtures are stable across runs.
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

const zSlicePayload = (points, overrides = {}) => ({
    id: 1,
    points,
    normal: [0, 0, 1],
    uVector: [1, 0, 0],
    vVector: [0, 1, 0],
    range: [0, 1],
    zCenter: 0.5,
    thickness: 0.2,
    bandwidth: 0.1,
    gridSize: 41,
    logScale: false,
    ...overrides
});

describe('localKdeWorker periodic boundaries', () => {
    it('keeps only images within the margin and maps them to their source atom', () => {
        const { augmented, sourceIndex } = augmentPeriodicImages(
            [{ x: 0.02, y: 0.5, z: 0.5 }],
            0.1
        );

        expect(augmented).toHaveLength(2);
        expect(sourceIndex).toEqual([0, 0]);
        expect(augmented.some((p) => Math.abs(p.x - 1.02) < 1e-12)).toBe(true);
    });

    it('wraps density across the cell boundary', async () => {
        const random = makeRandom(1);
        const points = [];
        for (let i = 0; i < 400; i += 1) {
            points.push({ x: random(), y: random(), z: 0.5 + 0.04 * (random() - 0.5) });
        }
        // Cluster hugging the x=0 face: its periodic image must contribute the
        // same density at the x=1 edge of the slice.
        for (let i = 0; i < 200; i += 1) {
            points.push({
                x: 0.01 + 0.02 * random(),
                y: 0.5 + 0.04 * (random() - 0.5),
                z: 0.5 + 0.04 * (random() - 0.5)
            });
        }

        const result = await computeKde(zSlicePayload(points));

        expect(result.slabCount).toBe(600);
        const mid = 20; // y = 0.5 row of the 41-point grid
        const leftEdge = result.density[mid][0];
        const rightEdge = result.density[mid][40];
        expect(leftEdge).toBeGreaterThan(0);
        // By periodicity the two edges see the cluster at the same distance.
        expect(Math.abs(rightEdge - leftEdge) / Math.max(rightEdge, leftEdge)).toBeLessThan(0.2);
        // Without wrapping the right edge only sees the uniform background; the
        // cluster image must raise it well above a cluster-free edge point.
        const farRight = result.density[36][40]; // y = 0.9, away from the cluster
        expect(rightEdge).toBeGreaterThan(2 * farRight);
    });

    it('finds atoms through depth wrapping when the slab crosses z=0', async () => {
        const random = makeRandom(2);
        // All atoms sit just below z=1; a slab centered at z=0 only finds them
        // through their periodic images at z-1.
        const points = [];
        for (let i = 0; i < 30; i += 1) {
            points.push({ x: random(), y: random(), z: 0.97 + 0.02 * random() });
        }

        const result = await computeKde(zSlicePayload(points, {
            zCenter: 0.0,
            thickness: 0.1,
            gridSize: 16
        }));

        expect(result.slabCount).toBe(30);
        expect(result.fitCount).toBeGreaterThan(0);
        expect(result.vmax).toBeGreaterThan(0);
    });
});
