// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The custom slice is a Miller-plane family (h k l), and is labelled as one.
// An atom's depth is its fractional coordinate dotted with h, so the slab is
// bounded by the lattice planes x . h = const, whose normal is h1 a* + h2 b* +
// h3 c*; the real-space direction [h k l] is a different vector in a
// non-orthogonal cell. The page used to print the input as "[1 1 0]" under the
// heading "Direction"; these checks pin the plane semantics and the notation.

import { describe, expect, it } from 'vitest';
import { computeKde } from '../localKdeWorker';
import { millerPlaneFileLabel, millerPlaneLabel } from '../slabSelection';

const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const angleDeg = (a, b) => (Math.acos(dot(a, b) / Math.sqrt(dot(a, a) * dot(b, b))) * 180) / Math.PI;

describe('custom slice = Miller plane (h k l)', () => {
    it('is labelled with parentheses on the canvas and in file names', () => {
        expect(millerPlaneLabel([1, 1, 0])).toBe('(1 1 0)');
        expect(millerPlaneFileLabel([1, 1, 0])).toBe('(1_1_0)');
        expect(millerPlaneFileLabel([-1, 0.5, 2])).toBe('(-1_0.5_2)');
        expect(millerPlaneFileLabel([Number.NaN, 1, 0])).toBe('(0_1_0)');
        expect(millerPlaneLabel([1, 0, 0])).not.toContain('[');
    });

    it('selects by x . h, so its normal is the reciprocal vector, not [h k l]', async () => {
        // Hexagonal cell (a = b = 5.66 A, gamma = 120 deg): two atoms with the
        // same fractional x lie on the same (1 0 0) plane although they sit at
        // different positions along the real-space direction a.
        const a = [5.66, 0, 0];
        const b = [-2.83, 5.66 * Math.sqrt(3) / 2, 0];
        const c = [0, 0, 4.53];
        const volume = dot(a, cross(b, c));
        const aStar = cross(b, c).map((value) => value / volume);
        expect(angleDeg(aStar, a)).toBeCloseTo(30, 10);

        const points = [
            ...Array.from({ length: 20 }, (_, i) => ({ x: 0.5, y: i / 20, z: (i * 7 % 20) / 20 })),
            ...Array.from({ length: 20 }, (_, i) => ({ x: 0.9, y: i / 20, z: (i * 3 % 20) / 20 }))
        ];
        const result = await computeKde({
            points,
            normal: [1, 0, 0],
            uVector: [0, 1, 0],
            vVector: [0, 0, 1],
            range: [0, 1],
            zCenter: 0.5,
            thickness: 0.08,
            bandwidth: 0.05,
            gridSize: 16,
            logScale: false
        });
        expect(result.slabCount).toBe(20);
    });
});
