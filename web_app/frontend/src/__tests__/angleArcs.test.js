// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The angle marks (components/angleArcs.js): an arc at the apex from one
// bond to the other, at a fraction of the shorter bond, clear of the site's
// cloud; all arcs one mesh, each as opaque as asked.

import { describe, expect, it } from 'vitest';
import { ARC_RADIUS_FRACTION, arcFrame, arcMesh, arcPoint } from '../components/angleArcs';
import { BOND_END_GAP } from '../components/bondTubes';

const GPU = { software: false };
const SOFTWARE = { software: true };
const distance = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
const direction = (from, to) => {
    const d = to.map((value, axis) => value - from[axis]);
    const length = Math.hypot(...d);
    return d.map((value) => value / length);
};

describe('an arc\'s frame', () => {
    // Apex at (1, 2, 3); A end 2.6 Å along x, C end 2.4 Å at 75° in the xy plane.
    const apex = [1, 2, 3];
    const theta = 75 * Math.PI / 180;
    const a = [apex[0] + 2.6, apex[1], apex[2]];
    const c = [apex[0] + 2.4 * Math.cos(theta), apex[1] + 2.4 * Math.sin(theta), apex[2]];
    const frame = arcFrame([...apex, ...a, ...c]);

    it('spans the angle between the bonds, at a fraction of the shorter one', () => {
        expect(frame.angle).toBeCloseTo(theta, 12);
        expect(frame.radius).toBeCloseTo(ARC_RADIUS_FRACTION * 2.4, 12);
        // Clear of the site's cloud (the bonds stop BOND_END_GAP short), well
        // short of the end atoms.
        expect(frame.radius).toBeGreaterThan(2 * BOND_END_GAP);
        expect(frame.radius).toBeLessThan(0.5 * 2.4);
    });

    it('starts on the A bond and ends on the C bond', () => {
        const start = arcPoint(frame, 0);
        const end = arcPoint(frame, frame.angle);
        direction(apex, start).forEach((value, axis) => expect(value).toBeCloseTo(direction(apex, a)[axis], 12));
        direction(apex, end).forEach((value, axis) => expect(value).toBeCloseTo(direction(apex, c)[axis], 12));
        // Every point at the same distance from the apex, in the bonds' plane.
        for (let k = 0; k <= 8; k += 1) {
            const point = arcPoint(frame, (k / 8) * frame.angle);
            expect(distance(point, apex)).toBeCloseTo(frame.radius, 12);
            expect(point[2]).toBeCloseTo(apex[2], 12);
        }
    });

    it('bends a straight angle the same way every time, and skips a bond of no length', () => {
        const straight = arcFrame([0, 0, 0, 2, 0, 0, -2, 0, 0]);
        expect(straight.angle).toBeCloseTo(Math.PI, 12);
        // Perpendicular to the bond, of unit length; its end on the C bond.
        expect(straight.w[0]).toBeCloseTo(0, 12);
        expect(Math.hypot(...straight.w)).toBeCloseTo(1, 12);
        expect(arcPoint(straight, Math.PI)[0]).toBeCloseTo(-straight.radius, 12);
        expect(arcFrame([0, 0, 0, 0, 0, 0, 1, 0, 0])).toBe(null);
    });
});

describe('the arcs\' mesh', () => {
    // Two arcs: 90° and 180°.
    const arcs = new Float64Array([0, 0, 0, 2, 0, 0, 0, 2, 0, 5, 5, 5, 7, 5, 5, 3, 5, 5]);

    it('draws every arc in one mesh, each with its own opacity', () => {
        const mesh = arcMesh(arcs, 0.05, { opacities: [1, 0.25], quality: GPU });
        expect(mesh.isMesh).toBe(true);
        const { position, color } = mesh.geometry.attributes;
        // 7.5° steps, 8 around: 90° → 12 steps, 180° → 24 steps.
        expect(position.count).toBe((13 + 25) * 8);
        expect(color.itemSize).toBe(4);
        expect(color.getW(0)).toBe(1);
        expect(color.getW(position.count - 1)).toBe(0.25);
        expect(mesh.material.vertexColors).toBe(true);
        expect(mesh.material.transparent).toBe(true);
    });

    it('keeps each tube on its arc, as thick as asked', () => {
        const mesh = arcMesh(arcs.subarray(0, 9), 0.05, { scales: [2], quality: GPU });
        const { position } = mesh.geometry.attributes;
        const radius = ARC_RADIUS_FRACTION * 2;
        for (let v = 0; v < position.count; v += 1) {
            const p = [position.getX(v), position.getY(v), position.getZ(v)];
            // Distance from the arc's circle (in the xy plane, about the apex).
            const off = Math.hypot(Math.hypot(p[0], p[1]) - radius, p[2]);
            expect(off).toBeCloseTo(0.1, 5);
        }
    });

    it('takes fewer vertices on a software renderer, and is null when nothing is drawn', () => {
        const soft = arcMesh(arcs, 0.05, { quality: SOFTWARE });
        expect(soft.geometry.attributes.position.count).toBe((7 + 13) * 5);
        expect(arcMesh(new Float64Array(0), 0.05, { quality: GPU })).toBe(null);
        expect(arcMesh(new Float64Array([0, 0, 0, 0, 0, 0, 1, 0, 0]), 0.05, { quality: GPU })).toBe(null);
    });
});
