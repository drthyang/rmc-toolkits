// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The angle marks (components/angleArcs.js): a flat wedge at the apex from
// one bond to the other, at a fraction of the shorter bond, clear of the
// site's cloud; all wedges one mesh, each as strong as asked.

import { describe, expect, it } from 'vitest';
import * as THREE from 'three';
import {
    ARC_RADIUS_FRACTION, arcFrame, arcPoint, WEDGE_COLOR, WEDGE_FILL_OPACITY, WEDGE_RIM_FRACTION, WEDGE_RIM_OPACITY, wedgeMesh
} from '../components/angleArcs';
import { PLOT_PALETTE } from '../plotPalette';
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

describe('the wedges\' mesh', () => {
    // Two angles: 90° and 180°.
    const arcs = new Float64Array([0, 0, 0, 2, 0, 0, 0, 2, 0, 5, 5, 5, 7, 5, 5, 3, 5, 5]);
    const radius = ARC_RADIUS_FRACTION * 2;
    const pointAt = (position, v) => [position.getX(v), position.getY(v), position.getZ(v)];

    it('draws every wedge in one mesh, each with its own strength', () => {
        const mesh = wedgeMesh(arcs, { strengths: [1, 0.25], quality: GPU });
        expect(mesh.isMesh).toBe(true);
        const { position, color } = mesh.geometry.attributes;
        // 4° steps: 90° → 23 steps, 180° → 45; per wedge the apex, the fan's
        // arc and the rim's two edges.
        expect(position.count).toBe((3 * 23 + 4) + (3 * 45 + 4));
        expect(mesh.geometry.index.count).toBe(9 * (23 + 45));
        expect(color.itemSize).toBe(4);
        // The fill is light, the rim nearly opaque; both scale with strength.
        expect(color.getW(0)).toBeCloseTo(WEDGE_FILL_OPACITY, 6);
        expect(color.getW(3 * 23 + 3)).toBeCloseTo(WEDGE_RIM_OPACITY, 6);
        expect(color.getW(3 * 23 + 4)).toBeCloseTo(0.25 * WEDGE_FILL_OPACITY, 6);
        expect(color.getW(position.count - 1)).toBeCloseTo(0.25 * WEDGE_RIM_OPACITY, 6);
    });

    it('is flat, unlit and seen from both faces, in the angle plot\'s blue with a darker rim', () => {
        const mesh = wedgeMesh(arcs.subarray(0, 9), { quality: GPU });
        expect(WEDGE_COLOR).toBe(PLOT_PALETTE[0]);
        expect(mesh.material.isMeshBasicMaterial).toBe(true);
        expect(mesh.material.vertexColors).toBe(true);
        expect(mesh.material.transparent).toBe(true);
        expect(mesh.material.depthWrite).toBe(false);
        expect(mesh.material.side).toBe(THREE.DoubleSide);
        const { color } = mesh.geometry.attributes;
        const blue = new THREE.Color(WEDGE_COLOR);
        expect([color.getX(0), color.getY(0), color.getZ(0)].map((v) => Number(v.toFixed(5))))
            .toEqual(blue.toArray().map((v) => Number(v.toFixed(5))));
        const last = color.count - 1;
        expect(color.getX(last) + color.getY(last) + color.getZ(last)).toBeLessThan(blue.r + blue.g + blue.b);
    });

    it('fills from the apex to the arc, with the rim just inside it', () => {
        const mesh = wedgeMesh(arcs.subarray(0, 9), { quality: GPU });
        const { position } = mesh.geometry.attributes;
        expect(pointAt(position, 0)).toEqual([0, 0, 0]);
        for (let v = 1; v < position.count; v += 1) {
            const p = pointAt(position, v);
            // In the bonds' plane (xy), between the two bonds.
            expect(p[2]).toBeCloseTo(0, 6);
            expect(p[0]).toBeGreaterThan(-1e-6);
            expect(p[1]).toBeGreaterThan(-1e-6);
            const r = Math.hypot(p[0], p[1]);
            const inner = v > 24 && (v - 25) % 2 === 0;
            expect(r).toBeCloseTo(inner ? radius * (1 - WEDGE_RIM_FRACTION) : radius, 6);
        }
    });

    it('takes fewer vertices on a software renderer, and skips what is not drawn', () => {
        const soft = wedgeMesh(arcs, { quality: SOFTWARE });
        // 10° steps: 9 and 18.
        expect(soft.geometry.attributes.position.count).toBe((3 * 9 + 4) + (3 * 18 + 4));
        // A strength of zero draws nothing.
        expect(wedgeMesh(arcs, { strengths: [1, 0], quality: GPU }).geometry.attributes.position.count).toBe(3 * 23 + 4);
        expect(wedgeMesh(arcs, { strengths: [0, 0], quality: GPU })).toBe(null);
        expect(wedgeMesh(new Float64Array(0), { quality: GPU })).toBe(null);
        expect(wedgeMesh(new Float64Array([0, 0, 0, 0, 0, 0, 1, 0, 0]), { quality: GPU })).toBe(null);
    });
});
