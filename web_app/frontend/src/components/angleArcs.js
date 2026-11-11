// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The angle marks of the folded cell's highlight (FoldedCellPanel): for each
// drawn triplet, a flat shaded wedge at its central atom between its two
// bonds, so the view shows which angle at a site is in the range, not only
// which bonds. angleHighlight.js highlightForBins gives the angles (`arcs`:
// apex, A end, C end) and their shares; a wedge is drawn as strongly as the
// bonds are (highlightStrength: the square root of the share).
//
// - Radius ARC_RADIUS_FRACTION of the shorter bond: past the site's thermal
//   cloud (and the bond's BOND_END_GAP), well short of the end atoms.
// - A light fill in the angle plot's blue, with a darker rim along the arc;
//   unlit and double-sided, so a wedge looks the same from either face.
// - All wedges are one merged mesh (one draw call), each one's strength a
//   vertex alpha; a software renderer takes fewer segments along the arc.
// - An angle near 180° has no plane of its own: its wedge bends towards the
//   cell axis least along its bonds, the same way every time.

import * as THREE from 'three';
import { PLOT_PALETTE } from '../plotPalette';
import { renderQuality } from './bondTubes';

// The arc's radius as a fraction of its shorter bond.
export const ARC_RADIUS_FRACTION = 0.35;
// The wedges' colour: the angle plot's line and range band.
export const WEDGE_COLOR = PLOT_PALETTE[0];
// A wedge of strength 1: its fill's and its rim's opacity.
export const WEDGE_FILL_OPACITY = 0.2;
export const WEDGE_RIM_OPACITY = 0.95;
// The rim's width as a fraction of the wedge's radius, and how far its
// colour goes from WEDGE_COLOR towards RIM_SHADE.
export const WEDGE_RIM_FRACTION = 0.043;
const RIM_SHADE = '#0b1b33';
const RIM_SHADE_MIX = 0.28;
// Degrees per segment along a wedge's arc (GPU, software renderer).
const ARC_STEP_DEGREES = { gpu: 4, software: 10 };
// Below this sine the two bonds are taken as collinear.
const COLLINEAR_SINE = 1e-6;

const sub = (p, q) => [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
const dot = (p, q) => p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
const scaled = (p, k) => [p[0] * k, p[1] * k, p[2] * k];
const unitOf = (p) => {
    const length = Math.hypot(...p);
    return length > 0 ? scaled(p, 1 / length) : null;
};

/**
 * An arc's frame from nine numbers (apex, A end, C end): its centre (the
 * apex), radius, angle (radians) and the in-plane unit vectors u (towards A)
 * and w (perpendicular, on C's side), so a point at φ ∈ [0, angle] is
 * apex + radius·(cos φ·u + sin φ·w). Null when a bond has no length.
 */
export const arcFrame = (arcs, index = 0) => {
    const at = index * 9;
    const apex = [arcs[at], arcs[at + 1], arcs[at + 2]];
    const toA = sub([arcs[at + 3], arcs[at + 4], arcs[at + 5]], apex);
    const toC = sub([arcs[at + 6], arcs[at + 7], arcs[at + 8]], apex);
    const u = unitOf(toA);
    const c = unitOf(toC);
    if (!u || !c) return null;
    const cosine = Math.min(1, Math.max(-1, dot(u, c)));
    const angle = Math.acos(cosine);
    let w = unitOf(sub(c, scaled(u, cosine)));
    if (!w || Math.sin(angle) < COLLINEAR_SINE) {
        // Collinear: bend towards the cell axis least along the bond.
        const axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        const axis = axes.reduce((best, next) => (Math.abs(dot(next, u)) < Math.abs(dot(best, u)) ? next : best));
        w = unitOf(sub(axis, scaled(u, dot(axis, u))));
    }
    const radius = ARC_RADIUS_FRACTION * Math.min(Math.hypot(...toA), Math.hypot(...toC));
    return { center: apex, radius, angle, u, w };
};

/** The point at φ (radians) along an arc's frame. */
export const arcPoint = ({ center, radius, u, w }, phi) => [0, 1, 2].map(
    (axis) => center[axis] + radius * (Math.cos(phi) * u[axis] + Math.sin(phi) * w[axis])
);

/**
 * One merged mesh for every angle (nine numbers each): a filled sector from
 * the apex to the arc, and a rim band along the arc, in `color`, each as
 * strong as `strengths[i]` (0..1; default 1). Angles are in the given order;
 * a degenerate one (a bond of no length) or one of no strength is skipped.
 * Null when nothing is drawn.
 */
export const wedgeMesh = (arcs, { strengths = null, quality = renderQuality(), color = WEDGE_COLOR } = {}) => {
    const count = Math.floor(arcs.length / 9);
    const stepDegrees = quality.software ? ARC_STEP_DEGREES.software : ARC_STEP_DEGREES.gpu;
    const fill = new THREE.Color(color);
    const rim = new THREE.Color(color).lerp(new THREE.Color(RIM_SHADE), RIM_SHADE_MIX);
    const positions = [];
    const colors = [];
    const indices = [];
    for (let i = 0; i < count; i += 1) {
        const frame = arcFrame(arcs, i);
        const strength = strengths ? strengths[i] : 1;
        if (!frame || !(frame.angle > 0) || !(strength > 0)) continue;
        const along = Math.max(2, Math.ceil(frame.angle * 180 / Math.PI / stepDegrees - 1e-9));
        // The fill: a fan from the apex.
        const apexAt = positions.length / 3;
        positions.push(...frame.center);
        colors.push(fill.r, fill.g, fill.b, WEDGE_FILL_OPACITY * strength);
        for (let k = 0; k <= along; k += 1) {
            positions.push(...arcPoint(frame, (k / along) * frame.angle));
            colors.push(fill.r, fill.g, fill.b, WEDGE_FILL_OPACITY * strength);
        }
        for (let k = 0; k < along; k += 1) indices.push(apexAt, apexAt + 1 + k, apexAt + 2 + k);
        // The rim: a thin band in the wedge's plane, ending on the arc.
        const inner = { ...frame, radius: frame.radius * (1 - WEDGE_RIM_FRACTION) };
        const rimAt = positions.length / 3;
        for (let k = 0; k <= along; k += 1) {
            const phi = (k / along) * frame.angle;
            positions.push(...arcPoint(inner, phi), ...arcPoint(frame, phi));
            for (let side = 0; side < 2; side += 1) colors.push(rim.r, rim.g, rim.b, WEDGE_RIM_OPACITY * strength);
        }
        for (let k = 0; k < along; k += 1) {
            const a = rimAt + 2 * k;
            indices.push(a, a + 1, a + 2, a + 2, a + 1, a + 3);
        }
    }
    if (!indices.length) return null;
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 4));
    geometry.setIndex(indices);
    geometry.computeBoundingSphere();
    return new THREE.Mesh(geometry, new THREE.MeshBasicMaterial({
        vertexColors: true,
        transparent: true,
        side: THREE.DoubleSide,
        depthWrite: false
    }));
};
