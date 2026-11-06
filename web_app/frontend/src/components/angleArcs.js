// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The angle marks of the folded cell's highlight (FoldedCellPanel): for each
// drawn triplet, a thin lit arc at its central atom from one bond to the
// other, so the view shows which angle at a site is in the range, not only
// which bonds. angleHighlight.js highlightForBins gives the arcs (apex, A end,
// C end) and their shares; the arc is drawn as strongly as the bonds are
// (highlightStrength: the square root of the share).
//
// - Radius ARC_RADIUS_FRACTION of the shorter bond: past the site's thermal
//   cloud (and the bond's BOND_END_GAP), well short of the end atoms.
// - All arcs are one merged tube mesh (one draw call), each arc's opacity a
//   vertex alpha; a software renderer takes fewer segments around and along.
// - An angle near 180° has no plane of its own: its arc bends towards the
//   cell axis least along its bonds, the same way every time.

import * as THREE from 'three';
import { renderQuality } from './bondTubes';

// The arc's radius as a fraction of its shorter bond.
export const ARC_RADIUS_FRACTION = 0.35;
// The arcs' colour: apart from the bond-role blue and orange and the
// elements' colours, readable on light and dark.
export const ARC_COLOR = '#c026d3';
// Degrees per segment along an arc (GPU, software renderer).
const ARC_STEP_DEGREES = { gpu: 7.5, software: 15 };
// Below this sine the two bonds are taken as collinear.
const COLLINEAR_SINE = 1e-6;

const sub = (p, q) => [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
const dot = (p, q) => p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
const cross = (p, q) => [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]];
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
 * One merged tube mesh for every arc (nine numbers each), `tubeRadius` thick
 * × `scales[i]`, opacity `opacities[i]` (both default 1), in ARC_COLOR.
 * Arcs are in the given order; a degenerate one (a bond of no length) is
 * skipped. Null when nothing is drawn.
 */
export const arcMesh = (arcs, tubeRadius, {
    scales = null, opacities = null, quality = renderQuality(), color = ARC_COLOR
} = {}) => {
    const count = Math.floor(arcs.length / 9);
    const step = (quality.software ? ARC_STEP_DEGREES.software : ARC_STEP_DEGREES.gpu) * Math.PI / 180;
    const around = quality.software ? 5 : 8;
    const frames = [];
    let vertexCount = 0;
    let indexCount = 0;
    for (let i = 0; i < count; i += 1) {
        const frame = arcFrame(arcs, i);
        if (!frame || !(frame.angle > 0)) continue;
        const along = Math.max(2, Math.ceil(frame.angle / step));
        frames.push({ frame, along, index: i });
        vertexCount += (along + 1) * around;
        indexCount += along * around * 6;
    }
    if (!frames.length) return null;
    const positions = new Float32Array(vertexCount * 3);
    const normals = new Float32Array(vertexCount * 3);
    const colors = new Float32Array(vertexCount * 4);
    const indices = vertexCount > 65535 ? new Uint32Array(indexCount) : new Uint16Array(indexCount);
    let v = 0;
    let n = 0;
    frames.forEach(({ frame, along, index }) => {
        const r = tubeRadius * (scales ? scales[index] : 1);
        const alpha = opacities ? opacities[index] : 1;
        // The tube's cross-section at φ: the in-plane radial direction and
        // the plane's normal.
        const normal = cross(frame.u, frame.w);
        const first = v;
        for (let k = 0; k <= along; k += 1) {
            const phi = (k / along) * frame.angle;
            const centre = arcPoint(frame, phi);
            const radial = [0, 1, 2].map((axis) => Math.cos(phi) * frame.u[axis] + Math.sin(phi) * frame.w[axis]);
            for (let j = 0; j < around; j += 1) {
                const theta = (j / around) * Math.PI * 2;
                const out = [0, 1, 2].map((axis) => Math.cos(theta) * radial[axis] + Math.sin(theta) * normal[axis]);
                positions.set([centre[0] + r * out[0], centre[1] + r * out[1], centre[2] + r * out[2]], v * 3);
                normals.set(out, v * 3);
                colors.set([1, 1, 1, alpha], v * 4);
                v += 1;
            }
        }
        for (let k = 0; k < along; k += 1) {
            for (let j = 0; j < around; j += 1) {
                const a = first + k * around + j;
                const b = first + k * around + ((j + 1) % around);
                indices.set([a, a + around, b, b, a + around, b + around], n);
                n += 6;
            }
        }
    });
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('normal', new THREE.BufferAttribute(normals, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 4));
    geometry.setIndex(new THREE.BufferAttribute(indices, 1));
    geometry.computeBoundingSphere();
    const material = new THREE.MeshPhongMaterial({
        color,
        emissive: color,
        emissiveIntensity: 0.35,
        specular: 0x3a3a3a,
        shininess: 40,
        vertexColors: true,
        transparent: true,
        depthWrite: false
    });
    return new THREE.Mesh(geometry, material);
};
