// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Shared axis palettes and rod builders for the Three.js panels (PCA Ellipsoid
// viewport, Site-ellipsoids unit-cell picker, Orientation sphere). One module
// so PC1/PC2/PC3 and a/b/c always read as the same objects across every panel.

import * as THREE from 'three';

export const TRIAD_COLORS = [0xd64545, 0x3fa34d, 0x3f7fd6]; // PC1 red, PC2 green, PC3 blue
// CSS twins of TRIAD_COLORS, used to tint PC rows/labels so they read as the
// same axes shown by the 3D triads.
export const PC_CSS_COLORS = ['#d64545', '#3fa34d', '#3f7fd6'];

// Crystallographic-axis palette — deliberately different from the PC tricolor
// so the two frames never read as the same thing. gold a, teal b, orchid c.
export const CELL_AXIS_COLORS = [0xe0a419, 0x18a3a0, 0xb15ad8];
export const CELL_AXIS_CSS = ['#e0a419', '#18a3a0', '#b15ad8'];
export const CELL_AXIS_LABELS = ['a', 'b', 'c'];

// Contrast highlight for the selected atom in the unit-cell view — a warm color
// that stands out against the cool element palette.
export const SELECTED_ATOM_COLOR = 0xff7a1a;

export const TRIAD_UP = new THREE.Vector3(0, 1, 0);

// PC1/PC2/PC3 triad as thin rods from `origin` along the principal axes (rows).
// Thin cylinders rather than lines so a little thickness reads at any zoom
// without dominating.
export const buildAxisTriad = (origin, axes, length, radius) => {
    const meshes = [];
    for (let a = 0; a < 3; a += 1) {
        const dir = new THREE.Vector3(axes[a][0], axes[a][1], axes[a][2]).normalize();
        const mesh = new THREE.Mesh(
            new THREE.CylinderGeometry(radius, radius, length, 10),
            new THREE.MeshBasicMaterial({ color: TRIAD_COLORS[a] })
        );
        mesh.position.copy(origin).addScaledVector(dir, length / 2);
        mesh.quaternion.setFromUnitVectors(TRIAD_UP, dir);
        meshes.push(mesh);
    }
    return meshes;
};

// Crystallographic a/b/c triad: thin rods from `origin` along the unit-cell
// edge directions (normalised, so all three read at one scale), each in its
// own color. The rods carry only the cell ORIENTATION in the shared Cartesian
// frame; a/b/c are identified by color via the panel legends.
export const buildCrystalAxes = (origin, unitCell, length, radius) => {
    const rods = [];
    for (let a = 0; a < 3; a += 1) {
        const raw = new THREE.Vector3(unitCell[a][0], unitCell[a][1], unitCell[a][2]);
        if (raw.lengthSq() < 1e-12) continue;
        const dir = raw.normalize();
        const rod = new THREE.Mesh(
            new THREE.CylinderGeometry(radius, radius, length, 10),
            new THREE.MeshBasicMaterial({ color: CELL_AXIS_COLORS[a] })
        );
        rod.position.copy(origin).addScaledVector(dir, length / 2);
        rod.quaternion.setFromUnitVectors(TRIAD_UP, dir);
        rods.push(rod);
    }
    return rods;
};
