// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import API_BASE_URL from '../api';
import { isStaticMode } from '../browserData';
import { COLORMAP_NAMES, getLut, sampleColormap } from '../colormaps';
import { buildElementColors, DEFAULT_ELEMENT_COLOR } from '../atomColors';
import { marchingCubes } from '../workers/marchingCubes';
import { downloadBlob, sanitizeFilename, saveCanvasAsPng } from '../figureExport';
import InfoBadge from './InfoBadge';
import SaveMenu from './SaveMenu';
import './PcaKdePage.css';

// The main viewport exports as PNG at native or 3× resolution, matching the
// KDE / 3D page's save options.
const SAVE_OPTIONS = [
    { id: 'png', label: 'Standard PNG', hint: '1×' },
    { id: 'png3x', label: 'High quality PNG', hint: '3×' }
];

const GRID_OPTIONS = [24, 32, 40, 48, 56, 64];
const BW_OPTIONS = [
    { value: 'scott', label: 'Scott' },
    { value: 'silverman', label: 'Silverman' }
];
const DEFAULTS = { grid: 40, bw: 'scott', extent: 4, probability: 0.5, isoPercent: 25, colormap: 'viridis', shellColormap: 'viridis', shellContrast: 1, clusterThreshold: 1.5 };

// Contrast highlight for the selected atom in the unit-cell view — a warm color
// that stands out against the cool element palette (teal Se, blue Ga/Ta).
const SELECTED_ATOM_COLOR = 0xff7a1a;

const numberFormat = (value, digits = 4) =>
    Number.isFinite(value) ? value.toFixed(digits) : '—';

// Cartesian point for a PCA-frame grid vertex (fractional grid indices) given the
// engine's per-axis coordinate arrays, principal axes (rows), and cloud mean.
const pcaVertexToCartesian = (fi, fj, fk, axisCoords, axes, mean, out) => {
    const px = axisCoords[0][0] + fi * (axisCoords[0][1] - axisCoords[0][0]);
    const py = axisCoords[1][0] + fj * (axisCoords[1][1] - axisCoords[1][0]);
    const pz = axisCoords[2][0] + fk * (axisCoords[2][1] - axisCoords[2][0]);
    out[0] = mean[0] + px * axes[0][0] + py * axes[1][0] + pz * axes[2][0];
    out[1] = mean[1] + px * axes[0][1] + py * axes[1][1] + pz * axes[2][1];
    out[2] = mean[2] + px * axes[0][2] + py * axes[1][2] + pz * axes[2][2];
    return out;
};

// Trilinear sample of the flat KDE density grid (C order over PC1, PC2, PC3:
// index = (i*grid + j)*grid + k) at continuous grid indices (fi, fj, fk).
const sampleDensityTrilinear = (density, grid, fi, fj, fk) => {
    const clamp = (v) => (v < 0 ? 0 : v > grid - 1 ? grid - 1 : v);
    const ci = clamp(fi);
    const cj = clamp(fj);
    const ck = clamp(fk);
    const i0 = Math.floor(ci);
    const j0 = Math.floor(cj);
    const k0 = Math.floor(ck);
    const i1 = Math.min(i0 + 1, grid - 1);
    const j1 = Math.min(j0 + 1, grid - 1);
    const k1 = Math.min(k0 + 1, grid - 1);
    const di = ci - i0;
    const dj = cj - j0;
    const dk = ck - k0;
    const at = (i, j, k) => density[(i * grid + j) * grid + k];
    const c00 = at(i0, j0, k0) * (1 - di) + at(i1, j0, k0) * di;
    const c01 = at(i0, j0, k1) * (1 - di) + at(i1, j0, k1) * di;
    const c10 = at(i0, j1, k0) * (1 - di) + at(i1, j1, k0) * di;
    const c11 = at(i0, j1, k1) * (1 - di) + at(i1, j1, k1) * di;
    const c0 = c00 * (1 - dj) + c10 * dj;
    const c1 = c01 * (1 - dj) + c11 * dj;
    return c0 * (1 - dk) + c1 * dk;
};

// Solid p% ellipsoid whose surface is colored by the KDE density sampled at each
// vertex -- "projecting" the density onto the harmonic reference shell. A perfectly
// Gaussian cloud gives a near-uniform color; hotter/colder patches mark where the
// real density departs from the ellipsoid (anharmonicity). The ellipsoid-surface
// point for a unit-sphere vertex is (semi .* vertex) in the PCA frame, which both
// maps to world coordinates (baked in, like the isosurface) and indexes the grid.
// Coloring is stretched to the shell's OWN density range (not the global 0..vmax),
// so the small variation across an iso-probability shell reads with full contrast.
const makeEllipsoidKdeSurface = (semi, axes, mean, kde, colormap, contrast = 1) => {
    const geometry = new THREE.SphereGeometry(1, 96, 64);
    const unit = geometry.attributes.position;
    const count = unit.count;
    const world = new Float32Array(count * 3);
    const colors = new Float32Array(count * 3);
    const values = new Float64Array(count);
    const grid = kde.grid;
    const axisCoords = kde.axisCoords;
    const step = [0, 1, 2].map((a) => (axisCoords[a][1] - axisCoords[a][0]) || 1);
    // First pass: bake world positions and sample the density at each vertex,
    // tracking the min/max that actually occur on this shell.
    let shellMin = Infinity;
    let shellMax = -Infinity;
    for (let v = 0; v < count; v += 1) {
        const p0 = semi[0] * unit.getX(v);
        const p1 = semi[1] * unit.getY(v);
        const p2 = semi[2] * unit.getZ(v);
        world[3 * v] = mean[0] + p0 * axes[0][0] + p1 * axes[1][0] + p2 * axes[2][0];
        world[3 * v + 1] = mean[1] + p0 * axes[0][1] + p1 * axes[1][1] + p2 * axes[2][1];
        world[3 * v + 2] = mean[2] + p0 * axes[0][2] + p1 * axes[1][2] + p2 * axes[2][2];
        const d = sampleDensityTrilinear(
            kde.density, grid,
            (p0 - axisCoords[0][0]) / step[0],
            (p1 - axisCoords[1][0]) / step[1],
            (p2 - axisCoords[2][0]) / step[2]
        );
        values[v] = d;
        if (d < shellMin) shellMin = d;
        if (d > shellMax) shellMax = d;
    }
    // Second pass: stretch the shell's density range across the full colormap, then
    // apply the user contrast as a symmetric gain about the mid-tone (0.5) — >1
    // narrows the effective vmin/vmax so faint departures stand out, <1 flattens it.
    const shellRange = shellMax - shellMin || 1;
    for (let v = 0; v < count; v += 1) {
        const t = 0.5 + ((values[v] - shellMin) / shellRange - 0.5) * contrast;
        const rgb = sampleColormap(colormap, t < 0 ? 0 : t > 1 ? 1 : t);
        colors[3 * v] = rgb[0] / 255;
        colors[3 * v + 1] = rgb[1] / 255;
        colors[3 * v + 2] = rgb[2] / 255;
    }
    geometry.setAttribute('position', new THREE.BufferAttribute(world, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    geometry.computeVertexNormals();
    const material = new THREE.MeshBasicMaterial({
        vertexColors: true,
        transparent: true,
        opacity: 0.92,
        side: THREE.FrontSide
    });
    return new THREE.Mesh(geometry, material);
};

// Colormap texture for a 2D KDE projection. density is [gridFirst][gridSecond]
// over the projection's two principal axes; the texel grid is laid out so u runs
// along the first axis and v along the second, matching the wall plane's basis.
const projectionTexture = (projection, colormap) => {
    const density = projection.density;
    const nFirst = density.length;
    const nSecond = density[0] ? density[0].length : 0;
    if (!nFirst || !nSecond) return null;
    const canvas = document.createElement('canvas');
    canvas.width = nFirst;
    canvas.height = nSecond;
    const ctx = canvas.getContext('2d');
    const image = ctx.createImageData(nFirst, nSecond);
    const lut = getLut(colormap);
    const stops = lut.length / 3 - 1;
    const vmax = projection.vmax || 1;
    for (let i = 0; i < nFirst; i += 1) {
        for (let j = 0; j < nSecond; j += 1) {
            const value = Math.max(0, Math.min(1, density[i][j] / vmax));
            const lutIndex = Math.round(value * stops) * 3;
            // CanvasTexture flips vertically on upload, so row (nSecond-1-j) maps
            // to v = j: second axis then increases along the plane's local +Y.
            const pixel = ((nSecond - 1 - j) * nFirst + i) * 4;
            image.data[pixel] = lut[lutIndex];
            image.data[pixel + 1] = lut[lutIndex + 1];
            image.data[pixel + 2] = lut[lutIndex + 2];
            image.data[pixel + 3] = 255;
        }
    }
    ctx.putImageData(image, 0, 0);
    const texture = new THREE.CanvasTexture(canvas);
    texture.colorSpace = THREE.SRGBColorSpace;
    texture.needsUpdate = true;
    return texture;
};

// Build the wall plane for one projection: a textured quad in the PCA frame on
// the far wall along the axis not spanned by the projection (Maksim Eremenko's
// shadow-box layout). first/second index the projection's axes; the plane's
// local X → axes[first], local Y → axes[second], placed at -halfWidth of the
// remaining axis and offset slightly outward so it sits just past the cloud.
const makeProjectionWall = (projection, axes, mean, halfWidths, colormap) => {
    const [first, second] = projection.axes;
    const third = 3 - first - second;
    const texture = projectionTexture(projection, colormap);
    if (!texture) return null;

    const geometry = new THREE.PlaneGeometry(2 * halfWidths[first], 2 * halfWidths[second]);
    const material = new THREE.MeshBasicMaterial({
        map: texture,
        side: THREE.DoubleSide,
        transparent: true,
        opacity: 0.96,
        depthWrite: false
    });
    const plane = new THREE.Mesh(geometry, material);

    const xAxis = new THREE.Vector3(axes[first][0], axes[first][1], axes[first][2]);
    const yAxis = new THREE.Vector3(axes[second][0], axes[second][1], axes[second][2]);
    const normal = new THREE.Vector3().crossVectors(xAxis, yAxis).normalize();
    const offset = -(halfWidths[third] + 0.06 * halfWidths[third]);
    const position = new THREE.Vector3(
        mean[0] + offset * axes[third][0],
        mean[1] + offset * axes[third][1],
        mean[2] + offset * axes[third][2]
    );
    plane.applyMatrix4(new THREE.Matrix4().makeBasis(xAxis, yAxis, normal).setPosition(position));
    return plane;
};

// Marching-squares line segments of a 2D field at one level. Each cell yields 0,
// 1, or 2 segments; endpoints come back in continuous grid coordinates (fi along
// the first axis, fj along the second).
const contourSegments = (density, level) => {
    const nFirst = density.length;
    const nSecond = density[0] ? density[0].length : 0;
    const segments = [];
    for (let i = 0; i < nFirst - 1; i += 1) {
        for (let j = 0; j < nSecond - 1; j += 1) {
            const corners = [
                { fi: i, fj: j, v: density[i][j] },
                { fi: i + 1, fj: j, v: density[i + 1][j] },
                { fi: i + 1, fj: j + 1, v: density[i + 1][j + 1] },
                { fi: i, fj: j + 1, v: density[i][j + 1] }
            ];
            const crossings = [];
            for (const [a, b] of [[0, 1], [1, 2], [2, 3], [3, 0]]) {
                const ca = corners[a];
                const cb = corners[b];
                if ((ca.v < level) !== (cb.v < level)) {
                    const t = (level - ca.v) / (cb.v - ca.v);
                    crossings.push([ca.fi + (cb.fi - ca.fi) * t, ca.fj + (cb.fj - ca.fj) * t]);
                }
            }
            if (crossings.length === 2) {
                segments.push(crossings[0], crossings[1]);
            } else if (crossings.length === 4) {
                segments.push(crossings[0], crossings[1], crossings[2], crossings[3]);
            }
        }
    }
    return segments;
};

// Contour lines for one projection, drawn just in front of its wall so they read
// on top of the colormap (like the projected contours in Maksim's Plotly figure).
const makeProjectionContours = (projection, axes, mean, halfWidths, color) => {
    const density = projection.density;
    const nFirst = density.length;
    const nSecond = density[0] ? density[0].length : 0;
    const vmax = projection.vmax || 0;
    if (!nFirst || !nSecond || vmax <= 0) return null;

    const [first, second] = projection.axes;
    const third = 3 - first - second;
    const stepFirst = (2 * halfWidths[first]) / Math.max(nFirst - 1, 1);
    const stepSecond = (2 * halfWidths[second]) / Math.max(nSecond - 1, 1);
    // Sit the lines just inside the wall (toward the box interior).
    const wallOffset = -(halfWidths[third] + 0.05 * halfWidths[third]);

    const points = [];
    [0.1, 0.25, 0.4, 0.55, 0.7, 0.85].forEach((fraction) => {
        contourSegments(density, fraction * vmax).forEach(([fi, fj]) => {
            const pFirst = -halfWidths[first] + fi * stepFirst;
            const pSecond = -halfWidths[second] + fj * stepSecond;
            points.push(new THREE.Vector3(
                mean[0] + pFirst * axes[first][0] + pSecond * axes[second][0] + wallOffset * axes[third][0],
                mean[1] + pFirst * axes[first][1] + pSecond * axes[second][1] + wallOffset * axes[third][1],
                mean[2] + pFirst * axes[first][2] + pSecond * axes[second][2] + wallOffset * axes[third][2]
            ));
        });
    });
    if (!points.length) return null;
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    const material = new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.55, depthWrite: false });
    const lines = new THREE.LineSegments(geometry, material);
    // The contour lines and their wall are transparent and nearly coplanar, so a
    // plain distance sort makes the wall randomly draw over the lines and hide
    // them. renderOrder forces the lines to paint after the walls; depthTest stays
    // on, so the central isosurface/ellipsoid still occlude them correctly.
    lines.renderOrder = 2;
    return lines;
};

// Wireframe box around the sampling volume (±halfWidths in the PCA frame), so
// the wall projections read as the faces of a shadow box.
const makeBoundingBox = (axes, mean, halfWidths, color) => {
    const corners = [];
    for (const sx of [-1, 1]) {
        for (const sy of [-1, 1]) {
            for (const sz of [-1, 1]) {
                const p = [sx * halfWidths[0], sy * halfWidths[1], sz * halfWidths[2]];
                corners.push(new THREE.Vector3(
                    mean[0] + p[0] * axes[0][0] + p[1] * axes[1][0] + p[2] * axes[2][0],
                    mean[1] + p[0] * axes[0][1] + p[1] * axes[1][1] + p[2] * axes[2][1],
                    mean[2] + p[0] * axes[0][2] + p[1] * axes[1][2] + p[2] * axes[2][2]
                ));
            }
        }
    }
    // Corner index = ((sx>0)<<2) | ((sy>0)<<1) | (sz>0); connect Hamming-1 pairs.
    const edges = [];
    for (let a = 0; a < 8; a += 1) {
        for (let b = a + 1; b < 8; b += 1) {
            const diff = a ^ b;
            if (diff === 1 || diff === 2 || diff === 4) edges.push(corners[a], corners[b]);
        }
    }
    const geometry = new THREE.BufferGeometry().setFromPoints(edges);
    return new THREE.LineSegments(geometry, new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.35 }));
};

const TRIAD_COLORS = [0xd64545, 0x3fa34d, 0x3f7fd6]; // PC1 red, PC2 green, PC3 blue
// CSS twins of TRIAD_COLORS, used to tint the PC rows in the statistics tables so
// they read as the same axes shown by the 3D triad and the viewport legend.
const PC_CSS_COLORS = ['#d64545', '#3fa34d', '#3f7fd6'];

// Wireframe color for the p% thermal-ellipsoid reference cage, chosen from the
// controls bar. The cage is always drawn transparent, so whatever color is picked
// it never crowds the KDE shell / isosurface it wraps. The value feeds the 3D
// material and the viewport legend swatch.
const ELLIPSOID_COLOR_OPTIONS = [
    { value: '#b892ff', label: 'Violet' },
    { value: '#ff7a1a', label: 'Amber' },
    { value: '#ffffff', label: 'White' },
    { value: '#111111', label: 'Black' },
    { value: '#c3d4e6', label: 'Silver' },
    { value: '#8fd4ef', label: 'Cyan' }
];
// Violet reads well against the default (viridis) KDE shell.
const ELLIPSOID_COLOR = ELLIPSOID_COLOR_OPTIONS[0].value;
const ELLIPSOID_OPACITY = 0.4;

const TRIAD_UP = new THREE.Vector3(0, 1, 0);

// PC1/PC2/PC3 triad as thin rods from `origin` along the principal axes (rows).
// Thin cylinders rather than lines so a little thickness reads at any zoom
// without dominating; used both around the volume and on the selected atom.
const buildAxisTriad = (origin, axes, length, radius) => {
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

const PROJECTION_META = [
    { key: 'pc12', label: 'PC1 – PC2' },
    { key: 'pc13', label: 'PC1 – PC3' },
    { key: 'pc23', label: 'PC2 – PC3' }
];

// Vertical field of view (deg) of the perspective camera. The orthographic
// camera's frustum is sized from it so switching projection preserves framing.
const CAMERA_FOV = 45;

// Size an orthographic camera so its view at `distance` matches a perspective
// camera of CAMERA_FOV at the same distance (same on-screen scale of the box).
const setOrthoFrustum = (camera, distance, aspect) => {
    const halfHeight = distance * Math.tan((CAMERA_FOV * Math.PI) / 360);
    const halfWidth = halfHeight * aspect;
    camera.left = -halfWidth;
    camera.right = halfWidth;
    camera.top = halfHeight;
    camera.bottom = -halfHeight;
};

// Place the main camera at `radius` along `dir` from the cloud mean, looking
// back at it, with `up` as the screen-up axis. Handles perspective and
// orthographic cameras alike (the latter gets a matched frustum + reset zoom).
const placeMainCamera = (camera, controls, kde, dir, up, aspect) => {
    const meanVec = new THREE.Vector3(kde.mean[0], kde.mean[1], kde.mean[2]);
    const axisLength = Math.max(...kde.halfWidths);
    const radius = axisLength * 4.3 || 1;
    // OrbitControls carries a damped orbit velocity after a drag, which its own
    // update() loop normally decays to zero over ~1 s. Flush it FIRST, with a
    // damping-off update() (that both applies and then zeroes the pending delta),
    // so none of that leftover rotation gets applied to the pose we are about to
    // set. Without this, a reframe fired before the glide settles — e.g. right
    // after loading a new run — lands rotated off the intended default.
    const damping = controls.enableDamping;
    controls.enableDamping = false;
    controls.update();
    // Now place the camera on the requested framing and commit it with zero velocity.
    controls.target.copy(meanVec);
    camera.up.set(up[0], up[1], up[2]);
    camera.position.copy(meanVec).addScaledVector(dir, radius);
    camera.near = radius / 100;
    camera.far = radius * 100;
    if (camera.isOrthographicCamera) {
        camera.zoom = 1;
        setOrthoFrustum(camera, radius, aspect);
    }
    camera.updateProjectionMatrix();
    controls.update();
    controls.enableDamping = damping;
};

// Point the main camera down the box body diagonal (+PC1/+PC2/+PC3) at a
// comfortable distance, so the three min-corner walls sit symmetrically at the
// back and the whole cube fits corner-on. This is the default framing, applied
// on a new site and re-applied by the panel's "Reset view" button.
const frameMainCamera = (camera, controls, kde, aspect = 1) => {
    const axes = kde.axes;
    const dir = new THREE.Vector3(
        axes[0][0] + axes[1][0] + axes[2][0],
        axes[0][1] + axes[1][1] + axes[2][1],
        axes[0][2] + axes[1][2] + axes[2][2]
    ).normalize();
    placeMainCamera(camera, controls, kde, dir, axes[2], aspect);
};

// Screen-up axis when looking straight down each principal axis: along PC1 or
// PC2, PC3 points up; along PC3, PC2 points up.
const VIEW_UP_AXIS = [2, 2, 1];

// Look straight down principal axis `axisIndex`, so its conjugate PC-plane faces
// the camera (the projection wall for that pair fills the frame).
const frameAlongAxis = (camera, controls, kde, axisIndex, aspect = 1) => {
    const axes = kde.axes;
    const dir = new THREE.Vector3(
        axes[axisIndex][0], axes[axisIndex][1], axes[axisIndex][2]
    ).normalize();
    placeMainCamera(camera, controls, kde, dir, axes[VIEW_UP_AXIS[axisIndex]], aspect);
};

export default function PcaKdePage({ directory, localRun, onSitesChange }) {
    // The loaded .rmc6f text, tagged with the file it came from, so a just-changed
    // dataset never runs against the previous model's text (see rmc6fText below).
    const [loadedText, setLoadedText] = useState({ file: null, text: null });
    const [sites, setSites] = useState(null);
    const [selectedRef, setSelectedRef] = useState(null);
    const [kde, setKde] = useState(null);
    const [sitesError, setSitesError] = useState(null);
    const [kdeError, setKdeError] = useState(null);
    const [loadingSites, setLoadingSites] = useState(false);
    const [loadingKde, setLoadingKde] = useState(false);

    const [grid, setGrid] = useState(DEFAULTS.grid);
    const [bw, setBw] = useState(DEFAULTS.bw);
    const [extent, setExtent] = useState(DEFAULTS.extent);
    // Fold-and-cluster distance (Å) used only when the loaded file has no reference-
    // site/cell columns and its sites must be reconstructed (see the PCA worker).
    const [clusterThreshold, setClusterThreshold] = useState(DEFAULTS.clusterThreshold);
    const [probability, setProbability] = useState(DEFAULTS.probability);
    const [isoPercent, setIsoPercent] = useState(DEFAULTS.isoPercent);
    const [colormap, setColormap] = useState(DEFAULTS.colormap);
    // The KDE shell has its own colormap so it can be read independently of the
    // wall projections / isosurface (which use `colormap`).
    const [shellColormap, setShellColormap] = useState(DEFAULTS.shellColormap);
    // Contrast gain for the shell colormap: a single knob that tightens (>1) or
    // loosens (<1) the effective vmin/vmax around the shell's mid-tone.
    const [shellContrast, setShellContrast] = useState(DEFAULTS.shellContrast);
    const [ellipsoidColor, setEllipsoidColor] = useState(ELLIPSOID_COLOR);
    const [showEllipsoid, setShowEllipsoid] = useState(true);
    // The KDE shell is the default density view; it and the isosurface are mutually
    // exclusive (the shell shows the same density from the outside), so the
    // isosurface starts off.
    const [showSurface, setShowSurface] = useState(false);
    const [showProjections, setShowProjections] = useState(true);
    const [showEllipsoidKde, setShowEllipsoidKde] = useState(true);
    const [perspective, setPerspective] = useState(true);
    // Which principal axis the camera is snapped to look down (null = free / the
    // default diagonal view); drives the highlight on the view-axis buttons.
    const [viewAxis, setViewAxis] = useState(null);

    const workerRef = useRef(null);
    const requestIdRef = useRef(0);
    const mountRef = useRef(null);
    const sceneRef = useRef(null);
    const framedRefRef = useRef(null);
    const framedDatasetRef = useRef(null);
    const structureMountRef = useRef(null);
    const structureSceneRef = useRef(null);
    const selectRef = useRef(null);

    const staticMode = isStaticMode();
    const structureFile = localRun?.structureFile || null;
    // A locally-loaded run (the Demo, or a picked folder) carries its .rmc6f as a
    // browser file, so it is parsed in the worker in BOTH runtimes; only a typed
    // backend directory goes through the Flask API. Static mode has no backend, so
    // it always relies on a local file.
    const localFile = structureFile?.sourceFile || null;
    // Only treat the text as current when it was read from the file in props right
    // now. On a dataset switch this is null until the new file's text loads, so the
    // effects below never fire against the previous model.
    const rmc6fText = loadedText.file === localFile ? loadedText.text : null;
    // Stable per-dataset identity: App's monotonic runId (bumped on every folder /
    // demo load, but KEPT across Live Data polls of the same run), falling back to
    // the backend directory. A change here means a genuinely different dataset was
    // loaded — a live refresh of the same run leaves it untouched.
    const datasetKey = localRun?.runId ?? directory ?? null;

    // --- Read the raw .rmc6f text once per local file. ------------------------
    useEffect(() => {
        let cancelled = false;
        if (localFile) {
            localFile.text().then((text) => {
                if (!cancelled) setLoadedText({ file: localFile, text });
            }).catch(() => {
                if (!cancelled) { setLoadedText({ file: null, text: null }); setSitesError('Could not read the structure file.'); }
            });
        } else {
            setLoadedText({ file: null, text: null });
        }
        return () => { cancelled = true; };
    }, [localFile]);

    // --- Worker lifecycle (used whenever a local file is the data source). -----
    useEffect(() => {
        const worker = new Worker(new URL('../workers/pcaKdeWorker.js', import.meta.url), { type: 'module' });
        workerRef.current = worker;
        return () => { worker.terminate(); workerRef.current = null; };
    }, []);

    // A single request path for both data sources: the worker when a local file
    // is loaded, otherwise the Flask API against the run directory.
    const requestPca = useCallback((kind, params) => {
        if (rmc6fText) {
            return new Promise((resolve, reject) => {
                const worker = workerRef.current;
                if (!worker) { reject(new Error('PCA-KDE worker unavailable')); return; }
                const id = requestIdRef.current + 1;
                requestIdRef.current = id;
                const handler = (event) => {
                    if (event.data.id !== id) return;
                    worker.removeEventListener('message', handler);
                    if (event.data.error) reject(new Error(event.data.error));
                    else resolve(event.data.result);
                };
                worker.addEventListener('message', handler);
                worker.postMessage({ id, kind, text: rmc6fText, ...params });
            });
        }
        const endpoint = kind === 'sites' ? '/api/pca/sites' : '/api/pca/kde';
        return axios
            .get(`${API_BASE_URL}${endpoint}`, { params: { dir: directory || '.', ...params } })
            .then((response) => response.data);
    }, [rmc6fText, directory]);

    // --- Load the per-site ellipsoid table. -----------------------------------
    useEffect(() => {
        let cancelled = false;
        const loadSites = async () => {
            // Wait for a local file's text before requesting; a backend directory
            // needs no text and proceeds immediately.
            if (localFile && !rmc6fText) { setSites(null); return; }
            setLoadingSites(true);
            setSitesError(null);
            try {
                const data = await requestPca('sites', { probability, clusterThreshold });
                if (cancelled) return;
                setSites(data);
                setSelectedRef((current) => {
                    const refs = data.sites.map((site) => site.referenceNumber);
                    if (current != null && refs.includes(current)) return current;
                    // Prefer a clean reconstructed site (members == one-per-cell) so the
                    // default view is a genuine thermal ellipsoid, not a disordered shell.
                    // Inert for files with real site columns (copiesPerCell is null there).
                    const clean = data.sites.find((site) => site.copiesPerCell && site.count === site.copiesPerCell);
                    return clean?.referenceNumber ?? refs[0] ?? null;
                });
            } catch (error) {
                if (!cancelled) { setSites(null); setSitesError(error.message); }
            } finally {
                if (!cancelled) setLoadingSites(false);
            }
        };
        loadSites();
        return () => { cancelled = true; };
    }, [requestPca, localFile, rmc6fText, probability, clusterThreshold]);

    // --- Load the KDE volume for the selected site. ---------------------------
    useEffect(() => {
        let cancelled = false;
        const loadKde = async () => {
            if (selectedRef == null) { setKde(null); return; }
            if (localFile && !rmc6fText) return;
            setLoadingKde(true);
            setKdeError(null);
            try {
                const data = await requestPca(
                    'kde',
                    // cubicBox keeps all three axes on one scale, so the wall
                    // projections come out the same size (Maksim's cubic layout).
                    { referenceNumber: selectedRef, grid, bw, extent, probability, clusterThreshold, cubicBox: true, projections: true }
                );
                if (cancelled) return;
                // Silent view reset on a genuinely new run: null the framing marker
                // right before storing THIS fresh volume, so the rebuild effect
                // reframes to the default against it — never a stale one — even when
                // the new run's first site shares a reference number with the old.
                // Site and slider changes keep datasetKey (and a Live Data refresh
                // reuses the same runId), so those never trigger this reset.
                if (framedDatasetRef.current !== datasetKey) {
                    framedDatasetRef.current = datasetKey;
                    framedRefRef.current = null;
                }
                setKde(data);
            } catch (error) {
                if (cancelled) return;
                // Changing the cluster distance re-clusters and re-numbers the sites, so
                // the KDE effect can momentarily request a reference number the new
                // clustering dropped. That self-heals when the reloaded sites re-pick a
                // valid site, so keep the current volume and don't surface the transient.
                if (/Unknown reference number/i.test(error.message || '')) return;
                setKde(null); setKdeError(error.message);
            } finally {
                if (!cancelled) setLoadingKde(false);
            }
        };
        loadKde();
        return () => { cancelled = true; };
    }, [requestPca, localFile, rmc6fText, selectedRef, grid, bw, extent, probability, clusterThreshold, datasetKey]);

    const selectedEllipsoid = useMemo(
        () => sites?.sites.find((site) => site.referenceNumber === selectedRef) || null,
        [sites, selectedRef]
    );

    const elementColors = useMemo(
        () => buildElementColors(sites?.elements ?? []),
        [sites]
    );

    // Current aspect ratio of the main canvas, needed to size the orthographic
    // frustum when (re)framing.
    const mainAspect = () => {
        const mount = mountRef.current;
        return mount && mount.clientHeight ? mount.clientWidth / mount.clientHeight : 1;
    };

    // Re-apply the default (diagonal) framing to the main panel on demand (the
    // user may have orbited/zoomed away). No-op until a volume is loaded.
    const resetMainView = useCallback(() => {
        const handle = sceneRef.current;
        if (handle && kde) {
            frameMainCamera(handle.camera, handle.controls, kde, mainAspect());
            setViewAxis(null);
        }
    }, [kde]);

    // Snap the camera to look straight down a principal axis (PC1/PC2/PC3).
    const lookAlong = useCallback((axisIndex) => {
        const handle = sceneRef.current;
        if (!handle || !kde) return;
        frameAlongAxis(handle.camera, handle.controls, kde, axisIndex, mainAspect());
        setViewAxis(axisIndex);
    }, [kde]);

    // Toggle the main viewport between perspective and orthographic projection,
    // preserving the current orientation (the scene handle swaps the camera).
    const applyPerspective = useCallback((value) => {
        setPerspective(value);
        sceneRef.current?.setProjection?.(value);
    }, []);

    // Export the current main-panel figure as PNG. Standard reads the live canvas
    // (preserveDrawingBuffer keeps it readable); high quality re-renders the same
    // frame at 3× pixel ratio, captures, then restores.
    const saveMainView = useCallback(async (format) => {
        const handle = sceneRef.current;
        if (!handle) return;
        const { renderer, scene, camera } = handle;
        const name = selectedEllipsoid
            ? `PCA_Ellipsoid_${selectedEllipsoid.element}_site${selectedEllipsoid.referenceNumber}`
            : 'PCA_Ellipsoid';
        if (format === 'png3x') {
            const size = renderer.getSize(new THREE.Vector2());
            const previousRatio = renderer.getPixelRatio();
            renderer.setPixelRatio(3);
            renderer.setSize(size.x, size.y, false);
            renderer.render(scene, camera);
            const blob = await new Promise((resolve) => renderer.domElement.toBlob(resolve, 'image/png'));
            renderer.setPixelRatio(previousRatio);
            renderer.setSize(size.x, size.y, false);
            renderer.render(scene, camera);
            if (blob) downloadBlob(blob, `${sanitizeFilename(name)}.png`);
        } else {
            renderer.render(scene, camera);
            await saveCanvasAsPng(renderer.domElement, name);
        }
    }, [selectedEllipsoid]);

    // Publish the per-site ellipsoid table upward (App → AI Assistant), so the
    // assistant can reason about the thermal displacements once this page has
    // computed them. Follows the current dataset, and clears on error/reset.
    useEffect(() => {
        onSitesChange?.(sites);
    }, [sites, onSitesChange]);

    // Keep the click handler pointing at the current setter without re-creating
    // the (once-mounted) structure scene.
    selectRef.current = setSelectedRef;

    // --- Three.js scene: build once, keep a handle for updates. ---------------
    useEffect(() => {
        const mount = mountRef.current;
        if (!mount) return undefined;
        const width = mount.clientWidth || 640;
        const height = mount.clientHeight || 520;

        const scene = new THREE.Scene();
        // Two cameras share the scene: a perspective one (default) and an
        // orthographic one (parallel projection, no foreshortening — useful for
        // comparing the KDE isosurface against the harmonic ellipsoid). The
        // projection toggle swaps which is active; `camera`/`controls` are `let`
        // so the animate loop and handle always read the live pair.
        const perspectiveCamera = new THREE.PerspectiveCamera(CAMERA_FOV, width / height, 0.001, 1000);
        perspectiveCamera.position.set(0.6, 0.5, 0.9);
        const orthographicCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.001, 1000);
        let camera = perspectiveCamera;
        // preserveDrawingBuffer keeps the rendered frame readable so the figure
        // can be captured for PNG export at any time.
        const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, preserveDrawingBuffer: true });
        renderer.setPixelRatio(window.devicePixelRatio || 1);
        renderer.setSize(width, height);
        mount.appendChild(renderer.domElement);

        let controls = new OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.dampingFactor = 0.12;

        scene.add(new THREE.AmbientLight(0xffffff, 0.65));
        const key = new THREE.DirectionalLight(0xffffff, 0.8);
        key.position.set(1, 1, 1);
        scene.add(key);
        const fill = new THREE.DirectionalLight(0xffffff, 0.35);
        fill.position.set(-1, -0.5, -1);
        scene.add(fill);

        const surfaceGroup = new THREE.Group();
        const ellipsoidGroup = new THREE.Group();
        const axesGroup = new THREE.Group();
        const wallsGroup = new THREE.Group();
        scene.add(wallsGroup, surfaceGroup, ellipsoidGroup, axesGroup);

        let animationId = 0;
        const animate = () => {
            controls.update();
            renderer.render(scene, camera);
            animationId = requestAnimationFrame(animate);
        };
        animate();

        const handleResize = () => {
            const w = mount.clientWidth || width;
            const h = mount.clientHeight || height;
            if (w === 0 || h === 0) return;
            const aspect = w / h;
            perspectiveCamera.aspect = aspect;
            perspectiveCamera.updateProjectionMatrix();
            // Keep the ortho camera's vertical extent; recompute the horizontal
            // half-width from the new aspect so a resize doesn't stretch the view.
            const halfHeight = orthographicCamera.top || 1;
            orthographicCamera.left = -halfHeight * aspect;
            orthographicCamera.right = halfHeight * aspect;
            orthographicCamera.updateProjectionMatrix();
            renderer.setSize(w, h);
        };
        window.addEventListener('resize', handleResize);
        // Track the mount itself so the canvas follows the panel's aspect-ratio
        // box and the responsive column layout, not just window resizes.
        const resizeObserver = new ResizeObserver(handleResize);
        resizeObserver.observe(mount);

        const handle = { scene, camera, renderer, controls, surfaceGroup, ellipsoidGroup, axesGroup, wallsGroup };

        // Swap the active camera between perspective and orthographic while keeping
        // the current orientation, distance, and pivot, so the toggle looks like a
        // pure projection change. Rebuilds OrbitControls on the new camera so its
        // up-vector bookkeeping is correct for the swapped type.
        handle.setProjection = (usePerspective) => {
            const next = usePerspective ? perspectiveCamera : orthographicCamera;
            if (next === camera) return;
            const target = controls.target.clone();
            const distance = camera.position.distanceTo(target);
            const w = mount.clientWidth || width;
            const h = mount.clientHeight || height;
            const aspect = h ? w / h : 1;
            next.up.copy(camera.up);
            next.position.copy(camera.position);
            next.near = camera.near;
            next.far = camera.far;
            if (next.isOrthographicCamera) {
                next.zoom = 1;
                setOrthoFrustum(next, distance, aspect);
            } else {
                next.aspect = aspect;
            }
            next.updateProjectionMatrix();
            next.lookAt(target);
            const previousControls = controls;
            const nextControls = new OrbitControls(next, renderer.domElement);
            nextControls.enableDamping = true;
            nextControls.dampingFactor = 0.12;
            nextControls.target.copy(target);
            previousControls.dispose();
            nextControls.update();
            camera = next;
            controls = nextControls;
            handle.camera = camera;
            handle.controls = controls;
        };

        sceneRef.current = handle;

        return () => {
            cancelAnimationFrame(animationId);
            window.removeEventListener('resize', handleResize);
            resizeObserver.disconnect();
            controls.dispose();
            renderer.dispose();
            if (renderer.domElement.parentNode === mount) mount.removeChild(renderer.domElement);
            sceneRef.current = null;
        };
    }, []);

    // --- Rebuild the isosurface, ellipsoid, and axis triad when data changes. --
    useEffect(() => {
        const handle = sceneRef.current;
        if (!handle || !kde) return;
        const { surfaceGroup, ellipsoidGroup, axesGroup, wallsGroup, camera, controls } = handle;

        const dispose = (group) => {
            while (group.children.length) {
                const child = group.children.pop();
                child.geometry?.dispose();
                // Wall planes carry a CanvasTexture that must be released too.
                child.material?.map?.dispose();
                child.material?.dispose();
                group.remove(child);
            }
        };
        dispose(surfaceGroup);
        dispose(ellipsoidGroup);
        dispose(axesGroup);
        dispose(wallsGroup);

        const axes = kde.axes;
        const mean = kde.mean;
        const axisCoords = kde.axisCoords;
        const gridN = kde.grid;
        const halfWidths = kde.halfWidths;

        // Density projections on the far walls + a wireframe cage (Maksim
        // Eremenko's shadow-box layout): the isosurface floats inside a box whose
        // three back walls show the PC-plane KDE projections.
        if (showProjections && kde.projections) {
            wallsGroup.add(makeBoundingBox(axes, mean, halfWidths, 0x8a97a8));
            PROJECTION_META.forEach(({ key }) => {
                const projection = kde.projections[key];
                const wall = makeProjectionWall(projection, axes, mean, halfWidths, colormap);
                if (wall) wallsGroup.add(wall);
                const contours = makeProjectionContours(projection, axes, mean, halfWidths, 0xeef3f8);
                if (contours) wallsGroup.add(contours);
            });
        }

        // Isosurface at the enclosed-mass threshold from the slider.
        const massLevel = kde.massLevels?.[isoPercent]?.level ?? (kde.vmax * (1 - isoPercent / 100));
        if (showSurface && massLevel > 0) {
            const field = kde.density instanceof Float64Array || Array.isArray(kde.density)
                ? kde.density
                : Float64Array.from(kde.density);
            const surface = marchingCubes(field, gridN, gridN, gridN, massLevel);
            if (surface.count > 0) {
                const cartesian = new Float32Array(surface.positions.length);
                const tmp = [0, 0, 0];
                for (let v = 0; v < surface.count; v += 1) {
                    pcaVertexToCartesian(
                        surface.positions[3 * v], surface.positions[3 * v + 1], surface.positions[3 * v + 2],
                        axisCoords, axes, mean, tmp
                    );
                    cartesian[3 * v] = tmp[0];
                    cartesian[3 * v + 1] = tmp[1];
                    cartesian[3 * v + 2] = tmp[2];
                }
                const geometry = new THREE.BufferGeometry();
                geometry.setAttribute('position', new THREE.BufferAttribute(cartesian, 3));
                geometry.computeVertexNormals();
                const rgb = sampleColormap(colormap, 0.72);
                const material = new THREE.MeshPhongMaterial({
                    color: new THREE.Color(rgb[0] / 255, rgb[1] / 255, rgb[2] / 255),
                    transparent: true,
                    opacity: 0.55,
                    side: THREE.DoubleSide,
                    shininess: 40
                });
                surfaceGroup.add(new THREE.Mesh(geometry, material));
            }
        }

        // Optional: paint the KDE density onto the p% ellipsoid surface. Drawn
        // before the wireframe so, when both are on, the cage sits on top.
        if (showEllipsoidKde && selectedEllipsoid) {
            ellipsoidGroup.add(
                makeEllipsoidKdeSurface(selectedEllipsoid.semiAxes, axes, mean, kde, shellColormap, shellContrast)
            );
        }

        // p% thermal ellipsoid: unit sphere scaled by the semi-axes, oriented by
        // the principal axes. This is the harmonic reference the KDE isosurface
        // is compared against -- where they differ, the motion is anharmonic.
        if (showEllipsoid && selectedEllipsoid) {
            const semi = selectedEllipsoid.semiAxes;
            const sphere = new THREE.SphereGeometry(1, 48, 32);
            const orientation = new THREE.Matrix4().set(
                axes[0][0], axes[1][0], axes[2][0], 0,
                axes[0][1], axes[1][1], axes[2][1], 0,
                axes[0][2], axes[1][2], axes[2][2], 0,
                0, 0, 0, 1
            );
            const scaleMatrix = new THREE.Matrix4().makeScale(semi[0], semi[1], semi[2]);
            const model = orientation.multiply(scaleMatrix);
            model.setPosition(mean[0], mean[1], mean[2]);
            // Always a translucent cage, so it never crowds the KDE shell / isosurface.
            const material = new THREE.MeshBasicMaterial({
                color: new THREE.Color(ellipsoidColor),
                wireframe: true,
                transparent: true,
                opacity: ELLIPSOID_OPACITY
            });
            const mesh = new THREE.Mesh(sphere, material);
            mesh.applyMatrix4(model);
            ellipsoidGroup.add(mesh);
        }

        // Principal-axis triad (PC1 red, PC2 green, PC3 blue) from the centre out
        // to the box faces. Thin rods so they read without dominating.
        const axisLength = Math.max(...kde.halfWidths);
        const meanVec = new THREE.Vector3(mean[0], mean[1], mean[2]);
        buildAxisTriad(meanVec, axes, axisLength, axisLength * 0.012)
            .forEach((rod) => axesGroup.add(rod));

        // Reframe to the default view only when the displayed volume is for a new
        // site or a new run, so sweeping a slider or toggling a layer doesn't yank
        // the view. framedRefRef tracks the site the camera is framed for; loadKde
        // resets it to null on a new run so this fires there too, even when the new
        // run's first site shares a reference number with the old one. frameMainCamera
        // is deterministic (it neutralizes orbit momentum), so it lands cleanly even
        // when this runs while the panel is off-screen after a load.
        const kdeRef = kde.referenceNumber ?? selectedRef;
        if (framedRefRef.current !== kdeRef) {
            const m = mountRef.current;
            frameMainCamera(camera, controls, kde, m && m.clientHeight ? m.clientWidth / m.clientHeight : 1);
            framedRefRef.current = kdeRef;
            // A fresh diagonal framing supersedes any snapped PC-axis view.
            setViewAxis(null);
        } else {
            // Same site: keep the pivot and clip planes synced to the current volume
            // without moving the user's camera.
            const radius = axisLength * 4.3 || 1;
            controls.target.copy(meanVec);
            camera.near = radius / 100;
            camera.far = radius * 100;
            camera.updateProjectionMatrix();
            controls.update();
        }
    }, [kde, isoPercent, colormap, shellColormap, shellContrast, ellipsoidColor, showEllipsoid, showSurface, showProjections, showEllipsoidKde, selectedEllipsoid, selectedRef]);

    // --- Unit-cell structure scene: one clickable marker per reference site. ---
    useEffect(() => {
        const mount = structureMountRef.current;
        if (!mount) return undefined;
        const width = mount.clientWidth || 260;
        const height = mount.clientHeight || 260;

        const scene = new THREE.Scene();
        const camera = new THREE.PerspectiveCamera(40, width / height, 0.01, 1000);
        const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
        renderer.setPixelRatio(window.devicePixelRatio || 1);
        renderer.setSize(width, height);
        mount.appendChild(renderer.domElement);

        const controls = new OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.dampingFactor = 0.12;
        controls.enablePan = false;

        scene.add(new THREE.AmbientLight(0xffffff, 0.7));
        const key = new THREE.DirectionalLight(0xffffff, 0.7);
        key.position.set(1, 1, 1);
        scene.add(key);

        const sitesGroup = new THREE.Group();
        const bondsGroup = new THREE.Group();
        const cellGroup = new THREE.Group();
        const highlightGroup = new THREE.Group();
        scene.add(cellGroup, bondsGroup, sitesGroup, highlightGroup);

        // Click (not drag) on a site marker selects that site. A small pointer
        // travel budget separates a click from an orbit drag.
        const raycaster = new THREE.Raycaster();
        const pointer = new THREE.Vector2();
        let downXY = null;
        const toNdc = (event) => {
            const rect = renderer.domElement.getBoundingClientRect();
            pointer.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
            pointer.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
        };
        const pick = () => {
            raycaster.setFromCamera(pointer, camera);
            return raycaster.intersectObjects(sitesGroup.children, false)[0]?.object ?? null;
        };
        const onPointerDown = (event) => { downXY = [event.clientX, event.clientY]; };
        const onPointerUp = (event) => {
            if (!downXY) return;
            const moved = Math.hypot(event.clientX - downXY[0], event.clientY - downXY[1]);
            downXY = null;
            if (moved > 5) return;
            toNdc(event);
            const hit = pick();
            if (hit) selectRef.current?.(hit.userData.referenceNumber);
        };
        const onPointerMove = (event) => {
            toNdc(event);
            renderer.domElement.style.cursor = pick() ? 'pointer' : 'grab';
        };
        renderer.domElement.addEventListener('pointerdown', onPointerDown);
        renderer.domElement.addEventListener('pointerup', onPointerUp);
        renderer.domElement.addEventListener('pointermove', onPointerMove);

        let animationId = 0;
        const animate = () => {
            controls.update();
            renderer.render(scene, camera);
            animationId = requestAnimationFrame(animate);
        };
        animate();

        const handleResize = () => {
            const w = mount.clientWidth || width;
            const h = mount.clientHeight || height;
            if (w === 0 || h === 0) return;
            camera.aspect = w / h;
            camera.updateProjectionMatrix();
            renderer.setSize(w, h);
        };
        const resizeObserver = new ResizeObserver(handleResize);
        resizeObserver.observe(mount);

        structureSceneRef.current = {
            scene, camera, renderer, controls, sitesGroup, bondsGroup, cellGroup, highlightGroup, framed: false
        };

        return () => {
            cancelAnimationFrame(animationId);
            resizeObserver.disconnect();
            renderer.domElement.removeEventListener('pointerdown', onPointerDown);
            renderer.domElement.removeEventListener('pointerup', onPointerUp);
            renderer.domElement.removeEventListener('pointermove', onPointerMove);
            controls.dispose();
            renderer.dispose();
            if (renderer.domElement.parentNode === mount) mount.removeChild(renderer.domElement);
            structureSceneRef.current = null;
        };
    }, []);

    // Rebuild the site markers + bonds + unit-cell box, and re-highlight selection.
    useEffect(() => {
        const handle = structureSceneRef.current;
        if (!handle || !sites) return;
        const { sitesGroup, bondsGroup, cellGroup, highlightGroup, camera, controls } = handle;

        const dispose = (group) => {
            while (group.children.length) {
                const child = group.children.pop();
                child.geometry?.dispose();
                child.material?.dispose();
                group.remove(child);
            }
        };
        dispose(sitesGroup);
        dispose(bondsGroup);
        dispose(cellGroup);
        dispose(highlightGroup);

        // Unit-cell vectors = supercell vectors / supercell counts.
        const lattice = sites.latticeVectors;
        const supercell = sites.supercell;
        const unit = lattice.map((row, i) => row.map((value) => value / supercell[i]));
        const toCartesian = (frac) => [0, 1, 2].map(
            (axis) => frac[0] * unit[0][axis] + frac[1] * unit[1][axis] + frac[2] * unit[2][axis]
        );
        const edgeLength = Math.min(...unit.map((row) => Math.hypot(...row)));
        const center = toCartesian([0.5, 0.5, 0.5]);

        // Unit-cell wireframe: the 12 edges of the parallelepiped.
        const cellCorners = [];
        for (const a of [0, 1]) {
            for (const b of [0, 1]) {
                for (const c of [0, 1]) cellCorners.push(toCartesian([a, b, c]));
            }
        }
        const cellPoints = [];
        for (let i = 0; i < 8; i += 1) {
            for (let j = i + 1; j < 8; j += 1) {
                const gi = [(i >> 2) & 1, (i >> 1) & 1, i & 1];
                const gj = [(j >> 2) & 1, (j >> 1) & 1, j & 1];
                if (gi.reduce((s, v, k) => s + (v !== gj[k] ? 1 : 0), 0) === 1) {
                    cellPoints.push(new THREE.Vector3(...cellCorners[i]), new THREE.Vector3(...cellCorners[j]));
                }
            }
        }
        cellGroup.add(new THREE.LineSegments(
            new THREE.BufferGeometry().setFromPoints(cellPoints),
            new THREE.LineBasicMaterial({ color: 0x8a97a8, transparent: true, opacity: 0.4 })
        ));

        const baseRadius = 0.05 * edgeLength;
        const positions = sites.sites.map((site) => ({
            ref: site.referenceNumber,
            el: site.element,
            rms: site.rms,
            axes: site.axes,
            pos: new THREE.Vector3(...toCartesian(site.siteFractional))
        }));

        // Draw each site as its calculated thermal ellipsoid (shape from the PCA
        // eigenvalues/axes) rather than a uniform sphere. The true RMS amplitudes
        // (~0.1 Å) are tiny next to the cell, so magnify them by a single global
        // factor tuned to the old sphere-marker size: the familiar footprint stays,
        // but each site now shows its anisotropy and its size relative to the others.
        const rmsValues = positions.flatMap((p) => p.rms || []);
        const meanRms = rmsValues.length
            ? rmsValues.reduce((sum, value) => sum + value, 0) / rmsValues.length
            : baseRadius;
        const ellipsoidScale = baseRadius / (meanRms || baseRadius);

        // --- Bonds: thin lines between in-cell atoms only (no periodic images, so
        // nothing is drawn outside the box). The cutoff tracks the nearest-
        // neighbour distance so it adapts to the material. -----------------------
        let nearest = Infinity;
        for (let i = 0; i < positions.length; i += 1) {
            for (let j = i + 1; j < positions.length; j += 1) {
                const d = positions[i].pos.distanceTo(positions[j].pos);
                if (d > 0.4 && d < nearest) nearest = d;
            }
        }
        const cutoff = Math.min(3.4, Math.max(2.0, nearest * 1.3));
        const bondPoints = [];
        for (let i = 0; i < positions.length; i += 1) {
            for (let j = i + 1; j < positions.length; j += 1) {
                const d = positions[i].pos.distanceTo(positions[j].pos);
                if (d > 0.4 && d <= cutoff) bondPoints.push(positions[i].pos, positions[j].pos);
            }
        }
        if (bondPoints.length) {
            bondsGroup.add(new THREE.LineSegments(
                new THREE.BufferGeometry().setFromPoints(bondPoints),
                new THREE.LineBasicMaterial({ color: 0x9aa3ad, transparent: true, opacity: 0.7 })
            ));
        }

        // --- Atoms: spheres colored by element. Unselected atoms are translucent
        // so they recede; the selected atom is opaque, bright, and wrapped in a
        // soft highlight cloud plus a correctly-oriented PC1/PC2/PC3 triad. -------
        const sphere = new THREE.SphereGeometry(1, 20, 16);
        const highlightColor = new THREE.Color(SELECTED_ATOM_COLOR);
        positions.forEach((site) => {
            const isSelected = site.ref === selectedRef;
            // The selected atom takes the contrast highlight color so it stands out
            // among the element-colored atoms; the rest keep their element color.
            const color = isSelected ? highlightColor : new THREE.Color(elementColors[site.el] || DEFAULT_ELEMENT_COLOR);
            const marker = new THREE.Mesh(sphere, new THREE.MeshPhongMaterial({
                color,
                emissive: isSelected ? highlightColor.clone().multiplyScalar(0.35) : new THREE.Color(0x000000),
                shininess: isSelected ? 70 : 20,
                transparent: !isSelected,
                opacity: isSelected ? 1 : 0.45
            }));
            marker.userData.referenceNumber = site.ref;
            // Semi-axes of this site's ellipsoid marker (magnified RMS, floored so a
            // near-degenerate third axis never collapses to an invisible sliver);
            // falls back to a uniform sphere when a site lacks PCA data.
            const hasEllipsoid = Boolean(site.rms && site.axes);
            const semi = hasEllipsoid
                ? site.rms.map((r) => Math.max(r * ellipsoidScale, baseRadius * 0.15))
                : [baseRadius, baseRadius, baseRadius];
            if (hasEllipsoid) {
                // Unit sphere scaled by the semi-axes and oriented by the site's
                // principal axes, centered on the site.
                const orientation = new THREE.Matrix4().set(
                    site.axes[0][0], site.axes[1][0], site.axes[2][0], 0,
                    site.axes[0][1], site.axes[1][1], site.axes[2][1], 0,
                    site.axes[0][2], site.axes[1][2], site.axes[2][2], 0,
                    0, 0, 0, 1
                );
                const model = orientation.multiply(new THREE.Matrix4().makeScale(semi[0], semi[1], semi[2]));
                model.setPosition(site.pos.x, site.pos.y, site.pos.z);
                marker.applyMatrix4(model);
            } else {
                marker.position.copy(site.pos);
                marker.scale.setScalar(baseRadius);
            }
            sitesGroup.add(marker);

            if (isSelected) {
                // Size the highlight to the selected marker's own extent so the glow
                // cloud always encloses it and the triad rods reach past it, however
                // large a soft / anisotropic site's ellipsoid gets.
                const extent = Math.max(...semi, baseRadius);
                // Soft concentric highlight shells in the same contrast color.
                [[1.9, 0.18], [1.35, 0.28]].forEach(([scale, opacity]) => {
                    const glow = new THREE.Mesh(sphere, new THREE.MeshBasicMaterial({
                        color: highlightColor, transparent: true, opacity, depthWrite: false
                    }));
                    glow.position.copy(site.pos);
                    glow.scale.setScalar(extent * scale);
                    highlightGroup.add(glow);
                });
                // Local principal-axis triad, using this site's PCA axes (same
                // Cartesian frame as the structure), so the orientation is exact.
                if (selectedEllipsoid?.axes) {
                    buildAxisTriad(site.pos, selectedEllipsoid.axes, extent * 2.3, extent * 0.1)
                        .forEach((rod) => highlightGroup.add(rod));
                }
            }
        });

        // Frame the cell once; keep the user's orientation afterwards.
        controls.target.set(...center);
        if (!handle.framed) {
            const span = Math.max(...unit.map((row) => Math.hypot(...row)));
            const radius = span * 1.7 || 1;
            camera.position.set(center[0] + radius, center[1] + radius * 0.7, center[2] + radius);
            camera.near = radius / 100;
            camera.far = radius * 100;
            camera.updateProjectionMatrix();
            handle.framed = true;
        }
        controls.update();
    }, [sites, selectedRef, selectedEllipsoid, elementColors]);

    // Re-fit the site-ellipsoids canvas to its mount whenever the Displacement-
    // statistics panel — which shares this column — changes height. Its content
    // lands in stages (the PCA result, then the KDE mass levels), each one
    // collapsing the shared row a little more; effects run after layout, so the
    // mount is settled here. In production a ResizeObserver also catches window
    // resizes, but this keeps the canvas exact without depending on its timing.
    useEffect(() => {
        const handle = structureSceneRef.current;
        const mount = structureMountRef.current;
        if (!handle || !mount || !mount.clientWidth || !mount.clientHeight) return;
        const { renderer, camera } = handle;
        const current = renderer.getSize(new THREE.Vector2());
        if (Math.abs(current.x - mount.clientWidth) < 1 && Math.abs(current.y - mount.clientHeight) < 1) return;
        camera.aspect = mount.clientWidth / mount.clientHeight;
        camera.updateProjectionMatrix();
        renderer.setSize(mount.clientWidth, mount.clientHeight);
    }, [sites, selectedRef, selectedEllipsoid, kde]);

    const isoMassLevel = kde?.massLevels?.[isoPercent]?.level;
    // For a reconstructed site, compare its member count against the one-per-cell
    // expectation: an exact match is a clean crystallographic site, a larger count
    // means atoms that don't separate at this distance (close/merged sites or an
    // orientationally-disordered shell), and a smaller one means it fragmented.
    const siteTag = selectedEllipsoid?.copiesPerCell
        ? {
            count: selectedEllipsoid.count,
            per: selectedEllipsoid.copiesPerCell,
            clean: selectedEllipsoid.count === selectedEllipsoid.copiesPerCell,
            status: selectedEllipsoid.count === selectedEllipsoid.copiesPerCell
                ? 'clean site'
                : selectedEllipsoid.count > selectedEllipsoid.copiesPerCell
                    ? 'merged / disordered'
                    : 'fragmented — raise the distance'
        }
        : null;
    const noRun = staticMode && !localFile;

    return (
        <div className="pca-page">
            <div className="pca-controls">
                {/* Site & KDE sampling */}
                <div className="control-group" role="group" aria-label="Site and sampling">
                    <span className="control-group-label">Sample</span>
                    <label className="control">
                        <span className="control-name">
                            Site
                            <InfoBadge label="About the site picker">
                                <p>
                                    Each reference site (an RMCProfile reference number) is one
                                    crystallographic position. Its cloud of per-atom displacements about
                                    the average structure is analysed by PCA to give the thermal ellipsoid.
                                </p>
                            </InfoBadge>
                        </span>
                        <select
                            value={selectedRef ?? ''}
                            onChange={(event) => setSelectedRef(Number(event.target.value))}
                            disabled={!sites}
                            aria-label="Site"
                        >
                            {sites?.sites.map((site) => (
                                <option key={site.referenceNumber} value={site.referenceNumber}>
                                    {`#${site.referenceNumber} ${site.element} — U=${numberFormat(site.uIso, 4)} Å²`}
                                    {site.copiesPerCell ? ` (${site.count}/${site.copiesPerCell})` : ''}
                                </option>
                            ))}
                        </select>
                    </label>
                    <label className="control">
                        <span className="control-name">Grid</span>
                        <select value={grid} onChange={(event) => setGrid(Number(event.target.value))} aria-label="Grid size">
                            {GRID_OPTIONS.map((value) => <option key={value} value={value}>{`${value}³`}</option>)}
                        </select>
                    </label>
                    <label className="control">
                        <span className="control-name">Bandwidth</span>
                        <select value={bw} onChange={(event) => setBw(event.target.value)} aria-label="Bandwidth rule">
                            {BW_OPTIONS.map((option) => <option key={option.value} value={option.value}>{option.label}</option>)}
                        </select>
                    </label>
                    <label className="control">
                        <span className="control-name">Box</span>
                        <input
                            type="range" min="2" max="5" step="0.5"
                            value={extent} onChange={(event) => setExtent(Number(event.target.value))}
                            aria-label="Box half-width in sigma"
                        />
                        <span className="control-value">{extent.toFixed(1)}σ</span>
                    </label>
                    {sites?.reconstructed && (
                        <label className="control">
                            <span className="control-name">
                                Cluster
                                <InfoBadge label="About site clustering">
                                    <p>
                                        This file carries no reference-site or cell columns, so sites are
                                        rebuilt by folding every atom into one unit cell and grouping atoms
                                        of the same element within this distance. Each site should gather one
                                        copy per supercell image; the count beside a site (e.g. 27/27) is its
                                        members against that expected number. Raise the distance to merge
                                        over-split sites, lower it to separate ones that ran together.
                                    </p>
                                </InfoBadge>
                            </span>
                            <input
                                type="range" min="0.4" max="2.5" step="0.1"
                                value={clusterThreshold} onChange={(event) => setClusterThreshold(Number(event.target.value))}
                                aria-label="Site clustering distance in Angstrom"
                            />
                            <span className="control-value">{clusterThreshold.toFixed(1)} Å</span>
                        </label>
                    )}
                </div>

                {/* Ellipsoid — the wireframe reference and the KDE density painted on it
                    (two aspects of the same surface, so they share one control group). */}
                <div className="control-group" role="group" aria-label="Ellipsoid and density shell">
                    <span className="control-group-label">Ellipsoid</span>
                    <label className="control switch">
                        <span className="control-name">Wireframe</span>
                        <input type="checkbox" checked={showEllipsoid} onChange={(event) => setShowEllipsoid(event.target.checked)} aria-label="Show ellipsoid wireframe" />
                        <i className="switch-track" aria-hidden="true" />
                    </label>
                    <label className="control">
                        <span className="control-name">
                            Level
                            <InfoBadge label="About the ellipsoid level">
                                <p>
                                    The enclosed-probability level drawn as the thermal-ellipsoid
                                    wireframe. 50% is the crystallographic convention.
                                </p>
                            </InfoBadge>
                        </span>
                        <input
                            type="range" min="0.1" max="0.99" step="0.01"
                            value={probability} onChange={(event) => setProbability(Number(event.target.value))}
                            aria-label="Ellipsoid probability"
                        />
                        <span className="control-value">{Math.round(probability * 100)}%</span>
                    </label>
                    <label className="control">
                        <span className="control-name">Color</span>
                        <i className="ellipsoid-color-swatch" style={{ background: ellipsoidColor }} aria-hidden="true" />
                        <select
                            value={ellipsoidColor}
                            onChange={(event) => setEllipsoidColor(event.target.value)}
                            aria-label="Ellipsoid wireframe color"
                        >
                            {ELLIPSOID_COLOR_OPTIONS.map((option) => (
                                <option key={option.value} value={option.value}>{option.label}</option>
                            ))}
                        </select>
                    </label>
                    {/* The InfoBadge is a sibling of the toggle label, not inside it, so the
                        "?" never sits inside the switch's clickable area (an interactive
                        element nested in a <label> makes the toggle click ambiguous). */}
                    <div className="switch-with-info">
                        <label className="control switch">
                            <span className="control-name">Shell</span>
                            <input
                                type="checkbox"
                                checked={showEllipsoidKde}
                                aria-label="Show KDE density shell"
                                onChange={(event) => {
                                    const on = event.target.checked;
                                    setShowEllipsoidKde(on);
                                    if (on) setShowSurface(false);
                                }}
                            />
                            <i className="switch-track" aria-hidden="true" />
                        </label>
                        <InfoBadge label="About the KDE shell" align="end">
                            <p>
                                Paints the KDE density onto the ellipsoid surface. A near-uniform
                                color means the motion is Gaussian; hotter and colder patches mark
                                where the real density departs from the harmonic ellipsoid.
                            </p>
                            <p>
                                It shows the same density as the isosurface from the outside, so the
                                two switch off each other.
                            </p>
                        </InfoBadge>
                    </div>
                    <label className="control">
                        <span className="control-name">Colormap</span>
                        <select value={shellColormap} onChange={(event) => setShellColormap(event.target.value)} aria-label="KDE shell colormap">
                            {COLORMAP_NAMES.map((name) => <option key={name} value={name}>{name}</option>)}
                        </select>
                    </label>
                    <label className="control">
                        <span className="control-name">
                            Contrast
                            <InfoBadge label="About the shell contrast" align="end">
                                <p>
                                    Stretches the shell colormap around its mid-tone: higher values
                                    push the effective range (vmin/vmax) together so faint departures
                                    from the ellipsoid stand out, lower values flatten it.
                                </p>
                            </InfoBadge>
                        </span>
                        <input
                            type="range" min="0.5" max="3" step="0.1"
                            value={shellContrast} onChange={(event) => setShellContrast(Number(event.target.value))}
                            aria-label="KDE shell contrast"
                        />
                        <span className="control-value">{shellContrast.toFixed(1)}×</span>
                    </label>
                </div>

                {/* Isosurface & wall projections — the volume density views (shared colormap) */}
                <div className="control-group" role="group" aria-label="Isosurface and wall projections">
                    <span className="control-group-label">Volume</span>
                    <label className="control switch">
                        <span className="control-name">
                            Isosurface
                            <InfoBadge label="About the isosurface">
                                <p>
                                    The KDE density isosurface enclosing this fraction of the cloud's
                                    mass. Compare its shape to the harmonic ellipsoid — where it sits
                                    inside, the motion is anharmonic (see Non-Gaussianity).
                                </p>
                            </InfoBadge>
                        </span>
                        <input
                            type="checkbox"
                            checked={showSurface}
                            aria-label="Show isosurface"
                            onChange={(event) => {
                                const on = event.target.checked;
                                setShowSurface(on);
                                // The isosurface is a clean standalone view of the density volume:
                                // turning it on clears the ellipsoid wireframe and its KDE shell
                                // (the shell shows the same density, drawn on the surface instead).
                                if (on) { setShowEllipsoidKde(false); setShowEllipsoid(false); }
                            }}
                        />
                        <i className="switch-track" aria-hidden="true" />
                    </label>
                    <label className="control">
                        <span className="control-name">Level</span>
                        <input
                            type="range" min="1" max="99" step="1"
                            value={isoPercent} onChange={(event) => setIsoPercent(Number(event.target.value))}
                            aria-label="Isosurface enclosed mass"
                        />
                        <span className="control-value">{isoPercent}%</span>
                    </label>
                    <label className="control switch">
                        <span className="control-name">Projections</span>
                        <input type="checkbox" checked={showProjections} onChange={(event) => setShowProjections(event.target.checked)} aria-label="Show wall projections" />
                        <i className="switch-track" aria-hidden="true" />
                    </label>
                    <label className="control">
                        <span className="control-name">Colormap</span>
                        <select value={colormap} onChange={(event) => setColormap(event.target.value)} aria-label="Isosurface and projections colormap">
                            {COLORMAP_NAMES.map((name) => <option key={name} value={name}>{name}</option>)}
                        </select>
                    </label>
                </div>
            </div>

            {noRun && (
                <p className="pca-hint">Open a run folder (with an <code>.rmc6f</code> file) to view thermal ellipsoids.</p>
            )}
            {sitesError && <p className="pca-error-banner">{sitesError}</p>}

            <div className="pca-layout">
                <div className="pca-panel pca-viewport">
                    <h3>
                        <span className="panel-title-label">
                            {selectedEllipsoid
                                ? `${selectedEllipsoid.element} site #${selectedEllipsoid.referenceNumber}`
                                : 'PCA ellipsoid'}
                        </span>
                        <span className="panel-title-actions">
                            {selectedEllipsoid && (
                                <span className="panel-title-count">{selectedEllipsoid.count.toLocaleString()} atoms</span>
                            )}
                            <button
                                type="button"
                                className="pca-reset-view"
                                onClick={resetMainView}
                                disabled={!kde}
                                title="Reset the camera to the default view"
                            >
                                <svg width="12" height="12" viewBox="0 0 24 24" aria-hidden="true" fill="none" stroke="currentColor" strokeWidth="2.2" strokeLinecap="round" strokeLinejoin="round">
                                    <path d="M3 12a9 9 0 1 0 9-9 9.75 9.75 0 0 0-6.74 2.74L3 8" />
                                    <path d="M3 3v5h5" />
                                </svg>
                                Reset view
                            </button>
                            <SaveMenu
                                onSave={saveMainView}
                                options={SAVE_OPTIONS}
                                label="Save"
                                align="right"
                                disabled={!kde}
                            />
                        </span>
                    </h3>
                    <div className="pca-canvas" ref={mountRef}>
                        {(loadingKde || loadingSites) && <div className="pca-badge">Computing…</div>}
                        {kdeError && <div className="pca-badge is-error">{kdeError}</div>}
                        {kde && (
                            <>
                                <div className="pca-view-controls pca-view-controls--left">
                                    <div className="pca-view-group" role="group" aria-label="Projection">
                                        <button
                                            type="button"
                                            className={`pca-view-btn ${perspective ? 'is-active' : ''}`}
                                            onClick={() => applyPerspective(true)}
                                            aria-pressed={perspective}
                                            title="Perspective projection"
                                        >
                                            Perspective
                                        </button>
                                        <button
                                            type="button"
                                            className={`pca-view-btn ${perspective ? '' : 'is-active'}`}
                                            onClick={() => applyPerspective(false)}
                                            aria-pressed={!perspective}
                                            title="Orthographic (parallel) projection"
                                        >
                                            Orthographic
                                        </button>
                                    </div>
                                </div>
                                <div className="pca-view-controls pca-view-controls--right">
                                    <div className="pca-view-group" role="group" aria-label="Camera orientation">
                                        <span className="pca-view-group-label">Along</span>
                                        {[0, 1, 2].map((axisIndex) => (
                                            <button
                                                key={axisIndex}
                                                type="button"
                                                className={`pca-view-btn ${viewAxis === axisIndex ? 'is-active' : ''}`}
                                                onClick={() => lookAlong(axisIndex)}
                                                aria-pressed={viewAxis === axisIndex}
                                                title={`Look down PC${axisIndex + 1}`}
                                            >
                                                {`PC${axisIndex + 1}`}
                                            </button>
                                        ))}
                                    </div>
                                </div>
                            </>
                        )}
                    </div>
                    <div className="pca-legend">
                        <span className="pca-legend-item"><i className="pca-legend-swatch" style={{ background: '#d64545' }} /> PC1</span>
                        <span className="pca-legend-item"><i className="pca-legend-swatch" style={{ background: '#3fa34d' }} /> PC2</span>
                        <span className="pca-legend-item"><i className="pca-legend-swatch" style={{ background: '#3f7fd6' }} /> PC3</span>
                        <span className="pca-legend-item"><i className="pca-legend-swatch" style={{ background: ellipsoidColor }} /> {Math.round(probability * 100)}% ellipsoid</span>
                        <span className="pca-legend-item pca-legend-note">walls: PC-plane density projections</span>
                        <a
                            className="pca-legend-credit"
                            href="https://github.com/MaximEremenko/Utilities/tree/main/RMCProfileUtilities/PCA_KDE"
                            target="_blank"
                            rel="noreferrer"
                        >
                            Analysis after Maksim Eremenko&rsquo;s PCA_KDE
                        </a>
                    </div>
                </div>

                <div className="pca-side">
                    <div className="pca-panel pca-stats-panel">
                        <h3>
                            <span className="panel-title-label">
                                Displacement statistics
                                <InfoBadge label="About the displacement statistics" align="end">
                                    <p>
                                        Anisotropic displacement parameters from the site's cloud
                                        covariance: U<sub>iso</sub>/B<sub>iso</sub> are the isotropic
                                        equivalents, and the tables below give the covariance tensor
                                        with its principal axes (PCA components) and per-axis amplitudes.
                                    </p>
                                </InfoBadge>
                            </span>
                        </h3>
                        <div className="pca-stats-body">
                        {selectedEllipsoid ? (
                            <>
                                {siteTag && (
                                    <p className={`pca-site-tag ${siteTag.clean ? 'is-clean' : 'is-flagged'}`}>
                                        <span className="pca-site-tag-count">{siteTag.count}/{siteTag.per}</span>
                                        <span>copies · {siteTag.status}</span>
                                        <InfoBadge label="About the reconstructed site" align="end">
                                            <p>
                                                Sites were rebuilt by folding + clustering because this file
                                                has no reference-site columns. The count is this site's members
                                                against the one-per-cell expectation from the supercell. A larger
                                                count means atoms that don't resolve into separate sites at the
                                                current distance — close sites, or an orientationally disordered
                                                group such as a rotor shell, whose “ellipsoid” is really a shell
                                                and is best read from the KDE rather than as a thermal ellipsoid.
                                            </p>
                                        </InfoBadge>
                                    </p>
                                )}
                                <dl className="pca-stats">
                                    <div className="pca-stat">
                                        <dt>U<sub>iso</sub> (Å²)</dt>
                                        <dd>{numberFormat(selectedEllipsoid.uIso)}</dd>
                                    </div>
                                    <div className="pca-stat">
                                        <dt>B<sub>iso</sub> (Å²)</dt>
                                        <dd>{numberFormat(selectedEllipsoid.bIso, 3)}</dd>
                                    </div>
                                    <div className="pca-stat">
                                        <dt>Anisotropy</dt>
                                        <dd>{numberFormat(selectedEllipsoid.anisotropy, 2)}{selectedEllipsoid.degenerate ? ' · degenerate' : ''}</dd>
                                    </div>
                                    <div className="pca-stat">
                                        <dt>
                                            Non-Gaussianity
                                            <InfoBadge label="About non-Gaussianity" align="end">
                                                <p>
                                                    Mean excess kurtosis of the displacement cloud
                                                    (0 = harmonic/Gaussian). Positive means a peaked,
                                                    fat-tailed distribution — the KDE isosurface then
                                                    sits inside the ellipsoid, signalling anharmonic
                                                    motion or split sites.
                                                </p>
                                            </InfoBadge>
                                        </dt>
                                        <dd>{numberFormat(selectedEllipsoid.nonGaussianity ?? kde?.nonGaussianity, 2)}</dd>
                                    </div>
                                </dl>
                                <div className="pca-matrices">
                                    <div className="pca-matrix-block">
                                        <div className="pca-matrix-title">
                                            Covariance U (Å²)
                                            <InfoBadge label="About the covariance matrix" align="end">
                                                <p>
                                                    The displacement covariance tensor in Cartesian
                                                    (x, y, z) axes — the anisotropic displacement
                                                    parameters. Its eigen-decomposition gives the
                                                    principal axes below.
                                                </p>
                                            </InfoBadge>
                                        </div>
                                        <div className="pca-matrix-scroll">
                                            <table className="pca-matrix">
                                                <thead>
                                                    <tr>
                                                        <th aria-hidden="true" />
                                                        <th scope="col">x</th>
                                                        <th scope="col">y</th>
                                                        <th scope="col">z</th>
                                                    </tr>
                                                </thead>
                                                <tbody>
                                                    {['x', 'y', 'z'].map((label, i) => (
                                                        <tr key={label}>
                                                            <th scope="row">{label}</th>
                                                            {selectedEllipsoid.covariance[i].map((value, j) => (
                                                                <td key={j} className={i === j ? 'is-diagonal' : ''}>
                                                                    {numberFormat(value, 4)}
                                                                </td>
                                                            ))}
                                                        </tr>
                                                    ))}
                                                </tbody>
                                            </table>
                                        </div>
                                    </div>

                                    <div className="pca-matrix-block">
                                        <div className="pca-matrix-title">
                                            Principal axes
                                            <InfoBadge label="About the principal axes" align="end">
                                                <p>
                                                    The PCA components (eigenvectors of the
                                                    covariance), one per row: each is a unit direction
                                                    in Cartesian (x, y, z). λ is its eigenvalue
                                                    (variance, Å²), RMS the amplitude (Å), and κ the
                                                    excess kurtosis along it (0 = Gaussian).
                                                </p>
                                            </InfoBadge>
                                        </div>
                                        <div className="pca-matrix-scroll">
                                            <table className="pca-matrix pca-matrix--axes">
                                                <thead>
                                                    <tr>
                                                        <th aria-hidden="true" />
                                                        <th scope="col">x</th>
                                                        <th scope="col">y</th>
                                                        <th scope="col">z</th>
                                                        <th scope="col">λ (Å²)</th>
                                                        <th scope="col">RMS (Å)</th>
                                                        <th scope="col">
                                                            <abbr title="Excess kurtosis along this axis (0 = Gaussian)">κ</abbr>
                                                        </th>
                                                    </tr>
                                                </thead>
                                                <tbody>
                                                    {selectedEllipsoid.axes.map((axis, i) => (
                                                        <tr key={i}>
                                                            <th scope="row">
                                                                <span
                                                                    className="pca-pc-dot"
                                                                    style={{ background: PC_CSS_COLORS[i] }}
                                                                    aria-hidden="true"
                                                                />
                                                                PC{i + 1}
                                                            </th>
                                                            {axis.map((value, j) => (
                                                                <td key={j}>{numberFormat(value, 3)}</td>
                                                            ))}
                                                            <td>{numberFormat(selectedEllipsoid.eigenvalues[i], 4)}</td>
                                                            <td>{numberFormat(selectedEllipsoid.rms[i], 3)}</td>
                                                            <td>{numberFormat(selectedEllipsoid.excessKurtosis?.[i], 2)}</td>
                                                        </tr>
                                                    ))}
                                                </tbody>
                                            </table>
                                        </div>
                                    </div>
                                </div>
                                {kde && (
                                    <p className="pca-meta">
                                        Volume {kde.grid}³ · fit {kde.fitCount.toLocaleString()}/{kde.count.toLocaleString()} ·
                                        captured mass {numberFormat(kde.mass * 100, 1)}%
                                        {Number.isFinite(isoMassLevel) ? ` · iso @ ${isoPercent}% mass` : ''}
                                        {kde.browserPcaKde ? ' · browser' : ' · server'}
                                    </p>
                                )}
                            </>
                        ) : (
                            <p className="pca-meta">{loadingSites ? 'Loading sites…' : 'No site selected.'}</p>
                        )}
                        </div>
                    </div>

                    <div className="pca-panel pca-unitcell-panel">
                        <h3>
                            <span className="panel-title-label">
                                Site ellipsoids
                                <InfoBadge label="About the site ellipsoids" align="end">
                                    <p>
                                        Every reference site in the average unit cell, drawn as its
                                        calculated thermal ellipsoid and colored by element. Click one to
                                        load its PCA-KDE in the main panel; the selected site is highlighted.
                                    </p>
                                </InfoBadge>
                            </span>
                        </h3>
                        <div className="pca-structure" ref={structureMountRef} />
                        {sites?.elements?.length > 0 && (
                            <div className="pca-legend">
                                {sites.elements.map((element) => (
                                    <span key={element} className="pca-legend-item">
                                        <i
                                            className="pca-legend-swatch"
                                            style={{ background: elementColors[element] || DEFAULT_ELEMENT_COLOR }}
                                        />
                                        {element}
                                    </span>
                                ))}
                            </div>
                        )}
                    </div>
                </div>
            </div>

            <footer className="app-footer">
                &copy; 2026 Tsung-Han Yang &middot;{' '}
                <a
                    href="https://github.com/drthyang/rmc-toolkits/blob/main/LICENSE"
                    target="_blank"
                    rel="noreferrer"
                >
                    AGPLv3
                </a>
                {' '}&middot;{' '}
                <a
                    href="https://github.com/drthyang/rmc-toolkits#readme"
                    target="_blank"
                    rel="noreferrer"
                >
                    About & documentation
                </a>
            </footer>
        </div>
    );
}
