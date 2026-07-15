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
import InfoBadge from './InfoBadge';
import './PcaKdePage.css';

const GRID_OPTIONS = [24, 32, 40, 48, 56, 64];
const BW_OPTIONS = [
    { value: 'scott', label: 'Scott' },
    { value: 'silverman', label: 'Silverman' }
];
const DEFAULTS = { grid: 40, bw: 'scott', extent: 4, probability: 0.5, isoPercent: 25, colormap: 'viridis' };

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

export default function PcaKdePage({ directory, localRun }) {
    const [rmc6fText, setRmc6fText] = useState(null);
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
    const [probability, setProbability] = useState(DEFAULTS.probability);
    const [isoPercent, setIsoPercent] = useState(DEFAULTS.isoPercent);
    const [colormap, setColormap] = useState(DEFAULTS.colormap);
    const [showEllipsoid, setShowEllipsoid] = useState(true);
    const [showSurface, setShowSurface] = useState(true);
    const [showProjections, setShowProjections] = useState(true);

    const workerRef = useRef(null);
    const requestIdRef = useRef(0);
    const mountRef = useRef(null);
    const sceneRef = useRef(null);
    const framedRefRef = useRef(null);
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

    // --- Read the raw .rmc6f text once per local file. ------------------------
    useEffect(() => {
        let cancelled = false;
        if (localFile) {
            localFile.text().then((text) => {
                if (!cancelled) setRmc6fText(text);
            }).catch(() => {
                if (!cancelled) { setRmc6fText(null); setSitesError('Could not read the structure file.'); }
            });
        } else {
            setRmc6fText(null);
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
                worker.postMessage({ id, kind, text: rmc6fText, cacheKey: structureFile?.path || 'run', ...params });
            });
        }
        const endpoint = kind === 'sites' ? '/api/pca/sites' : '/api/pca/kde';
        return axios
            .get(`${API_BASE_URL}${endpoint}`, { params: { dir: directory || '.', ...params } })
            .then((response) => response.data);
    }, [rmc6fText, structureFile, directory]);

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
                const data = await requestPca('sites', { probability });
                if (cancelled) return;
                setSites(data);
                setSelectedRef((current) => {
                    const refs = data.sites.map((site) => site.referenceNumber);
                    return current != null && refs.includes(current) ? current : refs[0] ?? null;
                });
            } catch (error) {
                if (!cancelled) { setSites(null); setSitesError(error.message); }
            } finally {
                if (!cancelled) setLoadingSites(false);
            }
        };
        loadSites();
        return () => { cancelled = true; };
    }, [requestPca, localFile, rmc6fText, probability]);

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
                    { referenceNumber: selectedRef, grid, bw, extent, probability, cubicBox: true, projections: true }
                );
                if (!cancelled) setKde(data);
            } catch (error) {
                if (!cancelled) { setKde(null); setKdeError(error.message); }
            } finally {
                if (!cancelled) setLoadingKde(false);
            }
        };
        loadKde();
        return () => { cancelled = true; };
    }, [requestPca, localFile, rmc6fText, selectedRef, grid, bw, extent, probability]);

    const selectedEllipsoid = useMemo(
        () => sites?.sites.find((site) => site.referenceNumber === selectedRef) || null,
        [sites, selectedRef]
    );

    const elementColors = useMemo(
        () => buildElementColors(sites?.elements ?? []),
        [sites]
    );

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
        const camera = new THREE.PerspectiveCamera(45, width / height, 0.001, 1000);
        camera.position.set(0.6, 0.5, 0.9);
        const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
        renderer.setPixelRatio(window.devicePixelRatio || 1);
        renderer.setSize(width, height);
        mount.appendChild(renderer.domElement);

        const controls = new OrbitControls(camera, renderer.domElement);
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
            camera.aspect = w / h;
            camera.updateProjectionMatrix();
            renderer.setSize(w, h);
        };
        window.addEventListener('resize', handleResize);
        // Track the mount itself so the canvas follows the panel's aspect-ratio
        // box and the responsive column layout, not just window resizes.
        const resizeObserver = new ResizeObserver(handleResize);
        resizeObserver.observe(mount);

        sceneRef.current = { scene, camera, renderer, controls, surfaceGroup, ellipsoidGroup, axesGroup, wallsGroup };

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
            const material = new THREE.MeshBasicMaterial({ color: 0xff5a5a, wireframe: true, transparent: true, opacity: 0.5 });
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

        // Reframe the camera only when the displayed volume is for a new site, so
        // toggling layers or sweeping a slider doesn't yank the view. Keying off
        // the volume's own reference number (not selectedRef) avoids a premature
        // reframe against the previous site's axes while its KDE is still loading —
        // that mismatch was what made the view snap to a wrong angle on click.
        // The default looks straight down the box body diagonal (the +PC1/+PC2/+PC3
        // direction), so the three min-corner walls sit symmetrically at the back;
        // the distance fits the whole cube corner-on with margin.
        const kdeRef = kde.referenceNumber ?? selectedRef;
        const radius = axisLength * 4.3 || 1;
        controls.target.copy(meanVec);
        if (framedRefRef.current !== kdeRef) {
            const dir = new THREE.Vector3(
                axes[0][0] + axes[1][0] + axes[2][0],
                axes[0][1] + axes[1][1] + axes[2][1],
                axes[0][2] + axes[1][2] + axes[2][2]
            ).normalize();
            camera.up.set(axes[2][0], axes[2][1], axes[2][2]);
            camera.position.copy(meanVec).addScaledVector(dir, radius);
            framedRefRef.current = kdeRef;
        }
        camera.near = radius / 100;
        camera.far = radius * 100;
        camera.updateProjectionMatrix();
        controls.update();
    }, [kde, isoPercent, colormap, showEllipsoid, showSurface, showProjections, selectedEllipsoid, selectedRef]);

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
            pos: new THREE.Vector3(...toCartesian(site.siteFractional))
        }));

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
            marker.position.copy(site.pos);
            marker.scale.setScalar(baseRadius);
            marker.userData.referenceNumber = site.ref;
            sitesGroup.add(marker);

            if (isSelected) {
                // Soft concentric highlight shells in the same contrast color.
                [[2.6, 0.18], [1.7, 0.28]].forEach(([scale, opacity]) => {
                    const glow = new THREE.Mesh(sphere, new THREE.MeshBasicMaterial({
                        color: highlightColor, transparent: true, opacity, depthWrite: false
                    }));
                    glow.position.copy(site.pos);
                    glow.scale.setScalar(baseRadius * scale);
                    highlightGroup.add(glow);
                });
                // Local principal-axis triad, using this site's PCA axes (same
                // Cartesian frame as the structure), so the orientation is exact.
                if (selectedEllipsoid?.axes) {
                    buildAxisTriad(site.pos, selectedEllipsoid.axes, baseRadius * 3.2, baseRadius * 0.14)
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

    const isoMassLevel = kde?.massLevels?.[isoPercent]?.level;
    const noRun = staticMode && !localFile;

    return (
        <div className="pca-page">
            <div className="pca-controls">
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
                <label className="control">
                    <span className="control-name">
                        Ellipsoid
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
                        type="range" min="1" max="99" step="1"
                        value={isoPercent} onChange={(event) => setIsoPercent(Number(event.target.value))}
                        aria-label="Isosurface enclosed mass"
                    />
                    <span className="control-value">{isoPercent}%</span>
                </label>
                <label className="control">
                    <span className="control-name">Colormap</span>
                    <select value={colormap} onChange={(event) => setColormap(event.target.value)} aria-label="Colormap">
                        {COLORMAP_NAMES.map((name) => <option key={name} value={name}>{name}</option>)}
                    </select>
                </label>
                <label className="control switch">
                    <span className="control-name">Isosurface</span>
                    <input type="checkbox" checked={showSurface} onChange={(event) => setShowSurface(event.target.checked)} />
                    <i className="switch-track" aria-hidden="true" />
                </label>
                <label className="control switch">
                    <span className="control-name">Ellipsoid</span>
                    <input type="checkbox" checked={showEllipsoid} onChange={(event) => setShowEllipsoid(event.target.checked)} />
                    <i className="switch-track" aria-hidden="true" />
                </label>
                <label className="control switch">
                    <span className="control-name">Projections</span>
                    <input type="checkbox" checked={showProjections} onChange={(event) => setShowProjections(event.target.checked)} />
                    <i className="switch-track" aria-hidden="true" />
                </label>
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
                        {selectedEllipsoid && (
                            <span className="panel-title-count">{selectedEllipsoid.count.toLocaleString()} atoms</span>
                        )}
                    </h3>
                    <div className="pca-canvas" ref={mountRef}>
                        {(loadingKde || loadingSites) && <div className="pca-badge">Computing…</div>}
                        {kdeError && <div className="pca-badge is-error">{kdeError}</div>}
                    </div>
                    <div className="pca-legend">
                        <span className="pca-legend-item"><i className="pca-legend-swatch" style={{ background: '#d64545' }} /> PC1</span>
                        <span className="pca-legend-item"><i className="pca-legend-swatch" style={{ background: '#3fa34d' }} /> PC2</span>
                        <span className="pca-legend-item"><i className="pca-legend-swatch" style={{ background: '#3f7fd6' }} /> PC3</span>
                        <span className="pca-legend-item"><i className="pca-legend-swatch" style={{ background: '#ff5a5a' }} /> {Math.round(probability * 100)}% ellipsoid</span>
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
                    <div className="pca-panel">
                        <h3>
                            <span className="panel-title-label">
                                Displacement statistics
                                <InfoBadge label="About the displacement statistics" align="end">
                                    <p>
                                        Anisotropic displacement parameters from the site's cloud
                                        covariance: U<sub>iso</sub>/B<sub>iso</sub> are the isotropic
                                        equivalents, RMS axes are the principal amplitudes.
                                    </p>
                                </InfoBadge>
                            </span>
                        </h3>
                        {selectedEllipsoid ? (
                            <>
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
                                        <dt>RMS axes (Å)</dt>
                                        <dd>{selectedEllipsoid.rms.map((value) => numberFormat(value, 3)).join(', ')}</dd>
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
                                    <div className="pca-stat">
                                        <dt>Eigenvalues (Å²)</dt>
                                        <dd>{selectedEllipsoid.eigenvalues.map((value) => numberFormat(value, 4)).join(', ')}</dd>
                                    </div>
                                </dl>
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

                    <div className="pca-panel">
                        <h3>
                            <span className="panel-title-label">
                                Unit cell
                                <InfoBadge label="About the unit cell view" align="end">
                                    <p>
                                        Every reference site in the average unit cell, colored by
                                        element. Click an atom to load its PCA-KDE in the main panel;
                                        the selected site is enlarged and highlighted.
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
