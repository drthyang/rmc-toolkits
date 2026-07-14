// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import API_BASE_URL from '../api';
import { isStaticMode } from '../browserData';
import { COLORMAP_NAMES, getLut, sampleColormap } from '../colormaps';
import { marchingCubes } from '../workers/marchingCubes';
import InfoBadge from './InfoBadge';
import './PcaKdePage.css';

const GRID_OPTIONS = [24, 32, 40, 48, 56, 64];
const BW_OPTIONS = [
    { value: 'scott', label: 'Scott' },
    { value: 'silverman', label: 'Silverman' }
];
const DEFAULTS = { grid: 40, bw: 'scott', extent: 3, probability: 0.5, isoPercent: 50, colormap: 'viridis' };

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

// Draw a 2D KDE projection (density[rows][cols]) into a canvas with a colormap.
const paintProjection = (canvas, projection, colormap) => {
    if (!canvas || !projection) return;
    const density = projection.density;
    const rows = density.length;
    const cols = density[0] ? density[0].length : 0;
    if (!rows || !cols) return;
    canvas.width = cols;
    canvas.height = rows;
    const ctx = canvas.getContext('2d');
    const image = ctx.createImageData(cols, rows);
    const lut = getLut(colormap);
    const vmax = projection.vmax || 1;
    for (let r = 0; r < rows; r += 1) {
        // Flip vertically so the second axis increases upward on screen.
        const srcRow = rows - 1 - r;
        for (let c = 0; c < cols; c += 1) {
            const value = Math.max(0, Math.min(1, density[srcRow][c] / vmax));
            const lutIndex = Math.round(value * (lut.length / 3 - 1)) * 3;
            const pixel = (r * cols + c) * 4;
            image.data[pixel] = lut[lutIndex];
            image.data[pixel + 1] = lut[lutIndex + 1];
            image.data[pixel + 2] = lut[lutIndex + 2];
            image.data[pixel + 3] = 255;
        }
    }
    ctx.putImageData(image, 0, 0);
};

const PROJECTION_META = [
    { key: 'pc12', label: 'PC1 – PC2' },
    { key: 'pc13', label: 'PC1 – PC3' },
    { key: 'pc23', label: 'PC2 – PC3' }
];

export default function PcaKdePage({ directory, localRun, theme = 'light' }) {
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

    const workerRef = useRef(null);
    const requestIdRef = useRef(0);
    const mountRef = useRef(null);
    const sceneRef = useRef(null);

    const staticMode = isStaticMode();
    const structureFile = localRun?.structureFile || null;

    // --- Load the raw .rmc6f text once per selected file (static mode). --------
    useEffect(() => {
        let cancelled = false;
        if (staticMode) {
            if (structureFile?.sourceFile) {
                structureFile.sourceFile.text().then((text) => {
                    if (!cancelled) setRmc6fText(text);
                }).catch(() => {
                    if (!cancelled) { setRmc6fText(null); setSitesError('Could not read the structure file.'); }
                });
            } else {
                setRmc6fText(null);
            }
        }
        return () => { cancelled = true; };
    }, [staticMode, structureFile]);

    // --- Static-mode worker lifecycle. ----------------------------------------
    useEffect(() => {
        if (!staticMode) return undefined;
        const worker = new Worker(new URL('../workers/pcaKdeWorker.js', import.meta.url), { type: 'module' });
        workerRef.current = worker;
        return () => { worker.terminate(); workerRef.current = null; };
    }, [staticMode]);

    // A single request path for both runtimes: Flask GET, or the static worker.
    const requestPca = useCallback((kind, params) => {
        if (staticMode) {
            return new Promise((resolve, reject) => {
                const worker = workerRef.current;
                if (!worker) { reject(new Error('PCA-KDE worker unavailable')); return; }
                if (!rmc6fText) { reject(new Error('Open a run folder to view thermal ellipsoids.')); return; }
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
    }, [staticMode, rmc6fText, structureFile, directory]);

    // --- Load the per-site ellipsoid table. -----------------------------------
    useEffect(() => {
        let cancelled = false;
        const loadSites = async () => {
            if (staticMode && !rmc6fText) { setSites(null); return; }
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
    }, [requestPca, staticMode, rmc6fText, probability]);

    // --- Load the KDE volume for the selected site. ---------------------------
    useEffect(() => {
        let cancelled = false;
        const loadKde = async () => {
            if (selectedRef == null) { setKde(null); return; }
            if (staticMode && !rmc6fText) return;
            setLoadingKde(true);
            setKdeError(null);
            try {
                const data = await requestPca(
                    'kde',
                    { referenceNumber: selectedRef, grid, bw, extent, probability, projections: true }
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
    }, [requestPca, staticMode, rmc6fText, selectedRef, grid, bw, extent, probability]);

    const selectedEllipsoid = useMemo(
        () => sites?.sites.find((site) => site.referenceNumber === selectedRef) || null,
        [sites, selectedRef]
    );

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
        scene.add(surfaceGroup, ellipsoidGroup, axesGroup);

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
            camera.aspect = w / h;
            camera.updateProjectionMatrix();
            renderer.setSize(w, h);
        };
        window.addEventListener('resize', handleResize);

        sceneRef.current = { scene, camera, renderer, controls, surfaceGroup, ellipsoidGroup, axesGroup };

        return () => {
            cancelAnimationFrame(animationId);
            window.removeEventListener('resize', handleResize);
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
        const { surfaceGroup, ellipsoidGroup, axesGroup, camera, controls } = handle;

        const dispose = (group) => {
            while (group.children.length) {
                const child = group.children.pop();
                child.geometry?.dispose();
                child.material?.dispose();
                group.remove(child);
            }
        };
        dispose(surfaceGroup);
        dispose(ellipsoidGroup);
        dispose(axesGroup);

        const axes = kde.axes;
        const mean = kde.mean;
        const axisCoords = kde.axisCoords;
        const gridN = kde.grid;

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

        // Principal-axis triad (PC1 red, PC2 green, PC3 blue), scaled to the box.
        const axisLength = Math.max(...kde.halfWidths);
        const triadColors = [0xd64545, 0x3fa34d, 0x3f7fd6];
        for (let a = 0; a < 3; a += 1) {
            const dir = new THREE.Vector3(axes[a][0], axes[a][1], axes[a][2]);
            const points = [
                new THREE.Vector3(mean[0], mean[1], mean[2]),
                new THREE.Vector3(
                    mean[0] + dir.x * axisLength,
                    mean[1] + dir.y * axisLength,
                    mean[2] + dir.z * axisLength
                )
            ];
            const geometry = new THREE.BufferGeometry().setFromPoints(points);
            axesGroup.add(new THREE.Line(geometry, new THREE.LineBasicMaterial({ color: triadColors[a] })));
        }

        // Frame the camera to the box on first load of a site.
        const radius = axisLength * 2.4 || 1;
        controls.target.set(mean[0], mean[1], mean[2]);
        const offset = new THREE.Vector3(0.8, 0.6, 1).normalize().multiplyScalar(radius);
        camera.position.copy(controls.target).add(offset);
        camera.near = radius / 100;
        camera.far = radius * 100;
        camera.updateProjectionMatrix();
        controls.update();
    }, [kde, isoPercent, colormap, showEllipsoid, showSurface, selectedEllipsoid]);

    // --- Paint the projection canvases. ---------------------------------------
    const projectionRefs = { pc12: useRef(null), pc13: useRef(null), pc23: useRef(null) };
    useEffect(() => {
        if (!kde?.projections) return;
        PROJECTION_META.forEach(({ key }) => {
            paintProjection(projectionRefs[key].current, kde.projections[key], colormap);
        });
    }, [kde, colormap]); // eslint-disable-line react-hooks/exhaustive-deps

    const isoMassLevel = kde?.massLevels?.[isoPercent]?.level;
    const noRun = staticMode && !structureFile;

    return (
        <div className="pca-page" data-theme={theme}>
            <div className="pca-controls">
                <div className="pca-control">
                    <label htmlFor="pca-site">Site
                        <InfoBadge text="Reference site (RMCProfile reference number). Its displacement cloud is analyzed by PCA." />
                    </label>
                    <select
                        id="pca-site"
                        value={selectedRef ?? ''}
                        onChange={(event) => setSelectedRef(Number(event.target.value))}
                        disabled={!sites}
                    >
                        {sites?.sites.map((site) => (
                            <option key={site.referenceNumber} value={site.referenceNumber}>
                                {`#${site.referenceNumber} ${site.element} — U=${numberFormat(site.uIso, 4)} Å²`}
                            </option>
                        ))}
                    </select>
                </div>
                <div className="pca-control">
                    <label htmlFor="pca-grid">Grid</label>
                    <select id="pca-grid" value={grid} onChange={(event) => setGrid(Number(event.target.value))}>
                        {GRID_OPTIONS.map((value) => <option key={value} value={value}>{`${value}³`}</option>)}
                    </select>
                </div>
                <div className="pca-control">
                    <label htmlFor="pca-bw">Bandwidth</label>
                    <select id="pca-bw" value={bw} onChange={(event) => setBw(event.target.value)}>
                        {BW_OPTIONS.map((option) => <option key={option.value} value={option.value}>{option.label}</option>)}
                    </select>
                </div>
                <div className="pca-control">
                    <label htmlFor="pca-extent">Box (σ): {extent.toFixed(1)}</label>
                    <input
                        id="pca-extent" type="range" min="2" max="5" step="0.5"
                        value={extent} onChange={(event) => setExtent(Number(event.target.value))}
                    />
                </div>
                <div className="pca-control">
                    <label htmlFor="pca-prob">Ellipsoid: {Math.round(probability * 100)}%
                        <InfoBadge text="Enclosed-probability level for the thermal ellipsoid (50% is the crystallographic convention)." />
                    </label>
                    <input
                        id="pca-prob" type="range" min="0.1" max="0.99" step="0.01"
                        value={probability} onChange={(event) => setProbability(Number(event.target.value))}
                    />
                </div>
                <div className="pca-control">
                    <label htmlFor="pca-iso">Isosurface mass: {isoPercent}%
                        <InfoBadge text="The density isosurface enclosing this fraction of the cloud. Compare its shape to the harmonic ellipsoid." />
                    </label>
                    <input
                        id="pca-iso" type="range" min="1" max="99" step="1"
                        value={isoPercent} onChange={(event) => setIsoPercent(Number(event.target.value))}
                    />
                </div>
                <div className="pca-control">
                    <label htmlFor="pca-cmap">Colormap</label>
                    <select id="pca-cmap" value={colormap} onChange={(event) => setColormap(event.target.value)}>
                        {COLORMAP_NAMES.map((name) => <option key={name} value={name}>{name}</option>)}
                    </select>
                </div>
                <div className="pca-control pca-toggles">
                    <label><input type="checkbox" checked={showSurface} onChange={(event) => setShowSurface(event.target.checked)} /> Isosurface</label>
                    <label><input type="checkbox" checked={showEllipsoid} onChange={(event) => setShowEllipsoid(event.target.checked)} /> Ellipsoid</label>
                </div>
            </div>

            {noRun && <p className="pca-hint">Open a run folder (with an <code>.rmc6f</code> file) to view thermal ellipsoids.</p>}
            {sitesError && <p className="pca-error">{sitesError}</p>}

            <div className="pca-body">
                <div className="pca-viewport">
                    <div className="pca-canvas" ref={mountRef} />
                    {(loadingKde || loadingSites) && <div className="pca-spinner">Computing…</div>}
                    {kdeError && <div className="pca-error pca-overlay">{kdeError}</div>}
                    <div className="pca-legend">
                        <span><i style={{ background: '#d64545' }} /> PC1</span>
                        <span><i style={{ background: '#3fa34d' }} /> PC2</span>
                        <span><i style={{ background: '#3f7fd6' }} /> PC3</span>
                        <span><i style={{ background: '#ff5a5a' }} /> {Math.round(probability * 100)}% ellipsoid</span>
                    </div>
                </div>

                <div className="pca-side">
                    {selectedEllipsoid && (
                        <div className="pca-stats">
                            <h3>{selectedEllipsoid.element} site #{selectedEllipsoid.referenceNumber}
                                <span className="pca-count">{selectedEllipsoid.count.toLocaleString()} atoms</span>
                            </h3>
                            <table>
                                <tbody>
                                    <tr><th>U<sub>iso</sub> (Å²)</th><td>{numberFormat(selectedEllipsoid.uIso)}</td></tr>
                                    <tr><th>B<sub>iso</sub> (Å²)</th><td>{numberFormat(selectedEllipsoid.bIso, 3)}</td></tr>
                                    <tr><th>RMS axes (Å)</th><td>{selectedEllipsoid.rms.map((value) => numberFormat(value, 3)).join(', ')}</td></tr>
                                    <tr><th>Anisotropy</th><td>{numberFormat(selectedEllipsoid.anisotropy, 2)}{selectedEllipsoid.degenerate ? ' (degenerate)' : ''}</td></tr>
                                    <tr>
                                        <th>Non-Gaussianity
                                            <InfoBadge text="Mean excess kurtosis of the displacement cloud (0 = harmonic/Gaussian). Positive means a peaked, fat-tailed distribution — the KDE isosurface sits inside the ellipsoid, signalling anharmonic motion or split sites." />
                                        </th>
                                        <td>{numberFormat(selectedEllipsoid.nonGaussianity ?? kde?.nonGaussianity, 2)}</td>
                                    </tr>
                                    <tr><th>Eigenvalues (Å²)</th><td>{selectedEllipsoid.eigenvalues.map((value) => numberFormat(value, 4)).join(', ')}</td></tr>
                                </tbody>
                            </table>
                            {kde && (
                                <p className="pca-meta">
                                    Volume {kde.grid}³ · fit {kde.fitCount.toLocaleString()}/{kde.count.toLocaleString()} ·
                                    captured mass {numberFormat(kde.mass * 100, 1)}%
                                    {Number.isFinite(isoMassLevel) ? ` · iso @ ${isoPercent}% mass` : ''}
                                    {kde.browserPcaKde ? ' · browser' : ' · server'}
                                </p>
                            )}
                        </div>
                    )}

                    <div className="pca-projections">
                        {PROJECTION_META.map(({ key, label }) => (
                            <figure key={key}>
                                <canvas ref={projectionRefs[key]} className="pca-projection" />
                                <figcaption>{label}</figcaption>
                            </figure>
                        ))}
                    </div>
                </div>
            </div>
        </div>
    );
}
