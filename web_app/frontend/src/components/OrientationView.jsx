// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Displacement-orientation sphere: the solid-angle histogram of a site's
// displacement *directions* (u = Δr/|Δr|), hex-binned on a Goldberg tiling and
// color-mapped on the unit sphere. Rendered as a tab of the PCA Ellipsoid main
// panel, sharing the page's site picker; the engine lives in
// workers/orientation.js (browser) and /api/pca/orientation (Flask), and this
// component only fetches, renders, and reads out.
//
// What to look for that the ellipsoid cannot show: discrete spots (hop sites),
// and a +u/−u imbalance (static off-centring / odd anharmonicity) — the map is
// never antipodally folded, and the asymmetry readout flags a real imbalance
// against its Poisson noise floor.

import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { COLORMAP_NAMES } from '../colormaps';
import {
    buildCellMesh,
    buildCellOutline,
    colorbarGradient,
    formatDirection,
    reliefFactors,
    vertexRadii
} from '../orientationSphere';
import { downloadBlob, sanitizeFilename, saveCanvasAsPng } from '../figureExport';
import InfoBadge from './InfoBadge';
import SaveMenu from './SaveMenu';
import './PcaKdePage.css';

const SAVE_OPTIONS = [
    { id: 'png', label: 'Standard PNG', hint: '1×' },
    { id: 'png3x', label: 'High quality PNG', hint: '3×' }
];

// Manual resolution choices (geodesic frequency ν → 10ν²+2 cells). 'auto' asks
// the engine for recommended_frequency, the ~12-points-per-cell guard against
// reading Poisson confetti as structure.
const FREQUENCY_OPTIONS = ['auto', 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24];

const WEIGHT_OPTIONS = [
    { value: 'count', label: 'Count' },
    { value: 'amplitude', label: '|Δr|' },
    { value: 'amplitude2', label: '|Δr|²' }
];

// Same axis palettes as the density view, so PC1/PC2/PC3 and a/b/c read as the
// same objects across the two tabs.
const TRIAD_COLORS = [0xd64545, 0x3fa34d, 0x3f7fd6];
const PC_CSS_COLORS = ['#d64545', '#3fa34d', '#3f7fd6'];
const CELL_AXIS_COLORS = [0xe0a419, 0x18a3a0, 0xb15ad8];
const CELL_AXIS_CSS = ['#e0a419', '#18a3a0', '#b15ad8'];
const CELL_AXIS_LABELS = ['a', 'b', 'c'];
const TRIAD_UP = new THREE.Vector3(0, 1, 0);

const numberFormat = (value, digits = 2) =>
    Number.isFinite(value) ? value.toFixed(digits) : '—';

// Thin axis rods through the sphere centre along ±dir (both signs, since a
// direction map has no preferred end of an axis), one mesh per sign.
const buildAxisRods = (axes, colors, radius, length) => {
    const rods = [];
    axes.forEach((axis, index) => {
        const dir = new THREE.Vector3(axis[0], axis[1], axis[2]);
        if (dir.lengthSq() < 1e-12) return;
        dir.normalize();
        for (const sign of [1, -1]) {
            const rod = new THREE.Mesh(
                new THREE.CylinderGeometry(radius, radius, length, 10),
                new THREE.MeshBasicMaterial({ color: colors[index], transparent: true, opacity: 0.9 })
            );
            rod.position.copy(dir.clone().multiplyScalar(sign * length / 2));
            rod.quaternion.setFromUnitVectors(TRIAD_UP, dir.clone().multiplyScalar(sign));
            rods.push(rod);
        }
    });
    return rods;
};

export default function OrientationView({
    requestPca,
    ready,
    selectedRef,
    selectedEllipsoid,
    clusterThreshold,
    unitCell
}) {
    const [result, setResult] = useState(null);
    const [error, setError] = useState(null);
    const [loading, setLoading] = useState(false);

    const [frequency, setFrequency] = useState('auto');
    const [weight, setWeight] = useState('count');
    const [frame, setFrame] = useState('cartesian');
    const [smoothing, setSmoothing] = useState(0);
    const [minQuantile, setMinQuantile] = useState(0);
    const [colormap, setColormap] = useState('viridis');
    // Amplitude relief: 0 keeps the unit sphere; higher values bulge each cell
    // radially by its mean |Δr| relative to the site average, so shape (how far
    // atoms move this way) and color (how often) carry independent information.
    const [relief, setRelief] = useState(0.5);
    const [showOutline, setShowOutline] = useState(true);
    const [showAxes, setShowAxes] = useState(true);
    // Hovered cell readout: { cell, x, y } in canvas-local coordinates.
    const [hover, setHover] = useState(null);

    const mountRef = useRef(null);
    const sceneRef = useRef(null);
    const resultRef = useRef(null);
    resultRef.current = result;

    // --- Fetch the histogram whenever the site or an option changes. ----------
    useEffect(() => {
        let cancelled = false;
        const load = async () => {
            if (!ready || selectedRef == null) { setResult(null); return; }
            setLoading(true);
            setError(null);
            try {
                const data = await requestPca('orientation', {
                    referenceNumber: selectedRef,
                    frequency: frequency === 'auto' ? null : frequency,
                    weight,
                    frame,
                    smoothing,
                    minAmplitudeQuantile: minQuantile,
                    clusterThreshold,
                    geometry: true
                });
                if (!cancelled) {
                    // A hover index from the previous tiling may be out of range
                    // for the new one (fewer cells) — clear it with the result.
                    setHover(null);
                    setResult(data);
                }
            } catch (requestError) {
                if (cancelled) return;
                // A cluster-distance change can transiently request a re-numbered
                // site (same self-healing race as the KDE volume); keep the current
                // sphere rather than surfacing it.
                if (/Unknown reference number/i.test(requestError.message || '')) return;
                setResult(null);
                setError(requestError.message);
            } finally {
                if (!cancelled) setLoading(false);
            }
        };
        load();
        return () => { cancelled = true; };
    }, [requestPca, ready, selectedRef, frequency, weight, frame, smoothing, minQuantile, clusterThreshold]);

    // --- Three.js scene: build once. ------------------------------------------
    useEffect(() => {
        const mount = mountRef.current;
        if (!mount) return undefined;
        const width = mount.clientWidth || 640;
        const height = mount.clientHeight || 520;

        const scene = new THREE.Scene();
        const camera = new THREE.PerspectiveCamera(45, width / height, 0.01, 100);
        camera.position.set(1.9, 1.4, 1.9);
        const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, preserveDrawingBuffer: true });
        renderer.setPixelRatio(window.devicePixelRatio || 1);
        renderer.setSize(width, height);
        mount.appendChild(renderer.domElement);

        const controls = new OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.dampingFactor = 0.12;
        controls.enablePan = false;

        const meshGroup = new THREE.Group();
        const outlineGroup = new THREE.Group();
        const axesGroup = new THREE.Group();
        scene.add(meshGroup, outlineGroup, axesGroup);

        // Raycast picking for the per-cell hover readout. Only re-renders React
        // when the hovered cell changes, so orbiting stays smooth.
        const raycaster = new THREE.Raycaster();
        const pointer = new THREE.Vector2();
        let lastCell = null;
        const onPointerMove = (event) => {
            const rect = renderer.domElement.getBoundingClientRect();
            pointer.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
            pointer.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
            raycaster.setFromCamera(pointer, camera);
            const hit = raycaster.intersectObjects(meshGroup.children, false)[0];
            const cell = hit ? hit.object.userData.triangleCell[hit.faceIndex] : null;
            if (cell !== lastCell) {
                lastCell = cell;
                setHover(cell == null ? null : { cell, x: event.clientX - rect.left, y: event.clientY - rect.top });
            } else if (cell != null) {
                setHover({ cell, x: event.clientX - rect.left, y: event.clientY - rect.top });
            }
        };
        const onPointerLeave = () => { lastCell = null; setHover(null); };
        renderer.domElement.addEventListener('pointermove', onPointerMove);
        renderer.domElement.addEventListener('pointerleave', onPointerLeave);

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

        sceneRef.current = { scene, camera, renderer, controls, meshGroup, outlineGroup, axesGroup };

        return () => {
            cancelAnimationFrame(animationId);
            resizeObserver.disconnect();
            renderer.domElement.removeEventListener('pointermove', onPointerMove);
            renderer.domElement.removeEventListener('pointerleave', onPointerLeave);
            controls.dispose();
            // Unlike the density canvas (kept mounted via CSS), this component
            // unmounts on every tab switch — release the last-built geometries
            // and force the GL context loss, or repeated Density↔Orientation
            // toggles pile up GPU buffers and live WebGL contexts until GC.
            [meshGroup, outlineGroup, axesGroup].forEach((group) => {
                while (group.children.length) {
                    const child = group.children.pop();
                    child.geometry?.dispose();
                    child.material?.dispose();
                    group.remove(child);
                }
            });
            renderer.dispose();
            renderer.forceContextLoss();
            if (renderer.domElement.parentNode === mount) mount.removeChild(renderer.domElement);
            sceneRef.current = null;
        };
    }, []);

    // --- Rebuild the sphere mesh when the histogram or display options change. -
    useEffect(() => {
        const handle = sceneRef.current;
        if (!handle) return;
        const { meshGroup, outlineGroup, axesGroup } = handle;
        const dispose = (group) => {
            while (group.children.length) {
                const child = group.children.pop();
                child.geometry?.dispose();
                child.material?.dispose();
                group.remove(child);
            }
        };
        dispose(meshGroup);
        dispose(outlineGroup);
        dispose(axesGroup);
        if (!result?.polygons) return;

        // Amplitude relief: per-cell radial factors from the mean |Δr| of each
        // cell's movers, averaged at shared polygon vertices so the surface
        // stays crack-free while every cell keeps its flat color.
        const radii = relief > 0 && result.cellMeanAmplitude
            ? vertexRadii(result.polygons, reliefFactors(result.cellMeanAmplitude, result.meanAmplitude, relief))
            : null;

        const mesh = buildCellMesh(result.polygons, result.enhancement, colormap, result.vmax, radii);
        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute('position', new THREE.BufferAttribute(mesh.positions, 3));
        geometry.setAttribute('color', new THREE.BufferAttribute(mesh.colors, 3));
        geometry.computeVertexNormals();
        const material = new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.FrontSide });
        const sphere = new THREE.Mesh(geometry, material);
        sphere.userData.triangleCell = mesh.triangleCell;
        meshGroup.add(sphere);

        if (showOutline) {
            const outline = buildCellOutline(result.polygons, 1.002, radii);
            const outlineGeometry = new THREE.BufferGeometry();
            outlineGeometry.setAttribute('position', new THREE.BufferAttribute(outline, 3));
            outlineGroup.add(new THREE.LineSegments(
                outlineGeometry,
                new THREE.LineBasicMaterial({ color: 0x10151c, transparent: true, opacity: 0.35 })
            ));
        }

        if (showAxes) {
            // Principal axes: the engine reports the PCA rotation it used, so in
            // the PCA frame the axes are the view's own x/y/z (identity).
            const pcAxes = frame === 'pca'
                ? [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                : result.pcaAxes || [];
            buildAxisRods(pcAxes, TRIAD_COLORS, 0.006, 2.7).forEach((rod) => axesGroup.add(rod));
            // Crystallographic a/b/c (normalized directions) in the Cartesian
            // frame only — in the PCA frame they would need the same rotation and
            // read as clutter over an already-rotated map.
            if (frame === 'cartesian' && unitCell) {
                buildAxisRods(unitCell, CELL_AXIS_COLORS, 0.0045, 2.35).forEach((rod) => axesGroup.add(rod));
            }
        }
    }, [result, colormap, relief, showOutline, showAxes, frame, unitCell]);

    const resetView = useCallback(() => {
        const handle = sceneRef.current;
        if (!handle) return;
        const { camera, controls } = handle;
        controls.target.set(0, 0, 0);
        camera.up.set(0, 1, 0);
        camera.position.set(1.9, 1.4, 1.9);
        camera.updateProjectionMatrix();
        controls.update();
    }, []);

    const saveView = useCallback(async (format) => {
        const handle = sceneRef.current;
        if (!handle) return;
        const { renderer, scene, camera } = handle;
        const name = selectedEllipsoid
            ? `Orientation_${selectedEllipsoid.element}_site${selectedEllipsoid.referenceNumber}`
            : 'Orientation_sphere';
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

    const gradient = useMemo(() => colorbarGradient(colormap), [colormap]);

    // A real inversion asymmetry: well above what Poisson noise alone produces.
    const asymmetrySignificant = result
        ? result.antipodalAsymmetry > 3 * result.antipodalAsymmetryNull
        : false;

    // Bounds-guarded: a raycast fired between a resolution change and the hover
    // reset could carry an index from the previous, larger tiling.
    const hoverCell = hover && result && hover.cell < result.cellCount ? hover.cell : null;

    return (
        <div className="orient-view">
            <div className="pca-controls orient-controls">
                <div className="control-group" role="group" aria-label="Orientation histogram options">
                    <label className="control">
                        <span className="control-name">
                            Resolution
                            <InfoBadge label="About the sphere resolution">
                                <p>
                                    Geodesic frequency ν of the hex tiling (10ν² + 2 cells — hexagons
                                    plus the 12 pentagons every hexagonal tiling of a sphere must
                                    contain). Auto targets ~12 displacements per cell, the guard
                                    against reading Poisson noise as structure.
                                </p>
                            </InfoBadge>
                        </span>
                        <select
                            value={frequency}
                            onChange={(event) => setFrequency(event.target.value === 'auto' ? 'auto' : Number(event.target.value))}
                            aria-label="Sphere resolution"
                        >
                            {FREQUENCY_OPTIONS.map((option) => (
                                <option key={option} value={option}>
                                    {option === 'auto'
                                        ? `Auto${result && frequency === 'auto' ? ` (ν=${result.frequency})` : ''}`
                                        : `ν=${option} (${10 * option * option + 2})`}
                                </option>
                            ))}
                        </select>
                    </label>
                    <label className="control">
                        <span className="control-name">
                            Weight
                            <InfoBadge label="About the weighting">
                                <p>
                                    Count: every atom votes once — the orientation distribution
                                    proper. |Δr| / |Δr|²: longer displacements vote more; the |Δr|²
                                    map is the angular decomposition of the mean-square displacement
                                    the U tensor summarizes.
                                </p>
                            </InfoBadge>
                        </span>
                        <select value={weight} onChange={(event) => setWeight(event.target.value)} aria-label="Cell weighting">
                            {WEIGHT_OPTIONS.map((option) => (
                                <option key={option.value} value={option.value}>{option.label}</option>
                            ))}
                        </select>
                    </label>
                    <label className="control">
                        <span className="control-name">
                            Min |Δr|
                            <InfoBadge label="About the amplitude cutoff">
                                <p>
                                    Drops the shortest displacements before binning (quantile of
                                    |Δr|). A near-zero displacement has a direction dominated by
                                    noise, which dilutes a real pattern toward uniform.
                                </p>
                            </InfoBadge>
                        </span>
                        <input
                            type="range" min="0" max="0.5" step="0.05"
                            value={minQuantile}
                            onChange={(event) => setMinQuantile(Number(event.target.value))}
                            aria-label="Minimum displacement quantile"
                        />
                        <span className="control-value">{Math.round(minQuantile * 100)}%</span>
                    </label>
                    <label className="control">
                        <span className="control-name">Smoothing</span>
                        <input
                            type="range" min="0" max="4" step="1"
                            value={smoothing}
                            onChange={(event) => setSmoothing(Number(event.target.value))}
                            aria-label="Neighbour smoothing passes"
                        />
                        <span className="control-value">{smoothing}×</span>
                    </label>
                    <label className="control">
                        <span className="control-name">
                            Relief
                            <InfoBadge label="About the amplitude relief">
                                <p>
                                    Bulges the sphere radially by each cell's mean |Δr| relative to
                                    the site average — directions where atoms move farther stick
                                    out, shorter ones dent in. Color (how often) and shape (how far)
                                    then carry independent information. 0% keeps a perfect sphere.
                                </p>
                            </InfoBadge>
                        </span>
                        <input
                            type="range" min="0" max="1" step="0.05"
                            value={relief}
                            onChange={(event) => setRelief(Number(event.target.value))}
                            aria-label="Amplitude relief strength"
                        />
                        <span className="control-value">{Math.round(relief * 100)}%</span>
                    </label>
                    <label className="control">
                        <span className="control-name">Colormap</span>
                        <select value={colormap} onChange={(event) => setColormap(event.target.value)} aria-label="Sphere colormap">
                            {COLORMAP_NAMES.map((name) => <option key={name} value={name}>{name}</option>)}
                        </select>
                    </label>
                    <label className="control switch">
                        <span className="control-name">Cell borders</span>
                        <input type="checkbox" checked={showOutline} onChange={(event) => setShowOutline(event.target.checked)} aria-label="Show cell borders" />
                        <i className="switch-track" aria-hidden="true" />
                    </label>
                    <label className="control switch">
                        <span className="control-name">Axes</span>
                        <input type="checkbox" checked={showAxes} onChange={(event) => setShowAxes(event.target.checked)} aria-label="Show axis rods" />
                        <i className="switch-track" aria-hidden="true" />
                    </label>
                </div>
            </div>

            <div className="pca-canvas orient-canvas" ref={mountRef}>
                {loading && <div className="pca-badge">Computing…</div>}
                {error && <div className="pca-badge is-error">{error}</div>}
                <div className="pca-view-controls pca-view-controls--left">
                    <div className="pca-view-group" role="group" aria-label="Direction frame">
                        <button
                            type="button"
                            className={`pca-view-btn ${frame === 'cartesian' ? 'is-active' : ''}`}
                            onClick={() => setFrame('cartesian')}
                            aria-pressed={frame === 'cartesian'}
                            title="Crystal Cartesian frame"
                        >
                            Crystal
                        </button>
                        <button
                            type="button"
                            className={`pca-view-btn ${frame === 'pca' ? 'is-active' : ''}`}
                            onClick={() => setFrame('pca')}
                            aria-pressed={frame === 'pca'}
                            title="This site's principal-axis frame (PC1 = x)"
                        >
                            PCA
                        </button>
                    </div>
                </div>
                <div className="pca-view-controls pca-view-controls--right">
                    <div className="pca-view-group" role="group" aria-label="View actions">
                        <button
                            type="button"
                            className="pca-view-btn"
                            onClick={resetView}
                            title="Reset the camera to the default view"
                        >
                            Reset
                        </button>
                    </div>
                    <SaveMenu
                        onSave={saveView}
                        options={SAVE_OPTIONS}
                        label="Save"
                        align="right"
                        disabled={!result}
                    />
                </div>
                {hoverCell != null && result && (
                    <div className="orient-tooltip" style={{ left: hover.x + 14, top: hover.y + 12 }}>
                        <div>{formatDirection(result.centers[hoverCell])}</div>
                        <div>{numberFormat(result.enhancement[hoverCell], 2)}× isotropic</div>
                        <div>
                            {result.counts[hoverCell]} atoms · z = {numberFormat(result.zScore[hoverCell], 1)}
                        </div>
                        {result.cellMeanAmplitude?.[hoverCell] > 0 && (
                            <div>⟨|Δr|⟩ = {numberFormat(result.cellMeanAmplitude[hoverCell], 3)} Å</div>
                        )}
                    </div>
                )}
            </div>

            {result && (
                <>
                    <div className="orient-colorbar-row">
                        <span className="orient-colorbar-label">0</span>
                        <div className="orient-colorbar" style={{ background: gradient }}>
                            {result.vmax > 1 && (
                                <i
                                    className="orient-colorbar-marker"
                                    style={{ left: `${(100 / result.vmax).toFixed(2)}%` }}
                                    title="1× = isotropic"
                                />
                            )}
                        </div>
                        <span className="orient-colorbar-label">{numberFormat(result.vmax, 1)}× isotropic</span>
                    </div>
                    <div className="orient-summary">
                        <span className="orient-stat">
                            <b>{result.cellCount}</b> cells (ν={result.frequency}, 12 pentagons)
                        </span>
                        <span className="orient-stat">
                            <b>{result.usedPoints.toLocaleString()}</b>/{result.totalPoints.toLocaleString()} vectors
                            {result.rejectedPoints > 0 ? ` (${result.rejectedPoints} below cutoff)` : ''}
                        </span>
                        <span className="orient-stat">
                            peak <b>{numberFormat(result.peakEnhancement, 2)}×</b> at {formatDirection(result.peakDirection)}
                            {' '}(z = {numberFormat(result.peakZScore, 1)})
                        </span>
                        <span className="orient-stat">
                            anisotropy <b>{numberFormat(result.orientationAnisotropy, 2)}</b>
                            <InfoBadge label="About the orientation anisotropy" align="end">
                                <p>
                                    3λ₁ − 1 of the orientation tensor ⟨u uᵀ⟩: 0 for an isotropic
                                    direction distribution, 2 for a perfect single axis. Resolution
                                    independent (computed from the vectors, not the bins).
                                </p>
                            </InfoBadge>
                        </span>
                        <span className={`orient-stat ${asymmetrySignificant ? 'is-flagged' : ''}`}>
                            ± asymmetry <b>{numberFormat(result.antipodalAsymmetry, 2)}</b>
                            {' '}<span className="orient-stat-null">(noise floor {numberFormat(result.antipodalAsymmetryNull, 2)})</span>
                            <InfoBadge label="About the antipodal asymmetry" align="end">
                                <p>
                                    The +u vs −u imbalance, Σ|n(u) − n(−u)| / N over antipodal cell
                                    pairs: 0 for an inversion-symmetric cloud, 1 for a fully
                                    one-sided one. The thermal ellipsoid is blind to this — a value
                                    well above the Poisson noise floor is real off-centring or
                                    odd-order anharmonicity.
                                </p>
                            </InfoBadge>
                        </span>
                        <span className="orient-stat">
                            map significance <b>{numberFormat(result.significance, 1)}σ</b>
                            <InfoBadge label="About the map significance" align="end">
                                <p>
                                    RMS of the per-cell Poisson z-scores against the isotropic null:
                                    ≈1 means the pattern is consistent with pure counting noise;
                                    well above 1 means real directional structure.
                                </p>
                            </InfoBadge>
                        </span>
                    </div>
                    <div className="pca-legend">
                        {frame === 'pca' ? (
                            <>
                                <span className="pca-legend-item"><i className="pca-legend-swatch" style={{ background: PC_CSS_COLORS[0] }} /> PC1</span>
                                <span className="pca-legend-item"><i className="pca-legend-swatch" style={{ background: PC_CSS_COLORS[1] }} /> PC2</span>
                                <span className="pca-legend-item"><i className="pca-legend-swatch" style={{ background: PC_CSS_COLORS[2] }} /> PC3</span>
                            </>
                        ) : (
                            <>
                                <span className="pca-legend-item"><i className="pca-legend-swatch" style={{ background: PC_CSS_COLORS[0] }} /> PC1</span>
                                <span className="pca-legend-item"><i className="pca-legend-swatch" style={{ background: PC_CSS_COLORS[1] }} /> PC2</span>
                                <span className="pca-legend-item"><i className="pca-legend-swatch" style={{ background: PC_CSS_COLORS[2] }} /> PC3</span>
                                {unitCell && CELL_AXIS_LABELS.map((label, i) => (
                                    <span key={label} className="pca-legend-item">
                                        <i className="pca-legend-swatch" style={{ background: CELL_AXIS_CSS[i] }} /> {label}
                                    </span>
                                ))}
                            </>
                        )}
                        <span className="pca-legend-note">
                            {result.weight === 'count' ? 'direction distribution' : `weighted by ${result.weight === 'amplitude' ? '|Δr|' : '|Δr|²'}`}
                            {result.smoothing > 0 ? ` · smoothed ${result.smoothing}×` : ''}
                            {result.browserOrientation ? ' · browser' : ' · server'}
                        </span>
                    </div>
                </>
            )}
        </div>
    );
}
