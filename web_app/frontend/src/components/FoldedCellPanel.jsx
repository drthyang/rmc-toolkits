// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Folded unit cell with the analysis' detected bonds drawn over it. The atom
// cloud is the same one the Atomic Density page shows — every atom in the
// supercell folded back into a single cell, so the spread around each site is
// the thermal cloud itself rather than a fitted ellipsoid. Bonds come from the
// triplet windows and are drawn as thin transparent lines so a full network
// reads as a framework instead of hiding the cloud behind it.

import React, { useCallback, useEffect, useRef, useState } from 'react';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { DEFAULT_ELEMENT_COLOR } from '../atomColors';
import { downloadBlob, sanitizeFilename, saveCanvasAsPng } from '../figureExport';
import InfoBadge from './InfoBadge';
import SaveMenu from './SaveMenu';
import { buildCrystalAxes } from './sceneAxes';
import './PcaKdePage.css';

const SAVE_OPTIONS = [
    { id: 'png', label: 'Standard PNG', hint: '1×' },
    { id: 'png3x', label: 'High quality PNG', hint: '3×' }
];

// The cloud is one point per atom; a full supercell is large enough that the
// draw call, not the parse, is the cost. Sampling keeps the panel interactive
// on 10×10×10 models while leaving the shape of every site's cloud intact.
const MAX_CLOUD_POINTS = 120000;

// Unit-cell vectors = supercell vectors divided by the supercell counts.
const unitVectorsOf = (source) => {
    if (!source?.latticeVectors || !source?.supercell) return null;
    return source.latticeVectors.map((row, i) => row.map((v) => v / (source.supercell[i] || 1)));
};

const FoldedCellPanel = ({
    structure,
    sites,
    bondSets = null,
    elementColors = {},
    loading = false,
    title = 'Folded unit cell'
}) => {
    const [showCellAxes, setShowCellAxes] = useState(true);
    const mountRef = useRef(null);
    const sceneRef = useRef(null);

    // Scene, camera and controls are built once and kept; only the contents of
    // the groups below are rebuilt when the data changes.
    useEffect(() => {
        const mount = mountRef.current;
        if (!mount) return undefined;

        const scene = new THREE.Scene();
        const camera = new THREE.PerspectiveCamera(45, 1, 0.01, 5000);
        const renderer = new THREE.WebGLRenderer({
            antialias: true,
            alpha: true,
            preserveDrawingBuffer: true
        });
        renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        mount.replaceChildren(renderer.domElement);

        const controls = new OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.dampingFactor = 0.08;

        const cloudGroup = new THREE.Group();
        const bondsGroup = new THREE.Group();
        const cellGroup = new THREE.Group();
        const cellAxesGroup = new THREE.Group();
        scene.add(cloudGroup, bondsGroup, cellGroup, cellAxesGroup);

        let frame = 0;
        const renderLoop = () => {
            controls.update();
            renderer.render(scene, camera);
            frame = requestAnimationFrame(renderLoop);
        };
        renderLoop();

        const resize = () => {
            const width = mount.clientWidth;
            const height = mount.clientHeight;
            if (!width || !height) return;
            camera.aspect = width / height;
            camera.updateProjectionMatrix();
            // Let setSize style the canvas too: .pca-structure gives it no CSS
            // size, so without the style it would lay out at attribute size
            // (CSS × pixelRatio) and overflow the panel on retina displays.
            renderer.setSize(width, height);
        };
        resize();
        const observer = new ResizeObserver(resize);
        observer.observe(mount);

        sceneRef.current = {
            scene, camera, renderer, controls, cloudGroup, bondsGroup, cellGroup, cellAxesGroup, resize
        };

        return () => {
            cancelAnimationFrame(frame);
            observer.disconnect();
            controls.dispose();
            renderer.dispose();
            mount.replaceChildren();
            sceneRef.current = null;
        };
    }, []);

    // Rebuild the cloud, the cell box and the bond network.
    useEffect(() => {
        const handle = sceneRef.current;
        if (!handle) return;
        const { cloudGroup, bondsGroup, cellGroup, cellAxesGroup, camera, controls } = handle;

        const dispose = (group) => {
            while (group.children.length) {
                const child = group.children.pop();
                child.geometry?.dispose();
                child.material?.map?.dispose();
                child.material?.dispose();
                group.remove(child);
            }
        };
        [cloudGroup, bondsGroup, cellGroup, cellAxesGroup].forEach(dispose);

        const unit = unitVectorsOf(structure) || unitVectorsOf(sites);
        if (!unit) return;
        const toCartesian = (frac) => [0, 1, 2].map(
            (axis) => frac[0] * unit[0][axis] + frac[1] * unit[1][axis] + frac[2] * unit[2][axis]
        );
        const edgeLength = Math.min(...unit.map((row) => Math.hypot(...row)));
        const center = new THREE.Vector3(...toCartesian([0.5, 0.5, 0.5]));

        // --- Unit-cell wireframe: the 12 edges of the parallelepiped.
        const corners = [];
        for (const a of [0, 1]) {
            for (const b of [0, 1]) {
                for (const c of [0, 1]) corners.push(toCartesian([a, b, c]));
            }
        }
        const cellPoints = [];
        for (let i = 0; i < 8; i += 1) {
            for (let j = i + 1; j < 8; j += 1) {
                const gi = [(i >> 2) & 1, (i >> 1) & 1, i & 1];
                const gj = [(j >> 2) & 1, (j >> 1) & 1, j & 1];
                if (gi.reduce((s, v, k) => s + (v !== gj[k] ? 1 : 0), 0) === 1) {
                    cellPoints.push(new THREE.Vector3(...corners[i]), new THREE.Vector3(...corners[j]));
                }
            }
        }
        cellGroup.add(new THREE.LineSegments(
            new THREE.BufferGeometry().setFromPoints(cellPoints),
            new THREE.LineBasicMaterial({ color: 0x8a97a8, transparent: true, opacity: 0.4 })
        ));

        if (showCellAxes) {
            const origin = new THREE.Vector3(...toCartesian([0, 0, 0]));
            const gizmoLength = 0.5 * edgeLength;
            buildCrystalAxes(origin, unit, gizmoLength, gizmoLength * 0.035)
                .forEach((obj) => cellAxesGroup.add(obj));
        }

        // --- Atom cloud, one THREE.Points per element so each keeps its color.
        const allPoints = structure?.points || [];
        const stride = Math.max(1, Math.ceil(allPoints.length / MAX_CLOUD_POINTS));
        const byElement = new Map();
        for (let i = 0; i < allPoints.length; i += stride) {
            const point = allPoints[i];
            if (!byElement.has(point.element)) byElement.set(point.element, []);
            byElement.get(point.element).push(point);
        }
        byElement.forEach((elementPoints, element) => {
            const positions = new Float32Array(elementPoints.length * 3);
            elementPoints.forEach((point, index) => {
                const cart = toCartesian([point.x, point.y, point.z]);
                positions[index * 3] = cart[0];
                positions[index * 3 + 1] = cart[1];
                positions[index * 3 + 2] = cart[2];
            });
            const geometry = new THREE.BufferGeometry();
            geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
            cloudGroup.add(new THREE.Points(geometry, new THREE.PointsMaterial({
                color: elementColors[element] || DEFAULT_ELEMENT_COLOR,
                size: 0.02 * edgeLength,
                sizeAttenuation: true,
                transparent: true,
                opacity: 0.55
            })));
        });

        // --- Detected bonds between the average sites, periodic images included
        // (a bond may reach a neighbouring cell's image just outside the box —
        // that is the real coordination). Thin transparent lines, so the network
        // sits over the cloud rather than hiding it.
        if (bondSets?.length && sites?.sites?.length) {
            const positions = sites.sites.map((site) => ({
                ref: site.referenceNumber,
                el: site.element,
                frac: site.siteFractional,
                pos: new THREE.Vector3(...toCartesian(site.siteFractional))
            }));
            bondSets.forEach((set) => {
                const [elementA, elementB] = set.elements;
                const sameElement = elementA === elementB;
                const [rmin, rmax] = set.window;
                const segments = [];
                positions.forEach((from) => {
                    if (from.el !== elementA) return;
                    positions.forEach((to) => {
                        if (to.el !== elementB) return;
                        for (let mx = -1; mx <= 1; mx += 1) {
                            for (let my = -1; my <= 1; my += 1) {
                                for (let mz = -1; mz <= 1; mz += 1) {
                                    const inCell = mx === 0 && my === 0 && mz === 0;
                                    if (inCell && to.ref === from.ref) continue;
                                    // Same-element in-cell pairs would double-draw.
                                    if (sameElement && inCell && to.ref < from.ref) continue;
                                    const delta = [
                                        to.frac[0] + mx - from.frac[0],
                                        to.frac[1] + my - from.frac[1],
                                        to.frac[2] + mz - from.frac[2]
                                    ];
                                    const cart = toCartesian(delta);
                                    const length = Math.hypot(...cart);
                                    if (length < rmin || length > rmax || length === 0) continue;
                                    segments.push(
                                        from.pos.clone(),
                                        from.pos.clone().add(new THREE.Vector3(...cart))
                                    );
                                }
                            }
                        }
                    });
                });
                if (segments.length) {
                    bondsGroup.add(new THREE.LineSegments(
                        new THREE.BufferGeometry().setFromPoints(segments),
                        new THREE.LineBasicMaterial({
                            color: set.color,
                            transparent: true,
                            opacity: 0.45
                        })
                    ));
                }
            });
        }

        // Frame the cell once per rebuild: the camera keeps whatever the user
        // dragged it to only while the data is unchanged.
        const radius = Math.max(...corners.map(
            (corner) => new THREE.Vector3(...corner).distanceTo(center)
        ), edgeLength);
        controls.target.copy(center);
        camera.position.copy(center).add(new THREE.Vector3(radius * 1.5, radius * 1.1, radius * 1.6));
        camera.near = radius / 100;
        camera.far = radius * 100;
        camera.updateProjectionMatrix();
        controls.update();
        // Snapshot this framing as what "Reset view" restores — without it,
        // controls.reset() would return the camera to its construction state
        // at the origin.
        controls.saveState();
        handle.resize();
    }, [structure, sites, bondSets, elementColors, showCellAxes]);

    const resetView = useCallback(() => {
        const handle = sceneRef.current;
        if (!handle) return;
        handle.controls.reset();
    }, []);

    // Same 1× / 3× options as the other 3D panels: 3× re-renders at a higher
    // pixel ratio, then restores the on-screen one.
    const saveFigure = useCallback(async (format) => {
        const handle = sceneRef.current;
        if (!handle) return;
        const { renderer, scene, camera } = handle;
        const name = 'Folded_Unit_Cell';
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
    }, []);

    return (
        <div className="pca-panel pca-unitcell-panel">
            <h3>
                <span className="panel-title-label">
                    {title}
                    <InfoBadge label="About the folded unit cell" align="start">
                        <p>
                            Every atom in the supercell folded back into one unit cell and
                            colored by element — the same view as the Atomic Density page.
                            The spread around each site is the thermal cloud itself, not a
                            fitted ellipsoid.
                        </p>
                        {bondSets?.length ? (
                            <p>
                                The lines are the analysis' <em>detected bonds</em>: average
                                site pairs whose distance falls inside the bond window,
                                periodic images included, so a bond may reach an image just
                                outside the box.
                            </p>
                        ) : null}
                    </InfoBadge>
                </span>
                <span className="panel-title-actions">
                    {loading && <span className="panel-title-count">Loading…</span>}
                    <button
                        type="button"
                        className={`pca-reset-view pca-cell-toggle ${showCellAxes ? 'is-active' : ''}`}
                        onClick={() => setShowCellAxes((current) => !current)}
                        aria-label="Show the crystallographic axes (a, b, c)"
                    >
                        a b c
                    </button>
                    <button
                        type="button"
                        className="pca-reset-view"
                        onClick={resetView}
                        aria-label="Reset the camera to the default view"
                    >
                        Reset view
                    </button>
                    <SaveMenu onSave={saveFigure} options={SAVE_OPTIONS} label="Save" align="right" />
                </span>
            </h3>
            <div className="pca-structure" ref={mountRef} />
            {Object.keys(elementColors).length > 0 && (
                <div className="atom-legend" aria-label="Atom colors by element">
                    {Object.entries(elementColors).map(([element, color]) => (
                        <span key={element} className="atom-legend-item">
                            <span className="atom-legend-swatch" style={{ background: color }} />
                            {element}
                        </span>
                    ))}
                </div>
            )}
        </div>
    );
};

export default FoldedCellPanel;
