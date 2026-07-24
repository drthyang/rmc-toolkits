// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Site-ellipsoids unit-cell picker: every reference site in the average unit
// cell drawn as its calculated thermal ellipsoid, colored by element, with the
// selected site highlighted — click a site to select it. Extracted from
// PcaKdePage so the Orientation page can reuse the same picker; the host page
// supplies the sites table and the selection callback.

import React, { useCallback, useEffect, useRef, useState } from 'react';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { DEFAULT_ELEMENT_COLOR } from '../atomColors';
import { downloadBlob, sanitizeFilename, saveCanvasAsPng } from '../figureExport';
import InfoBadge from './InfoBadge';
import SaveMenu from './SaveMenu';
import {
    CELL_AXIS_CSS,
    CELL_AXIS_LABELS,
    SELECTED_ATOM_COLOR,
    buildAxisTriad,
    buildCrystalAxes
} from './sceneAxes';
import './PcaKdePage.css';

const SAVE_OPTIONS = [
    { id: 'png', label: 'Standard PNG', hint: '1×' },
    { id: 'png3x', label: 'High quality PNG', hint: '3×' }
];

export default function SiteStructurePanel({ sites, selectedRef, onSelectSite, selectedEllipsoid, elementColors }) {
    // Crystallographic a/b/c gizmo at the cell origin — on by default, since
    // this is a structural view where axes belong.
    const [showCellAxes, setShowCellAxes] = useState(true);

    const structureMountRef = useRef(null);
    const structureSceneRef = useRef(null);
    const selectRef = useRef(null);
    // Keep the click handler pointing at the current callback without
    // re-creating the (once-mounted) scene.
    useEffect(() => { selectRef.current = onSelectSite; }, [onSelectSite]);

    // --- Scene: build once, keep a handle for updates. ------------------------
    useEffect(() => {
        const mount = structureMountRef.current;
        if (!mount) return undefined;
        const width = mount.clientWidth || 260;
        const height = mount.clientHeight || 260;

        const scene = new THREE.Scene();
        const camera = new THREE.PerspectiveCamera(40, width / height, 0.01, 1000);
        // preserveDrawingBuffer so the figure can be captured for PNG export at
        // any time (matches the main viewport's renderer).
        const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, preserveDrawingBuffer: true });
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
        const cellAxesGroup = new THREE.Group();
        const highlightGroup = new THREE.Group();
        scene.add(cellGroup, cellAxesGroup, bondsGroup, sitesGroup, highlightGroup);

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
            scene, camera, renderer, controls, sitesGroup, bondsGroup, cellGroup, cellAxesGroup, highlightGroup,
            framed: false, center: null, defaultRadius: 1
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
        const { sitesGroup, bondsGroup, cellGroup, cellAxesGroup, highlightGroup, camera, controls } = handle;

        const dispose = (group) => {
            while (group.children.length) {
                const child = group.children.pop();
                child.geometry?.dispose();
                // Axis-label sprites carry a CanvasTexture to release too.
                child.material?.map?.dispose();
                child.material?.dispose();
                group.remove(child);
            }
        };
        dispose(sitesGroup);
        dispose(bondsGroup);
        dispose(cellGroup);
        dispose(cellAxesGroup);
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

        // Crystallographic a/b/c gizmo at the cell origin corner (opt-in), sized
        // to a fraction of the shortest edge so it reads as a corner widget.
        if (showCellAxes) {
            const origin = new THREE.Vector3(...toCartesian([0, 0, 0]));
            const gizmoLength = 0.5 * edgeLength;
            buildCrystalAxes(origin, unit, gizmoLength, gizmoLength * 0.035)
                .forEach((obj) => cellAxesGroup.add(obj));
        }

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

        // Store the default framing so "Reset view" can restore it; frame the
        // cell once, then keep the user's orientation on later rebuilds.
        const span = Math.max(...unit.map((row) => Math.hypot(...row)));
        const defaultRadius = span * 1.7 || 1;
        handle.center = center;
        handle.defaultRadius = defaultRadius;
        controls.target.set(...center);
        if (!handle.framed) {
            camera.position.set(center[0] + defaultRadius, center[1] + defaultRadius * 0.7, center[2] + defaultRadius);
            camera.near = defaultRadius / 100;
            camera.far = defaultRadius * 100;
            camera.updateProjectionMatrix();
            handle.framed = true;
        }
        controls.update();
    }, [sites, selectedRef, selectedEllipsoid, elementColors, showCellAxes]);

    // Re-fit the canvas to its mount whenever sibling panels sharing the column
    // change height: content lands in stages, each collapsing the shared row a
    // little, and effects run after layout so the mount is settled here. The
    // ResizeObserver also catches window resizes, but this keeps the canvas
    // exact without depending on its timing.
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
    }, [sites, selectedRef, selectedEllipsoid]);

    // Re-apply the default framing (the rebuild stores the cell centre + a
    // fitting radius on the handle). No-op until framed.
    const resetStructureView = useCallback(() => {
        const handle = structureSceneRef.current;
        if (!handle || !handle.center) return;
        const { camera, controls, center, defaultRadius } = handle;
        const radius = defaultRadius || 1;
        controls.target.set(center[0], center[1], center[2]);
        camera.up.set(0, 1, 0);
        camera.position.set(center[0] + radius, center[1] + radius * 0.7, center[2] + radius);
        camera.near = radius / 100;
        camera.far = radius * 100;
        camera.updateProjectionMatrix();
        controls.update();
    }, []);

    // Export the figure as PNG — same 1× / 3× options as the main panels.
    const saveStructureView = useCallback(async (format) => {
        const handle = structureSceneRef.current;
        if (!handle) return;
        const { renderer, scene, camera } = handle;
        const name = 'PCA_Site_ellipsoids';
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
                    Site ellipsoids
                    <InfoBadge label="About the site ellipsoids" align="start">
                        <p>
                            Every reference site in the average unit cell, drawn as its
                            calculated thermal ellipsoid and colored by element. Click one to
                            select it; the selected site is highlighted. The <b>a b c</b>{' '}
                            button toggles the crystallographic-axis gizmo at the cell origin.
                        </p>
                    </InfoBadge>
                </span>
                <span className="panel-title-actions">
                    <button
                        type="button"
                        className={`pca-reset-view pca-cell-toggle ${showCellAxes ? 'is-active' : ''}`}
                        onClick={() => setShowCellAxes((value) => !value)}
                        aria-pressed={showCellAxes}
                        title="Show the crystallographic axes (a, b, c)"
                    >
                        <span style={{ color: CELL_AXIS_CSS[0] }}>a</span>
                        <span style={{ color: CELL_AXIS_CSS[1] }}>b</span>
                        <span style={{ color: CELL_AXIS_CSS[2] }}>c</span>
                    </button>
                    <button
                        type="button"
                        className="pca-reset-view"
                        onClick={resetStructureView}
                        disabled={!sites}
                        title="Reset the camera to the default view"
                    >
                        <svg width="12" height="12" viewBox="0 0 24 24" aria-hidden="true" fill="none" stroke="currentColor" strokeWidth="2.2" strokeLinecap="round" strokeLinejoin="round">
                            <path d="M3 12a9 9 0 1 0 9-9 9.75 9.75 0 0 0-6.74 2.74L3 8" />
                            <path d="M3 3v5h5" />
                        </svg>
                        Reset view
                    </button>
                    <SaveMenu
                        onSave={saveStructureView}
                        options={SAVE_OPTIONS}
                        label="Save"
                        align="right"
                        disabled={!sites}
                    />
                </span>
            </h3>
            <div className="pca-structure" ref={structureMountRef} />
            {sites?.elements?.length > 0 && (
                <div className="pca-legend">
                    {showCellAxes && (
                        <span className="pca-legend-item pca-legend-cellaxes">
                            {CELL_AXIS_LABELS.map((label, i) => (
                                <span key={label} className="pca-legend-cellaxis">
                                    <i className="pca-legend-swatch" style={{ background: CELL_AXIS_CSS[i] }} /> {label}
                                </span>
                            ))}
                            <span className="pca-legend-note">crystal axes</span>
                        </span>
                    )}
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
    );
}
