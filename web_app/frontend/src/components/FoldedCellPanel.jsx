// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Folded unit cell with the analysis' detected bonds drawn over it. The atom
// cloud is the same one the Atomic Density page shows — every atom in the
// supercell folded back into a single cell, so the spread around each site is
// the thermal cloud itself rather than a fitted ellipsoid. Bonds come from the
// triplet windows (bondNetwork.js) and are drawn as shaded tubes in their
// bond-role colours, one instanced mesh (one draw call) per window; the
// clouds of elements outside the bonds fade, so the selected triplet stands
// out. A large cell draws at most MAX_DRAWN_BONDS, those nearest the centre,
// and the legend says so. The element key sits inside the canvas, as a pill
// at the bottom left.
//
// `highlight` (angleHighlight.js highlightForBin: the triplets of the angle
// bin under the pointer) draws those triplets' bonds as thicker tubes of
// their own and dims the other bonds; the legend names the bin. It has its
// own group, so a hover never rebuilds the cloud.

import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { DEFAULT_ELEMENT_COLOR, elementColorNote } from '../atomColors';
import { downloadBlob, sanitizeFilename, saveCanvasAsPng } from '../figureExport';
import { Card, CardHeader, CardMeta, ToolButton } from '../ui';
import InfoBadge from '../ui/InfoBadge';
import SaveMenu from '../ui/SaveMenu';
import { buildCrystalAxes } from './sceneAxes';
import { bondRadius, cloudOpacity, foldedCellNetwork, MAX_DRAWN_BONDS, unitVectorsOf } from '../bondNetwork';
import { highlightLabel } from '../angleHighlight';
import useCopilotViewer from '../useCopilotViewer';
import './UnitCellPanel.css';

const SAVE_OPTIONS = [
    { id: 'png', label: 'Standard PNG', hint: '1×' },
    { id: 'png3x', label: 'High quality PNG', hint: '3×' }
];

// The cloud is one point per atom; a full supercell is large enough that the
// draw call, not the parse, is the cost. Sampling keeps the panel interactive
// on 10×10×10 models while leaving the shape of every site's cloud intact.
const MAX_CLOUD_POINTS = 120000;

// While a highlight is shown, the other bonds fade to this opacity and the
// highlighted ones are drawn this much thicker.
const DIMMED_BOND_OPACITY = 0.18;
const HIGHLIGHT_RADIUS = 1.7;

// One instanced tube mesh for a list of segments (six numbers each).
const tubeMesh = (segments, radius, material) => {
    const count = segments.length / 6;
    const mesh = new THREE.InstancedMesh(new THREE.CylinderGeometry(1, 1, 1, 8, 1, true), material, count);
    const up = new THREE.Vector3(0, 1, 0);
    const start = new THREE.Vector3();
    const end = new THREE.Vector3();
    const direction = new THREE.Vector3();
    const rotation = new THREE.Quaternion();
    const scale = new THREE.Vector3();
    const matrix = new THREE.Matrix4();
    for (let i = 0; i < count; i += 1) {
        start.fromArray(segments, i * 6);
        end.fromArray(segments, i * 6 + 3);
        direction.subVectors(end, start);
        const length = direction.length();
        rotation.setFromUnitVectors(up, direction.divideScalar(length));
        scale.set(radius, length, radius);
        matrix.compose(start.add(end).multiplyScalar(0.5), rotation, scale);
        mesh.setMatrixAt(i, matrix);
    }
    mesh.instanceMatrix.needsUpdate = true;
    mesh.computeBoundingSphere();
    return mesh;
};

const disposeGroup = (group) => {
    while (group.children.length) {
        const child = group.children.pop();
        child.geometry?.dispose();
        child.material?.map?.dispose();
        child.material?.dispose();
        if (child.isInstancedMesh) child.dispose();
        group.remove(child);
    }
};

const FoldedCellPanel = ({
    structure,
    sites,
    bondSets = null,
    // The bond network to draw (bondNetwork.js foldedCellNetwork); computed
    // here from `bondSets` when the caller does not pass it.
    network: networkProp,
    // The hovered angle bin's triplets (angleHighlight.js highlightForBin).
    highlight = null,
    elementColors = {},
    loading = false,
    title = 'Folded unit cell',
    // Elements to name in bold in the legend; the others read muted. Null:
    // every element is plain.
    legendEmphasis = null,
    // Extra legend items after the elements (e.g. the bond swatch, or a
    // "Compute to draw bonds" cue).
    legendExtras = null,
    className,
    // The AI Copilot's id for this card and 3D view (viewer3d.js).
    cardId = null
}) => {
    const [showCellAxes, setShowCellAxes] = useState(true);
    const mountRef = useRef(null);
    const sceneRef = useRef(null);

    // The bonds to draw, from the average sites of one unit cell.
    const ownNetwork = useMemo(
        () => (networkProp === undefined ? foldedCellNetwork(structure, sites, bondSets) : null),
        [networkProp, structure, sites, bondSets]
    );
    const network = networkProp === undefined ? ownNetwork : networkProp;

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
        const highlightGroup = new THREE.Group();
        const cellGroup = new THREE.Group();
        const cellAxesGroup = new THREE.Group();
        scene.add(cloudGroup, bondsGroup, highlightGroup, cellGroup, cellAxesGroup);
        // The bond tubes are shaded (the clouds are points, unlit): a soft
        // ambient plus a key light that rides with the camera, so the tubes
        // read as round from every side.
        const keyLight = new THREE.DirectionalLight(0xffffff, 1.6);
        keyLight.position.set(1, 1.5, 2);
        camera.add(keyLight);
        scene.add(camera, new THREE.AmbientLight(0xffffff, 1.1));

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
            scene, camera, renderer, controls, cloudGroup, bondsGroup, highlightGroup, cellGroup, cellAxesGroup, resize
        };

        return () => {
            cancelAnimationFrame(frame);
            observer.disconnect();
            controls.dispose();
            renderer.dispose();
            // dispose() frees GPU resources but not the WebGL context itself;
            // without this, scene rebuilds and remounts pile up contexts until the browser drops
            // the oldest ('Too many active WebGL contexts').
            renderer.forceContextLoss();
            mount.replaceChildren();
            sceneRef.current = null;
        };
    }, []);

    // Rebuild the cloud, the cell box and the bond network.
    useEffect(() => {
        const handle = sceneRef.current;
        if (!handle) return;
        const { cloudGroup, bondsGroup, cellGroup, cellAxesGroup, camera, controls } = handle;

        [cloudGroup, bondsGroup, cellGroup, cellAxesGroup].forEach(disposeGroup);

        const unit = unitVectorsOf(structure) || unitVectorsOf(sites);
        handle.unit = unit || null;
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
                opacity: cloudOpacity(element, network?.sets)
            })));
        });

        // --- Detected bonds (bondNetwork.js): one instanced tube mesh per
        // window, in its bond-role colour, opaque over the translucent cloud.
        network?.sets.forEach((set) => {
            if (!set.segments.length) return;
            bondsGroup.add(tubeMesh(
                set.segments,
                bondRadius(set.window, edgeLength),
                new THREE.MeshLambertMaterial({ color: set.color })
            ));
        });
        handle.edgeLength = edgeLength;

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
    }, [structure, sites, network, elementColors, showCellAxes]);

    // The hovered bin's triplets: their bonds as thicker tubes in their
    // bond-role colours (A–B, then B–C), the other bonds dimmed. Rebuilds only
    // this group; runs after the scene rebuild above, which it follows.
    useEffect(() => {
        const handle = sceneRef.current;
        if (!handle) return;
        const { bondsGroup, highlightGroup, edgeLength } = handle;
        disposeGroup(highlightGroup);
        const on = Boolean(highlight && network?.sets.length);
        bondsGroup.children.forEach((mesh) => {
            mesh.material.transparent = on;
            mesh.material.opacity = on ? DIMMED_BOND_OPACITY : 1;
            mesh.material.depthWrite = !on;
            mesh.material.needsUpdate = true;
        });
        if (!on || !edgeLength) return;
        const [first, second = first] = network.sets;
        [[highlight.ab, first], [highlight.bc, second]].forEach(([segments, set]) => {
            if (!segments?.length) return;
            highlightGroup.add(tubeMesh(
                segments,
                HIGHLIGHT_RADIUS * bondRadius(set.window, edgeLength),
                new THREE.MeshLambertMaterial({ color: set.color, emissive: set.color, emissiveIntensity: 0.25 })
            ));
        });
    }, [highlight, network, structure, sites, elementColors, showCellAxes]);

    const resetView = useCallback(() => {
        const handle = sceneRef.current;
        if (!handle) return;
        handle.controls.reset();
    }, []);

    // The AI Copilot's handle on this view (viewer3d.js): look along a crystal
    // direction, reset.
    useCopilotViewer(cardId, {
        orbit: {
            getHandle: () => sceneRef.current,
            getCell: () => sceneRef.current?.unit ?? null,
            onReset: () => resetView()
        },
        ready: () => Boolean(sceneRef.current?.unit)
    });

    // Same 1× / 3× options as the other 3D panels: 3× re-renders at a higher
    // pixel ratio, then restores the on-screen one.
    const saveFigure = useCallback(async (format) => {
        const handle = sceneRef.current;
        // A save that cannot happen says so (SaveMenu lists it), never passes as done.
        if (!handle) throw new Error('Nothing to save yet: the view has not been drawn');
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
            // A failed capture must say so (SaveMenu lists it), not do nothing.
            if (!blob) throw new Error('Could not capture the 3D view');
            downloadBlob(blob, `${sanitizeFilename(name)}.png`);
        } else {
            renderer.render(scene, camera);
            await saveCanvasAsPng(renderer.domElement, name);
        }
    }, []);

    return (
        <Card
            roundEnds
            className={className ? `pca-unitcell-panel ${className}` : 'pca-unitcell-panel'}
            data-copilot-card={cardId || undefined}
        >
            <CardHeader
                wrap
                title={title}
                help={(
                    <InfoBadge label="About the folded unit cell" align="start">
                        <p>
                            Every atom in the supercell folded back into one unit cell and
                            colored by element — the same view as the Atomic Density page.
                            The spread around each site is the thermal cloud itself, not a
                            fitted ellipsoid.
                        </p>
                        {bondSets?.length ? (
                            <p>
                                The tubes are the analysis' <em>detected bonds</em>, colored by
                                bond window: average site pairs whose distance falls inside the
                                window, periodic images included, so a bond may reach an image
                                just outside the box. Elements outside the bonds are faded. A
                                large cell draws at most {MAX_DRAWN_BONDS.toLocaleString()} bonds,
                                those nearest its centre.
                            </p>
                        ) : null}
                        {bondSets?.length ? (
                            <p>
                                Point at the angle plot (tap it on a touch screen) to pick out the
                                triplets of that bin: their two bonds, drawn between average sites,
                                so a drawn angle is the sites' average one. In a large cell only
                                those inside the drawn bonds are shown, and the legend says how
                                many. For a run opened in the browser (the Demo or a picked folder),
                                not yet a server folder.
                            </p>
                        ) : null}
                    </InfoBadge>
                )}
                actions={(
                    <>
                        {loading && <CardMeta>Loading…</CardMeta>}
                        <ToolButton
                            axes
                            active={showCellAxes}
                            onClick={() => setShowCellAxes((current) => !current)}
                            aria-label="Show the crystallographic axes (a, b, c)"
                        >
                            a b c
                        </ToolButton>
                        <ToolButton
                            onClick={resetView}
                            aria-label="Reset the camera to the default view"
                        >
                            Reset view
                        </ToolButton>
                        <SaveMenu onSave={saveFigure} options={SAVE_OPTIONS} label="Save" name="Bonds view" align="right" />
                    </>
                )}
            />
            <div className="ui-stage ui-stage--glow-soft ui-stage--orbit pca-structure pca-structure--overlaid">
                <div className="pca-structure__mount" ref={mountRef} />
                {(Object.keys(elementColors).length > 0 || legendExtras) && (
                    <div className="ui-legend ui-legend--overlay" aria-label="Atom colors by element">
                        {Object.entries(elementColors).map(([element, color]) => (
                            <span
                                key={element}
                                className={legendEmphasis && !legendEmphasis.includes(element)
                                    ? 'ui-legend__item is-muted'
                                    : 'ui-legend__item'}
                                title={elementColorNote(elementColors, element)}
                            >
                                <span className="ui-legend__swatch" style={{ background: color }} />
                                {element}
                            </span>
                        ))}
                        {legendExtras}
                        {highlight && (
                            <span
                                className="ui-legend__item is-highlight"
                                data-testid="angle-highlight"
                                title="The triplets of the hovered angle bin, on the average sites."
                            >
                                {highlightLabel(highlight)}
                            </span>
                        )}
                        {network?.limited && (
                            <span
                                className="ui-legend__item is-muted"
                                title={`The ${network.shown.toLocaleString()} bonds nearest the cell centre`
                                    + ` of ${network.total.toLocaleString()}: a large cell draws at most`
                                    + ` ${MAX_DRAWN_BONDS.toLocaleString()}.`}
                            >
                                {network.shown.toLocaleString()} of {network.total.toLocaleString()} bonds shown
                            </span>
                        )}
                    </div>
                )}
            </div>
        </Card>
    );
};

export default FoldedCellPanel;
