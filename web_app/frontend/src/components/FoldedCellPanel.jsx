// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Folded unit cell with the analysis' detected bonds drawn over it. The atom
// cloud is the same one the Atomic Density page shows — every atom in the
// supercell folded back into a single cell, so the spread around each site is
// the thermal cloud itself rather than a fitted ellipsoid. Bonds come from the
// triplet windows (bondNetwork.js) and are drawn as shaded tubes in their
// bond-role colours (bondTubes.js: lit, round-ended, stopping at the edge of
// each site's cloud), one instanced mesh (one draw call) per window; the
// clouds of elements outside the bonds fade, so the selected triplet stands
// out. A large cell draws at most MAX_DRAWN_BONDS, those nearest the centre,
// and the legend says so. The element key sits inside the canvas, as a pill
// at the bottom left.
//
// `highlight` (angleHighlight.js highlightForBins: the triplets of the angle
// range, or of the bin under the pointer) draws those triplets' bonds as
// tubes of their own, each as thick and as strong in colour as its share,
// marks each triplet's angle with a shaded wedge at its central atom
// (angleArcs.js), and draws the other bonds thin and grey, paler with
// distance; the legend names the range. It has its own group, so a hover
// never rebuilds the cloud.

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
import { highlightLabel, highlightStrength, MAX_HIGHLIGHT_TRIPLETS, MIN_HIGHLIGHT_SHARE } from '../angleHighlight';
import { bondMaterial, renderQuality, setFadeRange, tubeMesh } from './bondTubes';
import { wedgeMesh } from './angleArcs';
import useCopilotViewer from '../useCopilotViewer';
import './UnitCellPanel.css';
import { isDisplayed } from './pageVisibility';

const SAVE_OPTIONS = [
    { id: 'png', label: 'Standard PNG', hint: '1×' },
    { id: 'png3x', label: 'High quality PNG', hint: '3×' }
];

// The cloud is one point per atom; a full supercell is large enough that the
// draw call, not the parse, is the cost. Sampling keeps the panel interactive
// on 10×10×10 models while leaving the shape of every site's cloud intact.
const MAX_CLOUD_POINTS = 120000;

// The tubes are drawn this much thinner than bondNetwork's bondRadius.
const BOND_RADIUS_SCALE = 0.8;
// While a highlight is shown, the other bonds are drawn REST_RADIUS of their
// radius, in REST_COLOR, paling to REST_FADE_COLOR from REST_FADE_NEAR cell
// edges in front of the view's centre to REST_FADE_FAR behind it.
const REST_RADIUS = 0.42;
const REST_COLOR = '#aab3bf';
const REST_FADE_COLOR = '#e3e7ec';
const REST_FADE_NEAR = 0.5;
const REST_FADE_FAR = 0.8;
// A highlighted bond of strength s (highlightStrength: the square root of
// its share) is drawn HIGHLIGHT_RADIUS + HIGHLIGHT_RADIUS_GAIN·s of its
// radius, opaque, its colour from WEAK_MIX of the way towards WEAK_TINT
// (s = 0) to its bond-role colour (s = 1). Thinner than a plain bond: the
// wedges at the apex carry the detail, and thick bonds hid them.
const HIGHLIGHT_RADIUS = 0.6;
const HIGHLIGHT_RADIUS_GAIN = 0.3;
const HIGHLIGHT_SHININESS = 24;
const WEAK_TINT = '#dce4ec';
const WEAK_MIX = 0.78;

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
        const quality = renderQuality();
        const renderer = new THREE.WebGLRenderer({
            antialias: quality.antialias,
            alpha: true,
            preserveDrawingBuffer: true
        });
        renderer.setPixelRatio(Math.min(window.devicePixelRatio, quality.maxPixelRatio));
        mount.replaceChildren(renderer.domElement);

        const controls = new OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.dampingFactor = 0.08;

        const cloudGroup = new THREE.Group();
        const bondsGroup = new THREE.Group();
        // The bonds as the highlight's backdrop; shown only with a highlight.
        const restGroup = new THREE.Group();
        restGroup.visible = false;
        const highlightGroup = new THREE.Group();
        const cellGroup = new THREE.Group();
        const cellAxesGroup = new THREE.Group();
        scene.add(cloudGroup, bondsGroup, restGroup, highlightGroup, cellGroup, cellAxesGroup);
        // The bond tubes are shaded (the clouds are points, unlit): a soft
        // sky/ground fill plus a key light that rides with the camera, so the
        // tubes read as round from every side.
        const keyLight = new THREE.DirectionalLight(0xffffff, 1.7);
        keyLight.position.set(1, 1.5, 2);
        camera.add(keyLight);
        scene.add(camera, new THREE.HemisphereLight(0xffffff, 0x9aa3b0, 1.25));

        let frame = 0;
        const renderLoop = () => {
            // A page the app keeps mounted but hidden draws nothing.
            if (isDisplayed(renderer.domElement)) {
                controls.update();
                renderer.render(scene, camera);
            }
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
            scene, camera, renderer, controls, cloudGroup, bondsGroup, restGroup, highlightGroup, cellGroup, cellAxesGroup, resize
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
        const { cloudGroup, bondsGroup, restGroup, cellGroup, cellAxesGroup, camera, controls } = handle;

        [cloudGroup, bondsGroup, restGroup, cellGroup, cellAxesGroup].forEach(disposeGroup);

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
                BOND_RADIUS_SCALE * bondRadius(set.window, edgeLength),
                bondMaterial({ color: set.color })
            ));
            // The same bonds as the highlight draws them when out of range.
            const material = bondMaterial({ color: REST_COLOR, specular: 0x000000, fade: REST_FADE_COLOR });
            const rest = tubeMesh(
                set.segments,
                REST_RADIUS * BOND_RADIUS_SCALE * bondRadius(set.window, edgeLength),
                material
            );
            rest.onBeforeRender = () => {
                const distance = camera.position.distanceTo(controls.target);
                setFadeRange(material, distance - REST_FADE_NEAR * edgeLength, distance + REST_FADE_FAR * edgeLength);
            };
            restGroup.add(rest);
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

    // The range's (or hovered bin's) triplets: their bonds as opaque tubes in
    // their bond-role colours (A–B, then B–C), each as thick and as strong in
    // colour as its share, and their angles as wedges; the other bonds thin
    // and grey. Rebuilds only this group; runs after the scene rebuild above,
    // which it follows.
    useEffect(() => {
        const handle = sceneRef.current;
        if (!handle) return;
        const { bondsGroup, restGroup, highlightGroup, edgeLength } = handle;
        disposeGroup(highlightGroup);
        const on = Boolean(highlight && network?.sets.length);
        bondsGroup.visible = !on;
        restGroup.visible = on;
        if (!on || !edgeLength) return;
        const [first, second = first] = network.sets;
        [[highlight.ab, highlight.abShare, first], [highlight.bc, highlight.bcShare, second]].forEach(([segments, shares, set]) => {
            if (!segments?.length) return;
            const strengths = shares ? Array.from(shares, highlightStrength) : Array(segments.length / 6).fill(1);
            const strong = new THREE.Color(set.color);
            const weak = strong.clone().lerp(new THREE.Color(WEAK_TINT), WEAK_MIX);
            highlightGroup.add(tubeMesh(
                segments,
                BOND_RADIUS_SCALE * bondRadius(set.window, edgeLength),
                bondMaterial({ color: '#ffffff', shininess: HIGHLIGHT_SHININESS }),
                {
                    scales: strengths.map((s) => HIGHLIGHT_RADIUS + HIGHLIGHT_RADIUS_GAIN * s),
                    colors: strengths.map((s) => weak.clone().lerp(strong, s))
                }
            ));
        });
        // The angle of each drawn triplet, as a wedge at its apex.
        if (highlight.arcs?.length) {
            const wedges = wedgeMesh(highlight.arcs, { strengths: Array.from(highlight.arcShare, highlightStrength) });
            if (wedges) {
                wedges.renderOrder = 2;
                highlightGroup.add(wedges);
            }
        }
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
                                Drag across the angle plot to pick a range of angles: its triplets
                                are drawn here, their two bonds brought forward, between average sites (so
                                a drawn angle is the sites' average one). Drag the band or its
                                edges to move it; double-click or Clear removes it. Without a
                                range, pointing at a bin (tapping it on a touch screen) previews
                                that bin. The range moves in whole bins, so its count is the
                                histogram's. In a large cell only the triplets inside the drawn
                                bonds are shown, the largest shares first, and the legend says how
                                many. For a run opened
                                in the browser (the Demo or a picked folder), not yet a server
                                folder.
                            </p>
                        ) : null}
                        {bondSets?.length ? (
                            <p>
                                <b>Strength</b> — the share of a triplet's angles that fall in the
                                range: thick and full-coloured when nearly all of them do, thin and pale when only the tail
                                of its spread reaches it, not drawn under {Math.round(MIN_HIGHLIGHT_SHARE * 100)}%. A bond shared by
                                several triplets shows the strongest. The other bonds turn thin and grey.
                            </p>
                        ) : null}
                        {bondSets?.length ? (
                            <p>
                                <b>Wedges</b> — a shaded wedge at the central atom marks each angle
                                in the range, as strong as its triplet.
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
                                title={'The triplets of the angle range, on the average sites; each bond as strong as'
                                    + ' the share of its angles in the range. A large cell draws at most'
                                    + ` ${MAX_HIGHLIGHT_TRIPLETS.toLocaleString()} triplets, the largest shares first.`}
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
