import React, { useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import API_BASE_URL from '../api';
import { COLORMAP_NAMES, getLut } from '../colormaps';
import ModelSummary from './ModelSummary';
import './StructurePage.css';

const colors = {
    Ga: '#4f8cff',
    Nb: '#f2b84b',
    Se: '#62c084',
    default: '#d5d9df'
};

const vectorLength = (vector) => Math.sqrt(vector.reduce((sum, value) => sum + value * value, 0));
const dot = (a, b) => a.reduce((sum, value, index) => sum + value * b[index], 0);
const add = (a, b) => a.map((value, index) => value + b[index]);
const subtract = (a, b) => a.map((value, index) => value - b[index]);
const scale = (vector, factor) => vector.map((value) => value * factor);

const vectorFromFraction = (fraction, basis) => basis.reduce(
    (acc, vector, index) => add(acc, scale(vector, fraction[index])),
    [0, 0, 0]
);

const toVector3 = (vector) => new THREE.Vector3(vector[0], vector[1], vector[2]);

const makePlane = (uVector, vVector) => {
    const uLength = vectorLength(uVector);
    const vLength = vectorLength(vVector);
    const cosine = dot(uVector, vVector) / Math.max(uLength * vLength, 1e-12);
    const clamped = Math.max(-1, Math.min(1, cosine));
    const sine = Math.sqrt(Math.max(0, 1 - clamped * clamped));
    const u = [uLength, 0];
    const v = [vLength * clamped, vLength * sine];
    const corners = [
        [0, 0],
        u,
        add(u, v),
        v
    ];
    const xs = corners.map((corner) => corner[0]);
    const ys = corners.map((corner) => corner[1]);
    const bounds = {
        minX: Math.min(...xs),
        maxX: Math.max(...xs),
        minY: Math.min(...ys),
        maxY: Math.max(...ys)
    };
    return {
        u,
        v,
        bounds,
        aspect: (bounds.maxX - bounds.minX) / Math.max(bounds.maxY - bounds.minY, 1e-9)
    };
};

const makePlaneMapper = (plane, width, height, padding = 18) => {
    const spanX = plane.bounds.maxX - plane.bounds.minX || 1;
    const spanY = plane.bounds.maxY - plane.bounds.minY || 1;
    const scaleFactor = Math.min(
        (width - padding * 2) / spanX,
        (height - padding * 2) / spanY
    );
    const offsetX = (width - spanX * scaleFactor) / 2;
    const offsetY = (height - spanY * scaleFactor) / 2;

    const map = (uFraction, vFraction) => {
        const x = plane.u[0] * uFraction + plane.v[0] * vFraction;
        const y = plane.u[1] * uFraction + plane.v[1] * vFraction;
        return {
            x: offsetX + (x - plane.bounds.minX) * scaleFactor,
            y: offsetY + (plane.bounds.maxY - y) * scaleFactor
        };
    };

    return { map };
};

const drawPolygon = (ctx, points) => {
    ctx.beginPath();
    points.forEach((point, index) => {
        if (index === 0) ctx.moveTo(point.x, point.y);
        else ctx.lineTo(point.x, point.y);
    });
    ctx.closePath();
};

const cellCorners = (basis, center, zStart = 0, zEnd = 1) => [
    [0, 0, zStart],
    [1, 0, zStart],
    [1, 1, zStart],
    [0, 1, zStart],
    [0, 0, zEnd],
    [1, 0, zEnd],
    [1, 1, zEnd],
    [0, 1, zEnd]
].map((fraction) => toVector3(subtract(vectorFromFraction(fraction, basis), center)));

const edgeIndexPairs = [
    [0, 1], [1, 2], [2, 3], [3, 0],
    [4, 5], [5, 6], [6, 7], [7, 4],
    [0, 4], [1, 5], [2, 6], [3, 7]
];

const makeCellEdgeGeometry = (corners) => {
    const points = edgeIndexPairs.flatMap(([start, end]) => [corners[start], corners[end]]);
    return new THREE.BufferGeometry().setFromPoints(points);
};

const makeSlabGeometry = (corners) => {
    const vertices = new Float32Array(corners.flatMap((corner) => [corner.x, corner.y, corner.z]));
    const indices = [
        0, 1, 2, 0, 2, 3,
        4, 6, 5, 4, 7, 6,
        0, 4, 5, 0, 5, 1,
        1, 5, 6, 1, 6, 2,
        2, 6, 7, 2, 7, 3,
        3, 7, 4, 3, 4, 0
    ];
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(vertices, 3));
    geometry.setIndex(indices);
    geometry.computeVertexNormals();
    return geometry;
};

const StructurePage = ({ directory, theme }) => {
    const [structure, setStructure] = useState(null);
    const [error, setError] = useState(null);
    const [loading, setLoading] = useState(false);
    const [selectedElement, setSelectedElement] = useState('all');
    const [zCenter, setZCenter] = useState(0.5);
    const [thickness, setThickness] = useState(0.08);
    const [bandwidth, setBandwidth] = useState(0.03);
    const [gridSize, setGridSize] = useState(120);
    const [colormap, setColormap] = useState('viridis');
    const [showContours, setShowContours] = useState(true);
    const [logScale, setLogScale] = useState(false);
    const [kde, setKde] = useState(null);
    const [kdeLoading, setKdeLoading] = useState(false);
    const [kdeError, setKdeError] = useState(null);
    const canvasRef = useRef(null);
    const slabCanvasRef = useRef(null);
    const mountRef = useRef(null);
    const themeVars = useMemo(() => {
        if (typeof window === 'undefined') {
            return {
                canvasBg: '#101419',
                text: '#ecf1f7',
                muted: '#9ba7b5',
                border: '#2d3540',
                contour: 'rgba(40, 44, 48, 0.85)'
            };
        }
        const styles = getComputedStyle(document.documentElement);
        return {
            canvasBg: styles.getPropertyValue('--canvas-bg').trim() || '#101419',
            text: styles.getPropertyValue('--text').trim() || '#ecf1f7',
            muted: styles.getPropertyValue('--muted').trim() || '#9ba7b5',
            border: styles.getPropertyValue('--border-strong').trim() || '#2d3540',
            contour: theme === 'light' ? 'rgba(21, 34, 50, 0.72)' : 'rgba(230, 236, 244, 0.76)'
        };
    }, [theme]);

    useEffect(() => {
        const fetchStructure = async () => {
            setLoading(true);
            setError(null);
            try {
                const response = await axios.get(`${API_BASE_URL}/api/structure`, {
                    params: { dir: directory || '.', maxPoints: 75000 }
                });
                setStructure(response.data);
            } catch (err) {
                setStructure(null);
                setError(err.response?.data?.error || 'No structure data available in this folder');
            } finally {
                setLoading(false);
            }
        };

        fetchStructure();
    }, [directory]);

    const points = useMemo(() => {
        const allPoints = structure?.points || [];
        if (selectedElement === 'all') {
            return allPoints;
        }
        return allPoints.filter((point) => point.element === selectedElement);
    }, [structure, selectedElement]);

    const ranges = useMemo(() => {
        if (!points.length) {
            return { x: [0, 1], y: [0, 1], z: [0, 1] };
        }
        const xs = points.map((point) => point.x);
        const ys = points.map((point) => point.y);
        const zs = points.map((point) => point.z);
        return {
            x: [Math.min(...xs), Math.max(...xs)],
            y: [Math.min(...ys), Math.max(...ys)],
            z: [Math.min(...zs), Math.max(...zs)]
        };
    }, [points]);

    const unitCell = useMemo(() => {
        if (!structure?.latticeVectors || !structure?.supercell) {
            const unitVectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
            return {
                lengths: [1, 1, 1],
                basis: unitVectors,
                center: [0.5, 0.5, 0.5],
                abPlane: makePlane(unitVectors[0], unitVectors[1]),
                acPlane: makePlane(unitVectors[0], unitVectors[2])
            };
        }
        const unitVectors = structure.latticeVectors.map((vector, index) => (
            vector.map((value) => value / Math.max(structure.supercell[index], 1e-12))
        ));
        const lengths = unitVectors.map(vectorLength);
        const maxLength = Math.max(...lengths, 1e-9);
        const basis = unitVectors.map((vector) => scale(vector, 1 / maxLength));
        const center = scale(basis.reduce((acc, vector) => add(acc, vector), [0, 0, 0]), 0.5);
        return {
            lengths,
            basis,
            center,
            abPlane: makePlane(unitVectors[0], unitVectors[1]),
            acPlane: makePlane(unitVectors[0], unitVectors[2])
        };
    }, [structure]);

    // Default the z-slice to the densest band so the slice is populated on load
    // (the geometric midpoint can fall in a gap between atomic layers).
    useEffect(() => {
        if (!points.length) return;
        const bins = 50;
        const counts = new Array(bins).fill(0);
        points.forEach((point) => {
            const bin = Math.max(0, Math.min(bins - 1, Math.floor(point.z * bins)));
            counts[bin] += 1;
        });
        let best = 0;
        counts.forEach((count, index) => {
            if (count > counts[best]) best = index;
        });
        setZCenter((best + 0.5) / bins);
    }, [points]);

    // Fetch a real server-side gaussian_kde slice (debounced) whenever a
    // parameter that changes the density field updates.
    useEffect(() => {
        if (!structure) return undefined;
        const controller = new AbortController();
        const handle = setTimeout(async () => {
            setKdeLoading(true);
            setKdeError(null);
            try {
                const response = await axios.get(`${API_BASE_URL}/api/kde/slice`, {
                    params: {
                        dir: directory || '.',
                        element: selectedElement,
                        z: zCenter,
                        dz: thickness,
                        bw: bandwidth,
                        grid: gridSize,
                        log: logScale,
                        levels: 8
                    },
                    signal: controller.signal
                });
                setKde(response.data);
            } catch (err) {
                if (!axios.isCancel(err) && err.code !== 'ERR_CANCELED') {
                    setKde(null);
                    setKdeError(err.response?.data?.error || 'KDE computation failed');
                }
            } finally {
                setKdeLoading(false);
            }
        }, 160);

        return () => {
            controller.abort();
            clearTimeout(handle);
        };
    }, [structure, directory, selectedElement, zCenter, thickness, bandwidth, gridSize, logScale]);

    // Render the KDE density grid + contour overlay. Colormap and contour
    // visibility are pure client-side re-renders (no refetch).
    useEffect(() => {
        const canvas = canvasRef.current;
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const size = canvas.getBoundingClientRect();
        const width = Math.max(320, Math.floor(size.width));
        const height = Math.max(260, Math.floor(size.height));
        const dpr = window.devicePixelRatio || 1;
        canvas.width = width * dpr;
        canvas.height = height * dpr;
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, width, height);
        ctx.fillStyle = themeVars.canvasBg;
        ctx.fillRect(0, 0, width, height);

        const mapper = makePlaneMapper(unitCell.abPlane, width, height, 18);
        const cellOutline = [
            mapper.map(0, 0),
            mapper.map(1, 0),
            mapper.map(1, 1),
            mapper.map(0, 1)
        ];
        const density = kde?.density;
        const grid = kde?.grid || 0;
        if (density && grid > 0 && kde.vmax > kde.vmin) {
            const lut = getLut(colormap);
            const offscreen = document.createElement('canvas');
            offscreen.width = grid;
            offscreen.height = grid;
            const offCtx = offscreen.getContext('2d');
            const imageData = offCtx.createImageData(grid, grid);
            const span = kde.vmax - kde.vmin || 1;
            for (let py = 0; py < grid; py += 1) {
                const densityRow = density[py];
                for (let px = 0; px < grid; px += 1) {
                    const normalized = (densityRow[px] - kde.vmin) / span;
                    const lutIndex = Math.max(0, Math.min(255, Math.round(normalized * 255))) * 3;
                    const offset = (py * grid + px) * 4;
                    imageData.data[offset] = lut[lutIndex];
                    imageData.data[offset + 1] = lut[lutIndex + 1];
                    imageData.data[offset + 2] = lut[lutIndex + 2];
                    imageData.data[offset + 3] = 255;
                }
            }
            offCtx.putImageData(imageData, 0, 0);
            ctx.imageSmoothingEnabled = true;
            const origin = mapper.map(0, 0);
            const uEdge = mapper.map(1, 0);
            const vEdge = mapper.map(0, 1);
            ctx.save();
            drawPolygon(ctx, cellOutline);
            ctx.clip();
            ctx.transform(
                uEdge.x - origin.x,
                uEdge.y - origin.y,
                vEdge.x - origin.x,
                vEdge.y - origin.y,
                origin.x,
                origin.y
            );
            ctx.drawImage(offscreen, 0, 0, grid, grid, 0, 0, 1, 1);
            ctx.restore();

            if (showContours && kde.contours?.length) {
                const [xMin, xMax, yMin, yMax] = kde.extent;
                const spanX = xMax - xMin || 1;
                const spanY = yMax - yMin || 1;
                ctx.lineWidth = 1;
                ctx.strokeStyle = themeVars.contour;
                kde.contours.forEach((contour) => {
                    contour.lines.forEach((line) => {
                        ctx.beginPath();
                        line.forEach(([dataX, dataY], index) => {
                            const point = mapper.map((dataX - xMin) / spanX, (dataY - yMin) / spanY);
                            if (index === 0) ctx.moveTo(point.x, point.y);
                            else ctx.lineTo(point.x, point.y);
                        });
                        ctx.stroke();
                    });
                });
            }
        } else {
            ctx.fillStyle = themeVars.muted;
            ctx.font = '500 13px Inter, system-ui';
            ctx.fillText(kdeLoading ? 'Computing KDE...' : 'No atoms in this slab', 14, 28);
        }

        ctx.strokeStyle = themeVars.border;
        ctx.lineWidth = 1;
        drawPolygon(ctx, cellOutline);
        ctx.stroke();
        // Outlined text stays legible over any colormap.
        const drawOverlayText = (text, x, y) => {
            ctx.lineWidth = 3;
            ctx.lineJoin = 'round';
            ctx.strokeStyle = 'rgba(13, 18, 28, 0.62)';
            ctx.strokeText(text, x, y);
            ctx.fillStyle = '#fff';
            ctx.fillText(text, x, y);
        };
        ctx.font = '500 12px Inter, system-ui';
        if (kde) {
            drawOverlayText(`${kde.slabCount} atoms in slab (fit ${kde.fitCount})`, 12, 22);
            drawOverlayText(`c=${kde.z.toFixed(3)}  dc=${kde.dz.toFixed(3)}  bw=${kde.bw}`, 12, 40);
            if (kde.log) drawOverlayText('log10 density', 12, 58);
        }
    }, [kde, colormap, showContours, kdeLoading, themeVars, unitCell]);

    useEffect(() => {
        const canvas = slabCanvasRef.current;
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const size = canvas.getBoundingClientRect();
        const width = Math.max(220, Math.floor(size.width));
        const height = Math.max(260, Math.floor(size.height));
        const dpr = window.devicePixelRatio || 1;
        canvas.width = width * dpr;
        canvas.height = height * dpr;
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, width, height);
        ctx.fillStyle = themeVars.canvasBg;
        ctx.fillRect(0, 0, width, height);

        const mapper = makePlaneMapper(unitCell.acPlane, width, height, 18);
        const cellOutline = [
            mapper.map(0, 0),
            mapper.map(1, 0),
            mapper.map(1, 1),
            mapper.map(0, 1)
        ];
        const zStart = Math.max(0, zCenter - thickness / 2);
        const zEnd = Math.min(1, zCenter + thickness / 2);
        const slabOutline = [
            mapper.map(0, zStart),
            mapper.map(1, zStart),
            mapper.map(1, zEnd),
            mapper.map(0, zEnd)
        ];

        ctx.strokeStyle = themeVars.border;
        ctx.lineWidth = 1;
        drawPolygon(ctx, cellOutline);
        ctx.stroke();

        ctx.fillStyle = 'rgba(79, 140, 255, 0.18)';
        drawPolygon(ctx, slabOutline);
        ctx.fill();
        ctx.strokeStyle = '#74a7ff';
        drawPolygon(ctx, slabOutline);
        ctx.stroke();

        const sampleLimit = Math.min(points.length, 9000);
        const stride = Math.max(1, Math.floor(points.length / sampleLimit));
        for (let index = 0; index < points.length; index += stride) {
            const point = points[index];
            const projected = mapper.map(point.x, point.z);
            const inSlab = Math.abs(point.z - zCenter) <= thickness / 2;
            ctx.fillStyle = inSlab ? (colors[point.element] || colors.default) : 'rgba(166, 176, 188, 0.22)';
            ctx.fillRect(projected.x, projected.y, inSlab ? 2 : 1, inSlab ? 2 : 1);
        }

        ctx.fillStyle = themeVars.text;
        ctx.font = '500 12px Inter, system-ui';
        const xLabel = mapper.map(1, 0);
        const zLabel = mapper.map(0, 1);
        const slabLabel = mapper.map(0, zStart);
        ctx.fillText('a', Math.min(width - 16, xLabel.x + 4), Math.min(height - 8, xLabel.y + 14));
        ctx.fillText('c', Math.max(8, zLabel.x - 12), Math.max(16, zLabel.y - 6));
        ctx.fillText(`c=${zCenter.toFixed(3)}`, 10, Math.max(30, slabLabel.y - 6));
        ctx.fillText(`dc=${thickness.toFixed(3)}`, 10, Math.min(height - 16, slabLabel.y + 18));
    }, [points, zCenter, thickness, unitCell, themeVars]);

    useEffect(() => {
        const mount = mountRef.current;
        if (!mount || points.length === 0) return undefined;

        const width = mount.clientWidth;
        const height = mount.clientHeight;
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(themeVars.canvasBg);
        const camera = new THREE.PerspectiveCamera(45, width / height, 0.01, 20);

        const renderer = new THREE.WebGLRenderer({ antialias: true });
        renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        renderer.setSize(width, height);
        mount.replaceChildren(renderer.domElement);

        const controls = new OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.dampingFactor = 0.08;
        controls.enablePan = true;

        const group = new THREE.Group();
        scene.add(group);

        const byElement = points.reduce((acc, point) => {
            acc[point.element] = acc[point.element] || [];
            acc[point.element].push(point);
            return acc;
        }, {});

        Object.entries(byElement).forEach(([element, elementPoints]) => {
            const geometry = new THREE.BufferGeometry();
            const positions = new Float32Array(elementPoints.length * 3);
            elementPoints.forEach((point, index) => {
                const position = subtract(
                    vectorFromFraction([point.x, point.y, point.z], unitCell.basis),
                    unitCell.center
                );
                positions[index * 3] = position[0];
                positions[index * 3 + 1] = position[1];
                positions[index * 3 + 2] = position[2];
            });
            geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
            const material = new THREE.PointsMaterial({
                color: colors[element] || colors.default,
                size: 0.018,
                sizeAttenuation: true
            });
            group.add(new THREE.Points(geometry, material));
        });

        const corners = cellCorners(unitCell.basis, unitCell.center);
        const cellEdgeGeometry = makeCellEdgeGeometry(corners);
        const cellEdgeMaterial = new THREE.LineBasicMaterial({ color: '#737c86' });
        const cellEdges = new THREE.LineSegments(cellEdgeGeometry, cellEdgeMaterial);
        scene.add(cellEdges);

        const zStart = Math.max(0, zCenter - thickness / 2);
        const zEnd = Math.min(1, zCenter + thickness / 2);
        const slabCorners = cellCorners(unitCell.basis, unitCell.center, zStart, zEnd);
        const slabGeometry = makeSlabGeometry(slabCorners);
        const slabMaterial = new THREE.MeshBasicMaterial({
            color: '#4f8cff',
            transparent: true,
            opacity: 0.12,
            side: THREE.DoubleSide,
            depthWrite: false
        });
        const slabMesh = new THREE.Mesh(slabGeometry, slabMaterial);
        scene.add(slabMesh);

        const slabEdgeGeometry = makeCellEdgeGeometry(slabCorners);
        const slabEdgeMaterial = new THREE.LineBasicMaterial({
            color: '#8c96a3',
            transparent: true,
            opacity: 0.95
        });
        const slabEdges = new THREE.LineSegments(slabEdgeGeometry, slabEdgeMaterial);
        scene.add(slabEdges);

        const bounds = new THREE.Box3().setFromPoints(corners);
        const sphere = new THREE.Sphere();
        bounds.getBoundingSphere(sphere);
        const radius = Math.max(sphere.radius, 0.5);
        camera.near = radius / 100;
        camera.far = radius * 20;
        camera.position.copy(sphere.center).add(new THREE.Vector3(radius * 1.7, radius * 1.45, radius * 1.55));
        camera.lookAt(sphere.center);
        camera.updateProjectionMatrix();
        controls.target.copy(sphere.center);
        controls.minDistance = radius * 0.35;
        controls.maxDistance = radius * 8;

        let frameId;
        const animate = () => {
            controls.update();
            renderer.render(scene, camera);
            frameId = requestAnimationFrame(animate);
        };
        animate();

        return () => {
            cancelAnimationFrame(frameId);
            controls.dispose();
            renderer.dispose();
            slabGeometry.dispose();
            slabMaterial.dispose();
            cellEdgeGeometry.dispose();
            cellEdgeMaterial.dispose();
            slabEdgeGeometry.dispose();
            slabEdgeMaterial.dispose();
            group.traverse((object) => {
                object.geometry?.dispose?.();
                object.material?.dispose?.();
            });
        };
    }, [points, unitCell, zCenter, thickness, themeVars]);

    return (
        <section className="structure-page">
            <div className="structure-header">
                <div>
                    <h2>KDE And 3D Model</h2>
                    <p>{directory}</p>
                </div>
                {(loading || kdeLoading) && <span className="status-pill">{loading ? 'Loading' : 'KDE'}</span>}
            </div>

            {error && <div className="structure-error">{error}</div>}

            {structure && (
                <>
                    <ModelSummary structure={structure} />

                    <div className="structure-controls">
                        <label className="control">
                            <span className="control-name">Element</span>
                            <select value={selectedElement} onChange={(event) => setSelectedElement(event.target.value)}>
                                <option value="all">All</option>
                                {structure.elements.map((element) => (
                                    <option key={element} value={element}>{element}</option>
                                ))}
                            </select>
                        </label>
                        <label className="control">
                            <span className="control-name">Slice z</span>
                            <input
                                type="range"
                                min={ranges.z[0]}
                                max={ranges.z[1]}
                                step="0.001"
                                value={zCenter}
                                onChange={(event) => setZCenter(Number(event.target.value))}
                            />
                            <span className="control-value">{zCenter.toFixed(2)}</span>
                        </label>
                        <label className="control">
                            <span className="control-name">Thickness</span>
                            <input type="range" min="0.01" max="0.5" step="0.01" value={thickness} onChange={(event) => setThickness(Number(event.target.value))} />
                            <span className="control-value">{thickness.toFixed(2)}</span>
                        </label>
                        <label className="control">
                            <span className="control-name">Bandwidth</span>
                            <input type="range" min="0.005" max="0.15" step="0.005" value={bandwidth} onChange={(event) => setBandwidth(Number(event.target.value))} />
                            <span className="control-value">{bandwidth.toFixed(3)}</span>
                        </label>
                        <label className="control">
                            <span className="control-name">Colormap</span>
                            <select value={colormap} onChange={(event) => setColormap(event.target.value)}>
                                {COLORMAP_NAMES.map((name) => (
                                    <option key={name} value={name}>{name}</option>
                                ))}
                            </select>
                        </label>
                        <label className="control">
                            <span className="control-name">Grid</span>
                            <select value={gridSize} onChange={(event) => setGridSize(Number(event.target.value))}>
                                <option value={80}>80</option>
                                <option value={120}>120</option>
                                <option value={160}>160</option>
                                <option value={220}>220</option>
                            </select>
                        </label>
                        <label className="control switch">
                            <span className="control-name">Contours</span>
                            <input type="checkbox" checked={showContours} onChange={(event) => setShowContours(event.target.checked)} />
                            <i className="switch-track" aria-hidden="true" />
                        </label>
                        <label className="control switch">
                            <span className="control-name">Log scale</span>
                            <input type="checkbox" checked={logScale} onChange={(event) => setLogScale(event.target.checked)} />
                            <i className="switch-track" aria-hidden="true" />
                        </label>
                    </div>

                    {kdeError && <div className="structure-error">{kdeError}</div>}

                    <div className="analysis-layout">
                        <div
                            className="kde-panel"
                            style={{ '--panel-aspect': unitCell.abPlane.aspect }}
                        >
                            <h3>KDE Slice</h3>
                            <canvas ref={canvasRef} className="kde-canvas" />
                        </div>
                        <div
                            className="slab-panel"
                            style={{ '--panel-aspect': unitCell.acPlane.aspect }}
                        >
                            <div className="slab-panel-header">
                                <span>Slab In Cell</span>
                                <strong>{(Math.max(0, zCenter - thickness / 2)).toFixed(2)} - {(Math.min(1, zCenter + thickness / 2)).toFixed(2)}</strong>
                            </div>
                            <canvas ref={slabCanvasRef} />
                        </div>
                        <div
                            className="model-panel"
                            style={{ '--panel-aspect': Math.max(unitCell.abPlane.aspect, 1) }}
                        >
                            <h3>3D Model</h3>
                            <div ref={mountRef} className="three-mount" />
                        </div>
                    </div>
                </>
            )}
        </section>
    );
};

export default StructurePage;
