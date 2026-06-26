import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import API_BASE_URL from '../api';
import { isStaticMode } from '../browserData';
import { COLORMAP_NAMES, getLut } from '../colormaps';
import { buildElementColors, DEFAULT_ELEMENT_COLOR } from '../atomColors';
import { downloadBlob, sanitizeFilename, saveCanvasAsPng } from '../figureExport';
import ModelSummary from './ModelSummary';
import SaveMenu from './SaveMenu';
import './StructurePage.css';

const vectorLength = (vector) => Math.sqrt(vector.reduce((sum, value) => sum + value * value, 0));
const dot = (a, b) => a.reduce((sum, value, index) => sum + value * b[index], 0);
const add = (a, b) => a.map((value, index) => value + b[index]);
const subtract = (a, b) => a.map((value, index) => value - b[index]);
const scale = (vector, factor) => vector.map((value) => value * factor);
const cross = (a, b) => [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
];
const normalize = (vector, fallback = [0, 0, 1]) => {
    const length = vectorLength(vector);
    if (length <= 1e-9) return fallback;
    return vector.map((value) => value / length);
};

const CUBE_CORNERS = [
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
    [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]
];
const CUBE_EDGES = [
    [0, 1], [0, 2], [1, 3], [2, 3],
    [4, 5], [4, 6], [5, 7], [6, 7],
    [0, 4], [1, 5], [2, 6], [3, 7]
];
const STRUCTURE_MAX_POINTS = 1000000;
const SLAB_CANVAS_MAX_POINTS = 1000000;

const SLICE_PRESETS = {
    a: { label: 'a', normal: [1, 0, 0], u: [0, 1, 0], v: [0, 0, 1], uLabel: 'b', vLabel: 'c' },
    b: { label: 'b', normal: [0, 1, 0], u: [1, 0, 0], v: [0, 0, 1], uLabel: 'a', vLabel: 'c' },
    c: { label: 'c', normal: [0, 0, 1], u: [1, 0, 0], v: [0, 1, 0], uLabel: 'a', vLabel: 'b' }
};
const NORMAL_OPTIONS = [
    { value: 'a', label: 'a' },
    { value: 'b', label: 'b' },
    { value: 'c', label: 'c' },
    { value: 'custom', label: 'Custom' }
];
const CUSTOM_DIRECTION_LABELS = ['a', 'b', 'c'];

// KDE/3D panels are raster, so they export as PNG at native or 3x resolution.
const PANEL_SAVE_OPTIONS = [
    { id: 'png', label: 'PNG image', hint: '1×' },
    { id: 'png3x', label: 'High-res PNG', hint: '3×' },
];

const projectionRange = (normal) => {
    const values = CUBE_CORNERS.map((corner) => dot(corner, normal));
    return [Math.min(...values), Math.max(...values)];
};

const makeFreePlaneBasis = (normal) => {
    const reference = Math.abs(normal[0]) < 0.85 ? [1, 0, 0] : [0, 1, 0];
    const u = normalize(subtract(reference, scale(normal, dot(reference, normal))), [0, 1, 0]);
    const v = normalize(cross(normal, u), [0, 0, 1]);
    return { u, v };
};

const makeSliceConfig = (sliceDirection, customDirection) => {
    if (sliceDirection !== 'custom') {
        const preset = SLICE_PRESETS[sliceDirection] || SLICE_PRESETS.c;
        const normal = normalize(preset.normal);
        return {
            key: sliceDirection,
            label: preset.label,
            normal,
            u: normalize(preset.u),
            v: normalize(preset.v),
            uLabel: preset.uLabel,
            vLabel: preset.vLabel,
            range: projectionRange(normal)
        };
    }
    const normal = normalize(customDirection, [0, 0, 1]);
    const { u, v } = makeFreePlaneBasis(normal);
    return {
        key: 'custom',
        label: `[${customDirection.map((value) => Number(value).toLocaleString(undefined, { maximumFractionDigits: 2 })).join(' ')}]`,
        normal,
        u,
        v,
        uLabel: 'u',
        vLabel: 'v',
        range: projectionRange(normal)
    };
};

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

const makeProjectedPlane = (uVector, vVector, uvPoints) => {
    const base = makePlane(uVector, vVector);
    const projected = uvPoints.map(([uValue, vValue]) => (
        add(scale(base.u, uValue), scale(base.v, vValue))
    ));
    const xs = projected.map((point) => point[0]);
    const ys = projected.map((point) => point[1]);
    return {
        u: base.u,
        v: base.v,
        bounds: {
            minX: Math.min(...xs),
            maxX: Math.max(...xs),
            minY: Math.min(...ys),
            maxY: Math.max(...ys)
        },
        aspect: (Math.max(...xs) - Math.min(...xs)) / Math.max(Math.max(...ys) - Math.min(...ys), 1e-9)
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

    // Recover the (uFraction, vFraction) plane coordinates from a screen point.
    const invert = (screenX, screenY) => {
        const x = (screenX - offsetX) / scaleFactor + plane.bounds.minX;
        const y = plane.bounds.maxY - (screenY - offsetY) / scaleFactor;
        const det = plane.u[0] * plane.v[1] - plane.v[0] * plane.u[1];
        if (Math.abs(det) < 1e-12) return { uFraction: 0, vFraction: 0 };
        return {
            uFraction: (x * plane.v[1] - plane.v[0] * y) / det,
            vFraction: (plane.u[0] * y - x * plane.u[1]) / det
        };
    };

    return { map, invert };
};

const drawPolygon = (ctx, points) => {
    ctx.beginPath();
    points.forEach((point, index) => {
        if (index === 0) ctx.moveTo(point.x, point.y);
        else ctx.lineTo(point.x, point.y);
    });
    ctx.closePath();
};

const cellCorners = (basis, center) => CUBE_CORNERS.map((fraction) => (
    toVector3(subtract(vectorFromFraction(fraction, basis), center))
));

const makeCellEdgeGeometry = (corners) => {
    const points = CUBE_EDGES.flatMap(([start, end]) => [corners[start], corners[end]]);
    return new THREE.BufferGeometry().setFromPoints(points);
};

const planeSectionVertices = (normal, offset) => {
    const vertices = [];
    CUBE_EDGES.forEach(([start, end]) => {
        const p0 = CUBE_CORNERS[start];
        const p1 = CUBE_CORNERS[end];
        const d0 = dot(p0, normal) - offset;
        const d1 = dot(p1, normal) - offset;
        if (Math.abs(d0) <= 1e-9) vertices.push(p0);
        if (Math.abs(d1) <= 1e-9) vertices.push(p1);
        if (d0 * d1 < 0) {
            const t = d0 / (d0 - d1);
            vertices.push(add(p0, scale(subtract(p1, p0), t)));
        }
    });

    const unique = [];
    vertices.forEach((vertex) => {
        if (!unique.some((existing) => vectorLength(subtract(vertex, existing)) <= 1e-8)) {
            unique.push(vertex);
        }
    });
    if (unique.length < 3) return [];

    const { u, v } = makeFreePlaneBasis(normal);
    const center = scale(unique.reduce((acc, vertex) => add(acc, vertex), [0, 0, 0]), 1 / unique.length);
    return unique.sort((a, b) => {
        const aDelta = subtract(a, center);
        const bDelta = subtract(b, center);
        return Math.atan2(dot(aDelta, v), dot(aDelta, u)) - Math.atan2(dot(bDelta, v), dot(bDelta, u));
    });
};

const makeSlabGeometry = (sections) => {
    const corners = sections.flat();
    if (!corners.length) return new THREE.BufferGeometry();
    const vertices = new Float32Array(corners.flatMap((corner) => [corner.x, corner.y, corner.z]));
    const indices = [];
    let offset = 0;
    sections.forEach((section) => {
        if (section.length >= 3) {
            for (let index = 1; index < section.length - 1; index += 1) {
                indices.push(offset, offset + index, offset + index + 1);
            }
        }
        offset += section.length;
    });
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(vertices, 3));
    geometry.setIndex(indices);
    geometry.computeVertexNormals();
    return geometry;
};

const makeSectionEdgeGeometry = (sections) => {
    const points = [];
    sections.forEach((section) => {
        section.forEach((point, index) => {
            points.push(point, section[(index + 1) % section.length]);
        });
    });
    if (sections.length === 2 && sections[0].length === sections[1].length) {
        sections[0].forEach((point, index) => {
            points.push(point, sections[1][index]);
        });
    }
    return new THREE.BufferGeometry().setFromPoints(points);
};

const StructurePage = ({ directory, localRun, theme }) => {
    const [structure, setStructure] = useState(null);
    const [error, setError] = useState(null);
    const [loading, setLoading] = useState(false);
    const [selectedElement, setSelectedElement] = useState('all');
    const [sliceDirection, setSliceDirection] = useState('c');
    const [normalMenuOpen, setNormalMenuOpen] = useState(false);
    const [customDirection, setCustomDirection] = useState([1, 1, 0]);
    const [zCenter, setZCenter] = useState(0.5);
    const [thickness, setThickness] = useState(0.08);
    const [bandwidth, setBandwidth] = useState(0.03);
    const [gridSize, setGridSize] = useState(120);
    const [colormap, setColormap] = useState('viridis');
    const [showContours, setShowContours] = useState(false);
    const [logScale, setLogScale] = useState(false);
    const [kde, setKde] = useState(null);
    const [kdeLoading, setKdeLoading] = useState(false);
    const [kdeError, setKdeError] = useState(null);
    const canvasRef = useRef(null);
    const slabCanvasRef = useRef(null);
    const mountRef = useRef(null);
    // three.js render context (renderer/scene/camera + size), used to re-render
    // the model at a higher resolution for export.
    const modelExportRef = useRef(null);
    // Geometry of the last "Slab In Cell" render, used to drag the slice band.
    const slabGeometryRef = useRef(null);
    const slabDragRef = useRef(null);
    const zCenterRef = useRef(zCenter);
    const normalMenuRef = useRef(null);
    const cameraStateRef = useRef(null);
    const localStructureWorkerRef = useRef(null);
    const localStructureRequestRef = useRef(0);
    const localKdeWorkerRef = useRef(null);
    const localKdeRequestRef = useRef(0);
    const isLocalStructure = Boolean(localRun);
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
        if (!isLocalStructure) return undefined;
        const worker = new Worker(new URL('../workers/localStructureWorker.js', import.meta.url), {
            type: 'module'
        });
        localStructureWorkerRef.current = worker;

        worker.onmessage = (event) => {
            if (event.data.id !== localStructureRequestRef.current) return;
            setLoading(false);
            if (event.data.error) {
                setStructure(null);
                setError(event.data.error);
                return;
            }
            setStructure(event.data.result);
            setError(null);
        };

        worker.onerror = () => {
            setLoading(false);
            setStructure(null);
            setError('Browser structure parser failed');
        };

        return () => {
            worker.terminate();
            localStructureWorkerRef.current = null;
        };
    }, [isLocalStructure]);

    useEffect(() => {
        if (!isLocalStructure) return undefined;
        const worker = new Worker(new URL('../workers/localKdeWorker.js', import.meta.url), {
            type: 'module'
        });
        localKdeWorkerRef.current = worker;

        worker.onmessage = (event) => {
            if (event.data.id !== localKdeRequestRef.current) return;
            setKdeLoading(false);
            if (event.data.error) {
                setKde(null);
                setKdeError(event.data.error);
                return;
            }
            setKde(event.data.result);
            setKdeError(null);
        };

        worker.onerror = () => {
            setKdeLoading(false);
            setKde(null);
            setKdeError('Browser KDE worker failed');
        };

        return () => {
            worker.terminate();
            localKdeWorkerRef.current = null;
        };
    }, [isLocalStructure]);

    useEffect(() => {
        if (localRun) {
            setKde(null);
            setKdeError(null);
            if (localRun.structure) {
                setStructure(localRun.structure);
                setError(null);
                setLoading(false);
                return;
            }
            if (localRun.structureFile) {
                const worker = localStructureWorkerRef.current;
                if (!worker) {
                    setStructure(null);
                    setError('Browser structure parser is not available');
                    setLoading(false);
                    return;
                }
                const id = localStructureRequestRef.current + 1;
                localStructureRequestRef.current = id;
                setStructure(null);
                setError(null);
                setLoading(true);
                worker.postMessage({
                    id,
                    file: localRun.structureFile,
                    maxPoints: STRUCTURE_MAX_POINTS
                });
                return;
            }
            setStructure(null);
            setError(localRun.structureError || 'No structure data available in this folder');
            setLoading(false);
            return;
        }

        // Static mode has no Flask backend; structure data comes from a selected local run.
        if (isStaticMode()) {
            setStructure(null);
            setError('Open a run folder to view the structure.');
            setLoading(false);
            return;
        }

        const fetchStructure = async () => {
            setLoading(true);
            setError(null);
            try {
                const response = await axios.get(`${API_BASE_URL}/api/structure`, {
                    params: { dir: directory || '.', maxPoints: STRUCTURE_MAX_POINTS }
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
    }, [directory, localRun]);

    useEffect(() => {
        if (!normalMenuOpen) return undefined;
        const handlePointerDown = (event) => {
            if (!normalMenuRef.current?.contains(event.target)) {
                setNormalMenuOpen(false);
            }
        };
        const handleKeyDown = (event) => {
            if (event.key === 'Escape') {
                setNormalMenuOpen(false);
            }
        };
        document.addEventListener('pointerdown', handlePointerDown);
        document.addEventListener('keydown', handleKeyDown);
        return () => {
            document.removeEventListener('pointerdown', handlePointerDown);
            document.removeEventListener('keydown', handleKeyDown);
        };
    }, [normalMenuOpen]);

    const points = useMemo(() => {
        const allPoints = structure?.points || [];
        if (selectedElement === 'all') {
            return allPoints;
        }
        return allPoints.filter((point) => point.element === selectedElement);
    }, [structure, selectedElement]);

    // One distinct color per element, stable across element filtering and shared
    // by the slab-in-cell canvas and the folded-unit-cell 3D view.
    const elementColors = useMemo(() => {
        const elements = structure?.elements?.length
            ? structure.elements
            : (structure?.points || []).map((point) => point.element);
        return buildElementColors(elements);
    }, [structure]);

    const unitCell = useMemo(() => {
        if (!structure?.latticeVectors || !structure?.supercell) {
            const unitVectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
            return {
                lengths: [1, 1, 1],
                unitVectors,
                basis: unitVectors,
                center: [0.5, 0.5, 0.5]
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
            unitVectors,
            basis,
            center
        };
    }, [structure]);

    const sliceConfig = useMemo(() => makeSliceConfig(sliceDirection, customDirection), [sliceDirection, customDirection]);
    const slicePanelGeometry = useMemo(() => {
        const planePoints = CUBE_CORNERS.map((corner) => [dot(corner, sliceConfig.u), dot(corner, sliceConfig.v)]);
        const sidePoints = CUBE_CORNERS.map((corner) => [dot(corner, sliceConfig.u), dot(corner, sliceConfig.normal)]);
        return {
            planeAspect: makeProjectedPlane(
                vectorFromFraction(sliceConfig.u, unitCell.unitVectors),
                vectorFromFraction(sliceConfig.v, unitCell.unitVectors),
                planePoints
            ).aspect,
            sideAspect: makeProjectedPlane(
                vectorFromFraction(sliceConfig.u, unitCell.unitVectors),
                vectorFromFraction(sliceConfig.normal, unitCell.unitVectors),
                sidePoints
            ).aspect
        };
    }, [sliceConfig, unitCell]);

    const updateCustomDirection = (axisIndex, value) => {
        setCustomDirection((current) => current.map((axisValue, index) => (
            index === axisIndex ? Number(value) : axisValue
        )));
    };

    const chooseSliceDirection = (value) => {
        setSliceDirection(value);
        setNormalMenuOpen(false);
    };

    // Native PNG reads the live canvas; "png3x" re-renders the same drawing onto
    // a 3x offscreen canvas so the export is genuinely higher-resolution.
    const save2dPanel = async (canvas, drawFn, name, minW, minH, format) => {
        if (!canvas) return;
        if (format === 'png3x') {
            const rect = canvas.getBoundingClientRect();
            const width = Math.max(minW, Math.floor(rect.width));
            const height = Math.max(minH, Math.floor(rect.height));
            const scale = 3;
            const off = document.createElement('canvas');
            off.width = width * scale;
            off.height = height * scale;
            const ctx = off.getContext('2d');
            ctx.setTransform(scale, 0, 0, scale, 0, 0);
            drawFn(ctx, width, height);
            await saveCanvasAsPng(off, name);
        } else {
            await saveCanvasAsPng(canvas, name);
        }
    };

    // Re-render the three.js scene at a higher pixel ratio, capture, then restore.
    const captureModelBlob = (scale) => new Promise((resolve) => {
        const context = modelExportRef.current;
        if (!context) {
            resolve(null);
            return;
        }
        const { renderer, scene, camera, width, height } = context;
        const previousRatio = renderer.getPixelRatio();
        renderer.setPixelRatio(scale);
        renderer.setSize(width, height, false);
        renderer.render(scene, camera);
        renderer.domElement.toBlob((blob) => {
            renderer.setPixelRatio(previousRatio);
            renderer.setSize(width, height, false);
            renderer.render(scene, camera);
            resolve(blob);
        }, 'image/png');
    });

    const saveKdeSlice = (format) => save2dPanel(canvasRef.current, drawKdeSlice, `KDE_Slice_${sliceConfig.label}`, 320, 260, format);
    const saveSlab = (format) => save2dPanel(slabCanvasRef.current, drawSlab, `Slab_In_Cell_${sliceConfig.label}`, 220, 260, format);
    const saveModel = async (format) => {
        if (format === 'png3x') {
            const blob = await captureModelBlob(3);
            if (blob) downloadBlob(blob, `${sanitizeFilename('Folded_Unit_Cell')}.png`);
        } else {
            const canvas = mountRef.current?.querySelector('canvas');
            if (canvas) await saveCanvasAsPng(canvas, 'Folded_Unit_Cell');
        }
    };

    const pointDepth = useCallback((point) => {
        const rawDepth = dot([point.x, point.y, point.z], sliceConfig.normal);
        const span = sliceConfig.range[1] - sliceConfig.range[0] || 1;
        return (rawDepth - sliceConfig.range[0]) / span;
    }, [sliceConfig]);

    const inActiveSlab = useCallback(
        (point) => Math.abs(pointDepth(point) - zCenter) <= thickness / 2,
        [pointDepth, zCenter, thickness]
    );

    // Default the slice to the densest band so the view is populated on load
    // (the geometric midpoint can fall in a gap between atomic layers).
    useEffect(() => {
        if (!points.length) return;
        const bins = 50;
        const counts = new Array(bins).fill(0);
        points.forEach((point) => {
            const bin = Math.max(0, Math.min(bins - 1, Math.floor(pointDepth(point) * bins)));
            counts[bin] += 1;
        });
        let best = 0;
        counts.forEach((count, index) => {
            if (count > counts[best]) best = index;
        });
        setZCenter((best + 0.5) / bins);
    }, [points, sliceConfig, pointDepth]);

    // Fetch a real server-side gaussian_kde slice, or compute a lightweight
    // browser-side density field for uploaded static-mode data.
    useEffect(() => {
        if (!structure) return undefined;
        if (isLocalStructure) {
            const handle = setTimeout(() => {
                const worker = localKdeWorkerRef.current;
                if (!worker) return;
                const id = localKdeRequestRef.current + 1;
                localKdeRequestRef.current = id;
                setKdeLoading(true);
                setKdeError(null);
                worker.postMessage({
                    id,
                    points,
                    normal: sliceConfig.normal,
                    uVector: sliceConfig.u,
                    vVector: sliceConfig.v,
                    range: sliceConfig.range,
                    zCenter,
                    thickness,
                    bandwidth,
                    gridSize,
                    logScale
                });
            }, 80);

            return () => clearTimeout(handle);
        }

        const controller = new AbortController();
        const handle = setTimeout(async () => {
            setKdeLoading(true);
            setKdeError(null);
            try {
                const response = await axios.get(`${API_BASE_URL}/api/kde/slice`, {
                    params: {
                        dir: directory || '.',
                        element: selectedElement,
                        orientation: sliceDirection,
                        nx: sliceConfig.normal[0],
                        ny: sliceConfig.normal[1],
                        nz: sliceConfig.normal[2],
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
    }, [structure, isLocalStructure, points, directory, selectedElement, sliceDirection, sliceConfig, zCenter, thickness, bandwidth, gridSize, logScale]);

    // Render the KDE density grid + contour overlay into a 2D context sized to
    // (width, height) CSS units. The caller sets the backing resolution and the
    // matching transform, so the same draw serves the live canvas and a
    // higher-resolution offscreen canvas for export.
    const drawKdeSlice = useCallback((ctx, width, height) => {
        ctx.clearRect(0, 0, width, height);
        ctx.fillStyle = themeVars.canvasBg;
        ctx.fillRect(0, 0, width, height);

        const extent = kde?.extent || [-0.5, 0.5, -0.5, 0.5];
        const [xMin, xMax, yMin, yMax] = extent;
        const uVector = kde?.uVector || sliceConfig.u;
        const vVector = kde?.vVector || sliceConfig.v;
        const planePolygon = kde?.planePolygon?.length
            ? kde.planePolygon
            : [[xMin, yMin], [xMax, yMin], [xMax, yMax], [xMin, yMax]];
        const projectedPlane = makeProjectedPlane(
            vectorFromFraction(uVector, unitCell.unitVectors),
            vectorFromFraction(vVector, unitCell.unitVectors),
            planePolygon
        );
        const mapper = makePlaneMapper(projectedPlane, width, height, 18);
        const cellOutline = [
            ...planePolygon.map(([uValue, vValue]) => mapper.map(uValue, vValue))
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
            const origin = mapper.map(xMin, yMin);
            const uEdge = mapper.map(xMax, yMin);
            const vEdge = mapper.map(xMin, yMax);
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
                ctx.lineWidth = 1;
                ctx.strokeStyle = themeVars.contour;
                kde.contours.forEach((contour) => {
                    contour.lines.forEach((line) => {
                        ctx.beginPath();
                        line.forEach(([dataX, dataY], index) => {
                            const point = mapper.map(dataX, dataY);
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
            drawOverlayText(`${sliceConfig.label}=${kde.center.toFixed(3)}  d=${kde.thickness.toFixed(3)}  bw=${kde.bw}`, 12, 40);
            if (kde.log) drawOverlayText('log10 density', 12, 58);
        }
    }, [kde, colormap, showContours, kdeLoading, themeVars, unitCell, sliceConfig]);

    // Colormap and contour visibility are pure client-side re-renders (no refetch).
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
        drawKdeSlice(ctx, width, height);
    }, [drawKdeSlice]);

    // Draw the "Slab In Cell" side view into a 2D context sized to (width,
    // height) CSS units. Returns the band geometry the drag handler needs; the
    // live effect publishes it, export ignores it. Like drawKdeSlice, the caller
    // owns the backing resolution and transform.
    const drawSlab = useCallback((ctx, width, height) => {
        ctx.clearRect(0, 0, width, height);
        ctx.fillStyle = themeVars.canvasBg;
        ctx.fillRect(0, 0, width, height);

        const uVector = vectorFromFraction(sliceConfig.u, unitCell.unitVectors);
        const normalVector = vectorFromFraction(sliceConfig.normal, unitCell.unitVectors);
        const projectedCorners = CUBE_CORNERS.map((corner) => [dot(corner, sliceConfig.u), dot(corner, sliceConfig.normal)]);
        const depthSpan = sliceConfig.range[1] - sliceConfig.range[0] || 1;
        const centerDepth = sliceConfig.range[0] + zCenter * depthSpan;
        const depthThickness = thickness * depthSpan;
        const depthStart = Math.max(sliceConfig.range[0], centerDepth - depthThickness / 2);
        const depthEnd = Math.min(sliceConfig.range[1], centerDepth + depthThickness / 2);
        const uValues = projectedCorners.map(([uValue]) => uValue);
        const uMin = Math.min(...uValues);
        const uMax = Math.max(...uValues);
        const sidePlane = makeProjectedPlane(uVector, normalVector, [
            ...projectedCorners,
            [uMin, depthStart],
            [uMax, depthStart],
            [uMax, depthEnd],
            [uMin, depthEnd]
        ]);
        const mapper = makePlaneMapper(sidePlane, width, height, 18);
        // Band geometry so the drag handler can map cursor -> slice (returned below).
        const geometry = {
            invert: mapper.invert,
            uMin,
            uMax,
            depthStart,
            depthEnd,
            rangeMin: sliceConfig.range[0],
            depthSpan
        };
        const cellOutline = [
            mapper.map(uMin, sliceConfig.range[0]),
            mapper.map(uMax, sliceConfig.range[0]),
            mapper.map(uMax, sliceConfig.range[1]),
            mapper.map(uMin, sliceConfig.range[1])
        ];
        const slabOutline = [
            mapper.map(uMin, depthStart),
            mapper.map(uMax, depthStart),
            mapper.map(uMax, depthEnd),
            mapper.map(uMin, depthEnd)
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

        const sampleLimit = Math.min(points.length, SLAB_CANVAS_MAX_POINTS);
        const stride = Math.max(1, Math.floor(points.length / sampleLimit));
        for (let index = 0; index < points.length; index += stride) {
            const point = points[index];
            const fraction = [point.x, point.y, point.z];
            const projected = mapper.map(dot(fraction, sliceConfig.u), dot(fraction, sliceConfig.normal));
            const inSlab = inActiveSlab(point);
            ctx.fillStyle = inSlab ? (elementColors[point.element] || DEFAULT_ELEMENT_COLOR) : 'rgba(166, 176, 188, 0.22)';
            ctx.fillRect(projected.x, projected.y, inSlab ? 2 : 1, inSlab ? 2 : 1);
        }

        ctx.fillStyle = themeVars.text;
        ctx.font = '500 12px Inter, system-ui';
        const xLabel = mapper.map(uMax, sliceConfig.range[0]);
        const normalLabel = mapper.map(uMin, sliceConfig.range[1]);
        const slabLabel = mapper.map(uMin, depthStart);
        ctx.fillText(sliceConfig.uLabel, Math.min(width - 24, xLabel.x + 4), Math.min(height - 8, xLabel.y + 14));
        ctx.fillText(sliceConfig.label, Math.max(8, normalLabel.x - 12), Math.max(16, normalLabel.y - 6));
        ctx.fillText(`${sliceConfig.label}=${zCenter.toFixed(3)}`, 10, Math.max(30, slabLabel.y - 6));
        ctx.fillText(`d=${thickness.toFixed(3)}`, 10, Math.min(height - 16, slabLabel.y + 18));
        return geometry;
    }, [points, zCenter, thickness, unitCell, themeVars, sliceConfig, inActiveSlab, elementColors]);

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
        // Publish the band geometry so the drag handler can map cursor -> slice.
        slabGeometryRef.current = drawSlab(ctx, width, height);
    }, [drawSlab]);

    useEffect(() => {
        zCenterRef.current = zCenter;
    }, [zCenter]);

    // Let users drag the highlighted band in "Slab In Cell" to move the slice.
    useEffect(() => {
        const canvas = slabCanvasRef.current;
        if (!canvas) return undefined;

        const planeCoordsAt = (event) => {
            const geometry = slabGeometryRef.current;
            if (!geometry) return null;
            const rect = canvas.getBoundingClientRect();
            return geometry.invert(event.clientX - rect.left, event.clientY - rect.top);
        };

        const overBand = (coords) => {
            const geometry = slabGeometryRef.current;
            if (!geometry || !coords) return false;
            return (
                coords.uFraction >= geometry.uMin &&
                coords.uFraction <= geometry.uMax &&
                coords.vFraction >= geometry.depthStart &&
                coords.vFraction <= geometry.depthEnd
            );
        };

        const zCenterAt = (coords) => {
            const geometry = slabGeometryRef.current;
            return (coords.vFraction - geometry.rangeMin) / geometry.depthSpan;
        };

        const handlePointerDown = (event) => {
            const coords = planeCoordsAt(event);
            if (!overBand(coords)) return;
            slabDragRef.current = { offset: zCenterRef.current - zCenterAt(coords) };
            canvas.setPointerCapture(event.pointerId);
            canvas.style.cursor = 'grabbing';
            event.preventDefault();
        };

        const handlePointerMove = (event) => {
            const coords = planeCoordsAt(event);
            if (!coords) return;
            if (slabDragRef.current) {
                const next = zCenterAt(coords) + slabDragRef.current.offset;
                setZCenter(Math.min(1, Math.max(0, next)));
                event.preventDefault();
            } else {
                canvas.style.cursor = overBand(coords) ? 'grab' : 'default';
            }
        };

        const endDrag = (event) => {
            if (!slabDragRef.current) return;
            slabDragRef.current = null;
            if (canvas.hasPointerCapture?.(event.pointerId)) {
                canvas.releasePointerCapture(event.pointerId);
            }
            canvas.style.cursor = overBand(planeCoordsAt(event)) ? 'grab' : 'default';
        };

        const handlePointerLeave = () => {
            if (!slabDragRef.current) canvas.style.cursor = 'default';
        };

        canvas.addEventListener('pointerdown', handlePointerDown);
        canvas.addEventListener('pointermove', handlePointerMove);
        canvas.addEventListener('pointerup', endDrag);
        canvas.addEventListener('pointercancel', endDrag);
        canvas.addEventListener('pointerleave', handlePointerLeave);

        return () => {
            canvas.removeEventListener('pointerdown', handlePointerDown);
            canvas.removeEventListener('pointermove', handlePointerMove);
            canvas.removeEventListener('pointerup', endDrag);
            canvas.removeEventListener('pointercancel', endDrag);
            canvas.removeEventListener('pointerleave', handlePointerLeave);
        };
    // Re-attach once the slab canvas mounts (it only renders after a structure loads).
    }, [structure]);

    useEffect(() => {
        const mount = mountRef.current;
        if (!mount || points.length === 0) return undefined;

        const width = mount.clientWidth;
        const height = mount.clientHeight;
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(themeVars.canvasBg);
        const camera = new THREE.PerspectiveCamera(45, width / height, 0.01, 20);

        // preserveDrawingBuffer keeps the rendered frame readable so the figure
        // can be captured for PNG export at any time.
        const renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
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
                color: elementColors[element] || DEFAULT_ELEMENT_COLOR,
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

        const depthSpan = sliceConfig.range[1] - sliceConfig.range[0] || 1;
        const centerDepth = sliceConfig.range[0] + zCenter * depthSpan;
        const depthThickness = thickness * depthSpan;
        const depthStart = Math.max(sliceConfig.range[0], centerDepth - depthThickness / 2);
        const depthEnd = Math.min(sliceConfig.range[1], centerDepth + depthThickness / 2);
        const slabSections = [depthStart, depthEnd].map((depth) => (
            planeSectionVertices(sliceConfig.normal, depth).map((fraction) => (
                toVector3(subtract(vectorFromFraction(fraction, unitCell.basis), unitCell.center))
            ))
        ));
        const slabGeometry = makeSlabGeometry(slabSections);
        const slabMaterial = new THREE.MeshBasicMaterial({
            color: '#4f8cff',
            transparent: true,
            opacity: 0.12,
            side: THREE.DoubleSide,
            depthWrite: false
        });
        const slabMesh = new THREE.Mesh(slabGeometry, slabMaterial);
        scene.add(slabMesh);

        const slabEdgeGeometry = makeSectionEdgeGeometry(slabSections);
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
        controls.minDistance = radius * 0.35;
        controls.maxDistance = radius * 8;
        if (cameraStateRef.current) {
            camera.position.fromArray(cameraStateRef.current.position);
            controls.target.fromArray(cameraStateRef.current.target);
            camera.zoom = cameraStateRef.current.zoom;
        } else {
            camera.position.copy(sphere.center).add(new THREE.Vector3(radius * 1.7, radius * 1.45, radius * 1.55));
            controls.target.copy(sphere.center);
        }
        camera.lookAt(controls.target);
        camera.updateProjectionMatrix();

        // Expose the render context so the figure can be re-rendered at a higher
        // pixel ratio for export.
        modelExportRef.current = { renderer, scene, camera, width, height };

        let frameId;
        const animate = () => {
            controls.update();
            renderer.render(scene, camera);
            frameId = requestAnimationFrame(animate);
        };
        animate();

        return () => {
            cameraStateRef.current = {
                position: camera.position.toArray(),
                target: controls.target.toArray(),
                zoom: camera.zoom
            };
            modelExportRef.current = null;
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
    }, [points, unitCell, zCenter, thickness, themeVars, sliceConfig, elementColors]);

    return (
        <section className="structure-page">
            <div className="structure-header">
                <div>
                    <h2>KDE And Folded Unit Cell</h2>
                    <p>{localRun ? localRun.name : directory}</p>
                </div>
                {(loading || kdeLoading) && <span className="status-pill">{loading ? 'Loading' : isLocalStructure ? 'Density' : 'KDE'}</span>}
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
                            <span className="control-name">Normal</span>
                            <span className="normal-menu" ref={normalMenuRef}>
                                <button
                                    type="button"
                                    className="normal-menu-button"
                                    aria-haspopup="listbox"
                                    aria-expanded={normalMenuOpen}
                                    onClick={() => setNormalMenuOpen((open) => !open)}
                                >
                                    {NORMAL_OPTIONS.find((option) => option.value === sliceDirection)?.label || 'c'}
                                </button>
                                {normalMenuOpen && (
                                    <span className="normal-menu-list" role="listbox" aria-label="Normal">
                                        {NORMAL_OPTIONS.map((option) => (
                                            <button
                                                key={option.value}
                                                type="button"
                                                role="option"
                                                aria-selected={sliceDirection === option.value}
                                                className={sliceDirection === option.value ? 'is-selected' : ''}
                                                onClick={() => chooseSliceDirection(option.value)}
                                            >
                                                {option.label}
                                            </button>
                                        ))}
                                    </span>
                                )}
                            </span>
                        </label>
                        {sliceDirection === 'custom' && (
                            <label className="control custom-direction">
                                <span className="control-name">Direction</span>
                                {customDirection.map((value, index) => (
                                    <input
                                        key={CUSTOM_DIRECTION_LABELS[index]}
                                        type="number"
                                        step="0.1"
                                        value={value}
                                        aria-label={`Direction ${CUSTOM_DIRECTION_LABELS[index]}`}
                                        onChange={(event) => updateCustomDirection(index, event.target.value)}
                                    />
                                ))}
                            </label>
                        )}
                        <label className="control">
                            <span className="control-name">Slice</span>
                            <input
                                type="range"
                                min="0"
                                max="1"
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
                            style={{ '--panel-aspect': slicePanelGeometry.planeAspect }}
                        >
                            <h3>
                                KDE Slice
                                <SaveMenu onSave={saveKdeSlice} options={PANEL_SAVE_OPTIONS} label="Save" align="right" />
                            </h3>
                            <canvas ref={canvasRef} className="kde-canvas" />
                            {isLocalStructure && (
                                <div className="local-density-note">
                                    Browser-side Gaussian KDE. The Flask app uses SciPy KDE for reference-grade values.
                                </div>
                            )}
                        </div>
                        <div
                            className="slab-panel"
                            style={{ '--panel-aspect': slicePanelGeometry.sideAspect }}
                        >
                            <div className="slab-panel-header">
                                <span>Slab In Cell</span>
                                <div className="slab-panel-header-meta">
                                    <strong>{sliceConfig.label} {(Math.max(0, zCenter - thickness / 2)).toFixed(2)} - {(Math.min(1, zCenter + thickness / 2)).toFixed(2)}</strong>
                                    <SaveMenu onSave={saveSlab} options={PANEL_SAVE_OPTIONS} label="Save" align="right" />
                                </div>
                            </div>
                            <canvas ref={slabCanvasRef} />
                        </div>
                        <div
                            className="model-panel"
                            style={{ '--panel-aspect': Math.max(slicePanelGeometry.planeAspect, 1) }}
                        >
                            <h3>
                                Folded Unit Cell
                                <SaveMenu onSave={saveModel} options={PANEL_SAVE_OPTIONS} label="Save" align="right" />
                            </h3>
                            <div ref={mountRef} className="three-mount" />
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
                    </div>
                </>
            )}
            
            <footer className="app-footer">
                &copy; 2026 Tsung-Han Yang &middot;{' '}
                <a
                    href="https://github.com/drthyang/rmc-toolkits/blob/main/LICENSE"
                    target="_blank"
                    rel="noreferrer"
                >
                    MIT License
                </a>
            </footer>
        </section>
    );
};

export default StructurePage;
