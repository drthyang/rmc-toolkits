import React, { useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import API_BASE_URL from '../api';
import { COLORMAP_NAMES, getLut, sampleColormap } from '../colormaps';
import './StructurePage.css';

const colors = {
    Ga: '#4f8cff',
    Nb: '#f2b84b',
    Se: '#62c084',
    default: '#d5d9df'
};

const StructurePage = ({ directory }) => {
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
            return { lengths: [1, 1, 1], scale: [1, 1, 1], xyAspect: 1 };
        }
        const lengths = structure.latticeVectors.map((vector, index) => {
            const vectorLength = Math.sqrt(vector.reduce((sum, value) => sum + value * value, 0));
            return vectorLength / structure.supercell[index];
        });
        const maxLength = Math.max(...lengths, 1e-9);
        return {
            lengths,
            scale: lengths.map((length) => length / maxLength),
            xyAspect: lengths[0] / Math.max(lengths[1], 1e-9)
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
        ctx.fillStyle = '#0b0c0e';
        ctx.fillRect(0, 0, width, height);

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
                // Flip rows: image row 0 is the top (max y), density row 0 is min y.
                const densityRow = density[grid - 1 - py];
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
            ctx.drawImage(offscreen, 0, 0, grid, grid, 0, 0, width, height);

            if (showContours && kde.contours?.length) {
                const [xMin, xMax, yMin, yMax] = kde.extent;
                const spanX = xMax - xMin || 1;
                const spanY = yMax - yMin || 1;
                ctx.lineWidth = 1;
                ctx.strokeStyle = 'rgba(40, 44, 48, 0.85)';
                kde.contours.forEach((contour) => {
                    contour.lines.forEach((line) => {
                        ctx.beginPath();
                        line.forEach(([dataX, dataY], index) => {
                            const cx = ((dataX - xMin) / spanX) * width;
                            const cy = height - ((dataY - yMin) / spanY) * height;
                            if (index === 0) ctx.moveTo(cx, cy);
                            else ctx.lineTo(cx, cy);
                        });
                        ctx.stroke();
                    });
                });
            }
        } else {
            ctx.fillStyle = '#7c858f';
            ctx.font = '13px system-ui';
            ctx.fillText(kdeLoading ? 'Computing KDE...' : 'No atoms in this slab', 14, 28);
        }

        ctx.strokeStyle = '#384048';
        ctx.strokeRect(0.5, 0.5, width - 1, height - 1);
        ctx.fillStyle = '#dce3ea';
        ctx.font = '12px system-ui';
        if (kde) {
            ctx.fillText(`${kde.slabCount} atoms in slab (fit ${kde.fitCount})`, 12, 22);
            ctx.fillText(`z=${kde.z.toFixed(2)} A  dz=${kde.dz.toFixed(2)} A  bw=${kde.bw}`, 12, 40);
            if (kde.log) ctx.fillText('log10 density', 12, 58);
        }
    }, [kde, colormap, showContours, kdeLoading]);

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
        ctx.fillStyle = '#111316';
        ctx.fillRect(0, 0, width, height);

        // Size the plot box so equal distances in x and z map to equal pixels,
        // i.e. width:height matches the in-plane lattice parameters a:c.
        const aspect = unitCell.lengths[0] / Math.max(unitCell.lengths[2], 1e-9);
        const padX = 28;
        const padTop = 16;
        const padBottom = 24;
        const availW = width - padX * 2;
        const availH = height - padTop - padBottom;
        let plotW = availW;
        let plotH = plotW / aspect;
        if (plotH > availH) {
            plotH = availH;
            plotW = plotH * aspect;
        }
        const plotX = padX + (availW - plotW) / 2;
        const plotY = padTop + (availH - plotH) / 2;
        const zStart = Math.max(0, zCenter - thickness / 2);
        const zEnd = Math.min(1, zCenter + thickness / 2);
        const bandTop = plotY + (1 - zEnd) * plotH;
        const bandHeight = Math.max(2, (zEnd - zStart) * plotH);

        ctx.strokeStyle = '#48515b';
        ctx.lineWidth = 1;
        ctx.strokeRect(plotX, plotY, plotW, plotH);

        ctx.fillStyle = 'rgba(79, 140, 255, 0.18)';
        ctx.fillRect(plotX, bandTop, plotW, bandHeight);
        ctx.strokeStyle = '#74a7ff';
        ctx.strokeRect(plotX, bandTop, plotW, bandHeight);

        const sampleLimit = Math.min(points.length, 9000);
        const stride = Math.max(1, Math.floor(points.length / sampleLimit));
        for (let index = 0; index < points.length; index += stride) {
            const point = points[index];
            const x = plotX + point.x * plotW;
            const y = plotY + (1 - point.z) * plotH;
            const inSlab = Math.abs(point.z - zCenter) <= thickness / 2;
            ctx.fillStyle = inSlab ? (colors[point.element] || colors.default) : 'rgba(166, 176, 188, 0.22)';
            ctx.fillRect(x, y, inSlab ? 2 : 1, inSlab ? 2 : 1);
        }

        ctx.fillStyle = '#dce3ea';
        ctx.font = '12px system-ui';
        ctx.fillText('x', plotX + plotW - 4, height - 7);
        ctx.fillText('z', 8, plotY + 8);
        ctx.fillText(`z=${zCenter.toFixed(3)}`, plotX + 8, Math.max(30, bandTop - 6));
        ctx.fillText(`dz=${thickness.toFixed(3)}`, plotX + 8, Math.min(height - 16, bandTop + bandHeight + 16));
    }, [points, zCenter, thickness, unitCell]);

    useEffect(() => {
        const mount = mountRef.current;
        if (!mount || points.length === 0) return undefined;

        const width = mount.clientWidth;
        const height = mount.clientHeight;
        const scene = new THREE.Scene();
        scene.background = new THREE.Color('#101113');
        const camera = new THREE.PerspectiveCamera(45, width / height, 0.01, 10);
        camera.position.set(1.8, 1.6, 1.9);
        camera.lookAt(0, 0, 0);

        const renderer = new THREE.WebGLRenderer({ antialias: true });
        renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        renderer.setSize(width, height);
        mount.replaceChildren(renderer.domElement);

        const controls = new OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.dampingFactor = 0.08;
        controls.enablePan = true;
        controls.minDistance = 0.7;
        controls.maxDistance = 5;

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
                positions[index * 3] = (point.x - 0.5) * unitCell.scale[0];
                positions[index * 3 + 1] = (point.y - 0.5) * unitCell.scale[1];
                positions[index * 3 + 2] = (point.z - 0.5) * unitCell.scale[2];
            });
            geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
            const material = new THREE.PointsMaterial({
                color: colors[element] || colors.default,
                size: 0.018,
                sizeAttenuation: true
            });
            group.add(new THREE.Points(geometry, material));
        });

        const box = new THREE.Box3(
            new THREE.Vector3(-0.5 * unitCell.scale[0], -0.5 * unitCell.scale[1], -0.5 * unitCell.scale[2]),
            new THREE.Vector3(0.5 * unitCell.scale[0], 0.5 * unitCell.scale[1], 0.5 * unitCell.scale[2])
        );
        const helper = new THREE.Box3Helper(box, new THREE.Color('#737c86'));
        scene.add(helper);

        const slabZ = (zCenter - 0.5) * unitCell.scale[2];
        const slabThickness = thickness * unitCell.scale[2];
        const slabGeometry = new THREE.BoxGeometry(unitCell.scale[0], unitCell.scale[1], slabThickness);
        const slabMaterial = new THREE.MeshBasicMaterial({
            color: '#4f8cff',
            transparent: true,
            opacity: 0.12,
            depthWrite: false
        });
        const slabMesh = new THREE.Mesh(slabGeometry, slabMaterial);
        slabMesh.position.set(0, 0, slabZ);
        scene.add(slabMesh);

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
            group.traverse((object) => {
                object.geometry?.dispose?.();
                object.material?.dispose?.();
            });
        };
    }, [points, unitCell, zCenter, thickness]);

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
                    <div className="structure-summary">
                        <div><span>Total atoms</span><strong>{structure.totalAtoms}</strong></div>
                        <div><span>Rendered atoms</span><strong>{points.length}</strong></div>
                        <div><span>Unit cell</span><strong>{unitCell.lengths.map((value) => value.toFixed(3)).join(' x ')} A</strong></div>
                        <div><span>Source</span><strong>{structure.source.split('/').pop()}</strong></div>
                    </div>

                    <div className="structure-controls">
                        <label>
                            Element
                            <select value={selectedElement} onChange={(event) => setSelectedElement(event.target.value)}>
                                <option value="all">All</option>
                                {structure.elements.map((element) => (
                                    <option key={element} value={element}>{element}</option>
                                ))}
                            </select>
                        </label>
                        <label>
                            z
                            <input
                                type="range"
                                min={ranges.z[0]}
                                max={ranges.z[1]}
                                step="0.001"
                                value={zCenter}
                                onChange={(event) => setZCenter(Number(event.target.value))}
                            />
                            <span>{zCenter.toFixed(2)}</span>
                        </label>
                        <label>
                            dz
                            <input type="range" min="0.01" max="0.5" step="0.01" value={thickness} onChange={(event) => setThickness(Number(event.target.value))} />
                            <span>{thickness.toFixed(2)}</span>
                        </label>
                        <label>
                            bw
                            <input type="range" min="0.005" max="0.15" step="0.005" value={bandwidth} onChange={(event) => setBandwidth(Number(event.target.value))} />
                            <span>{bandwidth.toFixed(3)}</span>
                        </label>
                        <label>
                            cmap
                            <select value={colormap} onChange={(event) => setColormap(event.target.value)}>
                                {COLORMAP_NAMES.map((name) => (
                                    <option key={name} value={name}>{name}</option>
                                ))}
                            </select>
                        </label>
                        <label>
                            grid
                            <select value={gridSize} onChange={(event) => setGridSize(Number(event.target.value))}>
                                <option value={80}>80</option>
                                <option value={120}>120</option>
                                <option value={160}>160</option>
                                <option value={220}>220</option>
                            </select>
                        </label>
                        <label className="checkbox">
                            <input type="checkbox" checked={showContours} onChange={(event) => setShowContours(event.target.checked)} />
                            contours
                        </label>
                        <label className="checkbox">
                            <input type="checkbox" checked={logScale} onChange={(event) => setLogScale(event.target.checked)} />
                            log
                        </label>
                    </div>

                    {kdeError && <div className="structure-error">{kdeError}</div>}

                    <div className="analysis-layout">
                        <div className="kde-panel">
                            <h3>KDE Slice</h3>
                            <div className="kde-layout">
                                <canvas
                                    ref={canvasRef}
                                    className="kde-canvas"
                                    style={{ aspectRatio: `${unitCell.lengths[0]} / ${unitCell.lengths[1]}` }}
                                />
                                <div className="slab-panel">
                                    <div className="slab-panel-header">
                                        <span>Slab In Cell</span>
                                        <strong>{(Math.max(0, zCenter - thickness / 2)).toFixed(2)} - {(Math.min(1, zCenter + thickness / 2)).toFixed(2)}</strong>
                                    </div>
                                    <canvas
                                        ref={slabCanvasRef}
                                        style={{ aspectRatio: `${unitCell.lengths[0]} / ${unitCell.lengths[2]}` }}
                                    />
                                </div>
                            </div>
                        </div>
                        <div className="model-panel">
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
