import React, { useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import API_BASE_URL from '../api';
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

    useEffect(() => {
        if (!points.length) return;
        const zs = points.map((point) => point.z).sort((a, b) => a - b);
        setZCenter(zs[Math.floor(zs.length / 2)]);
    }, [points]);

    useEffect(() => {
        const canvas = canvasRef.current;
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const size = canvas.getBoundingClientRect();
        const width = Math.max(320, Math.floor(size.width));
        const height = Math.max(260, Math.floor(size.height));
        canvas.width = width * window.devicePixelRatio;
        canvas.height = height * window.devicePixelRatio;
        ctx.scale(window.devicePixelRatio, window.devicePixelRatio);
        ctx.clearRect(0, 0, width, height);
        ctx.fillStyle = '#111316';
        ctx.fillRect(0, 0, width, height);

        const slab = points.filter((point) => Math.abs(point.z - zCenter) <= thickness / 2);
        const grid = 220;
        const density = Array.from({ length: grid }, () => Array(grid).fill(0));
        slab.forEach((point) => {
            const gx = Math.max(0, Math.min(grid - 1, Math.floor(point.x * grid)));
            const gy = Math.max(0, Math.min(grid - 1, Math.floor(point.y * grid)));
            for (let dx = -5; dx <= 5; dx += 1) {
                for (let dy = -5; dy <= 5; dy += 1) {
                    const x = gx + dx;
                    const y = gy + dy;
                    if (x >= 0 && y >= 0 && x < grid && y < grid) {
                        const distanceSq = dx * dx + dy * dy;
                        density[y][x] += Math.exp(-distanceSq / 10);
                    }
                }
            }
        });

        const maxDensity = Math.max(1, ...density.flat());
        const cellW = width / grid;
        const cellH = height / grid;
        for (let y = 0; y < grid; y += 1) {
            for (let x = 0; x < grid; x += 1) {
                const value = density[y][x] / maxDensity;
                if (value > 0) {
                    ctx.fillStyle = `rgba(${Math.floor(65 + value * 190)}, ${Math.floor(125 + value * 80)}, ${Math.floor(135 + value * 70)}, ${Math.min(0.9, 0.15 + value)})`;
                    ctx.fillRect(x * cellW, height - (y + 1) * cellH, cellW + 1, cellH + 1);
                }
            }
        }
        ctx.strokeStyle = '#384048';
        ctx.strokeRect(0.5, 0.5, width - 1, height - 1);
        ctx.fillStyle = '#dce3ea';
        ctx.font = '12px system-ui';
        ctx.fillText(`${slab.length} atoms in slab`, 12, 22);
        ctx.fillText(`a:b = ${unitCell.lengths[0].toFixed(3)}:${unitCell.lengths[1].toFixed(3)}`, 12, 40);
    }, [points, zCenter, thickness, unitCell]);

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

        const pad = 28;
        const plotX = pad;
        const plotY = 16;
        const plotW = width - pad * 2;
        const plotH = height - pad * 2;
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
    }, [points, zCenter, thickness]);

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
                {loading && <span className="status-pill">Loading</span>}
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
                    </div>

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
                                    <canvas ref={slabCanvasRef} />
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
