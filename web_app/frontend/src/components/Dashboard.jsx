import React, { useEffect, useMemo, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import InteractivePlot from './InteractivePlot';
import './Dashboard.css';

const plotOrder = ['r_value', 'bragg', 'xray_sq', 'neutron_sq', 'xpdf', 'npdf', 'pdf_partials', 'stog'];

const vectorLength = (vector) => Math.sqrt(vector.reduce((sum, value) => sum + value * value, 0));

const angleBetween = (a, b) => {
    const denominator = Math.max(vectorLength(a) * vectorLength(b), 1e-12);
    const cosine = a.reduce((sum, value, index) => sum + value * b[index], 0) / denominator;
    const clamped = Math.max(-1, Math.min(1, cosine));
    return Math.acos(clamped) * (180 / Math.PI);
};

const formatNumber = (value, digits = 3) => Number(value).toLocaleString(undefined, {
    maximumFractionDigits: digits
});

const Dashboard = ({ directory }) => {
    const [files, setFiles] = useState([]);
    const [metadata, setMetadata] = useState({});
    const [structure, setStructure] = useState(null);
    const [structureError, setStructureError] = useState(null);
    const [error, setError] = useState(null);
    const [loading, setLoading] = useState(false);

    useEffect(() => {
        const fetchDashboard = async () => {
            setLoading(true);
            setError(null);
            try {
                const response = await axios.get(`${API_BASE_URL}/api/files`, {
                    params: { dir: directory || '.' }
                });
                const loadedFiles = response.data.files || [];
                setFiles(loadedFiles);
                const plotFiles = loadedFiles.filter((file) => file.plotKind);
                const metadataEntries = await Promise.all(
                    plotFiles.map(async (file) => {
                        try {
                            const meta = await axios.get(`${API_BASE_URL}/api/plot/metadata`, {
                                params: { path: file.path }
                            });
                            return [file.path, meta.data];
                        } catch {
                            return [file.path, null];
                        }
                    })
                );
                setMetadata(Object.fromEntries(metadataEntries));

                try {
                    const structureResponse = await axios.get(`${API_BASE_URL}/api/structure`, {
                        params: { dir: directory || '.', maxPoints: 100 }
                    });
                    setStructure(structureResponse.data);
                    setStructureError(null);
                } catch (structureErr) {
                    setStructure(null);
                    setStructureError(structureErr.response?.data?.error || 'No model structure detected');
                }
            } catch (err) {
                setError(err.response?.data?.error || 'Failed to load dashboard data');
                setStructure(null);
                setStructureError(null);
            } finally {
                setLoading(false);
            }
        };

        fetchDashboard();
    }, [directory]);

    const plotFiles = useMemo(() => {
        return files
            .filter((file) => file.plotKind)
            .sort((a, b) => plotOrder.indexOf(a.plotKind) - plotOrder.indexOf(b.plotKind));
    }, [files]);

    const modelSummary = useMemo(() => {
        if (!structure?.latticeVectors || !structure?.supercell) return null;

        const boxLengths = structure.latticeVectors.map(vectorLength);
        const cellLengths = boxLengths.map((length, index) => length / Math.max(structure.supercell[index], 1));
        const angles = [
            angleBetween(structure.latticeVectors[1], structure.latticeVectors[2]),
            angleBetween(structure.latticeVectors[0], structure.latticeVectors[2]),
            angleBetween(structure.latticeVectors[0], structure.latticeVectors[1])
        ];
        const elementEntries = Object.entries(structure.elementCounts || {})
            .sort(([a], [b]) => a.localeCompare(b))
            .map(([element, count]) => ({
                element,
                count,
                referenceSites: structure.atomIndices?.[element]?.length || 0
            }));
        const referenceSites = Object.values(structure.atomIndices || {}).reduce((sum, indices) => sum + indices.length, 0);

        return {
            source: structure.source?.split('/').pop() || 'Structure file',
            totalAtoms: structure.totalAtoms,
            referenceSites,
            supercell: structure.supercell,
            cellLengths,
            boxLengths,
            angles,
            elementEntries
        };
    }, [structure]);

    return (
        <section className="dashboard-page">
            <div className="dashboard-toolbar">
                <div>
                    <h2>Run Dashboard</h2>
                    <p>{directory}</p>
                </div>
                {loading && <span className="status-pill">Loading</span>}
            </div>

            {error && <div className="dashboard-error">{error}</div>}

            {modelSummary && (
                <section className="model-summary" aria-label="Model information">
                    <div className="model-summary-header">
                        <div>
                            <h2>Model information</h2>
                            <p>{modelSummary.source}</p>
                        </div>
                        <span>{formatNumber(modelSummary.totalAtoms, 0)} atoms</span>
                    </div>
                    <div className="model-stats">
                        <div className="stat">
                            <span className="stat-label">Lattice a · b · c (Å)</span>
                            <span className="stat-value">
                                {modelSummary.cellLengths.map((value) => formatNumber(value)).join('  ·  ')}
                            </span>
                        </div>
                        <div className="stat">
                            <span className="stat-label">Angles α · β · γ (°)</span>
                            <span className="stat-value">
                                {modelSummary.angles.map((value) => formatNumber(value, 1)).join('  ·  ')}
                            </span>
                        </div>
                        <div className="stat">
                            <span className="stat-label">Supercell</span>
                            <span className="stat-value">
                                {modelSummary.supercell.map((value) => formatNumber(value, 0)).join(' × ')}
                            </span>
                        </div>
                        <span className="stat-divider" aria-hidden="true" />
                        {modelSummary.elementEntries.map(({ element, referenceSites, count }) => (
                            <div className="stat" key={element}>
                                <span className="stat-label">{element}</span>
                                <span className="stat-value">{formatNumber(count, 0)}</span>
                                <span className="stat-sub">{formatNumber(referenceSites, 0)} sites</span>
                            </div>
                        ))}
                    </div>
                </section>
            )}

            {!modelSummary && structureError && (
                <div className="model-summary-empty">{structureError}</div>
            )}

            <div className="plot-grid">
                {plotFiles.map((file) => {
                    const meta = metadata[file.path];
                    return (
                        <article className="plot-card" key={file.path}>
                            <div className="plot-card-header">
                                <div>
                                    <h3>{meta?.title || file.name}</h3>
                                    <p>{file.plotKind}</p>
                                </div>
                                {meta?.metrics?.rwp !== undefined && (
                                    <span className="rwp-chip">Rwp {Number(meta.metrics.rwp).toPrecision(4)}</span>
                                )}
                            </div>
                            <InteractivePlot file={file} />
                        </article>
                    );
                })}
            </div>

            {!loading && plotFiles.length === 0 && (
                <div className="empty-state">Open a run folder such as data to populate the dashboard.</div>
            )}
        </section>
    );
};

export default Dashboard;
