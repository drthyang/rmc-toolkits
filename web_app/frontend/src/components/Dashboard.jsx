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
    const [showRValue, setShowRValue] = useState(false);

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
            angles,
            elementEntries
        };
    }, [structure]);

    const rValueFile = useMemo(
        () => plotFiles.find((file) => file.plotKind === 'r_value') || null,
        [plotFiles]
    );
    const gridFiles = useMemo(
        () => plotFiles.filter((file) => file.plotKind !== 'r_value'),
        [plotFiles]
    );

    const renderPlotCard = (file) => {
        const meta = metadata[file.path];
        return (
            <article className="plot-card" key={file.path}>
                <div className="plot-card-header">
                    <h3>{meta?.title || file.name}</h3>
                    {meta?.metrics?.rwp !== undefined && (
                        <span className="rwp-chip">Rwp {Number(meta.metrics.rwp).toPrecision(4)}</span>
                    )}
                </div>
                <InteractivePlot file={file} />
            </article>
        );
    };

    const renderRValuePanel = () => {
        if (!rValueFile) return null;
        const meta = metadata[rValueFile.path];
        return (
            <article className={`plot-card r-value-card${showRValue ? '' : ' is-collapsed'}`}>
                <div className="plot-card-header">
                    <h3>{meta?.title || rValueFile.name}</h3>
                    <div className="plot-card-header-actions">
                        {meta?.metrics?.rwp !== undefined && (
                            <span className="rwp-chip">Rwp {Number(meta.metrics.rwp).toPrecision(4)}</span>
                        )}
                        <button
                            type="button"
                            className="panel-toggle"
                            onClick={() => setShowRValue((value) => !value)}
                            aria-expanded={showRValue}
                        >
                            {showRValue ? 'Hide' : 'Show'}
                        </button>
                    </div>
                </div>
                {showRValue && <InteractivePlot file={rValueFile} variant="wide" />}
            </article>
        );
    };

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
                    <h2 className="model-summary-title" title={modelSummary.source}>Model information</h2>
                    <dl className="model-stats">
                        <div className="model-stat">
                            <dt>Cell (Å)</dt>
                            <dd>{modelSummary.cellLengths.map((value) => formatNumber(value)).join(' × ')}</dd>
                        </div>
                        <div className="model-stat">
                            <dt>Angles</dt>
                            <dd>{modelSummary.angles.map((value) => `${formatNumber(value, 1)}°`).join(' · ')}</dd>
                        </div>
                        <div className="model-stat">
                            <dt>Supercell</dt>
                            <dd>{modelSummary.supercell.map((value) => formatNumber(value, 0)).join(' × ')}</dd>
                        </div>
                        {modelSummary.elementEntries.map(({ element, count, referenceSites }) => (
                            <div className="model-stat" key={element}>
                                <dt>{element}</dt>
                                <dd>
                                    {formatNumber(count, 0)}
                                    <span className="model-stat-sub">{formatNumber(referenceSites, 0)} sites</span>
                                </dd>
                            </div>
                        ))}
                        <div className="model-stat">
                            <dt>Total atoms</dt>
                            <dd>
                                {formatNumber(modelSummary.totalAtoms, 0)}
                                <span className="model-stat-sub">{formatNumber(modelSummary.referenceSites, 0)} sites</span>
                            </dd>
                        </div>
                    </dl>
                </section>
            )}

            {!modelSummary && structureError && (
                <div className="model-summary-empty">{structureError}</div>
            )}

            {renderRValuePanel()}

            <div className="plot-grid">
                {gridFiles.map((file) => renderPlotCard(file))}
            </div>

            {!loading && plotFiles.length === 0 && (
                <div className="empty-state">Open a run folder such as data to populate the dashboard.</div>
            )}
        </section>
    );
};

export default Dashboard;
