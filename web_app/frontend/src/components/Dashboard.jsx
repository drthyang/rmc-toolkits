import React, { useEffect, useMemo, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import InteractivePlot from './InteractivePlot';
import ModelSummary from './ModelSummary';
import './Dashboard.css';

const plotOrder = ['r_value', 'bragg', 'xray_sq', 'neutron_sq', 'xpdf', 'npdf', 'pdf_partials', 'stog'];

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

            <ModelSummary structure={structure} />

            {!structure && structureError && (
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
