import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import { plotMetadataFromFile } from '../browserData';
import InteractivePlot from './InteractivePlot';
import ModelSummary from './ModelSummary';
import './Dashboard.css';

const plotOrder = ['r_value', 'bragg', 'xray_sq', 'neutron_sq', 'xpdf', 'npdf', 'pdf_partials', 'stog'];
const WATCH_INTERVAL_MS = 3000;

const fileSignature = (items) => items
    .map((file) => `${file.path}:${file.modified ?? ''}:${file.size ?? ''}:${file.plotKind ?? ''}`)
    .sort()
    .join('|');

const Dashboard = ({ directory, localRun, watchFiles = false }) => {
    const [files, setFiles] = useState([]);
    const [metadata, setMetadata] = useState({});
    const [structure, setStructure] = useState(null);
    const [structureError, setStructureError] = useState(null);
    const [error, setError] = useState(null);
    const [loading, setLoading] = useState(false);
    const [showRValue, setShowRValue] = useState(false);
    const [showLoadedFiles, setShowLoadedFiles] = useState(false);
    const [hiddenPlotPaths, setHiddenPlotPaths] = useState(() => new Set());
    const signatureRef = useRef('');
    const pollInFlightRef = useRef(false);

    const loadDashboard = useCallback(async ({ silent = false, loadedFiles: knownFiles = null } = {}) => {
        if (silent) {
            setError(null);
        } else {
            setLoading(true);
            setError(null);
        }

        try {
            let loadedFiles = knownFiles;
            if (!loadedFiles) {
                const response = await axios.get(`${API_BASE_URL}/api/files`, {
                    params: { dir: directory || '.' }
                });
                loadedFiles = response.data.files || [];
            }

            signatureRef.current = fileSignature(loadedFiles);
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
            if (!silent) setLoading(false);
        }
    }, [directory]);

    useEffect(() => {
        if (localRun) {
            const loadedFiles = localRun.files || [];
            signatureRef.current = fileSignature(loadedFiles);
            setFiles(loadedFiles);
            setMetadata(Object.fromEntries(
                loadedFiles
                    .filter((file) => file.plotKind)
                    .map((file) => [file.path, plotMetadataFromFile(file)])
            ));
            setStructure(localRun.structure || null);
            setStructureError(localRun.structure ? null : localRun.structureError || 'No model structure detected');
            setError(null);
            setLoading(false);
            return;
        }

        loadDashboard();
    }, [directory, loadDashboard, localRun]);

    useEffect(() => {
        setShowLoadedFiles(false);
        setHiddenPlotPaths(new Set());
    }, [directory, localRun]);

    useEffect(() => {
        if (!watchFiles || localRun) return undefined;

        const pollForUpdates = async () => {
            if (pollInFlightRef.current) return;
            pollInFlightRef.current = true;
            try {
                const response = await axios.get(`${API_BASE_URL}/api/files`, {
                    params: { dir: directory || '.' }
                });
                const loadedFiles = response.data.files || [];
                const nextSignature = fileSignature(loadedFiles);
                if (nextSignature !== signatureRef.current) {
                    await loadDashboard({ silent: true, loadedFiles });
                }
            } catch (err) {
                setError(err.response?.data?.error || 'Failed to monitor dashboard files');
            } finally {
                pollInFlightRef.current = false;
            }
        };

        const interval = window.setInterval(pollForUpdates, WATCH_INTERVAL_MS);
        return () => window.clearInterval(interval);
    }, [directory, loadDashboard, localRun, watchFiles]);

    const allPlotFiles = useMemo(() => {
        return files
            .filter((file) => file.plotKind)
            .sort((a, b) => plotOrder.indexOf(a.plotKind) - plotOrder.indexOf(b.plotKind));
    }, [files]);

    const plotFiles = useMemo(() => {
        return allPlotFiles.filter((file) => !hiddenPlotPaths.has(file.path));
    }, [allPlotFiles, hiddenPlotPaths]);

    const rValueFile = useMemo(
        () => plotFiles.find((file) => file.plotKind === 'r_value') || null,
        [plotFiles]
    );
    const gridFiles = useMemo(
        () => plotFiles.filter((file) => file.plotKind !== 'r_value'),
        [plotFiles]
    );

    const handleHidePlot = (path) => {
        setHiddenPlotPaths((current) => {
            const next = new Set(current);
            next.add(path);
            return next;
        });
    };

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
                <InteractivePlot
                    file={file}
                    plotData={file.plotData}
                    refreshKey={`${file.modified ?? ''}:${file.size ?? ''}`}
                />
                {file.parseError && <div className="dashboard-error">{file.parseError}</div>}
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
                {showRValue && (
                    <InteractivePlot
                        file={rValueFile}
                        plotData={rValueFile.plotData}
                        refreshKey={`${rValueFile.modified ?? ''}:${rValueFile.size ?? ''}`}
                        variant="wide"
                    />
                )}
            </article>
        );
    };

    const renderLoadedFilesPanel = () => {
        if (allPlotFiles.length === 0) {
            return structureError ? <div className="model-summary-empty">{structureError}</div> : null;
        }

        return (
            <article className={`plot-card loaded-files-card${showLoadedFiles ? '' : ' is-collapsed'}`}>
                <div className="plot-card-header">
                    <div>
                        <h3>
                            Loaded {plotFiles.length} plot {plotFiles.length === 1 ? 'file' : 'files'}
                        </h3>
                        {structureError && <p>{structureError}</p>}
                    </div>
                    <button
                        type="button"
                        className="panel-toggle"
                        onClick={() => setShowLoadedFiles((value) => !value)}
                        aria-expanded={showLoadedFiles}
                    >
                        {showLoadedFiles ? 'Hide' : 'Show'}
                    </button>
                </div>
                {showLoadedFiles && (
                    <ul className="loaded-files-list">
                        {plotFiles.map((file) => (
                            <li key={file.path}>
                                <span className="loaded-file-badge">
                                    <span className="loaded-file-kind">{file.plotKind}</span>
                                    <span className="loaded-file-name">{file.name}</span>
                                    <button
                                        type="button"
                                        className="loaded-file-hide"
                                        onClick={() => handleHidePlot(file.path)}
                                        aria-label={`Hide ${file.name} chart`}
                                        title="Hide chart"
                                    >
                                        &times;
                                    </button>
                                </span>
                            </li>
                        ))}
                    </ul>
                )}
            </article>
        );
    };

    return (
        <section className="dashboard-page">
            <div className="dashboard-toolbar">
                <div>
                    <h2>Run Dashboard</h2>
                    <p>{localRun ? localRun.name : directory}</p>
                </div>
                {loading && <span className="status-pill">Loading</span>}
            </div>

            {error && <div className="dashboard-error">{error}</div>}

            <ModelSummary structure={structure} />

            {renderLoadedFilesPanel()}

            {renderRValuePanel()}

            <div className="plot-grid">
                {gridFiles.map((file) => renderPlotCard(file))}
            </div>

            {!loading && allPlotFiles.length === 0 && (
                <div className="empty-state">Open a run folder such as data to populate the dashboard.</div>
            )}
        </section>
    );
};

export default Dashboard;
