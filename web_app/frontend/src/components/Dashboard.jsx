import React, { useEffect, useMemo, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import { plotMetadataFromFile, readAndParseLocalPlotFile } from '../browserData';
import InteractivePlot from './InteractivePlot';
import ModelSummary from './ModelSummary';
import './Dashboard.css';

const plotOrder = ['r_value', 'bragg', 'xray_sq', 'neutron_sq', 'xpdf', 'npdf', 'pdf_partials', 'stog'];

const Dashboard = ({ directory, localRun }) => {
    const [files, setFiles] = useState([]);
    const [metadata, setMetadata] = useState({});
    const [structure, setStructure] = useState(null);
    const [structureError, setStructureError] = useState(null);
    const [error, setError] = useState(null);
    const [loading, setLoading] = useState(false);
    const [localStatus, setLocalStatus] = useState(null);
    const [showRValue, setShowRValue] = useState(false);

    useEffect(() => {
        if (localRun) {
            let cancelled = false;
            const loadedFiles = localRun.files || [];
            const plotFiles = loadedFiles.filter((file) => file.plotKind);
            const diagnostics = localRun.diagnostics;
            setFiles(loadedFiles);
            setMetadata(Object.fromEntries(
                plotFiles
                    .map((file) => [file.path, plotMetadataFromFile(file)])
            ));
            setStructure(null);
            setStructureError(localRun.structureFile ? 'Loading structure summary...' : localRun.structureError || 'No model structure detected');
            setError(null);
            setLoading(true);
            setLocalStatus(
                diagnostics
                    ? `Indexed ${diagnostics.supportedFileCount} supported files from ${diagnostics.selectedFileCount} selected files`
                    : 'Indexed local run files'
            );

            const parsePlots = async () => {
                if (!plotFiles.length) {
                    if (!cancelled) {
                        setLoading(false);
                        setLocalStatus('No plot files detected');
                    }
                    return;
                }
                setLocalStatus(`Parsing ${plotFiles.length} plot files...`);
                const parsedEntries = await Promise.all(plotFiles.map(async (file) => {
                    try {
                        return [file.path, { ...file, plotData: await readAndParseLocalPlotFile(file) }];
                    } catch (plotError) {
                        return [file.path, { ...file, parseError: plotError.message || 'Could not parse plot file' }];
                    }
                }));
                if (cancelled) return;
                const parsedByPath = Object.fromEntries(parsedEntries);
                const nextFiles = loadedFiles.map((file) => parsedByPath[file.path] || file);
                setFiles(nextFiles);
                setMetadata(Object.fromEntries(
                    nextFiles
                        .filter((file) => file.plotKind)
                        .map((file) => [file.path, plotMetadataFromFile(file)])
                ));
                setLoading(false);
                setLocalStatus(`Loaded ${plotFiles.length} plot files`);
            };

            parsePlots();

            let structureWorker = null;
            if (localRun.structureFile) {
                structureWorker = new Worker(new URL('../workers/localStructureWorker.js', import.meta.url), {
                    type: 'module'
                });
                structureWorker.onmessage = (event) => {
                    if (cancelled) return;
                    if (event.data.error) {
                        setStructure(null);
                        setStructureError(event.data.error);
                        return;
                    }
                    setStructure(event.data.result);
                    setStructureError(null);
                };
                structureWorker.onerror = () => {
                    if (cancelled) return;
                    setStructure(null);
                    setStructureError('Browser structure summary parser failed');
                };
                structureWorker.postMessage({
                    id: 1,
                    file: localRun.structureFile,
                    maxPoints: 100
                });
            }

            return () => {
                cancelled = true;
                structureWorker?.terminate();
            };
        }

        const fetchDashboard = async () => {
            setLoading(true);
            setError(null);
            setLocalStatus(null);
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
                setError(err.response?.data?.error || null);
                setStructure(null);
                setStructureError(null);
            } finally {
                setLoading(false);
            }
        };

        fetchDashboard();
    }, [directory, localRun]);

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

    const renderPlotBody = (file, variant) => {
        if (file.sourceFile && !file.plotData && !file.parseError) {
            return <div className="plot-loading">Parsing plot file...</div>;
        }
        if (file.sourceFile && file.parseError) {
            return null;
        }
        return <InteractivePlot file={file} plotData={file.plotData} variant={variant} />;
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
                {renderPlotBody(file)}
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
                {showRValue && renderPlotBody(rValueFile, 'wide')}
                {rValueFile.parseError && <div className="dashboard-error">{rValueFile.parseError}</div>}
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

            {localStatus && <div className="dashboard-local-status">{localStatus}</div>}

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

            <footer className="app-footer">
                &copy; 2026 Tsung-Han Yang. All rights reserved.
            </footer>
        </section>
    );
};

export default Dashboard;
