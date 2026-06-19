import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import { fileSignature, isStaticMode, plotMetadataFromFile, readAndParseLocalPlotFile, WATCH_INTERVAL_MS } from '../browserData';
import InteractivePlot from './InteractivePlot';
import ModelSummary from './ModelSummary';
import './Dashboard.css';

const plotOrder = ['r_value', 'bragg', 'xray_sq', 'neutron_sq', 'xpdf', 'npdf', 'pdf_partials', 'stog'];

const defaultHiddenPlotPaths = (items) => new Set(
    items
        .filter((file) => file.plotKind === 'stog')
        .map((file) => file.path)
);

const rValueLogParts = (name) => {
    const match = name.match(/^(.+)-(\d{2,})\.log$/);
    return match ? { stem: match[1].toLowerCase(), sequence: Number(match[2]) } : null;
};

const comparePlotFiles = (a, b) => {
    const kindOrder = plotOrder.indexOf(a.plotKind) - plotOrder.indexOf(b.plotKind);
    if (kindOrder !== 0) return kindOrder;

    const aLog = rValueLogParts(a.name);
    const bLog = rValueLogParts(b.name);
    if (aLog && bLog) {
        const stemOrder = aLog.stem.localeCompare(bLog.stem);
        if (stemOrder !== 0) return stemOrder;
        if (aLog.sequence !== bLog.sequence) return aLog.sequence - bLog.sequence;
    }

    return a.name.localeCompare(b.name, undefined, { numeric: true, sensitivity: 'base' });
};

const combineRValueFiles = (rValueFiles) => {
    if (!rValueFiles.length) return null;
    if (
        rValueFiles.length === 1
        || !rValueFiles.some((file) => file.sourceFile || file.plotData || file.parseError)
        || rValueFiles.some((file) => file.sourceFile && !file.plotData && !file.parseError)
    ) {
        return rValueFiles[0];
    }

    const parsedFiles = rValueFiles.filter((file) => file.plotData?.series?.[0]?.y?.length);
    if (!parsedFiles.length) {
        return {
            ...rValueFiles[0],
            parseError: rValueFiles.map((file) => file.parseError).filter(Boolean).join('; ') || 'Could not parse R-value logs'
        };
    }

    const yValues = parsedFiles.flatMap((file) => file.plotData.series[0].y);
    const lastParsed = parsedFiles[parsedFiles.length - 1];
    const parseErrors = rValueFiles
        .filter((file) => file.parseError)
        .map((file) => `${file.name}: ${file.parseError}`);

    return {
        ...parsedFiles[0],
        name: 'R-value',
        path: `r-value:${parsedFiles.map((file) => file.path).join('|')}`,
        sourceFile: undefined,
        parseError: parseErrors.join('; '),
        plotData: {
            kind: 'r_value',
            title: 'R-value',
            metrics: { final_chi_r: lastParsed.plotData.metrics?.final_chi_r },
            xLabel: 'Time steps',
            yLabel: 'log(χ)',
            series: [{
                label: 'R',
                x: yValues.map((_, index) => index),
                y: yValues
            }]
        }
    };
};

const Dashboard = ({ directory, localRun, watchFiles = false }) => {
    const [files, setFiles] = useState([]);
    const [metadata, setMetadata] = useState({});
    const [structure, setStructure] = useState(null);
    const [structureError, setStructureError] = useState(null);
    const [error, setError] = useState(null);
    const [loading, setLoading] = useState(false);
    const [localStatus, setLocalStatus] = useState(null);
    const [showRValue, setShowRValue] = useState(false);
    const [showLoadedFiles, setShowLoadedFiles] = useState(false);
    const [hiddenPlotPaths, setHiddenPlotPaths] = useState(() => new Set());
    const [dismissedErrors, setDismissedErrors] = useState(() => new Set());
    const signatureRef = useRef('');
    const pollInFlightRef = useRef(false);
    const manuallyToggledPathsRef = useRef(new Set());
    const currentRunIdRef = useRef(null);
    const filesRef = useRef([]);

    useEffect(() => {
        filesRef.current = files;
    }, [files]);

    const loadServerDashboard = useCallback(async ({ silent = false, loadedFiles: knownFiles = null } = {}) => {
        if (!silent) setLoading(true);
        setError(null);
        setLocalStatus(null);
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
            setHiddenPlotPaths((current) => {
                if (!silent) return defaultHiddenPlotPaths(loadedFiles);
                const next = new Set(current);
                loadedFiles
                    .filter((file) => file.plotKind === 'stog' && !manuallyToggledPathsRef.current.has(file.path))
                    .forEach((file) => next.add(file.path));
                return next;
            });
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
            if (!silent) setLoading(false);
        }
    }, [directory]);

    useEffect(() => {
        if (localRun) {
            let cancelled = false;
            const loadedFiles = localRun.files || [];
            const plotFiles = loadedFiles.filter((file) => file.plotKind);
            const diagnostics = localRun.diagnostics;
            // Live Data sends a new localRun (same runId) when files change. Refresh the existing
            // view in place instead of tearing it down, mirroring the Flask silent poll.
            const sameRun = localRun.runId != null && localRun.runId === currentRunIdRef.current;
            currentRunIdRef.current = localRun.runId ?? null;
            const prevByPath = new Map(filesRef.current.map((file) => [file.path, file]));
            signatureRef.current = fileSignature(loadedFiles);

            if (!sameRun) {
                setFiles(loadedFiles);
                setHiddenPlotPaths(defaultHiddenPlotPaths(loadedFiles));
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
            }

            const parsePlots = async () => {
                if (!plotFiles.length) {
                    if (!cancelled && !sameRun) {
                        setLoading(false);
                        setLocalStatus('No plot files detected');
                    }
                    return;
                }
                if (!sameRun) setLocalStatus(`Parsing ${plotFiles.length} plot files...`);
                const parsedEntries = await Promise.all(plotFiles.map(async (file) => {
                    const prev = prevByPath.get(file.path);
                    // Reuse already-parsed data for files that did not change between polls.
                    if (sameRun && prev?.plotData && prev.modified === file.modified && prev.size === file.size) {
                        return [file.path, { ...file, plotData: prev.plotData }];
                    }
                    try {
                        return [file.path, { ...file, plotData: await readAndParseLocalPlotFile(file) }];
                    } catch (plotError) {
                        return [file.path, { ...file, parseError: plotError.message || 'Could not parse plot file' }];
                    }
                }));
                if (cancelled) return;
                const parsedByPath = Object.fromEntries(parsedEntries);
                const nextFiles = loadedFiles.map((file) => parsedByPath[file.path] || file);
                signatureRef.current = fileSignature(nextFiles);
                setFiles(nextFiles);
                setMetadata(Object.fromEntries(
                    nextFiles
                        .filter((file) => file.plotKind)
                        .map((file) => [file.path, plotMetadataFromFile(file)])
                ));
                if (sameRun) {
                    // Hide newly-appeared STOG plots by default without discarding manual toggles.
                    setHiddenPlotPaths((current) => {
                        const next = new Set(current);
                        nextFiles
                            .filter((file) => file.plotKind === 'stog' && !manuallyToggledPathsRef.current.has(file.path))
                            .forEach((file) => next.add(file.path));
                        return next;
                    });
                } else {
                    setLoading(false);
                    setLocalStatus(null);
                }
            };

            parsePlots();

            // On a live update only re-parse the structure if the .rmc6f actually changed, and keep
            // the current model summary visible until the new one is ready (no flash to empty).
            const prevStructure = prevByPath.get(localRun.structureFile?.path);
            const structureChanged = !sameRun || (
                !prevStructure
                || prevStructure.modified !== localRun.structureFile?.modified
                || prevStructure.size !== localRun.structureFile?.size
            );
            let structureWorker = null;
            if (localRun.structureFile && structureChanged) {
                structureWorker = new Worker(new URL('../workers/localStructureWorker.js', import.meta.url), {
                    type: 'module'
                });
                structureWorker.onmessage = (event) => {
                    if (cancelled) return;
                    if (event.data.error) {
                        if (!sameRun) setStructure(null);
                        setStructureError(event.data.error);
                        return;
                    }
                    setStructure(event.data.result);
                    setStructureError(null);
                };
                structureWorker.onerror = () => {
                    if (cancelled) return;
                    if (!sameRun) setStructure(null);
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

        // Static mode has no Flask backend; the dashboard is driven entirely by localRun.
        if (!isStaticMode()) {
            loadServerDashboard();
        }
    }, [directory, loadServerDashboard, localRun]);

    useEffect(() => {
        // Reset view state only when the folder changes, not on each Live Data refresh
        // (directory stays constant across live updates of the same run).
        setShowLoadedFiles(false);
        manuallyToggledPathsRef.current = new Set();
        setDismissedErrors(new Set());
    }, [directory]);

    useEffect(() => {
        // Server-side file watching is Flask-only; static mode watches via App-level handle polling.
        if (!watchFiles || localRun || isStaticMode()) return undefined;

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
                    await loadServerDashboard({ silent: true, loadedFiles });
                }
            } catch (err) {
                setError(err.response?.data?.error || 'Failed to monitor dashboard files');
            } finally {
                pollInFlightRef.current = false;
            }
        };

        const interval = window.setInterval(pollForUpdates, WATCH_INTERVAL_MS);
        return () => window.clearInterval(interval);
    }, [directory, loadServerDashboard, localRun, watchFiles]);

    const allPlotFiles = useMemo(() => {
        return files
            .filter((file) => file.plotKind)
            .sort(comparePlotFiles);
    }, [files]);

    const plotFiles = useMemo(() => {
        return allPlotFiles.filter((file) => !hiddenPlotPaths.has(file.path));
    }, [allPlotFiles, hiddenPlotPaths]);

    const rValueFiles = useMemo(
        () => plotFiles.filter((file) => file.plotKind === 'r_value'),
        [plotFiles]
    );
    const rValueFile = useMemo(() => combineRValueFiles(rValueFiles), [rValueFiles]);
    const gridFiles = useMemo(
        () => plotFiles.filter((file) => file.plotKind !== 'r_value'),
        [plotFiles]
    );

    const handleTogglePlotVisibility = (path) => {
        manuallyToggledPathsRef.current.add(path);
        setHiddenPlotPaths((current) => {
            const next = new Set(current);
            if (next.has(path)) {
                next.delete(path);
            } else {
                next.add(path);
            }
            return next;
        });
    };

    const dismissError = (key) => {
        setDismissedErrors((current) => {
            const next = new Set(current);
            next.add(key);
            return next;
        });
    };

    const renderDashboardError = (key, message) => {
        if (!message || dismissedErrors.has(key)) return null;
        return (
            <div className="dashboard-error" role="alert">
                <span>{message}</span>
                <button
                    type="button"
                    className="notification-close"
                    onClick={() => dismissError(key)}
                    aria-label="Close notification"
                    title="Close"
                >
                    &times;
                </button>
            </div>
        );
    };

    const renderPlotBody = (file, variant) => {
        if (file.sourceFile && !file.plotData && !file.parseError) {
            return <div className="plot-loading">Parsing plot file...</div>;
        }
        if (file.sourceFile && file.parseError) {
            return null;
        }
        return (
            <InteractivePlot
                file={file}
                plotData={file.plotData}
                refreshKey={`${file.modified ?? ''}:${file.size ?? ''}`}
                variant={variant}
            />
        );
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
                {renderDashboardError(`plot:${file.path}:${file.parseError}`, file.parseError)}
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
                {renderDashboardError(`r-value:${rValueFile.parseError}`, rValueFile.parseError)}
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
                            Loaded {allPlotFiles.length} plot {allPlotFiles.length === 1 ? 'file' : 'files'}
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
                        {allPlotFiles.map((file) => {
                            const isHidden = hiddenPlotPaths.has(file.path);
                            const kindClass = `kind-${file.plotKind}`;
                            return (
                                <li key={file.path}>
                                    <span className={`loaded-file-badge ${kindClass}${isHidden ? ' is-hidden' : ''}`}>
                                        <span className="loaded-file-kind">{file.plotKind}</span>
                                        <span className="loaded-file-name">{file.name}</span>
                                        <button
                                            type="button"
                                            className="loaded-file-hide"
                                            onClick={() => handleTogglePlotVisibility(file.path)}
                                            aria-label={isHidden ? `Show ${file.name} chart` : `Hide ${file.name} chart`}
                                            aria-pressed={isHidden}
                                            title={isHidden ? 'Show chart' : 'Hide chart'}
                                        >
                                            &times;
                                        </button>
                                    </span>
                                </li>
                            );
                        })}
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

            {renderDashboardError(`dashboard:${error}`, error)}

            {localStatus && <div className="dashboard-local-status">{localStatus}</div>}

            <ModelSummary structure={structure} />

            {renderLoadedFilesPanel()}

            {renderRValuePanel()}

            <div className="plot-grid">
                {gridFiles.map((file) => renderPlotCard(file))}
            </div>

            {!loading && allPlotFiles.length === 0 && (
                <div className="empty-state">Open a run folder to populate the dashboard.</div>
            )}

            <footer className="app-footer">
                &copy; 2026 Tsung-Han Yang. All rights reserved.
            </footer>
        </section>
    );
};

export default Dashboard;
