// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useCallback, useContext, useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import {
    combineRValueFiles,
    fileSignature,
    isStaticMode,
    parseRunSettings,
    plotMetadataFromFile,
    readAndParseLocalPlotFile,
    WATCH_INTERVAL_MS
} from '../browserData';
import { saveSvgFiguresAsZip } from '../figureExport';
import { WatchdogBadge } from '../llm';
import { describeSymmetry, toleranceLadder } from '../symmetryModel';
import { isIncompleteStructure } from '../structureReport';
import { SymTolContext } from '../symTolContext';
import InteractivePlot from './InteractivePlot';
import {
    Banner, Card, CardTitle, Chip, EmptyState, IconButton, Page, PageIssues, Pill,
    useIssue, useIssueSet, useReportIssue
} from '../ui';
import InfoBadge from '../ui/InfoBadge';
import SaveMenu from '../ui/SaveMenu';
import ModelSummary from './ModelSummary';
import AppFooter from './AppFooter';
import './Dashboard.css';

const plotOrder = ['r_value', 'bragg', 'xray_sq', 'neutron_sq', 'exafs_q', 'exafs_r', 'xpdf', 'npdf', 'pdf_partials'];
const isDashboardPlotFile = (file) => file.plotKind && file.plotKind !== 'stog';

const CHART_SAVE_OPTIONS = [
    { id: 'png', label: 'PNG image', hint: '.png' },
    { id: 'svg', label: 'SVG vector', hint: '.svg' },
];

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

const STRUCTURE_LOADING = 'Loading structure…';

// Structure messages that only say "there is none (yet)": with no plots the
// single empty-state line covers them. Anything else is a real failure.
const isMissingStructure = (message) => message === STRUCTURE_LOADING
    || message === 'No model structure detected'
    || /^No \.rmc6f file found/.test(message);

// The metric is present but null when Rwp is undefined for the data (an observed
// column with no finite values, or one that is entirely zero). Show a dash there:
// a number in that slot reads as a fit quality, and 0.000 reads as a perfect one.
const renderRwpChip = (meta) => {
    const value = meta?.metrics?.rwp;
    if (value === undefined) return null;
    return (
        <Chip tone="success" strong>
            Rwp {Number.isFinite(value) ? Number(value).toPrecision(4) : '—'}
        </Chip>
    );
};

const Dashboard = ({ directory, localRun, watchFiles = false, wantAssistantData = false, onRunContextChange }) => {
    const [files, setFiles] = useState([]);
    const [metadata, setMetadata] = useState({});
    const [structure, setStructure] = useState(null);
    const [structureError, setStructureError] = useState(null);
    // Set while a Live Data re-read of the structure came back incomplete and
    // the previous complete summary is being kept on screen (ModelSummary's
    // "previous read" chip).
    const [structureStale, setStructureStale] = useState(false);
    const structureRef = useRef(null);
    // Parsed <stem>.dat run-control settings (static mode) for the AI assistant.
    const [runSettings, setRunSettings] = useState(null);
    const settingsSigRef = useRef('');
    const [error, setError] = useState(null);
    const [loading, setLoading] = useState(false);
    const [showRValue, setShowRValue] = useState(false);
    const [showLoadedFiles, setShowLoadedFiles] = useState(false);
    const [hiddenPlotPaths, setHiddenPlotPaths] = useState(() => new Set());
    const [dismissedErrors, setDismissedErrors] = useState(() => new Set());
    const [savingAll, setSavingAll] = useState(false);
    // A failed "Save all figures": one line under the Loaded-files header,
    // cleared by the next save (the save menu does not await the handler).
    const [saveAllError, setSaveAllError] = useState(null);
    const pageRef = useRef(null);
    const signatureRef = useRef('');
    const pollInFlightRef = useRef(false);
    const manuallyToggledPathsRef = useRef(new Set());
    const currentRunIdRef = useRef(null);
    const filesRef = useRef([]);

    useEffect(() => {
        filesRef.current = files;
    }, [files]);

    useEffect(() => {
        structureRef.current = structure;
    }, [structure]);

    const loadServerDashboard = useCallback(async ({ silent = false, loadedFiles: knownFiles = null } = {}) => {
        if (!silent) setLoading(true);
        setError(null);
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
                if (!silent) return new Set();
                return current;
            });
            const plotFiles = loadedFiles.filter(isDashboardPlotFile);
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
                // A silent (Live Data) re-read that lands mid-write keeps the
                // previous complete model instead of flashing a short count.
                if (silent && structureRef.current && isIncompleteStructure(structureResponse.data)) {
                    setStructureStale(true);
                } else {
                    setStructure(structureResponse.data);
                    setStructureStale(false);
                }
                setStructureError(null);
            } catch (structureErr) {
                setStructure(null);
                setStructureStale(false);
                setStructureError(structureErr.response?.data?.error || 'No model structure detected');
            }
        } catch (err) {
            setError(err.response?.data?.error || 'Could not list the run folder');
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
            const plotFiles = loadedFiles.filter(isDashboardPlotFile);
            // Live Data sends a new localRun (same runId) when files change. Refresh the existing
            // view in place instead of tearing it down, mirroring the Flask silent poll.
            const sameRun = localRun.runId != null && localRun.runId === currentRunIdRef.current;
            currentRunIdRef.current = localRun.runId ?? null;
            const prevByPath = new Map(filesRef.current.map((file) => [file.path, file]));
            signatureRef.current = fileSignature(loadedFiles);

            if (!sameRun) {
                setFiles(loadedFiles);
                setHiddenPlotPaths(new Set());
                setMetadata(Object.fromEntries(
                    plotFiles
                        .map((file) => [file.path, plotMetadataFromFile(file)])
                ));
                setStructure(null);
                setStructureStale(false);
                setStructureError(localRun.structureFile ? STRUCTURE_LOADING : localRun.structureError || 'No model structure detected');
                setError(null);
                setLoading(true);
            }

            const parsePlots = async () => {
                if (!plotFiles.length) {
                    if (!cancelled && !sameRun) setLoading(false);
                    return;
                }
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
                        .filter(isDashboardPlotFile)
                        .map((file) => [file.path, plotMetadataFromFile(file)])
                ));
                if (!sameRun) setLoading(false);
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
                    // Live Data poll that caught the .rmc6f mid-write: keep the
                    // previous complete summary rather than a short composition.
                    if (sameRun && structureRef.current && isIncompleteStructure(event.data.result)) {
                        setStructureStale(true);
                        return;
                    }
                    setStructure(event.data.result);
                    setStructureStale(false);
                    setStructureError(null);
                };
                structureWorker.onerror = () => {
                    if (cancelled) return;
                    if (!sameRun) setStructure(null);
                    setStructureError('Structure parser failed');
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
            return undefined;
        }

        // No run loaded in static mode (e.g. the demo was toggled off): clear everything.
        currentRunIdRef.current = null;
        signatureRef.current = '';
        setFiles([]);
        setMetadata({});
        setStructure(null);
        setStructureStale(false);
        setStructureError(null);
        setError(null);
        setSaveAllError(null);
        setLoading(false);
        setHiddenPlotPaths(new Set());
        return undefined;
    }, [directory, loadServerDashboard, localRun]);

    // Parse the run-control settings (<stem>.dat) only once the assistant is in
    // use — the file is tiny but this keeps Dashboard/KDE-only startup clean and
    // re-reads only when the file identity/mtime changes.
    useEffect(() => {
        const settingsFile = wantAssistantData ? localRun?.settingsFile : null;
        const sig = settingsFile ? `${settingsFile.path}:${settingsFile.modified}` : '';
        if (sig === settingsSigRef.current) return undefined;
        settingsSigRef.current = sig;
        if (!settingsFile?.sourceFile) {
            setRunSettings(null);
            return undefined;
        }
        let cancelled = false;
        settingsFile.sourceFile.text()
            .then((text) => { if (!cancelled) setRunSettings(parseRunSettings(text)); })
            .catch(() => { if (!cancelled) setRunSettings(null); });
        return () => { cancelled = true; };
    }, [wantAssistantData, localRun]);

    useEffect(() => {
        // Reset view state only when the folder changes, not on each Live Data refresh
        // (directory stays constant across live updates of the same run).
        setShowLoadedFiles(false);
        manuallyToggledPathsRef.current = new Set();
        setDismissedErrors(new Set());
        // A failed "Save all figures" belongs to the previous run's figures.
        setSaveAllError(null);
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
            .filter(isDashboardPlotFile)
            .sort(comparePlotFiles);
    }, [files]);

    const plotFiles = useMemo(() => {
        return allPlotFiles.filter((file) => !hiddenPlotPaths.has(file.path));
    }, [allPlotFiles, hiddenPlotPaths]);

    const rValueFiles = useMemo(
        () => plotFiles.filter((file) => file.plotKind === 'r_value'),
        [plotFiles]
    );
    // The run the Model card describes picks which log group is charted.
    const structurePath = localRun ? localRun.structureFile?.path : structure?.source;
    const rValueFile = useMemo(
        () => combineRValueFiles(rValueFiles, structurePath),
        [rValueFiles, structurePath]
    );
    const gridFiles = useMemo(
        () => plotFiles.filter((file) => file.plotKind !== 'r_value'),
        [plotFiles]
    );

    // Every problem on this page is listed in its Problems section: the run
    // folder, a structure that failed to load, and each plot file or χ² log
    // that failed to parse (their cards keep their own line too).
    useIssue('run-folder', error, { source: 'Run folder' });
    useIssue('structure', structureError && !isMissingStructure(structureError) ? structureError : null, { source: 'Structure' });
    const fileIssues = useMemo(() => {
        const items = allPlotFiles
            .filter((file) => file.parseError && file.plotKind !== 'r_value')
            .map((file) => ({ key: file.path, message: file.parseError, source: file.name }));
        if (rValueFile?.parseErrors?.length) {
            rValueFile.parseErrors.forEach(({ name, message }) => items.push({ key: `log:${name}`, message, source: name }));
        } else if (rValueFile?.parseError) {
            items.push({ key: `log:${rValueFile.path}`, message: rValueFile.parseError, source: rValueFile.name });
        }
        return items;
    }, [allPlotFiles, rValueFile]);
    useIssueSet('files', fileIssues);
    const reportIssue = useReportIssue();

    // Detected space group for the AI assistant's run context, at the shared
    // tolerance — keeps symmetryModel out of the llm module's imports. The
    // ladder rides along so the context can express distortion magnitude.
    // Gated on wantAssistantData: this symmetry finder (+ ladder) is redundant
    // with ModelSummary's and only feeds the assistant, so it stays idle until
    // the user opens the AI Assistant page — Dashboard/KDE-only startup does no
    // extra work.
    const sharedSymTol = useContext(SymTolContext);
    const symTol = sharedSymTol ? sharedSymTol[0] : 0.2;
    const symmetry = useMemo(() => {
        if (!wantAssistantData) return null;
        const detected = describeSymmetry(structure, symTol);
        if (!detected) return null;
        return { ...detected, toleranceA: symTol, ladder: toleranceLadder(structure, 1.0) };
    }, [wantAssistantData, structure, symTol]);

    // Publish the parsed run context for App → AssistantPage, but only once the
    // assistant has been opened; the WatchdogBadge below stays on the dashboard.
    const assistantRun = useMemo(() => (wantAssistantData ? {
        runName: localRun ? localRun.name : directory,
        plotFiles: allPlotFiles,
        rValueFile,
        structure,
        symmetry,
        runSettings,
    } : null), [wantAssistantData, localRun, directory, allPlotFiles, rValueFile, structure, symmetry, runSettings]);

    useEffect(() => {
        onRunContextChange?.(assistantRun);
    }, [assistantRun, onRunContextChange]);

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

    // Rasterize every chart currently rendered in the dashboard to its own PNG.
    // We read the live SVG nodes (rather than holding refs to each plot) so the
    // R-value strip is included only when expanded, matching what the user sees.
    const handleSaveAllFigures = async (format) => {
        const root = pageRef.current;
        if (!root || savingAll) return;
        setSavingAll(true);
        setSaveAllError(null);
        try {
            const figures = [];
            let index = 0;
            root.querySelectorAll('[data-figure-card]').forEach((card) => {
                const svg = card.querySelector('.interactive-plot svg');
                if (!svg) return;
                index += 1;
                const title = card.querySelector('.ui-card__title')?.textContent?.trim();
                figures.push({ svgElement: svg, name: title || `figure-${index}` });
            });
            // Nothing drawn yet is not a successful save.
            if (!figures.length) throw new Error('No chart has been drawn yet');
            await saveSvgFiguresAsZip(figures, format, `figures-${format}.zip`);
        } catch (failure) {
            const message = failure?.message || 'Could not save the figures';
            // In a page, SaveMenu lists it in the Problems section.
            if (reportIssue) throw new Error(message);
            setSaveAllError(message);
        } finally {
            setSavingAll(false);
        }
    };

    const dismissError = (key) => {
        setDismissedErrors((current) => {
            const next = new Set(current);
            next.add(key);
            return next;
        });
    };

    // One line; `details` (optional) goes behind an "Error details" ? help.
    const renderDashboardError = (key, message, details = null) => {
        if (!message || dismissedErrors.has(key)) return null;
        return (
            <Banner tone="danger" flush role="alert" onDismiss={() => dismissError(key)}>
                {message}
                {details && (
                    <>
                        {' '}
                        <InfoBadge label="Error details" align="end">{details}</InfoBadge>
                    </>
                )}
            </Banner>
        );
    };

    const renderPlotBody = (file, variant) => {
        if (file.sourceFile && !file.plotData && !file.parseError) {
            return <div className="ui-loading ui-loading--sm">Parsing plot file…</div>;
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
        const title = meta?.title || file.name;
        return (
            <Card as="article" clip lift data-figure-card="" key={file.path}>
                <div className="ui-card__header-flush">
                    <div className="ui-card__heading">
                        <CardTitle>{title}</CardTitle>
                        {title !== file.name && (
                            <span className="ui-card__source" title={file.path}>{file.name}</span>
                        )}
                    </div>
                    {renderRwpChip(meta)}
                </div>
                {renderPlotBody(file)}
                {renderDashboardError(`plot:${file.path}:${file.parseError}`, file.parseError)}
            </Card>
        );
    };

    const renderRValuePanel = () => {
        if (!rValueFile) return null;
        const meta = metadata[rValueFile.path];
        const title = meta?.title || rValueFile.name;
        const sourceNames = rValueFile.sourceNames?.length ? rValueFile.sourceNames : [rValueFile.name];
        const sourceLabel = sourceNames.join(', ');
        const otherRuns = rValueFile.otherRuns || [];
        // Several failed logs: one counted line, the per-log list behind "?".
        const failedLogs = rValueFile.parseErrors || [];
        const errorLine = failedLogs.length > 1
            ? `${failedLogs.length} chi² logs could not be parsed`
            : rValueFile.parseError;
        const errorDetails = failedLogs.length > 1
            ? failedLogs.map(({ name, message }) => <p key={name}>{name}: {message}</p>)
            : null;
        const showErrorDetails = Boolean(errorDetails && errorLine
            && !dismissedErrors.has(`r-value:${rValueFile.parseError}`));
        return (
            <Card
                as="article"
                // A clipped card would hide the error-details popover; round
                // the ends instead while it is on screen.
                clip={!showErrorDetails}
                roundEnds={showErrorDetails}
                lift
                data-figure-card=""
                className={`r-value-card${showRValue ? '' : ' is-collapsed'}`}
            >
                <div className="ui-card__header-flush ui-card__header-flush--padded">
                    <div className="ui-card__heading">
                        <CardTitle>{title}</CardTitle>
                        {sourceLabel !== title && (
                            // One unit per log name, so a long list wraps between
                            // names (not at the hyphen inside one).
                            <span className="ui-card__source" title={sourceLabel}>
                                {sourceNames.map((name, index) => (
                                    <React.Fragment key={index}>
                                        {index > 0 && ', '}
                                        <span className="ui-card__source-item">{name}</span>
                                    </React.Fragment>
                                ))}
                            </span>
                        )}
                        {otherRuns.length > 0 && (
                            <Chip title={`Other runs in this folder, not charted: ${otherRuns.join(', ')}`}>
                                +{otherRuns.length} {otherRuns.length === 1 ? 'run' : 'runs'}
                            </Chip>
                        )}
                    </div>
                    <div className="ui-card__header-actions">
                        <WatchdogBadge rValueFile={rValueFile} />
                        {renderRwpChip(meta)}
                        <Pill
                            onClick={() => setShowRValue((value) => !value)}
                            aria-expanded={showRValue}
                        >
                            {showRValue ? 'Hide' : 'Show'}
                        </Pill>
                    </div>
                </div>
                {showRValue && renderPlotBody(rValueFile, 'wide')}
                {renderDashboardError(`r-value:${rValueFile.parseError}`, errorLine, errorDetails)}
            </Card>
        );
    };

    const renderLoadedFilesPanel = () => {
        if (allPlotFiles.length === 0) {
            // "Not found" / loading are covered by the empty-state line; a real
            // structure failure is listed in the Problems section.
            return null;
        }

        return (
            <Card as="article" lift data-figure-card="" className={`loaded-files-card${showLoadedFiles ? '' : ' is-collapsed'}`}>
                <div className="ui-card__header-flush ui-card__header-flush--padded">
                    <div>
                        <CardTitle>
                            Loaded {allPlotFiles.length} plot {allPlotFiles.length === 1 ? 'file' : 'files'}
                        </CardTitle>
                        {structureError && <p className="ui-card__subtitle">{structureError}</p>}
                    </div>
                    <div className="ui-card__header-actions">
                        {hasFigures && (
                            <SaveMenu
                                onSave={handleSaveAllFigures}
                                options={CHART_SAVE_OPTIONS}
                                label="Save all figures"
                                name="All figures"
                                align="right"
                                busy={savingAll}
                                className="ui-save--accent"
                            />
                        )}
                        <Pill
                            onClick={() => setShowLoadedFiles((value) => !value)}
                            aria-expanded={showLoadedFiles}
                        >
                            {showLoadedFiles ? 'Hide' : 'Show'}
                        </Pill>
                    </div>
                </div>
                {saveAllError && (
                    <Banner tone="danger" sm role="alert" className="save-all-error" onDismiss={() => setSaveAllError(null)}>
                        {saveAllError}
                    </Banner>
                )}
                {showLoadedFiles && (
                    <ul className="ui-card__section loaded-files-list">
                        {allPlotFiles.map((file) => {
                            const isHidden = hiddenPlotPaths.has(file.path);
                            const kindClass = `kind-${file.plotKind}`;
                            return (
                                <li key={file.path}>
                                    <span className={`ui-file-chip ${kindClass}${isHidden ? ' is-hidden' : ''}`}>
                                        <span className="ui-file-chip__kind">{file.plotKind}</span>
                                        <span className="ui-file-chip__name">{file.name}</span>
                                        <IconButton
                                            variant="remove"
                                            onClick={() => handleTogglePlotVisibility(file.path)}
                                            aria-label={isHidden ? `Show ${file.name} chart` : `Hide ${file.name} chart`}
                                            aria-pressed={isHidden}
                                            title={isHidden ? 'Show chart' : 'Hide chart'}
                                        >
                                            &times;
                                        </IconButton>
                                    </span>
                                </li>
                            );
                        })}
                    </ul>
                )}
            </Card>
        );
    };

    const hasFigures = gridFiles.length > 0 || (rValueFile && showRValue);
    const hasRun = Boolean(localRun) || files.length > 0;

    return (
        <Page wide ref={pageRef}>
            <PageIssues />

            <ModelSummary structure={structure} stale={structureStale} />

            {renderLoadedFilesPanel()}

            {renderRValuePanel()}

            <div className="plot-grid">
                {gridFiles.map((file) => renderPlotCard(file))}
            </div>

            {loading && allPlotFiles.length === 0 && <EmptyState>Loading…</EmptyState>}

            {!loading && allPlotFiles.length === 0 && (
                <EmptyState>{hasRun ? 'No plot files in this folder.' : 'Open a run folder.'}</EmptyState>
            )}

            <AppFooter />
        </Page>
    );
};

export default Dashboard;
