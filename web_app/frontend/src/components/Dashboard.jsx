// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import {
    chartedRValueFile,
    DASHBOARD_MAX_POINTS,
    dashboardPlotFiles,
    fileSignature,
    isDashboardPlotFile,
    isStaticMode,
    parseDashboardPlotFile,
    plotDataVersion,
    plotMetadataFromFile,
    rValueGroupKey,
    WATCH_INTERVAL_MS
} from '../browserData';
import { saveSvgFiguresAsZip } from '../figureExport';
import { attachFtUncertainty } from '../ftUncertaintyRun';
import { WatchdogBadge } from '../llm';
import { ribbonIssues, ribbonProp, ribbonSources } from '../grRibbon';
import { fetchPlotData } from '../plotDataCache';
import { RibbonNote, RibbonToggle } from './GrRibbon';
import { isIncompleteStructure } from '../structureReport';
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

const CHART_SAVE_OPTIONS = [
    { id: 'png', label: 'PNG image', hint: '.png' },
    { id: 'svg', label: 'SVG vector', hint: '.svg' },
];

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

// `cardMarks` ({ [cardId]: [{ x, label }] }): the x positions the AI Copilot
// marks on a plot card while it points at it (highlight_cards).
const NO_MARKS = {};

const NO_TARGETS = Object.freeze([]);
const NO_SERVER_PLOTS = Object.freeze({});
const NO_CHOICES = Object.freeze({});

// A plot file's version: what its chart's data is fetched for (InteractivePlot's refreshKey).
const plotVersion = (file) => `${file.modified ?? ''}:${file.size ?? ''}`;

// `wantAssistantData` (the AI Copilot has been opened) gates `onRunChange`:
// only then does the page publish its run ({ files, structure, structurePath,
// hiddenPlots }) for App's Copilot run module (useCopilotRun.js), so before
// that its parses and Live Data polls do not re-render App. In Flask mode it
// also fetches each plot file's data for that run (below), as a picked
// folder's files are parsed. `dataEpoch` is App's Flask-mode configuration
// epoch: the structure follows it (below).
const Dashboard = ({
    directory, localRun, watchFiles = false, wantAssistantData = false, onRunChange, cardMarks = NO_MARKS, dataEpoch = 0
}) => {
    const [files, setFiles] = useState([]);
    const [metadata, setMetadata] = useState({});
    const [structure, setStructure] = useState(null);
    const [structureError, setStructureError] = useState(null);
    // Set while a Live Data re-read of the structure came back incomplete and
    // the previous complete summary is being kept on screen (ModelSummary's
    // "previous read" chip).
    const [structureStale, setStructureStale] = useState(false);
    const structureRef = useRef(null);
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
    // Every Flask structure request is numbered and only the latest one's
    // answer is applied: a Live Data poll and a new configuration's re-read in
    // flight together, or a request for the folder or run that was open
    // before, never overwrite a newer structure.
    const structureRequestRef = useRef(0);
    // Flask mode, once the Copilot is open: the data of each plot file the
    // published run carries, `{ [path]: { plotData } | { parseError } }`, and
    // the version last asked for per path. A load of a folder starts a new
    // map, so an answer for the folder before, or for a version older than
    // the one last asked for, is dropped.
    const [serverPlots, setServerPlots] = useState(NO_SERVER_PLOTS);
    const serverPlotRequestsRef = useRef(new Map());
    // An xPDF plot's σ choice where its file carries its own σ beside the
    // propagated one ({ [path]: 'file' | 'propagated' }; the file's by default).
    const [ribbonChoices, setRibbonChoices] = useState(NO_CHOICES);
    // Flask mode: the plot data each chart fetched ({ [path]: plotData }), for
    // the ribbon's caption and the Problems section (static mode's files carry it).
    const [chartData, setChartData] = useState(NO_SERVER_PLOTS);

    useEffect(() => {
        filesRef.current = files;
    }, [files]);

    useEffect(() => {
        structureRef.current = structure;
    }, [structure]);

    // The run's structure from Flask: the Model information card, the Detected
    // SG card and its CIF (from the response's average-structure basis, as the
    // browser's structureFromRmc6f gives a picked folder) and, published below,
    // the Copilot's structure, symmetry and symmetry tools. A silent (Live
    // Data) re-read that lands mid-write keeps the previous complete model
    // instead of flashing a short count.
    const loadServerStructure = useCallback(async ({ silent = false } = {}) => {
        structureRequestRef.current += 1;
        const request = structureRequestRef.current;
        try {
            const structureResponse = await axios.get(`${API_BASE_URL}/api/structure`, {
                params: { dir: directory || '.', maxPoints: DASHBOARD_MAX_POINTS }
            });
            if (request !== structureRequestRef.current) return;
            if (silent && structureRef.current && isIncompleteStructure(structureResponse.data)) {
                setStructureStale(true);
            } else {
                setStructure(structureResponse.data);
                setStructureStale(false);
            }
            setStructureError(null);
        } catch (structureErr) {
            if (request !== structureRequestRef.current) return;
            setStructure(null);
            setStructureStale(false);
            setStructureError(structureErr.response?.data?.error || 'No model structure detected');
        }
    }, [directory]);

    const loadServerDashboard = useCallback(async ({ silent = false, loadedFiles: knownFiles = null } = {}) => {
        // A (non-silent) load of a folder: a structure or plot data still in
        // flight is the folder's before.
        if (!silent) {
            structureRequestRef.current += 1;
            serverPlotRequestsRef.current = new Map();
            setServerPlots(NO_SERVER_PLOTS);
            setLoading(true);
        }
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

            await loadServerStructure({ silent });
        } catch (err) {
            setError(err.response?.data?.error || 'Could not list the run folder');
            structureRequestRef.current += 1;
            setStructure(null);
            setStructureError(null);
        } finally {
            if (!silent) setLoading(false);
        }
    }, [directory, loadServerStructure]);

    // Flask mode: App's configuration epoch moves when the run's .rmc6f changes
    // on disk (seen by a Load of the same folder or a Live Data poll), and the
    // structure is re-read in place, as the analysis pages re-read theirs. So
    // the model card, the Detected SG card and the basis the Copilot's
    // symmetry tools read are the configuration the analysis pages and the
    // Copilot's Flask analyses read. The plots keep following Live Data.
    const epochRef = useRef(dataEpoch);
    useEffect(() => {
        if (epochRef.current === dataEpoch) return;
        epochRef.current = dataEpoch;
        if (localRun || isStaticMode()) return;
        loadServerStructure({ silent: true });
    }, [dataEpoch, localRun, loadServerStructure]);

    useEffect(() => {
        if (localRun) {
            let cancelled = false;
            // A server structure still in flight belongs to the folder before.
            structureRequestRef.current += 1;
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
                    return [file.path, await parseDashboardPlotFile(file)];
                }));
                if (cancelled) return;
                const parsedByPath = Object.fromEntries(parsedEntries);
                // Each xPDF transform's error ribbon, from the run's files
                // (ftUncertaintyRun.js; Flask sends the same in /api/plot/data).
                const nextFiles = await attachFtUncertainty(
                    loadedFiles.map((file) => parsedByPath[file.path] || file),
                    localRun.allFiles || loadedFiles
                );
                if (cancelled) return;
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
                    maxPoints: DASHBOARD_MAX_POINTS
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

    useEffect(() => {
        // Reset view state only when the folder changes, not on each Live Data refresh
        // (directory stays constant across live updates of the same run).
        setShowLoadedFiles(false);
        manuallyToggledPathsRef.current = new Set();
        setDismissedErrors(new Set());
        setRibbonChoices(NO_CHOICES);
        setChartData(NO_SERVER_PLOTS);
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

    // Every plot file in the Dashboard's order, hidden ones included (browserData:
    // the Copilot's run module reads the same list).
    const allPlotFiles = useMemo(() => dashboardPlotFiles(files), [files]);

    const plotFiles = useMemo(() => {
        return allPlotFiles.filter((file) => !hiddenPlotPaths.has(file.path));
    }, [allPlotFiles, hiddenPlotPaths]);

    // The run the Model card describes picks which log group is charted.
    const structurePath = localRun ? localRun.structureFile?.path : structure?.source;
    const rValueFile = useMemo(
        () => chartedRValueFile(allPlotFiles, hiddenPlotPaths, structurePath),
        [allPlotFiles, hiddenPlotPaths, structurePath]
    );
    const gridFiles = useMemo(
        () => plotFiles.filter((file) => file.plotKind !== 'r_value'),
        [plotFiles]
    );

    // Flask mode: the χ² chart's data is the server's concatenation of every
    // log of the charted run (related_r_value_logs, hidden logs included), so
    // its version is all of theirs: a later log a running job appends to
    // changes it, not only the first log's time and size.
    const rValueVersion = useMemo(() => {
        if (!rValueFile || localRun || isStaticMode()) return null;
        const run = rValueGroupKey(rValueFile.path);
        return fileSignature(allPlotFiles.filter((file) => file.plotKind === 'r_value' && rValueGroupKey(file.path) === run));
    }, [rValueFile, allPlotFiles, localRun]);

    // Flask mode, once the Copilot is open: what the published run fetches,
    // so the run module gets the plot data a picked folder's parsed files
    // carry. Every Dashboard plot file, hidden ones included (a picked folder
    // parses them all), at its own version; and the charted χ² log, whose
    // data is its run's whole history, at the run's version. Not the other
    // logs: each would bring the same whole history, which the run module
    // would concatenate again.
    const serverPlotTargets = useMemo(() => {
        if (!wantAssistantData || localRun || isStaticMode()) return NO_TARGETS;
        const targets = allPlotFiles
            .filter((file) => file.plotKind !== 'r_value')
            .map((file) => ({ file, version: plotDataVersion(file, files) }));
        if (rValueFile) targets.push({ file: rValueFile, version: rValueVersion });
        return targets;
    }, [wantAssistantData, localRun, allPlotFiles, files, rValueFile, rValueVersion]);

    // Flask mode: each xPDF chart's data (its own request, plotDataCache), for
    // its error ribbon, caption and Problems row. A chart that cannot load
    // lists that failure itself (InteractivePlot's plot-load issue).
    const xpdfTargets = useMemo(() => {
        if (localRun || isStaticMode()) return NO_TARGETS;
        return allPlotFiles
            .filter((file) => file.plotKind === 'xpdf')
            .map((file) => ({ file, version: plotDataVersion(file, files) }));
    }, [localRun, allPlotFiles, files]);

    useEffect(() => {
        let cancelled = false;
        xpdfTargets.forEach(({ file, version }) => {
            fetchPlotData(file.path, version).then(
                (plotData) => {
                    if (!cancelled) setChartData((current) => ({ ...current, [file.path]: { version, plotData } }));
                },
                () => {} // InteractivePlot lists the failure under this file's name.
            );
        });
        return () => { cancelled = true; };
    }, [xpdfTargets]);

    // The file an xPDF card's ribbon reads: with its data (a browser-loaded
    // run's own; in Flask mode the chart's latest answer), else null.
    const ribbonFile = useCallback((file) => {
        if (file.plotKind !== 'xpdf') return null;
        if (file.plotData) return file;
        const entry = chartData[file.path];
        return entry ? { ...file, plotData: entry.plotData } : null;
    }, [chartData]);

    // Each target is asked for once per version, through the cache the charts
    // share (plotDataCache.js): a chart already drawn costs no second request,
    // and a Live Data poll asks again only for the files whose version moved.
    // The answers of one round land together, as a picked folder's parses do.
    useEffect(() => {
        const requested = serverPlotRequestsRef.current;
        const asked = serverPlotTargets.filter(({ file, version }) => requested.get(file.path) !== version);
        if (!asked.length) return;
        asked.forEach(({ file, version }) => requested.set(file.path, version));
        Promise.all(asked.map(({ file, version }) => fetchPlotData(file.path, version).then(
            (plotData) => ({ path: file.path, version, entry: { plotData } }),
            (failure) => ({ path: file.path, version, entry: { parseError: failure.message || 'Could not load the plot data' } })
        ))).then((answers) => {
            const current = answers.filter(({ path, version }) => (
                requested === serverPlotRequestsRef.current && requested.get(path) === version
            ));
            if (!current.length) return;
            setServerPlots((plots) => ({
                ...plots,
                ...Object.fromEntries(current.map(({ path, entry }) => [path, entry]))
            }));
        });
    }, [serverPlotTargets]);

    // The files as the Copilot's run gets them: a target carries its latest
    // answer (the one before stays until a newer version's lands, as a picked
    // folder keeps a file's parse until its re-parse is done).
    const copilotFiles = useMemo(() => {
        if (!serverPlotTargets.length) return files;
        const targets = new Set(serverPlotTargets.map(({ file }) => file.path));
        return files.map((file) => {
            const entry = targets.has(file.path) ? serverPlots[file.path] : null;
            return entry ? { ...file, ...entry } : file;
        });
    }, [files, serverPlots, serverPlotTargets]);

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
    // Flask mode: a file whose data the Copilot's run could not get, hidden
    // or collapsed charts included (a chart on screen lists the same failure
    // under the same file name, and the section shows it once).
    const serverPlotIssues = useMemo(() => serverPlotTargets
        .filter(({ file }) => serverPlots[file.path]?.parseError)
        .map(({ file }) => ({ key: file.path, message: serverPlots[file.path].parseError, source: file.name })),
    [serverPlotTargets, serverPlots]);
    useIssueSet('plot-data', serverPlotIssues);
    // An xPDF plot whose error ribbon failed (a reproduction test, a check or
    // a read): the caption says why, and so does the Problems section. A ribbon
    // that is merely not possible (no σ, a missing file) is the caption's alone.
    const ribbonProblemRows = useMemo(
        () => ribbonIssues(allPlotFiles.map(ribbonFile).filter(Boolean)),
        [allPlotFiles, ribbonFile]
    );
    useIssueSet('ribbon', ribbonProblemRows);
    const reportIssue = useReportIssue();

    // Publish the run as this page holds it for App's Copilot run module
    // (useCopilotRun.js), which turns it into what the AI Copilot is handed:
    // the symmetry at the shared tolerance with its ladder, the χ² history,
    // the parsed .dat. Only once the Copilot has been opened, so Dashboard- or
    // KDE-only use does no extra work; the WatchdogBadge below stays on the
    // dashboard. The hidden list keeps its identity while the hidden set does.
    // In Flask mode the files carry the plot data fetched above.
    const hiddenPlots = useMemo(() => [...hiddenPlotPaths], [hiddenPlotPaths]);
    const publishedRun = useMemo(() => (wantAssistantData ? {
        files: copilotFiles,
        structure,
        structurePath: structurePath ?? null,
        hiddenPlots,
    } : null), [wantAssistantData, copilotFiles, structure, structurePath, hiddenPlots]);

    useEffect(() => {
        onRunChange?.(publishedRun);
    }, [publishedRun, onRunChange]);

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

    const renderPlotBody = (file, variant, marks = null, refreshKey = plotVersion(file), ribbon = null) => {
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
                refreshKey={refreshKey}
                variant={variant}
                marks={marks}
                ribbon={ribbon}
            />
        );
    };

    const renderPlotCard = (file) => {
        const meta = metadata[file.path];
        const title = meta?.title || file.name;
        // An xPDF plot's error ribbon: the file's σ, else the propagated one
        // (GrRibbon.jsx), with the toggle when both exist and its caption.
        const ribbon = file.plotKind === 'xpdf' ? ribbonSources(ribbonFile(file)) : null;
        const choice = ribbonChoices[file.path] || null;
        const chooseRibbon = (value) => setRibbonChoices((current) => ({ ...current, [file.path]: value }));
        return (
            <Card as="article" clip lift data-figure-card="" data-copilot-card={`plot:${file.name}`} key={file.path}>
                <div className="ui-card__header-flush">
                    <div className="ui-card__heading">
                        <CardTitle>{title}</CardTitle>
                        {title !== file.name && (
                            <span className="ui-card__source" title={file.path}>{file.name}</span>
                        )}
                    </div>
                    {ribbon?.fileSigma && ribbon?.propagated ? (
                        <div className="ui-card__header-actions">
                            <RibbonToggle sources={ribbon} choice={choice} onChoose={chooseRibbon} />
                            {renderRwpChip(meta)}
                        </div>
                    ) : renderRwpChip(meta)}
                </div>
                {renderPlotBody(
                    file, undefined, cardMarks[`plot:${file.name}`] || null, plotDataVersion(file, files), ribbonProp(ribbon, choice)
                )}
                {ribbon && <RibbonNote sources={ribbon} choice={choice} />}
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
                data-copilot-card="chi2-history"
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
                {showRValue && renderPlotBody(
                    rValueFile, 'wide', cardMarks['chi2-history'] || null, rValueVersion ?? plotVersion(rValueFile)
                )}
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
