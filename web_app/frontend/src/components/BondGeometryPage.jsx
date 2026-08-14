// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Bond Geometry page: bond-angle (triplet) distribution plus the bond-length
// and coordination statistics that fall out of the same neighbour search —
// the RMCProfile `triplets` workflow. Pick an A–B–C triplet with B central,
// bound the two bond lengths, and Compute histograms the angle at B.
//
// Design mirrors the PCA Ellipsoid / Orientation pages: all options in the
// top controls bar, results in pca-panel cards. The engine runs in the shared
// worker for browser-loaded runs (both runtimes) or via /api/triplets for a
// typed backend directory — identical payloads either way (source of truth:
// rmc_toolkits/triplets.py; port: workers/triplets.js).

import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import { isStaticMode, readAndParseLocalPlotFile } from '../browserData';
import { buildElementColors } from '../atomColors';
import { PLOT_PALETTE } from '../plotPalette';
import InfoBadge from './InfoBadge';
import InteractivePlot from './InteractivePlot';
import ModelSummary from './ModelSummary';
import FoldedCellPanel from './FoldedCellPanel';
import useSiteCloud from '../useSiteCloud';
import './PcaKdePage.css';
import './BondGeometryPage.css';

// Same cap as StructurePage/Dashboard: the Model information + Detected SG
// cards need the full basis and counts, and parsing is worker-side anyway.
const STRUCTURE_MAX_POINTS = 1000000;

const DEGREES = '°';
const ANGSTROM = 'Å';

const formatNumber = (value, digits = 2) =>
    Number.isFinite(value) ? value.toFixed(digits) : '—';

// Trailing-debounced copy of a value: the helper plot's window guides follow
// typing only after a pause, so each keystroke doesn't rebuild the plot (and
// reset its zoom/hover state, which is keyed on plotData identity).
const useDebounced = (value, delay = 400) => {
    const [debounced, setDebounced] = useState(value);
    useEffect(() => {
        const timer = setTimeout(() => setDebounced(value), delay);
        return () => clearTimeout(timer);
    }, [value, delay]);
    return debounced;
};

export default function BondGeometryPage({ directory, localRun }) {
    const {
        sites,
        sitesError,
        loadingSites,
        requestPca,
        localFile,
        ready,
        datasetKey
    } = useSiteCloud({ directory, localRun });

    const elements = useMemo(() => sites?.elements ?? [], [sites]);
    const elementColors = useMemo(() => buildElementColors(sites?.elements ?? []), [sites]);

    // Triplet selection: A–B–C with B central. Initialized per dataset once
    // the element list arrives; RMCProfile-style, ends default to the last
    // element (commonly the anion) around the first differing central one.
    const [end1, setEnd1] = useState('');
    const [apex, setApex] = useState('');
    const [end2, setEnd2] = useState('');
    // Seed once per *sites payload* (not per dataset key): on a dataset
    // switch the new key arrives while the previous run's sites are still in
    // state, so keying on the data itself waits for the real element list —
    // and a Live Data refresh keeps the user's picks when they still apply.
    const seededSites = useRef(null);
    useEffect(() => {
        if (!sites?.sites || sites === seededSites.current) return;
        seededSites.current = sites;
        const valid = elements.length
            && [end1, apex, end2].every((element) => elements.includes(element));
        if (valid) return;
        // Ends default to the most abundant element (in practice the anion),
        // the central atom to the next most abundant — Se–Nb–Se for GaNb4Se8.
        const totals = new Map(elements.map((element) => [element, 0]));
        sites.sites.forEach((site) => {
            totals.set(site.element, (totals.get(site.element) ?? 0) + site.count);
        });
        const ranked = [...totals.entries()].sort((a, b) => b[1] - a[1]).map(([element]) => element);
        const end = ranked[0];
        const center = ranked.find((element) => element !== end) ?? end;
        setEnd1(end);
        setApex(center);
        setEnd2(end);
    }, [sites, elements, end1, apex, end2]);

    // Bond windows (Å, inclusive). The B–C window follows A–B unless split.
    const [r12Min, setR12Min] = useState('2.0');
    const [r12Max, setR12Max] = useState('3.0');
    const [split23, setSplit23] = useState(false);
    const [r23Min, setR23Min] = useState('2.0');
    const [r23Max, setR23Max] = useState('3.0');
    const [binWidth, setBinWidth] = useState('1.0');

    const [result, setResult] = useState(null);
    const [resultError, setResultError] = useState(null);
    const [computing, setComputing] = useState(false);
    const [angleView, setAngleView] = useState('sin');

    const staticMode = isStaticMode();
    const noRun = staticMode && !localFile;

    // A dataset switch clears any previous result immediately and bumps the
    // epoch, so a compute that was in flight for the old run can never land
    // its (stale) payload on the new one.
    const runEpoch = useRef(0);
    useEffect(() => {
        runEpoch.current += 1;
        setResult(null);
        setResultError(null);
    }, [datasetKey]);

    const compute = useCallback(async () => {
        const epoch = runEpoch.current;
        setComputing(true);
        setResultError(null);
        try {
            const params = {
                end1,
                apex,
                end2,
                r12Min: Number(r12Min),
                r12Max: Number(r12Max),
                binWidth: Number(binWidth)
            };
            if (split23) {
                params.r23Min = Number(r23Min);
                params.r23Max = Number(r23Max);
            }
            const data = await requestPca('triplets', params);
            if (runEpoch.current === epoch) setResult(data);
        } catch (error) {
            if (runEpoch.current === epoch) {
                setResult(null);
                setResultError(error.message);
            }
        } finally {
            if (runEpoch.current === epoch) setComputing(false);
        }
    }, [requestPca, end1, apex, end2, r12Min, r12Max, split23, r23Min, r23Max, binWidth]);

    // --- Structure for the Model information / Detected SG cards. -----------
    // Same source as the Dashboard and Atomic Density pages: a local run's
    // .rmc6f parses in the structure worker (both runtimes); a typed backend
    // directory asks /api/structure. ModelSummary then renders the model card
    // and, when the payload carries a basis, the Detected SG card.
    const [structure, setStructure] = useState(null);
    const structureWorkerRef = useRef(null);
    const structureRequestRef = useRef(0);
    // Null the ref on terminate: StrictMode's dev double-mount runs this
    // cleanup between the two mounts, and a ref still holding the terminated
    // worker would make the second mount post requests into a dead worker.
    useEffect(() => () => {
        structureWorkerRef.current?.terminate();
        structureWorkerRef.current = null;
    }, []);
    useEffect(() => {
        let cancelled = false;
        if (localRun) {
            if (localRun.structure) {
                setStructure(localRun.structure);
                return undefined;
            }
            if (!localRun.structureFile) {
                setStructure(null);
                return undefined;
            }
            if (!structureWorkerRef.current) {
                structureWorkerRef.current = new Worker(
                    new URL('../workers/localStructureWorker.js', import.meta.url),
                    { type: 'module' }
                );
            }
            const worker = structureWorkerRef.current;
            const id = structureRequestRef.current + 1;
            structureRequestRef.current = id;
            setStructure(null);
            const onMessage = (event) => {
                if (event.data.id !== id) return;
                worker.removeEventListener('message', onMessage);
                if (!cancelled) setStructure(event.data.error ? null : event.data.result);
            };
            worker.addEventListener('message', onMessage);
            worker.postMessage({ id, file: localRun.structureFile, maxPoints: STRUCTURE_MAX_POINTS });
            return () => {
                cancelled = true;
                worker.removeEventListener('message', onMessage);
            };
        }
        if (isStaticMode()) {
            setStructure(null);
            return undefined;
        }
        axios
            .get(`${API_BASE_URL}/api/structure`, {
                params: { dir: directory || '.', maxPoints: STRUCTURE_MAX_POINTS }
            })
            .then((response) => { if (!cancelled) setStructure(response.data); })
            .catch(() => { if (!cancelled) setStructure(null); });
        return () => { cancelled = true; };
    }, [localRun, directory, datasetKey]);

    // --- Partial g(r) for the window helper. ---------------------------------
    // The run's PDFpartials.csv (when present) shows where the first
    // coordination shell ends, which is how the bond windows are chosen.
    const [partials, setPartials] = useState(null);
    const [partialsLoading, setPartialsLoading] = useState(false);
    useEffect(() => {
        let cancelled = false;
        setPartials(null);
        setPartialsLoading(true);
        const load = async () => {
            try {
                if (localRun) {
                    const file = localRun.files?.find((entry) => entry.plotKind === 'pdf_partials');
                    if (!file) return;
                    const parsed = file.plotData ?? (await readAndParseLocalPlotFile(file));
                    if (!cancelled) setPartials(parsed);
                    return;
                }
                // Static mode with no run loaded: there is no backend to ask.
                if (isStaticMode()) return;
                const listing = await axios.get(`${API_BASE_URL}/api/files`, {
                    params: { dir: directory || '.' }
                });
                const file = (listing.data.files ?? []).find(
                    (entry) => entry.plotKind === 'pdf_partials'
                );
                if (!file) return;
                const parsed = await axios.get(`${API_BASE_URL}/api/plot/data`, {
                    params: { path: file.path }
                });
                if (!cancelled) setPartials(parsed.data);
            } catch {
                if (!cancelled) setPartials(null);
            } finally {
                if (!cancelled) setPartialsLoading(false);
            }
        };
        load();
        return () => { cancelled = true; };
    }, [localRun, directory, datasetKey]);

    // PDFpartials.csv labels a pair in one order only, so try both.
    const findPartial = useCallback((a, b) => {
        if (!partials || !a || !b) return null;
        return (
            partials.series?.find((series) => series.label === `${a}-${b}`) ??
            partials.series?.find((series) => series.label === `${b}-${a}`) ??
            null
        );
    }, [partials]);

    const partialSeries = useMemo(() => findPartial(end1, apex), [findPartial, end1, apex]);
    // A second curve whenever the two bonds are different types: Ga–Ta–Se has
    // a Ga-Ta and a Ta-Se shell to bracket, Se–Ta–Se only ever has the one.
    // Independent of the window split — the shells differ either way.
    const partialCurve23 = useMemo(() => {
        const found = findPartial(apex, end2);
        return found && found !== partialSeries ? found : null;
    }, [findPartial, apex, end2, partialSeries]);

    // --- Plot data (memoized: a new object identity resets plot view state). --
    const anglePlot = useMemo(() => {
        if (!result) return null;
        const sin = angleView === 'sin';
        return {
            title: `${result.triplet.join('-')} bond angles`,
            xLabel: `angle at ${result.triplet[1]} (deg)`,
            yLabel: sin ? 'sin-corrected' : 'density (deg^{-1})',
            series: [
                {
                    label: sin ? 'sin-corrected' : 'density',
                    x: result.binCenters,
                    y: sin ? result.sinCorrected : result.density
                }
            ]
        };
    }, [result, angleView]);

    // Guides follow the inputs after a typing pause, so the plot (whose view
    // state resets on every plotData identity change) stays stable per key.
    const guideMin = useDebounced(r12Min);
    const guideMax = useDebounced(r12Max);
    const guide23Min = useDebounced(r23Min);
    const guide23Max = useDebounced(r23Max);
    const helperPlot = useMemo(() => {
        if (!partialSeries) return null;
        // The helper is about choosing the first-shell window, so show the
        // short-range part only — beyond it nothing informs a bond window.
        // With two windows on screen the crop has to clear the further one.
        const furthest = split23
            ? Math.max(Number(guideMax) || 0, Number(guide23Max) || 0)
            : Number(guideMax) || 0;
        const limit = Math.max(6, furthest * 2);
        const crop = (series) => {
            const cut = series.x.findIndex((value) => value > limit);
            const end = cut === -1 ? series.x.length : cut;
            return {
                label: series.label,
                x: series.x.slice(0, end),
                y: series.y.slice(0, end)
            };
        };
        // One curve per distinct bond type; same-type triplets draw it once, or
        // the series label would be duplicated.
        const shells = [crop(partialSeries)];
        if (partialCurve23) shells.push(crop(partialCurve23));
        const max = shells.reduce(
            (acc, shell) => shell.y.reduce((inner, value) => Math.max(inner, value), acc),
            0
        ) || 1;
        // The windows are what the split governs: off, both bonds use the A–B
        // bounds and one pair of guides covers them. Labelled by role, not by
        // pair — with the same element at both ends the pair names would match.
        // Two windows also get the two curve colors (guides take no palette
        // slot, so shell N is PLOT_PALETTE[N]): where the bonds are different
        // types each window is then the color of the shell it brackets. A lone
        // window keeps the neutral guide stroke — there is nothing to tell apart.
        const windows = split23
            ? [
                { prefix: 'A–B ', min: guideMin, max: guideMax, color: PLOT_PALETTE[0] },
                { prefix: 'B–C ', min: guide23Min, max: guide23Max, color: PLOT_PALETTE[1] }
            ]
            : [{ prefix: '', min: guideMin, max: guideMax, color: undefined }];
        const guides = windows
            .flatMap((window) => [
                { label: `${window.prefix}rmin ${window.min}`, value: Number(window.min), color: window.color },
                { label: `${window.prefix}rmax ${window.max}`, value: Number(window.max), color: window.color }
            ])
            .filter((guide) => Number.isFinite(guide.value))
            .map((guide) => ({
                label: guide.label,
                x: [guide.value, guide.value],
                y: [0, max],
                role: 'guide',
                color: guide.color
            }));
        return {
            title: partialCurve23
                ? `${partialSeries.label} and ${partialCurve23.label} partial g(r)`
                : `${partialSeries.label} partial g(r)`,
            xLabel: `r (${ANGSTROM})`,
            yLabel: 'g(r)',
            series: [...shells, ...guides]
        };
    }, [partialSeries, partialCurve23, split23, guideMin, guideMax, guide23Min, guide23Max]);

    // Coordination summary: "5.78 per B — 6-fold 82.8%".
    const coordinationSummary = useMemo(() => {
        if (!result) return null;
        const total = result.coordination.reduce((acc, value) => acc + value, 0);
        const bonds = result.coordination.reduce((acc, value, n) => acc + value * n, 0);
        if (!total) return null;
        let best = 0;
        result.coordination.forEach((value, n) => {
            if (value > result.coordination[best]) best = n;
        });
        return {
            mean: bonds / total,
            mode: best,
            modeShare: (100 * result.coordination[best]) / total
        };
    }, [result]);

    const windowLabel = (window) => `${formatNumber(window[0])}–${formatNumber(window[1])} ${ANGSTROM}`;

    // Detected-bond overlay for the unit-cell panel: the computed windows,
    // A–B in the app accent, B–C in amber when the windows are distinct.
    const bondSets = useMemo(() => {
        if (!result) return null;
        const sets = [
            { elements: [result.triplet[0], result.triplet[1]], window: result.bond12, color: 0x2563eb }
        ];
        if (!result.sharedEnds) {
            sets.push({ elements: [result.triplet[1], result.triplet[2]], window: result.bond23, color: 0xd97706 });
        }
        return sets;
    }, [result]);

    return (
        <div className="pca-page">
            {/* Model information first, as on the Dashboard; the Detected SG
                card stays on the Dashboard/Atomic Density pages only. */}
            {structure && (
                <div className="geom-model">
                    <ModelSummary structure={structure} showSymmetry={false} />
                </div>
            )}
            <div className="pca-controls">
                <div className="control-group" role="group" aria-label="Triplet">
                    <label className="control">
                        <span className="control-name">
                            A
                            <InfoBadge label="About the triplet">
                                <p>
                                    The RMCProfile <code>triplets</code> convention: the angle is
                                    measured at the <strong>central atom B</strong> between its bond
                                    to an A atom and its bond to a C atom. A and C may name the same
                                    element (each unordered pair of bonds then counts once).
                                </p>
                            </InfoBadge>
                        </span>
                        <select value={end1} onChange={(event) => setEnd1(event.target.value)} disabled={!elements.length} aria-label="End element A">
                            {elements.map((element) => <option key={element} value={element}>{element}</option>)}
                        </select>
                    </label>
                    <label className="control">
                        <span className="control-name">B (central)</span>
                        <select value={apex} onChange={(event) => setApex(event.target.value)} disabled={!elements.length} aria-label="Central element B" className="geom-central">
                            {elements.map((element) => <option key={element} value={element}>{element}</option>)}
                        </select>
                    </label>
                    <label className="control">
                        <span className="control-name">C</span>
                        <select value={end2} onChange={(event) => setEnd2(event.target.value)} disabled={!elements.length} aria-label="End element C">
                            {elements.map((element) => <option key={element} value={element}>{element}</option>)}
                        </select>
                    </label>
                </div>

                <div className="control-group" role="group" aria-label="Bond windows">
                    <label className="control">
                        <span className="control-name">
                            A{'–'}B window ({ANGSTROM})
                            <InfoBadge label="About the bond windows">
                                <p>
                                    Two atoms are bonded when their distance falls inside the window
                                    (inclusive). Read the window off the first-shell peak of the
                                    partial g(r) — the helper panel shades it once a partials file is
                                    in the run folder.
                                </p>
                            </InfoBadge>
                        </span>
                        <span className="geom-window">
                            <input type="number" step="0.05" min="0" max="15" value={r12Min} onChange={(event) => setR12Min(event.target.value)} aria-label="A-B window minimum" />
                            {'–'}
                            <input type="number" step="0.05" min="0" max="15" value={r12Max} onChange={(event) => setR12Max(event.target.value)} aria-label="A-B window maximum" />
                        </span>
                    </label>
                    <label className="control switch">
                        <span className="control-name">Distinct B{'–'}C</span>
                        <input type="checkbox" checked={split23} onChange={(event) => setSplit23(event.target.checked)} aria-label="Use a distinct B-C window" />
                        <i className="switch-track" aria-hidden="true" />
                    </label>
                    {split23 && (
                        <label className="control">
                            <span className="control-name">B{'–'}C window ({ANGSTROM})</span>
                            <span className="geom-window">
                                <input type="number" step="0.05" min="0" max="15" value={r23Min} onChange={(event) => setR23Min(event.target.value)} aria-label="B-C window minimum" />
                                {'–'}
                                <input type="number" step="0.05" min="0" max="15" value={r23Max} onChange={(event) => setR23Max(event.target.value)} aria-label="B-C window maximum" />
                            </span>
                        </label>
                    )}
                </div>

                <div className="control-group" role="group" aria-label="Histogram">
                    <label className="control">
                        <span className="control-name">Bin width</span>
                        <span className="geom-window">
                            <input type="number" step="0.5" min="0.1" max="45" value={binWidth} onChange={(event) => setBinWidth(event.target.value)} aria-label="Angle bin width in degrees" />
                            <span className="control-value">deg</span>
                        </span>
                    </label>
                    <button
                        type="button"
                        className="geom-compute"
                        onClick={compute}
                        disabled={computing || !ready || !elements.length || !end1 || !apex || !end2}
                    >
                        {computing ? 'Computing…' : 'Compute'}
                    </button>
                </div>
            </div>

            {noRun && (
                <p className="pca-hint">Open a run folder (with an <code>.rmc6f</code> file) to analyse bond angles.</p>
            )}
            {(sitesError || resultError) && <p className="pca-error-banner">{sitesError || resultError}</p>}
            {!noRun && !result && !resultError && !sitesError && (
                <p className="pca-hint">
                    Pick the A{'–'}B{'–'}C triplet (B central), bound the bond lengths, then
                    Compute. Angles are counted over the periodic configuration exactly, images included.
                </p>
            )}

            {result && (
                <div className="model-cards" role="status">
                    {/* Same presentation as the Model information card: labeled
                        columns, not badges. */}
                    <section className="model-summary" aria-label="Triplet result">
                        <h2 className="model-summary-title">
                            Triplet result
                            <span className="model-summary-source">
                                {`${result.triplet.join('–')} · ${windowLabel(result.bond12)}`}
                            </span>
                        </h2>
                        <dl className="model-stats">
                            <div className="model-stat">
                                <dt>Central atoms</dt>
                                <dd>
                                    {result.apexCount.toLocaleString()}
                                    <span className="model-stat-sub">{result.triplet[1]}</span>
                                </dd>
                            </div>
                            <div className="model-stat">
                                <dt>{result.sharedEnds ? 'Bonds' : 'Bonds A–B'}</dt>
                                <dd>
                                    {result.lengths12.count.toLocaleString()}
                                    {result.lengths12.meanLength != null && (
                                        <span className="model-stat-sub">
                                            mean {formatNumber(result.lengths12.meanLength, 3)} {ANGSTROM}
                                        </span>
                                    )}
                                </dd>
                            </div>
                            {!result.sharedEnds && result.lengths23 && (
                                <div className="model-stat">
                                    <dt>Bonds B–C</dt>
                                    <dd>
                                        {result.lengths23.count.toLocaleString()}
                                        {result.lengths23.meanLength != null && (
                                            <span className="model-stat-sub">
                                                mean {formatNumber(result.lengths23.meanLength, 3)} {ANGSTROM}
                                            </span>
                                        )}
                                    </dd>
                                </div>
                            )}
                            {coordinationSummary && (
                                <div className="model-stat">
                                    <dt>Coordination</dt>
                                    <dd>
                                        {formatNumber(coordinationSummary.mean)}
                                        <span className="model-stat-sub">
                                            {`per ${result.triplet[1]} · ${coordinationSummary.mode}-fold ${formatNumber(coordinationSummary.modeShare, 1)}%`}
                                        </span>
                                    </dd>
                                </div>
                            )}
                            <div className="model-stat">
                                <dt>Angles</dt>
                                <dd>
                                    {result.angleCount.toLocaleString()}
                                    {result.meanAngle != null && (
                                        <span className="model-stat-sub">
                                            {`mean ${formatNumber(result.meanAngle, 1)}${DEGREES} ± ${formatNumber(result.stdAngle, 1)}${DEGREES}`}
                                        </span>
                                    )}
                                </dd>
                            </div>
                        </dl>
                    </section>
                </div>
            )}

            <div className="geom-layout">
                <div className="pca-panel geom-main-panel">
                    <h3>
                        <span className="panel-title-label">
                            {result ? `${result.triplet.join('–')} bond-angle distribution` : 'Bond-angle distribution'}
                            <InfoBadge label="About the two normalizations">
                                <p>
                                    <b>Density</b> — the raw distribution: probability per degree,
                                    area 1. Note that randomly oriented bonds do <em>not</em> give a
                                    flat line here. There are simply fewer ways to form an angle near
                                    0{DEGREES} or 180{DEGREES} than near 90{DEGREES}, so the geometry
                                    alone bows the curve.
                                </p>
                                <p>
                                    <b>Sin-corrected</b> — divides that geometric factor out. Random
                                    bonds now read as a flat 1, anything above it is real structure,
                                    and a peak near 180{DEGREES} is no longer flattened. This is
                                    RMCProfile's <code>sinth</code> view.
                                </p>
                            </InfoBadge>
                        </span>
                        <span className="panel-title-actions">
                            {/* The angle count lives in the Triplet result card;
                                the header only carries the toggle. */}
                            <div className="pca-frame-toggle" role="group" aria-label="Angle normalization">
                                <button type="button" className={angleView === 'sin' ? 'is-active' : ''} onClick={() => setAngleView('sin')}>
                                    sin-corrected
                                </button>
                                <button type="button" className={angleView === 'raw' ? 'is-active' : ''} onClick={() => setAngleView('raw')}>
                                    density
                                </button>
                            </div>
                        </span>
                    </h3>
                    <div className="geom-plot">
                        {anglePlot
                            ? <InteractivePlot file={{ path: `geometry:angles:${datasetKey}`, name: anglePlot.title }} plotData={anglePlot} variant="fit" />
                            : <p className="geom-empty">No distribution yet.</p>}
                    </div>
                </div>

                <div className="pca-panel">
                        <h3>
                            <span className="panel-title-label">
                                Partial PDF
                                <InfoBadge label="About the partial PDF">
                                    <p>
                                        The A{'–'}B partial pair distribution from the run's{' '}
                                        <code>PDFpartials.csv</code>. Set the bond window to bracket
                                        the first-shell peak; the dashed guides track the current
                                        A{'–'}B bounds. With <b>distinct B{'–'}C</b> on, the B{'–'}C
                                        partial is plotted alongside it with its own pair of guides,
                                        so both windows can be set against their own shell.
                                    </p>
                                </InfoBadge>
                            </span>
                            {partialSeries && (
                                <span className="panel-title-count">
                                    {partialCurve23
                                        ? `${partialSeries.label} · ${partialCurve23.label}`
                                        : partialSeries.label}
                                </span>
                            )}
                        </h3>
                        <div className="geom-plot">
                            {helperPlot
                                ? <InteractivePlot file={{ path: `geometry:partial:${datasetKey}`, name: 'partial-gr' }} plotData={helperPlot} variant="fit" />
                                : (
                                    <p className="geom-empty">
                                        {partialsLoading
                                            ? 'Looking for a partials file…'
                                            : <>No <code>PDFpartials.csv</code> in this run.</>}
                                    </p>
                                )}
                        </div>
                </div>

                {/* The Atomic Density page's folded cell — the atom cloud, not
                    fitted ellipsoids — with the detected bonds over it. */}
                <FoldedCellPanel
                    title="Detected bonds"
                    structure={structure}
                    sites={sites}
                    elementColors={elementColors}
                    loading={loadingSites}
                    bondSets={bondSets}
                />
            </div>

            <footer className="app-footer">
                &copy; 2026 Tsung-Han Yang &middot;{' '}
                <a href="https://github.com/drthyang/rmc-toolkits/blob/main/LICENSE" target="_blank" rel="noreferrer">AGPLv3</a>
                {' '}&middot;{' '}
                <a href="https://github.com/drthyang/rmc-toolkits#readme" target="_blank" rel="noreferrer">About & documentation</a>
            </footer>
        </div>
    );
}
