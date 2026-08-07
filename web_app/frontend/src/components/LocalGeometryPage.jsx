// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Local Geometry page: bond-angle (triplet) distribution plus the bond-length
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
import InfoBadge from './InfoBadge';
import InteractivePlot from './InteractivePlot';
import useSiteCloud from '../useSiteCloud';
import './PcaKdePage.css';
import './LocalGeometryPage.css';

const DEGREES = '°';
const ANGSTROM = 'Å';

// Outline of a histogram as a step polyline (InteractivePlot draws lines
// only): duplicated x at every bin edge, level y across each bin.
const stepSeries = (binCenters, counts) => {
    if (!binCenters?.length) return { x: [], y: [] };
    const width = binCenters.length > 1 ? binCenters[1] - binCenters[0] : 1;
    const x = [];
    const y = [];
    binCenters.forEach((center, index) => {
        x.push(center - width / 2, center + width / 2);
        y.push(counts[index], counts[index]);
    });
    return { x, y };
};

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

export default function LocalGeometryPage({ directory, localRun }) {
    const { sites, sitesError, requestPca, localFile, ready, datasetKey } =
        useSiteCloud({ directory, localRun });

    const elements = useMemo(() => sites?.elements ?? [], [sites]);

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

    const partialSeries = useMemo(() => {
        if (!partials || !end1 || !apex) return null;
        return (
            partials.series?.find((series) => series.label === `${end1}-${apex}`) ??
            partials.series?.find((series) => series.label === `${apex}-${end1}`) ??
            null
        );
    }, [partials, end1, apex]);

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

    const lengthsPlot = useMemo(() => {
        if (!result) return null;
        const series = [];
        const first = stepSeries(result.lengths12.binCenters, result.lengths12.counts);
        series.push({
            label: result.sharedEnds
                ? `${result.triplet[0]}-${result.triplet[1]}`
                : `${result.triplet[0]}-${result.triplet[1]} (1-2)`,
            x: first.x,
            y: first.y
        });
        if (result.lengths23) {
            const second = stepSeries(result.lengths23.binCenters, result.lengths23.counts);
            series.push({
                label: `${result.triplet[1]}-${result.triplet[2]} (2-3)`,
                x: second.x,
                y: second.y
            });
        }
        return {
            title: `bond lengths in window`,
            xLabel: `length (${ANGSTROM})`,
            yLabel: 'bonds',
            series
        };
    }, [result]);

    // Guides follow the inputs after a typing pause, so the plot (whose view
    // state resets on every plotData identity change) stays stable per key.
    const guideMin = useDebounced(r12Min);
    const guideMax = useDebounced(r12Max);
    const helperPlot = useMemo(() => {
        if (!partialSeries) return null;
        // The helper is about choosing the first-shell window, so show the
        // short-range part only — beyond it nothing informs a bond window.
        const limit = Math.max(6, Number(guideMax) * 2 || 0);
        const cut = partialSeries.x.findIndex((value) => value > limit);
        const end = cut === -1 ? partialSeries.x.length : cut;
        const shell = {
            label: partialSeries.label,
            x: partialSeries.x.slice(0, end),
            y: partialSeries.y.slice(0, end)
        };
        const max = shell.y.reduce((acc, value) => Math.max(acc, value), 0) || 1;
        const guides = [
            { label: `rmin ${guideMin}`, value: Number(guideMin) },
            { label: `rmax ${guideMax}`, value: Number(guideMax) }
        ]
            .filter((guide) => Number.isFinite(guide.value))
            .map((guide) => ({
                label: guide.label,
                x: [guide.value, guide.value],
                y: [0, max],
                role: 'guide'
            }));
        return {
            title: `${partialSeries.label} partial g(r)`,
            xLabel: `r (${ANGSTROM})`,
            yLabel: 'g(r)',
            series: [shell, ...guides]
        };
    }, [partialSeries, guideMin, guideMax]);

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

    return (
        <div className="pca-page">
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
                <div className="geom-chips" role="status">
                    <span className="geom-chip">
                        <span className="geom-chip-name">Central atoms</span>
                        <strong>{result.apexCount.toLocaleString()}</strong> {result.triplet[1]}
                    </span>
                    <span className="geom-chip">
                        <span className="geom-chip-name">Bonds {windowLabel(result.bond12)}</span>
                        <strong>{result.lengths12.count.toLocaleString()}</strong>
                        {result.lengths12.meanLength != null && ` mean ${formatNumber(result.lengths12.meanLength, 3)} ${ANGSTROM}`}
                    </span>
                    {coordinationSummary && (
                        <span className="geom-chip">
                            <span className="geom-chip-name">Coordination</span>
                            <strong>{formatNumber(coordinationSummary.mean)}</strong>
                            {` per ${result.triplet[1]} · ${coordinationSummary.mode}-fold ${formatNumber(coordinationSummary.modeShare, 1)}%`}
                        </span>
                    )}
                    <span className="geom-chip">
                        <span className="geom-chip-name">Angles</span>
                        <strong>{result.angleCount.toLocaleString()}</strong>
                        {result.meanAngle != null &&
                            ` mean ${formatNumber(result.meanAngle, 1)}${DEGREES} ± ${formatNumber(result.stdAngle, 1)}${DEGREES}`}
                    </span>
                </div>
            )}

            <div className="geom-layout">
                <div className="pca-panel geom-main-panel">
                    <h3>
                        <span className="panel-title-label">
                            {result ? `${result.triplet.join('–')} bond-angle distribution` : 'Bond-angle distribution'}
                            <InfoBadge label="About the two normalizations">
                                <p>
                                    Sin-corrected divides each bin by the exact isotropic reference
                                    (flat 1 = random directions) — the RMCProfile{' '}
                                    <code>sinth</code> view, which keeps angles near 180{DEGREES}{' '}
                                    from being suppressed by geometry. Density is the raw per-degree
                                    probability with unit integral.
                                </p>
                            </InfoBadge>
                        </span>
                        <span className="panel-title-actions">
                            {result && <span className="panel-title-count">{result.angleCount.toLocaleString()} angles</span>}
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
                            ? <InteractivePlot file={{ path: `geometry:angles:${datasetKey}`, name: anglePlot.title }} plotData={anglePlot} />
                            : <p className="geom-empty">No distribution yet.</p>}
                    </div>
                </div>

                <div className="geom-side">
                    <div className="pca-panel">
                        <h3>
                            <span className="panel-title-label">Bond lengths in window</span>
                            {result?.lengths12 && (
                                <span className="panel-title-count">
                                    {result.sharedEnds
                                        ? `${result.lengths12.count.toLocaleString()} bonds`
                                        : `${result.lengths12.count.toLocaleString()} + ${result.lengths23.count.toLocaleString()} bonds`}
                                </span>
                            )}
                        </h3>
                        <div className="geom-plot">
                            {lengthsPlot
                                ? <InteractivePlot file={{ path: `geometry:lengths:${datasetKey}`, name: 'bond-lengths' }} plotData={lengthsPlot} />
                                : <p className="geom-empty">Compute to see the length histogram.</p>}
                        </div>
                    </div>
                    <div className="pca-panel">
                        <h3>
                            <span className="panel-title-label">
                                Window helper
                                <InfoBadge label="About the window helper">
                                    <p>
                                        The A{'–'}B partial pair distribution from the run's{' '}
                                        <code>PDFpartials.csv</code>. Set the window to bracket the
                                        first-shell peak; the dashed guides track the current
                                        A{'–'}B bounds.
                                    </p>
                                </InfoBadge>
                            </span>
                            {partialSeries && <span className="panel-title-count">{partialSeries.label}</span>}
                        </h3>
                        <div className="geom-plot">
                            {helperPlot
                                ? <InteractivePlot file={{ path: `geometry:partial:${datasetKey}`, name: 'partial-gr' }} plotData={helperPlot} />
                                : (
                                    <p className="geom-empty">
                                        {partialsLoading
                                            ? 'Looking for a partials file…'
                                            : <>No <code>PDFpartials.csv</code> in this run.</>}
                                    </p>
                                )}
                        </div>
                    </div>
                </div>
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
