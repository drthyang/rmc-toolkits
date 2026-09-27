// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Bond Geometry page: bond-angle (triplet) distribution plus the bond-length
// and coordination statistics that fall out of the same neighbour search —
// the RMCProfile `triplets` workflow. Pick an A–B–C triplet with B central,
// bound the two bond lengths, and Compute histograms the angle at B.
//
// Layout: the controls bar is a form (Enter in any field computes); below it
// the angle distribution is the hero card, with the headline results in a KPI
// rail under its header, and the folded cell and the partial g(r) stack on
// the right. The cards keep their skeleton before and after Compute: the
// hero shows a ghost of the angle axis with a prompt until the first result.
// One colour system runs through the three cards: element colours for atoms,
// BOND_COLORS for the A–B / B–C bond roles, neutral grey for references.
//
// The engine runs in the shared worker for browser-loaded runs (both
// runtimes) or via /api/triplets for a typed backend directory — identical
// payloads either way (source of truth: rmc_toolkits/triplets.py; port:
// workers/triplets.js).

import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import { isStaticMode, readAndParseLocalPlotFile } from '../browserData';
import { buildElementColors } from '../atomColors';
import { BOND_COLORS, PLOT_PALETTE } from '../plotPalette';
import {
    Banner, BondDash, Card, CardHeader, CardMeta, Chip, Control, ControlGroup, ControlsBar, ElementChip, Hint, Kpi,
    KpiRail, Page, PrimaryButton, Segmented, SegmentedButton, Switch, ToolButton, UnitField,
} from '../ui';
import InfoBadge from '../ui/InfoBadge';
import InteractivePlot from './InteractivePlot';
import ModelSummary from './ModelSummary';
import FoldedCellPanel from './FoldedCellPanel';
import useSiteCloud from '../useSiteCloud';
import { tripletRequestFromInputs } from '../workers/triplets';
import AppFooter from './AppFooter';
import './BondGeometryPage.css';

// Same cap as StructurePage/Dashboard: the Model information card needs the
// full counts, and parsing is worker-side anyway.
const STRUCTURE_MAX_POINTS = 1000000;

const DEGREES = '°';
const ANGSTROM = 'Å';

// The angle axis every angle plot uses (the ghost frame included, so nothing
// moves when the first result lands): 0–180° unpadded, labelled every 30°,
// minor marks every 10°, vertical grid at the labels, y from 0.
const ANGLE_AXIS = {
    xDomain: [0, 180],
    xTicks: [0, 30, 60, 90, 120, 150, 180],
    xMinorStep: 10,
    xGrid: true,
    yMin: 0
};

// The request fields a result depends on: a shown result whose fields differ
// from the current inputs is stale ("inputs changed", Compute reads Update).
const REQUEST_KEYS = ['end1', 'apex', 'end2', 'r12Min', 'r12Max', 'r23Min', 'r23Max', 'binWidth'];
const sameRequest = (one, two) => Boolean(one && two) && REQUEST_KEYS.every((key) => one[key] === two[key]);

// The input a validation error names (tripletRequestFromInputs starts its
// message with the field's label), so that field can be marked and focused.
const FIELD_LABELS = {
    r12Min: 'A–B window minimum',
    r12Max: 'A–B window maximum',
    r23Min: 'B–C window minimum',
    r23Max: 'B–C window maximum',
    binWidth: 'Bin width'
};
const fieldOfError = (message) => Object.keys(FIELD_LABELS)
    .find((key) => String(message ?? '').startsWith(FIELD_LABELS[key])) ?? null;

// A typed bound as a number; a blank box is NaN (no guide), never 0.
const typedNumber = (value) => (String(value ?? '').trim() === '' ? NaN : Number(value));

const formatNumber = (value, digits = 2) =>
    Number.isFinite(value) ? value.toFixed(digits) : '—';

// Realized bin width: 1 → "1.0", 0.5 → "0.5", 180/19 → "9.474".
const formatBinWidth = (width) => {
    const rounded = Math.round(width * 1000) / 1000;
    return Number.isInteger(rounded) ? rounded.toFixed(1) : String(rounded);
};

// Bounds at the precision they were given (2–4 decimals): a B–C window
// nudged to 3.4001 must not read as the A–B window's 3.40.
const formatBound = (value) => {
    if (!Number.isFinite(value)) return '—';
    const decimals = (String(value).split('.')[1] ?? '').length;
    return value.toFixed(Math.min(4, Math.max(2, decimals)));
};
const windowLabel = (window) => `${formatBound(window[0])}–${formatBound(window[1])} ${ANGSTROM}`;

// The isotropic reference in the density view: randomly oriented bonds put
// the fraction (cos θlo − cos θhi)/2 of their angles in each bin — the same
// exact bin integral the engine divides by for sin-corrected — so per degree
// it is (cos θlo − cos θhi)/(2·w). In the sin-corrected view it is exactly 1.
const isotropicDensity = (centers, width) => centers.map((center) => {
    const lo = ((center - width / 2) * Math.PI) / 180;
    const hi = ((center + width / 2) * Math.PI) / 180;
    return (Math.cos(lo) - Math.cos(hi)) / (2 * width);
});

const randomBondsGuide = (sin, centers, width) => (sin
    ? { label: 'random bonds', x: [0, 180], y: [1, 1], role: 'guide' }
    : { label: 'random bonds', x: centers, y: isotropicDensity(centers, width), role: 'guide', curve: 'step', binWidth: width });

// Element chips joined by bond dashes: "(Se)–(Ta)–(Se)" with the central
// atom ringed. Reads as plain text ("Se–Ta–Se") to a screen reader.
const TripletLabel = ({ elements, central, bonds, colors }) => (
    <span className="ui-element-chain">
        {elements.map((element, index) => (
            <React.Fragment key={`${element}-${index}`}>
                {index > 0 && <BondDash color={bonds[index - 1]} />}
                <ElementChip color={colors[element]} central={index === central}>{element}</ElementChip>
            </React.Fragment>
        ))}
    </span>
);

// The bonds a triplet draws, as a chain: one pair (B)–(A) when both bonds are
// the same type, else the whole triplet with the B–C dash in its own colour.
const bondChain = ([end1, apex, end2], shared) => (
    shared || end1 === end2
        ? { elements: [apex, end1], central: 0, bonds: [BOND_COLORS.ab] }
        : { elements: [end1, apex, end2], central: 1, bonds: [BOND_COLORS.ab, BOND_COLORS.bc] }
);

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

export default function BondGeometryPage({ directory, localRun, dataEpoch = 0 }) {
    const {
        sites,
        sitesError,
        loadingSites,
        requestPca,
        localFile,
        rmc6fText,
        ready,
        datasetKey
    } = useSiteCloud({ directory, localRun, dataEpoch });

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
    // String state: a cleared box stays empty and is an error, never 0.
    const [r12Min, setR12Min] = useState('2.00');
    const [r12Max, setR12Max] = useState('3.00');
    const [split23, setSplit23] = useState(false);
    const [r23Min, setR23Min] = useState('2.00');
    const [r23Max, setR23Max] = useState('3.00');
    const [binWidth, setBinWidth] = useState('1.0');

    const [result, setResult] = useState(null);
    // The request the shown result was computed from (stale check).
    const [resultRequest, setResultRequest] = useState(null);
    const [resultError, setResultError] = useState(null);
    // The input a validation error named: marked aria-invalid and focused.
    const [invalidField, setInvalidField] = useState(null);
    const [computing, setComputing] = useState(false);
    // Set when a new configuration of the SAME run replaced the one the
    // shown result was computed from (see below); cleared by Compute.
    const [configChanged, setConfigChanged] = useState(false);
    const [angleView, setAngleView] = useState('sin');

    const staticMode = isStaticMode();
    const noRun = staticMode && !localFile;
    // Where Compute runs: the shared worker for a browser-loaded run, the
    // Flask API for a typed backend directory (useSiteCloud's requestPca).
    const runtime = localFile ? 'browser' : 'server';

    const fieldRefs = {
        r12Min: useRef(null),
        r12Max: useRef(null),
        r23Min: useRef(null),
        r23Max: useRef(null),
        binWidth: useRef(null)
    };

    // A dataset switch clears any previous result immediately and bumps the
    // epoch, so a compute that was in flight for the old run can never land
    // its (stale) payload on the new one.
    const runEpoch = useRef(0);
    // Whether a computed angle distribution is on screen (read by the
    // configuration-change effect below without making it a dependency).
    const hasResult = useRef(false);
    // The configuration the page last saw: the Flask Live Data epoch and a
    // browser-loaded run's .rmc6f text (null while a new file is being read).
    const seenConfig = useRef({ epoch: dataEpoch, text: null });
    useEffect(() => {
        runEpoch.current += 1;
        seenConfig.current = { ...seenConfig.current, text: null };
        hasResult.current = false;
        setResult(null);
        setResultRequest(null);
        setResultError(null);
        // The in-flight compute (if any) can no longer land, so its
        // `finally` will not clear the busy state: clear it here.
        setComputing(false);
        setConfigChanged(false);
    }, [datasetKey]);
    // A new .rmc6f saved under the same run (Live Data) reloads the element
    // list, the Model information card and the partials in place, keeping
    // the triplet and the typed windows. The angle distribution is computed
    // on demand, so it is dropped rather than recomputed unasked: a result
    // from the previous configuration must never sit next to the new model.
    useEffect(() => {
        const seen = seenConfig.current;
        const epochChanged = dataEpoch !== seen.epoch;
        const textChanged = rmc6fText != null && seen.text != null && rmc6fText !== seen.text;
        seenConfig.current = { epoch: dataEpoch, text: rmc6fText ?? seen.text };
        if (!epochChanged && !textChanged) return;
        runEpoch.current += 1;
        if (hasResult.current) setConfigChanged(true);
        hasResult.current = false;
        setResult(null);
        setResultRequest(null);
        setResultError(null);
        setComputing(false);
    }, [dataEpoch, rmc6fText]);

    // The request the current inputs would send, or null while one is invalid.
    const currentRequest = useMemo(() => {
        try {
            return tripletRequestFromInputs({ end1, apex, end2, r12Min, r12Max, split23, r23Min, r23Max, binWidth });
        } catch {
            return null;
        }
    }, [end1, apex, end2, r12Min, r12Max, split23, r23Min, r23Max, binWidth]);
    // While a compute runs the badge says so; the stale cues wait for it.
    const stale = Boolean(result) && !computing && !sameRequest(resultRequest, currentRequest);

    const compute = useCallback(async () => {
        const epoch = runEpoch.current;
        setConfigChanged(false);
        // A cleared or non-numeric box is an error naming the field — never
        // sent as Number('') = 0, which silently widened the window to 0 Å.
        let params;
        try {
            params = tripletRequestFromInputs({
                end1, apex, end2, r12Min, r12Max, split23, r23Min, r23Max, binWidth
            });
        } catch (error) {
            hasResult.current = false;
            setResult(null);
            setResultRequest(null);
            setResultError(error.message);
            const field = fieldOfError(error.message);
            setInvalidField(field);
            fieldRefs[field]?.current?.focus();
            return;
        }
        setInvalidField(null);
        setComputing(true);
        setResultError(null);
        try {
            const data = await requestPca('triplets', params);
            if (runEpoch.current === epoch) {
                hasResult.current = true;
                setResult(data);
                setResultRequest(params);
            }
        } catch (error) {
            if (runEpoch.current === epoch) {
                hasResult.current = false;
                setResult(null);
                setResultRequest(null);
                setResultError(error.message);
            }
        } finally {
            if (runEpoch.current === epoch) setComputing(false);
        }
    // fieldRefs holds stable refs; listing it would re-create compute each render.
    // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [requestPca, end1, apex, end2, r12Min, r12Max, split23, r23Min, r23Max, binWidth]);

    const canCompute = !computing && ready && elements.length > 0 && Boolean(end1 && apex && end2);
    const submit = (event) => {
        event.preventDefault();
        if (canCompute) compute();
    };

    // Editing the field an error named clears the mark and the message.
    const edit = (key, setter) => (event) => {
        setter(event.target.value);
        if (invalidField === key) {
            setInvalidField(null);
            setResultError(null);
        }
    };

    // --- Structure for the Model information card. ---------------------------
    // Same source as the Dashboard and Atomic Density pages: a local run's
    // .rmc6f parses in the structure worker (both runtimes); a typed backend
    // directory asks /api/structure. ModelSummary renders the model card only
    // (showSymmetry={false}): the Detected SG card stays on the Dashboard and
    // Atomic Density pages.
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
    }, [localRun, directory, datasetKey, dataEpoch]);

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
    }, [localRun, directory, datasetKey, dataEpoch]);

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
    // No triplet fell inside the windows: the engine's curves are all zero,
    // which drawn under the random-bonds line would read as "far below
    // random". The hero keeps the empty axis and says so instead.
    const noAngles = Boolean(result) && result.angleCount === 0;
    // The angle distribution as a step curve (one step per bin) with a light
    // area, over the dashed isotropic reference.
    const anglePlot = useMemo(() => {
        if (!result || result.angleCount === 0) return null;
        const sin = angleView === 'sin';
        const width = result.binWidth;
        return {
            title: `${result.triplet.join('-')} bond angles`,
            xLabel: `angle at ${result.triplet[1]}, θ (${DEGREES})`,
            yLabel: sin ? 'sin-corrected (random = 1)' : 'density (deg^{-1})',
            ...ANGLE_AXIS,
            series: [
                {
                    label: sin ? 'sin-corrected' : 'density',
                    x: result.binCenters,
                    y: sin ? result.sinCorrected : result.density,
                    curve: 'step',
                    binWidth: width,
                    fill: true,
                    width: 1.75,
                    color: PLOT_PALETTE[0]
                },
                randomBondsGuide(sin, result.binCenters, width)
            ]
        };
    }, [result, angleView]);

    // Before the first result — or after one with no angles — the hero shows
    // the same axis, empty but for the random-bonds reference (in the density
    // view at the typed bin width, or the result's realized one).
    const ghostPlot = useMemo(() => {
        if (result && result.angleCount > 0) return null;
        const sin = angleView === 'sin';
        let width = result?.binWidth;
        let centers = result?.binCenters;
        if (!result) {
            const requested = Number(binWidth);
            const count = Number.isFinite(requested) && requested > 0 ? Math.max(1, Math.round(180 / requested)) : 180;
            width = 180 / count;
            centers = Array.from({ length: count }, (_, index) => (index + 0.5) * width);
        }
        return {
            title: result
                ? `${result.triplet.join('-')} bond angles (none in the windows)`
                : 'Bond-angle distribution (not computed)',
            xLabel: `angle at ${result ? result.triplet[1] : apex || 'B'}, θ (${DEGREES})`,
            yLabel: sin ? 'sin-corrected (random = 1)' : 'density (deg^{-1})',
            ...ANGLE_AXIS,
            series: [randomBondsGuide(sin, centers, width)]
        };
    }, [result, angleView, binWidth, apex]);

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
            ? Math.max(typedNumber(guideMax) || 0, typedNumber(guide23Max) || 0)
            : typedNumber(guideMax) || 0;
        const limit = Math.max(6, furthest * 2);
        const crop = (series, color) => {
            const cut = series.x.findIndex((value) => value > limit);
            const end = cut === -1 ? series.x.length : cut;
            return {
                label: series.label,
                x: series.x.slice(0, end),
                y: series.y.slice(0, end),
                color
            };
        };
        // One curve per distinct bond type; same-type triplets draw it once, or
        // the series label would be duplicated. Each curve wears its bond
        // role's colour (BOND_COLORS are the first two plot colours).
        const shells = [crop(partialSeries, BOND_COLORS.ab)];
        if (partialCurve23) shells.push(crop(partialCurve23, BOND_COLORS.bc));
        const max = shells.reduce(
            (acc, shell) => shell.y.reduce((inner, value) => Math.max(inner, value), acc),
            0
        ) || 1;
        // The windows are what the split governs: off, both bonds use the A–B
        // bounds and one pair of guides covers them. Labelled by role, not by
        // pair — with the same element at both ends the pair names would match.
        // Two windows also get the two bond-role colours: where the bonds are
        // different types each window is then the colour of the shell it
        // brackets. A lone window keeps the neutral guide stroke — there is
        // nothing to tell apart. Guides stay out of the legend: the window
        // chip in the card header names them.
        const windows = split23
            ? [
                { prefix: 'A–B ', min: guideMin, max: guideMax, color: BOND_COLORS.ab },
                { prefix: 'B–C ', min: guide23Min, max: guide23Max, color: BOND_COLORS.bc }
            ]
            : [{ prefix: '', min: guideMin, max: guideMax, color: undefined }];
        const guides = windows
            .flatMap((window) => [
                { label: `${window.prefix}rmin ${window.min}`, value: typedNumber(window.min), color: window.color },
                { label: `${window.prefix}rmax ${window.max}`, value: typedNumber(window.max), color: window.color }
            ])
            .filter((guide) => Number.isFinite(guide.value))
            .map((guide) => ({
                label: guide.label,
                x: [guide.value, guide.value],
                y: [0, max],
                role: 'guide',
                legend: false,
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

    // Coordination summary: mean bonds per B, the modal n and its share.
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

    // Bond tiles show physical bonds, each once (uniqueBonds). With the end
    // element equal to the central one every bond is found from both of its
    // ends, so the B-centred count (lengths.count, which the coordination
    // averages) is twice that; the tooltip says so.
    const bondsTitle = (lengths, end) => (
        end === result?.triplet[1]
            ? `Each ${end}–${end} bond counted once. Counted from the central atoms, as the `
                + `coordination is, there are ${lengths.count.toLocaleString()}: each bond is seen from both ends.`
            : `Each ${end}–${result?.triplet[1]} bond counted once, from its central ${result?.triplet[1]} atom.`
    );

    // Detected-bond overlay for the unit-cell panel: the computed windows in
    // their bond-role colours, B–C only when it is a bond of its own.
    const bondSets = useMemo(() => {
        if (!result) return null;
        const sets = [
            { elements: [result.triplet[0], result.triplet[1]], window: result.bond12, color: BOND_COLORS.ab }
        ];
        if (!result.sharedEnds) {
            sets.push({ elements: [result.triplet[1], result.triplet[2]], window: result.bond23, color: BOND_COLORS.bc });
        }
        return sets;
    }, [result]);

    // What the cards name: the computed triplet while a result is shown (the
    // plot and the 3D bonds are that result's), else the current picks.
    const hasTriplet = Boolean(end1 && apex && end2);
    const shownTriplet = result ? result.triplet : [end1, apex, end2];
    const [tripletA, tripletB, tripletC] = shownTriplet;
    const previewShared = end1 === end2
        && (!split23 || (Number(r12Min) === Number(r23Min) && Number(r12Max) === Number(r23Max)));
    const shownShared = result ? result.sharedEnds : previewShared;
    const tripletBonds = [BOND_COLORS.ab, shownShared ? BOND_COLORS.ab : BOND_COLORS.bc];
    const heroChain = hasTriplet || result
        ? <TripletLabel elements={shownTriplet} central={1} bonds={tripletBonds} colors={elementColors} />
        : null;
    const cellChain = hasTriplet || result ? bondChain(shownTriplet, shownShared) : null;
    // The g(r) card is a live helper: it follows the current picks.
    const pairChain = hasTriplet
        ? (partialCurve23
            ? { elements: [end1, apex, end2], central: 1, bonds: [BOND_COLORS.ab, BOND_COLORS.bc] }
            : { elements: [apex, end1], central: 0, bonds: [BOND_COLORS.ab] })
        : null;

    // Bond tiles: A–B always, B–C when it is a bond of its own. Same-pair
    // tiles (A = C with two windows) are told apart by role.
    const bondTiles = [
        { role: 'A–B', end: tripletA, lengths: result?.lengths12, window: result?.bond12, color: BOND_COLORS.ab }
    ];
    if (!shownShared) {
        bondTiles.push({ role: 'B–C', end: tripletC, lengths: result?.lengths23, window: result?.bond23, color: BOND_COLORS.bc });
    }
    const samePairTiles = bondTiles.length === 2 && tripletA === tripletC;

    // Before the first result: the central-atom count and the box, for the
    // prompt's chips.
    const apexAtoms = useMemo(() => (sites?.sites ?? [])
        .filter((site) => site.element === apex)
        .reduce((acc, site) => acc + (site.count ?? 0), 0), [sites, apex]);
    const box = sites?.supercell ? sites.supercell.join('×') : null;
    const typed = (value) => (String(value).trim() === '' ? '…' : value);
    const windowText = split23
        ? `${typed(r12Min)}–${typed(r12Max)} / ${typed(r23Min)}–${typed(r23Max)} ${ANGSTROM}`
        : `${typed(r12Min)}–${typed(r12Max)} ${ANGSTROM}`;
    // The windows a result used (resolved payload values), B–C when it differs.
    const resultWindowText = result
        ? (result.sharedEnds || windowLabel(result.bond23) === windowLabel(result.bond12)
            ? windowLabel(result.bond12)
            : `${windowLabel(result.bond12).replace(` ${ANGSTROM}`, '')} / ${windowLabel(result.bond23)}`)
        : null;

    const noHelper = !partialSeries && !partialsLoading;
    // The helper plot's Save (and Reset zoom) live in its card header: the
    // title's chips already name the curves, so the plot drops its legend
    // row and keeps that height for the g(r) (InteractivePlot actionsTarget).
    const [helperActions, setHelperActions] = useState(null);
    const ghost = !result;
    const dimmed = ghost || noAngles || computing;
    const showPrompt = ghost && !computing && elements.length > 0 && hasTriplet && !noRun;
    const showNoAngles = noAngles && !computing;

    const runButtonProps = {
        run: true,
        busy: computing,
        stale,
        disabled: !canCompute,
        'aria-busy': computing || undefined
    };

    return (
        <Page as="div" column>
            {/* Model information first, as on the Dashboard; the Detected SG
                card stays on the Dashboard/Atomic Density pages only. */}
            {structure && (
                <div className="geom-model">
                    <ModelSummary structure={structure} showSymmetry={false} />
                </div>
            )}
            {/* A form, so Enter in any field computes; validation stays the
                page's own (noValidate: a cleared box is a named error). */}
            <ControlsBar as="form" className="geom-controls" onSubmit={submit} noValidate aria-label="Triplet and bond windows">
                <ControlGroup label="Triplet">
                    <Control
                        label={(
                            <>
                                A
                                <InfoBadge label="About the triplet">
                                    <p>
                                        The RMCProfile <code>triplets</code> convention: the angle is
                                        measured at the <strong>central atom B</strong> between its bond
                                        to an A atom and its bond to a C atom. A and C may name the same
                                        element (each unordered pair of bonds then counts once).
                                    </p>
                                    <p>
                                        Angles are counted over the periodic configuration exactly,
                                        periodic images included.
                                    </p>
                                </InfoBadge>
                            </>
                        )}
                    >
                        <i className="ui-color-dot" style={{ background: elementColors[end1] }} aria-hidden="true" />
                        <select className="ui-select" value={end1} onChange={(event) => setEnd1(event.target.value)} disabled={!elements.length} aria-label="End element A">
                            {elements.map((element) => <option key={element} value={element}>{element}</option>)}
                        </select>
                    </Control>
                    <Control label="B (central)">
                        <i className="ui-color-dot" style={{ background: elementColors[apex] }} aria-hidden="true" />
                        <select className="ui-select ui-select--ring" value={apex} onChange={(event) => setApex(event.target.value)} disabled={!elements.length} aria-label="Central element B">
                            {elements.map((element) => <option key={element} value={element}>{element}</option>)}
                        </select>
                    </Control>
                    <Control label="C">
                        <i className="ui-color-dot" style={{ background: elementColors[end2] }} aria-hidden="true" />
                        <select className="ui-select" value={end2} onChange={(event) => setEnd2(event.target.value)} disabled={!elements.length} aria-label="End element C">
                            {elements.map((element) => <option key={element} value={element}>{element}</option>)}
                        </select>
                    </Control>
                    {end1 && end2 && end1 !== end2 && (
                        <ToolButton
                            aria-label="Swap A and C"
                            title="Swap A and C"
                            onClick={() => { setEnd1(end2); setEnd2(end1); }}
                        >
                            ⇄
                        </ToolButton>
                    )}
                </ControlGroup>

                <ControlGroup label="Bond windows">
                    <Control
                        label={(
                            <>
                                <i className="ui-role-bar" style={{ '--role': BOND_COLORS.ab }} aria-hidden="true" />
                                A{'–'}B
                                <InfoBadge label="About the bond windows">
                                    <p>
                                        Two atoms are bonded when their distance falls inside the window
                                        (inclusive). Read the window off the first-shell peak of the
                                        partial g(r) — the partial g(r) card marks the current bounds
                                        with dashed guides once a partials file is in the run folder.
                                    </p>
                                    <p>
                                        The angle count grows roughly as r_max⁶, so a wide window can
                                        exceed the per-request limit.
                                    </p>
                                </InfoBadge>
                            </>
                        )}
                    >
                        <span className="ui-pair">
                            <UnitField
                                unit={ANGSTROM}
                                invalid={invalidField === 'r12Min'}
                                inputProps={{
                                    ref: fieldRefs.r12Min, type: 'number', inputMode: 'decimal', step: '0.05', min: '0', max: '15',
                                    value: r12Min, onChange: edit('r12Min', setR12Min), 'aria-label': 'A-B window minimum',
                                    'aria-invalid': invalidField === 'r12Min' || undefined
                                }}
                            />
                            {'–'}
                            <UnitField
                                unit={ANGSTROM}
                                invalid={invalidField === 'r12Max'}
                                inputProps={{
                                    ref: fieldRefs.r12Max, type: 'number', inputMode: 'decimal', step: '0.05', min: '0', max: '15',
                                    value: r12Max, onChange: edit('r12Max', setR12Max), 'aria-label': 'A-B window maximum',
                                    'aria-invalid': invalidField === 'r12Max' || undefined
                                }}
                            />
                        </span>
                    </Control>
                    <Switch
                        label={<>Distinct B{'–'}C</>}
                        checked={split23}
                        onChange={(event) => setSplit23(event.target.checked)}
                        inputProps={{ 'aria-label': 'Use a distinct B-C window' }}
                    />
                    {split23 && (
                        <Control
                            className="ui-reveal"
                            label={(
                                <>
                                    <i className="ui-role-bar" style={{ '--role': BOND_COLORS.bc }} aria-hidden="true" />
                                    B{'–'}C
                                </>
                            )}
                        >
                            <span className="ui-pair">
                                <UnitField
                                    unit={ANGSTROM}
                                    invalid={invalidField === 'r23Min'}
                                    inputProps={{
                                        ref: fieldRefs.r23Min, type: 'number', inputMode: 'decimal', step: '0.05', min: '0', max: '15',
                                        value: r23Min, onChange: edit('r23Min', setR23Min), 'aria-label': 'B-C window minimum',
                                        'aria-invalid': invalidField === 'r23Min' || undefined
                                    }}
                                />
                                {'–'}
                                <UnitField
                                    unit={ANGSTROM}
                                    invalid={invalidField === 'r23Max'}
                                    inputProps={{
                                        ref: fieldRefs.r23Max, type: 'number', inputMode: 'decimal', step: '0.05', min: '0', max: '15',
                                        value: r23Max, onChange: edit('r23Max', setR23Max), 'aria-label': 'B-C window maximum',
                                        'aria-invalid': invalidField === 'r23Max' || undefined
                                    }}
                                />
                            </span>
                        </Control>
                    )}
                </ControlGroup>

                <ControlGroup label="Histogram">
                    <Control label="Bin width">
                        <UnitField
                            unit={DEGREES}
                            invalid={invalidField === 'binWidth'}
                            inputProps={{
                                ref: fieldRefs.binWidth, type: 'number', inputMode: 'decimal', step: '0.5', min: '0.1', max: '45',
                                value: binWidth, onChange: edit('binWidth', setBinWidth), 'aria-label': 'Angle bin width in degrees',
                                'aria-invalid': invalidField === 'binWidth' || undefined
                            }}
                        />
                    </Control>
                </ControlGroup>

                <div className="ui-cluster ui-cluster--end">
                    <kbd className="ui-kbd" title="Enter in any field computes">Enter</kbd>
                    <PrimaryButton type="submit" title="Or press Enter in any field" {...runButtonProps}>
                        {stale ? 'Update' : 'Compute'}
                    </PrimaryButton>
                </div>
            </ControlsBar>

            {noRun && <Hint>Open a run folder with an <code>.rmc6f</code> file.</Hint>}
            {sitesError && <Banner as="p" tone="danger" sm>{sitesError}</Banner>}

            <div className={noHelper ? 'geom-layout geom-layout--no-helper' : 'geom-layout'}>
                <Card roundEnds className="geom-hero">
                    <CardHeader
                        wrap
                        title={<>{heroChain}{heroChain && ' '}<span>{heroChain ? 'bond angles' : 'Bond-angle distribution'}</span></>}
                        help={(
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
                                    and a peak near 180{DEGREES} is no longer flattened.
                                </p>
                                <p>
                                    The dashed <b>random bonds</b> line is that reference in either
                                    view: 1 when sin-corrected, the exact isotropic fraction of each
                                    bin per degree in density.
                                </p>
                                <p>
                                    Same shape as the <code>norm/sin(theta)</code> column of
                                    RMCProfile's <code>triplets</code>, on another scale: that
                                    column is this curve × sin(w/2)/w for w-degree bins
                                    (≈ π/360 ≈ 0.00873), so compare shapes, or rescale.
                                </p>
                            </InfoBadge>
                        )}
                        actions={(
                            <>
                                {stale && (
                                    <Chip tone="warn" title="The plot is from other inputs: Update to recompute.">
                                        inputs changed
                                    </Chip>
                                )}
                                {configChanged && !result && (
                                    <Chip tone="warn" title="A new configuration of this run was saved.">
                                        new configuration
                                    </Chip>
                                )}
                                <Segmented role="group" aria-label="Angle normalization">
                                    <SegmentedButton type="button" active={angleView === 'sin'} onClick={() => setAngleView('sin')}>
                                        sin-corrected
                                    </SegmentedButton>
                                    <SegmentedButton type="button" active={angleView === 'raw'} onClick={() => setAngleView('raw')}>
                                        density
                                    </SegmentedButton>
                                </Segmented>
                            </>
                        )}
                    />
                    {/* The headline results, in place before Compute ("—"), so
                        nothing below moves when they arrive. */}
                    <KpiRail role="status" aria-live="polite" aria-label="Triplet result">
                        {/* The mean ± std angle is on the sub line (visible to
                            keyboard and touch users), not the headline: a
                            multimodal mean is not a bond angle. Compact, so it
                            fits a four-tile rail at 1280 px: "281,846 · mean
                            89.9 ± 42.8°"; the hover spells it out. */}
                        <Kpi
                            label="Angles"
                            value={result ? formatNumber(result.apexCount ? result.angleCount / result.apexCount : 0, 1) : null}
                            unit={`per ${tripletB}`}
                            sub={result && (result.meanAngle != null
                                ? `${result.angleCount.toLocaleString()} · mean ${formatNumber(result.meanAngle, 1)} ± ${formatNumber(result.stdAngle, 1)}${DEGREES}`
                                : `${result.angleCount.toLocaleString()} angles`)}
                            title={result
                                ? `${result.angleCount.toLocaleString()} angles`
                                    + (result.meanAngle != null
                                        ? `, mean ${formatNumber(result.meanAngle, 1)}${DEGREES} ± ${formatNumber(result.stdAngle, 1)}${DEGREES} (std)`
                                        : '')
                                    + `; ${formatBinWidth(result.binWidth)}${DEGREES} bins; computed in the ${runtime}.`
                                : undefined}
                        />
                        <Kpi
                            label="Coordination"
                            value={coordinationSummary ? formatNumber(coordinationSummary.mean) : null}
                            unit={`per ${tripletB}`}
                            sub={coordinationSummary
                                && `${coordinationSummary.mode}-fold ${formatNumber(coordinationSummary.modeShare, 1)}% · of ${result.apexCount.toLocaleString()}`}
                        />
                        {bondTiles.map((tile) => {
                            const mean = tile.lengths?.meanLength;
                            return (
                                <Kpi
                                    key={tile.role}
                                    label={(
                                        <>
                                            <BondDash color={tile.color} lead text="" />
                                            {tripletB && tile.end
                                                ? `${tripletB}–${tile.end} bond${samePairTiles ? ` (${tile.role})` : ''}`
                                                : 'Bond'}
                                        </>
                                    )}
                                    value={result ? (mean != null ? formatNumber(mean, 3) : '—') : null}
                                    unit={mean != null ? ANGSTROM : undefined}
                                    sub={result && tile.lengths
                                        && `${tile.lengths.uniqueBonds.toLocaleString()} bonds · ${windowLabel(tile.window)}`}
                                    title={result && tile.lengths ? bondsTitle(tile.lengths, tile.end) : undefined}
                                />
                            );
                        })}
                    </KpiRail>
                    <div className="geom-plot">
                        <div className={dimmed ? 'geom-plot__frame ui-dim' : 'geom-plot__frame'} inert={dimmed || undefined}>
                            <InteractivePlot
                                file={{ path: `geometry:angles:${datasetKey}`, name: (anglePlot ?? ghostPlot).title }}
                                plotData={anglePlot ?? ghostPlot}
                                variant="fit"
                            />
                        </div>
                        {computing && !result && <div className="ui-shimmer" aria-hidden="true" />}
                        {computing && (
                            <span className="ui-overlay-badge" role="status">Computing {end1}–{apex}–{end2}…</span>
                        )}
                        {showPrompt && (
                            <div className="ui-overlay-center">
                                <div className="ui-prompt">
                                    <div className="ui-prompt__row">
                                        {heroChain}
                                        <Chip title="Bond window">{windowText}</Chip>
                                    </div>
                                    {resultError
                                        ? <p className="ui-prompt__error" role="alert">{resultError}</p>
                                        : <p>{configChanged ? 'New configuration — Compute again.' : 'Pick a triplet, then Compute.'}</p>}
                                    <PrimaryButton type="button" onClick={compute} {...runButtonProps}>
                                        Compute
                                    </PrimaryButton>
                                    <div className="ui-prompt__row">
                                        {apexAtoms > 0 && (
                                            <Chip
                                                title={`The angle at ${apex} between two bonds, over all ${apexAtoms.toLocaleString()} ${apex}`
                                                    + `${box ? ` in the ${box} box` : ''}; periodic images included, exact.`}
                                            >
                                                {apexAtoms.toLocaleString()} {apex}
                                            </Chip>
                                        )}
                                        {box && <Chip title="Supercell">{box}</Chip>}
                                        <Chip title={runtime === 'browser' ? 'Runs in your browser' : 'Runs on the server'}>
                                            {runtime}
                                        </Chip>
                                    </div>
                                </div>
                            </div>
                        )}
                        {showNoAngles && (
                            <div className="ui-overlay-center">
                                <div className="ui-prompt">
                                    <div className="ui-prompt__row">
                                        {heroChain}
                                        <Chip title="The bond windows this result used">{resultWindowText}</Chip>
                                    </div>
                                    <p>No {result.triplet.join('–')} triplets in these windows.</p>
                                </div>
                            </div>
                        )}
                    </div>
                </Card>

                {/* The Atomic Density page's folded cell — the atom cloud, not
                    fitted ellipsoids — with the detected bonds over it. */}
                <FoldedCellPanel
                    className="geom-cell"
                    title={cellChain
                        ? <><TripletLabel {...cellChain} colors={elementColors} /> <span>bonds</span></>
                        : 'Folded unit cell'}
                    structure={structure}
                    sites={sites}
                    elementColors={elementColors}
                    loading={loadingSites}
                    bondSets={bondSets}
                    legendEmphasis={hasTriplet || result ? shownTriplet : null}
                    legendExtras={bondSets
                        ? bondSets.map((set) => (
                            <span key={set.color} className="ui-legend__item">
                                <BondDash color={set.color} text="" />
                                {`${tripletB}–${set.elements[1] === tripletB ? set.elements[0] : set.elements[1]} ${windowLabel(set.window)}`}
                            </span>
                        ))
                        : elements.length > 0 && (
                            <span className="ui-legend__item is-muted">Compute to draw bonds</span>
                        )}
                />

                <Card roundEnds className="geom-helper">
                    <CardHeader
                        wrap
                        title={pairChain
                            ? <><TripletLabel {...pairChain} colors={elementColors} /> <span>partial g(r)</span></>
                            : 'Partial g(r)'}
                        help={(
                            <InfoBadge label="About the partial g(r)" side="above">
                                <p>
                                    The A{'–'}B partial pair distribution from the run's{' '}
                                    <code>PDFpartials.csv</code>. Set the bond window to bracket
                                    the first-shell peak; the dashed guides track the current
                                    bounds.
                                </p>
                                <p>
                                    A second partial is drawn whenever A{'–'}B and B{'–'}C are
                                    different pair types (Ga{'–'}Nb{'–'}Se: Ga{'–'}Nb and Nb{'–'}Se),
                                    whether or not <b>Distinct B{'–'}C</b> is on; a same-type
                                    triplet such as Se{'–'}Nb{'–'}Se has one shell and one curve.
                                    The switch governs the guides: off, one neutral pair covers
                                    both bonds; on, each window gets its own labelled pair
                                    (A{'–'}B, B{'–'}C), colored like its shell's curve when the
                                    two bonds are different types.
                                </p>
                            </InfoBadge>
                        )}
                        actions={noHelper
                            ? (
                                <CardMeta>
                                    {noRun
                                        ? '—'
                                        : partials
                                            ? <>No {end1}{'–'}{apex} partial in <code>PDFpartials.csv</code>.</>
                                            : <>No <code>PDFpartials.csv</code> in this run.</>}
                                </CardMeta>
                            )
                            : (
                                <>
                                    {/* Split: one chip for both windows, each led by its
                                        bond-role dash (the guides' colours), so the
                                        header keeps one row with the Save beside it. */}
                                    {split23
                                        ? (
                                            <Chip title="A–B window (the blue guides), B–C window (the orange guides)">
                                                <BondDash color={BOND_COLORS.ab} lead text="A–B " />
                                                {typed(guideMin)}{'–'}{typed(guideMax)}{'\u2003'}
                                                <BondDash color={BOND_COLORS.bc} lead text="B–C " />
                                                {typed(guide23Min)}{'–'}{typed(guide23Max)} {ANGSTROM}
                                            </Chip>
                                        )
                                        : (
                                            <Chip title="Bond window (the dashed guides)">
                                                {typed(guideMin)}{'–'}{typed(guideMax)} {ANGSTROM}
                                            </Chip>
                                        )}
                                    <span className="ui-card__cluster" ref={setHelperActions} />
                                </>
                            )}
                    />
                    {!noHelper && (
                        <div className="geom-plot">
                            {helperPlot
                                ? (
                                    <InteractivePlot
                                        file={{ path: `geometry:partial:${datasetKey}`, name: 'partial-gr' }}
                                        plotData={helperPlot}
                                        variant="fit"
                                        legend={false}
                                        actionsTarget={helperActions}
                                    />
                                )
                                : <p className="ui-placeholder">Looking for a partials file…</p>}
                        </div>
                    )}
                </Card>
            </div>

            <AppFooter tight />
        </Page>
    );
}
