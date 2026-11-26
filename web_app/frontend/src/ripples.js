// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Fourier termination ripples: a line-for-line port of rmc_toolkits/ripples.py
// (the source of truth; its docstrings hold the reasoning, this file the
// arithmetic). Per fitted dataset of a run, from its files and first
// principles: the fit space and how RMCProfile makes the model it compares
// (whether the box is asked to fit ripples), the transform (Q_min, Q_max, ΔQ,
// window, where they come from) and the r-resolution (the kernel's FWHM and
// main lobe, the Nyquist step, the ripple period). Whether a given peak is a
// ripple is not judged (that reading is on the branch wip/joint-fit). Parity
// golden ripples_fixture.json (tests/generate_ripples_fixture.py →
// src/__tests__/ripplesParity.test.js). Pure: no I/O (ripplesRun.js reads a
// run's files).

import { FT_FILE_RE, RMC_FT_TOLERANCE, reproduceTransform, uniformStep } from './ftUncertainty';

// --- rules (rmc_toolkits/ripples.py) -----------------------------------------------

export const DEFAULT_DQ = 0.01;
export const INPUT_MAP_TOLERANCE = 1e-5;

export const RISKS = Object.freeze(['high', 'matched', 'unknown', 'none']);

/** Python's round(): halves to the even integer. */
const pyRound = (x) => {
    const floor = Math.floor(x);
    const fraction = x - floor;
    if (fraction > 0.5) return floor + 1;
    if (fraction < 0.5) return floor;
    return floor % 2 === 0 ? floor : floor + 1;
};

const sum = (values) => {
    let total = 0;
    for (let i = 0; i < values.length; i += 1) total += values[i];
    return total;
};
const dot = (a, b) => {
    let total = 0;
    for (let i = 0; i < a.length; i += 1) total += a[i] * b[i];
    return total;
};
const maxAbs = (values) => {
    let best = 0;
    for (let i = 0; i < values.length; i += 1) best = Math.max(best, Math.abs(values[i]));
    return best;
};
const pick = (values, mask) => {
    const out = [];
    for (let i = 0; i < values.length; i += 1) if (mask[i]) out.push(values[i]);
    return Float64Array.from(out);
};

// --- the transform's grid and kernel ---------------------------------------------------

/** Lorch's window sin(πQ/Q_max)/(πQ/Q_max). */
export const lorch = (q, qmax) => Float64Array.from(q, (value) => {
    const x = Math.PI * value / qmax;
    return x !== 0 ? Math.sin(x) / x : 1;
});

/** The transform's Q points q_min + k·ΔQ, k = 0 … n − 1. */
export const qGrid = (qMin, dq, n) => Float64Array.from({ length: n }, (_, k) => qMin + dq * k);

/** w_k: ΔQ (rectangle) or the trapezoid's, times the window. */
export const transformWeights = (q, dq, rule, window) => {
    let weights = new Float64Array(q.length).fill(dq);
    if (rule === 'trapezoid' && q.length > 1) {
        weights[0] *= 0.5;
        weights[weights.length - 1] *= 0.5;
    }
    if (window === 'lorch') {
        const values = lorch(q, q[q.length - 1]);
        weights = Float64Array.from(weights, (value, k) => value * values[k]);
    }
    return weights;
};

/** K(x) = Σ_k w_k cos(Q_k x). */
export const kernel = (q, weights, x) => {
    let total = 0;
    for (let k = 0; k < q.length; k += 1) total += weights[k] * Math.cos(q[k] * x);
    return total;
};

const firstCrossing = (q, weights, level) => {
    const qmax = q[q.length - 1];
    const step = 0.05 / qmax;
    let lo = 0;
    let hi = step;
    while (kernel(q, weights, hi) > level && hi < 100 / qmax) {
        lo = hi;
        hi += step;
    }
    for (let i = 0; i < 60; i += 1) {
        const mid = 0.5 * (lo + hi);
        if (kernel(q, weights, mid) > level) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
};

/** The Nyquist step, the kernel's FWHM and main lobe, the ripple period. */
export const resolution = (q, weights) => {
    const qmax = q[q.length - 1];
    return {
        nyquistStep: Math.PI / qmax,
        fwhm: 2 * firstCrossing(q, weights, 0.5 * kernel(q, weights, 0)),
        mainLobe: firstCrossing(q, weights, 0),
        ripplePeriod: 2 * Math.PI / qmax
    };
};

// --- the fit space, from the run-control file -----------------------------------------------

const DAT_KEY_RE = /^([A-Z_\d]+)\s*::\s*(.*)$/;
const DAT_SUB_RE = /^([A-Z_]+)\s*(?:::\s*(.*))?$/;
const LINE_BREAK_RE = /\r\n|\r|\n/;
const NUMBER_RE = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/;

const numbers = (value) => {
    const out = [];
    for (const token of value.split(/\s+/).filter(Boolean)) {
        if (!NUMBER_RE.test(token)) break;
        out.push(Number(token));
    }
    return out;
};
const firstNumber = (value) => {
    const tokens = value.split(/\s+/).filter(Boolean);
    if (!tokens.length) return null;
    return NUMBER_RE.test(tokens[0]) ? Number(tokens[0]) : null;
};

/** What the .dat says about each dataset's transform (ripple_settings). */
export const rippleSettings = (text) => {
    const datasets = [];
    let block = null;
    let minimum = null;
    for (const raw of text.split(LINE_BREAK_RE)) {
        const line = raw.trim();
        if (!line) continue;
        if (line.startsWith('>')) {
            if (block === null) continue;
            const match = DAT_SUB_RE.exec(line.slice(1).trim());
            if (!match) continue;
            const key = match[1];
            const value = (match[2] || '').trim();
            if (key === 'FILENAME') block.file = value;
            else if (key === 'FIT_TYPE') block.fitType = value;
            else if (key === 'DATA_TYPE') block.dataType = value;
            else if (key === 'REAL_SPACE_FIT') block.realSpaceFit = true;
            else if (key === 'CONVOLVE') { block.convolve = true; block.convolveQ = firstNumber(value); }
            else if (key === 'R_CUTOFF') block.rCutoff = firstNumber(value);
            else if (key === 'QMIN') block.qMin = firstNumber(value);
            else if (key === 'QMAX') block.qMax = firstNumber(value);
            else if (key === 'QSTEP') block.qStep = firstNumber(value);
            else if (key === 'NEUTRON_COEFFICIENTS') block.coefficients = numbers(value);
            continue;
        }
        const match = DAT_KEY_RE.exec(line);
        block = null;
        if (!match) continue;
        const [, key, value] = match;
        if (key === 'END') break;
        if (key.endsWith('_DATA')) {
            block = {
                block: key, file: null, fitType: null, dataType: null, realSpaceFit: false, convolve: false,
                convolveQ: null, rCutoff: null, qMin: null, qMax: null, qStep: null, coefficients: null
            };
            datasets.push(block);
        } else if (key === 'MINIMUM_DISTANCES') {
            minimum = numbers(value);
        }
    }
    return { minimumDistances: minimum, datasets };
};

const nthBlock = (settings, block, number) => {
    if (!settings) return null;
    const blocks = settings.datasets.filter((entry) => entry.block === block);
    return number > 0 && number <= blocks.length ? blocks[number - 1] : null;
};

/** How RMCProfile makes the model this dataset is compared with, and the risk (fit_space_risk). */
export const fitSpaceRisk = (kind, dataset, dataTransform) => {
    if (kind === 'reciprocal') {
        return { fitSpace: 'reciprocal', model: null, modelQ: null, rCutoff: null, risk: 'none', because: 'fitted-in-q' };
    }
    if (kind === 'xpdf') {
        return { fitSpace: 'real', model: 'same-transform', modelQ: null, rCutoff: null, risk: 'matched', because: 'same-transform' };
    }
    const known = dataTransform !== null && dataTransform !== undefined && ['provenance', 'stog'].includes(dataTransform.source);
    if (kind === 'npdf') {
        if (!dataset || !dataset.convolve) {
            return { fitSpace: 'real', model: 'histogram', modelQ: null, rCutoff: null, risk: 'high', because: 'untruncated-model' };
        }
        const modelQ = dataset.convolveQ ?? null;
        const out = { fitSpace: 'real', model: 'sinc', modelQ, rCutoff: dataset.rCutoff ?? null };
        if (modelQ === null) return { ...out, risk: 'unknown', because: 'sinc-no-q' };
        if (!known) return { ...out, risk: 'unknown', because: 'sinc-data-unverified' };
        if (dataTransform.window === 'lorch') return { ...out, risk: 'high', because: 'sinc-vs-lorch' };
        if (Math.abs(modelQ - dataTransform.qMax) > 0.02 * dataTransform.qMax) return { ...out, risk: 'high', because: 'sinc-q-differs' };
        return { ...out, risk: 'matched', because: 'sinc-matches' };
    }
    if (kind === 'xray-real') {
        const modelQ = dataset ? dataset.qMax ?? null : null;
        const out = { fitSpace: 'real', model: 'xray-transform', modelQ, rCutoff: null };
        if (modelQ === null || !known) return { ...out, risk: 'unknown', because: 'xray-real-unverified' };
        if (dataTransform.window === 'lorch') return { ...out, risk: 'high', because: 'transform-vs-lorch' };
        if (Math.abs(modelQ - dataTransform.qMax) > 0.02 * dataTransform.qMax) return { ...out, risk: 'high', because: 'transform-q-differs' };
        return { ...out, risk: 'matched', because: 'transform-matches' };
    }
    return { fitSpace: null, model: null, modelQ: null, rCutoff: null, risk: 'unknown', because: 'unknown-kind' };
};

const baseName = (path) => String(path).replace(/\\/g, '/').split('/').pop();

/** An Auto StoG provenance's transform, if one of its outputs is `dataName` (provenance_transform). */
export const provenanceTransform = (payload, dataName) => {
    if (!payload || typeof payload !== 'object' || payload.tool !== 'rmc-autoscale') return null;
    const names = new Set(Object.values(payload.outputs || {}).map((path) => baseName(path).toLowerCase()));
    if (!names.has(dataName.toLowerCase())) return null;
    const config = (payload.provenance || {}).config || {};
    const qMin = Number(config.qmin);
    const qMax = Number(config.qmax);
    if (config.qmin === undefined || config.qmax === undefined || config.qmin === null || config.qmax === null
        || !Number.isFinite(qMin) || !Number.isFinite(qMax)) return null;
    const nPoints = (payload.provenance || {}).n_q_points;
    const dq = typeof nPoints === 'number' && nPoints > 1 ? (qMax - qMin) / (nPoints - 1) : null;
    const enforcement = payload.enforcement || null;
    return {
        source: 'provenance',
        qMin,
        qMax,
        dQ: dq,
        window: config.lorch ? 'lorch' : 'none',
        filterCutoff: config.r_cutoff ? Number(config.r_cutoff) : null,
        enforceCutoff: enforcement && typeof enforcement === 'object' && enforcement.cutoff !== undefined && enforcement.cutoff !== null
            ? Number(enforcement.cutoff) : null
    };
};

/** A classic stog input's transform (autoScale.readStogInp), if one of its outputs is `dataName`. */
export const stogTransform = (inp, dataName) => {
    const names = new Set([inp.outRmcGr, inp.outRmcDr, inp.outRmcFq, inp.outFtGr, inp.outGr].map((name) => baseName(name).toLowerCase()));
    if (!names.has(dataName.toLowerCase())) return null;
    return {
        source: 'stog',
        qMin: inp.qmin,
        qMax: inp.qmax,
        dQ: null,
        window: inp.lorch ? 'lorch' : 'none',
        // readStogInp refuses a filter-disabled input, as read_stog_inp does.
        filterCutoff: inp.rCutoff > 0 ? inp.rCutoff : null,
        enforceCutoff: inp.peakCutoff > 0 ? inp.peakCutoff : null
    };
};

/** Whether a data file's title line names rmc-autoscale. */
export const autoStogTitle = (text) => text.split(LINE_BREAK_RE).slice(0, 3).some((line) => line.trim().startsWith('rmc-autoscale'));

// --- one dataset of a run -------------------------------------------------------------------

const XRAY_PDF_RE = /^(.+)_Xray_PDF(\d+)\.csv$/;
const PDF_RE = /^(.+)_PDF(\d+)\.csv$/;
const FQ_RE = /^(.+)_FQ(\d+)\.csv$/;
const SQ_RE = /^(.+)_SQ(\d+)\.csv$/;

export const BLOCKS = Object.freeze({
    xpdf: 'XRAY_RECIPROCAL_SPACE_DATA',
    npdf: 'NEUTRON_REAL_SPACE_DATA',
    'xray-real': 'XRAY_REAL_SPACE_DATA',
    'xray-reciprocal': 'XRAY_RECIPROCAL_SPACE_DATA',
    'neutron-reciprocal': 'NEUTRON_RECIPROCAL_SPACE_DATA'
});

/** What a run output file is: { kind, stem, number }, or null (dataset_file). */
export const datasetFile = (name) => {
    for (const [pattern, kind] of [[FT_FILE_RE, 'xpdf'], [XRAY_PDF_RE, 'xray-real'], [PDF_RE, 'npdf'],
        [FQ_RE, 'xray-reciprocal'], [SQ_RE, 'neutron-reciprocal']]) {
        const match = pattern.exec(name);
        if (match) return { kind, stem: match[1], number: Number(match[2]) };
    }
    return null;
};


const uniform = (r) => {
    if (r.length < 3) return false;
    const step = (r[r.length - 1] - r[0]) / (r.length - 1);
    if (!(step > 0)) return false;
    let worst = 0;
    for (let i = 0; i < r.length; i += 1) worst = Math.max(worst, Math.abs(r[i] - (r[0] + step * i)));
    return worst <= 1e-6 * step;
};

const searchLeft = (values, at) => {
    let lo = 0;
    let hi = values.length;
    while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (values[mid] < at) lo = mid + 1; else hi = mid;
    }
    return lo;
};

/** The input data file in the fitted file's units (map_input). */
export const mapInput = (rFile, yFile, columns) => {
    const rAll = columns[0];
    const yAll = columns[1];
    const finite = Array.from(rAll, (value, i) => Number.isFinite(value) && Number.isFinite(yAll[i]));
    const rIn = pick(rAll, finite);
    const yIn = pick(yAll, finite);
    if (rIn.length < 10) return null;
    const step = rFile.length > 1 ? rFile[1] - rFile[0] : 0;
    const matchedFile = [];
    const matchedInput = [];
    for (let k = 0; k < rFile.length; k += 1) {
        const i = searchLeft(rIn, rFile[k]);
        for (const j of [i - 1, i]) {
            if (j >= 0 && j < rIn.length && Math.abs(rIn[j] - rFile[k]) <= 1e-6 * Math.max(step, 1e-12)) {
                matchedFile.push(k);
                matchedInput.push(j);
                break;
            }
        }
    }
    if (matchedFile.length < 10) return null;
    const x = Float64Array.from(matchedInput, (j) => yIn[j]);
    const target = Float64Array.from(matchedFile, (k) => yFile[k]);
    const n = x.length;
    const sx = sum(x);
    const sy = sum(target);
    const sxx = dot(x, x);
    const sxy = dot(x, target);
    const det = n * sxx - sx * sx;
    if (det === 0) return null;
    const scale = (n * sxy - sx * sy) / det;
    const offset = (sxx * sy - sx * sxy) / det;
    const largest = maxAbs(target);
    let worst = 0;
    for (let i = 0; i < n; i += 1) worst = Math.max(worst, Math.abs(target[i] - scale * x[i] - offset));
    const residual = largest > 0 ? worst / largest : Infinity;
    if (!(residual <= INPUT_MAP_TOLERANCE) || !(scale > 0)) return null;
    const keep = Array.from(rIn, (value) => value > 0 && value <= rFile[rFile.length - 1] + 1e-9);
    return {
        r: pick(rIn, keep), y: Float64Array.from(pick(yIn, keep), (value) => scale * value + offset), scale, offset, residual
    };
};

/** A total function unless the NEUTRON_COEFFICIENTS keep one pair only. */
const blockSummary = (dataset) => (dataset
    ? Object.fromEntries(['block', 'file', 'fitType', 'realSpaceFit', 'convolve', 'convolveQ', 'rCutoff', 'qMin', 'qMax']
        .map((key) => [key, dataset[key] ?? null]))
    : null);


/** One fitted dataset of a run: how its model is made, its transform and its r-resolution (assess_dataset). */
export const assessDataset = ({
    name, kind, number, settings, x: xIn, data: dataIn, fq = null, inputColumns = null, dataTransform = null,
    override = null, processing = null
}) => {
    const dataset = nthBlock(settings, BLOCKS[kind] || '', number);
    const real = ['xpdf', 'npdf', 'xray-real'].includes(kind);
    const payload = {
        status: 'unavailable', reason: null, file: name, kind, number,
        block: blockSummary(dataset), fit: fitSpaceRisk(real ? kind : 'reciprocal', dataset, dataTransform),
        transform: null, resolution: null, processing
    };
    if (kind === 'xpdf' && dataset && !dataset.realSpaceFit) {
        payload.fit = { ...payload.fit, fitSpace: 'reciprocal', risk: 'none', because: 'fitted-in-q' };
    } else if (kind === 'xpdf' && !dataset) {
        // No .dat to say whether the transform is fitted: either way no ripple is fitted.
        payload.fit = { ...payload.fit, fitSpace: 'unknown', because: 'same-transform-no-dat' };
    }
    const x = Float64Array.from(xIn);
    const data = Float64Array.from(dataIn);
    if (!real) {
        let dq;
        try {
            dq = uniformStep(x);
        } catch (error) {
            return { ...payload, status: 'failed', reason: `the Q grid: ${error.message}` };
        }
        const weights = transformWeights(x, dq, 'rectangle', 'none');
        payload.transform = {
            qMin: x[0], qMax: x[x.length - 1], dQ: dq, nQ: x.length, rule: 'rectangle', window: 'none', source: 'q-file', windowKnown: false
        };
        payload.resolution = resolution(x, weights);
        payload.status = 'assessed';
        return payload;
    }
    let rr = x;
    let yy = data;
    let qMin;
    let qMax;
    let dq;
    let rule;
    let window;
    let source;
    let windowKnown;
    let transformExtra = {};
    if (kind === 'xpdf') {
        if (fq === null) return { ...payload, reason: 'no-fq-file' };
        const [q, , fExp] = fq.map((column) => Float64Array.from(column));
        try {
            dq = uniformStep(q);
        } catch (error) {
            return { ...payload, status: 'failed', reason: `the Q grid: ${error.message}` };
        }
        const reproduction = reproduceTransform(q, fExp, x, data);
        if (!reproduction.pass) {
            return {
                ...payload, status: 'failed', reason: 'transform-not-reproduced',
                transform: { reproduction: { ...reproduction, tolerance: RMC_FT_TOLERANCE } }
            };
        }
        qMin = q[0];
        qMax = q[q.length - 1];
        rule = 'rectangle';
        window = 'none';
        source = 'run-fq';
        windowKnown = true;
        transformExtra = { reproduction: { residual: reproduction.residual, offset: reproduction.offset, tolerance: RMC_FT_TOLERANCE } };
    } else {
        const given = override || {};
        if (given.qMax !== null && given.qMax !== undefined) {
            source = 'argument';
            qMax = Number(given.qMax);
            qMin = given.qMin ?? null;
            window = given.window || 'none';
            windowKnown = given.window !== null && given.window !== undefined;
            dq = DEFAULT_DQ;
        } else if (dataTransform !== null) {
            ({ source, qMax, qMin, window } = dataTransform);
            windowKnown = true;
            dq = dataTransform.dQ || DEFAULT_DQ;
        } else if (kind === 'npdf' && dataset && dataset.convolveQ) {
            source = 'convolve';
            qMax = dataset.convolveQ;
            qMin = null;
            window = 'none';
            windowKnown = false;
            dq = DEFAULT_DQ;
        } else if (kind === 'xray-real' && dataset && dataset.qMax) {
            source = 'xray-block';
            qMax = dataset.qMax;
            qMin = dataset.qMin;
            window = 'none';
            windowKnown = false;
            dq = dataset.qStep || DEFAULT_DQ;
        } else {
            return { ...payload, reason: 'qmax-unknown' };
        }
        qMin = qMin === null || qMin === undefined || qMin < dq ? dq : qMin;
        rule = source === 'xray-block' ? 'rectangle' : 'trapezoid';
        if (inputColumns !== null) {
            const mapped = mapInput(x, data, inputColumns);
            if (mapped !== null) {
                rr = mapped.r;
                yy = mapped.y;
                transformExtra = { input: { scale: mapped.scale, offset: mapped.offset, residual: mapped.residual } };
            }
        }
    }
    if (!uniform(pick(rr, Array.from(rr, (value, i) => Number.isFinite(value) && Number.isFinite(yy[i]))))) {
        return { ...payload, status: 'failed', reason: 'r-grid-not-uniform' };
    }
    const n = pyRound((qMax - qMin) / dq) + 1;
    const q = qGrid(qMin, dq, n);
    const weights = transformWeights(q, dq, rule, window);
    return {
        ...payload,
        status: 'assessed',
        transform: {
            qMin: q[0], qMax: q[q.length - 1], dQ: dq, nQ: q.length, rule, window, source, windowKnown, ...transformExtra
        },
        resolution: resolution(q, weights)
    };
};
