// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The AI Copilot's check_ripples (effect 'compute', llm copilot/effects.js):
// Fourier termination ripples in the run's G(r), from the run's files and
// first principles — RMC is a minimiser and fits whatever the data carry,
// ripples included. The engine (rmc_toolkits/ripples.py, its port
// src/ripples.js through ripplesRun.js, or Flask's /api/ripples) decides; this
// module renders: per fitted dataset how RMCProfile makes the model it
// compares (whether the box is asked to fit ripples), the Q range and window,
// and the r-resolution (the ripple period 2π/Q_max, where a peak's ripples
// fall, the kernel's FWHM); a plain verdict and the next step, which says how
// to test a suspect feature. It does not judge whether a given peak or
// shoulder is a ripple: that reading was not reliable enough to ship (it is
// kept on the branch wip/joint-fit). The run's ripple source is the adapter's
// (copilotRun.js LoadedRun `ripples`): (paths, args, { signal }) → the
// engine's result.
//
// Offered in the dock only for now: the validation's live arms leave it out
// by name (validation/harness.js LIVE_ARM_EXCLUDED) so protocol 0.8.0 and its
// recorded results stay valid; the arms gain it with a protocol bump.

// The run's fit files by plot kind, as compare_fit keys them.
export const RIPPLE_FIT_KINDS = Object.freeze({
    xpdf: 'X-ray G(r), RMCProfile\'s transform of its F(Q)',
    npdf: 'neutron G(r)',
    xray_sq: 'X-ray F(Q)',
    neutron_sq: 'neutron S(Q)'
});

const FIT_FILE_RE = /_(FT_XFQ|Xray_PDF|PDF|FQ|SQ)\d+\.csv$/;

const roundSig = (value, digits = 3) => (
    Number.isFinite(value) && value !== 0 ? Number(value.toPrecision(digits)) : (Number.isFinite(value) ? value : null)
);

// How RMCProfile makes the model, in a few words (fit_space_risk's `model`;
// the risk's sentence in the checks says what it means).
const MODEL_WORDS = {
    'same-transform': 'same transform as the data',
    histogram: 'untruncated histogram',
    sinc: 'histogram ⊗ sinc(CONVOLVE Q)',
    'xray-transform': 'histogram over QMIN-QMAX'
};

// Why a dataset carries the risk it does (fit_space_risk's `because`).
export const RISK_SENTENCES = {
    'fitted-in-q': 'fitted in Q: the box is not asked to make ripples (a G(r) transformed from these data still shows them)',
    'same-transform': 'fitted in r, but RMCProfile transforms the model F(Q) exactly as the data: both carry the same ripples, so the box is not asked to make them',
    'same-transform-no-dat': 'no .dat says whether this transform is fitted; if it is, RMCProfile transforms the model F(Q) exactly as the data, so either way the box is not asked to make its ripples',
    'untruncated-model': 'fitted in r against the configuration\'s untruncated histogram (no CONVOLVE): the box is asked to make the data\'s ripples',
    'sinc-no-q': 'CONVOLVE is set without its Q range: what the model is convolved with is unknown',
    'sinc-data-unverified': 'the model is sinc-convolved at the CONVOLVE Q, but the run\'s files do not say how the data were transformed (their Q_max, any window): whether the two are truncated alike is unknown',
    'sinc-vs-lorch': 'the model is sinc-convolved but the data were transformed with a Lorch window: truncated differently, so the box is asked to make the difference',
    'sinc-q-differs': 'the CONVOLVE Q differs from the data\'s Q_max: truncated differently, so the box is asked to make the difference',
    'sinc-matches': 'the model is sinc-convolved at the data\'s Q_max with no window on either: truncated alike',
    'xray-real-unverified': 'the model is transformed over the block\'s QMIN-QMAX, but the run\'s files do not say how the data were: whether the two are truncated alike is unknown',
    'transform-vs-lorch': 'the model is transformed with no window but the data with a Lorch window: truncated differently',
    'transform-q-differs': 'the block\'s QMAX differs from the data\'s Q_max: truncated differently',
    'transform-matches': 'the model is transformed over the data\'s Q range with no window on either: truncated alike',
    'unknown-kind': 'how RMCProfile makes this model is not established'
};

// Why a dataset could not be assessed (the engine's `reason` codes).
const REASON_SENTENCES = {
    'qmax-unknown': 'the run\'s files do not give the data\'s Q_max (no stog input or Auto StoG provenance names the file, no CONVOLVE Q): give qmax',
    'no-fq-file': 'no F(Q) file beside the transform to rebuild it from',
    'transform-not-reproduced': 'RMCProfile\'s transform of the run\'s F(Q) does not reproduce this G(r)',
    'r-grid-not-uniform': 'its r grid is not uniform'
};

// Where a dataset's transform comes from (the engine's transform `source`).
const SOURCE_WORDS = {
    'run-fq': 'the run\'s F(Q) file',
    'q-file': 'this file\'s Q grid',
    stog: 'the stog input',
    provenance: 'Auto StoG\'s provenance',
    convolve: 'the .dat\'s CONVOLVE Q',
    argument: 'given (qmax)',
    'xray-block': 'the .dat\'s QMIN-QMAX'
};

const keyFiles = (plotFiles) => {
    const keyed = new Map();
    (plotFiles || [])
        .filter((file) => RIPPLE_FIT_KINDS[file.plotKind] && FIT_FILE_RE.test(file.name || ''))
        .forEach((file) => {
            let key = file.plotKind;
            for (let n = 2; keyed.has(key); n += 1) key = `${file.plotKind}_${n}`;
            keyed.set(key, file);
        });
    return keyed;
};

/** check_ripples' result budget (characters): its datasets whole. */
export const RIPPLES_RESULT_BUDGET = 4500;

const datasetRow = (dataset, key) => {
    const row = { dataset: key, file: dataset.file, status: dataset.status };
    if (dataset.status !== 'assessed') {
        row.why = REASON_SENTENCES[dataset.reason] || dataset.reason || 'not assessed';
    }
    if (dataset.fit) {
        row.fitted_in = { reciprocal: 'Q', real: 'r' }[dataset.fit.fitSpace] || 'unknown';
        if (dataset.fit.model) row.model = MODEL_WORDS[dataset.fit.model] || dataset.fit.model;
        row.risk = dataset.fit.risk;
    }
    const transform = dataset.transform;
    if (transform && Number.isFinite(transform.qMax)) {
        row.Q_A = [roundSig(transform.qMin, 4), roundSig(transform.qMax, 4)];
        // A fit in Q has no transform, so no window; a G(r)'s is the files' or assumed.
        if (dataset.fit && dataset.fit.fitSpace !== 'reciprocal') {
            row.window = transform.windowKnown ? transform.window : `${transform.window} (assumed)`;
        }
        row.Q_from = SOURCE_WORDS[transform.source] || transform.source;
    }
    if (dataset.resolution) {
        row.resolution_A = {
            fwhm: roundSig(dataset.resolution.fwhm), nyquist: roundSig(dataset.resolution.nyquistStep),
            ripple_period: roundSig(dataset.resolution.ripplePeriod)
        };
    }
    return row;
};

const fmt = (value, digits = 3) => String(roundSig(value, digits));

/**
 * The data, checks and verdict of one check_ripples call from the engine's
 * result. `keys` maps a file name to its dataset key; `only` keeps one key.
 */
export const summarizeRipples = (result, { keys, only = null } = {}) => {
    const keyOf = (file) => keys.get(file) || file;
    // Real-space datasets first, the riskiest first (the result budget halves
    // lists from the end), then the fits in Q.
    const rank = (dataset) => {
        if (!dataset.fit || dataset.fit.fitSpace === 'reciprocal') return 4;
        return { high: 0, unknown: 1, matched: 2 }[dataset.fit.risk] ?? 3;
    };
    const all = (result.datasets || []).filter((dataset) => keys.has(dataset.file))
        .map((dataset, index) => ({ dataset, index }))
        .sort((a, b) => (rank(a.dataset) - rank(b.dataset)) || (a.index - b.index))
        .map(({ dataset }) => dataset);
    const datasets = only ? all.filter((dataset) => keyOf(dataset.file) === only) : all;
    const data = { datasets: datasets.map((dataset) => datasetRow(dataset, keyOf(dataset.file))) };
    const checks = [];
    const next = [];
    const real = datasets.filter((dataset) => dataset.fit && ['real', 'unknown'].includes(dataset.fit.fitSpace));
    const high = real.filter((dataset) => dataset.fit.risk === 'high');
    const unknown = real.filter((dataset) => dataset.fit.risk === 'unknown');
    const assessed = datasets.filter((dataset) => dataset.status === 'assessed');

    for (const dataset of datasets) {
        const key = keyOf(dataset.file);
        if (dataset.fit) {
            const level = dataset.fit.risk === 'high' ? 'warn' : 'note';
            checks.push({ level, text: `${key}: ${RISK_SENTENCES[dataset.fit.because] || dataset.fit.because}.` });
        }
        if (dataset.status !== 'assessed') {
            const reciprocal = all.find((other) => other.kind === 'neutron-reciprocal' && other.transform);
            const hint = dataset.reason === 'qmax-unknown' && reciprocal
                ? ` (the run's ${keyOf(reciprocal.file)} reaches Q = ${fmt(reciprocal.transform.qMax, 4)} Å⁻¹, if the G(r) came from the same S(Q))`
                : '';
            checks.push({ level: 'warn', text: `${key}: not assessed: ${REASON_SENTENCES[dataset.reason] || dataset.reason}${hint}.` });
        }
        if (dataset.status === 'assessed' && dataset.fit && ['real', 'unknown'].includes(dataset.fit.fitSpace)
            && dataset.resolution && dataset.transform) {
            // Where a peak's ripples fall: the truncation kernel sin(Q_max·Δr)/Δr
            // has its extrema where tan x = x (x = Q_max·Δr), x ≈ (n + ½)π for
            // n ≥ 1, so at Δr ≈ (n + ½)·π/Q_max. Their sign depends on the
            // peak's width (σ·Q_max), so no sign is stated.
            const { fwhm, ripplePeriod } = dataset.resolution;
            const lorch = dataset.transform.window === 'lorch' ? ' (Lorch window: ripples suppressed, resolution coarser)' : '';
            checks.push({
                level: 'note',
                text: `${key}: r-resolution: ripples repeat every ${fmt(ripplePeriod)} Å (2π/Q_max, Q_max = `
                    + `${fmt(dataset.transform.qMax, 4)} Å⁻¹)${lorch}; a peak at r₀ has its ripples near r₀ ± (n + ½)·`
                    + `${fmt(ripplePeriod / 2)} Å (n = 1, 2, …; π/Q_max), their sign set by the peak's width; features `
                    + `closer than the kernel's FWHM, ${fmt(fwhm)} Å, are not resolved.`
            });
        }
    }

    // The verdict, result first.
    let verdict = null;
    if (high.length) verdict = `the box is asked to fit termination ripples in ${high.map((d) => keyOf(d.file)).join(', ')} (fitted in r against an untruncated or differently truncated model)`;
    else if (unknown.length) verdict = `whether the box fits ripples in ${unknown.map((d) => keyOf(d.file)).join(', ')} cannot be established from the run's files`;
    else if (real.length) verdict = 'no real-space fit asks the box to make ripples (each model is truncated as its data are)';
    else if (datasets.length) verdict = 'every dataset is fitted in Q: the box is not asked to make ripples';
    data.verdict = verdict ? `${verdict.charAt(0).toUpperCase()}${verdict.slice(1)}.` : 'No fitted dataset to assess.';

    // The next step (the user runs RMC again only for a specific suspect finding).
    if (high.length) {
        next.push('Fit the data in Q (F(Q)/S(Q)) instead of r, or give RMCProfile the data\'s Q range with CONVOLVE :: Q_max on an unwindowed G(r), so the model is truncated as the data are.');
    }
    if (datasets.some((dataset) => dataset.reason === 'qmax-unknown')) {
        next.push('Give qmax (and window, if the G(r) was transformed with Lorch\'s) for the G(r) whose Q_max the files do not give.');
    }
    if (assessed.some((dataset) => ['real', 'unknown'].includes(dataset.fit?.fitSpace))) {
        next.push('This check does not say whether a given peak or shoulder is a ripple. To test one: a ripple moves with Q_max and a real feature stays, so compare a G(r) transformed with a higher Q_max, a Lorch-windowed transform (which suppresses ripples) or the other radiation\'s data (another Q_max).');
    }
    if (high.length) {
        next.push('If a feature matters: one run with the model truncated as the data are (or fitted in Q), from a new start; if the feature does not return, it was a ripple.');
    }
    if (next.length) data.next = next;
    checks.unshift({ level: assessed.length || datasets.some((d) => d.fit) ? 'ok' : 'fail', text: data.verdict });
    return { data, checks };
};

// The shape of check_ripples' data (dataShape.js).
const RIPPLES_SHAPE = {
    datasets: [{
        dataset: 'string',
        file: 'string',
        status: 'string',
        'why?': 'string',
        'fitted_in?': 'string',
        'model?': 'string',
        'risk?': 'string',
        'Q_A?': ['number|null', 'number|null'],
        'window?': 'string',
        'Q_from?': 'string',
        'resolution_A?': { fwhm: 'number|null', nyquist: 'number|null', ripple_period: 'number|null' }
    }],
    verdict: 'string',
    'next?': ['string']
};

/**
 * check_ripples for the run's fit files (`plotFiles`, the props' own) through
 * `ripples` (paths, args, { signal }) → the engine's result; null without a
 * source or a fit file.
 */
export const buildRippleTool = ({ plotFiles = null, ripples = null } = {}) => {
    const keyed = keyFiles(plotFiles);
    if (!ripples || !keyed.size) return null;
    const keys = new Map([...keyed.entries()].map(([key, file]) => [file.name, key]));
    const properties = {
        dataset: {
            type: 'string',
            enum: [...keyed.keys()],
            description: `${[...keyed.entries()].map(([key, file]) => `${key}: ${RIPPLE_FIT_KINDS[file.plotKind]}`).join('; ')}. Default: all.`
        },
        qmax: { type: 'number', minimum: 1, maximum: 200, description: 'Q_max (Å⁻¹) of a G(r) whose files lack it.' },
        window: { type: 'string', enum: ['none', 'lorch'], description: 'That G(r)\'s window.' }
    };
    return {
        name: 'check_ripples',
        description: 'Fourier termination ripples in the run\'s G(r), from the run\'s files and first principles: per '
            + 'fitted dataset, whether RMCProfile truncates the model as the data (whether the box is asked to fit '
            + 'ripples), the Q range and window, and the r-resolution (the ripple period 2π/Q_max and the kernel\'s '
            + 'FWHM). It does not judge whether a given peak or shoulder is a ripple.',
        // A short pointer to when to use it (method.js: sent with method 'full').
        method: 'Use it before reading a small peak, shoulder or split in a G(r) or partial g(r) as structure.',
        parameters: { type: 'object', properties },
        dataShape: RIPPLES_SHAPE,
        // Every dataset whole: within the default budget the rows and their
        // checks would be cut (as analyze_clusters' and analyze_modulations').
        resultBudget: RIPPLES_RESULT_BUDGET,
        label: ({ dataset }) => `Termination ripples${dataset ? ` ${dataset}` : ''}`,
        page: 'dashboard',
        cards: ({ dataset } = {}) => [...keyed.entries()]
            .filter(([key]) => !dataset || key === dataset)
            .map(([, file]) => `plot:${file.name}`),
        run: async (args = {}, { signal } = {}) => {
            const { dataset = null, qmax = null, window = null } = args;
            const paths = [...keyed.values()].map((file) => file.path || file.name);
            const result = await ripples(paths, { qmax, qmin: null, window }, { signal });
            return summarizeRipples(result, { keys, only: dataset });
        }
    };
};
