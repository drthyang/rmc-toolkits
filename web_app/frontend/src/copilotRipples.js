// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The AI Copilot's check_ripples (effect 'compute', llm copilot/effects.js):
// Fourier termination ripples in the run's G(r), so a small peak, a shoulder
// or a split is not read as structure when it is a ripple — RMC is a
// minimiser and fits whatever the data carry, ripples included. The engine
// (rmc_toolkits/ripples.py, its port src/ripples.js through ripplesRun.js, or
// Flask's /api/ripples) decides; this module renders: per fitted dataset how
// RMCProfile makes the model it compares (the risk that the box fits ripples),
// the transform and the r-resolution, the strong peaks' predicted ripples and
// the low-r level; the caller's features read as ripple-consistent, larger
// or smaller than the ripple there, below the resolution, across datasets of
// different Q_max; a plain verdict and the next step. The run's ripple source
// is the adapter's (copilotRun.js LoadedRun `ripples`): (paths, args,
// { signal }) → the engine's result.
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
const roundR = (value) => (Number.isFinite(value) ? Math.round(value * 1000) / 1000 : null);

// How RMCProfile makes the model, in words (fit_space_risk's `model`).
const MODEL_WORDS = {
    'same-transform': 'the same transform of the model F(Q) as the data\'s: truncated alike',
    histogram: 'the configuration\'s histogram, not truncated',
    sinc: 'the histogram convolved with the sinc of the CONVOLVE Q range',
    'xray-transform': 'the histogram transformed over the block\'s QMIN-QMAX and back'
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

// A feature's reading in one dataset (classify_feature's `class`).
export const FEATURE_SENTENCES = {
    'ripple-consistent': 'the size of the ripple predicted there: a termination ripple, not structure',
    'larger-than-ripple': 'larger than the ripple predicted there: more than termination (a weaker peak, the strong peak\'s own shape, or noise)',
    'smaller-than-ripple': 'smaller than the ripple predicted there: the ripple is not there in full (the data\'s transform differs from the one assumed, or something of the other sign fills it)',
    'no-ripple-here': 'no ripple of its sign is predicted there: not termination (a weaker peak, the strong peak\'s own shape, or noise)',
    'within-noise': 'within twice the data\'s noise level (from the low-r region): not resolved from noise',
    'resolved-peak': 'a peak in its own right (the larger peaks\' ripples do not explain it)'
};

// One feature across datasets of different Q_max (across_datasets' `verdict`).
export const ACROSS_SENTENCES = {
    'moves-with-qmax': 'it moves with Q_max: a ripple',
    'not-discriminating': 'the Q_max differ too little to tell a ripple from a feature here',
    'in-every-dataset': 'it is in every dataset and a ripple in none: not termination',
    'ripple-in-every-dataset': 'it is ripple-consistent in every dataset',
    mixed: 'the datasets disagree'
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

const LOW_R_WORDS = {
    measured: 'measured',
    enforced: 'set to the baseline (stog / Auto StoG): the low-r level cannot be read',
    filtered: 'cleaner than the predicted ripples: Fourier-filtered, it does not measure them',
    'too-short': 'too short below the first distance to read',
    unavailable: 'not read'
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

const STRONG_PEAKS_SHOWN = 4;

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
    if (dataset.peaks?.length) {
        row.strong_peaks = [...dataset.peaks]
            .sort((a, b) => (b.rippleScale - a.rippleScale) || (a.r - b.r))
            .slice(0, STRONG_PEAKS_SHOWN)
            .sort((a, b) => a.r - b.r)
            .map((peak) => {
                // The first ripple on each side.
                const left = peak.ripples.filter((ripple) => ripple.r < peak.r).sort((a, b) => b.r - a.r)[0];
                const right = peak.ripples.filter((ripple) => ripple.r > peak.r).sort((a, b) => a.r - b.r)[0];
                const nearest = [left, right].filter(Boolean);
                return {
                    r: roundR(peak.r),
                    sigma_A: roundSig(peak.sigma, 2),
                    first_ripples: nearest.map((ripple) => [roundR(ripple.r), roundSig(ripple.relative, 2)])
                };
            });
    }
    const low = dataset.lowR;
    if (low && low.state !== 'unavailable') {
        row.low_r = { state: low.state };
        if (low.range) row.low_r.range_A = [roundR(low.range[0]), roundR(low.range[1])];
        if (Number.isFinite(low.rms)) {
            const height = Math.max(0, ...(dataset.peaks || []).map((peak) => Math.abs(peak.height)));
            const relative = (value) => (height > 0 && Number.isFinite(value) ? roundSig(value / height, 2) : null);
            row.low_r.rms = relative(low.rms);
            row.low_r.predicted_rms = relative(low.predictedRms);
            if (Number.isFinite(low.modelRms)) row.low_r.model_rms = relative(low.modelRms);
            if (Number.isFinite(low.dataMinusModelRms)) row.low_r.data_minus_model_rms = relative(low.dataMinusModelRms);
            row.low_r.relative_to = 'the largest strong peak';
        }
        if (Number.isFinite(low.enforcedBelow)) row.low_r.enforced_below_A = roundR(low.enforcedBelow);
    }
    return row;
};

const readingRow = (reading, key) => {
    const row = { dataset: key, class: reading.class };
    if (reading.class === 'resolved-peak') {
        row.peak_A = roundR(reading.r);
        return row;
    }
    row.at_A = roundR(reading.r);
    row.size = roundSig(reading.size, 2);
    if (reading.ripple) {
        row.ripple_A = roundR(reading.ripple.r);
        row.ripple_size = roundSig(reading.ripple.amplitude, 2);
    }
    if (reading.fromPeak) row.from_peak_A = roundR(reading.fromPeak.distance);
    row.below_resolution = Boolean(reading.belowResolution);
    return row;
};

const fmt = (value, digits = 3) => String(roundSig(value, digits));

/**
 * The data, checks and verdict of one check_ripples call from the engine's
 * result. `keys` maps a file name to its dataset key; `only` keeps one key.
 */
export const summarizeRipples = (result, { keys, only = null, pair = null } = {}) => {
    const keyOf = (file) => keys.get(file) || file;
    const all = (result.datasets || []).filter((dataset) => keys.has(dataset.file));
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
        const low = dataset.lowR;
        if (low && low.state === 'measured' && Number.isFinite(low.rms) && low.predictedRms > 0) {
            const ratio = low.rms / low.predictedRms;
            // The model, where it is transformed as the data are, has only its own
            // ripples there: what the prediction should be.
            const model = Number.isFinite(low.modelRms)
                ? ` (the model, transformed alike, shows ${fmt(low.modelRms / low.predictedRms, 2)}× the prediction there)` : '';
            checks.push({
                level: 'note',
                text: `${key}: below the first distance (${fmt(low.range[0])}–${fmt(low.range[1])} Å) the data's rms about the baseline `
                    + `is ${fmt(ratio, 2)}× the ripples the strong peaks predict there${model}`
                    + (ratio > 2 ? ': noise or processing errors dominate the low r.' : ': ripple-dominated.')
            });
        } else if (low && ['enforced', 'filtered'].includes(low.state)) {
            checks.push({ level: 'note', text: `${key}: the low-r region is ${LOW_R_WORDS[low.state]}.` });
        }
    }

    // Features (in the datasets, or in the box's partial g(r)).
    if (result.features) {
        data.features = result.features.map((feature) => {
            const readings = feature.datasets.filter((reading) => keys.has(reading.file)
                && (!only || keyOf(reading.file) === only));
            const row = { r: roundR(feature.asked), in: readings.map((reading) => readingRow(reading, keyOf(reading.file))) };
            if (feature.across && !only) {
                row.across = { verdict: feature.across.verdict, discriminating: feature.across.discriminating };
            }
            return row;
        });
        data.found = Boolean(result.found);
        for (const feature of data.features) {
            const parts = feature.in.map((reading) => `in ${reading.dataset} ${FEATURE_SENTENCES[reading.class] || reading.class}`
                + (reading.below_resolution ? ` (and closer to its peak than the resolution)` : ''));
            const across = feature.across ? `; across the datasets ${ACROSS_SENTENCES[feature.across.verdict] || feature.across.verdict}` : '';
            const ripple = feature.in.some((reading) => reading.class === 'ripple-consistent')
                || feature.across?.verdict === 'moves-with-qmax';
            checks.push({ level: ripple ? 'warn' : 'note', text: `r = ${feature.r} Å: ${parts.join('; ')}${across}.` });
        }
    }
    if (result.pair) {
        data.pair = {
            pair: result.pair.label,
            features: result.pair.features.map((feature) => {
                const row = { r: roundR(feature.r), in_partial: feature.inPartial };
                if (Number.isFinite(feature.spacing)) row.spacing_A = roundR(feature.spacing);
                row.below_resolution = Boolean(feature.belowResolution);
                row.at_ripple = feature.atRipple.filter((ripple) => keys.has(ripple.file)).map((ripple) => ({
                    dataset: keyOf(ripple.file), r: roundR(ripple.r), relative: roundSig(ripple.relative, 2)
                }));
                return row;
            })
        };
        for (const feature of data.pair.features) {
            const at = feature.at_ripple.length
                ? `at a predicted ripple of ${feature.at_ripple.map((ripple) => ripple.dataset).join(', ')}, fitted in r against a model not truncated like the data: the box may be making it`
                : 'at no predicted ripple of a dataset whose box may fit ripples';
            checks.push({
                level: feature.at_ripple.length ? 'warn' : 'note',
                text: `${data.pair.pair} g(r) r = ${feature.r} Å: ${at}`
                    + (feature.below_resolution ? `; ${feature.spacing_A} Å from the next peak, below the data's resolution.` : '.')
            });
        }
    } else if (pair) {
        checks.push({ level: 'note', text: `No partial g(r) ${pair} in the run's partials file: the features were read in the datasets.` });
    }
    if (result.split) {
        const below = result.split.datasets.filter((entry) => keys.has(entry.file) && entry.belowResolution);
        data.split = {
            spacing_A: roundR(result.split.spacing),
            below_resolution_in: below.map((entry) => keyOf(entry.file)),
            fwhm_A: Object.fromEntries(result.split.datasets.filter((entry) => keys.has(entry.file))
                .map((entry) => [keyOf(entry.file), roundSig(entry.fwhm)]))
        };
        checks.push({
            level: below.length ? 'warn' : 'note',
            text: below.length
                ? `A split of ${data.split.spacing_A} Å is below the r-resolution of ${data.split.below_resolution_in.join(', ')}: those data cannot resolve it.`
                : `A split of ${data.split.spacing_A} Å is above every dataset's r-resolution.`
        });
    }

    // The verdict, result first.
    const featureRows = data.features || [];
    const ripples = featureRows.filter((row) => row.in.some((reading) => reading.class === 'ripple-consistent'));
    const notRipples = featureRows.filter((row) => row.in.length && row.in.every((reading) => (
        ['larger-than-ripple', 'no-ripple-here', 'resolved-peak'].includes(reading.class))));
    const parts = [];
    if (high.length) parts.push(`the box is asked to fit termination ripples in ${high.map((d) => keyOf(d.file)).join(', ')} (fitted in r against an untruncated or differently truncated model)`);
    else if (unknown.length) parts.push(`whether the box fits ripples in ${unknown.map((d) => keyOf(d.file)).join(', ')} cannot be established from the run's files`);
    else if (real.length) parts.push('no real-space fit asks the box to make ripples (each model is truncated as its data are)');
    else if (datasets.length) parts.push('every dataset is fitted in Q: the box is not asked to make ripples');
    if (ripples.length) parts.push(`ripple-consistent: r = ${ripples.map((row) => row.r).join(', ')} Å`);
    if (notRipples.length) parts.push(`not ripples: r = ${notRipples.map((row) => row.r).join(', ')} Å`);
    const sentence = parts.join('; ');
    data.verdict = parts.length ? `${sentence.charAt(0).toUpperCase()}${sentence.slice(1)}.` : 'No fitted dataset to assess.';

    // The next step (the user runs RMC again only for a specific suspect finding).
    if (high.length) {
        next.push('Fit the data in Q (F(Q)/S(Q)) instead of r, or give RMCProfile the data\'s Q range with CONVOLVE :: Q_max on an unwindowed G(r), so the model is truncated as the data are.');
    }
    if (datasets.some((dataset) => dataset.reason === 'qmax-unknown')) {
        next.push('Give qmax (and window, if the G(r) was transformed with Lorch\'s) for the G(r) whose Q_max the files do not give.');
    }
    if (ripples.length) {
        next.push(`Do not read r = ${ripples.map((row) => row.r).join(', ')} Å as structure; to test it, use data with a higher Q_max, a Lorch-windowed transform (which suppresses the ripples) or the other radiation's data (another Q_max).`);
        if (high.length) {
            next.push('If it matters: one run with the model truncated as the data are (or fitted in Q), from a new start, for that feature only — if it does not return, it was the ripple.');
        }
    }
    if (data.split?.below_resolution_in?.length || featureRows.some((row) => row.in.some((reading) => reading.below_resolution))) {
        next.push('A split or shoulder below the r-resolution is not constrained by the data: only a higher Q_max resolves it.');
    }
    if (featureRows.length && assessed.filter((dataset) => ['real', 'unknown'].includes(dataset.fit?.fitSpace) && dataset.total !== false).length < 2) {
        next.push('One G(r) only: add a dataset with a different Q_max (X-ray beside neutron) to tell a ripple, which moves with Q_max, from a feature, which stays.');
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
        'resolution_A?': { fwhm: 'number|null', nyquist: 'number|null', ripple_period: 'number|null' },
        'strong_peaks?': [{ r: 'number|null', sigma_A: 'number|null', first_ripples: [['number|null', 'number|null']] }],
        'low_r?': {
            state: 'string',
            'range_A?': ['number|null', 'number|null'],
            'rms?': 'number|null',
            'predicted_rms?': 'number|null',
            'model_rms?': 'number|null',
            'data_minus_model_rms?': 'number|null',
            'relative_to?': 'string',
            'enforced_below_A?': 'number|null'
        }
    }],
    'features?': [{
        r: 'number|null',
        in: [{
            dataset: 'string',
            class: 'string',
            'peak_A?': 'number|null',
            'at_A?': 'number|null',
            'size?': 'number|null',
            'ripple_A?': 'number|null',
            'ripple_size?': 'number|null',
            'from_peak_A?': 'number|null',
            'below_resolution?': 'boolean'
        }],
        'across?': { verdict: 'string', discriminating: 'boolean' }
    }],
    'found?': 'boolean',
    'pair?': {
        pair: 'string',
        features: [{
            r: 'number|null', in_partial: 'boolean', 'spacing_A?': 'number|null', below_resolution: 'boolean',
            at_ripple: [{ dataset: 'string', r: 'number|null', relative: 'number|null' }]
        }]
    },
    'split?': { spacing_A: 'number|null', below_resolution_in: ['string'], fwhm_A: { '*': 'number|null' } },
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
    const partials = (plotFiles || []).find((file) => file.plotKind === 'pdf_partials' && file.plotData?.series?.length);
    const pairs = partials
        ? partials.plotData.series.map((entry) => (entry.label || '').trim()).filter((label) => /^[A-Z][a-z]?\s*-\s*[A-Z][a-z]?$/.test(label))
        : [];
    const properties = {
        dataset: {
            type: 'string',
            enum: [...keyed.keys()],
            description: `${[...keyed.entries()].map(([key, file]) => `${key}: ${RIPPLE_FIT_KINDS[file.plotKind]}`).join('; ')}. Default: all.`
        },
        features: {
            type: 'array', items: { type: 'number' }, minItems: 1, maxItems: 8, description: 'Small peaks or shoulders to read (Å).'
        },
        r_min: { type: 'number', minimum: 0, maximum: 100, description: 'Or find them from r_min to r_max (Å).' },
        r_max: { type: 'number', minimum: 0, maximum: 100 },
        split: {
            type: 'array', items: { type: 'number' }, minItems: 2, maxItems: 2, description: 'A split peak\'s two positions (Å).'
        },
        qmax: { type: 'number', minimum: 1, maximum: 200, description: 'Q_max (Å⁻¹) of a G(r) whose files lack it.' },
        window: { type: 'string', enum: ['none', 'lorch'], description: 'That G(r)\'s window.' }
    };
    if (pairs.length) {
        properties.pair = { type: 'string', enum: pairs, description: 'Read the features in this partial g(r) of the box.' };
    }
    return {
        name: 'check_ripples',
        description: 'Fourier termination ripples in the run\'s G(r): per fitted dataset, whether RMCProfile truncates the '
            + 'model as the data (the risk that the box fits ripples), the Q range, window and r-resolution, the strong '
            + 'peaks\' predicted ripples, the low-r level; for features, whether each is a ripple, larger than one, below '
            + 'the resolution, or moves with Q_max across datasets; a verdict and the next step.',
        // A short pointer to when to use it (method.js: sent with method 'full').
        method: 'Use it before reading a small peak, shoulder or split in a G(r) or partial g(r) as structure.',
        parameters: { type: 'object', properties },
        dataShape: RIPPLES_SHAPE,
        label: ({ dataset, pair }) => `Termination ripples${dataset ? ` ${dataset}` : ''}${pair ? ` · ${pair}` : ''}`,
        page: 'dashboard',
        cards: ({ dataset } = {}) => [...keyed.entries()]
            .filter(([key]) => !dataset || key === dataset)
            .map(([, file]) => `plot:${file.name}`),
        run: async (args = {}, { signal } = {}) => {
            const { dataset = null, features = null, r_min: rMin = null, r_max: rMax = null, split = null, qmax = null, window = null } = args;
            const pair = args.pair ?? null;
            const checksOnly = [];
            if ((rMin === null) !== (rMax === null)) checksOnly.push('give r_min and r_max together');
            if (rMin !== null && rMax !== null && !(rMax > rMin)) checksOnly.push('r_max must exceed r_min');
            if (checksOnly.length) return { data: { datasets: [], verdict: 'Not run.' }, checks: [{ level: 'fail', text: `${checksOnly.join('; ')}.` }] };
            const paths = [...keyed.values()].map((file) => file.path || file.name);
            const result = await ripples(paths, {
                features: features && features.length ? features.map(Number) : null,
                rRange: rMin !== null ? [rMin, rMax] : null,
                split: split && split.length === 2 ? split.map(Number) : null,
                pair,
                qmax,
                qmin: null,
                window
            }, { signal });
            return summarizeRipples(result, { keys, only: dataset, pair });
        }
    };
};
