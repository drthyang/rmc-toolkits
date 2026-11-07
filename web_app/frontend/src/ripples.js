// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Fourier termination ripples: a line-for-line port of rmc_toolkits/ripples.py
// (the source of truth; its docstrings hold the reasoning, this file the
// arithmetic). What the Q range of a run's data puts into its G(r), per
// dataset: the fit space and how RMCProfile makes the model it compares (the
// risk that the box fits ripples), the transform (Q_min, Q_max, ΔQ, window),
// the r-resolution, the strong peaks as Gaussians fitted through the
// transform and their predicted ripples, the low-r level against them, and a
// feature's reading (ripple-consistent, larger or smaller than the ripple
// there, a resolved peak, below the resolution, across datasets of different
// Q_max). Parity golden ripples_fixture.json (tests/generate_ripples_fixture.py
// → src/__tests__/ripplesParity.test.js). Pure: no I/O (ripplesRun.js reads a
// run's files).

import { FT_FILE_RE, RMC_FT_TOLERANCE, reproduceTransform, uniformStep } from './ftUncertainty';

// --- rules (rmc_toolkits/ripples.py) -----------------------------------------------

export const RIPPLE_SIZE_FACTOR = 1.5;
export const RIPPLE_POSITION_FRACTION = 0.25;
export const NEGLIGIBLE_RIPPLE_FRACTION = 0.1;
export const STRONG_PEAK_FRACTION = 0.1;
export const MAX_PEAK_CANDIDATES = 24;
export const MAX_STRONG_PEAKS = 12;
export const STRONG_RIPPLE_FRACTION = 0.02;
export const PEAK_WINDOW_FRACTION = 0.3;
export const PEAK_WINDOW_MIN_POINTS = 2;
export const PEAK_WINDOW_MAX_PERIODS = 1.0;
export const RIPPLE_EXPLAINED_FRACTION = 0.5;
export const BACKFIT_SWEEPS = 2;
export const RIPPLE_REACH_PERIODS = 10;
export const RIPPLES_LISTED = 2;
export const LOW_R_MIN_PERIODS = 1.0;
export const ENFORCED_RTOL = 1e-9;
export const FILTERED_RATIO = 0.3;
export const PRESENT_FRACTION = 0.3;
export const BACKGROUND_PERIODS = 2.0;
export const MAX_FEATURES = 8;
export const FEATURE_MIN_FRACTION = 0.01;
export const SIGMA_MIN_QMAX = 0.1;
export const SIGMA_MAX = 1.0;
export const FIT_MAX_ITERATIONS = 60;
export const FIT_STEP_TOL = 1e-10;
export const FIT_POLISH_STEPS = 4;
export const DEFAULT_DQ = 0.01;
export const INPUT_MAP_TOLERANCE = 1e-5;

export const RISKS = Object.freeze(['high', 'matched', 'unknown', 'none']);
export const RULES = Object.freeze({
    rippleSizeFactor: RIPPLE_SIZE_FACTOR,
    ripplePositionFraction: RIPPLE_POSITION_FRACTION,
    strongPeakFraction: STRONG_PEAK_FRACTION,
    maxStrongPeaks: MAX_STRONG_PEAKS,
    lowRMinPeriods: LOW_R_MIN_PERIODS,
    filteredRatio: FILTERED_RATIO,
    presentFraction: PRESENT_FRACTION
});

const SQRT_2PI = Math.sqrt(2 * Math.PI);
const TWO_OVER_PI = 2 / Math.PI;

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
const mean = (values) => sum(values) / values.length;
const rmsOf = (values) => Math.sqrt(dot(values, values) / values.length);
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

/** np.interp(at, x, y) for increasing x. */
const interp = (at, x, y) => {
    const n = x.length;
    if (at <= x[0]) return y[0];
    if (at >= x[n - 1]) return y[n - 1];
    let lo = 0;
    let hi = n - 1;
    while (hi - lo > 1) {
        const mid = (lo + hi) >> 1;
        if (x[mid] <= at) lo = mid; else hi = mid;
    }
    if (x[hi] === x[lo]) return y[lo];
    return y[lo] + (at - x[lo]) * (y[hi] - y[lo]) / (x[hi] - x[lo]);
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

// --- a Gaussian peak through the transform -----------------------------------------------

/** A[e^{−(r−r₀)²/2σ²} − e^{−(r+r₀)²/2σ²}]. */
export const gaussianD = (r, peak) => Float64Array.from(r, (value) => {
    const s = peak.sigma;
    return peak.amplitude * (Math.exp(-((value - peak.r) ** 2) / (2 * s * s)) - Math.exp(-((value + peak.r) ** 2) / (2 * s * s)));
});

/** A√(2π)σ e^{−Q²σ²/2} sin(Q r₀). */
export const gaussianSine = (q, peak) => Float64Array.from(q, (value) => {
    const s = peak.sigma;
    return peak.amplitude * SQRT_2PI * s * Math.exp(-0.5 * value * value * s * s) * Math.sin(value * peak.r);
});

/** (2/π) Σ_k w_k S_k sin(Q_k r) at every r. */
export const applyTransform = (q, weights, spectrum, r) => {
    const coefficients = Float64Array.from(weights, (w, k) => TWO_OVER_PI * w * spectrum[k]);
    return Float64Array.from(r, (at) => {
        let total = 0;
        for (let k = 0; k < q.length; k += 1) total += Math.sin(at * q[k]) * coefficients[k];
        return total;
    });
};

/** Δ(r) = T[d](r) − T_[Q_min,∞)[d](r) (ripple_pattern). */
export const ripplePattern = (q, weights, dq, rule, peak, r) => {
    const truncated = applyTransform(q, weights, gaussianSine(q, peak), r);
    const belowCount = Math.floor(q[0] / dq + 1e-9);
    const below = [];
    for (let m = 1; m <= belowCount; m += 1) {
        const value = q[0] - dq * m;
        if (value > 0) below.push(value);
    }
    const belowQ = Float64Array.from(below);
    const low = belowQ.length
        ? applyTransform(belowQ, new Float64Array(belowQ.length).fill(dq), gaussianSine(belowQ, peak), r)
        : null;
    const firstWeight = dq * (rule === 'trapezoid' && q.length > 1 ? 0.5 : 1);
    const head = gaussianSine(q.subarray ? q.subarray(0, 1) : [q[0]], peak)[0];
    const d = gaussianD(r, peak);
    return Float64Array.from(r, (at, i) => {
        const edge = TWO_OVER_PI * (dq - firstWeight) * head * Math.sin(q[0] * at);
        return truncated[i] - d[i] + (low ? low[i] : 0) + edge;
    });
};

// --- strong peaks, fitted through the transform ------------------------------------------

/** Interior indices of strict local extrema (plateaus to their first point). */
export const localExtrema = (y) => {
    const out = [];
    for (let i = 1; i < y.length - 1; i += 1) {
        const left = y[i - 1];
        const here = y[i];
        const right = y[i + 1];
        if ((here > left && here >= right) || (here < left && here <= right)) out.push(i);
    }
    return out;
};

const prominences = (y, extrema) => extrema.map((i, position) => {
    const isMax = y[i] > y[i - 1];
    let left = null;
    let right = null;
    for (let j = position - 1; j >= 0; j -= 1) {
        const k = extrema[j];
        if ((y[k] < y[k - 1]) === isMax) { left = k; break; }
    }
    for (let j = position + 1; j < extrema.length; j += 1) {
        const k = extrema[j];
        if ((y[k] < y[k - 1]) === isMax) { right = k; break; }
    }
    const leftValue = left !== null ? y[left] : y[0];
    const rightValue = right !== null ? y[right] : y[y.length - 1];
    const base = 0.5 * (leftValue + rightValue);
    return [i, Math.abs(y[i] - base), base];
});

const copysignOne = (value) => (value < 0 || Object.is(value, -0) ? -1 : 1);

const fitWindow = (y, centre, base) => {
    const height = y[centre] - base;
    const limit = PEAK_WINDOW_FRACTION * Math.abs(height);
    const sign = copysignOne(height);
    let lo = centre;
    while (lo - 1 >= 0 && (y[lo - 1] - base) * sign >= limit) lo -= 1;
    let hi = centre;
    while (hi + 1 < y.length && (y[hi + 1] - base) * sign >= limit) hi += 1;
    lo = Math.max(0, Math.min(lo, centre - PEAK_WINDOW_MIN_POINTS));
    hi = Math.min(y.length - 1, Math.max(hi, centre + PEAK_WINDOW_MIN_POINTS));
    return [lo, hi + 1];
};

const det3 = (a) => (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
    - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
    + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]));

const solve3 = (m, v) => {
    const d = det3(m);
    if (d === 0 || !Number.isFinite(d)) return null;
    const out = [];
    for (let column = 0; column < 3; column += 1) {
        const a = [0, 1, 2].map((row) => [0, 1, 2].map((k) => (k === column ? v[row] : m[row][k])));
        out.push(det3(a) / d);
    }
    return out;
};

const sineRows = (r, q) => Array.from(r, (at) => Float64Array.from(q, (value) => Math.sin(at * value)));

const projected = (sines, q, weights, target, r, r0, sigma) => {
    const spectrum = gaussianSine(q, { amplitude: 1, r: r0, sigma });
    const coefficients = Float64Array.from(weights, (w, k) => TWO_OVER_PI * w * spectrum[k]);
    const phi = Float64Array.from(sines, (row) => dot(row, coefficients));
    const x = Float64Array.from(r, (value) => value - r0);
    const ones = new Float64Array(phi.length).fill(1);
    const columns = [phi, ones, x];
    const m = columns.map((a) => columns.map((b) => dot(a, b)));
    const v = columns.map((a) => dot(a, target));
    const solution = solve3(m, v);
    if (solution === null) return null;
    const [amplitude, b0, b1] = solution;
    const residual = Float64Array.from(target, (value, i) => value - amplitude * phi[i] - b0 - b1 * x[i]);
    return [amplitude, b0, b1, residual];
};

/** A Gaussian through the transform with a line under it, fitted on a window (fit_peak). */
export const fitPeak = (r, target, q, weights, r0, sigma, bounds) => {
    const qmax = q[q.length - 1];
    const sLo = Math.log(SIGMA_MIN_QMAX / qmax);
    const sHi = Math.log(SIGMA_MAX);
    const sines = sineRows(r, q);
    let theta = [r0, Math.min(Math.max(Math.log(sigma), sLo), sHi)];
    const [rLo, rHi] = bounds;
    const clamp = (values) => [Math.min(Math.max(values[0], rLo), rHi), Math.min(Math.max(values[1], sLo), sHi)];
    const evaluate = (values) => projected(sines, q, weights, target, r, values[0], Math.exp(values[1]));
    let current = evaluate(theta);
    if (current === null) return null;
    let cost = dot(current[3], current[3]);
    let lambda = 1e-3;
    const steps = [1e-5, 1e-5];
    let moved = 0;
    let previous = cost;
    for (let iteration = 0; iteration < FIT_MAX_ITERATIONS; iteration += 1) {
        const columns = [];
        for (let j = 0; j < 2; j += 1) {
            const plus = [...theta];
            const minus = [...theta];
            plus[j] += steps[j];
            minus[j] -= steps[j];
            const fp = evaluate(plus);
            const fm = evaluate(minus);
            if (fp === null || fm === null) return null;
            columns.push(Float64Array.from(fp[3], (value, i) => (value - fm[3][i]) / (2 * steps[j])));
        }
        const [j0, j1] = columns;
        const f = current[3];
        const h00 = dot(j0, j0);
        const h01 = dot(j0, j1);
        const h11 = dot(j1, j1);
        const g0 = dot(j0, f);
        const g1 = dot(j1, f);
        let accepted = false;
        while (lambda < 1e12) {
            const a00 = h00 * (1 + lambda);
            const a11 = h11 * (1 + lambda);
            const det = a00 * a11 - h01 * h01;
            if (det <= 0 || !Number.isFinite(det)) { lambda *= 4; continue; }
            const d0 = -(a11 * g0 - h01 * g1) / det;
            const d1 = -(a00 * g1 - h01 * g0) / det;
            const trial = clamp([theta[0] + d0, theta[1] + d1]);
            const result = evaluate(trial);
            if (result !== null) {
                const trialCost = dot(result[3], result[3]);
                if (trialCost < cost) {
                    moved = Math.max(Math.abs(trial[0] - theta[0]), Math.abs(trial[1] - theta[1]));
                    previous = cost;
                    theta = trial;
                    current = result;
                    cost = trialCost;
                    lambda = Math.max(lambda / 3, 1e-12);
                    accepted = true;
                    break;
                }
            }
            lambda *= 4;
        }
        if (!accepted || moved < FIT_STEP_TOL || previous - cost <= 1e-15 * previous) break;
    }
    // Polish: Gauss–Newton steps to the stationary point (fit_peak).
    for (let polish = 0; polish < FIT_POLISH_STEPS; polish += 1) {
        const columns = [];
        for (let j = 0; j < 2; j += 1) {
            const plus = [...theta];
            const minus = [...theta];
            plus[j] += steps[j];
            minus[j] -= steps[j];
            const fp = evaluate(plus);
            const fm = evaluate(minus);
            if (fp === null || fm === null) break;
            columns.push(Float64Array.from(fp[3], (value, i) => (value - fm[3][i]) / (2 * steps[j])));
        }
        if (columns.length < 2) break;
        const [j0, j1] = columns;
        const f = current[3];
        const h00 = dot(j0, j0);
        const h01 = dot(j0, j1);
        const h11 = dot(j1, j1);
        const g0 = dot(j0, f);
        const g1 = dot(j1, f);
        const det = h00 * h11 - h01 * h01;
        if (det <= 0 || !Number.isFinite(det)) break;
        const trial = clamp([theta[0] - (h11 * g0 - h01 * g1) / det, theta[1] - (h00 * g1 - h01 * g0) / det]);
        const result = evaluate(trial);
        if (result === null || dot(result[3], result[3]) > cost * (1 + 1e-9)) break;
        theta = trial;
        current = result;
        cost = dot(result[3], result[3]);
    }
    const [amplitude, background, slope, residual] = current;
    return {
        amplitude, r: theta[0], sigma: Math.exp(theta[1]), background, slope, rms: Math.sqrt(dot(residual, residual) / residual.length)
    };
};

/** |A| σ e^{−Q_max²σ²/2}. */
export const rippleScale = (peak, qmax) => {
    const s = peak.sigma;
    return Math.abs(peak.amplitude) * s * Math.exp(-0.5 * qmax * qmax * s * s);
};

const peakLike = (peak, sign, halfWidth) => (
    peak !== null && peak.amplitude * sign > 0 && peak.sigma < SIGMA_MAX && peak.sigma <= halfWidth
);

const truncatedModels = (r, q, weights, dq, rule, peaks) => {
    const total = new Float64Array(r.length);
    for (const peak of peaks) {
        const d = gaussianD(r, peak);
        const delta = ripplePattern(q, weights, dq, rule, peak, r);
        for (let i = 0; i < r.length; i += 1) total[i] += d[i] + delta[i];
    }
    return total;
};

const byCandidate = (a, b) => (b[1] - a[1]) || (a[0] - b[0]);

/** The strong peaks of d(r) (strong_peaks). */
export const strongPeaks = (r, d, q, weights, dq, rule, rFloor = null, signs = 'positive') => {
    const qmax = q[q.length - 1];
    const period = 2 * Math.PI / qmax;
    const step = r[1] - r[0];
    const scored = prominences(d, localExtrema(d)).filter((entry) => (
        (rFloor === null || r[entry[0]] >= rFloor) && (signs === 'both' || d[entry[0]] > entry[2])
    ));
    if (!scored.length) return [];
    const largest = Math.max(...scored.map((entry) => entry[1]));
    if (!(largest > 0)) return [];
    const candidates = scored.filter((entry) => entry[1] >= STRONG_PEAK_FRACTION * largest)
        .sort(byCandidate).slice(0, MAX_PEAK_CANDIDATES);
    const cap = Math.max(PEAK_WINDOW_MIN_POINTS, pyRound(PEAK_WINDOW_MAX_PERIODS * period / step));
    const accepted = [];
    for (const [index, prominence, base] of candidates) {
        let [lo, hi] = fitWindow(d, index, base);
        lo = Math.max(lo, index - cap);
        hi = Math.min(hi, index + cap + 1);
        const windowR = r.slice(lo, hi);
        const models = truncatedModels(windowR, q, weights, dq, rule, accepted.map((entry) => entry.peak));
        const target = Float64Array.from(windowR, (_, i) => d[lo + i] - models[i]);
        const halfWidth = Math.max(step, 0.5 * (r[hi - 1] - r[lo]));
        const sigma0 = Math.max(halfWidth / 1.1774, 0.5 / qmax);
        const peak = fitPeak(windowR, target, q, weights, r[index], sigma0, [r[lo], r[hi - 1]]);
        const sign = d[index] - base;
        if (!peakLike(peak, sign, halfWidth)) continue;
        const unit = { amplitude: 1, r: peak.r, sigma: peak.sigma };
        const height = Math.abs(peak.amplitude + 0) * Math.abs(applyTransform(q, weights, gaussianSine(q, unit), [peak.r])[0]);
        if (height < RIPPLE_EXPLAINED_FRACTION * prominence) continue;
        accepted.push({ peak, window: [lo, hi], sign, halfWidth });
    }
    accepted.sort((a, b) => a.peak.r - b.peak.r);
    for (let sweep = 0; sweep < BACKFIT_SWEEPS; sweep += 1) {
        accepted.forEach((entry, j) => {
            const [lo, hi] = entry.window;
            const windowR = r.slice(lo, hi);
            const others = accepted.filter((_, i) => i !== j).map((other) => other.peak);
            const models = truncatedModels(windowR, q, weights, dq, rule, others);
            const target = Float64Array.from(windowR, (_, i) => d[lo + i] - models[i]);
            const refit = fitPeak(windowR, target, q, weights, entry.peak.r, entry.peak.sigma, [r[lo], r[hi - 1]]);
            if (peakLike(refit, entry.sign, entry.halfWidth)) entry.peak = refit;
        });
    }
    let kept = accepted.map((entry) => ({ ...entry.peak, rippleScale: rippleScale(entry.peak, qmax) }));
    if (!kept.length) return [];
    const top = Math.max(...kept.map((peak) => peak.rippleScale));
    kept = kept.filter((peak) => peak.rippleScale >= STRONG_RIPPLE_FRACTION * top);
    kept.sort((a, b) => (b.rippleScale - a.rippleScale) || (a.r - b.r));
    kept = kept.slice(0, MAX_STRONG_PEAKS);
    kept.sort((a, b) => a.r - b.r);
    return kept;
};

/** Σ_j Δ_j(r). */
export const patternOn = (r, q, weights, dq, rule, peaks) => {
    const total = new Float64Array(r.length);
    for (const peak of peaks) {
        const delta = ripplePattern(q, weights, dq, rule, peak, r);
        for (let i = 0; i < r.length; i += 1) total[i] += delta[i];
    }
    return total;
};

/** Σ_j d_j(r). */
export const peaksOn = (r, peaks) => {
    const total = new Float64Array(r.length);
    for (const peak of peaks) {
        const d = gaussianD(r, peak);
        for (let i = 0; i < r.length; i += 1) total[i] += d[i];
    }
    return total;
};

/** Where a peak's ripples begin: beyond the kernel's main lobe and 2σ. */
export const rippleCore = (peak, mainLobe) => Math.max(mainLobe, 2 * peak.sigma);

const listedRipples = (r, q, weights, dq, rule, peak, mainLobe, toY) => {
    const reach = RIPPLE_REACH_PERIODS * 2 * Math.PI / q[q.length - 1];
    const rr = Float64Array.from(Array.from(r).filter((value) => Math.abs(value - peak.r) <= reach));
    if (rr.length < 3) return [];
    const delta = toY(rr, ripplePattern(q, weights, dq, rule, peak, rr));
    const height = toY([peak.r], [peak.amplitude])[0];
    const core = rippleCore(peak, mainLobe);
    const found = localExtrema(delta).filter((i) => Math.abs(rr[i] - peak.r) > core);
    const left = found.filter((i) => rr[i] < peak.r).sort((a, b) => rr[b] - rr[a]).slice(0, RIPPLES_LISTED);
    const right = found.filter((i) => rr[i] > peak.r).sort((a, b) => rr[a] - rr[b]).slice(0, RIPPLES_LISTED);
    return [...left, ...right].sort((a, b) => rr[a] - rr[b]).map((i) => ({
        r: rr[i], amplitude: delta[i], relative: height !== 0 ? delta[i] / height : null
    }));
};

// --- the low-r level ----------------------------------------------------------------------

const baselineOf = (r, y, form) => {
    if (form === 'constant') return mean(y);
    const denominator = dot(r, r);
    return denominator > 0 ? dot(r, y) / denominator : 0;
};
const baselineValues = (r, value, form) => Float64Array.from(r, (at) => (form === 'constant' ? value : value * at));

/** The largest r up to which y is exactly its baseline from the first point on, or null. */
export const enforcedBelow = (r, y, form, scale) => {
    if (r.length < 2 || !(scale > 0)) return null;
    let reference;
    if (form === 'constant') reference = y[0];
    else reference = r[0] !== 0 ? y[0] / r[0] : null;
    if (reference === null) return null;
    let last = null;
    for (let i = 0; i < r.length; i += 1) {
        const expected = form === 'constant' ? reference : reference * r[i];
        if (Math.abs(y[i] - expected) <= ENFORCED_RTOL * scale) last = i; else break;
    }
    return last !== null && last >= 1 ? r[last] : null;
};

/** The rms of y about its baseline below the first interatomic distance (low_r_level). */
export const lowRLevel = (r, y, predicted, { form, rFirst, fwhm, period, scale, model = null }) => {
    const out = {
        state: 'unavailable', range: null, baseline: null, rms: null, peak: null, predictedRms: null,
        residualRms: null, modelRms: null, dataMinusModelRms: null, enforcedBelow: null
    };
    const hi = rFirst - fwhm;
    const keep = Array.from(r, (value) => value > 0 && value <= hi);
    if (!keep.some(Boolean)) {
        out.state = 'too-short';
        return out;
    }
    const rKeep = pick(r, keep);
    let rr = rKeep;
    let yy = pick(y, keep);
    let pp = pick(predicted, keep);
    let m = model;
    const enforced = enforcedBelow(rr, yy, form, scale);
    out.enforcedBelow = enforced;
    if (enforced !== null) {
        const rest = Array.from(rr, (value) => value > enforced);
        rr = pick(rr, rest);
        yy = pick(yy, rest);
        pp = pick(pp, rest);
        if (m !== null) m = pick(pick(m, keep), rest);
        if (rr.length < 2 || (rr[rr.length - 1] - rr[0]) < LOW_R_MIN_PERIODS * period) {
            out.state = 'enforced';
            out.range = [rKeep[0], rKeep[rKeep.length - 1]];
            return out;
        }
    } else if (m !== null) {
        m = pick(m, keep);
    }
    if (rr.length < 2 || (rr[rr.length - 1] - rr[0]) < LOW_R_MIN_PERIODS * period) {
        out.state = 'too-short';
        out.range = rr.length ? [rr[0], rr[rr.length - 1]] : null;
        return out;
    }
    const value = baselineOf(rr, yy, form);
    const base = baselineValues(rr, value, form);
    const about = Float64Array.from(yy, (v, i) => v - base[i]);
    const rms = rmsOf(about);
    const predictedRms = rmsOf(pp);
    const residual = Float64Array.from(about, (v, i) => v - pp[i]);
    Object.assign(out, {
        range: [rr[0], rr[rr.length - 1]],
        baseline: { form, value },
        rms,
        peak: maxAbs(about),
        predictedRms,
        residualRms: rmsOf(residual),
        state: predictedRms > 0 && rms < FILTERED_RATIO * predictedRms ? 'filtered' : 'measured'
    });
    if (m !== null && m.length === rr.length) {
        const modelValue = baselineOf(rr, m, form);
        const modelBase = baselineValues(rr, modelValue, form);
        out.modelRms = rmsOf(Float64Array.from(m, (v, i) => v - modelBase[i]));
        out.baseline.model = modelValue;
        const difference = Float64Array.from(yy, (v, i) => v - m[i]);
        const average = mean(difference);
        out.dataMinusModelRms = rmsOf(Float64Array.from(difference, (v) => v - average));
    }
    return out;
};

// --- one series → the analysis ----------------------------------------------------------

const dSpace = (r, y, space) => Float64Array.from(y, (value, i) => value * (space === 'G' ? r[i] : 1));

const ySpace = (space) => (r, d) => Float64Array.from(d, (value, i) => {
    if (space !== 'G') return value;
    return r[i] !== 0 ? value / r[i] : 0;
});

/** The transform's numbers, the strong peaks and their ripples, the low-r level (analyse_series). */
export const analyseSeries = (rIn, yIn, {
    space, qMin, qMax, dq, rule, window, rFirst = null, model = null, signs = 'positive'
}) => {
    const r = Float64Array.from(rIn);
    const y = Float64Array.from(yIn);
    const n = pyRound((qMax - qMin) / dq) + 1;
    const q = qGrid(qMin, dq, n);
    const weights = transformWeights(q, dq, rule, window);
    const res = resolution(q, weights);
    const d = dSpace(r, y, space);
    const toY = ySpace(space);
    const floor = rFirst === null ? null : rFirst - res.fwhm;
    const peaks = strongPeaks(r, d, q, weights, dq, rule, floor, signs);
    const predictedD = patternOn(r, q, weights, dq, rule, peaks);
    const predicted = toY(r, predictedD);
    const listed = peaks.map((peak) => ({
        r: peak.r,
        sigma: peak.sigma,
        amplitude: peak.amplitude,
        height: toY([peak.r], [peak.amplitude])[0],
        rippleScale: peak.rippleScale,
        fitRms: peak.rms,
        ripples: listedRipples(r, q, weights, dq, rule, peak, res.mainLobe, toY)
    }));
    const scale = y.length ? maxAbs(y) : 0;
    const low = rFirst === null ? null : lowRLevel(r, y, predicted, {
        form: space === 'G' ? 'constant' : 'line', rFirst, fwhm: res.fwhm, period: res.ripplePeriod, scale,
        model: model === null ? null : Float64Array.from(model)
    });
    return {
        space,
        transform: { qMin: q[0], qMax: q[q.length - 1], dQ: dq, nQ: q.length, rule, window },
        resolution: res,
        peaks: listed,
        lowR: low,
        _arrays: { r, y, d, q, weights, predicted, predictedD, peaks }
    };
};

// --- features ----------------------------------------------------------------------------

const boxAverage = (y, half) => {
    const n = y.length;
    const cumulative = new Float64Array(n + 1);
    for (let i = 0; i < n; i += 1) cumulative[i + 1] = cumulative[i] + y[i];
    return Float64Array.from({ length: n }, (_, i) => {
        const lo = Math.max(0, i - half);
        const hi = Math.min(n, i + half + 1);
        return (cumulative[hi] - cumulative[lo]) / (hi - lo);
    });
};

/** The deviation h = d − P − B and the predicted ripples, in y's units (deviation). */
export const deviation = (analysis) => {
    const { r, d, predictedD, peaks, predicted } = analysis._arrays;
    const step = r[1] - r[0];
    const period = analysis.resolution.ripplePeriod;
    const peaksD = peaksOn(r, peaks);
    const rest = Float64Array.from(d, (value, i) => value - peaksD[i] - predictedD[i]);
    const half = Math.max(1, pyRound(0.5 * BACKGROUND_PERIODS * period / step));
    const background = boxAverage(rest, half);
    const h = Float64Array.from(d, (value, i) => value - peaksD[i] - background[i]);
    return { r, h: ySpace(analysis.space)(r, h), predicted };
};

const nearestExtremum = (r, values, at, tolerance, sign = null) => {
    let best = null;
    for (const i of localExtrema(values)) {
        if (Math.abs(r[i] - at) > tolerance) continue;
        if (sign !== null && values[i] * sign <= 0) continue;
        if (best === null || Math.abs(r[i] - at) < Math.abs(r[best] - at)) best = i;
    }
    return best;
};

/** Features in [rMin, rMax]: extrema of h beyond the strong peaks' cores, by size (find_features). */
export const findFeatures = (analysis, rMin, rMax) => {
    const { r, h } = deviation(analysis);
    const lobe = analysis.resolution.mainLobe;
    const { peaks } = analysis;
    const height = peaks.length ? Math.max(...peaks.map((peak) => Math.abs(peak.height))) : 0;
    const found = localExtrema(h).filter((i) => (
        rMin <= r[i] && r[i] <= rMax
        && !peaks.some((peak) => Math.abs(r[i] - peak.r) <= rippleCore(peak, lobe))
        && !(height > 0 && Math.abs(h[i]) < FEATURE_MIN_FRACTION * height)
    ));
    found.sort((a, b) => (Math.abs(h[b]) - Math.abs(h[a])) || (r[a] - r[b]));
    return found.slice(0, MAX_FEATURES).map((i) => r[i]).sort((a, b) => a - b);
};

/** One feature in one dataset (classify_feature). */
export const classifyFeature = (analysis, at, noise = null) => {
    const { r, h, predicted } = deviation(analysis);
    const res = analysis.resolution;
    const tolerance = RIPPLE_POSITION_FRACTION * res.ripplePeriod;
    let nearest = null;
    for (const peak of analysis.peaks) {
        const distance = Math.abs(at - peak.r);
        if (nearest === null || distance < nearest.distance) {
            nearest = { r: peak.r, distance, height: peak.height, sigma: peak.sigma };
        }
    }
    if (nearest !== null && nearest.distance <= tolerance) {
        return {
            r: nearest.r, asked: at, size: nearest.height, predictedHere: interp(nearest.r, r, predicted),
            ripple: null, class: 'resolved-peak', fromPeak: null, belowResolution: false
        };
    }
    const index = nearestExtremum(r, h, at, tolerance);
    const position = index !== null ? r[index] : at;
    const size = index !== null ? h[index] : interp(at, r, h);
    const sign = size >= 0 ? 1 : -1;
    const rippleIndex = nearestExtremum(r, predicted, position, tolerance, sign);
    const ripple = rippleIndex === null ? null : { r: r[rippleIndex], amplitude: predicted[rippleIndex] };
    nearest = null;
    for (const peak of analysis.peaks) {
        const distance = Math.abs(position - peak.r);
        if (nearest === null || distance < nearest.distance) nearest = { r: peak.r, distance };
    }
    let kind;
    if (noise !== null && Math.abs(size) < 2 * noise) kind = 'within-noise';
    else if (ripple === null || Math.abs(ripple.amplitude) < NEGLIGIBLE_RIPPLE_FRACTION * Math.abs(size)) kind = 'no-ripple-here';
    else if (Math.abs(size) > RIPPLE_SIZE_FACTOR * Math.abs(ripple.amplitude)) kind = 'larger-than-ripple';
    else if (Math.abs(size) < Math.abs(ripple.amplitude) / RIPPLE_SIZE_FACTOR) kind = 'smaller-than-ripple';
    else kind = 'ripple-consistent';
    return {
        r: position,
        asked: at,
        size,
        predictedHere: interp(position, r, predicted),
        ripple,
        class: kind,
        fromPeak: nearest,
        belowResolution: Boolean(nearest !== null && nearest.distance < res.fwhm)
    };
};

/** A split or shoulder of spacing |r₂ − r₁| against each dataset's FWHM (split_reading). */
export const splitReading = (analyses, r1, r2) => {
    const spacing = Math.abs(r2 - r1);
    return {
        spacing,
        datasets: analyses.map((a) => ({ fwhm: a.resolution.fwhm, belowResolution: spacing < a.resolution.fwhm }))
    };
};

/** One feature across datasets: a ripple moves with Q_max (across_datasets). */
export const acrossDatasets = (readings, analyses) => {
    const relative = readings.map((reading, k) => {
        const { peaks } = analyses[k];
        const height = peaks.length ? Math.max(...peaks.map((peak) => Math.abs(peak.height))) : 0;
        return height > 0 ? Math.abs(reading.size) / height : 0;
    });
    const largest = relative.length ? Math.max(...relative) : 0;
    const reference = readings.length ? readings[relative.indexOf(largest)] : null;
    const sign = reference === null || reference.size >= 0 ? 1 : -1;
    const present = readings.map((reading, k) => reading.class === 'resolved-peak'
        || (largest > 0 && relative[k] >= PRESENT_FRACTION * largest && reading.size * sign > 0
            && reading.class !== 'within-noise'));
    const ripple = readings.map((reading) => reading.class === 'ripple-consistent');
    const qmaxes = analyses.map((a) => a.transform.qMax);
    const distances = readings.filter((reading) => reading.fromPeak).map((reading) => reading.fromPeak.distance);
    const distance = distances.length ? Math.max(...distances) : 0;
    const tolerance = Math.min(...analyses.map((a) => RIPPLE_POSITION_FRACTION * a.resolution.ripplePeriod));
    const shift = qmaxes.length && Math.max(...qmaxes) > 0
        ? distance * Math.abs(1 - Math.min(...qmaxes) / Math.max(...qmaxes)) : 0;
    const discriminating = shift > tolerance;
    let verdict;
    if (ripple.some(Boolean) && !present.every(Boolean)) verdict = discriminating ? 'moves-with-qmax' : 'not-discriminating';
    else if (present.every(Boolean) && !ripple.some(Boolean)) verdict = 'in-every-dataset';
    else if (ripple.every(Boolean)) verdict = 'ripple-in-every-dataset';
    else verdict = 'mixed';
    return { present, ripple, verdict, discriminating };
};

/** Each feature in every dataset, and across the totals (classify_features). */
export const classifyFeatures = (analyses, features, noise = null, totals = null) => {
    const compared = analyses.map((_, k) => k).filter((k) => totals === null || totals[k]);
    return features.map((at) => {
        const readings = analyses.map((a, k) => classifyFeature(a, at, noise === null ? null : noise[k]));
        const entry = { asked: at, datasets: readings };
        if (compared.length > 1) {
            entry.across = {
                ...acrossDatasets(compared.map((k) => readings[k]), compared.map((k) => analyses[k])),
                compared
            };
        }
        return entry;
    });
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

/** D(r), T(r) and G(r)P carry the factor r. */
export const isDLike = (fitType) => {
    const text = (fitType || '').replace(/ /g, '').toUpperCase();
    return text.startsWith('D(R)') || text.startsWith('T(R)') || text.startsWith('G(R)P');
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

/** The first interatomic distance (closest_approach): `columns` [r, partial…]. */
export const closestApproach = (columns, minimumDistances) => {
    if (columns && columns.length >= 2) {
        const r = columns[0];
        for (let i = 0; i < r.length; i += 1) {
            for (let c = 1; c < columns.length; c += 1) {
                // np.nan_to_num(…) > 0: NaN is 0, +∞ counts.
                if (columns[c][i] > 0) return { value: r[i], source: 'partials' };
            }
        }
    }
    const positive = (minimumDistances || []).filter((value) => value > 0);
    return positive.length ? { value: Math.min(...positive), source: 'minimum-distances' } : null;
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
export const isTotal = (dataset) => {
    const coefficients = (dataset || {}).coefficients || [];
    return !(coefficients.length && coefficients.filter((value) => value !== 0).length === 1);
};

/** −Σ w_ij of the block's coefficients, a G(r)'s value below the closest approach. */
export const expectedBaseline = (dataset, space) => {
    const coefficients = (dataset || {}).coefficients;
    if (space !== 'G' || !coefficients || !coefficients.length) return null;
    return -sum(coefficients);
};

const blockSummary = (dataset) => (dataset
    ? Object.fromEntries(['block', 'file', 'fitType', 'realSpaceFit', 'convolve', 'convolveQ', 'rCutoff', 'qMin', 'qMax']
        .map((key) => [key, dataset[key] ?? null]))
    : null);

/** One fitted dataset of a run (assess_dataset). */
export const assessDataset = ({
    name, kind, number, settings, x: xIn, calc, data: dataIn, fq = null, inputColumns = null, dataTransform = null,
    override = null, rFirst = null, processing = null
}) => {
    const dataset = nthBlock(settings, BLOCKS[kind] || '', number);
    const real = ['xpdf', 'npdf', 'xray-real'].includes(kind);
    const payload = {
        status: 'unavailable', reason: null, file: name, kind, number, total: isTotal(dataset),
        block: blockSummary(dataset), fit: fitSpaceRisk(real ? kind : 'reciprocal', dataset, dataTransform),
        transform: null, resolution: null, peaks: [], lowR: null, processing, space: null, rFirst, _analysis: null
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
    let model = null;
    let rr;
    let yy;
    let qMin;
    let qMax;
    let dq;
    let rule;
    let window;
    let source;
    let windowKnown;
    let transformExtra = {};
    let space;
    if (kind === 'xpdf') {
        if (fq === null) return { ...payload, reason: 'no-fq-file' };
        const [q, fCalc, fExp] = fq.map((column) => Float64Array.from(column));
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
        const step = x[1] - x[0];
        const count = Math.floor(x[0] / step - 1e-9);
        const lowR = [];
        for (let k = count; k >= 1; k -= 1) {
            const value = x[0] - step * k;
            if (value > 0) lowR.push(value);
        }
        const coefficients = Float64Array.from(q, (value) => dq * value);
        const transformOf = (f) => lowR.map((at) => {
            let total = 0;
            for (let k = 0; k < q.length; k += 1) total += Math.sin(at * q[k]) * (coefficients[k] * f[k]);
            return total;
        });
        const shifted = Float64Array.from(fExp, (value) => value + reproduction.offset);
        rr = Float64Array.from([...lowR, ...x]);
        yy = Float64Array.from([...transformOf(shifted), ...data]);
        model = Float64Array.from([...transformOf(fCalc), ...Float64Array.from(calc)]);
        qMin = q[0];
        qMax = q[q.length - 1];
        rule = 'rectangle';
        window = 'none';
        source = 'run-fq';
        windowKnown = true;
        transformExtra = { reproduction: { residual: reproduction.residual, offset: reproduction.offset, tolerance: RMC_FT_TOLERANCE } };
        space = 'D';
    } else {
        space = isDLike(dataset ? dataset.fitType : null) ? 'D' : 'G';
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
        rr = x;
        yy = data;
        if (inputColumns !== null) {
            const mapped = mapInput(x, data, inputColumns);
            if (mapped !== null) {
                rr = mapped.r;
                yy = mapped.y;
                transformExtra = { input: { scale: mapped.scale, offset: mapped.offset, residual: mapped.residual } };
            }
        }
    }
    payload.space = space;
    const finite = Array.from(rr, (value, i) => Number.isFinite(value) && Number.isFinite(yy[i]));
    rr = pick(rr, finite);
    yy = pick(yy, finite);
    if (model !== null) model = pick(model, finite);
    if (!uniform(rr)) return { ...payload, status: 'failed', reason: 'r-grid-not-uniform' };
    const signs = dataset && (dataset.coefficients || []).some((value) => value < 0) ? 'both' : 'positive';
    const analysis = analyseSeries(rr, yy, {
        space, qMin, qMax, dq, rule, window, rFirst: rFirst === null ? null : rFirst.value, model, signs
    });
    const expected = expectedBaseline(dataset, space);
    if (expected !== null && analysis.lowR && analysis.lowR.baseline) analysis.lowR.baseline.expected = expected;
    return {
        ...payload,
        status: 'assessed',
        transform: { ...analysis.transform, source, windowKnown, ...transformExtra },
        resolution: analysis.resolution,
        peaks: analysis.peaks,
        lowR: analysis.lowR,
        range: [rr[0], rr[rr.length - 1]],
        _analysis: analysis
    };
};

/** A dataset's noise from its measured low-r level (noise_level). */
export const noiseLevel = (payload) => {
    const low = payload.lowR;
    if (!low || low.state !== 'measured') return null;
    if (low.dataMinusModelRms !== null && low.dataMinusModelRms !== undefined) return low.dataMinusModelRms;
    return low.residualRms ?? null;
};

/** Features of the box's partial g(r) against the risky datasets' ripples (pair_reading). */
export const pairReading = (r, g, features, payloads) => {
    const analysed = payloads.filter((p) => p.status === 'assessed' && p._analysis);
    const resolved = payloads.filter((p) => p.resolution);
    const fwhm = resolved.length ? Math.min(...resolved.map((p) => p.resolution.fwhm)) : null;
    const period = resolved.length ? Math.min(...resolved.map((p) => p.resolution.ripplePeriod)) : null;
    const tolerance = period ? RIPPLE_POSITION_FRACTION * period : 0.05;
    let largest = -Infinity;
    for (const value of g) if (Number.isFinite(value) && value > largest) largest = value;
    if (!g.length) largest = 0;
    const maxima = localExtrema(g).filter((i) => g[i] > g[i - 1] && g[i] >= 0.15 * largest);
    return features.map((at) => {
        const near = maxima.filter((i) => Math.abs(r[i] - at) <= tolerance);
        const position = near.length
            ? r[near.reduce((best, i) => (Math.abs(r[i] - at) < Math.abs(r[best] - at) ? i : best), near[0])]
            : at;
        const others = maxima.filter((i) => Math.abs(r[i] - position) > tolerance).map((i) => r[i]);
        const neighbour = others.length
            ? others.reduce((best, value) => (Math.abs(value - position) < Math.abs(best - position) ? value : best), others[0])
            : null;
        const spacing = neighbour !== null ? Math.abs(neighbour - position) : null;
        const ripples = [];
        for (const p of analysed) {
            if (!['high', 'unknown'].includes(p.fit.risk)) continue;
            const arrays = p._analysis._arrays;
            const index = nearestExtremum(arrays.r, arrays.predicted, position, tolerance);
            if (index === null) continue;
            const height = p.peaks.length ? Math.max(...p.peaks.map((peak) => Math.abs(peak.height))) : 0;
            ripples.push({ file: p.file, r: arrays.r[index], relative: height > 0 ? arrays.predicted[index] / height : null });
        }
        return {
            asked: at, r: position, inPartial: near.length > 0, neighbour, spacing,
            belowResolution: spacing !== null && fwhm !== null && spacing < fwhm,
            atRipple: ripples
        };
    });
};

/** The payload without its internal analysis. */
export const publicPayload = (payload) => Object.fromEntries(Object.entries(payload).filter(([key]) => key !== '_analysis'));

/** The run's datasets and the caller's features read in them (assess_run). */
export const assessRun = (payloads, { features = null, rRange = null, split = null, pair = null } = {}) => {
    const analysed = payloads.filter((p) => p.status === 'assessed' && p._analysis);
    const analyses = analysed.map((p) => p._analysis);
    let found = null;
    let wanted = features;
    if (wanted === null && rRange !== null && analyses.length) {
        const candidates = [];
        for (const analysis of analyses) {
            const height = analysis.peaks.length ? Math.max(...analysis.peaks.map((peak) => Math.abs(peak.height))) : 0;
            let dev = null;
            for (const at of findFeatures(analysis, rRange[0], rRange[1])) {
                if (dev === null) dev = deviation(analysis);
                const size = Math.abs(interp(at, dev.r, dev.h));
                candidates.push([height > 0 ? size / height : 0, at]);
            }
        }
        const tolerance = Math.min(...analyses.map((a) => RIPPLE_POSITION_FRACTION * a.resolution.ripplePeriod));
        const merged = [];
        for (const [relative, at] of candidates.sort((a, b) => (b[0] - a[0]) || (a[1] - b[1]))) {
            if (merged.every(([, other]) => Math.abs(at - other) > tolerance)) merged.push([relative, at]);
        }
        found = merged.slice(0, MAX_FEATURES).map(([, at]) => at).sort((a, b) => a - b);
        wanted = found;
    }
    const out = {
        datasets: payloads.map(publicPayload),
        features: null,
        found: found !== null,
        split: null,
        pair: null,
        rules: RULES
    };
    if (wanted && wanted.length) {
        if (pair !== null) {
            out.pair = { label: pair.label, features: pairReading(pair.r, pair.g, wanted, payloads) };
        } else if (analysed.length) {
            const readings = classifyFeatures(analyses, wanted, analysed.map(noiseLevel), analysed.map((p) => p.total !== false));
            for (const entry of readings) entry.datasets.forEach((reading, k) => { reading.file = analysed[k].file; });
            out.features = readings;
        }
    }
    if (split !== null && split.length === 2) {
        const resolved = payloads.filter((p) => p.resolution);
        const reading = splitReading(resolved.map((p) => ({ resolution: p.resolution })), split[0], split[1]);
        reading.datasets.forEach((entry, k) => { entry.file = resolved[k].file; });
        out.split = reading;
    }
    return out;
};
