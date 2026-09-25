// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { recentSlope, seriesStats } from '../context/runContext';

// Convergence heuristics over the ln(chi^2) history. These pure functions are
// the watchdog's source of truth — the LLM only narrates on top of them — so
// the badge keeps working with no model connected at all.

// Slope thresholds are expressed as the predicted ln(chi^2) change across the
// recent window (slope × window length): additive shifts in ln units are
// relative changes in chi^2, so 0.02 ≈ a 2% move regardless of magnitude.
const WINDOW_DELTA_EPSILON = 0.02;

// A run whose ln(chi^2) never dropped by at least this much is "stalled"
// rather than "converged" when it goes flat (0.1 ≈ a 10% chi^2 improvement).
const MIN_TOTAL_DROP = 0.1;

const recentWindowLength = (values) => (
    Math.max(Math.min(values.length, 5), Math.ceil(values.length * 0.2))
);

const windowDelta = (values) => recentSlope(values) * (recentWindowLength(values) - 1);

// A log row whose chi^2 is non-finite (NaN / Inf / Fortran overflow; JSON null
// from Flask) means the run produced non-finite values. The parsers keep such
// rows (they used to be dropped in the browser, which let a blown-up NaN tail
// read as "improving" on the last finite points), so: a non-finite LATEST value
// is divergence; earlier ones are ignored for the trend.
const lastIsNonFinite = (values) => !Number.isFinite(values[values.length - 1]);
const finiteOnly = (values) => values.filter(Number.isFinite);

export const detectDivergence = (values, epsilon = WINDOW_DELTA_EPSILON) => {
    if (!Array.isArray(values) || values.length < 2) return false;
    if (lastIsNonFinite(values)) return true;
    const finite = finiteOnly(values);
    return finite.length >= 2 && windowDelta(finite) > epsilon;
};

export const detectStall = (values, {
    epsilon = WINDOW_DELTA_EPSILON,
    minTotalDrop = MIN_TOTAL_DROP
} = {}) => {
    if (!Array.isArray(values) || values.length < 2 || lastIsNonFinite(values)) return false;
    const finite = finiteOnly(values);
    if (finite.length < 2) return false;
    if (Math.abs(windowDelta(finite)) > epsilon) return false;
    return finite[0] - finite[finite.length - 1] < minTotalDrop;
};

// Classify the run: 'improving' | 'converged' | 'stalled' | 'diverging' | 'unknown'.
export const classifyConvergence = (values, {
    epsilon = WINDOW_DELTA_EPSILON,
    minTotalDrop = MIN_TOTAL_DROP
} = {}) => {
    if (!Array.isArray(values) || values.length < 2) return 'unknown';
    if (lastIsNonFinite(values)) return 'diverging';
    const finite = finiteOnly(values);
    if (finite.length < 2) return 'unknown';
    const delta = windowDelta(finite);
    if (delta > epsilon) return 'diverging';
    if (delta < -epsilon) return 'improving';
    return finite[0] - finite[finite.length - 1] < minTotalDrop ? 'stalled' : 'converged';
};

// Recent-window statistics passed to the watchdog LLM prompt: small, rounded,
// and self-describing — never the full history.
const round = (value, digits) => (Number.isFinite(value) ? Number(value.toPrecision(digits)) : null);

export const watchdogStats = (values) => {
    const stats = seriesStats(values);
    if (!stats) return null;
    const finite = finiteOnly(values);
    const summary = {
        n_steps: stats.nSteps,
        first: round(stats.first, 3),
        last: round(stats.last, 3),
        min: round(stats.min, 3),
        recent_window_delta: finite.length >= 2 ? round(windowDelta(finite), 2) : null
    };
    if (stats.nonFiniteSteps) summary.non_finite_steps = stats.nonFiniteSteps;
    return summary;
};

// Has the history changed enough since the last LLM call to justify another?
// True on ≥ stepDelta new points, or when the last value moved by more than
// relativeDelta in chi^2 terms (ln values make relative change additive).
export const significantChange = (prevStats, nextStats, {
    relativeDelta = 0.02,
    stepDelta = 200
} = {}) => {
    if (!nextStats) return false;
    if (!prevStats) return true;
    if ((nextStats.n_steps || 0) - (prevStats.n_steps || 0) >= stepDelta) return true;
    if (!Number.isFinite(prevStats.last) || !Number.isFinite(nextStats.last)) {
        return prevStats.last !== nextStats.last;
    }
    return Math.abs(nextStats.last - prevStats.last) >= Math.log(1 + relativeDelta);
};
