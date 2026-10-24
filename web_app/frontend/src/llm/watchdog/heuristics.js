// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { seriesStats } from '../context/runContext';

// Convergence heuristics over the ln(chi^2) histories. These pure functions are
// the watchdog's source of truth — the LLM only narrates on top of them — so
// the badge keeps working with no model connected at all. The Copilot's
// assess_convergence, its run context (convergence.watchdog_status) and the
// answer check on "converged" read the same verdict.
//
// Every chi^2 term of the .log and their total are judged (the user's rule,
// 2026-10-24): a run is converged only when no term and not the total is
// still falling beyond its noise over the last tenth of the run. "Falling
// beyond its noise" (seriesTrend): the least-squares line through the
// series' last tenth (at least TREND_MIN_ROWS rows, never back past the last
// restart) drops by more than TREND_NOISE_FACTOR times the scatter of the
// points about that line. k = 3, the app's 3σ bar: for a flat series of
// independent noise of sd σ, the fitted change across m rows has sd
// σ·√(12(m−1)/(m(m+1))) — at most 0.99σ from m = 10 on, 0.23σ over the Demo's
// 227-row window — so a flat series passes for falling with probability
// ≤ 0.12 % with σ known, 0.8 % with σ estimated from the 10 rows of the
// shortest window (Student's t, 8 degrees of freedom), and never over 227
// rows. The rule compares the change with the scatter of the points, not a t
// statistic that grows with the row count, so a long window does not make a
// negligible drift "significant". Noise correlated from row to row wanders
// like a drift over a window not much longer than its correlation (lag-1
// ρ = 0.9: 0.25 % falsely falling over 227 rows, 8.6 % over 50, simulated):
// there the rule errs towards "still moving", never towards "converged".
// heuristics.test.js pins these rates.

// The window: the last tenth of the run …
export const TREND_WINDOW_FRACTION = 0.1;
// … and at least this many log rows (a shorter history: all of it).
export const TREND_MIN_ROWS = 10;
// Still falling (rising) when the fitted change across the window exceeds
// this many times the scatter about the fit (residual sd).
export const TREND_NOISE_FACTOR = 3;
// A flat run whose total and terms never dropped by this much in ln(chi^2)
// (0.1 ≈ a 10 % chi^2 improvement) is "stalled" rather than "converged".
export const MIN_TOTAL_DROP = 0.1;
// A fitted change below this (ln units) is round-off, never a trend.
const ROUND_OFF = 1e-12;
// The chart's clamp (browserData.js, plots.chi_history_ln): ln(1e-12), where a
// term whose chi^2 is 0 in every row (a constraint not in use) sits.
const ZERO_LN = Math.log(1e-12);

export const TRENDS = Object.freeze(['falling', 'rising', 'flat', 'unknown']);

const finiteCount = (values) => values.reduce((count, value) => count + (Number.isFinite(value) ? 1 : 0), 0);

/** The rows of the window over a history of n rows (TREND_WINDOW_FRACTION, TREND_MIN_ROWS). */
export const trendWindowRows = (n) => Math.max(Math.ceil(n * TREND_WINDOW_FRACTION), Math.min(n, TREND_MIN_ROWS));

/**
 * The trend of one ln(chi^2) series over the last tenth of the run: the last
 * trendWindowRows(n) rows of its n, from row `start` (the last restart: a
 * restart recomputes chi^2, so rows across it are not one series) at the
 * earliest; its finite values fitted by least squares against the row index.
 * Returns { trend, change, noiseSd, rows }: `change` the fitted line's change
 * across the rows it spans (negative falls), `noiseSd` the scatter about it
 * (residual sd, rows − 2 degrees of freedom), `rows` the finite values used;
 * `trend` falling when −change > TREND_NOISE_FACTOR × noiseSd, rising when
 * change > it, flat otherwise (a change below round-off is flat), unknown with
 * fewer than three finite values in the window (change and noiseSd null).
 */
export const seriesTrend = (values, { start = 0 } = {}) => {
    const unknown = { trend: 'unknown', change: null, noiseSd: null, rows: 0 };
    if (!Array.isArray(values) || !values.length) return unknown;
    const n = values.length;
    const from = Math.max(0, Math.min(start, n - 1), n - trendWindowRows(n));
    const x = [];
    const y = [];
    for (let i = from; i < n; i += 1) {
        if (Number.isFinite(values[i])) {
            x.push(i);
            y.push(values[i]);
        }
    }
    const m = x.length;
    if (m < 3) return { ...unknown, rows: m };
    const meanX = x.reduce((sum, value) => sum + value, 0) / m;
    const meanY = y.reduce((sum, value) => sum + value, 0) / m;
    let sxy = 0;
    let sxx = 0;
    for (let k = 0; k < m; k += 1) {
        sxy += (x[k] - meanX) * (y[k] - meanY);
        sxx += (x[k] - meanX) ** 2;
    }
    const slope = sxx > 0 ? sxy / sxx : 0;
    let rss = 0;
    for (let k = 0; k < m; k += 1) rss += (y[k] - meanY - slope * (x[k] - meanX)) ** 2;
    // A scatter below round-off is none (a noise-free or constant series).
    const sd = Math.sqrt(rss / (m - 2));
    const noiseSd = sd < ROUND_OFF ? 0 : sd;
    const change = slope * (x[m - 1] - x[0]);
    const beyond = Math.abs(change) > ROUND_OFF && Math.abs(change) > TREND_NOISE_FACTOR * noiseSd;
    let trend = 'flat';
    if (beyond) trend = change < 0 ? 'falling' : 'rising';
    return { trend, change, noiseSd, rows: m };
};

// The first and last finite values, the minimum and the drop first − last.
const extent = (values) => {
    const finite = values.filter(Number.isFinite);
    if (!finite.length) return { first: null, last: null, min: null, drop: null };
    let min = Infinity;
    finite.forEach((value) => { if (value < min) min = value; });
    return { first: finite[0], last: finite[finite.length - 1], min, drop: finite[0] - finite[finite.length - 1] };
};

/** A term whose chi^2 is 0 in every row (all at the clamp): a constraint not in use. */
export const isZeroTerm = (values) => Array.isArray(values) && finiteCount(values) > 0
    && values.every((value) => !Number.isFinite(value) || value <= ZERO_LN + 1e-9);

/** Header names made unique where they repeat: Curvature ×6 → Curvature #1 … #6. */
export const termLabels = (names) => {
    const total = new Map();
    names.forEach((name) => total.set(name, (total.get(name) || 0) + 1));
    const seen = new Map();
    return names.map((name, index) => {
        const label = name || `column ${index + 1 + 3}`;
        if (total.get(name) < 2) return label;
        seen.set(name, (seen.get(name) || 0) + 1);
        return `${label} #${seen.get(name)}`;
    });
};

const CHART_COLUMN_NAME = 'last .log column';

/**
 * The run's chi^2 series as the convergence rule judges them, from the
 * r-value plot data (browserData.js chiTermsPlotData / combineChiTerms, or
 * Flask's chiTerms): every term, labelled uniquely (termLabels), their total,
 * and the row the last log starts at (`start`, from the segments). Without
 * per-term histories (logs naming different terms) the charted column stands
 * for the terms beside the total; without terms at all (a log whose line 1
 * names no columns) the charted column alone, with no total. `column` is the
 * charted column's name. Null without a history.
 */
export const chiSeries = (plotData) => {
    const chart = plotData?.series?.[0]?.y;
    if (!Array.isArray(chart) || !chart.length) return null;
    const column = plotData.chiColumn || null;
    const terms = plotData.chiTerms || null;
    const charted = [{ name: column || CHART_COLUMN_NAME, values: chart }];
    const segments = Array.isArray(terms?.segments) ? terms.segments : [];
    const rows = segments.reduce((sum, count) => sum + count, 0);
    const start = segments.length > 1 && rows === chart.length ? rows - segments[segments.length - 1] : 0;
    const total = Array.isArray(terms?.total) && terms.total.length === chart.length ? terms.total : null;
    if (total && Array.isArray(terms.names) && Array.isArray(terms.series) && terms.series.length === terms.names.length
        && terms.series.every((series) => Array.isArray(series) && series.length === chart.length)) {
        const labels = termLabels(terms.names);
        return { terms: terms.series.map((values, k) => ({ name: labels[k], values })), total, start, column };
    }
    return { terms: charted, total, start, column };
};

/**
 * The convergence of a run's chi^2 series ({ terms: [{ name, values }],
 * total, start }, as chiSeries gives them; total may be null): every term's
 * and the total's trend over the last tenth (seriesTrend), and the run's
 * status, in this order —
 *   unknown    no series has three finite values in its window;
 *   diverging  a series' latest value is non-finite (a blown-up run), or the
 *              total rises beyond its noise;
 *   improving  a term or the total still falls beyond its noise;
 *   diverging  a term rises beyond its noise while none falls;
 *   stalled    every series flat, and none ever dropped by MIN_TOTAL_DROP
 *              from its first value to its last;
 *   converged  every series flat after such a drop.
 * Terms whose chi^2 is 0 in every row (isZeroTerm) are listed apart, judged
 * by no rule. Returns { status, windowRows, terms: [{ name, trend, change,
 * noiseSd, rows, first, last, min, drop }], total (the same, or null),
 * zeroTerms: [names], falling, rising: [names, 'total' for the total] }.
 */
export const judgeConvergence = (series) => {
    const empty = { status: 'unknown', windowRows: 0, terms: [], total: null, zeroTerms: [], falling: [], rising: [] };
    if (!series || !Array.isArray(series.terms) || !series.terms.length) return empty;
    const start = series.start || 0;
    const judged = (values) => ({ ...seriesTrend(values, { start }), ...extent(values) });
    const zeroTerms = series.terms.filter((term) => isZeroTerm(term.values)).map((term) => term.name);
    const terms = series.terms.filter((term) => !isZeroTerm(term.values)).map((term) => ({ name: term.name, ...judged(term.values) }));
    const total = Array.isArray(series.total) ? judged(series.total) : null;
    const all = [...terms, ...(total ? [{ ...total, name: 'total' }] : [])];
    const reference = series.total || series.terms[0].values;
    const windowRows = Math.min(trendWindowRows(reference.length), reference.length - Math.min(start, reference.length - 1));
    const falling = all.filter((entry) => entry.trend === 'falling').map((entry) => entry.name);
    const rising = all.filter((entry) => entry.trend === 'rising').map((entry) => entry.name);
    const result = { windowRows, terms, total, zeroTerms, falling, rising };
    const latestNonFinite = [...series.terms.map((term) => term.values), ...(series.total ? [series.total] : [])]
        .some((values) => values.length && !Number.isFinite(values[values.length - 1]));
    if (latestNonFinite) return { ...result, status: 'diverging' };
    if (!all.length || all.every((entry) => entry.trend === 'unknown')) return { ...result, status: 'unknown' };
    if (total?.trend === 'rising') return { ...result, status: 'diverging' };
    if (falling.length) return { ...result, status: 'improving' };
    if (rising.length) return { ...result, status: 'diverging' };
    const dropped = all.some((entry) => Number.isFinite(entry.drop) && entry.drop >= MIN_TOTAL_DROP);
    return { ...result, status: dropped ? 'converged' : 'stalled' };
};

/** The run's convergence from its r-value plot data (chiSeries → judgeConvergence). */
export const judgeRun = (plotData) => judgeConvergence(chiSeries(plotData));

/**
 * One ln(chi^2) series classified alone, by the same rule:
 * 'improving' | 'converged' | 'stalled' | 'diverging' | 'unknown'.
 */
export const classifyConvergence = (values) => (
    Array.isArray(values) && values.length
        ? judgeConvergence({ terms: [{ name: CHART_COLUMN_NAME, values }], total: null, start: 0 }).status
        : 'unknown'
);

/** A run classified from its r-value plot data: every term and their total. */
export const classifyRun = (plotData) => judgeRun(plotData).status;

/** One series alone is diverging (classifyConvergence). */
export const detectDivergence = (values) => classifyConvergence(values) === 'diverging';

/** One series alone is stalled (classifyConvergence). */
export const detectStall = (values) => classifyConvergence(values) === 'stalled';

// Statistics passed to the watchdog LLM prompt: small, rounded, and
// self-describing — never the full history.
const round = (value, digits) => (Number.isFinite(value) ? Number(value.toPrecision(digits)) : null);

/**
 * One series' summary for the watchdog prompt (watchdogRunStats builds the
 * run's); `start` is the last restart's row, as seriesTrend takes it.
 */
export const watchdogStats = (values, { start = 0 } = {}) => {
    const stats = seriesStats(values);
    if (!stats) return null;
    const { change } = seriesTrend(values, { start });
    const summary = {
        n_steps: stats.nSteps,
        first: round(stats.first, 3),
        last: round(stats.last, 3),
        min: round(stats.min, 3),
        change_last_10pct: round(change, 2)
    };
    if (stats.nonFiniteSteps) summary.non_finite_steps = stats.nonFiniteSteps;
    return summary;
};

/**
 * The run's summary for the watchdog prompt: the total's (or, without one,
 * the charted column's) watchdogStats, and the terms still falling or rising
 * beyond their noise over the last tenth. Null without a history.
 */
export const watchdogRunStats = (plotData) => {
    const series = chiSeries(plotData);
    if (!series) return null;
    const judgement = judgeConvergence(series);
    const summary = watchdogStats(series.total || series.terms[0].values, { start: series.start });
    if (!summary) return null;
    if (series.total) summary.of = 'total';
    if (judgement.falling.length) summary.still_falling = judgement.falling;
    if (judgement.rising.length) summary.rising = judgement.rising;
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
