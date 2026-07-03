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

export const detectDivergence = (values, epsilon = WINDOW_DELTA_EPSILON) => (
    Array.isArray(values) && values.length >= 2 && windowDelta(values) > epsilon
);

export const detectStall = (values, {
    epsilon = WINDOW_DELTA_EPSILON,
    minTotalDrop = MIN_TOTAL_DROP
} = {}) => {
    if (!Array.isArray(values) || values.length < 2) return false;
    if (Math.abs(windowDelta(values)) > epsilon) return false;
    return values[0] - values[values.length - 1] < minTotalDrop;
};

// Classify the run: 'improving' | 'converged' | 'stalled' | 'diverging' | 'unknown'.
export const classifyConvergence = (values, {
    epsilon = WINDOW_DELTA_EPSILON,
    minTotalDrop = MIN_TOTAL_DROP
} = {}) => {
    if (!Array.isArray(values) || values.length < 2) return 'unknown';
    const delta = windowDelta(values);
    if (delta > epsilon) return 'diverging';
    if (delta < -epsilon) return 'improving';
    return values[0] - values[values.length - 1] < minTotalDrop ? 'stalled' : 'converged';
};

// Recent-window statistics passed to the watchdog LLM prompt: small, rounded,
// and self-describing — never the full history.
export const watchdogStats = (values) => {
    const stats = seriesStats(values);
    if (!stats) return null;
    return {
        n_steps: stats.nSteps,
        first: Number(stats.first.toPrecision(3)),
        last: Number(stats.last.toPrecision(3)),
        min: Number(stats.min.toPrecision(3)),
        recent_window_delta: Number(windowDelta(values).toPrecision(2))
    };
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
