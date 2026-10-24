// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import {
    MIN_TOTAL_DROP,
    TREND_MIN_ROWS,
    TREND_NOISE_FACTOR,
    TREND_WINDOW_FRACTION,
    chiSeries,
    classifyConvergence,
    classifyRun,
    detectDivergence,
    detectStall,
    judgeConvergence,
    seriesTrend,
    significantChange,
    termLabels,
    trendWindowRows,
    watchdogRunStats,
    watchdogStats
} from '../watchdog/heuristics';

// A seeded normal deviate (mulberry32 + Box–Muller): the noise of the known answers.
const noise = (seed) => {
    let state = seed >>> 0;
    const uniform = () => {
        state = (state + 0x6D2B79F5) | 0;
        let t = Math.imul(state ^ (state >>> 15), 1 | state);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
    return () => Math.sqrt(-2 * Math.log(uniform() || 1e-300)) * Math.cos(2 * Math.PI * uniform());
};

// Synthetic ln(chi^2) histories for each regime.
const improving = Array.from({ length: 200 }, (_, index) => 2.3 - index * 0.01);
const converged = [
    ...Array.from({ length: 150 }, (_, index) => 2.3 - index * 0.01),
    ...Array.from({ length: 150 }, () => 0.8)
];
// Flat from the start: noise about one level, never a 10 % drop.
const flatNoise = noise(7);
const stalled = Array.from({ length: 200 }, () => 2.3 + 0.005 * flatNoise());
const diverging = [
    ...Array.from({ length: 100 }, (_, index) => 2.3 - index * 0.01),
    ...Array.from({ length: 100 }, (_, index) => 1.3 + index * 0.01)
];

describe('the rule: a series still falls when its last-tenth fit drops by more than 3 × its noise', () => {
    it('k = 3 over the last tenth of the run, at least 10 rows', () => {
        expect([TREND_NOISE_FACTOR, TREND_WINDOW_FRACTION, TREND_MIN_ROWS, MIN_TOTAL_DROP]).toEqual([3, 0.1, 10, 0.1]);
        expect([trendWindowRows(5), trendWindowRows(50), trendWindowRows(2262)]).toEqual([5, 10, 227]);
    });

    it('a flat noisy history after a drop → converged', () => {
        const draw = noise(11);
        const values = [
            ...Array.from({ length: 100 }, (_, t) => 3 - 0.02 * t + 0.01 * draw()),
            ...Array.from({ length: 200 }, () => 1 + 0.01 * draw())
        ];
        const trend = seriesTrend(values);
        expect(trend.trend).toBe('flat');
        expect(trend.rows).toBe(30);
        expect(Math.abs(trend.change)).toBeLessThan(TREND_NOISE_FACTOR * trend.noiseSd);
        expect(classifyConvergence(values)).toBe('converged');
    });

    it('a slowly decaying history → still improving, where the old fixed 2 % threshold called it flat', () => {
        // ln χ² = 1 + 0.3 e^(−t/1000) + noise 0.001: −0.0065 over the last
        // tenth (30 rows; its fit's sd from the noise 0.0006), 6 × its noise;
        // the old rule's window (the last 20 %) changed by −0.013 < 0.02 and
        // read it as flat ("stalled": 0.078 dropped over the run).
        const draw = noise(13);
        const values = Array.from({ length: 300 }, (_, t) => 1 + 0.3 * Math.exp(-t / 1000) + 0.001 * draw());
        const trend = seriesTrend(values);
        expect(trend.trend).toBe('falling');
        expect(trend.change).toBeGreaterThan(-0.0085);
        expect(trend.change).toBeLessThan(-0.0045);
        expect(-trend.change / trend.noiseSd).toBeGreaterThan(TREND_NOISE_FACTOR);
        expect(classifyConvergence(values)).toBe('improving');
    });

    it('a noise-free decline, however slow, still falls (its noise is zero)', () => {
        const values = Array.from({ length: 200 }, (_, index) => 2.3 - 0.0001 * index);
        expect(seriesTrend(values).trend).toBe('falling');
        expect(classifyConvergence(values)).toBe('improving');
        // A constant: round-off is never a trend.
        expect(seriesTrend(Array(50).fill(Math.log(0.123))).trend).toBe('flat');
    });

    it('one term decaying while the total looks flat → not converged (the old rule read the charted column only)', () => {
        // A small reciprocal-space term still falls (0.01 → 0.006 over the run)
        // while a large flat constraint term, charted last, dominates the total.
        const draw = noise(17);
        const n = 400;
        const small = Array.from({ length: n }, (_, t) => 0.01 * Math.exp(-t / 800) * Math.exp(0.002 * draw()));
        const large = Array.from({ length: n }, (_, t) => (t < 60 ? 20 - 0.15 * t : 11) * Math.exp(0.01 * draw()));
        const ln = (values) => values.map((value) => Math.log(value));
        const plotData = {
            chiColumn: 'Curvature',
            series: [{ y: ln(large) }],
            chiTerms: { names: ['F(Q)_1', 'Curvature'], series: [ln(small), ln(large)], total: ln(small.map((value, t) => value + large[t])), segments: [n] }
        };
        const judgement = judgeConvergence(chiSeries(plotData));
        expect(judgement.terms.map((term) => [term.name, term.trend])).toEqual([['F(Q)_1', 'falling'], ['Curvature', 'flat']]);
        expect(judgement.total.trend).toBe('flat');
        expect(judgement.status).toBe('improving');
        expect(judgement.falling).toEqual(['F(Q)_1']);
        // The charted column alone (the badge before 2026-10-24) looks converged.
        expect(classifyConvergence(plotData.series[0].y)).toBe('converged');
        expect(classifyRun(plotData)).toBe('improving');
    });

    it('never fits across the last restart: a jump at a restart inside the last tenth is no trend', () => {
        // 285 rows at ln χ² 2, a restart, 15 rows at 1 (both noisy): the last
        // tenth (30 rows) holds the restart, so only the 15 rows after it count.
        const draw = noise(19);
        const values = [...Array.from({ length: 285 }, () => 2 + 0.01 * draw()), ...Array.from({ length: 15 }, () => 1 + 0.01 * draw())];
        const plotData = { series: [{ y: values }], chiTerms: { names: ['X_ray_(R)1'], series: [values], total: values, segments: [285, 15] } };
        const judgement = judgeConvergence(chiSeries(plotData));
        expect(judgement.windowRows).toBe(15);
        expect(judgement.total.trend).toBe('flat');
        expect(judgement.status).toBe('converged');
        // Across the restart the same rows would read as falling.
        expect(seriesTrend(values).trend).toBe('falling');
    });

    it('a flat series passes for falling about as often as the 3σ bar says (seeded noise)', () => {
        const rate = (rows, rho, trials, seed) => {
            const draw = noise(seed);
            let falling = 0;
            for (let trial = 0; trial < trials; trial += 1) {
                let e = draw();
                const values = Array.from({ length: rows }, () => {
                    e = rho * e + Math.sqrt(1 - rho * rho) * draw();
                    return 1 + 0.01 * e;
                });
                if (seriesTrend(values).trend === 'falling') falling += 1;
            }
            return falling / trials;
        };
        // 10-row window (Student's t, 8 df: 0.82 %) and the Demo's 227 rows.
        expect(rate(100, 0, 4000, 23)).toBeLessThan(0.015);
        expect(rate(2270, 0, 300, 29)).toBe(0);
        // Correlated noise (lag-1 ρ = 0.9) over 227 rows stays rare.
        expect(rate(2270, 0.9, 1000, 31)).toBeLessThan(0.01);
    });
});

describe('classifyConvergence (one series alone)', () => {
    it('labels a steadily dropping series as improving', () => {
        expect(classifyConvergence(improving)).toBe('improving');
    });

    it('labels flat-at-minimum after a real drop as converged', () => {
        expect(classifyConvergence(converged)).toBe('converged');
    });

    it('labels flat-from-the-start as stalled', () => {
        expect(classifyConvergence(stalled)).toBe('stalled');
    });

    it('labels a rising tail as diverging', () => {
        expect(classifyConvergence(diverging)).toBe('diverging');
    });

    it('returns unknown for insufficient data', () => {
        expect(classifyConvergence([])).toBe('unknown');
        expect(classifyConvergence([1.0])).toBe('unknown');
        expect(classifyConvergence([1.0, 0.9])).toBe('unknown');
        expect(classifyConvergence(null)).toBe('unknown');
    });
});

describe('detectStall / detectDivergence', () => {
    it('match the classifier', () => {
        expect(detectStall(stalled)).toBe(true);
        expect(detectStall(converged)).toBe(false);
        expect(detectStall(improving)).toBe(false);
        expect(detectDivergence(diverging)).toBe(true);
        expect(detectDivergence(improving)).toBe(false);
    });
});

describe('every term and their total', () => {
    it('labels repeated header names apart', () => {
        expect(termLabels(['F(Q)_1', 'Curvature', 'Curvature', 'X_ray_(R)1'])).toEqual(['F(Q)_1', 'Curvature #1', 'Curvature #2', 'X_ray_(R)1']);
    });

    it('lists the terms that are zero throughout apart and judges the rest', () => {
        const zero = Array(100).fill(Math.log(1e-12));
        const plotData = {
            chiColumn: 'X_ray_(R)1',
            series: [{ y: improving.slice(0, 100) }],
            chiTerms: { names: ['Curvature', 'X_ray_(R)1'], series: [zero, improving.slice(0, 100)], total: improving.slice(0, 100), segments: [100] }
        };
        const judgement = judgeConvergence(chiSeries(plotData));
        expect(judgement.zeroTerms).toEqual(['Curvature']);
        expect(judgement.terms.map((term) => term.name)).toEqual(['X_ray_(R)1']);
        expect(judgement.status).toBe('improving');
    });

    it('a total that rises is diverging, even while a term falls; a term rising while nothing falls too', () => {
        const flat = Array(100).fill(0.5);
        const falling = improving.slice(0, 100);
        const rising = diverging.slice(100);
        const judge = (terms, total) => judgeConvergence({ terms: terms.map((values, k) => ({ name: `t${k}`, values })), total, start: 0 }).status;
        expect(judge([falling, rising], rising)).toBe('diverging');
        expect(judge([falling, rising], flat)).toBe('improving');
        expect(judge([flat, rising], flat)).toBe('diverging');
    });

    it('falls back on the charted column without terms, and on it beside the total when the logs name different terms', () => {
        expect(chiSeries({ chiColumn: 'X', series: [{ y: improving }] })).toMatchObject({ terms: [{ name: 'X' }], total: null, start: 0 });
        const mixed = chiSeries({ series: [{ y: improving }], chiTerms: { names: null, series: null, total: converged.slice(0, 200), segments: [150, 50] } });
        expect(mixed.terms.map((term) => term.name)).toEqual(['last .log column']);
        expect(mixed.start).toBe(150);
        expect(chiSeries(null)).toBeNull();
    });
});

describe('watchdogStats', () => {
    it('produces compact rounded stats', () => {
        const stats = watchdogStats(improving);
        expect(stats.n_steps).toBe(200);
        expect(stats.first).toBe(2.3);
        expect(stats.last).toBe(Number((2.3 - 199 * 0.01).toPrecision(3)));
        expect(stats.change_last_10pct).toBeCloseTo(-0.19, 6);
    });

    it('returns null with no history', () => {
        expect(watchdogStats([])).toBeNull();
    });

    it('summarises the total and names what still falls', () => {
        const plotData = {
            series: [{ y: converged }],
            chiTerms: { names: ['F(Q)_1', 'X_ray_(R)1'], series: [improving.concat(improving).slice(0, 300), converged], total: converged, segments: [300] }
        };
        expect(watchdogRunStats(plotData)).toMatchObject({ of: 'total', n_steps: 300, still_falling: ['F(Q)_1'] });
        expect(watchdogRunStats({ series: [{ y: converged }] })).not.toHaveProperty('of');
    });
});

describe('significantChange', () => {
    const base = { n_steps: 1000, last: 0.9 };

    it('is true with no previous stats and false with no next stats', () => {
        expect(significantChange(null, base)).toBe(true);
        expect(significantChange(base, null)).toBe(false);
    });

    it('triggers on enough new steps', () => {
        expect(significantChange(base, { n_steps: 1200, last: 0.9 })).toBe(true);
        expect(significantChange(base, { n_steps: 1050, last: 0.9 })).toBe(false);
    });

    it('triggers on a >2% chi^2 move (ln shift beyond ln(1.02))', () => {
        expect(significantChange(base, { n_steps: 1001, last: 0.9 - 0.03 })).toBe(true);
        expect(significantChange(base, { n_steps: 1001, last: 0.9 - 0.01 })).toBe(false);
    });
});

// A log row whose chi^2 is non-finite (NaN / Inf / overflow; null from Flask)
// means the run blew up; the parsers now keep those rows instead of dropping them.
describe('non-finite values', () => {
    const blownUp = [...improving.slice(0, 120), NaN, NaN, NaN];

    it('a non-finite latest value is divergence, never improvement', () => {
        expect(classifyConvergence(blownUp)).toBe('diverging');
        expect(detectDivergence(blownUp)).toBe(true);
        expect(detectStall(blownUp)).toBe(false);
        expect(classifyConvergence([...improving, null])).toBe('diverging');
    });

    it('earlier non-finite points do not break the trend', () => {
        const glitch = [...improving];
        glitch[50] = NaN;
        glitch[195] = NaN;
        expect(classifyConvergence(glitch)).toBe('improving');
    });

    it('watchdog stats stay numeric where defined and count the gaps', () => {
        const stats = watchdogStats(blownUp);
        expect(stats.last).toBeNull();
        expect(stats.non_finite_steps).toBe(3);
        expect(Number.isFinite(stats.min)).toBe(true);
    });
});
