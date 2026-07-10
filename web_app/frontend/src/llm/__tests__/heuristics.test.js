// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import {
    classifyConvergence,
    detectDivergence,
    detectStall,
    significantChange,
    watchdogStats
} from '../watchdog/heuristics';

// Synthetic ln(chi^2) histories for each regime.
const improving = Array.from({ length: 200 }, (_, index) => 2.3 - index * 0.01);
const converged = [
    ...Array.from({ length: 150 }, (_, index) => 2.3 - index * 0.01),
    ...Array.from({ length: 150 }, () => 0.8)
];
const stalled = Array.from({ length: 200 }, (_, index) => 2.3 - 0.0001 * index);
const diverging = [
    ...Array.from({ length: 100 }, (_, index) => 2.3 - index * 0.01),
    ...Array.from({ length: 100 }, (_, index) => 1.3 + index * 0.01)
];

describe('classifyConvergence', () => {
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

describe('watchdogStats', () => {
    it('produces compact rounded stats', () => {
        const stats = watchdogStats(improving);
        expect(stats.n_steps).toBe(200);
        expect(stats.first).toBe(2.3);
        expect(stats.last).toBe(Number((2.3 - 199 * 0.01).toPrecision(3)));
        expect(stats.recent_window_delta).toBeLessThan(0);
    });

    it('returns null with no history', () => {
        expect(watchdogStats([])).toBeNull();
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
