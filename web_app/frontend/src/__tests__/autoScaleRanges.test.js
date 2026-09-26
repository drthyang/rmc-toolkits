// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Q-range and low-r enforcement checks, the same rules as ScalingConfig and
// validate_enforcement in scaling.py (tests/test_scaling_range_validation.py).
// A NaN or negative qmin, a non-finite qmax, a cutoff at or beyond rmax (the
// whole G(r) replaced), a NaN cutoff (reported, never applied) and a reversed
// first-peak window used to pass.

import { afterEach, describe, expect, it, vi } from 'vitest';
import { makeConfig, validateEnforcement } from '../workers/autoScale';

const base = { qmin: 0.6, qmax: 29, rho0: 0.05, bAvgSq: 0.02, rmax: 25, nr: 1000 };

describe('makeConfig Q range', () => {
    it('requires a finite qmin >= 0 and a finite qmax', () => {
        expect(() => makeConfig({ ...base, qmin: NaN })).toThrow('qmin must be finite and >= 0');
        expect(() => makeConfig({ ...base, qmin: -3 })).toThrow('qmin must be finite and >= 0');
        expect(() => makeConfig({ ...base, qmax: NaN })).toThrow('qmax must be finite');
        expect(() => makeConfig({ ...base, qmax: Infinity })).toThrow('qmax must be finite');
        expect(makeConfig({ ...base, qmin: 0 }).qmin).toBe(0);
    });
});

describe('validateEnforcement', () => {
    it('accepts a real triple and refuses the rest', () => {
        expect(() => validateEnforcement({ cutoff: 2.48, peakRmin: 2.65, peakRmax: 3.1 }, 25)).not.toThrow();
        expect(() => validateEnforcement({ cutoff: 0, peakRmin: 0, peakRmax: 0 }, 25)).not.toThrow();
        for (const [triple, message] of [
            [{ cutoff: NaN, peakRmin: 2, peakRmax: 2 }, 'enforcement cutoff must be finite and >= 0'],
            [{ cutoff: -1, peakRmin: -1, peakRmax: -1 }, 'enforcement cutoff must be finite and >= 0'],
            [{ cutoff: 25, peakRmin: 25, peakRmax: 25 }, 'must be below rmax'],
            [{ cutoff: 1000, peakRmin: 1000, peakRmax: 1000 }, 'must be below rmax'],
            [{ cutoff: 2, peakRmin: 2.5, peakRmax: 1.5 }, 'first-peak window must be finite with rmin <= rmax'],
            [{ cutoff: 2, peakRmin: NaN, peakRmax: 3 }, 'first-peak window must be finite with rmin <= rmax'],
        ]) {
            expect(() => validateEnforcement(triple, 25)).toThrow(message);
        }
    });
});

describe('the Auto StoG worker checks an explicit enforcement before computing', () => {
    afterEach(() => vi.unstubAllGlobals());

    it('answers an error for a cutoff beyond rmax', async () => {
        vi.stubGlobal('self', { postMessage: () => {} });
        vi.resetModules();
        const { runAutoScaleJob } = await import('../workers/autoScaleWorker.js');
        const q = Array.from({ length: 1000 }, (_, i) => 0.6 + i * 0.0284);
        const sq = q.map((x) => 1 + 0.3 * Math.sin(2.7 * x) * Math.exp(-0.05 * x * x));
        expect(() => runAutoScaleJob({
            id: 1, config: makeConfig(base), q, sq, mode: 'manual', a: 1, b: 0,
            enforcement: { cutoff: 1000, peakRmin: 1000, peakRmax: 1000 },
        })).toThrow('must be below rmax');
    });
});
