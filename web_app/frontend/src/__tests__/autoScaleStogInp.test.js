// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Classic stog.inp first-peak line (line 22: "cutoff rmin rmax") semantics in
// the browser engine — must agree with scaling_cli.stog_inp_closest_approach.

import { describe, expect, it } from 'vitest';
import {
  autoscale, makeConfig, resolveEnforcementDescriptor, stogInpClosestApproach,
} from '../workers/autoScale';

const inp = (peakCutoff, peakRmin, peakRmax) => ({ peakCutoff, peakRmin, peakRmax });

describe('stogInpClosestApproach (r0 from the stog.inp peak line)', () => {
  it('uses the first-peak window start when it lies inside the cleanup cutoff', () => {
    // Pre-fix the page took max(cutoff, rmin) = 2.3 and fitted across the peak.
    expect(stogInpClosestApproach(inp(2.3, 1.6, 2.2), 1.0)).toBe(1.6);
  });

  it('uses the cutoff when the window lies outside [0, cutoff] or is empty', () => {
    expect(stogInpClosestApproach(inp(2.48, 2.65, 3.1), 1.0)).toBe(2.48);
    expect(stogInpClosestApproach(inp(2.7, 2.3, 2.2), 1.0)).toBe(2.7);
  });

  it('defers to detection when the line leaves no fit window (FeCoSn "1.0 0 0")', () => {
    expect(stogInpClosestApproach(inp(1.0, 0, 0), 1.0)).toBeNull();
  });

  it('defers to detection when the line leaves a sliver (< MIN_AUTO_WINDOW)', () => {
    // '1.46 0 0' at r-cut 1.0 pinned [1.2, 1.21] (FeCoSn a 15 % low); '1.0 0 0'
    // at r-cut 0.5 pinned [0.7, 0.75] (a 43 % low) — both reported converged.
    expect(stogInpClosestApproach(inp(1.46, 0, 0), 1.0)).toBeNull();
    expect(stogInpClosestApproach(inp(1.5, 0, 0), 1.0)).toBeNull();
    expect(stogInpClosestApproach(inp(1.0, 0, 0), 0.5)).toBeNull();
    expect(stogInpClosestApproach(inp(1.0, 0, 0), 0.54)).toBeNull();
    expect(stogInpClosestApproach(inp(1.55, 0, 0), 1.0)).toBe(1.55);
    expect(stogInpClosestApproach(inp(1.0, 0, 0), 0.4)).toBe(1.0);
  });
});

describe('low-r window validation (scaling.ScalingConfig / autoscale parity)', () => {
  const base = { qmin: 0.5, qmax: 30, rho0: 0.05, bAvgSq: 0.02 };

  it('rCutoff must be finite and >= 0', () => {
    [-1, -1e-9, NaN, Infinity].forEach((rCutoff) => {
      expect(() => makeConfig({ ...base, rCutoff })).toThrow(/rCutoff/);
    });
    expect(() => makeConfig({ ...base, rCutoff: 0 })).not.toThrow();
    ['r0', 'rFitMin', 'rFitMax'].forEach((name) => {
      expect(() => makeConfig({ ...base, [name]: NaN })).toThrow(new RegExp(name));
    });
    expect(() => makeConfig({ ...base, rFitMin: -0.5 })).toThrow(/rFitMin/);
  });

  it('autoscale refuses a pinned density-limit window narrower than 0.1 A', () => {
    const q = Float64Array.from({ length: 100 }, (_, i) => 0.5 + 0.2 * i);
    const sq = new Float64Array(100).fill(1);
    [{ r0: 1.2 }, { rFitMax: 0.95 }, { rFitMin: 1.2, rFitMax: 1.25 }].forEach((pins) => {
      const config = makeConfig({ ...base, rCutoff: 0.7, ...pins });
      expect(() => autoscale(q, sq, config)).toThrow(/narrower than 0.1/);
      // The page's one-line summary rides beside the unchanged message.
      let summary = null;
      try {
        autoscale(q, sq, config);
      } catch (error) {
        summary = error.summary;
      }
      expect(summary).toBe('Low-r fit window too narrow — widen it under Advanced → Low-r region.');
    });
  });
});

describe('resolveEnforcementDescriptor (page enforcement, CLI precedence)', () => {
  const inp59438 = { peakCutoff: 2.7, peakRmin: 2.3, peakRmax: 3.1 };

  it('keeps the stog.inp first-peak window for the pre-filled cutoff', () => {
    // selectSource pre-fills Cutoff with String(inp.peakCutoff); pre-fix the
    // page then flattened [2.3, 2.7] to -<b>^2 where the CLI keeps it.
    const prefilled = Number(String(inp59438.peakCutoff));
    expect(resolveEnforcementDescriptor({ enforce: true, cutoff: prefilled }, inp59438))
      .toEqual({ cutoff: 2.7, peakRmin: 2.3, peakRmax: 3.1 });
    expect(resolveEnforcementDescriptor({ enforce: true, cutoff: undefined }, inp59438))
      .toEqual({ cutoff: 2.7, peakRmin: 2.3, peakRmax: 3.1 });
  });

  it('a different typed cutoff is a flat replacement (like --enforce-cutoff)', () => {
    expect(resolveEnforcementDescriptor({ enforce: true, cutoff: 2.5 }, inp59438))
      .toEqual({ cutoff: 2.5, peakRmin: 2.5, peakRmax: 2.5 });
  });

  it('auto without a cutoff source, off when unchecked', () => {
    expect(resolveEnforcementDescriptor({ enforce: true, cutoff: undefined }, null)).toBe('auto');
    expect(resolveEnforcementDescriptor({ enforce: false, cutoff: 2.5 }, inp59438)).toBeNull();
  });
});
