// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Classic stog.inp first-peak line (line 22: "cutoff rmin rmax") semantics in
// the browser engine — must agree with scaling_cli.stog_inp_closest_approach.

import { describe, expect, it } from 'vitest';
import { stogInpClosestApproach } from '../workers/autoScale';

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
});
