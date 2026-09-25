// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// estimateRho0's `extrapolated` flag follows the first measured Q, not
// config.qmin (1.0 audit, stog-b; mirrors tests/test_stog_b_extrapolated.py).

import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import { FZ_FIT_WIDTH, estimateRho0, makeConfig } from '../workers/autoScale';

const Q = Float64Array.from(fixture.q);

describe('estimateRho0 extrapolated flag', () => {
  it('is set by NaN-padded data even when qmin lies below it', () => {
    const padded = Float64Array.from(fixture.sqMeas, (value, i) => (Q[i] < 1.3 ? NaN : value));
    const config = (qmin) => makeConfig({
      ...fixture.config, qmin, rho0: 0.02, bSqAvg: fixture.fzBSqAvg,
    });
    const low = estimateRho0(Q, padded, config(0.6));
    const atData = estimateRho0(Q, padded, config(1.3));
    expect(low.qFirst).toBeCloseTo(1.32, 9);
    // The same estimate (qmin only moves the C1 tail window edge) ...
    expect(Math.abs(low.rho0 - atData.rho0) / atData.rho0).toBeLessThan(1e-6);
    expect(low.extrapolated).toBe(true);
    expect(atData.extrapolated).toBe(true);
  }, 60000); // two unconverged 8-pass estimates: slow under a parallel run

  it('is clear when the data start below the fit width', () => {
    const estimate = estimateRho0(Q, Float64Array.from(fixture.sqMeas), makeConfig({
      ...fixture.config, qmin: 0.3, rho0: 0.02, bSqAvg: fixture.fzBSqAvg,
    }));
    expect(estimate.qFirst).toBeLessThan(FZ_FIT_WIDTH);
    expect(estimate.extrapolated).toBe(false);
  }, 60000);
});
