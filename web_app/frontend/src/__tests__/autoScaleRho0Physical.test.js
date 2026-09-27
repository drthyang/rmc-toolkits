// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// estimateRho0 never accepts a physically impossible root (0.6.0 audit, stog-b;
// mirrors tests/test_stog_b_rho0.py): the iterate stays in RHO0_PHYSICAL_RANGE
// and a concordant root must satisfy the density limit, else converged false
// with a reason the page reports.

import fs from 'node:fs';
import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import {
  RHO0_PHYSICAL_RANGE,
  estimateRho0,
  faberZiman,
  makeConfig,
  readStogXy,
  rho0NonConvergenceMessage,
} from '../workers/autoScale';

const Q = Float64Array.from(fixture.q);
const SQ = Float64Array.from(fixture.sqMeas);
const seeded = () => makeConfig({ ...fixture.config, rho0: 0.02, bSqAvg: fixture.fzBSqAvg });

describe('estimateRho0 physical range', () => {
  it('uses the Python default range', () => {
    expect(RHO0_PHYSICAL_RANGE).toEqual([0.005, 0.25]);
  });

  it('still converges on a genuine root', () => {
    const estimate = estimateRho0(Q, SQ, seeded());
    expect(estimate.converged).toBe(true);
    expect(estimate.reason).toBeNull();
  }, 60000);

  it('stops with a reason when the step leaves the range', () => {
    const estimate = estimateRho0(Q, SQ, seeded(), null, { rhoMax: 0.03 });
    expect(estimate.converged).toBe(false);
    expect(estimate.reason).toMatch(/physical density range/);
    estimate.history.forEach((row) => expect(row[0]).toBeLessThanOrEqual(0.03));
    expect(rho0NonConvergenceMessage(estimate)).toContain(estimate.reason);
    expect(() => estimateRho0(Q, SQ, seeded(), null, { rhoMin: 0.3 })).toThrow(/rhoMin < rhoMax/);
  }, 60000);
});

// Real missing-low-Q data (local-only, ~15 s, so opt-in like the other real-data
// JS checks: RMC_TOOLKITS_FULL_SWEEP=1 npm test; the Python suite runs it by
// default): with r0 pinned at 2.6 Å the 300 K run used to "converge" at
// 0.428 Å⁻³ from the true density as seed.
const RUN_300K = new URL(
  '../../../../data/stog_tests/stog_300K/PG3_55526_SQ_rebin.sq', import.meta.url,
);
const runRealData = Boolean(globalThis.process?.env?.RMC_TOOLKITS_FULL_SWEEP)
  && fs.existsSync(RUN_300K);

describe.skipIf(!runRealData)('Mn3Sn 300 K (local sample data)', () => {
  it('does not adopt the spurious high-density root', () => {
    const [q, sq] = readStogXy(fs.readFileSync(RUN_300K, 'utf8'));
    const fz = faberZiman('Mn3Sn');
    const config = makeConfig({
      qmin: 0.82, qmax: 28, rho0: 0.063049, r0: 2.6, bAvgSq: fz.bAvgSqBarn, bSqAvg: fz.bSqAvgBarn,
    });
    const estimate = estimateRho0(q, sq, config);
    expect(estimate.converged).toBe(false);
    expect(estimate.reason).toMatch(/physical density range/);
    estimate.history.forEach((row) => expect(row[0]).toBeLessThanOrEqual(0.25));
  }, 60000);
});
