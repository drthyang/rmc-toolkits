// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Low-r window placement loop of the static-mode Auto StoG engine
// (placeLowRWindow, the trial / confirm loop inside autoscale): replays the
// scripted-pass scenarios of tests/test_stog_a_placement.py and checks parity
// with scaling._place_low_r_window on the golden outcomes in the fixture.
// A refit with a <= 0 is never returned, an onset its own refit no longer
// shows is dropped, and r0Detected is the onset the window was built from.

import fs from 'node:fs';
import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import {
  MAX_WINDOW_REFITS,
  autoscale,
  diagnosticsSummary,
  faberZiman,
  firstShellCandidates,
  makeConfig,
  placeLowRWindow,
  readStogXy,
  rFitWindow,
} from '../workers/autoScale';

// Numbers in an error message, not the digit of "r0" (generate_autoscale_fixture.NUMBER).
const NUMBER = /(?<![A-Za-z_\d.])-?\d+(?:\.\d+)?/g;
const numbersIn = (message) => (message.match(NUMBER) || []).map(Number);

const { r, profiles, config: baseConfig, cases } = fixture.expected.placement;
const rArr = Float64Array.from(r);
const profileArrays = Object.fromEntries(
  Object.entries(profiles).map(([name, g]) => [name, Float64Array.from(g)]),
);

/** The scripted fit pass of tests/test_stog_a_placement.scripted_pass. */
const scriptedPass = (scenario, calls) => (config) => {
  const [lo, hi] = rFitWindow(config);
  if (hi - lo < 0.03) throw new Error('fit windows contain fewer than 2 points');
  let row;
  if (config.r0 == null) {
    row = scenario.trials.find((item) => Math.abs(hi - lo - item[0]) < 1e-6).slice(1);
  } else {
    const match = scenario.refits.find((item) => Math.abs(item[0] - config.r0) <= 0.05);
    row = match ? match.slice(1) : scenario.default;
    calls.push(config.r0);
  }
  const [a, name] = row;
  return { a, r: rArr, gFiltered: profileArrays[name], rFitWindowUsed: [lo, hi] };
};

describe('window placement parity with scaling._place_low_r_window', () => {
  cases.forEach((scenario) => {
    it(`${scenario.name}: ${scenario.expected ? 'confirmed below the first shell' : scenario.error.kind}`, () => {
      const config = makeConfig({ ...baseConfig, qmax: scenario.qmax });
      const calls = [];
      const run = () => placeLowRWindow(scriptedPass(scenario, calls), config);
      if (scenario.error) {
        let message = null;
        try {
          run();
        } catch (error) {
          message = error.message;
        }
        expect(message).not.toBeNull();
        expect(message).toContain(scenario.error.kind);
        expect(numbersIn(message)).toEqual(scenario.error.numbers);
      } else {
        const result = run();
        expect(result.a).toBe(scenario.expected.a);
        expect(result.a).toBeGreaterThan(0);
        expect(result.r0Detected).toBe(scenario.expected.r0Detected);
        expect(result.windowRefined).toBe(true);
        // r0Detected is the onset the returned window was built from.
        expect(result.rFitWindowUsed[1]).toBeCloseTo(result.r0Detected - 0.25, 12);
        expect(result.rFitWindowUsed[1]).toBeCloseTo(scenario.expected.rFitWindow[1], 12);
      }
      expect(calls).toEqual(scenario.refitOnsets);
      expect(calls.length).toBeLessThanOrEqual(MAX_WINDOW_REFITS);
    });
  });
});

describe('first-shell candidate lists match scaling.first_shell_candidates', () => {
  const { r: detR, searchMin, profiles: detProfiles, cases: detCases } = fixture.expected.detector;
  const detRArr = Float64Array.from(detR);
  detCases.forEach(({ name, qmax, candidates }) => {
    it(`${name} (qmax ${qmax})`, () => {
      expect(firstShellCandidates(detRArr, Float64Array.from(detProfiles[name]), { searchMin, qmax }))
        .toEqual(candidates);
    });
  });
});

// Real missing-low-Q data (local-only, ~25 s, so opt-in like the Python full
// sweep: RMC_TOOLKITS_FULL_SWEEP=1 npm test): the Mn3Sn 59438 run at Qmin 1.0 /
// Qmax 29 returned a = -0.584 on the window [1.2, 1.33] before the review fix.
const RUN_59438 = new URL(
  '../../../../data/stog_tests/stog_59438/PG3_59438_SQ_rebin.dat', import.meta.url,
);
const runRealData = Boolean(globalThis.process?.env?.RMC_TOOLKITS_FULL_SWEEP)
  && fs.existsSync(RUN_59438);

describe.skipIf(!runRealData)('Mn3Sn 59438 (local sample data)', () => {
  it('returns a positive scale with the window below the inverted first shell', () => {
    const [q, sq] = readStogXy(fs.readFileSync(RUN_59438, 'utf8'));
    const fz = faberZiman('Mn3Sn');
    const config = makeConfig({
      qmin: 1.0, qmax: 29.0, rho0: 0.063049, bAvgSq: fz.bAvgSqBarn, bSqAvg: fz.bSqAvgBarn,
    });
    const result = autoscale(q, sq, config);
    const summary = diagnosticsSummary(result, config);
    expect(result.a).toBeGreaterThan(0);
    expect(result.rFitWindowUsed[1]).toBeLessThan(2.65);
    expect(result.rFitWindowUsed[1]).toBeCloseTo(result.r0Detected - 0.25, 12);
    expect(result.r0Detected).toBeGreaterThan(2.4);
    expect(result.r0Detected).toBeLessThan(2.9);
    expect(summary.density_limit_satisfied).toBe(false);
  }, 120000);
});
