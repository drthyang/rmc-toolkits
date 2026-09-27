// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The main (density-limit) auto path with a composition (0.6.0 audit, stog-b):
// with <b^2> known the omitted-low-Q extrapolation targets
// S(0) = 1 - <b^2>/<b>^2, and the self-consistent loop's Fourier filter must use
// that target exactly like scaling._pipeline. The JS loop used S(0) = 0, so the
// page's (a, b) differed from the CLI/API by ~2 % on the Mn3Sn runs. Golden
// numbers from tests/generate_autoscale_fixture.py; the engines are the same
// deterministic algorithm, so agreement is at round-off.

import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import { autoscale, effectiveS0Target, makeConfig } from '../workers/autoScale';

const relError = (value, reference) => Math.abs(value - reference) / Math.abs(reference);

const Q = Float64Array.from(fixture.q);
const SQ = Float64Array.from(fixture.sqMeas);

describe('auto-scale with a composition (S(0) target inside the loop)', () => {
  fixture.expected.autoComposition.forEach((expected) => {
    it(`${expected.name}: (a, b) match Python to round-off`, () => {
      const config = makeConfig({ ...fixture.config, ...expected.config });
      expect(effectiveS0Target(config)).toBeCloseTo(expected.s0Target, 14);
      const result = autoscale(Q, SQ, config);
      expect(result.converged).toBe(expected.converged);
      expect(result.iterations).toBe(expected.iterations);
      expect(relError(result.a, expected.a)).toBeLessThan(1e-10);
      expect(relError(result.b, expected.b)).toBeLessThan(1e-10);
      expect(relError(result.lowRRms, expected.lowRRms)).toBeLessThan(1e-9);
      expect(relError(result.c1TailMean, expected.c1TailMean)).toBeLessThan(1e-12);
      expect(relError(result.aFz, expected.aFz)).toBeLessThan(1e-10);
      if (expected.config.r0 === null) {
        expect(result.r0Detected).toBeCloseTo(expected.r0Detected, 12);
      }
    });
  });
});
