// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Iteration count of an unconverged auto-scale (1.0 audit, stog-b): Python's
// `for iterations in range(1, max_iter + 1)` leaves iterations == max_iter, and
// the JS loop must report the same (it reported maxIter + 1, one more than its
// own history) — the count is shown on the page and written to the provenance.

import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import { autoscale, diagnosticsSummary, makeConfig } from '../workers/autoScale';

const Q = Float64Array.from(fixture.q);
const SQ = Float64Array.from(fixture.sqMeas);

describe('auto-scale iteration count', () => {
  [1, 2].forEach((maxIter) => {
    it(`an unconverged run with maxIter = ${maxIter} reports ${maxIter} iterations`, () => {
      const config = makeConfig({ ...fixture.config, maxIter });
      const result = autoscale(Q, SQ, config);
      expect(result.converged).toBe(false);
      expect(result.iterations).toBe(maxIter);
      expect(result.history.length).toBe(maxIter);
      expect(diagnosticsSummary(result, config).iterations).toBe(maxIter);
    });
  });

  it('a converged run reports the iteration it converged on', () => {
    const result = autoscale(Q, SQ, makeConfig({ ...fixture.config }));
    expect(result.converged).toBe(true);
    expect(result.iterations).toBe(fixture.expected.auto.iterations);
    expect(result.history.length).toBe(result.iterations);
  });
});
