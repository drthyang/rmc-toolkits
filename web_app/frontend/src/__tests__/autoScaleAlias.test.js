// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The r range the discrete sine transform resolves, r < π/ΔQ (1.0 audit,
// stog-b; mirrors tests/test_stog_b_alias.py): the engines report
// r_alias_limit and flag an r_max beyond it.

import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import { aliasLimit, diagnosticsSummary, fqToGpdf, makeConfig, scalePipeline } from '../workers/autoScale';

describe('aliasing limit π/ΔQ', () => {
  it('the transform folds a 20 Å shell to 2π/ΔQ − 20 Å, inverted', () => {
    const q = Float64Array.from({ length: 296 }, (_, i) => 0.1 * (i + 5));
    const fq = Float64Array.from(q, (value) => Math.exp(-0.5 * (0.1 * value) ** 2) * Math.sin(20 * value));
    const r = Float64Array.from({ length: 5000 }, (_, i) => 0.01 * (i + 1));
    const g = fqToGpdf(q, fq, r);
    let low = 0;
    let high = 0;
    for (let i = 1; i < g.length; i += 1) {
      if (g[i] < g[low]) low = i;
      if (g[i] > g[high]) high = i;
    }
    expect(r[low]).toBeCloseTo((2 * Math.PI) / 0.1 - 20, 1);
    expect(-g[low]).toBeGreaterThan(0.99 * g[high]);
    expect(aliasLimit(q)).toBeCloseTo(Math.PI / 0.1, 9);
  });

  it('diagnosticsSummary flags r_max beyond the limit', () => {
    const q = Float64Array.from(fixture.q); // ΔQ = 0.03 → π/ΔQ = 104.7 Å
    const sq = Float64Array.from(fixture.sqMeas);
    [[25, false], [120, true]].forEach(([rmax, beyond]) => {
      const config = makeConfig({ ...fixture.config, rmax, nr: rmax * 40 });
      const summary = diagnosticsSummary(scalePipeline(q, sq, config, fixture.aTrue, fixture.bTrue), config);
      expect(summary.r_alias_limit).toBeCloseTo(Math.PI / 0.03, 6);
      expect(summary.rmax_beyond_alias_limit).toBe(beyond);
    });
  });
});
