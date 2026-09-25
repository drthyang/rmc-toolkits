// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// fourierFilter on an r grid that starts at r = 0 (1.0 audit, stog-b; mirrors
// tests/test_stog_b_rzero.py): the section integrand 4π ρ0 r g was 0/0 at r = 0
// and NaN reached every output. Now it is G_PDF + 4π ρ0 r, and g(0) is the
// continuous extension 1 + G_PDF'(0)/(4π ρ0) (gpdfSlopeAtZero).

import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import { fourierFilter, fqToGpdf, gpdfSlopeAtZero } from '../workers/autoScale';

const Q = Float64Array.from(fixture.q);
const SQ = Float64Array.from(fixture.sqMeas, (value) => 10 * value - 9); // the true S(Q)
const RHO0 = 0.05;
const OPTIONS = [];
[false, true].forEach((lorch) => [false, true].forEach((lowQCorrection) => [0, -12.06]
  .forEach((s0Target) => OPTIONS.push({ lorch, lowQCorrection, s0Target }))));

describe('fourierFilter with r = 0 in the grid', () => {
  const r0 = Float64Array.from({ length: 1001 }, (_, i) => 0.01 * i);

  OPTIONS.forEach((options) => {
    it(`finite outputs and a continuous g(0): ${JSON.stringify(options)}`, () => {
      const out = fourierFilter(Q, SQ, r0, { rho0: RHO0, cutoff: 1, ...options });
      [out.sqFiltered, out.sqFt, out.gFiltered].forEach((array) => {
        expect(array.every(Number.isFinite)).toBe(true);
      });
      const extrapolated = (4 * out.gFiltered[1] - out.gFiltered[2]) / 3;
      expect(Math.abs(out.gFiltered[0] - extrapolated)).toBeLessThan(1e-3);
      const ref = fourierFilter(Q, SQ, r0.subarray(1), { rho0: RHO0, cutoff: 1, ...options });
      let worst = 0;
      for (let i = 0; i < Q.length; i += 1) worst = Math.max(worst, Math.abs(out.sqFiltered[i] - ref.sqFiltered[i]));
      expect(worst).toBeLessThan(1e-5);
    });
  });

  it('gpdfSlopeAtZero is the derivative of fqToGpdf at 0 (Richardson)', () => {
    const fq = Float64Array.from(Q, (q, i) => q * (SQ[i] - 1));
    OPTIONS.forEach((options) => {
      const h = 0.002;
      const g = fqToGpdf(Q, fq, Float64Array.from([h, 2 * h, 4 * h]), options);
      const s1 = (8 * g[0] - g[1]) / (6 * h);
      const s2 = (8 * g[1] - g[2]) / (12 * h);
      const richardson = (16 * s1 - s2) / 15;
      const slope = gpdfSlopeAtZero(Q, fq, options);
      expect(Math.abs(slope - richardson) / Math.abs(richardson)).toBeLessThan(1e-9);
    });
  });

  it('rejects a negative r grid', () => {
    expect(() => fourierFilter(Q, SQ, Float64Array.from([-1, 0, 1, 2]), { rho0: RHO0, cutoff: 1 }))
      .toThrow(/non-negative/);
  });
});
