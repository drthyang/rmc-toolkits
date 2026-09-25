// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Unwindowed omitted-low-Q basis at small v = Q0 r (stog-b; mirrors
// tests/test_stog_b_lowq_series.py): the closed forms cancel O(1) terms down to
// O(v^3), so below |v| = 0.5 the Taylor series is used. Checked against
// composite Simpson quadrature of the defining moments.

import { describe, expect, it } from 'vitest';
import { lowQCorrectionBasis } from '../workers/autoScale';

const simpson = (f, lo, hi, n = 2000) => {
  const h = (hi - lo) / n;
  let total = f(lo) + f(hi);
  for (let k = 1; k < n; k += 1) total += (k % 2 ? 4 : 2) * f(lo + k * h);
  return (total * h) / 3;
};

describe('lowQCorrectionBasis (unwindowed) at small Q0 r', () => {
  [0.01, 0.5, 1.0, 2.0].forEach((q0) => {
    it(`q0 ${q0}: moments match quadrature`, () => {
      const q = Float64Array.from({ length: 300 }, (_, i) => q0 + ((28 - q0) * i) / 299);
      const r = [1e-6, 1e-4, 1e-3, 0.002, 0.01, 0.1, 0.49 / q0, 0.51 / q0, 1.0, 5.0];
      const { coef, constant } = lowQCorrectionBasis(q, Float64Array.from(r));
      r.forEach((ri, i) => {
        const refCoef = (2 / Math.PI) * simpson((x) => x * (x / q0) * Math.sin(x * ri), 0, q0);
        const refConst = (2 / Math.PI) * simpson((x) => x * Math.sin(x * ri), 0, q0);
        expect(Math.abs(coef[i] - refCoef) / Math.abs(refCoef)).toBeLessThan(1e-10);
        expect(Math.abs(constant[i] - refConst) / Math.abs(refConst)).toBeLessThan(1e-10);
      });
    });
  });
});
