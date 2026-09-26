// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Lorch omitted-low-Q basis without catastrophic cancellation (1.0 audit,
// stog-b; mirrors tests/test_stog_b_lorch_basis.py): the old
// (cos v - 1)/(r - a)^2 form lost every digit within ~1e-6 Å of r = π/Qmax
// outside its 1e-9 patch. The band |r - a| = 1e-9 .. 1e-5 (and r = a, r = 0,
// the default 0.01 Å grid) is checked against composite Simpson quadrature of
// the defining integrals.

import { describe, expect, it } from 'vitest';
import { lowQCorrectionBasis } from '../workers/autoScale';

const simpson = (f, lo, hi, n = 4000) => {
  const h = (hi - lo) / n;
  let total = f(lo) + f(hi);
  for (let k = 1; k < n; k += 1) total += (k % 2 ? 4 : 2) * f(lo + k * h);
  return (total * h) / 3;
};

const reference = (q0, qmax, r) => {
  const a = Math.PI / qmax;
  const window = (q) => (q === 0 ? 1 : Math.sin(a * q) / (a * q));
  return {
    coef: (2 / Math.PI) * simpson((q) => q * (q / q0) * window(q) * Math.sin(q * r), 0, q0),
    constant: (2 / Math.PI) * simpson((q) => q * window(q) * Math.sin(q * r), 0, q0),
  };
};

const grid = (lo, hi, n) => Float64Array.from({ length: n }, (_, i) => lo + ((hi - lo) * i) / (n - 1));

describe('lowQCorrectionBasis (Lorch) near r = π/Qmax', () => {
  [[1.0, 28.0], [0.5, 28.56], [1.0, 52.36], [0.5, 26.18]].forEach(([q0, qmax]) => {
    it(`q0 ${q0}, Qmax ${qmax}: coefficients match quadrature across the band`, () => {
      const a = Math.PI / qmax;
      const r = [0, a];
      for (let e = -9; e <= -5; e += 0.5) r.push(a - 10 ** e, a + 10 ** e);
      for (let k = 1; k <= 30; k += 1) r.push(0.01 * k);
      const { coef, constant } = lowQCorrectionBasis(
        grid(q0, qmax, 400), Float64Array.from(r), { lorch: true },
      );
      r.forEach((ri, i) => {
        const ref = reference(q0, qmax, ri);
        if (ri === 0) {
          expect(coef[i]).toBe(0);
          expect(constant[i]).toBe(0);
          return;
        }
        expect(Math.abs(coef[i] - ref.coef) / Math.abs(ref.coef)).toBeLessThan(1e-9);
        expect(Math.abs(constant[i] - ref.constant) / Math.abs(ref.constant)).toBeLessThan(1e-9);
      });
    });
  });
});
