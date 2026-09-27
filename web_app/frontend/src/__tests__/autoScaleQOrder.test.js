// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Q-grid order (0.6.0 audit, stog-b; mirrors tests/test_stog_b_qorder.py). With
// descending Q the trapezoid panels are negative and autoscale "converged" to a
// negative scale in both engines. Now cropSq sorts ascending and rejects
// duplicate / overlapping Q, the transforms throw on non-increasing grids, and
// a <= 0 is a failed fit (converged false, fitFailure set).

import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import {
  autoscale,
  cropSq,
  diagnosticsSummary,
  fourierFilter,
  fqToGpdf,
  lowQCorrectionBasis,
  makeConfig,
  sineTransform,
} from '../workers/autoScale';

const Q = Float64Array.from(fixture.q);
const SQ = Float64Array.from(fixture.sqMeas);
const reversed = (values) => Float64Array.from(values).reverse();
const config = () => makeConfig({ ...fixture.config });

describe('cropSq orders Q ascending', () => {
  it('sorts a descending file together with its sigma', () => {
    const sigma = Float64Array.from(Q, (q) => 1e-3 * (1 + q));
    const asc = cropSq(Q, SQ, config(), sigma);
    const desc = cropSq(reversed(Q), reversed(SQ), config(), reversed(sigma));
    expect(Array.from(desc.q)).toEqual(Array.from(asc.q));
    expect(Array.from(desc.sq)).toEqual(Array.from(asc.sq));
    expect(Array.from(desc.sigma)).toEqual(Array.from(asc.sigma));
  });

  it('sorts non-overlapping segments written high-Q first', () => {
    const half = Math.floor(Q.length / 2);
    const order = [...Array.from({ length: Q.length - half }, (_, i) => half + i),
      ...Array.from({ length: half }, (_, i) => i)];
    const cropped = cropSq(order.map((i) => Q[i]), order.map((i) => SQ[i]), config());
    expect(Array.from(cropped.q)).toEqual(Array.from(Q));
  });

  it('rejects overlapping banks and duplicate Q', () => {
    const bank1 = Array.from(Q).filter((q) => q <= 16);
    const bank2 = Array.from(Q).filter((q) => q >= 14).map((q) => q + 0.005);
    const sq1 = Array.from(SQ).slice(0, bank1.length);
    const sq2 = Array.from(SQ).slice(Q.length - bank2.length);
    expect(() => cropSq([...bank1, ...bank2], [...sq1, ...sq2], config())).toThrow(/overlap/);
    const qDup = [...Array.from(Q).slice(0, 100), ...Array.from(Q).slice(99)];
    const sqDup = [...Array.from(SQ).slice(0, 100), ...Array.from(SQ).slice(99)];
    expect(() => cropSq(qDup, sqDup, config())).toThrow(/duplicate/);
  });
});

describe('auto-scale on descending data', () => {
  it('gives exactly the ascending result', () => {
    const asc = autoscale(Q, SQ, config());
    const desc = autoscale(reversed(Q), reversed(SQ), config());
    expect(asc.a).toBeGreaterThan(0);
    expect(desc.a).toBe(asc.a);
    expect(desc.b).toBe(asc.b);
    expect(desc.converged).toBe(true);
  });

  it('reports a non-positive scale as a failed fit', () => {
    // A sign-inverted measurement: the fit needs a < 0.
    const inverted = Float64Array.from(SQ, (value) => 2 - (10 * value - 9));
    const cfg = config();
    const result = autoscale(Q, inverted, cfg);
    expect(result.a).toBeLessThanOrEqual(0);
    expect(result.converged).toBe(false);
    expect(result.fitFailure).toMatch(/a <= 0/);
    expect(diagnosticsSummary(result, cfg).fit_failure).toMatch(/non-physical scale/);
  });
});

describe('transforms require an increasing grid', () => {
  const r = Float64Array.from({ length: 200 }, (_, i) => 0.05 * (i + 1));
  const fq = Float64Array.from(Q, (q, i) => q * (SQ[i] - 1));

  it('throws on descending Q or r', () => {
    expect(() => fqToGpdf(reversed(Q), reversed(fq), r)).toThrow(/strictly increasing/);
    expect(() => lowQCorrectionBasis(reversed(Q), r)).toThrow(/strictly increasing/);
    expect(() => fourierFilter(reversed(Q), reversed(SQ), r, { rho0: 0.05, cutoff: 1 }))
      .toThrow(/strictly increasing/);
    expect(() => fourierFilter(Q, SQ, reversed(r), { rho0: 0.05, cutoff: 1 }))
      .toThrow(/strictly increasing/);
  });

  it('throws on a repeated or NaN point, and integrates short sections to 0', () => {
    const q = Float64Array.from(Q);
    q[10] = q[9];
    expect(() => sineTransform(q, fq, r)).toThrow(/strictly increasing/);
    q[10] = NaN;
    expect(() => sineTransform(q, fq, r)).toThrow(/strictly increasing/);
    expect(Math.max(...sineTransform([], [], Q))).toBe(0);
    expect(Math.max(...sineTransform([0.5], [1], Q))).toBe(0);
  });
});
