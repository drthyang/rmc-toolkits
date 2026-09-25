// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The browser applies the CLI/API σ-column guard (1.0 audit, stog-b): one zero
// σ on a usable row got a 1e12 weight (a negative scale on real Mn3Sn data) and
// one NaN σ made every output NaN, where rmc-autoscale drops the column and
// fits unweighted. usableSigma is that guard (scaling_cli.usable_sigma).

import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import { autoscale, makeConfig, usableSigma } from '../workers/autoScale';

const Q = Float64Array.from(fixture.q);
const SQ = Float64Array.from(fixture.sqMeas);
const clean = () => Float64Array.from(Q, (q) => 1e-3 * (1 + 0.05 * q));

describe('usableSigma (CLI/API guard)', () => {
  it('keeps a clean column', () => {
    const sigma = clean();
    expect(usableSigma(Q, SQ, sigma)).toEqual({ sigma, nBad: 0 });
    expect(usableSigma(Q, SQ, null)).toEqual({ sigma: null, nBad: 0 });
  });

  it('drops the whole column for a zero, negative or NaN sigma on a usable row', () => {
    [0, -1e-3, NaN, Infinity].forEach((bad) => {
      const sigma = clean();
      sigma[900] = bad; // in the high-Q tail, where sigma weights the C1 rows
      expect(usableSigma(Q, SQ, sigma)).toEqual({ sigma: null, nBad: 1 });
    });
  });

  it('ignores sigma on rows that are not usable anyway', () => {
    const sq = Float64Array.from(SQ);
    sq[5] = NaN;
    const sigma = clean();
    sigma[5] = NaN;
    expect(usableSigma(Q, sq, sigma).sigma).toBe(sigma);
  });

  it('a guarded bad column fits exactly like no sigma (CLI parity)', () => {
    const config = makeConfig({ ...fixture.config });
    const sigma = clean();
    sigma[900] = 0;
    const guarded = autoscale(Q, SQ, config, usableSigma(Q, SQ, sigma).sigma);
    const unweighted = autoscale(Q, SQ, config, null);
    expect(guarded.a).toBe(unweighted.a);
    expect(guarded.b).toBe(unweighted.b);
    const nanSigma = clean();
    nanSigma[900] = NaN;
    // What the page used to send (a NaN fit never converges: cap the loop).
    const raw = autoscale(Q, SQ, makeConfig({ ...fixture.config, maxIter: 3 }), nanSigma);
    expect(raw.converged).toBe(false);
    expect(raw.fitFailure).toMatch(/non-physical scale a = NaN/);
  }, 60000);
});
