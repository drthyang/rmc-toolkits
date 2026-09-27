// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Q->0 head fit is Huber's M-estimator (twin of tests/test_huber_irls.py):
// IRLS rows scaled by sqrt(w), so each pass minimises sum w r^2. The reference
// solves the weighted normal equations X^T W X beta = X^T W y with W = diag(w)
// for the engine's pass count from the unweighted start; before 0.6.0 the rows
// were scaled by w (effective weight w^2), which this head tells apart.

import { describe, expect, it } from 'vitest';
import { fzLimitFit, makeConfig } from '../workers/autoScale';

const HUBER_C = 1.345;

const median = (values) => {
  const sorted = [...values].sort((x, y) => x - y);
  const mid = sorted.length >> 1;
  return sorted.length % 2 ? sorted[mid] : 0.5 * (sorted[mid - 1] + sorted[mid]);
};

const huberWeights = (residuals) => {
  const centre = median(residuals);
  const scale = 1.4826 * median(residuals.map((value) => Math.abs(value - centre)));
  if (scale <= 1e-14) return residuals.map(() => 1);
  return residuals.map((value) => Math.min(1, (HUBER_C * scale) / Math.max(Math.abs(value), 1e-14 * scale)));
};

// Weighted straight line y = c0 + c1 x by the 2 x 2 normal equations.
const weightedLine = (x, y, w) => {
  let s = 0; let sx = 0; let sxx = 0; let sy = 0; let sxy = 0;
  for (let i = 0; i < x.length; i += 1) {
    s += w[i]; sx += w[i] * x[i]; sxx += w[i] * x[i] * x[i];
    sy += w[i] * y[i]; sxy += w[i] * x[i] * y[i];
  }
  const det = s * sxx - sx * sx;
  return [(sxx * sy - sx * sxy) / det, (s * sxy - sx * sy) / det];
};

const referenceIntercept = (q, y, passes, rowPower = 1) => {
  const qMean = q.reduce((total, value) => total + value, 0) / q.length;
  const x = q.map((value) => value - qMean);
  let beta = weightedLine(x, y, x.map(() => 1));
  for (let pass = 0; pass < passes; pass += 1) {
    const weights = huberWeights(x.map((value, i) => beta[0] + beta[1] * value - y[i]));
    beta = weightedLine(x, y, weights.map((w) => w ** rowPower));
  }
  return beta[0] - beta[1] * qMean;
};

// Deterministic contaminated head: noise plus one-sided spikes on every eighth row.
const contaminatedHead = () => {
  let state = 5;
  const random = () => {
    state = (state * 16807) % 2147483647;
    return state / 2147483647;
  };
  const q = [];
  const y = [];
  for (let i = 0; i < 100; i += 1) {
    const value = 0.8 + 0.01 * i;
    const noise = 0.01 * Math.sqrt(-2 * Math.log(random())) * Math.cos(2 * Math.PI * random());
    const spike = i % 8 === 3 ? 0.1 + 0.1 * random() : 0;
    q.push(value);
    y.push(0.4 + 0.3 * value + noise + spike);
  }
  return { q, y };
};

describe('fzLimitFit Huber IRLS', () => {
  it('matches the reference IRLS pass for pass, not the old w^2 weighting', () => {
    const { q, y } = contaminatedHead();
    const config = makeConfig({ qmin: 0.8, qmax: 30, rho0: 0.05, bAvgSq: 1, bSqAvg: 2 });
    const fit = fzLimitFit(Float64Array.from([...q, 30]), Float64Array.from([...y, 1]), 1, config);
    const reference = referenceIntercept(q, y, 3);
    expect(Math.abs(fit.sMeas0 - reference)).toBeLessThan(1e-12);
    expect(Math.abs(referenceIntercept(q, y, 3, 2) - reference)).toBeGreaterThan(1e-3);
  });
});
