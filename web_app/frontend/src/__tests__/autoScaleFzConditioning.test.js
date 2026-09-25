// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Conditioning of the Q->0 Faber-Ziman amplitude (1.0 audit, stog-b; mirrors
// tests/test_stog_b_fz_conditioning.py): fzLimitFit reports the standard error
// of S_meas(0) - level and flags an ill-conditioned aFz; parity with the Python
// goldens (fixture expected.fzLimit).

import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import {
  FZ_REL_SE_MAX, autoscale, diagnosticsSummary, fzLimitFit, levelSweep, makeConfig,
} from '../workers/autoScale';

const relError = (value, reference) => Math.abs(value - reference) / Math.abs(reference);
const Q = Float64Array.from(fixture.q);
const golden = fixture.expected.fzLimit;
const config = (extra = {}) => makeConfig({ ...fixture.config, bSqAvg: fixture.fzBSqAvg, ...extra });

const compare = (fit, expected) => {
  expect(fit.reliable).toBe(expected.reliable);
  expect(relError(fit.aFz, expected.a_fz)).toBeLessThan(1e-10);
  expect(relError(fit.sMeas0, expected.s_meas_0)).toBeLessThan(1e-10);
  expect(relError(fit.sMeas0Se, expected.s_meas_0_se)).toBeLessThan(1e-10);
  expect(relError(fit.aFzRelSe, expected.a_fz_rel_se)).toBeLessThan(1e-10);
};

describe('fzLimitFit (scaling.fz_limit_fit port)', () => {
  it('matches Python on the clean model (reliable)', () => {
    const sq = Float64Array.from(fixture.sqMeas);
    const sweep = levelSweep(Q, sq);
    const fit = fzLimitFit(Q, sq, sweep.level, config(), { levelUncertainty: sweep.levelUncertainty });
    compare(fit, golden.good);
    expect(fit.aFzRelSe).toBeLessThan(FZ_REL_SE_MAX);
  });

  it('matches Python on a head within noise of the level (flagged)', () => {
    const sq = Float64Array.from(golden.badSqMeas);
    const fit = fzLimitFit(Q, sq, golden.badLevel, config(), {
      levelUncertainty: golden.badLevelUncertainty,
    });
    compare(fit, golden.bad);
    expect(fit.reliable).toBe(false);
  });

  it('the summary reports it in fz and density mode', () => {
    const sq = Float64Array.from(golden.badSqMeas);
    const fzConfig = config({ amplitudeCriterion: 'fz' });
    const fz = diagnosticsSummary(autoscale(Q, sq, fzConfig), fzConfig);
    expect(fz.a_fz_reliable).toBe(false);
    expect(fz.a_fz_rel_se).toBeGreaterThan(0.5);
    const clean = diagnosticsSummary(autoscale(Q, Float64Array.from(fixture.sqMeas), config()), config());
    expect(clean.a_fz_reliable).toBe(true);
  });
});
