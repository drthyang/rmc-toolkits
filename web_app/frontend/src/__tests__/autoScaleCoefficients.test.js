// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// <b>^2 and <b^2> from one consistent source (0.6.0 audit, stog-b; mirrors
// tests/test_stog_b_coefficients.py): makeConfig rejects S(0) = 1 - <b^2>/<b>^2
// > 0, and resolveCoefficients never pairs a formula's <b^2> with a <b>^2 from
// another source (the page's composition + x-ray <b>^2 = 1 case).

import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import { faberZiman, makeConfig, resolveCoefficients } from '../workers/autoScale';

describe('makeConfig coefficient validation', () => {
  it('rejects <b^2> < <b>^2 (Cauchy-Schwarz)', () => {
    expect(() => makeConfig({ ...fixture.config, bAvgSq: 1, bSqAvg: 0.447511 }))
      .toThrow(/Cauchy-Schwarz/);
    expect(() => makeConfig({ ...fixture.config, bSqAvg: 0 })).toThrow(/bSqAvg/);
    expect(() => makeConfig({ ...fixture.config, bSqAvg: NaN })).toThrow(/bSqAvg/);
  });

  it('accepts S(0) = 0 and polyatomic pairs', () => {
    expect(() => makeConfig({ ...fixture.config, bAvgSq: 1, bSqAvg: 1 })).not.toThrow();
    expect(() => makeConfig({ ...fixture.config, bAvgSq: 1, bSqAvg: 1.10426 })).not.toThrow();
  });
});

describe('resolveCoefficients (scaling_cli.resolve_coefficients port)', () => {
  it('takes the formula pair when nothing else is given', () => {
    const fz = faberZiman('Mn3Sn');
    const out = resolveCoefficients({ formula: 'Mn3Sn' });
    expect(out.bAvgSq).toBe(fz.bAvgSqBarn);
    expect(out.bSqAvg).toBeCloseTo(fz.bSqAvgBarn, 15);
    expect(out.dropped).toBeNull();
  });

  it('keeps the formula ratio on an agreeing <b>^2', () => {
    const fz = faberZiman('Mn3Sn');
    const out = resolveCoefficients({ bAvgSq: 0.015407, formula: 'Mn3Sn' });
    expect(out.bSqAvg / out.bAvgSq).toBeCloseTo(fz.bSqAvgBarn / fz.bAvgSqBarn, 12);
  });

  it('never pairs the formula <b^2> with another source', () => {
    const out = resolveCoefficients({ bAvgSq: 1, formula: 'FeCoSn' });
    expect(out.bAvgSq).toBe(1);
    expect(out.bSqAvg).toBeUndefined();
    expect(out.dropped.bSqAvg).toBeCloseTo(0.447511, 6);
    expect(resolveCoefficients({ bAvgSq: 1, bSqAvg: 1.10426, formula: 'FeCoSn' }).bSqAvg)
      .toBe(1.10426);
  });
});
