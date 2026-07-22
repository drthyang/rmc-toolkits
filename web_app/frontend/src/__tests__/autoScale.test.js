// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Parity tests for the static-mode Auto StoG engine against Python-generated
// golden numbers (tests/generate_autoscale_fixture.py). Tolerances are loose
// enough for summation-order float noise (numpy pairwise vs JS sequential)
// and tight enough that any real math drift fails.

import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import {
  autoscale,
  detectFirstPeakOnset,
  diagnosticsSummary,
  estimateRho0,
  faberZiman,
  firstPeakZero,
  levelSweep,
  makeConfig,
  massDensityFromNumberDensity,
  numberDensityFromMassDensity,
  parseFormula,
  readDatHeader,
  readStogInp,
  readStogXy,
  scalePipeline,
  writeStogXy,
} from '../workers/autoScale';

const relError = (value, reference) => Math.abs(value - reference) / Math.abs(reference);

const baseConfig = () => makeConfig({ ...fixture.config });

const Q = Float64Array.from(fixture.q);
const SQ = Float64Array.from(fixture.sqMeas);

describe('autoScale engine parity with the Python engine', () => {
  it('level sweep reproduces the Python window and level', () => {
    const sweep = levelSweep(Q, SQ);
    const expected = fixture.expected.sweep;
    expect(sweep.asymptoteFound).toBe(true);
    expect(sweep.nAdmissible).toBe(expected.nAdmissible);
    expect(relError(sweep.level, expected.level)).toBeLessThan(1e-9);
    expect(sweep.qLo).toBeCloseTo(expected.qLo, 9);
    expect(sweep.qHi).toBeCloseTo(expected.qHi, 9);
    expect(relError(sweep.levelUncertainty, expected.levelUncertainty)).toBeLessThan(1e-6);
  });

  it('auto-scale (sweep + density) matches Python and recovers the truth', () => {
    const result = autoscale(Q, SQ, baseConfig());
    const expected = fixture.expected.auto;
    expect(result.converged).toBe(true);
    expect(result.iterations).toBe(expected.iterations);
    expect(relError(result.a, expected.a)).toBeLessThan(1e-6);
    expect(relError(result.b, expected.b)).toBeLessThan(1e-6);
    expect(relError(result.lowRRms, expected.lowRRms)).toBeLessThan(1e-5);
    expect(relError(result.c1TailMean, expected.c1TailMean)).toBeLessThan(1e-8);
    // The physics check the whole feature stands on:
    expect(relError(result.a, fixture.aTrue)).toBeLessThan(0.02);
    expect(Math.abs(result.b - fixture.bTrue) / Math.abs(fixture.bTrue)).toBeLessThan(0.02);
  });

  it('FZ amplitude mode matches Python (closed form, no loop)', () => {
    const config = makeConfig({
      ...fixture.config,
      bSqAvg: fixture.fzBSqAvg,
      amplitudeCriterion: 'fz',
    });
    const result = autoscale(Q, SQ, config);
    expect(result.iterations).toBe(0);
    expect(relError(result.a, fixture.expected.fz.a)).toBeLessThan(1e-6);
    expect(relError(result.b, fixture.expected.fz.b)).toBeLessThan(1e-6);
    expect(result.b).toBeCloseTo(1 - result.a * result.sweep.level, 10);
  });

  it('rho0 self-consistency matches Python and recovers the density', () => {
    const config = makeConfig({
      ...fixture.config,
      rho0: 0.02, // seeded away from the truth, like the Python golden
      bSqAvg: fixture.fzBSqAvg,
    });
    const estimate = estimateRho0(Q, SQ, config);
    const expected = fixture.expected.rho0Estimate;
    expect(estimate.converged).toBe(expected.converged);
    expect(estimate.iterations).toBe(expected.iterations);
    // The fixed-point iteration compounds summation-order float noise, so
    // cross-engine agreement is bounded by the rtol=1e-3 stopping rule, not
    // by single-pass transform precision.
    expect(relError(estimate.rho0, expected.rho0)).toBeLessThan(1e-4);
    expect(Math.abs(estimate.concordance - 1)).toBeLessThan(1.5e-3);
    expect(relError(estimate.rho0, fixture.config.rho0)).toBeLessThan(0.05); // truth
    expect(() => estimateRho0(Q, SQ, baseConfig())).toThrow(/bSqAvg/);
  });

  it('manual pipeline matches Python sampled outputs', () => {
    const result = scalePipeline(Q, SQ, baseConfig(), fixture.aTrue, fixture.bTrue);
    const expected = fixture.expected.manual;
    expect(relError(result.lowRRms, expected.lowRRms)).toBeLessThan(1e-6);
    expect(relError(result.c1TailMean, expected.c1TailMean)).toBeLessThan(1e-8);
    expected.rSampleIdx.forEach((index, position) => {
      expect(result.gk[index]).toBeCloseTo(expected.gkSamples[position], 9);
    });
    expected.qSampleIdx.forEach((index, position) => {
      expect(result.sqFiltered[index]).toBeCloseTo(expected.sqFilteredSamples[position], 9);
    });
    const summary = diagnosticsSummary(result, baseConfig());
    expect(summary).toHaveProperty('low_r_rms_pre_enforcement');
    expect(summary).toHaveProperty('density_limit_satisfied');
  });

  it('first-peak zeroing enforces the flat level below the cutoff', () => {
    const config = baseConfig();
    const result = scalePipeline(Q, SQ, config, fixture.aTrue, fixture.bTrue);
    const g = firstPeakZero(result.r, result.gFiltered, {
      cutoff: 2.48, peakRmin: 2.65, peakRmax: 3.1,
    });
    for (let i = 0; i < result.r.length; i += 1) {
      if (result.r[i] <= 2.48) expect(g[i]).toBe(0);
    }
  });

  it('two-pass first-shell detection matches Python', () => {
    const config = makeConfig({ ...fixture.config, r0: null, rFitMax: null });
    const result = autoscale(Q, SQ, config);
    const expected = fixture.expected.autoDetect;
    expect(result.r0Detected).toBeCloseTo(expected.r0Detected, 9);
    expect(result.windowRefined).toBe(expected.windowRefined);
    expect(relError(result.a, expected.a)).toBeLessThan(1e-6);
    expect(result.rFitWindowUsed[1]).toBeCloseTo(expected.rFitWindow[1], 9);
    expect(detectFirstPeakOnset(result.r, result.gFiltered, { searchMin: 1.3 }))
      .toBeCloseTo(expected.r0Detected, 9);
  });

  it('density converters follow the ADDIE convention (Mn3Sn check)', () => {
    expect(faberZiman('Mn3Sn').bAvgSqBarn).toBeCloseTo(0.015407, 6);
    const massDensity = massDensityFromNumberDensity('Mn3Sn', 0.063049);
    expect(massDensity).toBeCloseTo(7.4209, 3);
    expect(numberDensityFromMassDensity('Mn3Sn', massDensity)).toBeCloseTo(0.063049, 9);
  });

  it('rejects invalid configurations like the Python engine', () => {
    expect(() => makeConfig({ ...fixture.config, rho0: -1 })).toThrow(/rho0/);
    expect(() => makeConfig({ ...fixture.config, qmin: 40 })).toThrow(/qmax/);
    expect(() => makeConfig({ ...fixture.config, amplitudeCriterion: 'fz' })).toThrow(/bSqAvg/);
    expect(() => makeConfig({
      ...fixture.config, amplitudeCriterion: 'fz', bSqAvg: 0.02, c1Mode: 'joint',
    })).toThrow(/sweep/);
  });
});

describe('stog file parsing (JS ports)', () => {
  const INP = (
    '1\ndata.dat\n1.0 28.0\n-9 0.1\n0\nscale.fq\nscale.gr\n50\n5000\nN\n'
    + '0.063049\n0\nN\nY\n1.0\nscale_ft.sq\nscale_ft.gr\n0.015407\n'
    + 'rmc.fq\nrmc.gr\nrmc.dr\n2.48 2.65 3.1\n'
  );

  it('readStogInp decodes the classic layout', () => {
    const inp = readStogInp(INP);
    expect(inp.a).toBeCloseTo(10, 12);
    expect(inp.b).toBeCloseTo(-9, 12);
    expect(inp.rho0).toBeCloseTo(0.063049, 12);
    expect(inp.peakRmin).toBeCloseTo(2.65, 12);
    expect(inp.lorch).toBe(false);
    expect(() => readStogInp(INP.replace('\nN\nY\n', '\nY\nY\n'))).toThrow(/try again/);
  });

  it('readStogXy keeps NaN padding and picks the modal column count', () => {
    const columns = readStogXy('        3401\n0.956998\n 0.50 NaN\n 0.51 1.25\n 0.52 1.30\ntitle\n');
    expect(columns.length).toBe(2);
    expect(Number.isNaN(columns[1][0])).toBe(true);
    expect(columns[1][2]).toBeCloseTo(1.3, 12);
  });

  it('readDatHeader parses :: metadata', () => {
    const header = readDatHeader(
      'TITLE :: FeCoSn 199K\nNUMBER_DENSITY :: 0.057329 Angstrom^(-3)\nMINIMUM_DISTANCES :: 2.4 2.2\n'
    );
    expect(header.title).toBe('FeCoSn 199K');
    expect(header.numberDensity).toBeCloseTo(0.057329, 12);
    expect(header.minDistance).toBeCloseTo(2.2, 12);
  });

  it('writeStogXy round-trips through readStogXy', () => {
    const x = Float64Array.from([1, 2, 3]);
    const y = Float64Array.from([0.5, -1.25, 3e-4]);
    const text = writeStogXy(x, y, { title: 'test' });
    const columns = readStogXy(text);
    for (let i = 0; i < 3; i += 1) {
      expect(columns[0][i]).toBeCloseTo(x[i], 12);
      expect(columns[1][i]).toBeCloseTo(y[i], 12);
    }
  });

  it('faberZiman matches the pystog argon cross-check', () => {
    expect(faberZiman('Ar').bAvgSqFm2).toBeCloseTo(3.644, 2);
    expect(parseFormula('Sr0.5Ba0.5TiO3').Ti).toBe(1);
    expect(parseFormula('Al2(SO4)3').O).toBe(12);
    expect(() => faberZiman('Xx')).toThrow(/scattering length/);
  });
});
