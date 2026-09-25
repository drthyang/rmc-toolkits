// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// First-coordination-shell handling of the static-mode Auto StoG engine:
// detection (the FIRST shell, either sign), parity with the Python engine on
// golden cases from tests/generate_autoscale_fixture.py.

import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import { autoscale, detectFirstPeakOnset, makeConfig } from '../workers/autoScale';

const relError = (value, reference) => Math.abs(value - reference) / Math.abs(reference);

describe('first-shell detector parity with scaling.detect_first_peak_onset', () => {
  const { r, searchMin, cases } = fixture.expected.detector;
  const rArr = Float64Array.from(r);

  cases.forEach(({ name, qmax, g, onset }) => {
    it(`${name} (qmax ${qmax})`, () => {
      const got = detectFirstPeakOnset(rArr, Float64Array.from(g), { searchMin, qmax });
      if (onset === null) expect(got).toBeNull();
      else expect(got).toBe(onset);
    });
  });

  it('finds an inverted first shell weaker than the second (not the strongest |g|)', () => {
    const inverted = cases.find((item) => item.name === 'invertedFirst' && item.qmax === 28);
    expect(inverted.onset).toBeGreaterThan(1.75);
    expect(inverted.onset).toBeLessThan(1.95);
  });
});

describe('low-r window placement parity with scaling.autoscale', () => {
  const { q, aTrue, cases } = fixture.expected.window;
  const qArr = Float64Array.from(q);

  cases.forEach(({ name, config, sqMeas, expected, error }) => {
    it(`${name}: ${expected ? 'window below the first shell' : 'fails loudly'}`, () => {
      const run = () => autoscale(qArr, Float64Array.from(sqMeas), makeConfig(config));
      if (error) {
        expect(run).toThrow(/first coordination shell/);
        return;
      }
      const result = run();
      expect(relError(result.a, expected.a)).toBeLessThan(1e-6);
      expect(relError(result.b, expected.b)).toBeLessThan(1e-6);
      expect(result.r0Detected).toBe(expected.r0Detected);
      expect(result.windowRefined).toBe(true);
      expect(result.rFitWindowUsed[0]).toBeCloseTo(expected.rFitWindow[0], 12);
      expect(result.rFitWindowUsed[1]).toBeCloseTo(expected.rFitWindow[1], 12);
      // Ti-O at 1.95 A: the window ends below it and the scale is recovered.
      expect(result.rFitWindowUsed[1]).toBeLessThan(1.7);
      expect(relError(result.a, aTrue)).toBeLessThan(0.02);
    });
  });
});
