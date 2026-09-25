// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// First-coordination-shell handling of the static-mode Auto StoG engine:
// detection (the FIRST shell, either sign), parity with the Python engine on
// golden cases from tests/generate_autoscale_fixture.py.

import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import { detectFirstPeakOnset } from '../workers/autoScale';

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
