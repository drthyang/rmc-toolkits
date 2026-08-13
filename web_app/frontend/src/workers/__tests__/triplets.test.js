// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Parity of the static-mode triplets engine against the Python reference
// (rmc_toolkits/triplets.py) via the golden fixture written by
// tests/generate_triplets_fixture.py, plus self-contained geometry checks.

import { describe, expect, it } from 'vitest';
import { bondAngleSummary } from '../triplets';
import fixture from '../../__tests__/fixtures/triplets_fixture.json';

const runSpec = (testCase, spec) =>
  bondAngleSummary(testCase.fractional, testCase.elements, testCase.lattice, {
    triplet: spec.triplet,
    bond12: spec.bond12,
    bond23: spec.bond23,
    binWidth: spec.binWidth,
    collectAngles: true
  });

describe('Python parity (triplets_fixture.json)', () => {
  for (const testCase of fixture.cases) {
    for (const [index, spec] of testCase.specs.entries()) {
      it(`${testCase.name} ${spec.triplet.join('-')} #${index}`, () => {
        const expected = testCase.expected[index];
        const result = runSpec(testCase, spec);

        // Histogram counts are integers: exact equality, no tolerance.
        expect(result.counts).toEqual(expected.counts);
        expect(result.angleCount).toBe(expected.angleCount);
        expect(result.apexCount).toBe(expected.apexCount);
        expect(result.lengths12.count).toBe(expected.bond12Count);
        const count23 = result.sharedEnds ? result.lengths12.count : result.lengths23.count;
        expect(count23).toBe(expected.bond23Count);

        // Single-pass float math: 1e-9 covers summation-order noise.
        for (const [key, values] of [['density', expected.density], ['sinCorrected', expected.sinCorrected]]) {
          values.forEach((value, bin) => expect(result[key][bin]).toBeCloseTo(value, 9));
        }
        expect(result.meanAngle).toBeCloseTo(expected.meanAngle, 7);
        expect(result.stdAngle).toBeCloseTo(expected.stdAngle, 7);
        expect(result.lengths12.meanLength).toBeCloseTo(expected.meanLength12, 7);

        // Sorted angles: full list for small cases, ends for large ones.
        // The fixture rounds to 1e-6 deg, so compare at 1e-5.
        const sorted = result.sortedAngles;
        const near = (actual, reference) =>
          reference.forEach((value, i) => expect(actual[i]).toBeCloseTo(value, 5));
        if (expected.sortedAngles) {
          expect(sorted.length).toBe(expected.sortedAngles.length);
          near(sorted, expected.sortedAngles);
        } else {
          near(sorted.slice(0, expected.sortedAnglesHead.length), expected.sortedAnglesHead);
          near(sorted.slice(-expected.sortedAnglesTail.length), expected.sortedAnglesTail);
        }

        // The full app payload: lengths + coordination parity.
        const summary = expected.summary;
        expect(result.coordination).toEqual(summary.coordination);
        expect(result.lengths12.counts).toEqual(summary.lengths12.counts);
        if (summary.lengths23) {
          expect(result.lengths23.counts).toEqual(summary.lengths23.counts);
        } else {
          expect(result.lengths23).toBeNull();
        }
        summary.binCenters.forEach((value, bin) =>
          expect(result.binCenters[bin]).toBeCloseTo(value, 9));
        expect(result.sharedEnds).toBe(summary.sharedEnds);
      });
    }
  }
});

describe('geometry invariants', () => {
  const cubic = [[10, 0, 0], [0, 10, 0], [0, 0, 10]];

  it('octahedron gives 12 right angles and 3 straight ones', () => {
    const fractional = [
      [0.5, 0.5, 0.5],
      [0.7, 0.5, 0.5], [0.3, 0.5, 0.5],
      [0.5, 0.7, 0.5], [0.5, 0.3, 0.5],
      [0.5, 0.5, 0.7], [0.5, 0.5, 0.3]
    ];
    const elements = ['Nb', 'O', 'O', 'O', 'O', 'O', 'O'];
    const result = bondAngleSummary(fractional, elements, cubic, {
      triplet: ['O', 'Nb', 'O'],
      bond12: [1, 3],
      collectAngles: true
    });
    expect(result.angleCount).toBe(15);
    expect(result.coordination).toEqual([0, 0, 0, 0, 0, 0, 1]);
    const sorted = result.sortedAngles;
    sorted.slice(0, 12).forEach((angle) => expect(angle).toBeCloseTo(90, 9));
    sorted.slice(12).forEach((angle) => expect(angle).toBeCloseTo(180, 9));
  });

  it('bonds cross the periodic wall and images of one atom both count', () => {
    const result = bondAngleSummary(
      [[0, 0, 0], [0.25, 0, 0]],
      ['Nb', 'Se'],
      [[4, 0, 0], [0, 4, 0], [0, 0, 4]],
      { triplet: ['Se', 'Nb', 'Se'], bond12: [0.5, 3.5], collectAngles: true }
    );
    expect(result.lengths12.count).toBe(2);
    expect(result.angleCount).toBe(1);
    expect(result.sortedAngles[0]).toBeCloseTo(180, 9);
  });

  it('rejects an unknown element with the available list', () => {
    expect(() =>
      bondAngleSummary([[0.5, 0.5, 0.5]], ['Nb'], cubic, {
        triplet: ['Se', 'Nb', 'Se'],
        bond12: [1, 2]
      })
    ).toThrow(/available: Nb/);
  });

  it('normalization: density integrates to 1, counts partition the angles', () => {
    const rng = (() => {
      let state = 42;
      return () => {
        state = (state * 1664525 + 1013904223) % 4294967296;
        return state / 4294967296;
      };
    })();
    const fractional = Array.from({ length: 30 }, () => [rng(), rng(), rng()]);
    const elements = fractional.map((_, index) => (index % 3 ? 'Se' : 'Nb'));
    const result = bondAngleSummary(fractional, elements, cubic, {
      triplet: ['Se', 'Nb', 'Se'],
      bond12: [1, 4],
      binWidth: 2
    });
    const total = result.counts.reduce((acc, value) => acc + value, 0);
    expect(total).toBe(result.angleCount);
    const integral = result.density.reduce((acc, value) => acc + value * result.binWidth, 0);
    expect(integral).toBeCloseTo(total ? 1 : 0, 9);
    const apexTotal = result.coordination.reduce((acc, value) => acc + value, 0);
    expect(apexTotal).toBe(result.apexCount);
  });
});
