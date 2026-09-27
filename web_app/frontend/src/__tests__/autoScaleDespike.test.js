// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Opt-in despike parity (0.6.0 audit, stog-b): Python used to despike twice in
// the auto path (fit on one point set, write and count another); both engines
// now despike once, so (a, b), the written grid and nDespiked agree.

import { describe, expect, it } from 'vitest';
import fixture from './fixtures/autoscale_fixture.json';
import { autoscale, cropSq, makeConfig } from '../workers/autoScale';

const relError = (value, reference) => Math.abs(value - reference) / Math.abs(reference);

describe('despike (single pass) matches Python', () => {
  it('fits, writes and counts the same point set', () => {
    const expected = fixture.expected.despike;
    const q = Float64Array.from(fixture.q);
    const sq = Float64Array.from(expected.sqMeas);
    const config = makeConfig({ ...fixture.config, despike: true });
    const result = autoscale(q, sq, config);
    expect(result.q.length).toBe(expected.nq);
    expect(result.nDespiked).toBe(expected.nDespiked);
    expect(result.nDespiked).toBe(12); // exactly the injected glitches
    expect(result.q.length).toBe(cropSq(q, sq, config).q.length);
    expect(result.iterations).toBe(expected.iterations);
    expect(relError(result.a, expected.a)).toBeLessThan(1e-10);
    expect(relError(result.b, expected.b)).toBeLessThan(1e-10);
    expect(relError(result.lowRRms, expected.lowRRms)).toBeLessThan(1e-10);
  });
});
