// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Parity of the static-mode modulation engine against the Python reference
// (rmc_toolkits/modulations.py) via the golden fixture written by
// tests/generate_modulations_fixture.py, plus the worker's request.

import { describe, expect, it } from 'vitest';
import {
    LONG_WAVE_Q, MODE_KINDS, MODULATION_BECAUSE, MODULATION_RULES, SIGMA_MAX,
    displacementModulations, generalModeSigma, isLongWave, zoneBoundaryPoints, zoneBoundarySigma
} from '../modulations';
import { siteDisplacementsFromRmc6f } from '../pcaKde';
import { buildCase, buildConfiguration, DEV_CASES } from '../../../validation/cases';
import fixture from '../../__tests__/fixtures/modulations_fixture.json';

// Integers, strings and booleans exactly; other floats to 1e-9.
const expectClose = (actual, expected, path) => {
    if (expected === null || typeof expected !== 'object') {
        if (typeof expected === 'number' && !Number.isInteger(expected)) expect(actual, path).toBeCloseTo(expected, 9);
        else expect(actual, path).toEqual(expected);
        return;
    }
    if (Array.isArray(expected)) {
        expect(actual, path).toHaveLength(expected.length);
        expected.forEach((value, index) => expectClose(actual[index], value, `${path}[${index}]`));
        return;
    }
    expect(Object.keys(actual).sort(), path).toEqual(Object.keys(expected).sort());
    Object.entries(expected).forEach(([key, value]) => expectClose(actual[key], value, `${path}.${key}`));
};

describe('Python parity (modulations_fixture.json)', () => {
    it('halves only even axes', () => {
        expect(zoneBoundaryPoints([4, 4, 4])).toHaveLength(7);
        expect(zoneBoundaryPoints([3, 4, 5])).toEqual([[0, 1, 0]]);
    });

    for (const testCase of fixture.cases) {
        const parsed = siteDisplacementsFromRmc6f(testCase.text);
        testCase.elements.forEach((element, index) => {
            it(`${testCase.name}: ${element ?? 'every element'}`, () => {
                expectClose(displacementModulations(parsed, element), testCase.expected[index], 'result');
            });
        });
    }
});

describe('the verdict (the engine decides; modulations.py is the reference)', () => {
    it('pins the closed vocabularies and the rules to the Python engine\'s (the fixture\'s enums)', () => {
        expect(MODULATION_BECAUSE).toEqual(fixture.enums.because);
        expect([...MODE_KINDS]).toEqual(fixture.enums.kinds);
        expect({ ...MODULATION_RULES }).toEqual(fixture.enums.rules);
    });

    it('every fixture point: ordered recomputable from its ratio, σ and bar; because from its closed set', () => {
        fixture.cases.forEach((testCase) => testCase.expected.forEach((result) => result.elements.forEach((entry) => {
            Object.entries(entry.verdict).forEach(([name, finding]) => expect(MODULATION_BECAUSE[name]).toContain(finding.because));
            [...(entry.zoneBoundary || []), ...(entry.strongest || [])].forEach((point) => {
                expect(MODULATION_BECAUSE.point).toContain(point.because);
                if (point.bar == null || point.because === 'undefined_statistic') {
                    expect(point.ordered).toBeNull();
                    return;
                }
                expect(point.ordered).toBe(point.ratio >= point.bar.ratio && point.sigma != null && point.sigma >= point.bar.sigma);
            });
        })));
    });

    it('zone-boundary σ: calibrated against independent displacements, corrected over the points tested', () => {
        // Power 1 is the independent-displacement mean: nothing.
        expect(zoneBoundarySigma(1, 1, 7)).toBeLessThan(0);
        // One site, one point: P(χ²₁ ≥ 9) = 0.0027, 2.78σ; over seven points 1 − (1 − p)⁷ = 0.019, 2.08σ.
        expect(zoneBoundarySigma(9, 1, 1)).toBeCloseTo(2.78, 2);
        expect(zoneBoundarySigma(9, 1, 7)).toBeCloseTo(2.08, 2);
        // More sites average the noise down: the same power is stronger evidence.
        expect(zoneBoundarySigma(3, 16, 7)).toBeGreaterThan(zoneBoundarySigma(3, 1, 7));
        expect(zoneBoundarySigma(250, 1, 21)).toBeGreaterThan(8);
        expect(zoneBoundarySigma(null, 1, 7)).toBeNull();
        expect(zoneBoundarySigma(0, 1, 7)).toBeNull();
    });

    it('general-mode σ: P ≥ χ²_2m / 2m, corrected over the compared modes', () => {
        // One site, one mode: P(χ²₂ ≥ 18) = e⁻⁹ = 1.23e-4, 3.67σ; corrected over 255 modes: 0.031, 1.87σ.
        expect(generalModeSigma(9, 1, 1)).toBeCloseTo(3.666, 2);
        expect(generalModeSigma(9, 1, 255)).toBeCloseTo(1.867, 2);
        expect(generalModeSigma(1, 1, 7)).toBeLessThan(0);
        // Power 1 over 255 modes: 1 − (1 − e⁻¹)²⁵⁵ rounds to 1, σ = −∞: reported as null.
        expect(generalModeSigma(1, 1, 255)).toBeNull();
        expect(generalModeSigma(9, 1, 0)).toBeNull();
    });

    it('σ is reported as the tool showed it: +∞ is SIGMA_MAX (shown >8), −∞ is null', () => {
        expect(SIGMA_MAX).toBe(40);
        expect(zoneBoundarySigma(1e5, 4, 7)).toBe(SIGMA_MAX);
        expect(generalModeSigma(1e5, 4, 99)).toBe(SIGMA_MAX);
        expect(zoneBoundarySigma(650, 2, 1)).toBeLessThan(SIGMA_MAX);
        expect(zoneBoundarySigma(0.01, 32, 21)).toBeNull();
    });

    it('one long-wave predicate: on every grid point it is the tool\'s |q_i| ≤ 0.25 and the engine\'s |q_i| ≤ 0.25 + 1e-12', () => {
        for (let size = 1; size <= 256; size += 1) {
            for (let index = 0; index < size; index += 1) {
                const value = index / size > 0.5 ? index / size - 1 : index / size;
                expect(Math.abs(value) <= LONG_WAVE_Q + 1e-12).toBe(Math.abs(value) <= LONG_WAVE_Q);
                expect(isLongWave([value, 0, 0])).toBe(Math.abs(value) <= LONG_WAVE_Q);
                // A zone-boundary component by its index (the engine) is the tool's 1e-9 test on the value.
                expect(index === 0 || 2 * index === size).toBe(Math.abs(value) < 1e-9 || Math.abs(Math.abs(value) - 0.5) < 1e-9);
            }
        }
    });

    it('order at a general wavevector (3/8, 0, 0) on an 8-cell axis: generalOrder, not a cell doubling (review L4)', () => {
        const config = buildConfiguration({
            cell: [4, 4, 4], supercell: [8, 8, 8], seed: 11,
            sites: [{ el: 'Ti', frac: [0.5, 0.5, 0.5] }, { el: 'O', frac: [0.5, 0.5, 0] }],
            shift: (index, _, rng, n) => {
                const thermal = [0, 1, 2].map(() => 0.06 * Math.sqrt(-2 * Math.log(rng() || 1e-12)) * Math.cos(2 * Math.PI * rng()));
                return index === 0 ? [thermal[0], thermal[1], thermal[2] + 0.08 * Math.cos(2 * Math.PI * 0.375 * n[0])] : thermal;
            }
        });
        const result = displacementModulations(siteDisplacementsFromRmc6f(config.text));
        const ti = result.elements.find((entry) => entry.element === 'Ti');
        expect(ti.verdict.generalOrder).toMatchObject({ answer: true, because: 'ordered_mode', decidedBy: [[expect.closeTo(0.375, 6), 0, 0]] });
        expect(ti.verdict.cellDoubling.answer).toBe(false);
        expect(ti.strongest[0]).toMatchObject({ kind: 'general', ordered: true, because: 'ordered' });
        expect(result.elements.find((entry) => entry.element === 'O').verdict.generalOrder.answer).toBe(false);
    });
});

describe('the verdict on the validation\'s dev cases (the decisions the tool used to make)', () => {
    const verdictsOf = (id) => {
        const built = buildCase(DEV_CASES.find((item) => item.id === id), { partials: false });
        return displacementModulations(siteDisplacementsFromRmc6f(built.rmc6f));
    };

    it('independent or randomly split displacements are not order (the earlier rule called them ordered)', () => {
        // Ba at (½, 0, 0) in pure thermal motion, Nb in random two-site
        // disorder, Cu in randomly oriented JT octahedra, Mo in breathing triangles.
        for (const id of ['dev-thermal-null', 'dev-offcentre', 'dev-jahn-teller', 'dev-breathing']) {
            verdictsOf(id).elements.filter((entry) => entry.sites).forEach((entry) => {
                expect(entry.verdict.cellDoubling.answer, `${id} ${entry.element}`).not.toBe(true);
                expect(entry.verdict.generalOrder.answer, `${id} ${entry.element}`).not.toBe(true);
            });
        }
    });

    it('dev-zone-boundary: the planted Ti order at (0, 0, ½) is found, and only it; the decoy Sr wave is a long wave', () => {
        const result = verdictsOf('dev-zone-boundary');
        const ordered = result.elements.flatMap((entry) => entry.verdict.cellDoubling.decidedBy.map((q) => [entry.element, q]));
        expect(ordered).toEqual([['Ti', [0, 0, 0.5]]]);
        expect(result.elements.find((entry) => entry.element === 'Sr').verdict.longWave.answer).toBe(true);
        expect(result.elements.every((entry) => entry.verdict.generalOrder.answer !== true)).toBe(true);
    });
});

// The worker's modulations request (every element by default, 'all' in any
// case, one element, an unknown one) is the request contract,
// fixtures/analysis_requests.json, replayed by src/__tests__/analysisRequests.test.js.
