// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Parity of the static-mode modulation engine against the Python reference
// (rmc_toolkits/modulations.py) via the golden fixture written by
// tests/generate_modulations_fixture.py, plus the worker's request.

import { describe, expect, it } from 'vitest';
import {
    AXIS_TIE_RTOL, DIFFUSE_CORRELATION, DIFFUSE_PHASES, DIFFUSE_RANK_RTOL, DIFFUSE_SIGMA, LONG_WAVE_Q, MODE_KINDS, MODULATION_BECAUSE,
    MODULATION_RULES, MODULATION_SIGMA, NULL_EQUAL_RTOL, SIGMA_MAX, diffuseOrder, displacementModulations, elementSpectrum,
    generalModeSigma, independentCellNull, isLongWave, modePolarisation, nullTail, tensorAt, zoneBoundaryPoints, zoneBoundarySigma
} from '../modulations';
import { evidenceClears as evidenceClearsOf } from '../clusters';
import { siteDisplacementsFromRmc6f } from '../pcaKde';
import { buildCase, buildConfiguration, DEV_CASES } from '../../../validation/cases';
import { TOLERANCES } from '../../../validation/scoring';
import fixture from '../../__tests__/fixtures/modulations_fixture.json';

// Integers, strings and booleans exactly; other floats to 1e-9. A
// polarisation holds floats only: a share or a component that is 1 or 0 in
// Python (a single site's zone-boundary tensor has rank 1: share 1) is a
// float there too, compared to 1e-9.
const expectClose = (actual, expected, path) => {
    if (expected === null || typeof expected !== 'object') {
        const float = typeof expected === 'number' && (!Number.isInteger(expected) || path.includes('.polarisation.'));
        if (float) expect(actual, path).toBeCloseTo(expected, 9);
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

// --- plan item 4.3: the polarisation of a mode ------------------------------------

const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
// The angle between two axes (degrees; an axis has no sign).
const axisAngle = (a, b) => (Math.acos(Math.min(Math.abs(dot(a, b)) / Math.sqrt(dot(a, a) * dot(b, b)), 1)) * 180) / Math.PI;
const cellOf = (parsed) => parsed.latticeVectors.map((row, axis) => row.map((value) => value / parsed.supercell[axis]));
const cartesian = (uvw, cell) => [0, 1, 2].map((k) => uvw[0] * cell[0][k] + uvw[1] * cell[1][k] + uvw[2] * cell[2][k]);
const outer = (v, scale = 1) => v.map((a) => v.map((b) => a * b * scale));
// The stated known answer (tests/test_modulations.py POLARISATION_DEG, _SHARE).
const POLARISATION_DEG = 5;
const POLARISATION_SHARE = 0.95;
// The fixture's planted mode: both sites of its first element along this
// direct-lattice direction of a monoclinic cell, at q = (½, 0, 0).
const PLANTED_UVW = [1, 0, 1];
const POLARISED = fixture.cases.find((entry) => entry.name.startsWith('polarised'));

describe('the polarisation of a mode (plan item 4.3; mode_polarisation is the reference)', () => {
    it('pins the axis-sense tie tolerance to the Python engine\'s', () => {
        expect(AXIS_TIE_RTOL).toBe(fixture.enums.axisTieRtol);
    });

    it('the leading axis as a crystal direction, in the sense that makes the closest crystal axis acute', () => {
        const cell = cellOf(siteDisplacementsFromRmc6f(POLARISED.text));
        const along = cartesian(PLANTED_UVW, cell).map((value, _, v) => value / Math.hypot(...v));
        // [1 0 1] of the monoclinic cell is not the Cartesian (1, 0, 1).
        expect(axisAngle(along, [1, 0, 1])).toBeGreaterThan(13);
        for (const sign of [1, -1]) {
            const pol = modePolarisation(outer(along.map((value) => sign * value), 2.5), cell);
            pol.direction.forEach((value, k) => expect(value).toBeCloseTo(along[k], 12));
            pol.uvw.forEach((value, k) => expect(value).toBeCloseTo(PLANTED_UVW[k], 12));
            expect(pol.share).toBeCloseTo(1, 12);
        }
    });

    it('share: the fraction of the power along the axis; a circular mode has no axis (½); no power, no polarisation', () => {
        const cubic = [[4, 0, 0], [0, 4, 0], [0, 0, 4]];
        const pol = modePolarisation([[0.2, 0, 0], [0, 1, 0], [0, 0, 0.05]], cubic);
        pol.uvw.forEach((value, k) => expect(value).toBeCloseTo([0, 1, 0][k], 12));
        expect(pol.share).toBeCloseTo(1 / 1.25, 12);
        // Re[(a + ib)(a + ib)ᴴ] = aaᵀ + bbᵀ for orthogonal unit a, b.
        expect(modePolarisation([[1, 0, 0], [0, 1, 0], [0, 0, 0]], cubic).share).toBeCloseTo(0.5, 12);
        expect(modePolarisation([[0, 0, 0], [0, 0, 0], [0, 0, 0]], cubic)).toBeNull();
        // Equally close to a and b: the first, a, is made acute.
        const d = [-1, 1, 0].map((value) => value / Math.SQRT2);
        modePolarisation(outer(d), cubic).uvw.forEach((value, k) => expect(value).toBeCloseTo([1, -1, 0][k], 12));
    });

    it('a planted mode in a monoclinic cell is found within the stated angle and share (the fixture\'s case)', () => {
        const parsed = siteDisplacementsFromRmc6f(POLARISED.text);
        const cell = cellOf(parsed);
        const want = cartesian(PLANTED_UVW, cell);
        const result = displacementModulations(parsed);
        const entry = result.elements.find((item) => item.element === 'Zr');
        const point = entry.zoneBoundary.find((item) => item.q.join() === '0.5,0,0');
        expect(point.ordered).toBe(true);
        expect(axisAngle(point.polarisation.direction, want)).toBeLessThan(POLARISATION_DEG);
        expect(axisAngle(cartesian(point.polarisation.uvw, cell), want)).toBeLessThan(POLARISATION_DEG);
        expect(point.polarisation.share).toBeGreaterThan(POLARISATION_SHARE);
        expect(dot(point.polarisation.direction, want)).toBeGreaterThan(0);
        // The tensor's trace is the power: the share is a part of the power at q.
        const { power, tensor } = elementSpectrum(parsed, entry.element);
        const q = (parsed.supercell[0] / 2) * parsed.supercell[1] * parsed.supercell[2];
        const t = tensorAt(tensor, q);
        expect(t[0][0] + t[1][1] + t[2][2]).toBeCloseTo(power[q], 9);
    });

    it('dev-zone-boundary: the planted order along its axis, the decoy long wave along the decoy\'s, within the scoring tolerance', () => {
        // One site over 512 cells: the sampling error is 0.06 / (0.1·√512) rad,
        // 1.5° per component, and this seed's draw is 5.5° off (the PC1 that
        // arm A+ read until now, an estimator of the same noise: 6.0°); the
        // validation scores a direction within 15° (scoring.js).
        const item = DEV_CASES.find((entry) => entry.id === 'dev-zone-boundary');
        const built = buildCase(item, { partials: false });
        const result = displacementModulations(siteDisplacementsFromRmc6f(built.rmc6f));
        const [decoyElement, orderedElement] = item.params.elements;
        const ordered = result.elements.find((entry) => entry.element === orderedElement);
        const point = ordered.zoneBoundary.find((entry) => entry.q.join() === item.params.q.join());
        expect(point.ordered).toBe(true);
        expect(axisAngle(point.polarisation.direction, item.params.axis)).toBeCloseTo(5.52, 2);
        expect(axisAngle(point.polarisation.direction, item.params.axis)).toBeLessThan(TOLERANCES.DIRECTION_DEG);
        expect(point.polarisation.share).toBeGreaterThan(POLARISATION_SHARE);
        const decoy = result.elements.find((entry) => entry.element === decoyElement).strongest[0];
        expect(decoy.kind).toBe('long_wave');
        expect(axisAngle(decoy.polarisation.direction, item.params.decoy.axis)).toBeLessThan(POLARISATION_DEG);
    });
});

// --- F3 (2026-10-19): the family is every element, whichever the call names ------

describe('the Šidák family is every element of the configuration (F3)', () => {
    const atQ = (result, element, q) => result.elements.find((entry) => entry.element === element).zoneBoundary.find((point) => point.q.join() === q.join());

    it('the fixture\'s weak doubling: 2.99σ over every element\'s points, the same point and verdict with Nb named', () => {
        const weak = fixture.cases.find((entry) => entry.name.startsWith('weak doubling'));
        const parsed = siteDisplacementsFromRmc6f(weak.text);
        const every = displacementModulations(parsed);
        const named = displacementModulations(parsed, 'Nb');
        expect(named.family).toEqual(every.family);
        expect(every.family.elements).toEqual(['Ga', 'Nb', 'Se']);
        const point = atQ(every, 'Nb', [0, 0, 0.5]);
        expect(atQ(named, 'Nb', [0, 0, 0.5])).toEqual(point);
        expect(point.sigma).toBeGreaterThan(2.95);
        expect(point.sigma).toBeLessThan(MODULATION_SIGMA);
        expect(point).toMatchObject({ ordered: false, because: 'below_sigma' });
        expect(named.elements[0].verdict.cellDoubling.answer).toBe(false);
        // Nb's own 7 points (the family before F3) would have called it ordered.
        const nb = named.elements[0];
        expect(zoneBoundarySigma(point.ratio * nb.medianPower, independentCellNull(parsed, 'Nb'), 7)).toBeGreaterThanOrEqual(MODULATION_SIGMA);
    });

    it('dev-weak-order-b (F3\'s case): the same point and verdict whether or not the call names its element', () => {
        // The pair re-set at protocol 2026-10-19.19 (its amplitude from the
        // independent-cell null): the weak order reads 2.63σ, below the bar,
        // whether or not the call names Nb. (At F3 it read 2.96σ under the
        // null before F1, at F1 5.6σ; its own 7 points now give 2.99σ.)
        const item = DEV_CASES.find((entry) => entry.id === 'dev-weak-order-b');
        const parsed = siteDisplacementsFromRmc6f(buildCase(item, { partials: false }).rmc6f);
        const element = item.params.elements[1];
        const every = displacementModulations(parsed);
        const named = displacementModulations(parsed, element);
        expect(named.family).toEqual(every.family);
        expect(every.family.zoneBoundaryPoints).toBe(21);
        const point = atQ(named, element, item.params.q);
        expect(point).toEqual(atQ(every, element, item.params.q));
        expect(point.sigma).toBeCloseTo(2.633, 3);
        expect(point).toMatchObject({ ordered: false, because: 'below_sigma' });
        expect(named.elements[0].verdict.cellDoubling).toMatchObject({ answer: false, because: 'no_ordered_point' });
        const power = point.ratio * named.elements[0].medianPower;
        expect(zoneBoundarySigma(power, independentCellNull(parsed, element), 7)).toBeCloseTo(2.987, 3);
        expect(zoneBoundarySigma(power, 1, 21)).toBeCloseTo(0.860, 3);
    });
});

// --- F1 (2026-10-19): the independent-cell null -----------------------------------

describe('the independent-cell null (F1; independent_cell_null and null_tail are the reference)', () => {
    const builtCase = (id) => siteDisplacementsFromRmc6f(buildCase(DEV_CASES.find((entry) => entry.id === id), { partials: false }).rmc6f);

    it('pins the null\'s tolerances to the Python engine\'s', () => {
        expect(NULL_EQUAL_RTOL).toBe(fixture.enums.nullEqualRtol);
        expect(DIFFUSE_RANK_RTOL).toBe(fixture.enums.diffuseRankRtol);
        expect([...DIFFUSE_PHASES]).toEqual(fixture.enums.diffusePhases);
    });

    it('equal weights are the exact gamma; an integer m is the bound the null took until F1', () => {
        expect(nullTail(9, 1, true)).toBeCloseTo(0.0026997960632601866, 15);
        expect(nullTail(9, [0.5, 0.5], false)).toBeCloseTo(Math.exp(-18) * 19, 15);
        expect(nullTail(9, 2, false)).toBe(nullTail(9, [0.5, 0.5], false));
        expect(nullTail(0, 3, true)).toBeNull();
        expect(nullTail(1, [], true)).toBeNull();
    });

    it('the weights sum to the mean power over q ≠ 0, N/(N − 1)', () => {
        const parsed = builtCase('dev-jahn-teller');
        const nQ = parsed.supercell.reduce((product, n) => product * n, 1);
        ['K', 'Cu', 'F'].forEach((name) => {
            const { power } = elementSpectrum(parsed, name);
            const mean = Array.from(power.subarray(1)).reduce((sum, value) => sum + value, 0) / (nQ - 1);
            expect(mean, name).toBeCloseTo(nQ / (nQ - 1), 12);
            expect(independentCellNull(parsed, name).reduce((sum, value) => sum + value, 0), name).toBeCloseTo(nQ / (nQ - 1), 12);
        });
    });

    it('sites that share a cell-level variable read no order, and carry fewer independent components', () => {
        // The six F of a Jahn-Teller octahedron share its long axis, the six
        // W of an octahedral cluster its elongation, the three Mo of a
        // triangle its breathing: their power has about 3, 2.5 and 4.3
        // independent components where independent isotropic sites give 18,
        // 18 and 9, and the null of independent sites assumed 6, 6 and 3 (it
        // found order in such configurations).
        for (const [id, element, sites, most] of [['dev-jahn-teller', 'F', 6, 4], ['dev-jt-off-centre', 'F', 6, 4], ['dev-octahedral-clusters', 'W', 6, 4], ['dev-breathing', 'Mo', 3, 5]]) {
            const result = displacementModulations(builtCase(id));
            const entry = result.elements.find((item) => item.element === element);
            expect(entry.sites, id).toBe(sites);
            expect(entry.nullDof, id).toBeLessThan(most);
            result.elements.forEach((item) => ['cellDoubling', 'generalOrder', 'diffuseOrder'].forEach((name) => {
                expect(item.verdict[name].answer, `${id} ${item.element} ${name}`).not.toBe(true);
            }));
        }
    });

    it('dev-weak-order-a and -b: one isotropic site has three components, so the weak order reads 6.4σ and 5.6σ', () => {
        // F1's known answer, on the pair as it was then (protocol
        // 2026-10-19.16): its amplitude and seeds before the re-set of
        // protocol 2026-10-19.19, which took the amplitude down to an expected
        // 3σ against this null (validation/__tests__/judgment.test.js).
        const atF1 = { 'dev-weak-order-a': { delta: 0.018, seed: 132 }, 'dev-weak-order-b': { delta: 0.018, seed: 133 } };
        for (const [id, sigma] of [['dev-weak-order-a', 6.413], ['dev-weak-order-b', 5.647]]) {
            const current = DEV_CASES.find((entry) => entry.id === id);
            const item = { ...current, params: { ...current.params, ...atF1[id] } };
            const result = displacementModulations(siteDisplacementsFromRmc6f(buildCase(item, { partials: false }).rmc6f));
            const entry = result.elements.find((candidate) => candidate.element === item.params.elements[1]);
            const point = entry.zoneBoundary.find((candidate) => candidate.q.join() === item.params.q.join());
            expect(point.sigma, id).toBeCloseTo(sigma, 2);
            expect(entry.verdict.cellDoubling.decidedBy, id).toEqual([item.params.q]);
            expect(entry.nullDof, id).toBeCloseTo(3, 1);
        }
    });
});

// --- plan item 4.4: diffuse order between cells -------------------------------------

describe('diffuse order between cells (plan item 4.4; diffuse_order is the reference)', () => {
    const fixtureCase = (prefix) => siteDisplacementsFromRmc6f(fixture.cases.find((entry) => entry.name.startsWith(prefix)).text);

    it('the planted short-range order is found, in antiphase along a and in phase along c, with no single mode ordered', () => {
        for (const [prefix, lag, phase] of [['diffuse order in antiphase', [1, 0, 0], 'antiphase'], ['diffuse order in phase', [0, 0, 1], 'in_phase']]) {
            const result = displacementModulations(fixtureCase(prefix));
            const nb = result.elements.find((entry) => entry.element === 'Nb');
            expect(nb.verdict.diffuseOrder).toMatchObject({ answer: true, because: `correlated_${phase}` });
            expect(nb.diffuse).toMatchObject({ lag, phase, df: 13 });
            expect(Math.abs(nb.diffuse.correlation)).toBeGreaterThan(DIFFUSE_CORRELATION);
            expect(nb.diffuse.sigma).toBeGreaterThan(DIFFUSE_SIGMA);
            expect(nb.verdict.cellDoubling.answer).toBe(false);
            expect(nb.verdict.generalOrder.answer).toBe(false);
            const [record] = nb.verdict.diffuseOrder.evidence;
            expect(record.clears).toBe(evidenceClearsOf(record));
        }
    });

    it('no dev case has order between cells, and none reads as diffuse order', () => {
        for (const item of DEV_CASES) {
            const result = displacementModulations(siteDisplacementsFromRmc6f(buildCase(item, { partials: false }).rmc6f));
            result.elements.filter((entry) => entry.sites).forEach((entry) => {
                expect(entry.verdict.diffuseOrder.answer, `${item.id} ${entry.element}`).not.toBe(true);
            });
        }
    });

    it('the statistic reads ranks: any increasing function of the power, and one spiked mode, give the same χ²', () => {
        const q = [];
        for (let h = 0; h < 8; h += 1) for (let k = -3; k < 5; k += 1) for (let l = -3; l < 5; l += 1) q.push([h / 8, k / 8, l / 8]);
        const modes = q.slice(0, 150);
        let state = 7;
        const power = modes.map(() => {
            state = (state * 16807) % 2147483647;
            return -Math.log(state / 2147483647);
        });
        const base = diffuseOrder(power, modes, 3);
        expect(diffuseOrder(power.map((value) => Math.log(value) + 10), modes, 3).chiSquare).toBeCloseTo(base.chiSquare, 9);
        const spiked = [...power];
        spiked[power.indexOf(Math.max(...power))] = 1e6;
        expect(diffuseOrder(spiked, modes, 3).chiSquare).toBeCloseTo(base.chiSquare, 9);
    });
});

// The worker's modulations request (every element by default, 'all' in any
// case, one element, an unknown one) is the request contract,
// fixtures/analysis_requests.json, replayed by src/__tests__/analysisRequests.test.js.
