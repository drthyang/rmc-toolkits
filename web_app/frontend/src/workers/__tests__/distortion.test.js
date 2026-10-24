// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Parity of the static-mode polyhedral-distortion engine (workers/distortion.js)
// against the Python reference (rmc_toolkits/distortion.py) via the golden
// written by tests/generate_distortion_fixture.py, a few known answers of its
// own, and the worker's 'distortion' request.

import { describe, expect, it } from 'vitest';
import {
    APP_MAX_BONDS, OCTAHEDRON_MATCHINGS, PAIRS, QUANTILES, TETRAHEDRAL_ANGLE_DEG, TRANS_TIE_DEG,
    polyhedralDistortion, polyhedronIndices, siteSpecies, summarize
} from '../distortion';
import { WINDOW_TOL } from '../triplets';
import { siteDisplacementsFromRmc6f } from '../pcaKde';
import { handlePcaMessage } from '../pcaKdeWorker';
import { structureFromRmc6f } from '../../browserData';
import { expectGolden } from '../../__tests__/goldenCompare';
import fixture from '../../__tests__/fixtures/distortion_fixture.json';

// The engine's input as the worker builds it from the .rmc6f text: every
// atom (siteDisplacementsFromRmc6f's atomList, both layouts) and the
// average-structure basis (structureFromRmc6f, as /api/structure's).
const configFromText = (text) => {
    const parsed = siteDisplacementsFromRmc6f(text);
    return {
        fractional: parsed.atomList.fractional,
        elements: parsed.atomList.elements,
        latticeVectors: parsed.latticeVectors,
        supercell: parsed.supercell,
        basis: structureFromRmc6f({ text, path: 'fixture.rmc6f' }, 1).basis
    };
};

const outcome = (config, { centre, ligand, rMin, rMax, maxBonds = null }) => {
    try {
        return polyhedralDistortion(config, { centre, ligand, rMin, rMax, maxBonds });
    } catch (error) {
        return { error: error.message };
    }
};

describe('Python parity (distortion_fixture.json)', () => {
    it('shares the engine\'s constants', () => {
        const { constants } = fixture;
        expect(APP_MAX_BONDS).toBe(constants.APP_MAX_BONDS);
        expect(WINDOW_TOL).toBe(constants.WINDOW_TOL);
        expect(TRANS_TIE_DEG).toBe(constants.TRANS_TIE_DEG);
        expectGolden(TETRAHEDRAL_ANGLE_DEG, constants.TETRAHEDRAL_ANGLE_DEG, 'TETRAHEDRAL_ANGLE_DEG');
        expect([...QUANTILES]).toEqual(constants.QUANTILES);
        expect({ 4: PAIRS[4], 6: PAIRS[6] }).toEqual(constants.PAIRS);
        expect(OCTAHEDRON_MATCHINGS).toEqual(constants.OCTAHEDRON_MATCHINGS);
    });

    for (const testCase of fixture.cases) {
        const config = configFromText(testCase.text);
        testCase.queries.forEach((query, index) => {
            const label = `${testCase.name}: ${query.centre}–${query.ligand} ${query.rMin}–${query.rMax} Å${query.maxBonds ? ` (≤ ${query.maxBonds} bonds)` : ''}`;
            it(label, () => {
                const expected = testCase.expected[index];
                const got = outcome(config, query);
                if (expected.error) expect(got).toEqual(expected);
                else expectGolden(got, expected, 'result');
            });
        });
    }
});

// Bonds as the neighbour search gives them, from Cartesian vectors.
const bondsOf = (vectors) => vectors.map(([vx, vy, vz]) => ({ vx, vy, vz, length: Math.hypot(vx, vy, vz) }));
const AXES = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]];

describe('known answers', () => {
    it('a regular octahedron and a regular tetrahedron are 0, 1 and 0', () => {
        const octahedron = polyhedronIndices(bondsOf(AXES.map((axis) => axis.map((value) => 2 * value))));
        expect(octahedron.distortionIndex).toBeLessThan(1e-14);
        expect(octahedron.quadraticElongation).toBeCloseTo(1, 13);
        expect(octahedron.angleVariance).toBeLessThan(1e-20);
        expect(octahedron.transAngle).toBeCloseTo(180, 10);
        const s = 2 / Math.sqrt(3);
        const tetrahedron = polyhedronIndices(bondsOf([[s, s, s], [s, -s, -s], [-s, s, -s], [-s, -s, s]]));
        expect(tetrahedron.distortionIndex).toBeLessThan(1e-14);
        expect(tetrahedron.quadraticElongation).toBeCloseTo(1, 13);
        expect(tetrahedron.angleVariance).toBeLessThan(1e-20);
        expect(tetrahedron.transAngle).toBeNull();
    });

    it('a 4 + 2 Jahn–Teller octahedron at 1.93 / 2.25 Å: DI = 0.0698, ⟨λ⟩ in closed form, σ²θ = 0', () => {
        const lengths = [1.93, 1.93, 1.93, 1.93, 2.25, 2.25];
        const indices = polyhedronIndices(bondsOf(AXES.map((axis, k) => axis.map((value) => value * lengths[k]))));
        const mean = 12.22 / 6;
        expect(indices.distortionIndex).toBeCloseTo((4 * Math.abs(1.93 - mean) + 2 * Math.abs(2.25 - mean)) / (6 * mean), 12);
        expect(indices.distortionIndex).toBeCloseTo(0.0698308783, 9);
        expect(indices.quadraticElongation).toBeCloseTo((4 * 1.93 ** 2 + 2 * 2.25 ** 2) / (6 * (1.93 ** 2 * 2.25) ** (2 / 3)), 12);
        expect(indices.angleVariance).toBeLessThan(1e-20);
    });

    it('a flat polyhedron has no quadratic elongation', () => {
        const square = polyhedronIndices(bondsOf([[2, 0, 0], [0, 2, 0], [-2, 0, 0], [0, -2, 0]]));
        expect(Number.isNaN(square.quadraticElongation)).toBe(true);
        expect(summarize([square.quadraticElongation])).toBeNull();
    });

    it('summaries and site labels follow the Python rules', () => {
        expect(summarize([4, 1, 3, 2, 5])).toEqual({ count: 5, mean: 3, sd: Math.sqrt(2.5), p10: expect.closeTo(1.4, 14), median: 3, p90: expect.closeTo(4.6, 14) });
        expect(summarize([1]).sd).toBeNull();
        expect(siteSpecies({ Zr: 3, Ti: 7 })).toBe('Ti');
        expect(siteSpecies({ Zr: 5, Ti: 5 })).toBe('Ti');
    });
});

describe('the worker\'s distortion request', () => {
    const text = fixture.cases.find((entry) => entry.name === 'perovskite-thermal').text;
    const ask = (params) => handlePcaMessage({ kind: 'distortion', ...params }, async () => text);

    it('answers the engine\'s payload with the parse warning, under the request caps', async () => {
        const result = await ask({ centre: 'ti', ligand: 'o', rMin: 1.5, rMax: 2.5 });
        const { parseWarning, ...payload } = result;
        expect(parseWarning).toBeNull();
        expectGolden(payload, fixture.cases.find((entry) => entry.name === 'perovskite-thermal').expected[0], 'result');
        await expect(ask({ centre: 'Ti', ligand: 'O', rMin: 1.5, rMax: 8.5 })).rejects.toThrow('rMax must be <= 8');
        await expect(ask({ centre: 'Ti', ligand: 'O', rMin: 1.5 })).rejects.toThrow('rMax required');
        await expect(ask({ centre: 'Ti', ligand: 'O', rMin: 3, rMax: 2 })).rejects.toThrow('window must be');
    });
});
