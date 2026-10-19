// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The browser displacement-direction engine (workers/orientation.js) against
// the Python reference (rmc_toolkits/orientation.py), through the golden that
// tests/generate_orientation_fixture.py writes -- one written interface for
// both engines instead of values copied into two suites. Python checks the
// committed file against its own engine (tests/test_orientation_fixture.py).
//
// Tolerances (goldenCompare.js): integers (counts, cells, frequencies, tie
// counts), booleans, strings and nulls exactly; every other number to 1e-9
// relative + 1e-12 absolute, p-values to 1e-9 relative alone; the axis sets
// (pcaAxes, orientationAxes) to 1e-9 eigenspace by eigenspace (an axis up to
// its sign, a degenerate pair -- the 8-well model's -- by its projector), with
// the eigenvalues the fixture's `spectra` give.

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';

import {
    assignCells,
    goldbergTiling,
    orientationHistogram,
    pearsonNullMoments,
    recommendedFrequency,
    siteOrientationHistogram
} from '../workers/orientation.js';
import { siteDisplacementsFromRmc6f } from '../workers/pcaKde.js';
import fixture from './fixtures/orientation_fixture.json';
import { expectGolden } from './goldenCompare.js';

const TILING_KEYS = ['centers', 'areas', 'sizes', 'antipode'];
const PER_CELL_KEYS = ['mass', 'density', 'enhancement', 'expected', 'zScore', 'cellMeanAmplitude'];
// The browser engine marks its payload; Python has no such key.
const JS_ONLY_KEYS = ['browserOrientation'];
const RELATIVE = new Set(['peakLocalPValue', 'peakPValue', 'mapPValue', 'orientationBinghamPValue']);

const demoCache = new Map();
const demoSites = (name) => {
    if (!demoCache.has(name)) {
        const text = readFileSync(new URL(`../../public/demo/${name}`, import.meta.url), 'utf8');
        demoCache.set(name, siteDisplacementsFromRmc6f(text));
    }
    return demoCache.get(name);
};

const vectorsOf = (spec) => spec.vectors.flatMap(([name, count]) => (
    count === null ? fixture.inputs[name] : fixture.inputs[name].slice(0, count)
));

const run = ({ input, options }) => {
    const merged = { geometry: false, ...options };
    if (input.demo) {
        return siteOrientationHistogram(demoSites(input.demo), {
            referenceNumber: input.referenceNumber ?? null,
            element: input.element ?? null,
            ...merged
        });
    }
    return orientationHistogram(vectorsOf(input), merged);
};

describe('orientation payload parity (orientation_fixture.json)', () => {
    it('covers every kind of case the fixture promises', () => {
        const names = fixture.histograms.map((entry) => entry.name);
        expect(names.length).toBeGreaterThanOrEqual(20);
        expect(fixture.histograms.some((entry) => entry.expected.mapSignificance === null)).toBe(true);
        expect(fixture.histograms.some((entry) => entry.expected.antipodalAsymmetryZ === null)).toBe(true);
        expect(fixture.histograms.some((entry) => entry.expected.antipodalAsymmetrySignificant)).toBe(true);
        expect(fixture.histograms.some((entry) => entry.input.demo)).toBe(true);
        expect(fixture.histograms.some((entry) => entry.options.geometry)).toBe(true);
    });

    for (const entry of fixture.histograms) {
        it(entry.name, () => {
            const result = run(entry);
            JS_ONLY_KEYS.forEach((key) => expect(result[key], key).toBe(true));
            const keys = Object.keys(result).filter((key) => !JS_ONLY_KEYS.includes(key)).sort();
            expect(keys).toEqual(entry.keys);

            const tiling = fixture.tilings[String(entry.expected.frequency)];
            const options = { relative: RELATIVE, bases: entry.spectra };
            entry.keys.forEach((key) => {
                if (key in entry.expected) {
                    expectGolden(result[key], entry.expected[key], key, { ...options, key });
                } else if (TILING_KEYS.includes(key) && entry.cells) {
                    expectGolden(result[key], tiling[key], key, { ...options, key });
                } else {
                    // Stored once per frequency or only for the `cells` cases:
                    // here only its shape.
                    expect(TILING_KEYS.includes(key) || PER_CELL_KEYS.includes(key), key).toBe(true);
                    expect(result[key], key).toHaveLength(entry.expected.cellCount);
                }
            });
        });
    }
});

describe('tiling parity (orientation_fixture.json)', () => {
    it('assigns the exact tie families to the same cells', () => {
        fixture.assignments.forEach(({ name, directions, cells }) => {
            Object.entries(cells).forEach(([frequency, expected]) => {
                expect(assignCells(goldbergTiling(Number(frequency)), directions), `${name} nu=${frequency}`)
                    .toEqual(expected);
            });
        });
    });

    it('orders every neighbour row as Python does (checksums)', () => {
        // sum over cells c and slots k of (k + 1)(neighbors[c][k] + 1)(c % 97 + 1)
        const checksum = (neighbors) => neighbors.reduce((total, row, cell) => total + row.reduce(
            (sum, neighbor, slot) => sum + (slot + 1) * (neighbor + 1) * ((cell % 97) + 1), 0
        ), 0);
        Object.entries(fixture.tilingChecks.neighborChecksums).forEach(([frequency, expected]) => {
            expect(checksum(goldbergTiling(Number(frequency)).neighbors), `nu=${frequency}`).toBe(expected);
        });
    });

    it('starts the polygons on the same vertex', () => {
        Object.entries(fixture.tilingChecks.polygonStarts).forEach(([frequency, starts]) => {
            const tiling = goldbergTiling(Number(frequency));
            Object.entries(starts).forEach(([cell, vertex]) => {
                expectGolden(tiling.polygons[Number(cell)][0], vertex, `nu=${frequency} cell ${cell}`);
            });
        });
    });
});

describe('statistics helpers parity (orientation_fixture.json)', () => {
    it('recommends the same frequency at every boundary', () => {
        fixture.recommendedFrequency.forEach(({ n, targetPerCell, maxFrequency, frequency }) => {
            expect(recommendedFrequency(n, { targetPerCell, maxFrequency }), `N=${n}`).toBe(frequency);
        });
    });

    it('has the same exact null moments of Pearson\'s X^2', () => {
        fixture.pearsonMoments.forEach(({ used, probabilities, moments }) => {
            expectGolden(pearsonNullMoments(used, probabilities), moments, `N=${used} C=${probabilities.length}`);
        });
    });
});
