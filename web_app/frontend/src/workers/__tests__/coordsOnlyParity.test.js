// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Cross-runtime parity on a legacy coords-only .rmc6f (`id element [label]
// x y z`: no reference-site or cell-index columns). The golden
// (rmc6f_coords_only_fixture.json, written by tests/generate_coords_only_fixture.py)
// holds the configuration and what the reference-grade Python FILE loaders
// make of it: kde.load_unit_cell_positions (the Flask /api/kde/slice input)
// and triplets.bond_angle_summary_from_file (the /api/triplets body).
// Here the browser path reads the same text: structureFromRmc6f (whose
// points feed the KDE worker) and the shared worker's 'triplets' request.
// Both runtimes must see every atom -- the Flask KDE used to see none.

import { describe, expect, it } from 'vitest';
import { structureFromRmc6f } from '../../browserData.js';
import { handlePcaMessage } from '../pcaKdeWorker.js';
import fixture from '../../__tests__/fixtures/rmc6f_coords_only_fixture.json';

describe('coords-only .rmc6f: browser vs Python file loaders', () => {
    it('folds the same atoms to the same unit-cell positions as the Flask KDE loader', () => {
        const structure = structureFromRmc6f({ text: fixture.rmc6f, path: 'coords_only.rmc6f' }, 1_000_000);
        const expected = fixture.kdePositions.fractional;
        expect(expected.length).toBeGreaterThan(0);
        expect(structure.parseReport.coordsOnlyAtoms).toBe(expected.length);
        expect(structure.parseWarning).toBeNull();
        expect(structure.points).toHaveLength(expected.length);
        let worst = 0;
        structure.points.forEach((point, index) => {
            [point.x, point.y, point.z].forEach((value, axis) => {
                worst = Math.max(worst, Math.abs(value - expected[index][axis]));
            });
        });
        expect(worst).toBeLessThan(1e-12);
    });

    for (const { request, summary } of fixture.triplets) {
        const label = `${request.end1}-${request.apex}-${request.end2}`
            + (request.r23Min != null ? ' (distinct windows)' : '');
        it(`computes the same ${label} bond-angle summary as /api/triplets`, async () => {
            const result = await handlePcaMessage({ kind: 'triplets', ...request }, async () => fixture.rmc6f);
            expect(summary.angleCount).toBeGreaterThan(0);
            expect(result.angleCount).toBe(summary.angleCount);
            expect(result.apexCount).toBe(summary.apexCount);
            expect(result.counts).toEqual(summary.counts);
            expect(result.coordination).toEqual(summary.coordination);
            expect(result.lengths12.counts).toEqual(summary.lengths12.counts);
            expect(result.lengths12.uniqueBonds).toBe(summary.lengths12.uniqueBonds);
            if (summary.lengths23) {
                expect(result.lengths23.counts).toEqual(summary.lengths23.counts);
                expect(result.lengths23.uniqueBonds).toBe(summary.lengths23.uniqueBonds);
            }
            expect(result.sharedEnds).toBe(summary.sharedEnds);
            expect(result.meanAngle).toBeCloseTo(summary.meanAngle, 9);
            expect(result.stdAngle).toBeCloseTo(summary.stdAngle, 9);
            summary.sinCorrected.forEach((value, bin) => {
                if (value === null) expect(result.sinCorrected[bin]).toBeNull();
                else expect(result.sinCorrected[bin]).toBeCloseTo(value, 9);
            });
            expect(result.parseWarning).toBe(summary.parseWarning);
        });
    }
});
