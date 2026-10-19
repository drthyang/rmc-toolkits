// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The browser per-site displacement statistics (workers/pcaKde.js,
// siteEllipsoids: the table of the PCA Ellipsoid and Displacement Directions
// pages) against the Python reference (rmc_toolkits/pca_kde.py,
// site_ellipsoids) through the golden tests/generate_pca_sites_fixture.py
// writes. Python checks the committed file against its own engine
// (tests/test_pca_sites_fixture.py).
//
// Tolerances (goldenCompare.js): integers (counts, element counts, reference
// numbers), booleans (axisResolved, mixed, zeroSpread, degenerate), strings
// and nulls (the zero-spread site's axes, kurtosis and anisotropy) exactly;
// every other number to 1e-9 relative + 1e-12 absolute (A^2 for the
// covariance); each site's axes to 1e-9 eigenspace by eigenspace with its own
// eigenvalues (an axis up to its sign).

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';

import { siteDisplacementsFromRmc6f, siteEllipsoids, siteSpecies } from '../workers/pcaKde.js';
import fixture from './fixtures/pca_sites_fixture.json';
import { expectGolden } from './goldenCompare.js';

// The browser row also carries the reconstruction path's copy count (null on
// the reference-number path every case here takes); Python has no such key.
const JS_ONLY_KEYS = ['copiesPerCell'];

const textOf = (source) => (source.demo
    ? readFileSync(new URL(`../../public/demo/${source.demo}`, import.meta.url), 'utf8')
    : fixture.synthetic);

describe('site statistics parity (pca_sites_fixture.json)', () => {
    it('covers the site kinds the fixture promises', () => {
        const rows = fixture.cases.flatMap((entry) => entry.expected);
        expect(rows.some((row) => row.zeroSpread)).toBe(true);
        expect(rows.some((row) => row.mixed)).toBe(true);
        expect(rows.some((row) => row.nonGaussianity !== null && row.nonGaussianity < 0)).toBe(true);
        expect(rows.some((row) => row.axisResolved.every(Boolean))).toBe(true);
        expect(rows.some((row) => !row.zeroSpread && row.axisResolved.every((value) => !value))).toBe(true);
        expect(Math.max(...rows.map((row) => row.count))).toBeLessThan(20000);
    });

    for (const entry of fixture.cases) {
        it(entry.name, () => {
            const parsed = siteDisplacementsFromRmc6f(textOf(entry.source));
            expect(parsed.reconstructed).toBe(false);
            expect(parsed.referenceNumbers).toEqual(entry.referenceNumbers);
            expect(siteSpecies(parsed.sites)).toEqual(entry.species);
            expect(parsed.parseWarning).toBe(entry.parseWarning);

            const rows = siteEllipsoids(parsed.sites, entry.probability);
            expect(rows).toHaveLength(entry.expected.length);
            rows.forEach((row, index) => {
                const expected = entry.expected[index];
                const path = `site ${expected.referenceNumber}`;
                JS_ONLY_KEYS.forEach((key) => expect(row[key], `${path}.${key}`).toBeNull());
                const keys = Object.keys(row).filter((key) => !JS_ONLY_KEYS.includes(key)).sort();
                expect(keys, path).toEqual(Object.keys(expected).sort());
                const trimmed = Object.fromEntries(keys.map((key) => [key, row[key]]));
                expectGolden(trimmed, expected, path, { bases: { axes: expected.eigenvalues } });
            });
        });
    }
});
