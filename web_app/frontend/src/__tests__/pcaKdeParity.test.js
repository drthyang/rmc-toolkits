// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The browser PCA-KDE (workers/pcaKde.js) above the fit cap against the Python
// reference (rmc_toolkits/pca_kde.py) through the golden
// tests/generate_pca_kde_fixture.py writes; tests/test_pca_kde_fixture.py checks
// the committed file against the Python engine.
//
// Above MAX_PCA_FIT_POINTS both engines fit the KDE on the rows one draw picks
// (fitSubsampleIndices ⟷ fit_subsample_indices: mulberry32 + a partial
// Fisher-Yates), so the draw is held exactly -- the uint32 stream, every drawn
// index, its order -- and the volume that follows to the goldenCompare.js rule:
// integers exactly, other numbers to 1e-9 relative + 1e-12 absolute, the axes
// eigenspace by eigenspace.

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';

import {
    MAX_PCA_FIT_POINTS,
    displacementCloud,
    fitSubsampleIndices,
    mulberry32,
    pcaKdeVolume,
    siteDisplacementsFromRmc6f,
    sitePcaKde
} from '../workers/pcaKde.js';
import fixture from './fixtures/pca_kde_fixture.json';
import { expectGolden } from './goldenCompare.js';

// The browser payload's transport flag; Python has no such key.
const JS_ONLY_KEYS = ['browserPcaKde'];

const indexChecksum = (indices) => indices.reduce((sum, index, i) => sum + (i + 1) * index, 0);

const cloudChecksum = (cloud) => cloud.reduce(
    (sum, [x, y, z], i) => sum + ((i % 1000) + 1) * (x + 2 * y + 3 * z), 0
);

// The synthetic site of generate_pca_kde_fixture.synthetic_cloud: the same
// operations in the same order, so the same doubles.
const syntheticCloud = ({ seed, count, transform: m }) => {
    const random = mulberry32(seed);
    const points = [];
    for (let i = 0; i < count; i += 1) {
        const e = [0, 1, 2].map(() => {
            const u0 = random();
            const u1 = random();
            const u2 = random();
            return u0 + u1 + u2 - 1.5;
        });
        points.push([0, 1, 2].map((j) => e[0] * m[0][j] + e[1] * m[1][j] + e[2] * m[2][j]));
    }
    return points;
};

const demoText = (name) => readFileSync(new URL(`../../public/demo/${name}`, import.meta.url), 'utf8');

describe('PCA-KDE fit subsample parity (pca_kde_fixture.json)', () => {
    it('shares the fit cap', () => {
        expect(MAX_PCA_FIT_POINTS).toBe(fixture.maxFitPoints);
    });

    it('draws the same mulberry32 stream, seed for seed', () => {
        fixture.stream.forEach(({ seed, uint32 }) => {
            const random = mulberry32(seed);
            expect(uint32.map(() => random() * 4294967296), `seed ${seed}`).toEqual(uint32);
        });
    });

    it('draws the same indices in the same order', () => {
        fixture.draws.forEach((draw) => {
            const label = `${draw.count} -> ${draw.limit}, seed ${draw.seed}`;
            const indices = fitSubsampleIndices(draw.count, draw.limit, draw.seed);
            expect(indices.length, label).toBe(draw.size);
            expect(indexChecksum(indices), label).toBe(draw.checksum);
            if (draw.indices) expect(indices, label).toEqual(draw.indices);
            if (draw.excluded) {
                const chosen = new Set(indices);
                const excluded = [];
                for (let index = 0; index < draw.count; index += 1) if (!chosen.has(index)) excluded.push(index);
                expect(excluded, label).toEqual(draw.excluded);
                for (let i = 1; i < indices.length; i += 1) expect(indices[i] > indices[i - 1], label).toBe(true);
            }
        });
    });

    for (const entry of fixture.volumes) {
        it(entry.name, () => {
            const { source, options, cloud: want } = entry;
            let cloud;
            let result;
            if (source.synthetic) {
                cloud = syntheticCloud(fixture.synthetic);
                result = pcaKdeVolume(cloud, options);
            } else {
                const parsed = siteDisplacementsFromRmc6f(demoText(source.demo));
                const selection = { referenceNumber: source.referenceNumber ?? null, element: source.element ?? null };
                cloud = displacementCloud(parsed, selection).cloud;
                result = sitePcaKde(parsed, { ...selection, ...options });
            }
            // The same rows in the same order: the draw picks rows by position.
            expect(cloud.length).toBe(want.count);
            // The synthetic points are built with exactly rounded operations only.
            if (source.synthetic) expect([...cloud.slice(0, 3), ...cloud.slice(-3)]).toEqual([...want.head, ...want.tail]);
            expectGolden(cloud.slice(0, 3), want.head, 'cloud head');
            expectGolden(cloud.slice(-3), want.tail, 'cloud tail');
            expectGolden(cloudChecksum(cloud), want.checksum, 'cloud checksum');
            expect(indexChecksum(fitSubsampleIndices(cloud.length, MAX_PCA_FIT_POINTS, 0))).toBe(entry.fitChecksum);

            const expected = entry.expected;
            expect(result.count).toBe(want.count);
            expect(result.fitCount).toBe(MAX_PCA_FIT_POINTS);
            JS_ONLY_KEYS.forEach((key) => expect(result[key]).toBe(true));
            const keys = Object.keys(result).filter((key) => !JS_ONLY_KEYS.includes(key)).sort();
            expect(keys).toEqual(Object.keys(expected).sort());
            // The walls' `axes` are index pairs, not an eigenbasis: compared apart.
            const { projections: wantWalls, ...wantRest } = expected;
            const trimmed = Object.fromEntries(keys.filter((key) => key !== 'projections').map((key) => [key, result[key]]));
            expectGolden(trimmed, wantRest, entry.name, { bases: { axes: expected.eigenvalues } });
            expectGolden(result.projections, wantWalls, `${entry.name} projections`);
        });
    }
});
