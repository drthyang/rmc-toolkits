// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Cross-runtime parity of the Structure-page KDE: the browser worker
// (computeKde, CPU path) against the reference-grade Python engine
// (rmc_toolkits/kde.py -> oriented_kde_slice) via the golden written by
// tests/generate_kde_fixture.py. Every case keeps the slab below the 6000-point
// fit cap, so both runtimes sum the same rows with the same kernel; the maps must
// agree to 1e-6 of the peak, and declined slabs must be declined for the same
// reason. Real-data cases read the .rmc6f through the browser parser, exactly as
// the app does; the GaNb4Se8 run under data/ is gitignored, so those cases skip
// when it is absent (CI), while the committed demo run always runs.

import { existsSync, readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { structureFromRmc6f } from '../../browserData.js';
import { computeKde } from '../localKdeWorker';
import fixture from '../../__tests__/fixtures/kde_parity_fixture.json';

const REPO_ROOT = fileURLToPath(new URL('../../../../../', import.meta.url));
const CUBE_CORNERS = [
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
    [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]
];
const dot = (a, b) => a.reduce((sum, value, index) => sum + value * b[index], 0);
// StructurePage.projectionRange()
const projectionRange = (normal) => {
    const values = CUBE_CORNERS.map((corner) => dot(corner, normal));
    return [Math.min(...values), Math.max(...values)];
};

const structures = new Map();
const pointsFor = (testCase) => {
    if (testCase.points) {
        return testCase.points.map(([x, y, z]) => ({ x, y, z, element: 'X' }));
    }
    const path = `${REPO_ROOT}${testCase.rmc6f}`;
    if (!structures.has(path)) {
        const text = readFileSync(path, 'utf8');
        structures.set(path, structureFromRmc6f({ text, path }, 1_000_000).points);
    }
    const all = structures.get(path);
    return testCase.element === 'all' ? all : all.filter((point) => point.element === testCase.element);
};

const available = (testCase) => !testCase.rmc6f || existsSync(`${REPO_ROOT}${testCase.rmc6f}`);

describe('browser KDE vs Python reference (kde_parity_fixture.json)', () => {
    for (const testCase of fixture.cases) {
        const run = available(testCase) ? it : it.skip;
        run(testCase.name, async () => {
            const expected = testCase.expected;
            const result = await computeKde({
                points: pointsFor(testCase),
                normal: testCase.normal,
                uVector: testCase.u,
                vVector: testCase.v,
                range: projectionRange(testCase.normal),
                zCenter: testCase.z,
                thickness: testCase.dz,
                bandwidth: testCase.bw,
                gridSize: testCase.grid,
                logScale: false
            });

            expect(result.slabCount).toBe(expected.slabCount);
            expect(result.fitCount).toBe(expected.fitCount);
            if (!expected.kernel) {
                // Declined by both runtimes, for the same reason.
                expect(result.vmax).toBe(0);
                expect(result.message).toBe(expected.message);
                expect(result.messageCode).toBe(expected.messageCode);
                expect(result.kernel).toBeNull();
                return;
            }

            // Same density grid to 1e-6 of the peak.
            let worst = 0;
            expected.densityOverPeak.forEach((row, y) => row.forEach((value, x) => {
                worst = Math.max(worst, Math.abs(result.density[y][x] / expected.vmax - value));
            }));
            expect(worst).toBeLessThan(1e-6);
            expect(Math.abs(result.vmax / expected.vmax - 1)).toBeLessThan(1e-6);

            // Because it is the same kernel H = f^2 C, entry by entry, with the
            // same diagnostics.
            expect(result.message).toBeNull();
            expect(result.messageCode).toBeNull();
            expect(result.warnings.map((warning) => warning.code)).toEqual(expected.warnings);
            const scale = expected.kernel.sigmaMajor ** 2;
            expected.kernel.covariance.flat().forEach((value, index) => {
                expect(Math.abs(result.kernel.covariance.flat()[index] - value) / scale).toBeLessThan(1e-9);
            });
            expect(result.kernel.sigmaMinor / expected.kernel.sigmaMinor).toBeCloseTo(1, 7);
            expect(result.kernel.sigmaMajor / expected.kernel.sigmaMajor).toBeCloseTo(1, 7);
        });
    }
});
