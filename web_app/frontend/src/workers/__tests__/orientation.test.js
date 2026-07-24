// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, it, expect } from 'vitest';

import {
    MIN_FREQUENCY,
    MAX_FREQUENCY,
    assignCells,
    goldbergTiling,
    orientationHistogram,
    recommendedFrequency,
    siteOrientationHistogram
} from '../orientation.js';
import { siteDisplacementsFromRmc6f } from '../pcaKde.js';
import { handlePcaMessage } from '../pcaKdeWorker.js';

// Deterministic Gaussian sampler so tests are reproducible without a dependency.
const makeRng = (seed) => {
    let value = seed >>> 0;
    const next = () => {
        value += 0x6d2b79f5;
        let mixed = value;
        mixed = Math.imul(mixed ^ (mixed >>> 15), mixed | 1);
        mixed ^= mixed + Math.imul(mixed ^ (mixed >>> 7), mixed | 61);
        return ((mixed ^ (mixed >>> 14)) >>> 0) / 4294967296;
    };
    return () => {
        let u = 0;
        let v = 0;
        while (u === 0) u = next();
        while (v === 0) v = next();
        return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
    };
};

const gaussianCloud = (n, seed, sigma) => {
    const gauss = makeRng(seed);
    const points = [];
    for (let i = 0; i < n; i += 1) {
        points.push([gauss() * sigma[0], gauss() * sigma[1], gauss() * sigma[2]]);
    }
    return points;
};

const norm = (u) => Math.hypot(u[0], u[1], u[2]);

describe('goldbergTiling', () => {
    it('has 10*nu^2 + 2 cells with exactly 12 pentagons', () => {
        for (const frequency of [1, 2, 3, 5, 8]) {
            const tiling = goldbergTiling(frequency);
            expect(tiling.cellCount).toBe(10 * frequency * frequency + 2);
            expect(tiling.sizes.filter((size) => size === 5)).toHaveLength(12);
            expect(tiling.sizes.every((size) => size === 5 || size === 6)).toBe(true);
        }
    });

    it('has unit cell centers and areas summing to 4*pi', () => {
        const tiling = goldbergTiling(6);
        tiling.centers.forEach((center) => expect(norm(center)).toBeCloseTo(1, 12));
        const total = tiling.areas.reduce((sum, area) => sum + area, 0);
        expect(total).toBeCloseTo(4 * Math.PI, 9);
        expect(Math.min(...tiling.areas)).toBeGreaterThan(0);
    });

    it('keeps cell areas nearly equal (no polar slivers)', () => {
        const tiling = goldbergTiling(8);
        expect(Math.max(...tiling.areas) / Math.min(...tiling.areas)).toBeLessThan(2);
    });

    it('has symmetric neighbor lists', () => {
        const tiling = goldbergTiling(4);
        tiling.neighbors.forEach((row, cell) => {
            row.forEach((neighbor) => {
                if (neighbor < 0) return;
                expect(tiling.neighbors[neighbor]).toContain(cell);
            });
        });
    });

    it('has an exact, involutive antipode map', () => {
        const tiling = goldbergTiling(5);
        tiling.antipode.forEach((other, cell) => {
            expect(tiling.centers[other][0]).toBeCloseTo(-tiling.centers[cell][0], 12);
            expect(tiling.centers[other][1]).toBeCloseTo(-tiling.centers[cell][1], 12);
            expect(tiling.centers[other][2]).toBeCloseTo(-tiling.centers[cell][2], 12);
            expect(tiling.antipode[other]).toBe(cell);
        });
    });

    it('rejects out-of-range frequencies', () => {
        expect(() => goldbergTiling(MIN_FREQUENCY - 1)).toThrow();
        expect(() => goldbergTiling(MAX_FREQUENCY + 1)).toThrow();
    });
});

describe('assignCells', () => {
    it('matches brute-force nearest center', () => {
        const gauss = makeRng(3);
        const directions = [];
        for (let i = 0; i < 2000; i += 1) {
            const u = [gauss(), gauss(), gauss()];
            const length = norm(u);
            directions.push([u[0] / length, u[1] / length, u[2] / length]);
        }
        for (const frequency of [2, 6, 11]) {
            const tiling = goldbergTiling(frequency);
            const assigned = assignCells(tiling, directions);
            let agree = 0;
            directions.forEach((direction, index) => {
                let best = 0;
                let bestDot = -Infinity;
                tiling.centers.forEach((center, cell) => {
                    const value = center[0] * direction[0] + center[1] * direction[1] + center[2] * direction[2];
                    if (value > bestDot) { bestDot = value; best = cell; }
                });
                if (assigned[index] === best) agree += 1;
            });
            expect(agree / directions.length).toBeGreaterThan(0.9999);
        }
    });

    it('assigns each cell center to itself', () => {
        const tiling = goldbergTiling(7);
        const assigned = assignCells(tiling, tiling.centers);
        assigned.forEach((cell, index) => expect(cell).toBe(index));
    });
});

describe('recommendedFrequency', () => {
    it('scales with point count and clamps at the floor', () => {
        expect(recommendedFrequency(100)).toBeLessThanOrEqual(recommendedFrequency(10000));
        expect(recommendedFrequency(0)).toBe(MIN_FREQUENCY);
    });
});

describe('orientationHistogram', () => {
    it('integrates the density to one over the sphere', () => {
        const cloud = gaussianCloud(20000, 1, [0.1, 0.1, 0.1]);
        const result = orientationHistogram(cloud, { frequency: 8 });
        const integral = result.density.reduce((sum, value, cell) => sum + value * result.areas[cell], 0);
        expect(integral).toBeCloseTo(1, 9);
    });

    it('is flat for an isotropic cloud', () => {
        const cloud = gaussianCloud(20000, 2, [0.1, 0.1, 0.1]);
        const result = orientationHistogram(cloud, { frequency: 8 });
        const mean = result.enhancement.reduce((sum, value) => sum + value, 0) / result.cellCount;
        expect(mean).toBeCloseTo(1, 1);
        expect(Math.max(...result.enhancement.map((value) => Math.abs(value - 1)))).toBeLessThan(1);
        expect(result.significance).toBeLessThan(1.5);
        expect(Math.abs(result.orientationAnisotropy)).toBeLessThan(0.05);
    });

    it('peaks on the long axis of a lobed cloud', () => {
        const cloud = gaussianCloud(20000, 3, [0.5, 0.1, 0.1]);
        const result = orientationHistogram(cloud, { frequency: 8 });
        expect(Math.abs(result.peakDirection[0])).toBeGreaterThan(0.9);
        expect(result.peakEnhancement).toBeGreaterThan(2);
        expect(result.significance).toBeGreaterThan(3);
        expect(result.orientationAnisotropy).toBeGreaterThan(0.5);
    });

    it('sees antipodal asymmetry in a one-sided cloud', () => {
        const gauss = makeRng(5);
        const cloud = [];
        for (let i = 0; i < 8000; i += 1) {
            cloud.push([Math.abs(gauss() * 0.05) + 0.05, gauss() * 0.05, gauss() * 0.05]);
        }
        const result = orientationHistogram(cloud, { frequency: 6 });
        const tiling = goldbergTiling(6);
        const peak = result.peakCell;
        expect(result.enhancement[peak]).toBeGreaterThan(4 * result.enhancement[tiling.antipode[peak]]);
        expect(result.antipodalAsymmetry).toBeGreaterThan(0.5);
        expect(result.antipodalAsymmetry).toBeGreaterThan(3 * result.antipodalAsymmetryNull);
    });

    it('keeps a symmetric cloud near the asymmetry noise floor', () => {
        const result = orientationHistogram(gaussianCloud(20000, 13, [0.5, 0.1, 0.1]), { frequency: 6 });
        expect(result.antipodalAsymmetry).toBeLessThan(2 * result.antipodalAsymmetryNull);
        expect(result.antipodalAsymmetry).toBeGreaterThan(0.2 * result.antipodalAsymmetryNull);
    });

    it('conserves mass under smoothing and keeps the lobe in place', () => {
        const cloud = gaussianCloud(5000, 6, [0.5, 0.1, 0.1]);
        const raw = orientationHistogram(cloud, { frequency: 8, smoothing: 0 });
        const smoothed = orientationHistogram(cloud, { frequency: 8, smoothing: 3 });
        const total = (values) => values.reduce((sum, value) => sum + value, 0);
        expect(total(smoothed.mass)).toBeCloseTo(total(raw.mass), 9);
        expect(smoothed.peakEnhancement).toBeLessThan(raw.peakEnhancement);
        expect(Math.abs(smoothed.peakDirection[0])).toBeGreaterThan(0.85);
    });

    it('changes the winner under amplitude^2 weighting', () => {
        const gauss = makeRng(9);
        const cloud = [];
        for (let i = 0; i < 9000; i += 1) {
            cloud.push([gauss() * 0.08, gauss() * 0.01, gauss() * 0.01]);
        }
        for (let i = 0; i < 1000; i += 1) {
            cloud.push([gauss() * 0.02, gauss() * 0.02, gauss() * 0.8]);
        }
        const byCount = orientationHistogram(cloud, { frequency: 6, weight: 'count' });
        const byAmp2 = orientationHistogram(cloud, { frequency: 6, weight: 'amplitude2' });
        expect(Math.abs(byCount.peakDirection[0])).toBeGreaterThan(0.9);
        expect(Math.abs(byAmp2.peakDirection[2])).toBeGreaterThan(0.9);
    });

    it('applies amplitude cutoffs', () => {
        const cloud = gaussianCloud(1000, 7, [0.1, 0.1, 0.1]);
        const result = orientationHistogram(cloud, { frequency: 3, minAmplitudeQuantile: 0.25 });
        expect(result.totalPoints).toBe(1000);
        expect(result.usedPoints).toBeGreaterThan(745);
        expect(result.usedPoints).toBeLessThan(755);
        expect(result.usedPoints + result.rejectedPoints).toBe(1000);
    });

    it('aligns an oblique lobe with PC1 in the pca frame', () => {
        const base = gaussianCloud(15000, 8, [0.5, 0.08, 0.08]);
        const angle = 0.7;
        const cloud = base.map(([x, y, z]) => [
            x * Math.cos(angle) - y * Math.sin(angle),
            x * Math.sin(angle) + y * Math.cos(angle),
            z
        ]);
        const result = orientationHistogram(cloud, { frequency: 8, frame: 'pca' });
        expect(Math.abs(result.peakDirection[0])).toBeGreaterThan(0.9);
    });

    it('matches the direct orientation tensor', () => {
        const cloud = gaussianCloud(3000, 10, [0.5, 0.1, 0.1]);
        const result = orientationHistogram(cloud, { frequency: 4 });
        const tensor = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        cloud.forEach((row) => {
            const length = norm(row);
            const u = [row[0] / length, row[1] / length, row[2] / length];
            for (let a = 0; a < 3; a += 1) {
                for (let b = 0; b < 3; b += 1) tensor[a][b] += u[a] * u[b] / cloud.length;
            }
        });
        for (let a = 0; a < 3; a += 1) {
            for (let b = 0; b < 3; b += 1) {
                expect(result.orientationTensor[a][b]).toBeCloseTo(tensor[a][b], 9);
            }
        }
        const trace = result.orientationEigenvalues.reduce((sum, value) => sum + value, 0);
        expect(trace).toBeCloseTo(1, 9);
    });

    it('honors the geometry toggle and pentagon vertex counts', () => {
        const cloud = gaussianCloud(500, 11, [0.1, 0.1, 0.1]);
        const withGeometry = orientationHistogram(cloud, { frequency: 3, geometry: true });
        const without = orientationHistogram(cloud, { frequency: 3, geometry: false });
        expect(withGeometry.polygons).toBeDefined();
        expect(without.polygons).toBeUndefined();
        withGeometry.polygons.forEach((polygon, cell) => {
            expect(polygon).toHaveLength(withGeometry.sizes[cell]);
        });
    });

    it('rejects bad input', () => {
        expect(() => orientationHistogram([[0, 0]])).toThrow();
        expect(() => orientationHistogram([[0, 0, 0], [0, 0, 0]])).toThrow();
        expect(() => orientationHistogram([[1, 1, 1]], { weight: 'nope' })).toThrow();
        expect(() => orientationHistogram([[1, 1, 1]], { frame: 'lab' })).toThrow();
        expect(() => orientationHistogram([[1, 1, 1]], { minAmplitudeQuantile: 1 })).toThrow();
    });
});

// --- site-level extraction + worker routing -----------------------------------

const makeCloudRmc6f = (sigmaFrac, { supercell = [6, 6, 6], seed = 1, element = 'Ga' } = {}) => {
    const gauss = makeRng(seed);
    const lattice = [
        [8 * supercell[0], 0, 0],
        [0, 8 * supercell[1], 0],
        [0, 0, 8 * supercell[2]]
    ];
    const lines = [
        `Supercell dimensions ${supercell[0]} ${supercell[1]} ${supercell[2]}`,
        'Lattice vectors (Ang):',
        lattice[0].join(' '),
        lattice[1].join(' '),
        lattice[2].join(' '),
        'Atoms:'
    ];
    let atom = 0;
    for (let ix = 0; ix < supercell[0]; ix += 1) {
        for (let iy = 0; iy < supercell[1]; iy += 1) {
            for (let iz = 0; iz < supercell[2]; iz += 1) {
                atom += 1;
                const coord = [
                    ix / supercell[0] + gauss() * sigmaFrac[0],
                    iy / supercell[1] + gauss() * sigmaFrac[1],
                    iz / supercell[2] + gauss() * sigmaFrac[2]
                ];
                lines.push(
                    `${atom} ${element} [1] ${coord[0].toFixed(10)} ${coord[1].toFixed(10)} ${coord[2].toFixed(10)} 1 ${ix} ${iy} ${iz}`
                );
            }
        }
    }
    return lines.join('\n');
};

describe('siteOrientationHistogram', () => {
    it('tags reference number and element and finds the drawn lobe', () => {
        const text = makeCloudRmc6f([0.03, 0.01, 0.01], { supercell: [8, 8, 8], seed: 2, element: 'Nb' });
        const parsed = siteDisplacementsFromRmc6f(text);
        const result = siteOrientationHistogram(parsed, {
            referenceNumber: 1,
            frequency: 4,
            geometry: false
        });
        expect(result.referenceNumber).toBe(1);
        expect(result.element).toBe('Nb');
        expect(result.totalPoints).toBe(512);
        expect(Math.abs(result.peakDirection[0])).toBeGreaterThan(0.8);
    });
});

describe('worker routing', () => {
    it("answers kind: 'orientation' with the histogram shape", async () => {
        const text = makeCloudRmc6f([0.03, 0.01, 0.01], { supercell: [6, 6, 6], seed: 4 });
        const result = await handlePcaMessage(
            { kind: 'orientation', referenceNumber: 1, frequency: 3, geometry: false },
            async () => text
        );
        expect(result.cellCount).toBe(92);
        expect(result.referenceNumber).toBe(1);
        expect(result.enhancement).toHaveLength(92);
        expect(result.browserOrientation).toBe(true);
    });
});
