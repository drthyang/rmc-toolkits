// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Regression tests for the browser PCA-ellipsoid engine (workers/pcaKde.js),
// one describe block per defect found by the 1.0 audit. They mirror
// tests/test_pca_regressions.py so both engines are pinned to the same rule.

import { existsSync, readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, it, expect } from 'vitest';

import {
    displacementCloud,
    eigenDecomposition,
    pcaKdeVolume,
    siteDisplacementsFromRmc6f,
    siteEllipsoids,
    sitePcaKde
} from '../pcaKde.js';
import { handlePcaMessage } from '../pcaKdeWorker.js';

const AVERAGE_RMC6F = fileURLToPath(new URL('../../../../../data/5K_try1/GaNb4Se8_5KAVERAGE.rmc6f', import.meta.url));

// Deterministic standard-normal sampler (mulberry32 + Box-Muller).
const makeGauss = (seed) => {
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

const header = (supercell, cellEdge) => [
    `Supercell dimensions ${supercell[0]} ${supercell[1]} ${supercell[2]}`,
    'Lattice vectors (Ang):',
    `${supercell[0] * cellEdge} 0 0`,
    `0 ${supercell[1] * cellEdge} 0`,
    `0 0 ${supercell[2] * cellEdge}`,
    'Atoms:'
];

// One site at unit-cell fraction `site` with isotropic Gaussian spread; every box
// copy gets one atom whose supercell fraction is wrapped into [0, 1), exactly as
// RMCProfile writes it.
const wrappedSiteLines = (site, sigma, { supercell, cellEdge, seed, element = 'Se', reference = 1, start = 1 }) => {
    const gauss = makeGauss(seed);
    const lines = [];
    let atom = start;
    for (let ix = 0; ix < supercell[0]; ix += 1) {
        for (let iy = 0; iy < supercell[1]; iy += 1) {
            for (let iz = 0; iz < supercell[2]; iz += 1) {
                const cell = [ix, iy, iz];
                const coord = [0, 1, 2].map((i) => {
                    const f = (cell[i] + site[i] + (gauss() * sigma) / cellEdge) / supercell[i];
                    return ((f % 1) + 1) % 1;
                });
                lines.push(`${atom} ${element} [1] ${coord.map((v) => v.toFixed(12)).join(' ')} `
                    + `${reference} ${ix} ${iy} ${iz}`);
                atom += 1;
            }
        }
    }
    return lines;
};

describe('site-centred supercell fold (pca.parity.1 / physics.9 / numerics.13 / parity.25 / numerics.31 / physics.38)', () => {
    const intact = (supercell, site, sigma = 0.08) => {
        const cellEdge = 10;
        const text = [...header(supercell, cellEdge),
            ...wrappedSiteLines(site, sigma, { supercell, cellEdge, seed: 3 })].join('\n');
        const parsed = siteDisplacementsFromRmc6f(text);
        const [entry] = siteEllipsoids(parsed.sites);
        entry.rms.forEach((value) => {
            expect(value).toBeGreaterThan(0.7 * sigma);
            expect(value).toBeLessThan(1.3 * sigma);
        });
        const maxDisp = Math.max(...parsed.sites[0].displacements.flat().map(Math.abs));
        expect(maxDisp).toBeLessThan(8 * sigma);
        entry.siteFractional.forEach((value, i) => {
            const delta = value - site[i];
            expect(Math.abs(delta - Math.round(delta))).toBeLessThan(0.01);
        });
    };

    it('keeps a site at x = 1/2 intact along a single-cell axis', () => {
        intact([6, 6, 1], [0.3, 0.3, 0.5]);
    });

    it('keeps a site near x = 1 intact along a two-cell axis', () => {
        intact([6, 6, 2], [0.3, 0.3, 0.99]);
    });

    it('keeps a large-spread site at x = 1/2 intact in a 1xNxN box', () => {
        intact([1, 8, 8], [0.5, 0.5, 0.5], 0.2);
    });

    it('leaves regular boxes unchanged', () => {
        intact([6, 6, 6], [0.5, 0.99, 0.0]);
    });
});

describe('zero-spread sites (pca.parity.3 / parity.24 / parity.30)', () => {
    const frozenText = () => {
        const supercell = [6, 6, 6];
        const cellEdge = 10;
        const frozen = wrappedSiteLines([0.25, 0.5, 0.75], 0, { supercell, cellEdge, seed: 1, element: 'Ga', reference: 1 });
        const moving = wrappedSiteLines([0.6, 0.1, 0.3], 0.07, {
            supercell, cellEdge, seed: 2, element: 'Se', reference: 2, start: frozen.length + 1
        });
        return [...header(supercell, cellEdge), ...frozen, ...moving].join('\n');
    };

    it('flags a frozen site and suppresses its round-off anisotropy, kurtosis and axes', () => {
        const parsed = siteDisplacementsFromRmc6f(frozenText());
        const [frozen, moving] = siteEllipsoids(parsed.sites);
        expect(frozen.zeroSpread).toBe(true);
        expect(frozen.degenerate).toBe(true);
        expect(frozen.anisotropy).toBeNull();
        expect(frozen.nonGaussianity).toBeNull();
        expect(frozen.excessKurtosis).toEqual([null, null, null]);
        expect(frozen.axes).toBeNull();
        expect(frozen.uIso).toBeLessThan(1e-12);
        expect(moving.zeroSpread).toBe(false);
        expect(moving.degenerate).toBe(false);
        moving.excessKurtosis.forEach((value) => expect(Number.isFinite(value)).toBe(true));
    });

    it('refuses a KDE of a zero-spread cloud', () => {
        const parsed = siteDisplacementsFromRmc6f(frozenText());
        expect(() => sitePcaKde(parsed, { referenceNumber: 1, grid: 12, projections: false })).toThrow(/zero spread/);
        const femto = Array.from({ length: 200 }, (_, i) => [1e-14 * Math.sin(i), 1e-14 * Math.cos(3 * i), 1e-14 * Math.sin(7 * i)]);
        expect(() => pcaKdeVolume(femto, { grid: 8 })).toThrow(/zero spread/);
    });

    it('reports no kurtosis along a collapsed axis', () => {
        const gauss = makeGauss(5);
        const cloud = Array.from({ length: 2000 }, () => [0.1 * gauss(), 0.07 * gauss(), 0]);
        const result = pcaKdeVolume(cloud, { grid: 16, projections: false });
        expect(result.degenerate).toBe(true);
        expect(result.excessKurtosis[2]).toBeNull();
        expect(Number.isFinite(result.excessKurtosis[0])).toBe(true);
        expect(Number.isFinite(result.nonGaussianity)).toBe(true);
    });

    it('Jacobi diagonalises a matrix of any scale (relative stopping test)', () => {
        const theta = 0.6;
        const rot = [[Math.cos(theta), -Math.sin(theta), 0], [Math.sin(theta), Math.cos(theta), 0], [0, 0, 1]];
        const diag = [4, 2, 1];
        for (const scale of [1, 1e-4, 1e-28]) {
            const M = [0, 1, 2].map((i) => [0, 1, 2].map((j) => scale * [0, 1, 2]
                .reduce((sum, k) => sum + rot[i][k] * diag[k] * rot[j][k], 0)));
            const { eigenvalues, axes } = eigenDecomposition(M);
            expect(eigenvalues[0] / scale).toBeCloseTo(4, 10);
            expect(eigenvalues[1] / scale).toBeCloseTo(2, 10);
            expect(eigenvalues[2] / scale).toBeCloseTo(1, 10);
            // PC1 is the first column of the rotation, up to sign.
            expect(Math.abs(axes[0][0] * rot[0][0] + axes[0][1] * rot[1][0])).toBeCloseTo(1, 10);
        }
    });

    it.skipIf(!existsSync(AVERAGE_RMC6F))('flags every site of the real GaNb4Se8 AVERAGE configuration', () => {
        const parsed = siteDisplacementsFromRmc6f(readFileSync(AVERAGE_RMC6F, 'utf8'));
        const entries = siteEllipsoids(parsed.sites);
        expect(entries).toHaveLength(52);
        entries.forEach((entry) => {
            expect(entry.zeroSpread).toBe(true);
            expect(entry.degenerate).toBe(true);
            expect(entry.nonGaussianity).toBeNull();
        });
    });
});

describe('rotation-invariant non-Gaussianity and axis resolution (pca.physics.8 / physics.40)', () => {
    // Spherical scale mixture of normals (elliptical): scale 1 (80%) or 2 (20%);
    // its excess kurtosis is 3 E[s^4] / E[s^2]^2 - 3 = 1.6875 along every direction.
    const scaleMixture = (n, seed, scale = [1, 1, 1]) => {
        const gauss = makeGauss(seed);
        let state = (seed * 2654435761) >>> 0;
        const uniform = () => {
            state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
            return state / 4294967296;
        };
        return Array.from({ length: n }, () => {
            const s = uniform() < 0.8 ? 1 : 2;
            return [gauss() * s * scale[0], gauss() * s * scale[1], gauss() * s * scale[2]];
        });
    };
    const transform = (cloud, m) => cloud.map((p) => [0, 1, 2].map((j) => p[0] * m[0][j] + p[1] * m[1][j] + p[2] * m[2][j]));

    it('is affine invariant, so a noise-oriented frame cannot move it', () => {
        const cloud = scaleMixture(3000, 1).map((p) => p.map((v) => 0.05 * v));
        const base = pcaKdeVolume(cloud, { grid: 8, projections: false }).nonGaussianity;
        const matrices = [
            [[1.0, 0.4, -0.3], [0.0, 0.8, 0.5], [0.2, 0.0, 0.6]],
            [[1, 0, 0], [0, 1 + 1e-7, 0], [0, 0, 1 + 2e-7]],
            [[1 + 2e-7, 0, 0], [0, 1, 0], [0, 0, 1 + 1e-7]]
        ];
        matrices.forEach((m) => {
            const moved = pcaKdeVolume(transform(cloud, m), { grid: 8, projections: false }).nonGaussianity;
            expect(Math.abs(moved - base)).toBeLessThan(1e-9);
        });
    });

    it('equals the marginal excess kurtosis of an elliptical cloud', () => {
        const cloud = scaleMixture(40000, 2, [0.12, 0.08, 0.05]);
        const result = pcaKdeVolume(cloud, { grid: 8, projections: false, maxFitPoints: 40000 });
        expect(Math.abs(result.nonGaussianity - 1.6875)).toBeLessThan(0.12);
    });

    it('flags axes inside a degenerate eigenvalue pair as unresolved', () => {
        const gauss = makeGauss(4);
        const cloudOf = (sigma) => Array.from({ length: 1000 }, () => sigma.map((s) => s * gauss()));
        expect(pcaKdeVolume(cloudOf([0.1, 0.1, 0.1]), { grid: 8, projections: false }).axisResolved)
            .toEqual([false, false, false]);
        expect(pcaKdeVolume(cloudOf([0.2, 0.1, 0.1]), { grid: 8, projections: false }).axisResolved)
            .toEqual([true, false, false]);
        expect(pcaKdeVolume(cloudOf([0.2, 0.1, 0.05]), { grid: 8, projections: false }).axisResolved)
            .toEqual([true, true, true]);
    });

    it('site table and volume report the same statistics', () => {
        const supercell = [10, 10, 10];
        const cellEdge = 8;
        const text = [...header(supercell, cellEdge),
            ...wrappedSiteLines([0.25, 0.25, 0.25], 0.08, { supercell, cellEdge, seed: 9 })].join('\n');
        const parsed = siteDisplacementsFromRmc6f(text);
        const [entry] = siteEllipsoids(parsed.sites);
        const volume = pcaKdeVolume(parsed.sites[0].displacements, { grid: 8, projections: false });
        expect(entry.axisResolved).toEqual([false, false, false]);
        expect(Math.abs(entry.nonGaussianity - volume.nonGaussianity)).toBeLessThan(1e-9);
    });

    it('reads a symmetric split site as platykurtic', () => {
        const gauss = makeGauss(6);
        const [d, s] = [0.15, 0.08];
        const cloud = Array.from({ length: 8000 }, (_, i) => [(i % 2 ? d : -d) + s * gauss(), s * gauss(), s * gauss()]);
        const result = pcaKdeVolume(cloud, { grid: 8, projections: false });
        const analytic = (-2 * d ** 4) / (s * s + d * d) ** 2;
        expect(result.axisResolved[0]).toBe(true);
        expect(Math.abs(result.excessKurtosis[0] - analytic)).toBeLessThan(0.08);
        expect(result.nonGaussianity).toBeLessThan(0);
    });
});

describe('mixed-occupancy sites (pca.parity.4 / numerics.17 / parity.20 / parity.26 / numerics.33 / physics.39)', () => {
    // Reference 1 is mixed (majority first in the file, minority LAST, as RMCProfile
    // groups atoms by type); reference 2 is one species written with an upper-case token.
    const mixedText = ({ minorityEvery = 4, majority = 'Ga', minority = 'In', otherToken = 'SE' } = {}) => {
        const supercell = [4, 4, 4];
        const cellEdge = 8;
        const lines = wrappedSiteLines([0.25, 0.25, 0.25], 0.06, { supercell, cellEdge, seed: 11, element: majority, reference: 1 });
        const major = lines.filter((_, k) => k % minorityEvery);
        const minor = lines.filter((_, k) => !(k % minorityEvery)).map((line) => line.replace(` ${majority} `, ` ${minority} `));
        const other = wrappedSiteLines([0.6, 0.1, 0.3], 0.07, { supercell, cellEdge, seed: 12, element: otherToken, reference: 2 });
        return { text: [...header(supercell, cellEdge), ...major, ...other, ...minor].join('\n'), nMajor: major.length, nMinor: minor.length };
    };

    it('labels a mixed site by its majority species and reports its composition', async () => {
        const { text, nMajor, nMinor } = mixedText();
        const parsed = siteDisplacementsFromRmc6f(text);
        const [mixed, pure] = siteEllipsoids(parsed.sites);
        expect(mixed.element).toBe('Ga');
        expect(mixed.mixed).toBe(true);
        expect(mixed.elementCounts).toEqual({ Ga: nMajor, In: nMinor });
        expect(pure.element).toBe('Se');
        expect(pure.mixed).toBe(false);
        expect(pure.elementCounts).toEqual({ Se: 64 });
        const summary = await handlePcaMessage({ kind: 'sites' }, async () => text);
        expect(summary.elements).toEqual(['Ga', 'In', 'Se']);
    });

    it('pools atoms by their own element', () => {
        const { text, nMajor, nMinor } = mixedText();
        const parsed = siteDisplacementsFromRmc6f(text);
        expect(displacementCloud(parsed, { element: 'In' }).cloud).toHaveLength(nMinor);
        expect(displacementCloud(parsed, { element: 'ga' }).cloud).toHaveLength(nMajor);
        expect(displacementCloud(parsed, { element: 'Se' }).cloud).toHaveLength(64);
        expect(() => displacementCloud(parsed, { element: 'Nb' })).toThrow(/Unknown element/);
        expect(sitePcaKde(parsed, { element: 'In', grid: 8, projections: false }).count).toBe(nMinor);
    });

    it('breaks a tie toward the alphabetically first species', () => {
        const { text } = mixedText({ minorityEvery: 2, majority: 'Fe', minority: 'Co' });
        const [entry] = siteEllipsoids(siteDisplacementsFromRmc6f(text).sites);
        expect(entry.elementCounts).toEqual({ Co: 32, Fe: 32 });
        expect(entry.element).toBe('Co');
    });
});

describe('non-finite inputs (pca.parity.6 / parity.23 / parity.27 / numerics.37)', () => {
    const gauss = makeGauss(3);
    const cloud = () => Array.from({ length: 300 }, () => [0.1 * gauss(), 0.1 * gauss(), 0.1 * gauss()]);

    it('rejects non-finite points with a clear message', () => {
        for (const bad of [Number.NaN, Infinity]) {
            const points = cloud();
            points[17][1] = bad;
            expect(() => pcaKdeVolume(points, { grid: 8 })).toThrow(/non-finite/);
        }
    });

    it('rejects non-finite parameters', () => {
        const points = cloud();
        const cases = [{ extent: Number.NaN }, { extent: Infinity }, { bwScale: Number.NaN },
            { bwScale: Infinity }, { bw: Number.NaN }, { bw: Infinity }, { grid: Number.NaN }];
        cases.forEach((options) => {
            expect(() => pcaKdeVolume(points, { grid: 8, ...options })).toThrow();
        });
    });
});

describe('cubic box is a display box only (pca.numerics.14 / numerics.32 / physics.44)', () => {
    it('keeps unit mass for a planar cloud on even and odd grids', () => {
        const gauss = makeGauss(8);
        const cloud = Array.from({ length: 1000 }, () => [0.1 * gauss(), 0.07 * gauss(), 0]);
        for (const grid of [40, 41]) {
            const result = pcaKdeVolume(cloud, { grid, extent: 4, cubicBox: true, projections: false });
            expect(result.degenerate).toBe(true);
            expect(result.mass).toBeGreaterThan(0.99);
            expect(result.mass).toBeLessThan(1.01);
            expect(result.massLevels).toHaveLength(101);
        }
    });

    it('samples the same volume whatever the display box', () => {
        const gauss = makeGauss(9);
        const cloud = Array.from({ length: 1000 }, () => [0.3 * gauss(), 0.1 * gauss(), 0.015 * gauss()]);
        const cubic = pcaKdeVolume(cloud, { grid: 40, extent: 4, cubicBox: true, projections: false });
        const plain = pcaKdeVolume(cloud, { grid: 40, extent: 4, cubicBox: false, projections: false });
        expect(cubic.mass).toBeGreaterThan(0.99);
        expect(Array.from(cubic.density)).toEqual(Array.from(plain.density));
        expect(new Set(cubic.boxHalfWidths).size).toBe(1);
        expect(plain.boxHalfWidths).toEqual(plain.halfWidths);
    });
});
