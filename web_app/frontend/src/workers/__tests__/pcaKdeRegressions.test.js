// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Regression tests for the browser PCA-ellipsoid engine (workers/pcaKde.js),
// one describe block per defect found by the 1.0 audit. They mirror
// tests/test_pca_regressions.py so both engines are pinned to the same rule.

import { existsSync, readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, it, expect } from 'vitest';

import {
    eigenDecomposition,
    pcaKdeVolume,
    siteDisplacementsFromRmc6f,
    siteEllipsoids,
    sitePcaKde
} from '../pcaKde.js';

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
