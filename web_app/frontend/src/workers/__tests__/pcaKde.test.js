// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, it, expect } from 'vitest';

import {
    pcaKdeVolume,
    eigenDecomposition,
    probabilityScale,
    siteDisplacementsFromRmc6f,
    siteEllipsoids,
    sitePcaKde
} from '../pcaKde.js';

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

const anisotropicCloud = (n, seed, transform) => {
    const gauss = makeRng(seed);
    const points = [];
    for (let i = 0; i < n; i += 1) {
        const g = [gauss(), gauss(), gauss()];
        points.push([
            g[0] * transform[0][0] + g[1] * transform[1][0] + g[2] * transform[2][0],
            g[0] * transform[0][1] + g[1] * transform[1][1] + g[2] * transform[2][1],
            g[0] * transform[0][2] + g[1] * transform[1][2] + g[2] * transform[2][2]
        ]);
    }
    return points;
};

// Brute-force reference: the full anisotropic 3D Gaussian KDE with scipy's
// bandwidth matrix H = factor^2 * C, evaluated by direct summation. This is the
// ground truth the separable engine must reproduce.
const bruteForceDensity = (cloud, result) => {
    const n = cloud.length;
    const mean = result.mean;
    const cov = result.covariance;
    const factor = result.factor;
    const H = cov.map((row) => row.map((value) => value * factor * factor));

    // Invert the symmetric 3x3 bandwidth matrix and get its determinant.
    const [[a, b, c], [, d, e], [, , f]] = H;
    const det = a * (d * f - e * e) - b * (b * f - e * c) + c * (b * e - d * c);
    const inv = [
        [(d * f - e * e) / det, (c * e - b * f) / det, (b * e - c * d) / det],
        [(c * e - b * f) / det, (a * f - c * c) / det, (b * c - a * e) / det],
        [(b * e - c * d) / det, (b * c - a * e) / det, (a * d - b * b) / det]
    ];
    const norm = 1 / (n * (2 * Math.PI) ** 1.5 * Math.sqrt(det));

    const axes = result.axes;
    const coords = result.axisCoords;
    const grid = result.grid;
    const volume = new Float64Array(grid * grid * grid);
    for (let i = 0; i < grid; i += 1) {
        for (let j = 0; j < grid; j += 1) {
            for (let k = 0; k < grid; k += 1) {
                // PCA-frame grid point -> Cartesian.
                const p = [coords[0][i], coords[1][j], coords[2][k]];
                const x = mean[0] + p[0] * axes[0][0] + p[1] * axes[1][0] + p[2] * axes[2][0];
                const y = mean[1] + p[0] * axes[0][1] + p[1] * axes[1][1] + p[2] * axes[2][1];
                const z = mean[2] + p[0] * axes[0][2] + p[1] * axes[1][2] + p[2] * axes[2][2];
                let sum = 0;
                for (let m = 0; m < n; m += 1) {
                    const dx = x - cloud[m][0];
                    const dy = y - cloud[m][1];
                    const dz = z - cloud[m][2];
                    const q = inv[0][0] * dx * dx + inv[1][1] * dy * dy + inv[2][2] * dz * dz
                        + 2 * (inv[0][1] * dx * dy + inv[0][2] * dx * dz + inv[1][2] * dy * dz);
                    sum += Math.exp(-0.5 * q);
                }
                volume[(i * grid + j) * grid + k] = sum * norm;
            }
        }
    }
    return volume;
};

describe('probabilityScale', () => {
    it('matches the crystallographic 50% constant', () => {
        expect(probabilityScale(0.5)).toBeCloseTo(1.5381722, 3);
    });

    it('is monotonic and rejects out-of-range probabilities', () => {
        expect(probabilityScale(0.5)).toBeLessThan(probabilityScale(0.99));
        for (const bad of [0, 1, -0.1, 1.2]) {
            expect(() => probabilityScale(bad)).toThrow();
        }
    });
});

describe('eigenDecomposition', () => {
    it('recovers a known symmetric matrix, descending and right-handed', () => {
        // Diagonal in a rotated frame: eigenvalues 4, 2, 1.
        const cov = [
            [3.25, 0.4330127, -0.75],
            [0.4330127, 2.75, -0.4330127],
            [-0.75, -0.4330127, 3.0]
        ];
        // Adjust to a cleanly constructed matrix instead: rotate diag(4,2,1).
        const diag = [4, 2, 1];
        const theta = 0.6;
        const rot = [
            [Math.cos(theta), -Math.sin(theta), 0],
            [Math.sin(theta), Math.cos(theta), 0],
            [0, 0, 1]
        ];
        const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        for (let i = 0; i < 3; i += 1) {
            for (let j = 0; j < 3; j += 1) {
                let s = 0;
                for (let k = 0; k < 3; k += 1) s += rot[i][k] * diag[k] * rot[j][k];
                M[i][j] = s;
            }
        }
        const { eigenvalues, axes } = eigenDecomposition(M);
        expect(eigenvalues[0]).toBeCloseTo(4, 6);
        expect(eigenvalues[1]).toBeCloseTo(2, 6);
        expect(eigenvalues[2]).toBeCloseTo(1, 6);
        // Orthonormal, right-handed.
        const det = axes[0][0] * (axes[1][1] * axes[2][2] - axes[1][2] * axes[2][1])
            - axes[0][1] * (axes[1][0] * axes[2][2] - axes[1][2] * axes[2][0])
            + axes[0][2] * (axes[1][0] * axes[2][1] - axes[1][1] * axes[2][0]);
        expect(det).toBeCloseTo(1, 6);
        expect(cov).toBeDefined();
    });
});

describe('pcaKdeVolume separability', () => {
    const transform = [[1.4, 0.3, -0.2], [0.3, 0.7, 0.1], [-0.2, 0.1, 0.35]];

    it('matches the brute-force anisotropic 3D KDE to round-off', () => {
        const cloud = anisotropicCloud(300, 1, transform);
        for (const bw of ['scott', 'silverman', 0.35]) {
            const result = pcaKdeVolume(cloud, { bw, grid: 14, projections: false });
            const expected = bruteForceDensity(cloud, result);
            let maxErr = 0;
            let maxRel = 0;
            for (let i = 0; i < expected.length; i += 1) {
                const err = Math.abs(result.density[i] - expected[i]);
                if (err > maxErr) maxErr = err;
                if (expected[i] > 1e-6) maxRel = Math.max(maxRel, err / expected[i]);
            }
            expect(maxErr).toBeLessThan(1e-9);
            expect(maxRel).toBeLessThan(1e-9);
        }
    });

    it('projection equals a fine Riemann sum of the volume over the dropped axis', () => {
        const cloud = anisotropicCloud(300, 3, transform);
        const result = pcaKdeVolume(cloud, { bw: 'scott', grid: 80, extent: 5 });
        const grid = result.grid;
        const dz = result.axisCoords[2][1] - result.axisCoords[2][0];
        const proj = result.projections.pc12.density;
        let maxRel = 0;
        for (let i = 0; i < grid; i += 1) {
            for (let j = 0; j < grid; j += 1) {
                let integral = 0;
                for (let k = 0; k < grid; k += 1) integral += result.density[(i * grid + j) * grid + k];
                integral *= dz;
                if (integral > 1e-5) maxRel = Math.max(maxRel, Math.abs(proj[i][j] - integral) / integral);
            }
        }
        expect(maxRel).toBeLessThan(3e-3);
    });

    it('captures nearly all the mass in a wide box', () => {
        const cloud = anisotropicCloud(1500, 5, [[0.4, 0, 0], [0, 0.2, 0], [0, 0, 0.1]]);
        const result = pcaKdeVolume(cloud, { bw: 'scott', grid: 56, extent: 4, projections: false });
        expect(result.mass).toBeGreaterThan(0.99);
    });

    it('recovers a known anisotropy', () => {
        const cloud = anisotropicCloud(4000, 7, [[0.5, 0, 0], [0, 0.25, 0], [0, 0, 0.1]]);
        const result = pcaKdeVolume(cloud, { bw: 'scott', grid: 14, projections: false });
        expect(result.rms[0]).toBeCloseTo(0.5, 1);
        expect(result.rms[1]).toBeCloseTo(0.25, 1);
        expect(result.rms[2]).toBeCloseTo(0.1, 1);
    });

    it('cubic box uses one uniform half-width', () => {
        const cloud = anisotropicCloud(400, 9, [[0.6, 0, 0], [0, 0.2, 0], [0, 0, 0.05]]);
        const result = pcaKdeVolume(cloud, { bw: 'scott', grid: 10, cubicBox: true, projections: false });
        expect(result.halfWidths[0]).toBeCloseTo(result.halfWidths[1], 9);
        expect(result.halfWidths[1]).toBeCloseTo(result.halfWidths[2], 9);
    });

    it('mass levels bracket the cloud (higher probability -> lower threshold)', () => {
        const cloud = anisotropicCloud(1500, 12, [[0.3, 0, 0], [0, 0.3, 0], [0, 0, 0.3]]);
        const result = pcaKdeVolume(cloud, { bw: 'scott', grid: 40, projections: false });
        const levels = Object.fromEntries(result.massLevels.map((entry) => [entry.p, entry.level]));
        expect(levels[0.1]).toBeGreaterThan(levels[0.9]);
    });

    it('excess kurtosis is ~0 for a Gaussian cloud', () => {
        const cloud = anisotropicCloud(4000, 21, [[0.3, 0, 0], [0, 0.2, 0], [0, 0, 0.1]]);
        const result = pcaKdeVolume(cloud, { bw: 'scott', grid: 12, projections: false });
        expect(result.excessKurtosis).toHaveLength(3);
        expect(Math.abs(result.nonGaussianity)).toBeLessThan(0.4);
    });

    it('rejects degenerate input', () => {
        expect(() => pcaKdeVolume(Array.from({ length: 10 }, () => [0, 0, 0]), { grid: 8 })).toThrow();
        expect(() => pcaKdeVolume([[0, 0, 0], [1, 1, 1]], { grid: 8 })).toThrow();
    });
});

const buildRmc6f = (atomLines, supercell = [6, 6, 6], edge = 8) => {
    const L = edge;
    return [
        `Supercell dimensions ${supercell[0]} ${supercell[1]} ${supercell[2]}`,
        'Lattice vectors (Ang):',
        `${L * supercell[0]} 0 0`,
        `0 ${L * supercell[1]} 0`,
        `0 0 ${L * supercell[2]}`,
        'Atoms:',
        ...atomLines
    ].join('\n');
};

const cloudRmc6f = (sigmaFrac, { supercell = [6, 6, 6], seed = 1, element = 'Ga', edge = 8 } = {}) => {
    const gauss = makeRng(seed);
    const lines = [];
    let atom = 0;
    for (let ix = 0; ix < supercell[0]; ix += 1) {
        for (let iy = 0; iy < supercell[1]; iy += 1) {
            for (let iz = 0; iz < supercell[2]; iz += 1) {
                atom += 1;
                const base = [ix / supercell[0], iy / supercell[1], iz / supercell[2]];
                const disp = [gauss() * sigmaFrac[0], gauss() * sigmaFrac[1], gauss() * sigmaFrac[2]];
                const coord = base.map((value, i) => value + disp[i]);
                lines.push(`${atom} ${element} [1] ${coord[0].toFixed(10)} ${coord[1].toFixed(10)} `
                    + `${coord[2].toFixed(10)} 1 ${ix} ${iy} ${iz}`);
            }
        }
    }
    return buildRmc6f(lines, supercell, edge);
};

describe('siteDisplacementsFromRmc6f', () => {
    it('extracts a centered cloud with the expected Cartesian sigma', () => {
        const sigmaFrac = [0.02, 0.01, 0.005];
        const edge = 8;
        const text = cloudRmc6f(sigmaFrac, { supercell: [8, 8, 8], seed: 4, edge });
        const parsed = siteDisplacementsFromRmc6f(text);
        expect(parsed.referenceNumbers).toEqual([1]);
        const site = parsed.sites[0];
        expect(site.count).toBe(512);
        // Centered.
        const mean = site.displacements.reduce(
            (acc, d) => [acc[0] + d[0], acc[1] + d[1], acc[2] + d[2]], [0, 0, 0]
        ).map((value) => value / site.count);
        mean.forEach((value) => expect(Math.abs(value)).toBeLessThan(1e-9));
        // Displacements are supercell fractions, so Cartesian sigma = fractional
        // sigma * supercell edge (edge per cell * number of cells).
        const supercellEdge = edge * 8;
        const variance = site.displacements.reduce(
            (acc, d) => [acc[0] + d[0] * d[0], acc[1] + d[1] * d[1], acc[2] + d[2] * d[2]], [0, 0, 0]
        ).map((value) => value / (site.count - 1));
        expect(Math.sqrt(variance[0])).toBeCloseTo(sigmaFrac[0] * supercellEdge, 1);
    });

    it('folds an atom that drifted across the supercell boundary', () => {
        const text = buildRmc6f([
            '1 Ga [1] 0.999 0.001 0.500 1 0 0 0',
            '2 Ga [1] 0.001 0.999 0.500 1 0 0 0',
            '3 Ga [1] 0.000 0.000 0.500 1 0 0 0',
            '4 Ga [1] 0.998 0.002 0.500 1 0 0 0'
        ], [1, 1, 1], 8);
        const parsed = siteDisplacementsFromRmc6f(text);
        const maxDisp = Math.max(...parsed.sites[0].displacements.flat().map(Math.abs));
        expect(maxDisp).toBeLessThan(0.1);
    });

    it('site ellipsoids match a direct covariance and sort descending', () => {
        const text = cloudRmc6f([0.03, 0.015, 0.008], { supercell: [6, 6, 6], seed: 8 });
        const parsed = siteDisplacementsFromRmc6f(text);
        const ellipsoids = siteEllipsoids(parsed.sites);
        expect(ellipsoids).toHaveLength(1);
        const eig = ellipsoids[0].eigenvalues;
        expect(eig[0]).toBeGreaterThanOrEqual(eig[1]);
        expect(eig[1]).toBeGreaterThanOrEqual(eig[2]);
    });

    it('sitePcaKde tags reference number and element', () => {
        const text = cloudRmc6f([0.02, 0.02, 0.02], { supercell: [6, 6, 6], seed: 2, element: 'Nb' });
        const parsed = siteDisplacementsFromRmc6f(text);
        const result = sitePcaKde(parsed, { referenceNumber: 1, grid: 12, projections: false });
        expect(result.referenceNumber).toBe(1);
        expect(result.element).toBe('Nb');
        expect(result.count).toBe(216);
    });
});
