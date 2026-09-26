// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Crystal-frame wall projections (pca.numerics.16). The page used to splat every
// grid node bilinearly into 2D bins of the node spacing; when the PCA axes are
// rotated against the frame the projected node lattice beats against the bins,
// a moiré of up to ~21% (6-9% max error on the per-axis volume). The walls are now
// the line-integral marginal of the same volume (projectVolumeOntoFrame); this
// pins its error against the exact analytic 2D marginal of the KDE.

import { describe, it, expect } from 'vitest';

import { projectVolumeOntoFrame } from '../pcaCrystalFrame.js';
import { pcaKdeVolume } from '../workers/pcaKde.js';

const makeRandom = (seed) => {
    let s = seed >>> 0;
    const uniform = () => { s = (Math.imul(s, 1664525) + 1013904223) >>> 0; return s / 4294967296; };
    return () => Math.sqrt(-2 * Math.log(uniform() + 1e-300)) * Math.cos(2 * Math.PI * uniform());
};

const rotatedCloud = (degrees, sigma, n = 1000, seed = 7) => {
    const normal = makeRandom(seed);
    const th = (degrees * Math.PI) / 180;
    return Array.from({ length: n }, () => {
        const x = normal() * sigma[0];
        const y = normal() * sigma[1];
        const z = normal() * sigma[2];
        return [x * Math.cos(th) - y * Math.sin(th), x * Math.sin(th) + y * Math.cos(th), z];
    });
};

// Exact marginal of the KDE on the x-y plane: a 2D Gaussian KDE whose bandwidth is
// the x-y block of the Cartesian bandwidth P^T diag(h^2) P.
const exactXyMarginal = (cloud, kde, half, nBins) => {
    const P = kde.axes;
    const h2 = kde.bandwidth.map((b) => b * b);
    const H = [[0, 0], [0, 0]];
    for (let a = 0; a < 3; a += 1) {
        for (let r = 0; r < 2; r += 1) for (let c = 0; c < 2; c += 1) H[r][c] += P[a][r] * h2[a] * P[a][c];
    }
    const det = H[0][0] * H[1][1] - H[0][1] * H[1][0];
    const inv = [[H[1][1] / det, -H[0][1] / det], [-H[1][0] / det, H[0][0] / det]];
    const step = (2 * half) / (nBins - 1);
    const out = [];
    for (let i = 0; i < nBins; i += 1) {
        const row = [];
        for (let j = 0; j < nBins; j += 1) {
            const u = -half + i * step;
            const v = -half + j * step;
            let acc = 0;
            cloud.forEach((pt) => {
                const du = u - (pt[0] - kde.mean[0]);
                const dv = v - (pt[1] - kde.mean[1]);
                acc += Math.exp(-0.5 * (du * (inv[0][0] * du + inv[0][1] * dv) + dv * (inv[1][0] * du + inv[1][1] * dv)));
            });
            row.push(acc / (cloud.length * 2 * Math.PI * Math.sqrt(det)));
        }
        out.push(row);
    }
    return out;
};

const IDENTITY_FRAME = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

const errors = (wall, exact) => {
    let peak = 0;
    exact.forEach((row) => row.forEach((v) => { peak = Math.max(peak, v); }));
    let maxRel = 0;
    let maxAbs = 0;
    let sumAbs = 0;
    let sumRef = 0;
    exact.forEach((row, i) => row.forEach((e, j) => {
        const w = wall[i][j];
        sumAbs += Math.abs(w - e);
        sumRef += e;
        maxAbs = Math.max(maxAbs, Math.abs(w - e) / peak);
        if (e > 0.1 * peak) maxRel = Math.max(maxRel, Math.abs(w / e - 1));
    }));
    return { maxRel, maxAbs, l1: sumAbs / sumRef };
};

describe('projectVolumeOntoFrame', () => {
    [0, 30, 45].forEach((degrees) => {
        it(`matches the exact marginal with PCs rotated ${degrees} deg against the frame (grid 40)`, () => {
            const cloud = rotatedCloud(degrees, [0.12, 0.07, 0.05]);
            const kde = pcaKdeVolume(cloud, { grid: 40, extent: 4, cubicBox: true, projections: false });
            const half = Math.max(...kde.halfWidths);
            const walls = projectVolumeOntoFrame(kde, IDENTITY_FRAME, half, 40);
            const { maxRel, maxAbs, l1 } = errors(walls.pc12.density, exactXyMarginal(cloud, kde, half, 40));
            expect(l1).toBeLessThan(0.015);
            expect(maxAbs).toBeLessThan(0.025);
            expect(maxRel).toBeLessThan(0.04);
        });
    });

    it('stays within a few percent on the coarsest UI grid and for a thin cloud', () => {
        const coarse = rotatedCloud(45, [0.12, 0.07, 0.05]);
        const kdeCoarse = pcaKdeVolume(coarse, { grid: 24, extent: 4, cubicBox: true, projections: false });
        const halfCoarse = Math.max(...kdeCoarse.halfWidths);
        const coarseErr = errors(projectVolumeOntoFrame(kdeCoarse, IDENTITY_FRAME, halfCoarse, 24).pc12.density,
            exactXyMarginal(coarse, kdeCoarse, halfCoarse, 24));
        expect(coarseErr.l1).toBeLessThan(0.03);
        expect(coarseErr.maxAbs).toBeLessThan(0.04);

        const thin = rotatedCloud(45, [0.3, 0.1, 0.02], 1000, 11);
        const kdeThin = pcaKdeVolume(thin, { grid: 40, extent: 4, cubicBox: true, projections: false });
        const halfThin = Math.max(...kdeThin.halfWidths);
        const thinErr = errors(projectVolumeOntoFrame(kdeThin, IDENTITY_FRAME, halfThin, 40).pc12.density,
            exactXyMarginal(thin, kdeThin, halfThin, 40));
        expect(thinErr.l1).toBeLessThan(0.015);
        expect(thinErr.maxRel).toBeLessThan(0.04);
    });

    it('conserves the captured mass on every plane', () => {
        const cloud = rotatedCloud(30, [0.12, 0.07, 0.05]);
        const kde = pcaKdeVolume(cloud, { grid: 40, extent: 4, cubicBox: true, projections: false });
        const half = Math.max(...kde.halfWidths);
        const walls = projectVolumeOntoFrame(kde, IDENTITY_FRAME, half, 40);
        const texel = (2 * half) / 39;
        ['pc12', 'pc13', 'pc23'].forEach((key) => {
            let mass = 0;
            walls[key].density.forEach((row) => row.forEach((v) => { mass += v * texel * texel; }));
            expect(Math.abs(mass - kde.mass)).toBeLessThan(0.02);
        });
        expect(walls.pc13.axes).toEqual([0, 2]);
        expect(walls.pc23.extent).toEqual([-half, half, -half, half]);
    });
});
