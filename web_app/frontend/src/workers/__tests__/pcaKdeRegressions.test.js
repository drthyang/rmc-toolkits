// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Regression tests for the browser PCA-ellipsoid engine (workers/pcaKde.js),
// one describe block per defect found by the 1.0 audit. They mirror
// tests/test_pca_regressions.py so both engines are pinned to the same rule.

import { describe, it, expect } from 'vitest';

import { siteDisplacementsFromRmc6f, siteEllipsoids } from '../pcaKde.js';

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
