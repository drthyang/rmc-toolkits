// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Regression tests for the 1.0 symmetry-finder audit: the Detected SG card and its
// tolerance ladder must never show a wrong space-group symbol or number.

import { describe, it, expect } from 'vitest';

import { spaceGroupAtTolerance, symmetryLadder, findSpaceGroupOps } from '../symmetry.js';
import { conventionalCell } from '../symmetryModel.js';
import { demoStructure, shuffled, withNoise, STRUCTURES, at } from './fixtures/symmetryStructures.js';
import { cellVectors } from './fixtures/spaceGroups.js';

const ladderKey = (ladder) => ladder.map((b) => `${b.spaceGroup}[${b.from.toFixed(6)},${b.to.toFixed(6)}]`).join(' > ');

describe('operation residuals do not depend on the order of the sites', () => {
    const demo = demoStructure();
    const A = conventionalCell(demo);

    it('gives the same ladder for every ordering of the demo basis', { timeout: 60000 }, () => {
        const reference = ladderKey(symmetryLadder(A, demo.basis, 1.0));
        for (const seed of [1, 2, 3, 5]) {
            expect(ladderKey(symmetryLadder(A, shuffled(demo.basis, seed), 1.0)), `seed ${seed}`).toBe(reference);
        }
    });

    it('gives the same headline group for every ordering', { timeout: 60000 }, () => {
        for (const tol of [0.015, 0.02, 0.03, 0.05]) {
            const reference = spaceGroupAtTolerance(A, demo.basis, tol);
            for (const seed of [2, 3]) {
                const other = spaceGroupAtTolerance(A, shuffled(demo.basis, seed), tol);
                expect(other.spaceGroup, `tol ${tol} seed ${seed}`).toBe(reference.spaceGroup);
                expect(other.nSpace, `tol ${tol} seed ${seed}`).toBe(reference.nSpace);
            }
        }
    });

    it('refines each translation over all matched sites (least squares)', () => {
        // With noise on every site the translation of {R|t} must be the least-squares
        // one: the mean offset between each site's image and its matched partner is
        // zero. A translation read off one reference atom carries that atom's noise.
        const { A, basis } = withNoise(STRUCTURES.perovskite(), 0.02, 7);
        const { ops } = findSpaceGroupOps(A, basis, 0.2);
        expect(ops.length).toBe(48);
        for (const { R, t } of ops) {
            const mean = [0, 0, 0];
            for (const s of basis) {
                const img = [0, 1, 2].map((i) => R[i][0] * s.frac[0] + R[i][1] * s.frac[1] + R[i][2] * s.frac[2] + t[i]);
                let best = Infinity;
                let bestD = null;
                for (const o of basis) {
                    if (o.el !== s.el) continue;
                    const d = [0, 1, 2].map((i) => { const u = o.frac[i] - img[i]; return u - Math.round(u); });
                    const c = [0, 1, 2].map((j) => d[0] * A[0][j] + d[1] * A[1][j] + d[2] * A[2][j]);
                    const dist = Math.hypot(...c);
                    if (dist < best) { best = dist; bestD = d; }
                }
                for (let i = 0; i < 3; i++) mean[i] += bestD[i] / basis.length;
            }
            for (let i = 0; i < 3; i++) expect(Math.abs(mean[i])).toBeLessThan(1e-9);
        }
    });
});

describe('lattice strain is measured on the atomic-position scale', () => {
    const perovskiteIn = (A) => ({ A, basis: STRUCTURES.perovskite().basis });

    it('sees a 0.4 % tetragonal strain at a tolerance tighter than the strain', () => {
        // Fractional coordinates stay those of cubic Pm-3m, so ONLY the cell is strained:
        // c − a = 0.0156 Å. A tolerance below that must not report cubic symmetry.
        const a = 3.905;
        const { A, basis } = perovskiteIn(cellVectors(a, a, a * 1.004, 90, 90, 90));
        expect(spaceGroupAtTolerance(A, basis, 0.005)).toMatchObject({ spaceGroup: 'P4/mmm', spaceGroupNumber: 123 });
        expect(spaceGroupAtTolerance(A, basis, 0.05)).toMatchObject({ spaceGroup: 'Pm-3m', spaceGroupNumber: 221 });
        const ladder = symmetryLadder(A, basis, 1.0);
        expect(ladder.map((b) => b.spaceGroup)).toEqual(['P4/mmm', 'Pm-3m']);
        // the cubic rung starts where the strain is taken up: ≈ c − a
        expect(ladder[1].from).toBeGreaterThan(0.01);
        expect(ladder[1].from).toBeLessThan(0.02);
    });

    it('does not loosen the a/b test because c is long', () => {
        // 2.8 % a/b splitting (0.14 Å) in a cell with a 20 Å c axis.
        const basis = [at('X', 0, 0, 0), at('Y', 0.5, 0.5, 0.5)];
        const A = cellVectors(5.0, 5.14, 20.0, 90, 90, 90);
        expect(spaceGroupAtTolerance(A, basis, 0.01)).toMatchObject({ spaceGroup: 'Pmmm', spaceGroupNumber: 47 });
        expect(spaceGroupAtTolerance(A, basis, 0.3)).toMatchObject({ spaceGroup: 'P4/mmm', spaceGroupNumber: 123 });
    });

    it('keeps a slightly sheared cell a group at every tolerance', () => {
        // β = 90.3° (a 0.02 Å shear of the cell edges): monoclinic at tight tolerance,
        // cubic once the shear is within it — never P1 at every tolerance because the
        // accepted lattice operations failed to compose.
        const { A, basis } = perovskiteIn(cellVectors(3.9, 3.9, 3.91, 90, 90.3, 90));
        expect(spaceGroupAtTolerance(A, basis, 0.003)).toMatchObject({ spaceGroup: 'P2/m', spaceGroupNumber: 10 });
        expect(spaceGroupAtTolerance(A, basis, 0.2)).toMatchObject({ spaceGroup: 'Pm-3m', spaceGroupNumber: 221 });
        const ladder = symmetryLadder(A, basis, 1.0);
        expect(ladder[0].spaceGroup).toBe('P2/m');
        expect(ladder[ladder.length - 1].spaceGroup).toBe('Pm-3m');
    });
});
