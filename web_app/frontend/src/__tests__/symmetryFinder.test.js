// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Regression tests for the 1.0 symmetry-finder audit: the Detected SG card and its
// tolerance ladder must never show a wrong space-group symbol or number.

import { describe, it, expect } from 'vitest';

import { spaceGroupAtTolerance, symmetryLadder, findSpaceGroupOps } from '../symmetry.js';
import { conventionalCell } from '../symmetryModel.js';
import {
    demoStructure, shuffled, withNoise, STRUCTURES, at, closureDefects, lacunarSpinel,
} from './fixtures/symmetryStructures.js';
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

describe('every reported operation set is a group', () => {
    // A ladder rung or headline is a set of operations {R|t} kept because each one's
    // residual is within the tolerance. With noise those residuals differ, so the kept
    // set is an arbitrary subset of the true group unless closure is enforced.
    const midpoints = (ladder) => ladder.map((b) => (b.from + b.to) / 2);

    it('closes every rung of the bundled demo ladder', { timeout: 60000 }, () => {
        const demo = demoStructure();
        const A = conventionalCell(demo);
        const ladder = symmetryLadder(A, demo.basis, 1.0);
        expect(ladder[ladder.length - 1]).toMatchObject({ spaceGroup: 'F-43m', spaceGroupNumber: 216 });
        for (const tol of midpoints(ladder)) {
            const found = spaceGroupAtTolerance(A, demo.basis, tol);
            expect(closureDefects(found.ops), `${found.spaceGroup} at ${tol.toFixed(4)} Å`).toHaveLength(0);
        }
    });

    it('closes the headline of a noisy lacunar spinel at every tolerance', { timeout: 60000 }, () => {
        const { A, basis } = withNoise(lacunarSpinel(), 0.02, 3);
        for (const tol of [0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2]) {
            const found = spaceGroupAtTolerance(A, basis, tol);
            expect(closureDefects(found.ops), `${found.spaceGroup} at ${tol} Å`).toHaveLength(0);
            expect(found.nSpace).toBe(found.ops.length);
        }
        expect(spaceGroupAtTolerance(A, basis, 0.2)).toMatchObject({ spaceGroup: 'F-43m', spaceGroupNumber: 216 });
    });

    it('reports the worst residual of the group it returns', () => {
        const { A, basis } = withNoise(STRUCTURES.perovskite(), 0.01, 11);
        const found = spaceGroupAtTolerance(A, basis, 0.2);
        expect(found.maxResidual).toBeCloseTo(Math.max(...found.ops.map((o) => o.residual)), 12);
    });

    it('never lets the operation count fall as the tolerance loosens', { timeout: 60000 }, () => {
        const { A, basis } = withNoise(lacunarSpinel(), 0.02, 3);
        const ladder = symmetryLadder(A, basis, 1.0);
        expect(ladder[0].from).toBe(0);
        for (let i = 1; i < ladder.length; i += 1) {
            expect(ladder[i].from).toBeCloseTo(ladder[i - 1].to, 12);
            expect(ladder[i].nSpace).toBeGreaterThan(ladder[i - 1].nSpace);
        }
    });
});
