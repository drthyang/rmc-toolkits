// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, it, expect } from 'vitest';

import {
    transpose3,
    multiply3,
    matVec3,
    invert3,
    unitCellVectors,
    crystalPcaTransforms,
    principalAxisOrientation
} from '../pcaCrystalFrame.js';

const IDENTITY = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

const closeMatrix = (actual, expected, digits = 9) => {
    for (let i = 0; i < 3; i += 1) {
        for (let j = 0; j < 3; j += 1) {
            expect(actual[i][j]).toBeCloseTo(expected[i][j], digits);
        }
    }
};

describe('3×3 linear algebra', () => {
    it('inverts a matrix and round-trips through its inverse', () => {
        const m = [[2, 1, 0], [0, 3, 1], [1, 0, 4]];
        const inv = invert3(m);
        closeMatrix(multiply3(m, inv), IDENTITY);
        closeMatrix(multiply3(inv, m), IDENTITY);
    });

    it('returns null for a singular matrix', () => {
        expect(invert3([[1, 2, 3], [2, 4, 6], [0, 1, 0]])).toBeNull();
    });

    it('transpose and matVec agree with the definition', () => {
        const m = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
        closeMatrix(transpose3(m), [[1, 4, 7], [2, 5, 8], [3, 6, 9]]);
        expect(matVec3(m, [1, 0, -1])).toEqual([-2, -2, -2]);
    });
});

describe('unitCellVectors', () => {
    it('divides the supercell lattice by the repeat counts', () => {
        const lattice = [[16, 0, 0], [0, 24, 0], [0, 0, 40]];
        const unit = unitCellVectors(lattice, [2, 3, 5]);
        closeMatrix(unit, [[8, 0, 0], [0, 8, 0], [0, 0, 8]]);
    });
});

describe('crystalPcaTransforms', () => {
    // Orthogonal cubic cell aligned with Cartesian axes: M = 8·I, so the
    // fractional→PCA matrix is just 8·P and the direction cosines equal P itself.
    const unit = [[8, 0, 0], [0, 8, 0], [0, 0, 8]];

    it('fracToPca = P·Mᵀ and pcaToFrac is its exact inverse', () => {
        // A right-handed rotation about z by 30° as the PCA axes (rows).
        const t = Math.PI / 6;
        const axes = [
            [Math.cos(t), Math.sin(t), 0],
            [-Math.sin(t), Math.cos(t), 0],
            [0, 0, 1]
        ];
        const { fracToPca, pcaToFrac } = crystalPcaTransforms(axes, unit);
        // Manually: P·Mᵀ = 8·P for a cubic cell.
        closeMatrix(fracToPca, axes.map((row) => row.map((v) => 8 * v)));
        // Round-trip is identity.
        closeMatrix(multiply3(fracToPca, pcaToFrac), IDENTITY);
        closeMatrix(multiply3(pcaToFrac, fracToPca), IDENTITY);
    });

    it('round-trips for an oblique (monoclinic) cell too', () => {
        // β = 105° monoclinic cell: a along x, b along y, c in the x–z plane.
        const beta = (105 * Math.PI) / 180;
        const cLen = 6;
        const unitMono = [
            [5, 0, 0],
            [0, 7, 0],
            [cLen * Math.cos(beta), 0, cLen * Math.sin(beta)]
        ];
        // A clean orthonormal PCA triad (columns of a rotation matrix).
        const axes = [
            [1 / Math.sqrt(3), 1 / Math.sqrt(3), 1 / Math.sqrt(3)],
            [-1 / Math.sqrt(2), 1 / Math.sqrt(2), 0],
            [1 / Math.sqrt(6), 1 / Math.sqrt(6), -2 / Math.sqrt(6)]
        ];
        const { fracToPca, pcaToFrac } = crystalPcaTransforms(axes, unitMono);
        closeMatrix(multiply3(fracToPca, pcaToFrac), IDENTITY);
        closeMatrix(multiply3(pcaToFrac, fracToPca), IDENTITY);
    });

    it('returns null for a collapsed cell', () => {
        const flat = [[8, 0, 0], [0, 8, 0], [0, 0, 0]];
        expect(crystalPcaTransforms(IDENTITY, flat)).toBeNull();
    });
});

describe('principalAxisOrientation', () => {
    const cubic = [[8, 0, 0], [0, 8, 0], [0, 0, 8]];

    it('reports 0°/90° and the right [u v w] when PCA axes align with a, b, c', () => {
        const axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        const { perAxis } = principalAxisOrientation(axes, cubic);
        // PC1 ∥ a
        expect(perAxis[0].anglesDeg[0]).toBeCloseTo(0, 6);
        expect(perAxis[0].anglesDeg[1]).toBeCloseTo(90, 6);
        expect(perAxis[0].anglesDeg[2]).toBeCloseTo(90, 6);
        expect(perAxis[0].dominant.label).toBe('a');
        expect(perAxis[0].dominant.angleDeg).toBeCloseTo(0, 6);
        expect(perAxis[0].crystalDirection).toEqual([1, 0, 0]);
        // PC3 ∥ c
        expect(perAxis[2].dominant.label).toBe('c');
        expect(perAxis[2].crystalDirection).toEqual([0, 0, 1]);
    });

    it('gives the direction cosine for a 30° tilt about z', () => {
        const t = Math.PI / 6;
        const axes = [
            [Math.cos(t), Math.sin(t), 0],
            [-Math.sin(t), Math.cos(t), 0],
            [0, 0, 1]
        ];
        const { perAxis } = principalAxisOrientation(axes, cubic);
        // PC1 makes 30° with a and 60° with b.
        expect(perAxis[0].anglesDeg[0]).toBeCloseTo(30, 6);
        expect(perAxis[0].anglesDeg[1]).toBeCloseTo(60, 6);
        expect(perAxis[0].cosines[0]).toBeCloseTo(Math.cos(t), 6);
    });

    it('crystalDirection of a PCA axis is a true lattice direction (oblique cell)', () => {
        // Monoclinic cell; take a PCA axis exactly along the c edge and confirm the
        // recovered [u v w] is [0 0 1].
        const beta = (110 * Math.PI) / 180;
        const cLen = 6.2;
        const unitMono = [
            [5.1, 0, 0],
            [0, 7.3, 0],
            [cLen * Math.cos(beta), 0, cLen * Math.sin(beta)]
        ];
        const cHat = unitMono[2].map((v) => v / Math.hypot(...unitMono[2]));
        // Build an orthonormal PCA set whose first axis is ĉ.
        const aHat = [1, 0, 0];
        const bHat = [
            aHat[1] * cHat[2] - aHat[2] * cHat[1],
            aHat[2] * cHat[0] - aHat[0] * cHat[2],
            aHat[0] * cHat[1] - aHat[1] * cHat[0]
        ];
        const bLen = Math.hypot(...bHat);
        const axes = [cHat, bHat.map((v) => v / bLen), [
            cHat[1] * bHat[2] - cHat[2] * bHat[1],
            cHat[2] * bHat[0] - cHat[0] * bHat[2],
            cHat[0] * bHat[1] - cHat[1] * bHat[0]
        ].map((v) => v / bLen)];
        const { perAxis } = principalAxisOrientation(axes, unitMono);
        // PC1 points along the crystallographic c direction → [0 0 1].
        expect(perAxis[0].crystalDirection[0]).toBeCloseTo(0, 6);
        expect(perAxis[0].crystalDirection[1]).toBeCloseTo(0, 6);
        expect(perAxis[0].crystalDirection[2]).toBeCloseTo(1, 6);
    });
});
