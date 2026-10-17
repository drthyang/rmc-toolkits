// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import {
    extractPeaks,
    firstShell,
    IDEAL_COORDINATION,
    SHELL_TROUGH_FRACTION,
    SPLIT_SHELL_REACH_A,
    nearestNeighborDistances,
    pairCorrelationsContext
} from '../context/pairCorrelations';

// A 2×2×2 supercell of a 4 Å cubic cell with Na at the origin and Cl at
// [0.5, 0, 0]: Na–Cl nearest neighbour 2 Å, like pairs 4 Å (via images).
const rockSaltish = () => ({
    latticeVectors: [[8, 0, 0], [0, 8, 0], [0, 0, 8]],
    supercell: [2, 2, 2],
    basis: [
        { el: 'Na', frac: [0, 0, 0] },
        { el: 'Cl', frac: [0.5, 0, 0] }
    ]
});

const gaussian = (r, center, sigma, height) => (
    height * Math.exp(-((r - center) ** 2) / (2 * sigma ** 2))
);

const grid = (n, step) => Array.from({ length: n }, (_, i) => i * step);

describe('nearestNeighborDistances', () => {
    it('finds pair distances across periodic images', () => {
        const nn = nearestNeighborDistances(rockSaltish());
        expect(nn.get('Cl-Na')).toBeCloseTo(2, 6);
        expect(nn.get('Na-Na')).toBeCloseTo(4, 6);
        expect(nn.get('Cl-Cl')).toBeCloseTo(4, 6);
    });

    it('returns null without a basis', () => {
        expect(nearestNeighborDistances({ latticeVectors: [[1, 0, 0]], supercell: [1, 1, 1] })).toBeNull();
        expect(nearestNeighborDistances(null)).toBeNull();
    });
});

describe('extractPeaks', () => {
    it('finds the first peaks below rMax with plausible widths', () => {
        const x = grid(400, 0.02);   // 0 … 8 Å
        const y = x.map((r) => (
            gaussian(r, 2, 0.1, 5) + gaussian(r, 3.5, 0.15, 3) + gaussian(r, 7, 0.1, 4)
        ));
        const peaks = extractPeaks(x, y);
        expect(peaks).toHaveLength(2);
        expect(peaks[0].r).toBeCloseTo(2, 1);
        expect(peaks[1].r).toBeCloseTo(3.5, 1);
        // Gaussian FWHM = 2.355σ ≈ 0.24; smoothing broadens it slightly.
        expect(peaks[0].fwhm).toBeGreaterThan(0.15);
        expect(peaks[0].fwhm).toBeLessThan(0.4);
    });

    it('ignores noise below the floor and handles flat input', () => {
        const x = grid(300, 0.02);
        const y = x.map((r) => gaussian(r, 3, 0.1, 10) + 0.3 * Math.sin(r * 40));
        const peaks = extractPeaks(x, y);
        expect(peaks.length).toBeGreaterThanOrEqual(1);
        expect(peaks[0].r).toBeCloseTo(3, 1);
        expect(extractPeaks(x, x.map(() => 0))).toEqual([]);
        expect(extractPeaks([], [])).toEqual([]);
    });
});

describe('pairCorrelationsContext', () => {
    const partialsFile = () => {
        const x = grid(300, 0.02);
        return {
            plotKind: 'pdf_partials',
            plotData: {
                series: [
                    { label: 'Na-Cl', x, y: x.map((r) => gaussian(r, 2.05, 0.1, 6)) },
                    { label: 'Na-Na', x, y: x.map((r) => gaussian(r, 3.4, 0.12, 4)) },
                    { label: 'not a pair', x, y: x.map(() => 1) }
                ]
            }
        };
    };

    it('pairs average-structure distances with g(r) peaks', () => {
        const pairs = pairCorrelationsContext(rockSaltish(), [partialsFile()]);
        expect(pairs).toHaveLength(2);
        expect(pairs[0].pair).toBe('Na-Cl');
        expect(pairs[0].avg_structure_d_A).toBe(2);
        expect(pairs[0].gr_peaks_A[0].r).toBeCloseTo(2.05, 1);
        // Na-Na g(r) peak at 3.4 Å vs 4 Å average distance — the short-range signal.
        expect(pairs[1].avg_structure_d_A).toBe(4);
        expect(pairs[1].gr_peaks_A[0].r).toBeCloseTo(3.4, 1);
    });

    it('still reports peaks without a basis, and null without partials', () => {
        const pairs = pairCorrelationsContext(null, [partialsFile()]);
        expect(pairs[0].avg_structure_d_A).toBeUndefined();
        expect(pairs[0].gr_peaks_A.length).toBeGreaterThan(0);
        expect(pairCorrelationsContext(rockSaltish(), [])).toBeNull();
        expect(pairCorrelationsContext(rockSaltish(), [{ plotKind: 'npdf' }])).toBeNull();
    });
});

// Shells of Gaussian peaks, each holding `cn` partner atoms at density rho
// (CN = 4π ρ ∫ r² g dr; a narrow Gaussian at r0 holds 4π ρ h σ √(2π) r0²).
const shellG = (shells, rho, step = 0.005) => {
    const x = grid(1400, step);
    const y = x.map((r) => shells.reduce((sum, [r0, sigma, cn]) => (
        sum + gaussian(r, r0, sigma, cn / (4 * Math.PI * rho * sigma * Math.sqrt(2 * Math.PI) * r0 * r0))
    ), 0) + (r > 5 ? 1 : 0));
    return { x, y };
};

describe('firstShell: a split shell cut at a deep trough (review L1)', () => {
    const rho = 0.05;

    it('joins a Jahn–Teller octahedron (4 short + 2 long, deep trough) to CN 6, and keeps the first peak\'s own window', () => {
        const { x, y } = shellG([[1.93, 0.05, 4], [2.25, 0.05, 2]], rho);
        const shell = firstShell(x, y, { partnerDensity: rho });
        expect(shell.peaksInShell).toEqual([expect.closeTo(1.93, 1), expect.closeTo(2.25, 1)]);
        expect(shell.cn).toBeCloseTo(6, 1);
        expect(shell.boundaryR).toBeGreaterThan(2.35);
        expect(shell.split.firstCn).toBeCloseTo(4, 1);
        expect(shell.split.firstBoundaryR).toBeGreaterThan(2.0);
        expect(shell.split.firstBoundaryR).toBeLessThan(2.2);
        expect(shell.split.firstTroughRatio).toBeLessThan(SHELL_TROUGH_FRACTION);
        // Without a density there is no CN to judge by: the old cut stands.
        const blind = firstShell(x, y);
        expect(blind.split).toBeUndefined();
        expect(blind.boundaryR).toBeCloseTo(shell.split.firstBoundaryR, 5);
    });

    it('never joins a peak beyond the reach, nor one whose summed CN is no polyhedron', () => {
        const far = shellG([[2.0, 0.05, 6], [2.0 + SPLIT_SHELL_REACH_A + 0.1, 0.05, 2]], rho);
        expect(firstShell(far.x, far.y, { partnerDensity: rho }).cn).toBeCloseTo(6, 1);
        const odd = shellG([[2.0, 0.05, 4], [2.35, 0.05, 3]], rho);
        const shell = firstShell(odd.x, odd.y, { partnerDensity: rho });
        expect(shell.split).toBeUndefined();
        expect(shell.cn).toBeCloseTo(4, 1);
    });

    it('takes the furthest peak that completes a polyhedron (2 + 2 + 2), and counts either direction', () => {
        const ortho = shellG([[1.95, 0.03, 2], [2.12, 0.03, 2], [2.32, 0.03, 2]], rho);
        const shell = firstShell(ortho.x, ortho.y, { partnerDensity: rho });
        expect(shell.cn).toBeCloseTo(6, 1);
        expect(shell.split.firstCn).toBeCloseTo(2, 1);
        // 1 + 0.5 B per A is no polyhedron, but A is four times as dense: 6 A per B.
        const reverse = shellG([[1.93, 0.05, 1], [2.25, 0.05, 0.5]], rho);
        expect(firstShell(reverse.x, reverse.y, { partnerDensity: rho }).split).toBeUndefined();
        const joined = firstShell(reverse.x, reverse.y, { partnerDensity: rho, centreDensity: 4 * rho });
        expect(joined.cn).toBeCloseTo(1.5, 1);
        expect(joined.split.firstCn).toBeCloseTo(1, 1);
        expect(IDEAL_COORDINATION).toEqual([4, 6, 8, 12]);
    });
});
