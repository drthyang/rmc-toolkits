// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { nearestNeighborDistances, pairCorrelationsContext } from '../context/pairCorrelations';
import { shellCatalogue } from '../context/shells';

// The shell rules themselves (peaks, windows, CN, the first peak, joins) are
// tested with the module that reads them: shells.test.js.

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

// Shells of Gaussian peaks in r²g, each holding `cn` partner atoms at density
// rho (g = Σ cn N(r; r₀, σ) / (4π ρ r²)), 1 beyond 5 Å.
const shellG = (shells, rho, step = 0.005) => {
    const x = grid(1400, step);
    const y = x.map((r) => (r > 0 ? shells.reduce((sum, [r0, sigma, cn]) => (
        sum + cn * gaussian(r, r0, sigma, 1 / (sigma * Math.sqrt(2 * Math.PI)))
    ), 0) / (4 * Math.PI * rho * r * r) : 0) + (r > 5 ? 1 : 0));
    return { x, y };
};

describe('pairCorrelationsContext reads the run\'s one shell reading (shells.js)', () => {
    // 100 Ta and 200 Se in a 20 Å box: ρ(Se) = 0.025 Å⁻³.
    const structure = { latticeVectors: [[20, 0, 0], [0, 20, 0], [0, 0, 20]], supercell: [1, 1, 1], elementCounts: { Ta: 100, Se: 200 } };
    const split = () => {
        const { x, y } = shellG([[2.52, 0.06, 3], [2.74, 0.1, 3]], 0.025);
        return [{ plotKind: 'pdf_partials', plotData: { series: [{ label: 'Ta-Se', x, y }] } }];
    };

    it('gives the shell, its CN and, for two resolved peaks, where the first peak ends and what it holds', () => {
        const plotFiles = split();
        const [entry] = pairCorrelationsContext(structure, plotFiles);
        const shell = shellCatalogue({ plotFiles, structure }).shell('Ta', 'Se');
        expect(entry.gr_peaks_A.map((peak) => peak.r)).toEqual(shell.peaks.slice(0, 2).map((peak) => Number(peak.r.toPrecision(3))));
        expect(entry.first_shell.to_A).toBe(Number(shell.window[1].toPrecision(3)));
        expect(entry.first_shell.cn_Se_per_Ta).toBe(Number(shell.cn.toPrecision(3)));
        expect(entry.first_shell.cn_Se_per_Ta).toBeCloseTo(6, 1);
        expect(entry.first_shell.first_peak_to_A).toBe(Number(shell.firstPeak.window[1].toPrecision(3)));
        expect(entry.first_shell.first_peak_cn).toBe(Number(shell.firstPeak.cn.toPrecision(3)));
        // Separated from the next shell: no gap flag; no settings: no constraint flag.
        expect(entry.first_shell.gap_g_fraction).toBeUndefined();
        expect(entry.cut_at_min_distance).toBeUndefined();
    });

    it('flags the closest-approach pile-up from the run settings, in either pair order', () => {
        const x = grid(1000, 0.01);
        const plotFiles = [{
            plotKind: 'pdf_partials',
            plotData: { series: [{ label: 'Ta-Se', x, y: x.map((r) => (r < 2.4 ? 0 : gaussian(r, 2.5, 0.08, 8))) }] }
        }];
        const [entry] = pairCorrelationsContext(structure, plotFiles, { runSettings: { atoms: ['Se', 'Ta'], minimumDistancesA: [3, 2.4, 3] } });
        expect(entry.cut_at_min_distance.at_A).toBe(2.4);
        expect(entry.cut_at_min_distance.g_fraction).toBeGreaterThan(0.2);
        expect(entry.first_shell.first_peak_to_A).toBeUndefined();
    });
});
