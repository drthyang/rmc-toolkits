// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The run's one shell reading (context/shells.js, architecture candidate #6):
// a known answer for every rule — the search range, the rounding, the CN
// integral, the first-peak window (the inter-peak minimum, approved
// 2026-10-19), separation, the closest-approach pile-up (the FWHM edge clause
// dropped), the CN match, the one cut rule (windowCut: the CN margin where a
// window meets a shell, the counts where there is none to compare; approved
// by the user, 2026-10-19), the split-shell joins of the review's L1 (moved
// here from pairCorrelations) and a joined shell's split by number (plan item
// 4.2): its parts, long − short and its sampling error.

import { describe, expect, it } from 'vitest';
import {
    BOUNDARY_REACH_A,
    CN_CHANGE_ABS,
    CN_CHANGE_REL,
    CN_MATCH_TOLERANCE,
    CONSTRAINT_LIMIT_FRACTION,
    CUT_SHIFT_A,
    IDEAL_COORDINATION,
    MODAL_COORDINATION_SHARE,
    PEAK_SEARCH_MAX_A,
    SEPARATION_G_FRACTION,
    SHELL_TROUGH_FRACTION,
    SPLIT_SHELL_REACH_A,
    UNBONDED_APEX_FRACTION,
    AVERAGE_SHELL_REACH_A,
    averageShellCheck,
    averageStructureBonds,
    cnChanges,
    cnMarginCut,
    cnMatches,
    countsCut,
    elementDensity,
    extractPeaks,
    isIdealCoordination,
    labelledMinimumDistances,
    shellCatalogue,
    windowCut
} from '../context/shells';

const grid = (n, step) => Array.from({ length: n }, (_, i) => Number((i * step).toFixed(6)));
const gaussian = (r, center, sigma, height) => height * Math.exp(-((r - center) ** 2) / (2 * sigma ** 2));

// A cubic L Å box holding these atom counts (density count / L³).
const box = (elementCounts, L = 10) => ({ latticeVectors: [[L, 0, 0], [0, L, 0], [0, 0, L]], elementCounts });
const partials = (series) => [{ name: 'p.csv', plotKind: 'pdf_partials', plotData: { series } }];

// A partial g(r) of Gaussian shells in r²g: shell i holds exactly cn_i partner
// atoms at density rho (g = Σ cn N(r; r₀, σ) / (4π ρ r²)), plus 1 beyond `tail`.
const shellSeries = (label, shells, rho, { step = 0.01, rMax = 10, tail = 5.5 } = {}) => {
    const x = grid(Math.round(rMax / step), step);
    const y = x.map((r) => (r > 0
        ? shells.reduce((sum, [r0, sigma, cn]) => sum + cn * gaussian(r, r0, sigma, 1 / (sigma * Math.sqrt(2 * Math.PI))), 0) / (4 * Math.PI * rho * r * r)
        : 0) + (r > tail ? 1 : 0));
    return { label, x, y };
};

// 50 partner atoms in 1000 Å³: ρ = 0.05 Å⁻³.
const RHO = 0.05;
const catalogueOf = (series, { structure = box({ Ta: 50, Se: 50 }), runSettings = null } = {}) => shellCatalogue({
    plotFiles: partials(Array.isArray(series) ? series : [series]), structure, runSettings
});

describe('extractPeaks (the display peaks; the module\'s search below 8 Å)', () => {
    it('finds the first peaks below rMax with plausible widths', () => {
        const x = grid(400, 0.02); // 0 … 8 Å
        const y = x.map((r) => gaussian(r, 2, 0.1, 5) + gaussian(r, 3.5, 0.15, 3) + gaussian(r, 7, 0.1, 4));
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

describe('the first-peak window: the inter-peak minimum (approved 2026-10-19)', () => {
    // The known answer: an octahedron split 3 + 3, at 2.5 and 2.75 Å, each
    // shell σ = 0.07 Å, so the two overlap (the trough between them stays at
    // ~40 % of the first peak and no deep minimum separates them). The first
    // peak alone holds 3 by construction.
    const SPLIT = [[2.5, 0.07, 3], [2.75, 0.07, 3]];
    const halfWidthWindow = (series) => {
        const [first] = extractPeaks(series.x, series.y, { rMax: PEAK_SEARCH_MAX_A, maxPeaks: 3 });
        return [Math.max(0, first.r - first.fwhm), first.r + first.fwhm];
    };

    it.each([0.005, 0.01, 0.02])('a 3 + 3 split with overlapping widths (grid %s Å): the first peak holds 3 ± 0.1, r ± FWHM overcounts', (step) => {
        const series = shellSeries('Ta-Se', SPLIT, RHO, { step });
        const catalogue = catalogueOf(series);
        const shell = catalogue.shell('Ta', 'Se');
        expect(shell.peaksInShell).toHaveLength(2);
        expect(shell.joined).toBe(false);
        expect(Math.abs(shell.cn - 6)).toBeLessThanOrEqual(0.1);
        const { firstPeak } = shell;
        expect(Math.abs(firstPeak.cn - 3)).toBeLessThanOrEqual(0.1);
        // From the shell's start to the minimum between the two peaks.
        expect(firstPeak.window[0]).toBe(shell.window[0]);
        expect(firstPeak.window[1]).toBeGreaterThan(shell.peaksInShell[0]);
        expect(firstPeak.window[1]).toBeLessThan(shell.peaksInShell[1]);
        expect(firstPeak.troughRatio).toBeGreaterThan(SHELL_TROUGH_FRACTION);
        expect(firstPeak.changesCn).toBe(true);
        // The rule it replaces: the half-height walk reaches into the second
        // peak (3.40–4.00 on these grids).
        const old = catalogue.window('Ta', 'Se', halfWidthWindow(series)).cn;
        expect(old - 3).toBeGreaterThan(0.3);
    });

    it('where the trough stays above half the first peak, r ± FWHM holds the whole shell; the minimum still finds the first', () => {
        const series = shellSeries('Ta-Se', [[2.5, 0.08, 3], [2.75, 0.08, 3]], RHO);
        const catalogue = catalogueOf(series);
        const shell = catalogue.shell('Ta', 'Se');
        expect(shell.firstPeak.troughRatio).toBeGreaterThan(0.5);
        expect(Math.abs(shell.firstPeak.cn - 3)).toBeLessThanOrEqual(0.1);
        expect(catalogue.window('Ta', 'Se', halfWidthWindow(series)).cn).toBeGreaterThan(5.9);
    });

    it('a sharp short shell and a broad long one (the Demo\'s shape): the minimum comes closer to 3 than r ± FWHM', () => {
        const series = shellSeries('Ta-Se', [[2.52, 0.06, 3], [2.74, 0.1, 3]], RHO);
        const catalogue = catalogueOf(series);
        const { firstPeak } = catalogue.shell('Ta', 'Se');
        const old = catalogue.window('Ta', 'Se', halfWidthWindow(series)).cn;
        expect(Math.abs(firstPeak.cn - 3)).toBeLessThan(Math.abs(old - 3));
    });

    it('a shell of one resolved peak has no first peak (no r ± FWHM fallback)', () => {
        const shell = catalogueOf(shellSeries('Ta-Se', [[2.6, 0.08, 6]], RHO)).shell('Ta', 'Se');
        expect(shell.peaksInShell).toHaveLength(1);
        expect(shell.firstPeak).toBeNull();
    });

    it('three resolved peaks: the first peak ends at the minimum between the first two', () => {
        const shell = catalogueOf(shellSeries('Se-Se', [[3.3, 0.07, 4], [3.55, 0.07, 4], [3.8, 0.07, 4]], RHO)).shell('Se', 'Se');
        expect(shell.peaksInShell).toHaveLength(3);
        expect(shell.joined).toBe(false);
        expect(Math.abs(shell.cn - 12)).toBeLessThanOrEqual(0.1);
        expect(shell.firstPeak.window[1]).toBeGreaterThan(3.3);
        expect(shell.firstPeak.window[1]).toBeLessThan(3.55);
        expect(Math.abs(shell.firstPeak.cn - 4)).toBeLessThanOrEqual(0.1);
    });
});

describe('one rounding rule and one CN integral', () => {
    it('rounds every window to 0.01 Å, then integrates exactly over it: each CN is window().cn of its window', () => {
        // A 0.005 Å grid, where rounding moves the edges.
        const catalogue = catalogueOf(shellSeries('Ta-Se', [[2.52, 0.06, 3], [2.74, 0.1, 3]], RHO, { step: 0.005 }));
        const shell = catalogue.shell('Ta', 'Se');
        const onGrid = (value) => Math.abs(value * 100 - Math.round(value * 100)) < 1e-9;
        const readings = [
            [shell.window, shell.cn],
            [shell.firstPeak.window, shell.firstPeak.cn],
            [shell.cut.earlier.window, shell.cut.earlier.cn],
            [shell.cut.later.window, shell.cut.later.cn]
        ];
        for (const [window, cn] of readings) {
            expect(window.every(onGrid)).toBe(true);
            expect(catalogue.window('Ta', 'Se', window).cn).toBe(cn);
        }
        expect(shell.cut.earlier.window[1]).toBeCloseTo(shell.window[1] - CUT_SHIFT_A, 9);
        expect(shell.cut.later.window[1]).toBeCloseTo(shell.window[1] + CUT_SHIFT_A, 9);
    });

    it('a shell narrower than the cut: the earlier cut stops at the shell\'s start and holds no atoms', () => {
        // A dimer bond (1 partner atom) at 2.26 Å, σ 0.02 Å: the shell is ~0.15 Å wide.
        const shell = catalogueOf(shellSeries('V-V', [[2.26, 0.02, 1]], RHO, { step: 0.005 }), { structure: box({ V: 50 }) }).shell('V', 'V');
        expect(shell.window[1] - shell.window[0]).toBeLessThan(CUT_SHIFT_A);
        expect(shell.cut.earlier.window).toEqual([shell.window[0], shell.window[0]]);
        expect(shell.cut.earlier.cn).toBe(0);
        expect(shell.cut.earlier.changesCn).toBe(true);
        expect(Math.abs(shell.cut.later.cn - 1)).toBeLessThanOrEqual(0.05);
    });

    it('integrates 4π ρ r² max(g, 0) over exactly the window, the edges interpolated', () => {
        // g = 1 on [2, 4] (a step), so CN = 4π ρ (b³ − a³)/3 for any [a, b] inside.
        const x = grid(1000, 0.01);
        const y = x.map((r) => (r >= 2 && r <= 4 ? 1 : 0) + gaussian(r, 3, 0.05, 20) - gaussian(r, 3, 0.05, 20));
        const catalogue = catalogueOf({ label: 'Ta-Se', x, y });
        const [a, b] = [2.333, 3.777];
        const cn = catalogue.window('Ta', 'Se', [a, b]).cn;
        expect(cn).toBeCloseTo(4 * Math.PI * RHO * (b ** 3 - a ** 3) / 3, 3);
        // A negative g counts as no atoms.
        const negative = catalogueOf({ label: 'Ta-Se', x, y: y.map((value, i) => (x[i] > 3 && x[i] < 3.5 ? -1 : value)) });
        expect(negative.window('Ta', 'Se', [3.05, 3.45]).cn).toBe(0);
    });

    it('reads the mean distance and spread from r²g over the quoted window, not from the g peak', () => {
        const shell = catalogueOf(shellSeries('Ta-Se', [[2.6, 0.12, 6]], RHO)).shell('Ta', 'Se');
        expect(shell.peaks[0].r).toBeLessThan(2.6);
        expect(shell.meanR).toBeCloseTo(2.6, 2);
        expect(shell.sdR).toBeCloseTo(0.12, 2);
        // The window starts one grid step before the onset (g at 1 % of its
        // peak): the tail below it, ~0.2 % of the shell, is outside.
        expect(shell.cn).toBeGreaterThan(5.98);
        expect(shell.cn).toBeLessThan(6);
    });
});

describe('one search range: peaks below 8 Å, the closing minimum to max(8, r₁ + 3) Å', () => {
    it('closes a shell beyond 8 Å when its first peak lies above 5 Å, whatever range a caller would pick', () => {
        // A broad shell at 6.8 Å (σ 0.45 Å): g falls to ~0 near 8.4 Å.
        const shell = catalogueOf(shellSeries('Ga-Ga', [[6.8, 0.45, 12]], RHO, { rMax: 12, tail: 10.5 }), { structure: box({ Ga: 50 }) })
            .shell('Ga', 'Ga');
        expect(PEAK_SEARCH_MAX_A).toBe(8);
        expect(BOUNDARY_REACH_A).toBe(3);
        expect(shell.window[1]).toBeGreaterThan(PEAK_SEARCH_MAX_A);
        expect(shell.window[1]).toBeLessThanOrEqual(6.8 + BOUNDARY_REACH_A);
        expect(Math.abs(shell.cn - 12)).toBeLessThanOrEqual(0.1);
    });

    it('reads no shell for a pair whose first peak lies beyond 8 Å', () => {
        const catalogue = catalogueOf(shellSeries('Ga-Ga', [[8.6, 0.1, 12]], RHO, { rMax: 12, tail: 11 }));
        expect(catalogue.shell('Ga', 'Ga')).toBeNull();
        expect(catalogue.pairs).toEqual([]);
        expect(catalogue.window('Ga', 'Ga', [8, 9])).toBeNull();
    });
});

describe('the constraint-limited rule: g above 5 % of the peak at the cut (the user\'s rule, 2026-10-24)', () => {
    // A Ga–Se peak at 2.42 Å (σ 0.1 Å) cut at its closest-approach
    // constraint: g is 0 below it, and just above it the share of the first
    // peak that the cut's position sets.
    const settings = (limit) => ({ atoms: ['Ga', 'Se'], minimumDistancesA: [3.0, limit, 3.0] });
    const cutAt = (limit) => {
        const x = grid(1000, 0.01);
        const series = { label: 'Ga-Se', x, y: x.map((r) => (r <= limit ? 0 : gaussian(r, 2.42, 0.1, 8))) };
        return catalogueOf(series, { runSettings: settings(limit) }).shell('Ga', 'Se').constraint;
    };

    it('its own constant, 5 %; the shell separation keeps 2 %', () => {
        expect(CONSTRAINT_LIMIT_FRACTION).toBe(0.05);
        expect(SEPARATION_G_FRACTION).toBe(0.02);
    });

    it('3 % of the peak at the cut is not constraint-limited (it was under the old 2 %)', () => {
        const constraint = cutAt(2.13);
        expect(Math.round(constraint.gFraction * 100)).toBe(3);
        expect(constraint.gFraction).toBeGreaterThan(SEPARATION_G_FRACTION);
        expect(constraint).toMatchObject({ limitA: 2.13, limited: false });
    });

    it('9 % of the peak at the cut is constraint-limited, like the Demo\'s Ga–Se (8.7 % at 2.2 Å)', () => {
        const constraint = cutAt(2.18);
        expect(Math.round(constraint.gFraction * 100)).toBe(9);
        expect(constraint).toMatchObject({ limitA: 2.18, limited: true });
    });
});

describe('the shell separation (SEPARATION_G_FRACTION, 2 %) and the closest-approach pile-up', () => {
    it('a shell is separated when g at its boundary is at most 2 % of its peak', () => {
        expect(SEPARATION_G_FRACTION).toBe(0.02);
        const x = grid(700, 0.01);
        const filled = (level) => catalogueOf({
            label: 'Ta-Ta', x, y: x.map((r) => gaussian(r, 3.0, 0.1, 5) + gaussian(r, 4.3, 0.1, 4) + (r > 3.2 && r < 4.1 ? level : 0))
        }, { structure: box({ Ta: 50 }) }).shell('Ta', 'Ta');
        const gap = filled(0.35);
        expect(gap.troughRatio).toBeGreaterThan(SEPARATION_G_FRACTION);
        expect(gap.separated).toBe(false);
        expect(gap.cut.earlier.cn).toBeLessThan(gap.cn);
        expect(gap.cut.later.cn).toBeGreaterThan(gap.cn);
        const clean = filled(0);
        expect(clean.troughRatio).toBeLessThanOrEqual(SEPARATION_G_FRACTION);
        expect(clean.separated).toBe(true);
    });

    // The title keeps the 2 % of the shell-catalogue decision that cites it
    // (2026-10-19); since 2026-10-24 the constraint has its own fraction, 5 %
    // (CONSTRAINT_LIMIT_FRACTION, the describe above). Both cases hold under it.
    it('a constraint limits a pair when raw g just above it is above 2 % of the peak — the FWHM edge clause is gone', () => {
        const settings = (limit) => ({ atoms: ['Ga', 'Se'], minimumDistancesA: [3.0, limit, 3.0] });
        // Cut at the constraint: g is 0 below 2.2 Å, ~9 % of the peak just above it.
        const x = grid(1000, 0.01);
        const cutSeries = { label: 'Ga-Se', x, y: x.map((r) => (r < 2.2 ? 0 : gaussian(r, 2.42, 0.1, 8))) };
        const cut = catalogueOf(cutSeries, { runSettings: settings(2.2) }).shell('Ga', 'Se').constraint;
        expect(cut).toMatchObject({ limitA: 2.2, limited: true });
        expect(cut.gFraction).toBeGreaterThan(0.05);
        // A sharp peak on a fine grid whose half maximum lies within 0.05 Å
        // above the constraint, while g just above it is ~1 % of the peak:
        // the old r − FWHM/2 ≤ d_min + 0.05 Å clause called it cut; it is not.
        const fine = grid(2000, 0.005);
        const sharp = { label: 'Ga-Se', x: fine, y: fine.map((r) => gaussian(r, 2.5, 0.02, 8)) };
        const [peak] = extractPeaks(sharp.x, sharp.y, { rMax: 8 });
        expect(peak.r - peak.fwhm / 2).toBeLessThanOrEqual(2.43 + 0.05);
        const clear = catalogueOf(sharp, { runSettings: settings(2.43) }).shell('Ga', 'Se').constraint;
        expect(clear.gFraction).toBeLessThan(SEPARATION_G_FRACTION);
        expect(clear.gFraction).toBeLessThan(CONSTRAINT_LIMIT_FRACTION);
        expect(clear.limited).toBe(false);
        // No constraint for the pair: none.
        expect(catalogueOf(cutSeries).shell('Ga', 'Se').constraint).toBeNull();
    });

    it('labels the run settings\' minimum distances in RMCProfile\'s pair order and matches them in either order', () => {
        expect(labelledMinimumDistances({ atoms: ['Ga', 'Se', 'Ta'], minimumDistancesA: [3, 2.38, 2.8, 3, 2.3, 3] }))
            .toEqual({ 'Ga-Ga': 3, 'Ga-Se': 2.38, 'Ga-Ta': 2.8, 'Se-Se': 3, 'Se-Ta': 2.3, 'Ta-Ta': 3 });
        // Counts that do not line up: the plain list, which labels no pair.
        expect(labelledMinimumDistances({ atoms: ['Ga'], minimumDistancesA: [2.2, 2.3] })).toEqual([2.2, 2.3]);
        expect(labelledMinimumDistances(null)).toBeNull();
        expect(labelledMinimumDistances({ minimumDistancesA: [] })).toBeNull();
        const x = grid(600, 0.01);
        const series = { label: 'Ta-Se', x, y: x.map((r) => gaussian(r, 2.6, 0.06, 8)) };
        const runSettings = { atoms: ['Ga', 'Se', 'Ta'], minimumDistancesA: [3, 2.38, 2.8, 3, 2.3, 3] };
        expect(catalogueOf(series, { runSettings }).shell('Ta', 'Se').constraint.limitA).toBe(2.3);
        expect(catalogueOf(series, { runSettings: { ...runSettings, atoms: ['Ga'] } }).shell('Ta', 'Se').constraint).toBeNull();
    });
});

describe('a split shell cut at a deep trough (the review\'s L1)', () => {
    it('joins a Jahn–Teller octahedron (4 short + 2 long, deep trough) to CN 6; the first peak alone ends in the trough', () => {
        const series = shellSeries('Cu-F', [[1.93, 0.05, 4], [2.25, 0.05, 2]], RHO, { step: 0.005 });
        const shell = catalogueOf(series, { structure: box({ Cu: 50, F: 50 }) }).shell('Cu', 'F');
        expect(shell.peaksInShell).toEqual([expect.closeTo(1.93, 1), expect.closeTo(2.25, 1)]);
        expect(shell.joined).toBe(true);
        expect(shell.cn).toBeCloseTo(6, 1);
        expect(shell.ideal).toBe(6);
        expect(shell.window[1]).toBeGreaterThan(2.35);
        expect(shell.firstPeak.cn).toBeCloseTo(4, 1);
        expect(shell.firstPeak.window[1]).toBeGreaterThan(2.0);
        expect(shell.firstPeak.window[1]).toBeLessThan(2.2);
        expect(shell.firstPeak.troughRatio).toBeLessThan(SHELL_TROUGH_FRACTION);
        // Without a density there is no CN to judge by: the old cut stands.
        const blind = catalogueOf(series, { structure: null }).shell('Cu', 'F');
        expect(blind.joined).toBe(false);
        expect(blind.cn).toBeNull();
        expect(blind.window[1]).toBeLessThan(2.2);
    });

    it('never joins a peak beyond the reach, nor one whose summed CN is no polyhedron', () => {
        const far = catalogueOf(shellSeries('Cu-F', [[2.0, 0.05, 6], [2.0 + SPLIT_SHELL_REACH_A + 0.1, 0.05, 2]], RHO, { step: 0.005 }),
            { structure: box({ Cu: 50, F: 50 }) }).shell('Cu', 'F');
        expect(far.cn).toBeCloseTo(6, 1);
        expect(far.joined).toBe(false);
        const odd = catalogueOf(shellSeries('Cu-F', [[2.0, 0.05, 4], [2.35, 0.05, 3]], RHO, { step: 0.005 }),
            { structure: box({ Cu: 50, F: 50 }) }).shell('Cu', 'F');
        expect(odd.joined).toBe(false);
        expect(odd.cn).toBeCloseTo(4, 1);
    });

    it('never joins a split M3/M4 cluster shell (1 + 1, 2 + 1): CN 2 and 3 are no ideal count (F5, a known limit)', () => {
        const triangle = catalogueOf(shellSeries('Mo-Mo', [[2.7, 0.04, 1], [2.95, 0.04, 1]], RHO, { step: 0.005 }),
            { structure: box({ Mo: 50, Se: 50 }) }).shell('Mo', 'Mo');
        expect(triangle.peaksInShell).toEqual([expect.closeTo(2.7, 1)]);
        expect(triangle.joined).toBe(false);
        expect(triangle.cn).toBeCloseTo(1, 1);
        const tetrahedron = catalogueOf(shellSeries('Ta-Ta', [[2.9, 0.04, 2], [3.15, 0.04, 1]], RHO, { step: 0.005 }),
            { structure: box({ Ta: 50, Se: 50 }) }).shell('Ta', 'Ta');
        expect(tetrahedron.joined).toBe(false);
        expect(tetrahedron.cn).toBeCloseTo(2, 1);
        expect(IDEAL_COORDINATION).not.toContain(2);
        expect(IDEAL_COORDINATION).not.toContain(3);
    });

    it('takes the furthest peak that completes a polyhedron (2 + 2 + 2), and counts either direction', () => {
        const ortho = catalogueOf(shellSeries('Cu-F', [[1.95, 0.03, 2], [2.12, 0.03, 2], [2.32, 0.03, 2]], RHO, { step: 0.005 }),
            { structure: box({ Cu: 50, F: 50 }) }).shell('Cu', 'F');
        expect(ortho.cn).toBeCloseTo(6, 1);
        expect(ortho.firstPeak.cn).toBeCloseTo(2, 1);
        // 1 + 0.5 F per Cu is no polyhedron, but Cu is four times as dense: 6 Cu per F.
        const reverse = shellSeries('Cu-F', [[1.93, 0.05, 1], [2.25, 0.05, 0.5]], RHO, { step: 0.005 });
        expect(catalogueOf(reverse, { structure: box({ Cu: 50, F: 50 }) }).shell('Cu', 'F').joined).toBe(false);
        const joined = catalogueOf(reverse, { structure: box({ Cu: 200, F: 50 }) }).shell('Cu', 'F');
        expect(joined.joined).toBe(true);
        expect(joined.cn).toBeCloseTo(1.5, 1);
        expect(joined.reverseCn).toBeCloseTo(6, 1);
        expect(joined.ideal).toBeNull();
        expect(joined.firstPeak.cn).toBeCloseTo(1, 1);
        expect(IDEAL_COORDINATION).toEqual([4, 6, 8, 12]);
    });
});

describe('a joined shell\'s split by number (plan item 4.2)', () => {
    // A Jahn–Teller CuF₆: r²g holds exactly 4 F at 1.93 Å and 2 at 2.25 Å
    // (σ 0.05 Å each), 50 Cu in the box. The tolerance on each mean and on
    // long − short, 0.002 Å: the means are moments of r²g, not grid
    // positions, so the grid moves them only through the quadrature (< 1e-4
    // Å at h/σ ≤ 0.4) and the window edges, and the tails cut off beyond the
    // inter-peak minimum and below the shell's onset lie 3σ out (~0.1 % of a
    // part), moving a part's mean by < 5e-4 Å; 0.002 Å is four times that,
    // and a fifth of the 0.01 Å every window is rounded to.
    const JT = [[1.93, 0.05, 4], [2.25, 0.05, 2]];
    const CU_F = box({ Cu: 50, F: 50 });
    const TOLERANCE_A = 0.002;

    it.each([0.005, 0.01, 0.02])('a 4 + 2 Jahn–Teller octahedron (grid %s Å): 4 short at 1.93 + 2 long at 2.25 Å, long − short 0.32 Å, shares 2/3 and 1/3, the parts summing to the shell', (step) => {
        const shell = catalogueOf(shellSeries('Cu-F', JT, RHO, { step }), { structure: CU_F }).shell('Cu', 'F');
        expect(shell.joined).toBe(true);
        const { split } = shell;
        expect(split.direction).toBe('forward');
        expect(split.parts).toHaveLength(2);
        const [short, long] = split.parts;
        // The short part is the first peak alone; the long part runs from its
        // end to the shell's boundary.
        expect(short.window).toEqual(shell.firstPeak.window);
        expect(short.cn).toBe(shell.firstPeak.cn);
        expect(long.window).toEqual([shell.firstPeak.window[1], shell.window[1]]);
        // Each CN over exactly its window: the parts sum to the shell's.
        expect(short.cn + long.cn).toBeCloseTo(shell.cn, 12);
        expect(short.share + long.share).toBeCloseTo(1, 12);
        expect(Math.abs(short.cn - 4)).toBeLessThanOrEqual(0.02);
        expect(Math.abs(long.cn - 2)).toBeLessThanOrEqual(0.02);
        expect(Math.abs(short.share - 2 / 3)).toBeLessThanOrEqual(0.005);
        expect(Math.abs(long.share - 1 / 3)).toBeLessThanOrEqual(0.005);
        expect(Math.abs(short.meanR - 1.93)).toBeLessThanOrEqual(TOLERANCE_A);
        expect(Math.abs(long.meanR - 2.25)).toBeLessThanOrEqual(TOLERANCE_A);
        expect(Math.abs(split.longMinusShort - 0.32)).toBeLessThanOrEqual(TOLERANCE_A);
        expect(short.sdR).toBeCloseTo(0.05, 2);
        expect(long.sdR).toBeCloseTo(0.05, 2);
        // The shell's mean distance is the parts', weighted by their shares.
        expect(short.share * short.meanR + long.share * long.meanR).toBeCloseTo(shell.meanR, 10);
        // Its sampling error: each part's spread over √ its bonds in the box
        // (50 Cu × 4 and × 2), the bonds taken as independent.
        expect(split.longMinusShortSe).toBeCloseTo(Math.hypot(short.sdR / Math.sqrt(50 * short.cn), long.sdR / Math.sqrt(50 * long.cn)), 12);
        expect(split.longMinusShortSe).toBeCloseTo(0.05 * Math.sqrt(1 / 200 + 1 / 100), 3);
    });

    it('2 + 2 + 2 (an orthorhombic octahedron): every part, a third of the shell each, long − short the longest minus the shortest', () => {
        // Short, medium and long pairs, as a Jahn–Teller octahedron with both
        // the Q2 and the Q3 mode; lumping the medium pair with either end
        // would read it as a 4 + 2 elongation or a 2 + 4 compression.
        const shell = catalogueOf(shellSeries('Cu-F', [[1.95, 0.03, 2], [2.12, 0.03, 2], [2.32, 0.03, 2]], RHO, { step: 0.005 }), { structure: CU_F })
            .shell('Cu', 'F');
        expect(shell.joined).toBe(true);
        const { parts, longMinusShort } = shell.split;
        expect(parts).toHaveLength(3);
        expect(parts[0].window).toEqual(shell.firstPeak.window);
        // Contiguous, from the shell's start to its boundary.
        expect(parts.map((part) => part.window[0])).toEqual([shell.window[0], parts[0].window[1], parts[1].window[1]]);
        expect(parts[2].window[1]).toBe(shell.window[1]);
        parts.forEach((part, k) => {
            expect(Math.abs(part.cn - 2), `part ${k}`).toBeLessThanOrEqual(0.03);
            expect(Math.abs(part.share - 1 / 3), `part ${k}`).toBeLessThanOrEqual(0.005);
            expect(Math.abs(part.meanR - [1.95, 2.12, 2.32][k]), `part ${k}`).toBeLessThanOrEqual(TOLERANCE_A);
        });
        expect(parts.reduce((sum, part) => sum + part.cn, 0)).toBeCloseTo(shell.cn, 12);
        expect(Math.abs(longMinusShort - 0.37)).toBeLessThanOrEqual(TOLERANCE_A);
    });

    it('says the counts in the direction the join decided: reverse when only the centre\'s count is ideal', () => {
        // 1 + 0.5 F per Cu is no polyhedron; with Cu four times as dense,
        // 4 + 2 Cu per F is.
        const series = shellSeries('Cu-F', [[1.93, 0.05, 1], [2.25, 0.05, 0.5]], RHO, { step: 0.005 });
        const structure = box({ Cu: 200, F: 50 });
        const shell = catalogueOf(series, { structure }).shell('Cu', 'F');
        expect(shell.joined).toBe(true);
        expect(shell.split.direction).toBe('reverse');
        expect(shell.split.parts.map((part) => part.reverseCn)).toEqual([expect.closeTo(4, 1), expect.closeTo(2, 1)]);
        expect(shell.split.parts.map((part) => part.share)).toEqual([expect.closeTo(2 / 3, 2), expect.closeTo(1 / 3, 2)]);
        expect(Math.abs(shell.split.longMinusShort - 0.32)).toBeLessThanOrEqual(TOLERANCE_A);
        // The same shell from the F side, whose own count (6 Cu per F) is
        // ideal: forward, the same split, the same bonds and so the same error.
        const fromF = catalogueOf(series, { structure }).shell('F', 'Cu');
        expect(fromF.split.direction).toBe('forward');
        expect(fromF.split.parts.map((part) => part.cn)).toEqual(shell.split.parts.map((part) => expect.closeTo(part.reverseCn, 9)));
        expect(fromF.split.longMinusShort).toBe(shell.split.longMinusShort);
        expect(fromF.split.longMinusShortSe).toBeCloseTo(shell.split.longMinusShortSe, 12);
    });

    it('no split for a shell that is not joined: one peak, two overlapping peaks with no deep trough, no density to join by', () => {
        expect(catalogueOf(shellSeries('Ta-Se', [[2.6, 0.08, 6]], RHO)).shell('Ta', 'Se').split).toBeNull();
        // 3 + 3 at 2.5 and 2.75 Å, the trough at ~40 % of the first peak:
        // the first peak alone is read, but the minimum cuts through both
        // bond groups, so no split by number.
        const overlapping = catalogueOf(shellSeries('Ta-Se', [[2.5, 0.07, 3], [2.75, 0.07, 3]], RHO)).shell('Ta', 'Se');
        expect(overlapping.firstPeak).not.toBeNull();
        expect(overlapping.joined).toBe(false);
        expect(overlapping.split).toBeNull();
        expect(catalogueOf(shellSeries('Cu-F', JT, RHO, { step: 0.005 }), { structure: null }).shell('Cu', 'F').split).toBeNull();
    });

    it('the cut at the minimum: a small positive bias as the parts overlap, below 1 % of the split at a 0.07 Å bond spread and 2 % at the edge of the join', () => {
        // Noise-free r²g, so only the cut acts: each part's tail beyond the
        // minimum moves its mean toward the other part's, the long part's
        // more (its tail is the larger share of it).
        for (const step of [0.005, 0.01, 0.02]) {
            const shell = catalogueOf(shellSeries('Cu-F', [[1.93, 0.07, 4], [2.25, 0.07, 2]], RHO, { step }), { structure: CU_F }).shell('Cu', 'F');
            const bias = shell.split.longMinusShort - 0.32;
            expect(bias, `grid ${step} Å`).toBeGreaterThan(0);
            expect(bias, `grid ${step} Å`).toBeLessThan(0.0032);
        }
        // σ 0.085 Å: the trough at 20 % of the first peak, the deepest the join takes.
        const edge = catalogueOf(shellSeries('Cu-F', [[1.93, 0.085, 4], [2.25, 0.085, 2]], RHO, { step: 0.005 }), { structure: CU_F }).shell('Cu', 'F');
        expect(edge.firstPeak.troughRatio).toBeGreaterThan(0.19);
        expect(edge.split.longMinusShort - 0.32).toBeGreaterThan(0);
        expect(edge.split.longMinusShort - 0.32).toBeLessThan(0.0064);
    });
});

describe('the split\'s sampling error: the configuration\'s own, the bonds independent', () => {
    // Seeded draws (mulberry32, Box–Muller), so the test is deterministic.
    const generator = (seed) => () => {
        seed = (seed + 0x6d2b79f5) | 0;
        let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
    const normal = (random) => {
        let u = 0;
        while (u === 0) u = random();
        return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * random());
    };
    // `centres` octahedra of independent bond lengths, 4 ~ N(1.93, σ) and 2 ~
    // N(2.25, σ), histogrammed as the dev cases' partials are
    // (validation/cases.js partialGr: 0.01 Å bins, each distance spread as a
    // 0.015 Å Gaussian), g normalised so that CN = 4π ρ ∫ r² g dr.
    const configuration = (seed, { centres = 216, sigma = 0.05, dr = 0.01, kernel = 0.015 } = {}) => {
        const random = generator(seed);
        const bins = Math.round(6 / dr);
        const counts = new Float64Array(bins);
        const spread = Math.ceil((3 * kernel) / dr);
        const weights = Array.from({ length: 2 * spread + 1 }, (_, k) => Math.exp(-0.5 * (((k - spread) * dr) / kernel) ** 2));
        const total = weights.reduce((sum, w) => sum + w, 0);
        for (let c = 0; c < centres; c += 1) {
            for (const [r0, n] of [[1.93, 4], [2.25, 2]]) {
                for (let j = 0; j < n; j += 1) {
                    const bin = Math.floor((r0 + sigma * normal(random)) / dr);
                    weights.forEach((w, k) => { if (bin + k - spread >= 0 && bin + k - spread < bins) counts[bin + k - spread] += w / total; });
                }
            }
        }
        const L = 8 * Math.cbrt(centres);
        const rho = (6 * centres) / L ** 3;
        const x = Array.from({ length: bins }, (_, k) => (k + 0.5) * dr);
        const y = x.map((r, k) => counts[k] / (centres * 4 * Math.PI * r * r * dr * rho));
        return { series: { label: 'Cu-F', x, y }, structure: box({ Cu: centres, F: 6 * centres }, L) };
    };

    it('(long − short − 0.32) / its sampling error is standard normal over 300 configurations of 216 octahedra', () => {
        const z = [];
        for (let seed = 1; seed <= 300; seed += 1) {
            const { series, structure } = configuration(seed);
            const { split } = catalogueOf(series, { structure }).shell('Cu', 'F');
            expect(split?.parts, `seed ${seed}`).toHaveLength(2);
            z.push((split.longMinusShort - 0.32) / split.longMinusShortSe);
        }
        const mean = z.reduce((sum, value) => sum + value, 0) / z.length;
        const sd = Math.sqrt(z.reduce((sum, value) => sum + (value - mean) ** 2, 0) / (z.length - 1));
        const covered = z.filter((value) => Math.abs(value) <= 1.96).length / z.length;
        // Each bound ~3.5 standard errors of its statistic at 300 draws.
        expect(Math.abs(mean)).toBeLessThan(0.2);
        expect(sd).toBeGreaterThan(0.85);
        expect(sd).toBeLessThan(1.15);
        expect(covered).toBeGreaterThan(0.91);
        expect(covered).toBeLessThan(0.99);
    });
});

describe('the CN match and the two cut readings', () => {
    it('cnMatches: within half an atom, inclusive; the ideal counts are built on it', () => {
        expect(CN_MATCH_TOLERANCE).toBe(0.5);
        expect(cnMatches(5.5, 6)).toBe(true);
        expect(cnMatches(6.5, 6)).toBe(true);
        expect(cnMatches(5.49, 6)).toBe(false);
        expect(cnMatches(NaN, 6)).toBe(false);
        expect(cnMatches(null, 0)).toBe(false);
        expect(isIdealCoordination(3.6)).toBe(true);
        expect(isIdealCoordination(5)).toBe(false);
        expect(isIdealCoordination(12.4)).toBe(true);
    });

    it('cnChanges: the Workbench\'s margin, max(0.5, 15 % of the shell\'s CN)', () => {
        expect([CN_CHANGE_ABS, CN_CHANGE_REL]).toEqual([0.5, 0.15]);
        expect(cnChanges(5.1, 6)).toBe(true); // 0.9 = 15 % of 6
        expect(cnChanges(5.2, 6)).toBe(false);
        expect(cnChanges(2.5, 3)).toBe(true); // 0.5 > 15 % of 3
        expect(cnChanges(2.6, 3)).toBe(false);
        expect(cnChanges(null, 3)).toBe(false);
    });

    it('window(): what a window holds, and its CN-margin reading', () => {
        const catalogue = catalogueOf(shellSeries('Ta-Se', [[2.52, 0.06, 3], [2.74, 0.1, 3]], RHO));
        const shell = catalogue.shell('Ta', 'Se');
        const [lo, hi] = shell.window;
        expect(catalogue.window('Ta', 'Se', [lo, hi])).toMatchObject({ holds: 'whole', cnMargin: { cut: false }, shell });
        expect(catalogue.window('Ta', 'Se', [lo, hi])).not.toHaveProperty('counts');
        expect(catalogue.window('Ta', 'Se', [lo - 0.3, hi])).toMatchObject({ holds: 'whole', cnMargin: { cut: false } });
        expect(catalogue.window('Ta', 'Se', [lo, hi + 0.3])).toMatchObject({ holds: 'more', cnMargin: { cut: false } });
        expect(catalogue.window('Ta', 'Se', shell.firstPeak.window)).toMatchObject({ holds: 'part', cnMargin: { cut: true } });
        // Part of the shell but within the margin: not a cut.
        expect(catalogue.window('Ta', 'Se', [lo + 0.1, hi])).toMatchObject({ holds: 'part', cnMargin: { cut: false } });
        expect(catalogue.window('Ta', 'Se', [hi + 0.5, hi + 1])).toMatchObject({ holds: 'apart', cnMargin: { cut: false } });
        expect(catalogue.window('Ta', 'Se', [2.6, 2.5])).toBeNull();
        expect(cnMarginCut(1, shell, [hi + 0.5, hi + 1])).toBe(false);
        expect(cnMarginCut(shell.firstPeak.cn, shell, shell.firstPeak.window)).toBe(true);
    });

    // Since 2026-10-19 compute_bond_angles' rule only where there is no shell
    // to compare a window with (windowCut's fallback, below). The title is the
    // one the shell-catalogue decision record names.
    it('countsCut: compute_bond_angles\' rule, pinned (60 % modal share, more than half unbonded)', () => {
        expect([MODAL_COORDINATION_SHARE, UNBONDED_APEX_FRACTION]).toEqual([0.6, 0.5]);
        const summary = (modalShare, unbondedShare) => ({ mean: 3, modal: 3, modalShare, unbondedShare });
        expect(countsCut(summary(0.59, 0))).toEqual({ cut: true, by: 'count-spread' });
        expect(countsCut(summary(0.61, 0))).toEqual({ cut: false, by: null });
        expect(countsCut(summary(0.6, 0))).toEqual({ cut: false, by: null });
        expect(countsCut(summary(0.51, 0.51))).toEqual({ cut: true, by: 'unbonded' });
        expect(countsCut(summary(0.5, 0.5))).toEqual({ cut: true, by: 'count-spread' });
        expect(countsCut(summary(0.7, 0.5))).toEqual({ cut: false, by: null });
        expect(countsCut({ mean: null, modal: null, modalShare: null, unbondedShare: null })).toEqual({ cut: false, by: null });
    });

    // The counts the bond-angle engine would give: a share with the modal
    // count, a share with none.
    const counts = (modalShare, unbondedShare = 0) => ({ mean: 3, modal: 3, modalShare, unbondedShare });

    it('windowCut, the one rule: where a window meets a shell with a CN the CN margin decides, whatever the counts', () => {
        // 3 + 3 partners at 2.52 and 2.74 Å: one split shell of 6.
        const catalogue = catalogueOf(shellSeries('Ta-Se', [[2.52, 0.06, 3], [2.74, 0.1, 3]], RHO));
        const shell = catalogue.shell('Ta', 'Se');
        const [lo, hi] = shell.window;
        const read = (window, given = null) => windowCut(catalogue.window('Ta', 'Se', window), given);
        // The whole shell: not cut, even with counts that spread (a shell
        // running into the next) — the counts no longer decide here.
        expect(read([lo, hi], counts(0.4))).toEqual({ cut: false, by: null, rule: 'cn-margin', whole: true });
        // The first peak alone (3 of 6): cut, even with uniform counts (a
        // sub-shell every centre has whole); and without counts (the Workbench).
        expect(read(shell.firstPeak.window, counts(0.95))).toEqual({ cut: true, by: 'cn-margin', rule: 'cn-margin', whole: false });
        expect(read(shell.firstPeak.window)).toEqual({ cut: true, by: 'cn-margin', rule: 'cn-margin', whole: false });
        // The top 0.2 Å early, within the margin: not cut, and its count is the shell's.
        expect(shell.cut.earlier.changesCn).toBe(false);
        expect(read(shell.cut.earlier.window, counts(0.4))).toEqual({ cut: false, by: null, rule: 'cn-margin', whole: true });
        // Beyond the shell by more than the margin (g = 1 from 5.5 Å on: the
        // next atoms too): not cut, and not the shell's count.
        const wide = [lo, 6.5];
        expect(catalogue.window('Ta', 'Se', wide).cn - shell.cn).toBeGreaterThan(Math.max(CN_CHANGE_ABS, CN_CHANGE_REL * shell.cn));
        expect(read(wide, counts(0.95))).toEqual({ cut: false, by: null, rule: 'cn-margin', whole: false });
    });

    it('windowCut: with no shell to compare, the counts decide (the fallback); a window apart from the shell never holds it', () => {
        const catalogue = catalogueOf(shellSeries('Ta-Se', [[2.52, 0.06, 3], [2.74, 0.1, 3]], RHO));
        const [, hi] = catalogue.shell('Ta', 'Se').window;
        // Apart from the first shell (the module reads first shells only).
        const apart = catalogue.window('Ta', 'Se', [hi + 0.5, hi + 1]);
        expect(apart.holds).toBe('apart');
        expect(windowCut(apart, counts(0.4))).toEqual({ cut: true, by: 'count-spread', rule: 'counts', whole: false });
        expect(windowCut(apart, counts(0.3, 0.7))).toEqual({ cut: true, by: 'unbonded', rule: 'counts', whole: false });
        expect(windowCut(apart, counts(0.95))).toEqual({ cut: false, by: null, rule: 'counts', whole: false });
        // No closed shell for the pair, or no partials: as before, the counts' verdict is the whole story.
        expect(windowCut(null, counts(0.47))).toEqual({ cut: true, by: 'count-spread', rule: 'counts', whole: false });
        expect(windowCut(null, counts(0.95))).toEqual({ cut: false, by: null, rule: 'counts', whole: true });
        // No density, so no CN to compare: the counts.
        const noDensity = catalogueOf(shellSeries('Ta-Se', [[2.52, 0.06, 3], [2.74, 0.1, 3]], RHO), { structure: null });
        const reading = noDensity.window('Ta', 'Se', noDensity.shell('Ta', 'Se').firstPeak.window);
        expect(reading.cn).toBeNull();
        expect(windowCut(reading, counts(0.95))).toEqual({ cut: false, by: null, rule: 'counts', whole: true });
        // Nothing to decide with (the Workbench, before a call): not cut.
        expect(windowCut(apart)).toEqual({ cut: false, by: null, rule: null, whole: false });
        expect(windowCut(null)).toEqual({ cut: false, by: null, rule: null, whole: false });
    });
});

describe('the average structure\'s first shell: find_gr_peaks\' cross-check (the user\'s rule, 2026-10-24)', () => {
    // A 10 Å cubic cell (supercell 2×2×2: the conventional cell is the
    // lattice over it): a centre M at the cell's centre with four X at 1.93 Å
    // in a square and a fifth at 2.40 Å above it, a 4 + 1 pyramid, and a
    // second X site with no M within the reach.
    const a = 10;
    const site = (el, x, y, z) => ({ el, frac: [x / a, y / a, z / a] });
    const structure = {
        latticeVectors: [[2 * a, 0, 0], [0, 2 * a, 0], [0, 0, 2 * a]],
        supercell: [2, 2, 2],
        basis: [
            site('M', 5, 5, 5),
            site('X', 6.93, 5, 5), site('X', 3.07, 5, 5), site('X', 5, 6.93, 5), site('X', 5, 3.07, 5), site('X', 5, 5, 7.4),
            site('X', 0.5, 0.5, 0.5)
        ]
    };

    it('reads every bond between site means from the shortest to the join\'s reach beyond it, either way round', () => {
        expect(AVERAGE_SHELL_REACH_A).toBe(SPLIT_SHELL_REACH_A);
        expect(averageStructureBonds(structure, 'M', 'X')).toEqual([1.93, 2.4]);
        expect(averageStructureBonds(structure, 'X', 'M')).toEqual([1.93, 2.4]);
        // Not each site's nearest partner: the lone X's nearest M (7.79 Å) is no bond.
        expect(Math.max(...averageStructureBonds(structure, 'X', 'M'))).toBeLessThan(1.93 + AVERAGE_SHELL_REACH_A + 1e-9);
        // Periodic, through the cell metric: two M a cell apart along a are 10 Å apart.
        expect(averageStructureBonds(structure, 'M', 'M')).toEqual([10]);
        expect(averageStructureBonds({ ...structure, basis: [] }, 'M', 'X')).toEqual([]);
        expect(averageStructureBonds(structure, 'M', 'Y')).toEqual([]);
        expect(averageStructureBonds(null, 'M', 'X')).toEqual([]);
    });

    it('disagrees only when a bond lies outside the window by more than the first peak\'s FWHM', () => {
        const shell = { window: [1.7, 2.17], peaks: [{ r: 1.93, fwhm: 0.12 }] };
        expect(averageShellCheck(shell, [1.93])).toEqual({ bonds: [1.93, 1.93], outsideA: 0, widthA: 0.12, disagrees: false });
        // The pyramid's fifth bond, 0.23 Å past the closing minimum.
        expect(averageShellCheck(shell, [1.93, 2.4])).toMatchObject({ bonds: [1.93, 2.4], disagrees: true });
        expect(averageShellCheck(shell, [1.93, 2.4]).outsideA).toBeCloseTo(0.23, 12);
        // Just within one width beyond the window, and just past it.
        expect(averageShellCheck(shell, [2.28]).disagrees).toBe(false);
        expect(averageShellCheck(shell, [2.3]).disagrees).toBe(true);
        // Below the window's onset too.
        expect(averageShellCheck(shell, [1.5]).disagrees).toBe(true);
        expect(averageShellCheck({ window: null, peaks: [] }, [1.93])).toBeNull();
        expect(averageShellCheck(shell, [])).toBeNull();
    });
});

describe('the catalogue', () => {
    const run = () => ({
        plotFiles: partials([
            shellSeries('Ta-Se', [[2.52, 0.06, 3], [2.74, 0.1, 3]], RHO),
            shellSeries('Ta-Ta', [[3.0, 0.08, 3]], RHO),
            { label: 'total', x: [0, 1], y: [0, 1] }
        ]),
        structure: box({ Ta: 50, Se: 100 }),
        runSettings: null
    });

    it('reads one shell per labelled pair, in the file\'s order, and the reverse with the other density', () => {
        const props = run();
        const catalogue = shellCatalogue(props);
        expect(catalogue.elements).toEqual(['Se', 'Ta']);
        expect(catalogue.pairs.map((shell) => [shell.label, shell.centre, shell.partner])).toEqual([['Ta-Se', 'Ta', 'Se'], ['Ta-Ta', 'Ta', 'Ta']]);
        const forward = catalogue.shell('Ta', 'Se');
        const reverse = catalogue.shell('Se', 'Ta');
        expect(reverse.label).toBe('Ta-Se');
        expect(reverse.window).toEqual(forward.window);
        // Se is twice as dense as Ta: twice the Se around each Ta.
        expect(forward.cn).toBeCloseTo(2 * reverse.cn, 9);
        expect(forward.reverseCn).toBeCloseTo(reverse.cn, 9);
        expect(catalogue.shell('Se', 'Se')).toBeNull();
        expect(elementDensity(props.structure, 'Se')).toBeCloseTo(0.1, 12);
    });

    it('is memoised per partials data, structure and run settings, and its shells are frozen', () => {
        const props = run();
        const catalogue = shellCatalogue(props);
        expect(shellCatalogue({ ...props })).toBe(catalogue);
        expect(shellCatalogue({ ...props, structure: box({ Ta: 50, Se: 50 }) })).not.toBe(catalogue);
        expect(shellCatalogue({ ...props, runSettings: { atoms: ['Se', 'Ta'], minimumDistancesA: [2, 2.2, 2] } })).not.toBe(catalogue);
        const shell = catalogue.shell('Ta', 'Se');
        expect(Object.isFrozen(shell)).toBe(true);
        expect(Object.isFrozen(shell.window)).toBe(true);
        expect(Object.isFrozen(shell.firstPeak)).toBe(true);
        expect(() => { shell.window[1] = 9; }).toThrow(TypeError);
    });

    it('is null without a partials file of element pairs', () => {
        expect(shellCatalogue({})).toBeNull();
        expect(shellCatalogue({ plotFiles: [{ plotKind: 'npdf', plotData: { series: [{ label: 'x' }] } }] })).toBeNull();
        expect(shellCatalogue({ plotFiles: partials([{ label: 'total', x: [0, 1], y: [0, 1] }]) })).toBeNull();
    });
});
