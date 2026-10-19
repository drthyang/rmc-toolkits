// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Pair-correlation evidence for the LLM context: the nearest-neighbour distance
// each element pair has in the AVERAGE structure (basis sites in the
// conventional cell, over periodic images) next to the first peaks of that
// pair's partial-PDF g(r) and its first coordination shell, as the run's one
// shell reading gives them (shells.js). A peak displaced from — or split
// around — the average-structure distance is the short-range-correlation
// signal the average structure cannot express. Pure functions of plain data,
// like runContext.

import { shellCatalogue } from './shells';

// 3-sig-digit rounding, duplicated from runContext to avoid a circular import.
const roundSig = (value, digits = 3) => (
    Number.isFinite(value) && value !== 0 ? Number(value.toPrecision(digits)) : value
);

// Conventional unit cell (rows, Å) = supercell lattice / supercell dims.
// Deliberately duplicated from symmetryModel.js: the llm module must not
// import host modules (see src/llm/README.md import boundary).
const conventionalCell = (structure) => structure.latticeVectors.map(
    (row, i) => row.map((v) => v / Math.max(structure.supercell[i], 1))
);

const cartDistance = (df, A) => {
    let sum = 0;
    for (let k = 0; k < 3; k++) {
        const c = df[0] * A[0][k] + df[1] * A[1][k] + df[2] * A[2][k];
        sum += c * c;
    }
    return Math.sqrt(sum);
};

const pairKey = (a, b) => [a, b].sort().join('-');

/**
 * Nearest-neighbour distance (Å) for every element pair of the average
 * structure: minimum over basis-site pairs and ±1 periodic images. Returns a
 * Map "El1-El2" (elements alphabetical) → distance, or null without a basis.
 */
export const nearestNeighborDistances = (structure) => {
    if (!structure?.basis?.length || !structure?.latticeVectors || !structure?.supercell) return null;
    const A = conventionalCell(structure);
    const sites = structure.basis;
    const best = new Map();
    for (let i = 0; i < sites.length; i++) {
        for (let j = i; j < sites.length; j++) {
            const key = pairKey(sites[i].el, sites[j].el);
            let min = best.get(key) ?? Infinity;
            for (let ix = -1; ix <= 1; ix++) {
                for (let iy = -1; iy <= 1; iy++) {
                    for (let iz = -1; iz <= 1; iz++) {
                        if (i === j && ix === 0 && iy === 0 && iz === 0) continue;
                        const d = cartDistance([
                            sites[j].frac[0] + ix - sites[i].frac[0],
                            sites[j].frac[1] + iy - sites[i].frac[1],
                            sites[j].frac[2] + iz - sites[i].frac[2]
                        ], A);
                        // Ignore numerically-coincident sites (duplicated basis entries).
                        if (d > 1e-3 && d < min) min = d;
                    }
                }
            }
            if (Number.isFinite(min)) best.set(key, min);
        }
    }
    return best;
};

/**
 * The pair_correlations context section: one entry per partial-PDF series
 * (labels like "Nb-Nb"), carrying the average-structure nearest-neighbour
 * distance when a basis is available, and from the run's shell reading
 * (shells.js shellCatalogue) the first two g(r) peaks and the first shell. A
 * pair whose distribution is cut at its closest-approach constraint (the run
 * settings' MINIMUM_DISTANCES, any pair order), or whose first shell runs into
 * the next, carries a flag. Returns null when the run has no parsed partials.
 */
export const pairCorrelationsContext = (structure, plotFiles, { maxPairs = 10, runSettings = null } = {}) => {
    const partials = (plotFiles || []).find(
        (file) => file.plotKind === 'pdf_partials' && file.plotData?.series?.length
    );
    if (!partials) return null;
    const catalogue = shellCatalogue({ plotFiles, structure, runSettings });
    const expected = nearestNeighborDistances(structure);
    const pairs = [];
    for (const series of partials.plotData.series) {
        if (pairs.length >= maxPairs) break;
        const label = (series.label || '').trim();
        const match = /^([A-Z][a-z]?)\s*-\s*([A-Z][a-z]?)$/.exec(label);
        if (!match) continue;
        const entry = { pair: label };
        const expectedD = expected?.get(pairKey(match[1], match[2]));
        if (Number.isFinite(expectedD)) entry.avg_structure_d_A = roundSig(expectedD);
        const shell = catalogue?.shell(match[1], match[2]) ?? null;
        if (shell?.peaks.length) {
            entry.gr_peaks_A = shell.peaks.slice(0, 2).map((peak) => ({
                r: roundSig(peak.r),
                fwhm: roundSig(peak.fwhm, 2)
            }));
        }
        // The first shell's mean distance (from r²g, not the g peak) and, with
        // the cell, its coordination number: B atoms around each A.
        if (shell?.window && Number.isFinite(shell.meanR)) {
            entry.first_shell = { mean_d_A: roundSig(shell.meanR), to_A: roundSig(shell.window[1]) };
            if (Number.isFinite(shell.cn)) entry.first_shell[`cn_${match[2]}_per_${match[1]}`] = roundSig(shell.cn);
            // A shell of two or more resolved peaks: where its first peak ends
            // (the minimum between the first two), and what it holds.
            if (shell.firstPeak) {
                entry.first_shell.first_peak_to_A = roundSig(shell.firstPeak.window[1]);
                if (Number.isFinite(shell.firstPeak.cn)) entry.first_shell.first_peak_cn = roundSig(shell.firstPeak.cn);
            }
            // Distances filling the gap to the next shell.
            if (!shell.separated) entry.first_shell.gap_g_fraction = roundSig(shell.troughRatio, 2);
        }
        // Pairs pile up against the hard constraint: it shapes the peak.
        if (shell?.constraint?.limited) {
            entry.cut_at_min_distance = { at_A: shell.constraint.limitA, g_fraction: roundSig(shell.constraint.gFraction, 2) };
        }
        if (entry.avg_structure_d_A !== undefined || entry.gr_peaks_A) pairs.push(entry);
    }
    return pairs.length ? pairs : null;
};
