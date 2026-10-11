// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Pair-correlation evidence for the LLM context: the nearest-neighbour distance
// each element pair has in the AVERAGE structure (basis sites in the
// conventional cell, over periodic images) next to the first peaks of that
// pair's partial-PDF g(r). A peak displaced from — or split around — the
// average-structure distance is the short-range-correlation signal the average
// structure cannot express. Pure functions of plain data, like runContext.

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
 * First peaks of a g(r)-like curve below rMax: local maxima of a lightly
 * smoothed copy, above a floor relative to the window maximum, with an
 * approximate full width at half maximum from a half-height walk (overlapping
 * peaks make the width an overestimate — good enough for LLM evidence).
 * Returns [{ r, fwhm }] sorted by r, at most maxPeaks entries.
 */
export const extractPeaks = (x, y, { rMax = 6, maxPeaks = 2 } = {}) => {
    if (!Array.isArray(x) || !Array.isArray(y) || x.length !== y.length) return [];
    let end = x.findIndex((value) => value > rMax);
    if (end < 0) end = x.length;
    if (end < 5) return [];
    // 5-point moving average so grid noise does not mint peaks.
    const smooth = new Array(end);
    for (let i = 0; i < end; i++) {
        let sum = 0;
        let n = 0;
        for (let k = -2; k <= 2; k++) {
            const j = i + k;
            if (j >= 0 && j < end) { sum += y[j]; n += 1; }
        }
        smooth[i] = sum / n;
    }
    const windowMax = Math.max(...smooth);
    if (!(windowMax > 0)) return [];
    const floor = 0.15 * windowMax;
    const peaks = [];
    for (let i = 2; i < end - 2; i++) {
        const value = smooth[i];
        if (value < floor) continue;
        if (value >= smooth[i - 1] && value >= smooth[i - 2] && value > smooth[i + 1] && value > smooth[i + 2]) {
            const half = value / 2;
            let left = i;
            while (left > 0 && smooth[left] > half) left -= 1;
            let right = i;
            while (right < end - 1 && smooth[right] > half) right += 1;
            peaks.push({ r: x[i], fwhm: x[right] - x[left] });
        }
    }
    return peaks.slice(0, maxPeaks);
};

// A dip shallower than this fraction of the first peak lies inside the shell
// (a split shell, like the demo's two Ta–Se distances of one TaSe6
// octahedron), not between shells.
export const SHELL_TROUGH_FRACTION = 0.2;
// The shell starts where g(r) first exceeds this fraction of its peak.
const SHELL_ONSET_FRACTION = 0.01;
// A shell whose tail decays into an empty gap ends where g falls below this
// fraction of its peak (3.7σ out on a Gaussian shell: < 0.1 % of its atoms).
const SHELL_EMPTY_FRACTION = 0.001;

const trapezoid = (x, f, from, to) => {
    let sum = 0;
    for (let i = from; i < to; i++) sum += 0.5 * (f[i] + f[i + 1]) * (x[i + 1] - x[i]);
    return sum;
};

const movingAverage5 = (y, end) => {
    const smooth = new Array(end);
    for (let i = 0; i < end; i++) {
        let sum = 0;
        let n = 0;
        for (let k = -2; k <= 2; k++) {
            const j = i + k;
            if (j >= 0 && j < end) { sum += y[j]; n += 1; }
        }
        smooth[i] = sum / n;
    }
    return smooth;
};

/**
 * The first coordination shell of a partial g(r) (normalised to 1 at large r,
 * as RMCProfile writes it): from the onset, where the smoothed g first exceeds
 * 1 % of the first peak, to the boundary, the first minimum after that peak
 * below SHELL_TROUGH_FRACTION of its height. The shell's distance statistics
 * come from the radial distribution r²·g(r): the g(r) peak itself sits below
 * the mean distance of a broad or skewed shell, because g divides by r².
 * With the partner's number density (atoms/Å³) it also integrates the
 * coordination number CN = 4π ρ ∫ r² g(r) dr, and the CN 0.2 Å either side of
 * the boundary — how much the count depends on where the cut is put.
 * Returns null when no peak is found below rMax.
 */
export const firstShell = (x, y, { rMax = 8, partnerDensity = null } = {}) => {
    if (!Array.isArray(x) || !Array.isArray(y) || x.length !== y.length) return null;
    let end = x.findIndex((value) => value > rMax);
    if (end < 0) end = x.length;
    if (end < 5) return null;
    const peaks = extractPeaks(x.slice(0, end), y.slice(0, end), { rMax, maxPeaks: 8 });
    if (!peaks.length) return null;
    const smooth = movingAverage5(y, end);
    const peakIndex = x.indexOf(peaks[0].r);
    const height = smooth[peakIndex];
    let onset = peakIndex;
    while (onset > 0 && smooth[onset - 1] > SHELL_ONSET_FRACTION * height) onset -= 1;
    // The boundary: a local minimum in a deep trough, or the point where g
    // has fallen back to ~0 (a tail that decays into an empty gap).
    let boundary = -1;
    for (let i = peakIndex + 1; i < end - 1; i++) {
        const deepMinimum = smooth[i] <= smooth[i - 1] && smooth[i] <= smooth[i + 1]
            && smooth[i] < SHELL_TROUGH_FRACTION * height;
        if (deepMinimum || smooth[i] < SHELL_EMPTY_FRACTION * height) {
            boundary = i;
            break;
        }
    }
    const shell = {
        onsetR: x[onset],
        peakR: peaks[0].r,
        peakG: height,
        boundaryR: boundary >= 0 ? x[boundary] : null,
        troughRatio: boundary >= 0 ? Math.max(0, smooth[boundary]) / height : null,
        peaksInShell: peaks.filter((peak) => boundary < 0 || peak.r < x[boundary]).map((peak) => peak.r)
    };
    if (boundary < 0) return shell;
    const radial = y.map((value, i) => Math.max(0, value) * x[i] * x[i]);
    const weight = trapezoid(x, radial, onset, boundary);
    if (weight > 0) {
        const mean = trapezoid(x, radial.map((value, i) => value * x[i]), onset, boundary) / weight;
        const second = trapezoid(x, radial.map((value, i) => value * (x[i] - mean) ** 2), onset, boundary) / weight;
        shell.meanR = mean;
        shell.sdR = Math.sqrt(Math.max(0, second));
    }
    if (Number.isFinite(partnerDensity) && partnerDensity > 0) {
        const cnUpTo = (r) => {
            let to = x.findIndex((value) => value > r);
            if (to < 0) to = x.length - 1;
            return 4 * Math.PI * partnerDensity * trapezoid(x, radial, 0, Math.max(0, to - 1));
        };
        shell.cn = 4 * Math.PI * partnerDensity * trapezoid(x, radial, 0, boundary);
        shell.cnRange = [cnUpTo(x[boundary] - 0.2), cnUpTo(x[boundary] + 0.2)];
    }
    return shell;
};

/** Atoms of `element` per Å³ in the configuration (null without counts or cell). */
export const elementDensity = (structure, element) => {
    const count = structure?.elementCounts?.[element];
    const L = structure?.latticeVectors;
    if (!Number.isFinite(count) || !Array.isArray(L) || L.length !== 3) return null;
    const volume = Math.abs(
        L[0][0] * (L[1][1] * L[2][2] - L[1][2] * L[2][1])
        - L[0][1] * (L[1][0] * L[2][2] - L[1][2] * L[2][0])
        + L[0][2] * (L[1][0] * L[2][1] - L[1][1] * L[2][0])
    );
    return volume > 0 ? count / volume : null;
};

// g(r) at a closest-approach constraint, or at the first shell's boundary,
// above this fraction of the first peak: flagged in the context (the
// find_gr_peaks tool says the same, with the numbers).
export const PAIR_FLAG_FRACTION = 0.02;

/** The largest raw g just above `limit` (its first 2.5 grid steps), or 0. */
export const gJustAbove = (x, y, limit) => {
    const step = x.length > 1 ? Math.abs(x[1] - x[0]) : 0;
    return x.reduce((best, r, i) => (r > limit && r <= limit + 2.5 * step ? Math.max(best, y[i]) : best), 0);
};

/**
 * The pair_correlations context section: one entry per partial-PDF series
 * (labels like "Nb-Nb"), carrying the average-structure nearest-neighbour
 * distance when a basis is available, the extracted g(r) peaks and the first
 * shell. `constraints` is the run settings' labelled min_distances_A map (any
 * pair order); a pair whose distribution is cut at its constraint, or whose
 * first shell runs into the next, carries a flag. Returns null when the run
 * has no parsed partials.
 */
export const pairCorrelationsContext = (structure, plotFiles, { maxPairs = 10, constraints = null } = {}) => {
    const partials = (plotFiles || []).find(
        (file) => file.plotKind === 'pdf_partials' && file.plotData?.series?.length
    );
    if (!partials) return null;
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
        const peaks = extractPeaks(series.x, series.y);
        if (peaks.length) {
            entry.gr_peaks_A = peaks.map((peak) => ({
                r: roundSig(peak.r),
                fwhm: roundSig(peak.fwhm, 2)
            }));
        }
        // The first shell's mean distance (from r²g, not the g peak) and, with
        // the cell, its coordination number: B atoms around each A.
        const shell = firstShell(series.x, series.y, { partnerDensity: elementDensity(structure, match[2]) });
        if (shell?.boundaryR != null && Number.isFinite(shell.meanR)) {
            entry.first_shell = { mean_d_A: roundSig(shell.meanR), to_A: roundSig(shell.boundaryR) };
            if (Number.isFinite(shell.cn)) entry.first_shell[`cn_${match[2]}_per_${match[1]}`] = roundSig(shell.cn);
            // Distances filling the gap to the next shell.
            if (shell.troughRatio > PAIR_FLAG_FRACTION) entry.first_shell.gap_g_fraction = roundSig(shell.troughRatio, 2);
        }
        const limitEntry = constraints && !Array.isArray(constraints)
            ? Object.entries(constraints).find(([key]) => {
                const parts = key.split('-').map((part) => part.trim());
                return parts.length === 2 && pairKey(parts[0], parts[1]) === pairKey(match[1], match[2]);
            })
            : null;
        if (shell && limitEntry && Number.isFinite(limitEntry[1])) {
            const fraction = gJustAbove(series.x, series.y, limitEntry[1]) / shell.peakG;
            // Pairs pile up against the hard constraint: it shapes the peak.
            if (fraction > PAIR_FLAG_FRACTION) entry.cut_at_min_distance = { at_A: limitEntry[1], g_fraction: roundSig(fraction, 2) };
        }
        if (entry.avg_structure_d_A !== undefined || entry.gr_peaks_A) pairs.push(entry);
    }
    return pairs.length ? pairs : null;
};
