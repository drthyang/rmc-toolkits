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

/**
 * The pair_correlations context section: one entry per partial-PDF series
 * (labels like "Nb-Nb"), carrying the average-structure nearest-neighbour
 * distance when a basis is available and the extracted g(r) peaks. Returns
 * null when the run has no parsed partials.
 */
export const pairCorrelationsContext = (structure, plotFiles, { maxPairs = 10 } = {}) => {
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
        if (entry.avg_structure_d_A !== undefined || entry.gr_peaks_A) pairs.push(entry);
    }
    return pairs.length ? pairs : null;
};
