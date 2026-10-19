// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// One shell reading per pair, per run (architecture candidate #6). The first
// coordination shell of every element pair in the run's partial-PDF file is
// read here once, and every reader takes it from here: the run context's
// pair_correlations (pairCorrelations.js), find_gr_peaks (runTools.js), the
// Workbench's window questions (clarify.js), the evidence gate (evidence.js),
// the answer checks (answerAudit.js: cnMatches), compute_bond_angles' window
// refusal (copilotTools.js: countsCut) and arm A+ (validation/armA.js).
// Before, five readers each searched their own r range, rounded their own way
// and integrated the coordination number their own way, and one shell could
// read three ways (the Demo's Ta–Se first peak: 3.9 Se per Ta in one, 3.86 in
// another).
//
// The rules, one per quantity (docs/algorithms/ai-assistant.md, Step 9):
//
//   search      peaks below PEAK_SEARCH_MAX_A (8 Å); the minimum that closes
//               the shell is sought to max(8, r₁ + 3) Å. No caller picks a
//               range.
//   shell       from the onset (where the smoothed g first exceeds 1 % of the
//               first peak) to the first deep minimum (below 20 % of that
//               peak) or the point where g falls back to ~0; a split shell cut
//               at a deep trough is joined when a later peak within 0.5 Å
//               brings the summed CN to an ideal polyhedral count either way
//               (the review's L1).
//   rounding    every window the module returns is rounded to 0.01 Å (the
//               grid), and its CN is integrated over the rounded window:
//               every CN returned equals window().cn for the window returned
//               with it.
//   CN          4π ρ ∫ r² max(g, 0) dr, the trapezoid of the linear
//               interpolant over exactly the quoted window (ρ the partner's
//               number density).
//   first peak  only with at least two resolved peaks in the shell: from the
//               shell's start to the minimum of the smoothed g between the
//               first two (approved 2026-10-19). Never r ± FWHM: a half-height
//               walk on overlapping peaks reaches into the second one (the
//               Demo's Ta–Se: 3.86 Se per Ta in r ± FWHM, 3.25 to the
//               inter-peak minimum, 3 in the ideal split).
//   separation  a shell is separated from the next when the smoothed g at its
//               boundary is at most SEPARATION_G_FRACTION (2 %) of its peak;
//               a closest-approach constraint limits the pair when the raw g
//               just above it is above the same fraction.
//   CN match    two counts agree within half an atom (cnMatches).
//   cut         two readings of whether a window cuts a shell, both today's:
//               the Workbench asks on the CN margin (cnMarginCut: the window
//               overlaps the shell and holds fewer atoms by at least
//               max(0.5, 15 % of the shell's)); compute_bond_angles refuses on
//               the coordination counts (countsCut: more than half the
//               central atoms unbonded, or fewer than 60 % with the most
//               common count). Which rule should decide is the user's
//               decision; both readings live here so they can be compared.
//
// Pure functions of plain data (no React), inside the llm module (its import
// boundary: README.md), so the context, the copilot and the host all reach it.

// --- constants ---------------------------------------------------------------------

// Peaks are sought below this r (Å).
export const PEAK_SEARCH_MAX_A = 8;
// The minimum closing the shell is sought to max(PEAK_SEARCH_MAX_A, r₁ + this) Å.
export const BOUNDARY_REACH_A = 3;
// A local maximum of the smoothed g below this fraction of the largest one in
// the search range is noise, not a peak.
const PEAK_FLOOR_FRACTION = 0.15;
// A dip shallower than this fraction of the first peak lies inside the shell
// (a split shell, like the Demo's two Ta–Se distances of one TaSe6
// octahedron), not between shells.
export const SHELL_TROUGH_FRACTION = 0.2;
// The shell starts where g(r) first exceeds this fraction of its peak.
const SHELL_ONSET_FRACTION = 0.01;
// A shell whose tail decays into an empty gap ends where g falls below this
// fraction of its peak (3.7σ out on a Gaussian shell: < 0.1 % of its atoms).
const SHELL_EMPTY_FRACTION = 0.001;
// A split shell whose two bond lengths are well resolved has a DEEP trough
// between them (a Jahn–Teller octahedron, 4 × 1.93 + 2 × 2.25 Å: g falls to
// 12 % of the first peak), so the first deep minimum cuts the polyhedron. A
// later peak within this distance of the first (Å; bond-length splits of
// distorted polyhedra are 0.1–0.4 Å, the next coordination shell lies a
// sizeable fraction of the bond further out) joins the shell when the summed
// coordination number reaches an ideal polyhedral count (the review's L1,
// approved 2026-10-14).
export const SPLIT_SHELL_REACH_A = 0.5;
export const IDEAL_COORDINATION = Object.freeze([4, 6, 8, 12]);
// Two coordination numbers agree within this many atoms (half an atom).
export const CN_MATCH_TOLERANCE = 0.5;
// g at a shell's boundary (smoothed), or just above a closest-approach
// constraint (raw), above this fraction of the first peak: the shell runs
// into the next, or the constraint shapes the peak. One constant for both.
export const SEPARATION_G_FRACTION = 0.02;
// The cut sensitivity: the shell's top moved this far (Å) either way.
export const CUT_SHIFT_A = 0.2;
// The Workbench's CN margin: a window changes the answer when its CN differs
// from the whole shell's by at least the larger of these (atoms; fraction of
// the shell's CN).
export const CN_CHANGE_ABS = 0.5;
export const CN_CHANGE_REL = 0.15;
// compute_bond_angles' counts rule (pinned as it was in copilotTools.js):
// central atoms with no A–B bond in the window above this fraction — the
// window misses them …
export const UNBONDED_APEX_FRACTION = 0.5;
// … or fewer central atoms than this share with the most common A–B bond
// count — the coordination is spread, as when the window edge cuts through a
// shell (Demo Se–Ta–Se at 2.3–2.7 Å: 47 %; whole shells there: 83–96 %).
export const MODAL_COORDINATION_SHARE = 0.6;
// Windows are rounded to this many decimals (0.01 Å, the grid).
const WINDOW_DECIMALS = 2;
// A caller's window covers a shell edge within half a rounding unit of it.
const EDGE_TOLERANCE_A = 0.005 + 1e-9;

// --- small helpers -----------------------------------------------------------------

/** Whether two coordination numbers agree: |cn − count| ≤ CN_MATCH_TOLERANCE (inclusive). */
export const cnMatches = (cn, count) => Number.isFinite(cn) && Number.isFinite(count) && Math.abs(cn - count) <= CN_MATCH_TOLERANCE;

/** Whether a coordination number is an ideal polyhedral count (IDEAL_COORDINATION) by cnMatches. */
export const isIdealCoordination = (cn) => IDEAL_COORDINATION.some((count) => cnMatches(cn, count));

/** The Workbench's CN margin: |cn − reference| ≥ max(CN_CHANGE_ABS, CN_CHANGE_REL · reference). */
export const cnChanges = (cn, reference) => Number.isFinite(cn) && Number.isFinite(reference)
    && Math.abs(cn - reference) >= Math.max(CN_CHANGE_ABS, CN_CHANGE_REL * reference);

const roundWindowEdge = (value) => Number(value.toFixed(WINDOW_DECIMALS)) + 0;

const ELEMENT_PAIR = /^([A-Z][a-z]?)\s*-\s*([A-Z][a-z]?)$/;
const pairOf = (label) => {
    const match = ELEMENT_PAIR.exec(String(label || '').trim());
    return match ? [match[1], match[2]] : null;
};
const unorderedKey = (a, b) => [a, b].sort().join('-');

const deepFreeze = (value) => {
    if (value && typeof value === 'object' && !Object.isFrozen(value)) {
        Object.values(value).forEach(deepFreeze);
        Object.freeze(value);
    }
    return value;
};

// The first index whose r exceeds rMax (x.length when none does).
const endIndex = (x, rMax) => {
    const end = x.findIndex((value) => value > rMax);
    return end < 0 ? x.length : end;
};

// 5-point moving average of y[0 … end) (fewer points at the ends).
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

// ∫ f dr over exactly [lo, hi]: the trapezoid of f's linear interpolant
// between grid points, the edges interpolated.
const integrate = (x, f, lo, hi) => {
    let sum = 0;
    for (let i = 0; i < x.length - 1; i += 1) {
        const a = Math.max(lo, x[i]);
        const b = Math.min(hi, x[i + 1]);
        if (b <= a) continue;
        const slope = (f[i + 1] - f[i]) / (x[i + 1] - x[i]);
        const fa = f[i] + slope * (a - x[i]);
        const fb = f[i] + slope * (b - x[i]);
        sum += 0.5 * (fa + fb) * (b - a);
    }
    return sum;
};

/**
 * Peaks of a g(r)-like curve below rMax: local maxima of a lightly smoothed
 * copy, above a floor relative to the window maximum, with an approximate full
 * width at half maximum from a half-height walk (overlapping peaks make the
 * width an overestimate — a display value, never a window). Returns
 * [{ r, fwhm }] sorted by r, at most maxPeaks entries.
 */
export const extractPeaks = (x, y, { rMax = 6, maxPeaks = 2 } = {}) => {
    if (!Array.isArray(x) || !Array.isArray(y) || x.length !== y.length) return [];
    const end = endIndex(x, rMax);
    if (end < 5) return [];
    // 5-point moving average so grid noise does not mint peaks.
    const smooth = movingAverage5(y, end);
    const windowMax = Math.max(...smooth);
    if (!(windowMax > 0)) return [];
    const floor = PEAK_FLOOR_FRACTION * windowMax;
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

// The largest raw g just above `limit` (its first 2.5 grid steps), or 0.
const gJustAbove = (x, y, limit) => {
    const step = x.length > 1 ? Math.abs(x[1] - x[0]) : 0;
    return x.reduce((best, r, i) => (r > limit && r <= limit + 2.5 * step ? Math.max(best, y[i]) : best), 0);
};

// Element-pair labels in RMCProfile's upper-triangle row-wise order (Ga-Ga,
// Ga-Nb, Ga-Se, Nb-Nb, …) — the order MINIMUM_DISTANCES values are given in,
// matching the partial-PDF column order.
const pairLabels = (atoms) => {
    const labels = [];
    for (let i = 0; i < atoms.length; i++) {
        for (let j = i; j < atoms.length; j++) labels.push(`${atoms[i]}-${atoms[j]}`);
    }
    return labels;
};

/**
 * The run settings' closest-approach constraints (parsed from <stem>.dat,
 * browserData.parseRunSettings), labelled per element pair when the counts
 * line up ({ 'Ga-Se': 2.2, … }), else the plain list as given; null without
 * any. The run context's run_settings.min_distances_A and the shells'
 * constraint read this one labelling.
 */
export const labelledMinimumDistances = (settings) => {
    const distances = settings?.minimumDistancesA;
    if (!Array.isArray(distances) || !distances.length) return null;
    const labels = pairLabels(settings.atoms || []);
    return labels.length === distances.length
        ? Object.fromEntries(distances.map((d, i) => [labels[i], d]))
        : distances;
};

// The labelled constraint of one element pair, matched in either order.
const constraintOf = (constraints, a, b) => {
    if (!constraints || Array.isArray(constraints)) return null;
    const entry = Object.entries(constraints).find(([key]) => {
        const parts = key.split('-').map((part) => part.trim());
        return parts.length === 2 && unorderedKey(parts[0], parts[1]) === unorderedKey(a, b);
    });
    return entry && Number.isFinite(entry[1]) ? entry[1] : null;
};

// --- the shell of one partial ------------------------------------------------------

// A partial as the shells read it: r, g and the radial distribution r²·max(g, 0).
const prepare = (entry) => ({
    x: entry.x,
    y: entry.y,
    radial: entry.y.map((value, i) => Math.max(0, value) * entry.x[i] * entry.x[i])
});

// The CN of [lo, hi]: 4π ρ ∫ r² max(g, 0) dr over exactly that window; null
// without a density or for an empty window.
const cnOf = ({ x, radial }, density, [lo, hi]) => (
    Number.isFinite(density) && density > 0 && hi > lo ? 4 * Math.PI * density * integrate(x, radial, lo, hi) : null
);

/**
 * The geometry of the first shell of one partial g(r) (normalised to 1 at
 * large r, as RMCProfile writes it), and the CN of any window over it. With
 * the partner's number density the CN counts partner atoms around each
 * centre; with the centre's too, a split shell can be joined by its count in
 * the reverse direction. Null when no peak lies below PEAK_SEARCH_MAX_A.
 */
const readShell = (data, { partnerDensity = null, centreDensity = null } = {}) => {
    const { x, y, radial } = data;
    const peaks = extractPeaks(x, y, { rMax: PEAK_SEARCH_MAX_A, maxPeaks: 8 });
    if (!peaks.length) return null;
    const end = endIndex(x, Math.max(PEAK_SEARCH_MAX_A, peaks[0].r + BOUNDARY_REACH_A));
    const smooth = movingAverage5(y, end);
    const step = x.length > 1 ? Math.abs(x[1] - x[0]) : 0;
    const hasDensity = Number.isFinite(partnerDensity) && partnerDensity > 0;
    const cnOver = (range) => cnOf(data, partnerDensity, range);

    const peakIndex = x.indexOf(peaks[0].r);
    const height = smooth[peakIndex];
    let onset = peakIndex;
    while (onset > 0 && smooth[onset - 1] > SHELL_ONSET_FRACTION * height) onset -= 1;
    const lo = roundWindowEdge(Math.max(0, x[onset] - step));
    // The boundary: a local minimum in a deep trough, or the point where g
    // has fallen back to ~0 (a tail that decays into an empty gap). Measured
    // against the height of the peak that opens the shell part (`from`).
    const boundaryAfter = (from, top) => {
        for (let i = from + 1; i < end - 1; i++) {
            const deepMinimum = smooth[i] <= smooth[i - 1] && smooth[i] <= smooth[i + 1]
                && smooth[i] < SHELL_TROUGH_FRACTION * top;
            if (deepMinimum || smooth[i] < SHELL_EMPTY_FRACTION * top) return i;
        }
        return -1;
    };
    let boundary = boundaryAfter(peakIndex, height);
    let joined = false;
    if (boundary >= 0 && hasDensity) {
        // Later peaks of the same split shell: each closes at its own deep
        // minimum; keep the furthest whose summed CN (over the rounded
        // window, as it is reported) is an ideal count either way.
        const idealEitherWay = (cn) => isIdealCoordination(cn)
            || (Number.isFinite(centreDensity) && centreDensity > 0 && isIdealCoordination((cn * centreDensity) / partnerDensity));
        let reach = boundary;
        let joinedAt = -1;
        for (const peak of peaks) {
            if (peak.r - peaks[0].r > SPLIT_SHELL_REACH_A) break;
            const index = x.indexOf(peak.r);
            if (index <= reach) continue;
            const closes = boundaryAfter(index, smooth[index]);
            if (closes < 0) break;
            reach = closes;
            if (idealEitherWay(cnOver([lo, roundWindowEdge(x[closes])]))) joinedAt = closes;
        }
        if (joinedAt > boundary) {
            boundary = joinedAt;
            joined = true;
        }
    }
    const peaksInShell = peaks.filter((peak) => boundary < 0 || peak.r < x[boundary]).map((peak) => peak.r);
    const geometry = { peaks, peaksInShell, peakG: height, joined, cnOver, hasDensity };
    if (boundary < 0) return { ...geometry, window: null };

    const window = [lo, roundWindowEdge(x[boundary])];
    const weight = integrate(x, radial, window[0], window[1]);
    let meanR = null;
    let sdR = null;
    if (weight > 0) {
        meanR = integrate(x, radial.map((value, i) => value * x[i]), window[0], window[1]) / weight;
        const second = integrate(x, radial.map((value, i) => value * (x[i] - meanR) ** 2), window[0], window[1]) / weight;
        sdR = Math.sqrt(Math.max(0, second));
    }
    // The first peak alone: only with two resolved peaks in the shell, to the
    // minimum of the smoothed g between the first two (the first lowest point
    // on a flat bottom).
    let firstPeak = null;
    if (peaksInShell.length >= 2) {
        const from = x.indexOf(peaksInShell[0]);
        const to = x.indexOf(peaksInShell[1]);
        let lowest = from;
        for (let i = from + 1; i <= to; i += 1) if (smooth[i] < smooth[lowest]) lowest = i;
        firstPeak = { window: [lo, roundWindowEdge(x[lowest])], troughRatio: Math.max(0, smooth[lowest]) / height };
    }
    return {
        ...geometry,
        window,
        meanR,
        sdR,
        troughRatio: Math.max(0, smooth[boundary]) / height,
        firstPeak
    };
};

// A window [lo, hi] measured against a shell: what it holds of the shell.
const holdsOf = ([lo, hi], [shellLo, shellHi]) => {
    if (!(hi > shellLo && lo < shellHi)) return 'apart';
    if (lo > shellLo + EDGE_TOLERANCE_A || hi < shellHi - EDGE_TOLERANCE_A) return 'part';
    return hi > shellHi + EDGE_TOLERANCE_A ? 'more' : 'whole';
};

/**
 * The Workbench's rule (today's, unchanged): a window cuts the shell when it
 * overlaps it and holds fewer partner atoms than the whole shell by the CN
 * margin (cnChanges).
 */
export const cnMarginCut = (cn, shell, [lo, hi]) => Boolean(shell?.window) && Number.isFinite(shell.cn)
    && hi > shell.window[0] && lo < shell.window[1]
    && Number.isFinite(cn) && cn < shell.cn && cnChanges(cn, shell.cn);

/**
 * compute_bond_angles' rule (today's, pinned): a window cuts the shell when
 * more than UNBONDED_APEX_FRACTION of the central atoms have no bond in it
 * ('unbonded'), or fewer than MODAL_COORDINATION_SHARE share the most common
 * count ('count-spread'). `summary` is the bond-angle engine's
 * coordinationSummary { mean, modal, modalShare, unbondedShare }
 * (triplets.py bond_angle_summary, workers/triplets.js bondAngleSummary).
 * Returns { cut, by }.
 */
export const countsCut = (summary) => {
    if (summary?.unbondedShare > UNBONDED_APEX_FRACTION) return { cut: true, by: 'unbonded' };
    if (Number.isFinite(summary?.modalShare) && summary.modalShare < MODAL_COORDINATION_SHARE) return { cut: true, by: 'count-spread' };
    return { cut: false, by: null };
};

// One Shell (frozen) of a centre–partner pair, from its partial's geometry.
const buildShell = (centre, partner, label, geometry, { partnerDensity, centreDensity, constraint, x, y }) => {
    const ratio = Number.isFinite(centreDensity) && centreDensity > 0 && Number.isFinite(partnerDensity) && partnerDensity > 0
        ? centreDensity / partnerDensity
        : null;
    const { window } = geometry;
    const cn = window ? geometry.cnOver(window) : null;
    // A window and its CN; an empty one (a shell narrower than the cut) holds no atoms.
    const reading = (range) => {
        const value = range[1] > range[0] ? geometry.cnOver(range) : (geometry.hasDensity ? 0 : null);
        return { window: range, cn: value, changesCn: cnChanges(value, cn) };
    };
    const limit = constraint;
    const gFraction = Number.isFinite(limit) && geometry.peakG > 0 ? gJustAbove(x, y, limit) / geometry.peakG : null;
    return deepFreeze({
        centre,
        partner,
        label,
        peaks: geometry.peaks.map((peak) => ({ r: peak.r, fwhm: peak.fwhm })),
        window: window ? [...window] : null,
        peaksInShell: [...geometry.peaksInShell],
        meanR: window ? geometry.meanR : null,
        sdR: window ? geometry.sdR : null,
        peakG: geometry.peakG,
        cn,
        reverseCn: Number.isFinite(cn) && ratio != null ? (centre === partner ? cn : cn * ratio) : null,
        ideal: Number.isFinite(cn) ? (IDEAL_COORDINATION.find((count) => cnMatches(cn, count)) ?? null) : null,
        joined: Boolean(window) && geometry.joined,
        troughRatio: window ? geometry.troughRatio : null,
        separated: window ? geometry.troughRatio <= SEPARATION_G_FRACTION : null,
        firstPeak: window && geometry.firstPeak ? { ...reading([...geometry.firstPeak.window]), troughRatio: geometry.firstPeak.troughRatio } : null,
        cut: window
            ? {
                // Never below the shell's start.
                earlier: reading([window[0], Math.max(window[0], roundWindowEdge(window[1] - CUT_SHIFT_A))]),
                later: reading([window[0], roundWindowEdge(window[1] + CUT_SHIFT_A)])
            }
            : null,
        constraint: Number.isFinite(limit) ? { limitA: limit, gFraction, limited: gFraction > SEPARATION_G_FRACTION } : null
    });
};

// --- the catalogue -----------------------------------------------------------------

const partialsOf = (plotFiles) => (plotFiles || []).find(
    (file) => file?.plotKind === 'pdf_partials' && file.plotData?.series?.length
) || null;

const buildCatalogue = (partials, structure, runSettings) => {
    const constraints = labelledMinimumDistances(runSettings);
    // Each ordered pair's partial: the file's own label first; the reverse of
    // a pair the file does not list reads the same g.
    const sources = new Map();
    const labelled = [];
    for (const entry of partials.plotData.series) {
        const pair = pairOf(entry?.label);
        if (!pair || !Array.isArray(entry.x) || !Array.isArray(entry.y) || entry.x.length !== entry.y.length) continue;
        const key = `${pair[0]}|${pair[1]}`;
        if (sources.has(key)) continue;
        sources.set(key, { data: prepare(entry), label: String(entry.label).trim() });
        labelled.push(key);
    }
    for (const key of labelled) {
        const [a, b] = key.split('|');
        if (!sources.has(`${b}|${a}`)) sources.set(`${b}|${a}`, sources.get(key));
    }
    if (!sources.size) return null;
    const shells = new Map();
    for (const [key, { data, label }] of sources) {
        const [centre, partner] = key.split('|');
        const partnerDensity = elementDensity(structure, partner);
        const centreDensity = elementDensity(structure, centre);
        const geometry = readShell(data, { partnerDensity, centreDensity });
        shells.set(key, geometry
            ? buildShell(centre, partner, label, geometry, {
                partnerDensity, centreDensity, constraint: constraintOf(constraints, centre, partner), x: data.x, y: data.y
            })
            : null);
    }
    const shell = (centre, partner) => shells.get(`${centre}|${partner}`) ?? null;
    const elements = [...new Set([...sources.keys()].flatMap((key) => key.split('|')))].sort();
    return Object.freeze({
        elements: Object.freeze(elements),
        // One shell per labelled series, in the file's order and orientation.
        pairs: Object.freeze(labelled.map((key) => shells.get(key)).filter(Boolean)),
        shell,
        /**
         * A window [lo, hi] (Å, as given: not rounded) against the
         * centre–partner shell: its CN (partner atoms per centre, integrated
         * over exactly that window), what it holds of the shell ('whole',
         * 'part', 'more' — the whole shell and beyond its boundary — or
         * 'apart'), the Workbench's CN-margin reading and, given the window's
         * coordination counts (the bond-angle engine's coordinationSummary),
         * compute_bond_angles' counts reading. Null when the pair has no
         * closed shell.
         */
        window: (centre, partner, [lo, hi] = [], { counts = null } = {}) => {
            const found = shell(centre, partner);
            if (!found?.window || !(hi > lo)) return null;
            const cn = cnOf(sources.get(`${centre}|${partner}`).data, elementDensity(structure, partner), [lo, hi]);
            return {
                cn,
                holds: holdsOf([lo, hi], found.window),
                cnMargin: { cut: cnMarginCut(cn, found, [lo, hi]) },
                counts: counts ? countsCut(counts) : null,
                shell: found
            };
        }
    });
};

// Memoised per partials data, structure and run settings (object identity):
// the context, the tools, the Workbench's questions and the evidence gate of
// one turn read one catalogue.
const NONE = Object.freeze({});
const CATALOGUES = new WeakMap();

/**
 * The shell catalogue of a run: `run` holds the Copilot's run props
 * ({ plotFiles, structure, runSettings }; only the partial-PDF file is
 * required — without the structure's counts and cell there is no CN, without
 * the run settings no constraint). Returns { elements, pairs, shell(centre,
 * partner), window(centre, partner, [lo, hi], { counts }) }, or null without a
 * partials file with element-pair series.
 */
export const shellCatalogue = (run) => {
    const partials = partialsOf(run?.plotFiles);
    if (!partials) return null;
    const structure = run.structure ?? null;
    const runSettings = run.runSettings ?? null;
    let byStructure = CATALOGUES.get(partials.plotData.series);
    if (!byStructure) {
        byStructure = new WeakMap();
        CATALOGUES.set(partials.plotData.series, byStructure);
    }
    let bySettings = byStructure.get(structure || NONE);
    if (!bySettings) {
        bySettings = new WeakMap();
        byStructure.set(structure || NONE, bySettings);
    }
    const key = runSettings || NONE;
    if (!bySettings.has(key)) bySettings.set(key, buildCatalogue(partials, structure, runSettings));
    return bySettings.get(key);
};
