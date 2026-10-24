// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// One shell reading per pair, per run (architecture candidate #6). The first
// coordination shell of every element pair in the run's partial-PDF file is
// read here once, and every reader takes it from here: the run context's
// pair_correlations (pairCorrelations.js), find_gr_peaks (runTools.js), the
// Workbench's window questions (clarify.js), the verification
// (copilot/verification/: the evidence gate's finite clusters,
// clusterElementsOf; the answer checks, cnMatches), compute_bond_angles'
// window refusal (copilotTools.js: windowCut) and arm A+ (validation/armA.js).
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
//   split       only for a joined shell (plan item 4.2, approved 2026-10-19):
//               its parts, cut at the minimum of the smoothed g between each
//               two successive resolved peaks in the shell (the first cut is
//               the first peak's), the last part ending at the shell's
//               boundary; each part's CN (the CN integral over exactly its
//               window, so the parts sum to the shell's), its share of the
//               shell's CN and its mean distance and spread (r²g over the
//               part, as the shell's); long − short = the last part's mean −
//               the first's, with its sampling error (each part's spread over
//               √ its bonds in the box, the bonds taken as independent); the
//               counts said in the direction the join decided (`direction`).
//   separation  a shell is separated from the next when the smoothed g at its
//               boundary is at most SEPARATION_G_FRACTION (2 %) of its peak.
//   constraint  a closest-approach constraint limits the pair when the raw g
//               just above it is above CONSTRAINT_LIMIT_FRACTION (5 %) of the
//               first peak: the peak is constraint-limited, and nothing is
//               interpreted on its low-r side (the user's rule, 2026-10-24;
//               2 %, the separation's fraction, before).
//   CN match    two counts agree within half an atom (cnMatches).
//   clusters    an element whose own first shell holds 0.5–4.5 neighbours
//               may form finite clusters (clusterElementsOf: the evidence
//               gate's "no order" and arm A+'s cluster sweep).
//   cut         one rule (windowCut; approved by the user, 2026-10-19): a
//               window that meets a shell cuts it when it holds fewer partner
//               atoms than the whole shell by at least max(0.5, 15 % of the
//               shell's CN) (cnMarginCut). The Workbench asks on it before an
//               analysis runs (clarify.js) and compute_bond_angles refuses on
//               it unless the user chose the window. Only where the run reads
//               no shell to compare the window with — no partials, no closed
//               shell for the pair, no density for a CN, or a window apart
//               from the first shell — does the tool fall back on the
//               coordination counts, as before (countsCut: more than half the
//               central atoms unbonded, or fewer than 60 % with the most
//               common count).
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
// An element whose own first shell holds at most this many neighbours may
// form finite clusters (a dimer 1, a trimer 2, a tetrahedron 3, an
// octahedron 4 per member): the evidence gate's "no order" needs
// analyze_clusters on it, and arm A+ tries its clusters (clusterElementsOf,
// one rule since architecture candidate 8).
export const CLUSTER_CN_MAX = 4.5;
// g at a shell's boundary (smoothed) above this fraction of the first peak:
// the shell runs into the next.
export const SEPARATION_G_FRACTION = 0.02;
// g just above a closest-approach constraint (raw, its first 2.5 grid steps)
// above this fraction of the first peak: the pairs pile up against the
// constraint, which, not the data, shapes the peak's low-r side — the peak is
// constraint-limited, and nothing on its low-r side is interpreted (the
// user's rule, 2026-10-24: "g(r) above 5 % of the peak at the cut"; until
// then the separation's 2 %). The Demo's Ga–Se at its 2.2 Å constraint: 8.7 %.
export const CONSTRAINT_LIMIT_FRACTION = 0.05;
// The cut sensitivity: the shell's top moved this far (Å) either way.
export const CUT_SHIFT_A = 0.2;
// The CN margin: a window changes the answer when its CN differs from the
// whole shell's by at least the larger of these (atoms; fraction of the
// shell's CN). One that meets the shell with fewer atoms by it cuts the shell
// (cnMarginCut): the cut rule of the Workbench's question and of
// compute_bond_angles' refusal (windowCut).
export const CN_CHANGE_ABS = 0.5;
export const CN_CHANGE_REL = 0.15;
// The counts fallback (countsCut), for a window the run reads no shell to
// compare with (windowCut); before 2026-10-19 compute_bond_angles' rule for
// every window. Central atoms with no A–B bond in the window above this
// fraction — the window misses them …
export const UNBONDED_APEX_FRACTION = 0.5;
// … or fewer central atoms than this share with the most common A–B bond
// count — the coordination is spread, as when the window edge cuts through a
// shell (Demo Se–Ta–Se at 2.3–2.7 Å: 47 %; whole shells there: 83–96 %).
export const MODAL_COORDINATION_SHARE = 0.6;
// A split part's mean distance has a sampling error only with at least this
// many bonds in the box (its spread over √ the bonds).
const SPLIT_MIN_BONDS = 2;
// Windows are rounded to this many decimals (0.01 Å, the grid).
const WINDOW_DECIMALS = 2;
// A caller's window covers a shell edge within half a rounding unit of it.
const EDGE_TOLERANCE_A = 0.005 + 1e-9;

// --- small helpers -----------------------------------------------------------------

/** Whether two coordination numbers agree: |cn − count| ≤ CN_MATCH_TOLERANCE (inclusive). */
export const cnMatches = (cn, count) => Number.isFinite(cn) && Number.isFinite(count) && Math.abs(cn - count) <= CN_MATCH_TOLERANCE;

/** Whether a coordination number is an ideal polyhedral count (IDEAL_COORDINATION) by cnMatches. */
export const isIdealCoordination = (cn) => IDEAL_COORDINATION.some((count) => cnMatches(cn, count));

/**
 * The elements of a shell catalogue (shellCatalogue, below) that may form
 * finite clusters: those whose own first shell (the element around itself)
 * holds from CN_MATCH_TOLERANCE to CLUSTER_CN_MAX neighbours, both
 * inclusive — a closed shell, since only a closed shell has a CN. In the
 * catalogue's (alphabetical) order; [] without a catalogue. The evidence
 * gate and arm A+ (validation/armA.js) read this one rule.
 */
export const clusterElementsOf = (catalogue) => (catalogue?.elements || []).filter((element) => {
    const cn = catalogue.shell(element, element)?.cn;
    return Number.isFinite(cn) && cn >= CN_MATCH_TOLERANCE && cn <= CLUSTER_CN_MAX;
});

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

// The mean distance and spread over exactly [lo, hi], from r² max(g, 0): the
// moments of the bond lengths there. Nulls for a window holding nothing.
const momentsOver = ({ x, radial }, [lo, hi]) => {
    const weight = hi > lo ? integrate(x, radial, lo, hi) : 0;
    if (!(weight > 0)) return { meanR: null, sdR: null };
    const meanR = integrate(x, radial.map((value, i) => value * x[i]), lo, hi) / weight;
    const second = integrate(x, radial.map((value, i) => value * (x[i] - meanR) ** 2), lo, hi) / weight;
    return { meanR, sdR: Math.sqrt(Math.max(0, second)) };
};

/**
 * The geometry of the first shell of one partial g(r) (normalised to 1 at
 * large r, as RMCProfile writes it), and the CN of any window over it. With
 * the partner's number density the CN counts partner atoms around each
 * centre; with the centre's too, a split shell can be joined by its count in
 * the reverse direction. Null when no peak lies below PEAK_SEARCH_MAX_A.
 */
const readShell = (data, { partnerDensity = null, centreDensity = null } = {}) => {
    const { x, y } = data;
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
    const moments = (range) => momentsOver(data, range);
    const geometry = { peaks, peaksInShell, peakG: height, joined, cnOver, moments, hasDensity };
    if (boundary < 0) return { ...geometry, window: null };

    const window = [lo, roundWindowEdge(x[boundary])];
    const { meanR, sdR } = moments(window);
    // The minimum of the smoothed g between two resolved peaks (their grid
    // r): the first lowest point on a flat bottom.
    const minimumBetween = (fromR, toR) => {
        const to = x.indexOf(toR);
        let lowest = x.indexOf(fromR);
        for (let i = lowest + 1; i <= to; i += 1) if (smooth[i] < smooth[lowest]) lowest = i;
        return lowest;
    };
    // The first peak alone: only with two resolved peaks in the shell, to the
    // minimum of the smoothed g between the first two.
    let firstPeak = null;
    if (peaksInShell.length >= 2) {
        const lowest = minimumBetween(peaksInShell[0], peaksInShell[1]);
        firstPeak = { window: [lo, roundWindowEdge(x[lowest])], troughRatio: Math.max(0, smooth[lowest]) / height };
    }
    // A joined shell's parts: cut at the minimum between each two successive
    // resolved peaks in it (the first cut is the first peak's), the last part
    // ending at the shell's boundary. Rounded as every window is.
    const splitEdges = joined && peaksInShell.length >= 2
        ? [lo, ...peaksInShell.slice(1).map((r, k) => roundWindowEdge(x[minimumBetween(peaksInShell[k], r)])), window[1]]
        : null;
    return {
        ...geometry,
        window,
        meanR,
        sdR,
        troughRatio: Math.max(0, smooth[boundary]) / height,
        firstPeak,
        splitEdges
    };
};

// A window [lo, hi] measured against a shell: what it holds of the shell.
const holdsOf = ([lo, hi], [shellLo, shellHi]) => {
    if (!(hi > shellLo && lo < shellHi)) return 'apart';
    if (lo > shellLo + EDGE_TOLERANCE_A || hi < shellHi - EDGE_TOLERANCE_A) return 'part';
    return hi > shellHi + EDGE_TOLERANCE_A ? 'more' : 'whole';
};

/**
 * The CN margin's cut: a window cuts the shell when it overlaps it and holds
 * fewer partner atoms than the whole shell by the CN margin (cnChanges). The
 * rule windowCut applies wherever the window meets a shell.
 */
export const cnMarginCut = (cn, shell, [lo, hi]) => Boolean(shell?.window) && Number.isFinite(shell.cn)
    && hi > shell.window[0] && lo < shell.window[1]
    && Number.isFinite(cn) && cn < shell.cn && cnChanges(cn, shell.cn);

/**
 * The counts fallback (windowCut's, for a window with no shell to compare):
 * a window cuts a shell when more than UNBONDED_APEX_FRACTION of the central
 * atoms have no bond in it ('unbonded'), or fewer than
 * MODAL_COORDINATION_SHARE share the most common count ('count-spread').
 * `summary` is the bond-angle engine's coordinationSummary { mean, modal,
 * modalShare, unbondedShare } (triplets.py bond_angle_summary,
 * workers/triplets.js bondAngleSummary). Returns { cut, by }.
 */
export const countsCut = (summary) => {
    if (summary?.unbondedShare > UNBONDED_APEX_FRACTION) return { cut: true, by: 'unbonded' };
    if (Number.isFinite(summary?.modalShare) && summary.modalShare < MODAL_COORDINATION_SHARE) return { cut: true, by: 'count-spread' };
    return { cut: false, by: null };
};

/**
 * Whether a bond window cuts its coordination shell: the one rule (approved
 * by the user, 2026-10-19). The Workbench asks on it before an analysis runs
 * (clarify.js beforeCall, without counts) and compute_bond_angles refuses on
 * it unless the user chose the window (copilotTools.js, with the bond-angle
 * engine's counts). `reading` is the catalogue's window() for the window (null
 * when the run reads no closed shell for the pair, or has no partials);
 * `counts` the engine's coordinationSummary for the window, or null.
 *
 *   - The window meets a shell with a CN (it overlaps the shell, and the
 *     structure gives the partner's density): the CN margin decides
 *     (cnMarginCut; rule 'cn-margin'). `whole`: its CN is the shell's within
 *     the margin, either way (a window that holds more than the shell by the
 *     margin is not cut, and does not count the shell either).
 *   - Otherwise — no partials, no closed shell for the pair, no density for a
 *     CN, or a window apart from the shell (the module reads first shells
 *     only, so it cannot say whether such a window cuts a later one) — the
 *     coordination counts decide, as they did for every window before
 *     2026-10-19 (countsCut; rule 'counts'). `whole` is "not cut", except that
 *     a window apart from the first shell never holds it. Without counts
 *     nothing decides (rule null: not cut, not whole).
 *
 * Returns { cut, by, rule, whole }; `by` is 'cn-margin', 'unbonded',
 * 'count-spread' or null.
 */
export const windowCut = (reading, counts = null) => {
    if (reading && reading.holds !== 'apart' && Number.isFinite(reading.cn) && Number.isFinite(reading.shell?.cn)) {
        const cut = reading.cnMargin.cut;
        return { cut, by: cut ? 'cn-margin' : null, rule: 'cn-margin', whole: !cnChanges(reading.cn, reading.shell.cn) };
    }
    if (!counts) return { cut: false, by: null, rule: null, whole: false };
    const { cut, by } = countsCut(counts);
    return { cut, by, rule: 'counts', whole: !cut && reading?.holds !== 'apart' };
};

/**
 * A joined shell's split (plan item 4.2): `parts` [{ window, cn, reverseCn,
 * share, meanR, sdR }], cut where readShell cut them (splitEdges), each CN
 * integrated over exactly its window (an empty window holds 0), so the parts
 * sum to the shell's CN and their shares to 1; `longMinusShort` = the last
 * part's mean distance − the first's; `longMinusShortSe` its sampling error,
 * √(se_first² + se_last²) with se = the part's spread / √(its bonds in the
 * box: the centres' count × its CN, halved for a same-element pair), the
 * bonds taken as independent — the configuration's own sampling, not the
 * fit's uncertainty (null without the centres' count); and `direction`, the
 * count the join decided by: 'forward' (partner atoms per centre: the
 * shell's own count is ideal, or there is no other) or 'reverse' (centre
 * atoms per partner). Null for a shell that is not joined.
 */
const splitOf = (geometry, { cn, ratio, ideal, centreCount, sameElement }) => {
    const edges = geometry.splitEdges;
    if (!edges || !(cn > 0)) return null;
    const parts = edges.slice(1).map((hi, k) => {
        const window = [edges[k], hi];
        const partCn = hi > edges[k] ? geometry.cnOver(window) : 0;
        const { meanR, sdR } = geometry.moments(window);
        const bonds = centreCount > 0 ? (centreCount * partCn) / (sameElement ? 2 : 1) : 0;
        return {
            window,
            cn: partCn,
            reverseCn: ratio == null ? null : (sameElement ? partCn : partCn * ratio),
            share: partCn / cn,
            meanR,
            sdR,
            se: bonds >= SPLIT_MIN_BONDS && Number.isFinite(sdR) ? sdR / Math.sqrt(bonds) : null
        };
    });
    const first = parts[0];
    const last = parts.at(-1);
    return {
        parts: parts.map(({ se, ...part }) => part),
        longMinusShort: Number.isFinite(first.meanR) && Number.isFinite(last.meanR) ? last.meanR - first.meanR : null,
        longMinusShortSe: first.se != null && last.se != null ? Math.hypot(first.se, last.se) : null,
        direction: ideal != null || ratio == null ? 'forward' : 'reverse'
    };
};

// One Shell (frozen) of a centre–partner pair, from its partial's geometry.
const buildShell = (centre, partner, label, geometry, { partnerDensity, centreDensity, centreCount, constraint, x, y }) => {
    const ratio = Number.isFinite(centreDensity) && centreDensity > 0 && Number.isFinite(partnerDensity) && partnerDensity > 0
        ? centreDensity / partnerDensity
        : null;
    const { window } = geometry;
    const cn = window ? geometry.cnOver(window) : null;
    const ideal = Number.isFinite(cn) ? (IDEAL_COORDINATION.find((count) => cnMatches(cn, count)) ?? null) : null;
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
        ideal,
        joined: Boolean(window) && geometry.joined,
        troughRatio: window ? geometry.troughRatio : null,
        separated: window ? geometry.troughRatio <= SEPARATION_G_FRACTION : null,
        firstPeak: window && geometry.firstPeak ? { ...reading([...geometry.firstPeak.window]), troughRatio: geometry.firstPeak.troughRatio } : null,
        split: window ? splitOf(geometry, { cn, ratio, ideal, centreCount, sameElement: centre === partner }) : null,
        cut: window
            ? {
                // Never below the shell's start.
                earlier: reading([window[0], Math.max(window[0], roundWindowEdge(window[1] - CUT_SHIFT_A))]),
                later: reading([window[0], roundWindowEdge(window[1] + CUT_SHIFT_A)])
            }
            : null,
        constraint: Number.isFinite(limit) ? { limitA: limit, gFraction, limited: gFraction > CONSTRAINT_LIMIT_FRACTION } : null
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
                partnerDensity,
                centreDensity,
                centreCount: structure?.elementCounts?.[centre] ?? null,
                constraint: constraintOf(constraints, centre, partner),
                x: data.x,
                y: data.y
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
         * 'apart') and its CN-margin reading (cnMarginCut), which the one cut
         * rule reads (windowCut). Null when the pair has no closed shell.
         */
        window: (centre, partner, [lo, hi] = []) => {
            const found = shell(centre, partner);
            if (!found?.window || !(hi > lo)) return null;
            const cn = cnOf(sources.get(`${centre}|${partner}`).data, elementDensity(structure, partner), [lo, hi]);
            return {
                cn,
                holds: holdsOf([lo, hi], found.window),
                cnMargin: { cut: cnMarginCut(cn, found, [lo, hi]) },
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
 * partner), window(centre, partner, [lo, hi]) }, or null without a partials
 * file with element-pair series.
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
