// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Displacement modulations — static-mode port of rmc_toolkits/modulations.py
// (the source of truth; keep the two in sync, the parity fixture pins them).
//
// The average structure folds every cell into one, so it cannot show order
// that enlarges the cell. The configuration can: an ordered pattern with
// wavevector q puts power into the Fourier component, at q, of every site's
// displacement field over the supercell's cell grid. Per site the power
// |U_s(q)|² / (n_s tr C_s) is 1 on average for independent displacements;
// an element's sites are averaged; the zone-boundary points and the
// strongest modes are reported against the empirical distribution over all
// q ≠ 0 (ratio to its median), and ranked among the compared modes: one per
// ±q pair, no long-wavelength wave (every |q_i| ≤ LONG_WAVE_Q) — percentile
// and rank. numpy's fftn is replaced by a separable DFT — the same numbers to
// floating-point noise. The engine also decides what these show (the verdict,
// modulations.py's): each tested point's calibrated σ, kind, ordered, because
// and bar, each element's cellDoubling / generalOrder / longWave, the family
// and the rules — which the Copilot's analyze_modulations only renders.

import { normalQuantile, regularizedGamma } from './orientation.js';

export const TOP_MODES = 5;
export const PERCENTILE_RTOL = 1e-9;
// A mode with every component within this of the zone centre is a
// long-wavelength wave, left out of the comparison (modulations.py).
export const LONG_WAVE_Q = 0.25;
// The long-wave test's tolerance: a grid component k/N is 1/4 exactly or at
// least 1/(4N) from it, so on any grid it changes nothing.
export const LONG_WAVE_TOL = 1e-12;

// --- the verdict's rules (modulations.py; keep the two in sync) ------------------
//
// A zone-boundary point is ordered displacement when its power is at least
// MODULATION_RATIO times the median over all wavevectors AND it is significant
// at MODULATION_SIGMA after a family-wise (Šidák) correction over every
// zone-boundary point and element tested: for independent displacements a
// site's normalised power at a half-integer q is Σ_c w_c χ²₁ (Σ w_c = 1), so an
// element's mean over m sites has at least m degrees of freedom, P ≥ χ²_m/m. At
// a general q (neither a zone-boundary point nor a long wave) the sum is
// complex: P ≥ χ²_2m/2m, corrected over the compared modes. These were the
// Copilot tool's rules until 2026-10-19, ported unchanged.

/** Ratio to the median power an ordered point needs (the effect size). */
export const MODULATION_RATIO = 3;
/** Family-wise one-sided significance (σ) an ordered point needs. */
export const MODULATION_SIGMA = 3;
/** An element's strongest mode this many times the median, and a long wave, is a box-scale wave. */
export const LONG_WAVE_RATIO = 5;
/** σ = +∞ (an upper tail that underflows) is reported as this; no finite double reaches it (≈ 38.5 at p = 5e-324). σ = −∞ is null. */
export const SIGMA_MAX = 40;
const P_CAP = 1 - 1e-16;

export const MODULATION_RULES = Object.freeze({
    ratio: MODULATION_RATIO,
    sigma: MODULATION_SIGMA,
    longWaveRatio: LONG_WAVE_RATIO,
    longWaveQ: LONG_WAVE_Q,
    sigmaMax: SIGMA_MAX
});
/** Every `because` a point and each element finding can give (closed sets, pinned by the parity fixture). */
export const MODULATION_BECAUSE = Object.freeze({
    point: ['ordered', 'below_ratio', 'below_sigma', 'undefined_statistic', 'zone_boundary', 'long_wave'],
    cellDoubling: ['ordered_point', 'no_ordered_point', 'undefined_statistic', 'untestable', 'no_displacements'],
    generalOrder: ['ordered_mode', 'no_ordered_mode', 'undefined_statistic', 'no_general_mode', 'no_displacements'],
    longWave: ['long_wave_dominates', 'below_ratio', 'strongest_not_long_wave', 'undefined_statistic', 'no_displacements']
});
export const MODE_KINDS = Object.freeze(['zone_boundary', 'long_wave', 'general']);

/** A long-wavelength mode: every component within LONG_WAVE_Q of the zone centre. The one predicate. */
export const isLongWave = (q) => q.every((value) => Math.abs(value) <= LONG_WAVE_Q + LONG_WAVE_TOL);

// One-sided σ of 1 − (1 − p)^family, without cancellation for small p: +∞
// where it underflows to 0, −∞ where it rounds to 1 (modulations.py _sidak_sigma).
const sidakSigma = (p, family) => {
    const familyP = -Math.expm1(family * Math.log1p(-Math.min(p, P_CAP)));
    if (familyP <= 0) return Infinity;
    if (familyP >= 1) return -Infinity;
    return -normalQuantile(familyP);
};
// σ as reported: −∞ as null, +∞ as SIGMA_MAX.
const reported = (sigma) => (sigma === -Infinity ? null : Math.min(sigma, SIGMA_MAX));

/**
 * The family-wise significance (one-sided σ) of an element's zone-boundary
 * power (normalised: 1 = independent displacements) over `sites` sites,
 * among `family` points tested. Null where undefined (power ≤ 0) or −∞;
 * +∞ is SIGMA_MAX. Mirrors zone_boundary_sigma.
 */
export const zoneBoundarySigma = (power, sites, family) => {
    if (!(power > 0) || !(sites > 0)) return null;
    return reported(sidakSigma(regularizedGamma(sites / 2, (sites * power) / 2).upper, family));
};

/** The same at a general wavevector (a complex sum, P ≥ χ²_2m/2m), `family` = the compared modes. Mirrors general_mode_sigma. */
export const generalModeSigma = (power, sites, family) => {
    if (!(power > 0) || !(sites > 0) || !(family >= 1)) return null;
    return reported(sidakSigma(regularizedGamma(sites, sites * power).upper, family));
};

/** The half-integer points (0/1 flags per axis) an even supercell can hold. */
export const zoneBoundaryPoints = (supercell) => {
    const points = [];
    for (const a of [0, 1]) {
        for (const b of [0, 1]) {
            for (const c of [0, 1]) {
                const flags = [a, b, c];
                if (!flags.some(Boolean)) continue;
                if (flags.every((flag, axis) => !flag || supercell[axis] % 2 === 0)) points.push(flags);
            }
        }
    }
    return points;
};

const folded = (index, size) => {
    const value = index / size;
    return value > 0.5 ? value - 1 : value;
};

// In-place DFT of `re`/`im` (length n·stride… strided lines) along one axis
// of a C-ordered [N0][N1][N2] grid: X_k = Σ_n x_n e^{−2πi kn/N}.
const dftAxis = (re, im, dims, axis) => {
    const n = dims[axis];
    const strides = [dims[1] * dims[2], dims[2], 1];
    const stride = strides[axis];
    const cos = Float64Array.from({ length: n * n }, (_, i) => Math.cos(-2 * Math.PI * ((Math.floor(i / n) * (i % n)) % n) / n));
    const sin = Float64Array.from({ length: n * n }, (_, i) => Math.sin(-2 * Math.PI * ((Math.floor(i / n) * (i % n)) % n) / n));
    const lineRe = new Float64Array(n);
    const lineIm = new Float64Array(n);
    const total = dims[0] * dims[1] * dims[2];
    for (let start = 0; start < total; start += 1) {
        // A line starts where this axis's index is 0.
        if (Math.floor(start / stride) % n !== 0) continue;
        for (let k = 0; k < n; k += 1) {
            lineRe[k] = re[start + k * stride];
            lineIm[k] = im[start + k * stride];
        }
        for (let k = 0; k < n; k += 1) {
            let sumRe = 0;
            let sumIm = 0;
            for (let j = 0; j < n; j += 1) {
                const c = cos[k * n + j];
                const s = sin[k * n + j];
                sumRe += lineRe[j] * c - lineIm[j] * s;
                sumIm += lineRe[j] * s + lineIm[j] * c;
            }
            re[start + k * stride] = sumRe;
            im[start + k * stride] = sumIm;
        }
    }
};

const capitalize = (symbol) => (symbol ? symbol[0].toUpperCase() + symbol.slice(1).toLowerCase() : symbol);

/**
 * Per-element displacement power at every wavevector of the supercell grid,
 * from a parsed configuration (siteDisplacementsFromRmc6f). Mirrors
 * rmc_toolkits.modulations.displacement_modulations; camelCase payload.
 */
export const displacementModulations = (parsed, element = null, { top = TOP_MODES } = {}) => {
    const sites = parsed.sites;
    if (parsed.reconstructed || sites.some((site) => !Array.isArray(site.cells))) {
        throw new Error('modulation analysis needs the reference-site and cell columns of the .rmc6f file');
    }
    const supercell = parsed.supercell.map(Number);
    const nQ = supercell[0] * supercell[1] * supercell[2];
    if (nQ < 8) throw new Error(`a [${supercell.join(', ')}] supercell has too few wavevectors to compare`);
    const labels = sites.map((site) => site.element);
    let elements;
    if (element != null) {
        const symbol = capitalize(String(element).trim());
        if (!labels.includes(symbol)) {
            const species = [...new Set(sites.flatMap((site) => Object.keys(site.elementCounts || { [site.element]: 1 })))].sort();
            throw new Error(`no site holds ${symbol}; elements: ${species.join(', ')}`);
        }
        elements = [symbol];
    } else {
        elements = [...new Set(labels)].sort();
    }
    const boundary = zoneBoundaryPoints(supercell);
    const untestable = [];
    for (const a of [0, 1]) for (const b of [0, 1]) for (const c of [0, 1]) {
        const flags = [a, b, c];
        if (flags.some(Boolean) && !boundary.some((p) => p.every((v, i) => v === flags[i]))) untestable.push(flags);
    }
    const result = {
        supercell,
        wavevectors: nQ - 1,
        untestable,
        elements: [],
        parseWarning: parsed.parseWarning ?? null
    };
    const linear = ([i, j, k]) => {
        const w = [i, j, k].map((v, axis) => ((v % supercell[axis]) + supercell[axis]) % supercell[axis]);
        return (w[0] * supercell[1] + w[1]) * supercell[2] + w[2];
    };
    const unravel = (q) => [Math.floor(q / (supercell[1] * supercell[2])), Math.floor(q / supercell[2]) % supercell[1], q % supercell[2]];
    // The compared modes: one per ±q pair (the smaller grid index), q ≠ 0, no long wave.
    const compared = [];
    for (let q = 1; q < nQ; q += 1) {
        const index = unravel(q);
        const partner = linear(index.map((v) => -v));
        if (q <= partner && !isLongWave(index.map((v, axis) => folded(v, supercell[axis])))) compared.push(q);
    }
    for (const name of elements) {
        const power = new Float64Array(nQ);
        let missing = 0;
        let used = 0;
        sites.forEach((site) => {
            if (site.element !== name) return;
            const u = site.displacements;
            const n = u.length;
            const mean = [0, 1, 2].map((k) => u.reduce((sum, row) => sum + row[k], 0) / n);
            const variance = [0, 1, 2].reduce((sum, k) => sum + u.reduce((acc, row) => acc + (row[k] - mean[k]) ** 2, 0) / n, 0);
            if (!(variance > 0)) return;
            missing += nQ - n;
            for (let component = 0; component < 3; component += 1) {
                const re = new Float64Array(nQ);
                const im = new Float64Array(nQ);
                site.cells.forEach((cell, row) => { re[linear(cell)] = u[row][component]; });
                for (let axis = 0; axis < 3; axis += 1) dftAxis(re, im, supercell, axis);
                for (let q = 0; q < nQ; q += 1) power[q] += (re[q] * re[q] + im[q] * im[q]) / (n * variance);
            }
            used += 1;
        });
        if (!used) {
            result.elements.push({ element: name, sites: 0, verdict: untestedVerdict() });
            continue;
        }
        for (let q = 0; q < nQ; q += 1) power[q] /= used;
        const others = Array.from(power.subarray(1));
        const sorted = [...others].sort((a, b) => a - b);
        const mid = sorted.length >> 1;
        const median = sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
        const reference = compared.map((q) => power[q]);
        const describe = (q) => {
            const value = power[q];
            const index = unravel(q);
            let lower = 0;
            let higher = 0;
            reference.forEach((other) => {
                if (other < value * (1 - PERCENTILE_RTOL)) lower += 1;
                if (other > value * (1 + PERCENTILE_RTOL)) higher += 1;
            });
            return {
                q: index.map((v, axis) => folded(v, supercell[axis])),
                ratio: median > 0 ? value / median : null,
                percentile: reference.length ? (lower / reference.length) * 100 : null,
                rank: reference.length ? higher + 1 : null
            };
        };
        const zoneBoundary = boundary.map((flags) => {
            const entry = describe(linear(flags.map((flag, axis) => (flag ? supercell[axis] / 2 : 0))));
            entry.q = flags.map((flag) => 0.5 * flag);
            return entry;
        });
        const order = Array.from({ length: nQ }, (_, q) => q).sort((a, b) => (power[b] - power[a]) || (a - b));
        const strongest = [];
        const seen = new Set();
        for (const q of order) {
            if (q === 0) continue;
            const index = unravel(q);
            const partner = index.map((v, axis) => ((-v % supercell[axis]) + supercell[axis]) % supercell[axis]);
            // +q and −q carry the same power (to rounding): report the pair
            // once, as its member with the smaller grid index, in both engines.
            const key = linear(index) < linear(partner) ? linear(index) : linear(partner);
            if (seen.has(key)) continue;
            seen.add(key);
            const mode = describe(key);
            // A zone-boundary point (every index 0 or N/2) is tested in
            // zoneBoundary; a long wave is not tested; the rest are general.
            const keyIndex = unravel(key);
            const boundaryMode = keyIndex.every((v, axis) => v === 0 || 2 * v === supercell[axis]);
            mode.kind = boundaryMode ? 'zone_boundary' : (isLongWave(mode.q) ? 'long_wave' : 'general');
            strongest.push(mode);
            if (strongest.length === top) break;
        }
        result.elements.push({
            element: name, sites: used, missingCopies: missing, medianPower: median, comparedModes: reference.length, zoneBoundary, strongest
        });
    }
    decide(result);
    return result;
};

// --- the verdict (modulations.py _decide; keep the two in sync) ------------------

// One point against its bar. The power is ratio × median, as the tool formed
// it: σ far below the background is ill-conditioned, and this keeps it bit for
// bit the value the tool reported.
const testPoint = (point, test, entry, family) => {
    const { ratio } = point;
    const power = ratio != null ? ratio * entry.medianPower : null;
    const { sites } = entry;
    const isDefined = power != null && power > 0 && sites > 0 && (test === 'zone' || family >= 1);
    let sigma = null;
    let ordered;
    let because;
    if (!isDefined) {
        ordered = null;
        because = 'undefined_statistic';
    } else {
        sigma = test === 'zone' ? zoneBoundarySigma(power, sites, family) : generalModeSigma(power, sites, family);
        if (ratio < MODULATION_RATIO) [ordered, because] = [false, 'below_ratio'];
        else if (sigma == null || sigma < MODULATION_SIGMA) [ordered, because] = [false, 'below_sigma'];
        else [ordered, because] = [true, 'ordered'];
    }
    Object.assign(point, { sigma, ordered, because, bar: { ratio: MODULATION_RATIO, sigma: MODULATION_SIGMA } });
};

// The point with the largest σ (null as −∞), the first on a tie: the nearest miss.
const nearestPoint = (points) => {
    let best = null;
    let bestKey = -Infinity;
    points.forEach((point) => {
        const key = point.sigma ?? -Infinity;
        if (best == null || key > bestKey) {
            best = point;
            bestKey = key;
        }
    });
    return best == null ? null : { q: [...best.q], ratio: best.ratio, sigma: best.sigma, because: best.because };
};

const modulationFinding = (answer, because, decided = [], nearest = null) => ({
    answer, because, decidedBy: decided.map((point) => [...point.q]), nearest
});

const elementVerdict = (entry) => {
    const zone = entry.zoneBoundary;
    let cellDoubling;
    if (!zone.length) cellDoubling = modulationFinding(null, 'untestable');
    else {
        const ordered = zone.filter((point) => point.ordered === true);
        if (ordered.length) cellDoubling = modulationFinding(true, 'ordered_point', ordered);
        else if (zone.every((point) => point.because === 'undefined_statistic')) cellDoubling = modulationFinding(null, 'undefined_statistic', [], nearestPoint(zone));
        else cellDoubling = modulationFinding(false, 'no_ordered_point', [], nearestPoint(zone));
    }
    const general = entry.strongest.filter((mode) => mode.kind === 'general');
    const orderedModes = general.filter((mode) => mode.ordered === true);
    let generalOrder;
    if (!general.length) generalOrder = modulationFinding(null, 'no_general_mode');
    else if (orderedModes.length) generalOrder = modulationFinding(true, 'ordered_mode', orderedModes);
    else if (general.every((mode) => mode.because === 'undefined_statistic')) generalOrder = modulationFinding(null, 'undefined_statistic', [], nearestPoint(general));
    else generalOrder = modulationFinding(false, 'no_ordered_mode', [], nearestPoint(general));
    const first = entry.strongest[0];
    let longWave;
    if (!first || first.ratio == null) longWave = modulationFinding(null, 'undefined_statistic');
    else if (first.kind !== 'long_wave') longWave = modulationFinding(false, 'strongest_not_long_wave');
    else if (first.ratio >= LONG_WAVE_RATIO) longWave = modulationFinding(true, 'long_wave_dominates', [first]);
    else longWave = modulationFinding(false, 'below_ratio');
    return { cellDoubling, generalOrder, longWave };
};

// An element with no displacement to analyse (zero spread everywhere).
const untestedVerdict = () => Object.fromEntries(
    ['cellDoubling', 'generalOrder', 'longWave'].map((name) => [name, modulationFinding(null, 'no_displacements')])
);

// Every tested point's σ, ordered, because and bar, and each element's
// verdict, over the family of this call: the zone-boundary points and the
// compared modes of every element tested.
const decide = (result) => {
    const tested = result.elements.filter((entry) => entry.sites);
    const family = {
        zoneBoundaryPoints: tested.reduce((sum, entry) => sum + entry.zoneBoundary.length, 0),
        comparedModes: tested.reduce((sum, entry) => sum + entry.comparedModes, 0)
    };
    tested.forEach((entry) => {
        entry.zoneBoundary.forEach((point) => {
            point.kind = 'zone_boundary';
            testPoint(point, 'zone', entry, family.zoneBoundaryPoints);
        });
        entry.strongest.forEach((mode) => {
            if (mode.kind === 'general') testPoint(mode, 'general', entry, family.comparedModes);
            else Object.assign(mode, { sigma: null, ordered: null, because: mode.kind, bar: null });
        });
        entry.verdict = elementVerdict(entry);
    });
    result.family = family;
    result.rules = { ...MODULATION_RULES };
};
