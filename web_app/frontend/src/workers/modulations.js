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
// floating-point noise. Every reported point carries its polarisation (plan
// item 4.3): the leading axis of the element's power tensor at q, as a
// Cartesian direction and a direct-lattice [u v w], with its share of the
// power (modePolarisation). The engine also decides what these show (the
// verdict, modulations.py's): each tested point's calibrated σ, kind,
// ordered, because and bar, each element's cellDoubling / generalOrder /
// longWave / diffuseOrder, the family (every element of the configuration,
// whichever the call names: F3) and the rules — which the Copilot's
// analyze_modulations only renders.
//
// The null (F1, 2026-10-19): the configuration's own cells placed at random,
// every correlation inside a cell kept (independentCellNull); the power is a
// weighted sum of χ²₁ (zone boundary) or χ²₂/2 (general q) terms, its tail a
// saddlepoint (nullTail). diffuseOrder (plan item 4.4): the ranks of the
// power over the compared modes regressed on the 13 neighbour-cell lags
// (diffuseOrder): short-range order between cells spread over many q.

import { normalQuantile, regularizedGamma } from './orientation.js';
import { eigenDecomposition } from './pcaKde.js';
import { evidenceClears } from './clusters.js';

export const TOP_MODES = 5;
export const PERCENTILE_RTOL = 1e-9;
// A mode with every component within this of the zone centre is a
// long-wavelength wave, left out of the comparison (modulations.py).
export const LONG_WAVE_Q = 0.25;
// The long-wave test's tolerance: a grid component k/N is 1/4 exactly or at
// least 1/(4N) from it, so on any grid it changes nothing.
export const LONG_WAVE_TOL = 1e-12;
// A polarisation is an axis: it is given in the sense that makes the closest
// crystal axis acute (crystalOrientationRows' rule); axes whose |cos| lies
// within this relative tolerance of the largest tie, and the first decides.
export const AXIS_TIE_RTOL = 1e-9;

// --- the verdict's rules (modulations.py; keep the two in sync) ------------------
//
// A zone-boundary point is ordered displacement when its power is at least
// MODULATION_RATIO times the median over all wavevectors AND it is significant
// at MODULATION_SIGMA after a family-wise (Šidák) correction over every
// zone-boundary point of every element the configuration has, whichever the
// call names (F3, 2026-10-19: choosing the element after seeing the data is
// selective inference, so the family never shrinks with the request). The
// null is the independent-cell null (F1, 2026-10-19): the element's power at q
// is a weighted sum of χ²₁ (a half-integer q, a real Fourier sum) or χ²₂/2 (a
// general q) terms whose weights are the eigenvalues of its cells' covariance
// N/(N − 1) Σ_c y(c) y(c)ᵀ normalised as the power (independentCellNull), with
// a Lugannani–Rice tail (nullTail). Until F1 the sites were taken as
// independent (P ≥ χ²_m/m, χ²_2m/2m at a general q): anti-conservative where
// they share a cell-level variable, conservative where they do not. General
// modes are corrected over the compared modes.

/** Ratio to the median power an ordered point needs (the effect size). */
export const MODULATION_RATIO = 3;
/** Family-wise one-sided significance (σ) an ordered point needs. */
export const MODULATION_SIGMA = 3;
/** An element's strongest mode this many times the median, and a long wave, is a box-scale wave. */
export const LONG_WAVE_RATIO = 5;
/** σ = +∞ (an upper tail that underflows) is reported as this; no finite double reaches it (≈ 38.5 at p = 5e-324). σ = −∞ is null. */
export const SIGMA_MAX = 40;
const P_CAP = 1 - 1e-16;

// --- the independent-cell null (F1; modulations.py's constants) -------------------

/** A null weight below this fraction of the largest is a direction no site moves in (or rounding): left out. */
export const NULL_WEIGHT_RTOL = 1e-12;
/** Weights within this relative spread are one weight: the null is then exactly a gamma. */
export const NULL_EQUAL_RTOL = 1e-12;
/** Bisection steps for the saddlepoint equation K'(t) = x (the same in both engines). */
export const SADDLEPOINT_STEPS = 200;
/** Within this |w| of the null mean the Lugannani–Rice formula cancels: its limit there. */
export const SADDLEPOINT_W_MIN = 1e-6;

// --- diffuse order between cells (plan item 4.4) ----------------------------------

/** The neighbour-cell lags (in cells) the diffuse test reads: one of each ±R pair of the 26 cells around a cell. */
export const DIFFUSE_LAGS = Object.freeze([
    [1, 0, 0], [0, 1, 0], [0, 0, 1],
    [1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1],
    [1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1]
].map((lag) => Object.freeze(lag)));
/** Family-wise one-sided σ diffuse order needs (Šidák over every element). */
export const DIFFUSE_SIGMA = 3;
/** The effect size it needs: the strongest neighbour-cell correlation |Γ(R)|. */
export const DIFFUSE_CORRELATION = 0.1;
/** Compared modes needed per term of the regression (the constant and every lag the grid resolves). */
export const DIFFUSE_MODES_PER_TERM = 2;
/** A lag whose template is, within this relative norm, a combination of the earlier ones is left out. */
export const DIFFUSE_RANK_RTOL = 1e-9;

export const MODULATION_RULES = Object.freeze({
    ratio: MODULATION_RATIO,
    sigma: MODULATION_SIGMA,
    longWaveRatio: LONG_WAVE_RATIO,
    longWaveQ: LONG_WAVE_Q,
    sigmaMax: SIGMA_MAX,
    diffuseSigma: DIFFUSE_SIGMA,
    diffuseCorrelation: DIFFUSE_CORRELATION,
    diffuseLags: DIFFUSE_LAGS.map((lag) => [...lag]),
    diffuseModesPerTerm: DIFFUSE_MODES_PER_TERM
});
/** Every `because` a point and each element finding can give (closed sets, pinned by the parity fixture). */
export const MODULATION_BECAUSE = Object.freeze({
    point: ['ordered', 'below_ratio', 'below_sigma', 'undefined_statistic', 'zone_boundary', 'long_wave'],
    cellDoubling: ['ordered_point', 'no_ordered_point', 'undefined_statistic', 'untestable', 'no_displacements'],
    generalOrder: ['ordered_mode', 'no_ordered_mode', 'undefined_statistic', 'no_general_mode', 'no_displacements'],
    longWave: ['long_wave_dominates', 'below_ratio', 'strongest_not_long_wave', 'undefined_statistic', 'no_displacements'],
    diffuseOrder: ['correlated_in_phase', 'correlated_antiphase', 'below_correlation', 'below_sigma', 'undefined_statistic',
        'too_few_modes', 'no_displacements']
});
export const MODE_KINDS = Object.freeze(['zone_boundary', 'long_wave', 'general']);
/** The phase of the strongest neighbour-cell correlation: in phase (Γ > 0) or in antiphase (Γ < 0). */
export const DIFFUSE_PHASES = Object.freeze(['in_phase', 'antiphase']);

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

// The null's weights as an array: an integer m stands for m equal weights 1/m
// (the bound P ≥ χ²_m/m used until F1). Mirrors _null_weights.
const nullWeightsOf = (weights) => {
    if (typeof weights === 'number') return Number.isInteger(weights) && weights > 0 ? Array(weights).fill(1 / weights) : [];
    return Array.from(weights ?? [], Number);
};

// P(Z > w) for a standard normal Z, accurate in both tails (modulations.py _normal_upper).
const normalUpper = (w) => {
    const half = 0.5 * regularizedGamma(0.5, 0.5 * w * w).upper;
    return w >= 0 ? half : 1 - half;
};

// Lugannani–Rice upper tail of X with K(t) = −a Σ log(1 − b_k t) (b_k > 0):
// a = ½, b = 2λ for Σ λ χ²₁; a = 1, b = λ for Σ λ χ²₂/2 (modulations.py
// _saddlepoint_tail, the same bisection and the same limit near the mean).
const saddlepointTail = (x, b, a) => {
    let sum = 0;
    let largest = 0;
    b.forEach((value) => { sum += value; largest = Math.max(largest, value); });
    const mean = a * sum;
    let [lo, hi] = x > mean ? [0, 1 / largest] : [(-a * b.length) / x, 0];
    for (let step = 0; step < SADDLEPOINT_STEPS; step += 1) {
        const mid = 0.5 * (lo + hi);
        if (!(lo < mid && mid < hi)) break;
        let slope = 0;
        b.forEach((value) => { slope += value / (1 - value * mid); });
        if (a * slope < x) lo = mid;
        else hi = mid;
    }
    const t = x > mean ? lo : 0.5 * (lo + hi);
    let cgf = 0;
    b.forEach((value) => { cgf -= Math.log1p(-value * t); });
    cgf *= a;
    const w = Math.sign(t) * Math.sqrt(Math.max(2 * (t * x - cgf), 0));
    if (Math.abs(w) < SADDLEPOINT_W_MIN) {
        let k2 = 0;
        let k3 = 0;
        b.forEach((value) => { k2 += value * value; k3 += value * value * value; });
        return 0.5 - (2 * a * k3) / (6 * Math.sqrt(2 * Math.PI) * (a * k2) ** 1.5);
    }
    let curvature = 0;
    b.forEach((value) => { curvature += (value / (1 - value * t)) ** 2; });
    const u = t * Math.sqrt(a * curvature);
    const tail = normalUpper(w) + (Math.exp(-0.5 * w * w) / Math.sqrt(2 * Math.PI)) * (1 / u - 1 / w);
    return Math.min(Math.max(tail, 0), 1);
};

/**
 * P(P_null ≥ power) under the independent-cell null with `weights` (the
 * eigenvalues of independentCellNull, or an integer m for m equal weights
 * 1/m); `real` at a zone-boundary point (Σ λ χ²₁), else a general q
 * (Σ λ χ²₂/2). One weight value: the exact gamma; several: the saddlepoint.
 * Null where undefined. Mirrors null_tail.
 */
export const nullTail = (power, weights, real) => {
    const lam = nullWeightsOf(weights).filter((value) => value > 0);
    if (!(power > 0) || !lam.length) return null;
    const [a, scale] = real ? [0.5, 2] : [1, 1];
    const largest = Math.max(...lam);
    const smallest = Math.min(...lam);
    if (largest - smallest <= NULL_EQUAL_RTOL * largest) {
        const mean = lam.reduce((sum, value) => sum + value, 0) / lam.length;
        return regularizedGamma(a * lam.length, power / (scale * mean)).upper;
    }
    return saddlepointTail(power, lam.map((value) => scale * value), a);
};

/**
 * The family-wise significance (one-sided σ) of an element's zone-boundary
 * power against the independent-cell null with `weights` (nullTail, real),
 * among `family` points tested. Null where undefined (power ≤ 0) or −∞; +∞
 * is SIGMA_MAX. Mirrors zone_boundary_sigma.
 */
export const zoneBoundarySigma = (power, weights, family) => {
    const p = nullTail(power, weights, true);
    return p == null ? null : reported(sidakSigma(p, family));
};

/** The same at a general wavevector (a complex sum), `family` = the compared modes. Mirrors general_mode_sigma. */
export const generalModeSigma = (power, weights, family) => {
    if (!(family >= 1)) return null;
    const p = nullTail(power, weights, false);
    return p == null ? null : reported(sidakSigma(p, family));
};

// Cyclic Jacobi eigenvalues of a symmetric n × n matrix (Float64Array, row
// major), largest first: the role numpy's eigvalsh has in the engine.
const symmetricEigenvalues = (matrix, n) => {
    const a = Float64Array.from(matrix);
    let scale = 0;
    for (let i = 0; i < n * n; i += 1) scale += a[i] * a[i];
    for (let sweep = 0; sweep < 100; sweep += 1) {
        let off = 0;
        for (let p = 0; p < n; p += 1) for (let q = p + 1; q < n; q += 1) off += a[p * n + q] * a[p * n + q];
        if (!(off > 1e-30 * scale)) break;
        for (let p = 0; p < n; p += 1) {
            for (let q = p + 1; q < n; q += 1) {
                const apq = a[p * n + q];
                if (apq === 0) continue;
                const theta = (a[q * n + q] - a[p * n + p]) / (2 * apq);
                const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
                const c = 1 / Math.sqrt(t * t + 1);
                const s = t * c;
                for (let k = 0; k < n; k += 1) {
                    const akp = a[k * n + p];
                    const akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for (let k = 0; k < n; k += 1) {
                    const apk = a[p * n + k];
                    const aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
            }
        }
    }
    return Array.from({ length: n }, (_, i) => a[i * n + i]).sort((x, y) => y - x);
};

// Mᵀx = d for the rows a, b, c of `cell` (Cramer's rule): the direct-lattice
// components of a Cartesian vector (np.linalg.solve(cell.T, d) in the engine).
const solveTransposed = (cell, d) => {
    const m = [0, 1, 2].map((i) => [0, 1, 2].map((j) => cell[j][i]));
    const det3 = (a) => a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    const det = det3(m);
    return [0, 1, 2].map((column) => det3(m.map((row, i) => row.map((value, j) => (j === column ? d[i] : value)))) / det);
};

/**
 * The displacement direction of a mode from its power tensor at q (a 3×3
 * array): the real part of the element's cross-spectral tensor, normalised
 * as the power (its trace). Returns { direction (the leading eigenvector, a
 * Cartesian unit vector), uvw (the same as a direct-lattice [u v w], M⁻ᵀd
 * scaled so its largest component is 1), share (λ₁ / tr) }, both in the sense
 * that makes the closest crystal axis acute; null without power. Mirrors
 * mode_polarisation.
 */
export const modePolarisation = (tensor, cell) => {
    const trace = tensor[0][0] + tensor[1][1] + tensor[2][2];
    if (!(Number.isFinite(trace) && trace > 0)) return null;
    const { eigenvalues, axes } = eigenDecomposition(tensor);
    const length = Math.hypot(...axes[0]);
    let direction = axes[0].map((value) => value / length);
    const cosines = cell.map((edge) => (edge[0] * direction[0] + edge[1] * direction[1] + edge[2] * direction[2]) / Math.hypot(...edge));
    const largest = Math.max(...cosines.map(Math.abs));
    const dominant = cosines.findIndex((value) => Math.abs(value) >= largest * (1 - AXIS_TIE_RTOL));
    if (cosines[dominant] < 0) direction = direction.map((value) => -value);
    const raw = solveTransposed(cell, direction);
    const scale = Math.max(...raw.map(Math.abs));
    // + 0: a −0 component reads 0, as the engine's + 0.0.
    return {
        direction: direction.map((value) => value + 0),
        uvw: raw.map((value) => value / scale + 0),
        share: Math.min(eigenvalues[0] / trace, 1)
    };
};

// A site's total displacement variance tr C (population); 0 for a site that does not move.
const siteVariance = (u) => {
    const n = u.length;
    const mean = [0, 1, 2].map((k) => u.reduce((sum, row) => sum + row[k], 0) / n);
    return [0, 1, 2].reduce((sum, k) => sum + u.reduce((acc, row) => acc + (row[k] - mean[k]) ** 2, 0) / n, 0);
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

// A grid index [i, j, k] (any integers) → its C-order linear index on the supercell grid.
const gridLinear = (supercell) => ([i, j, k]) => {
    const w = [i, j, k].map((v, axis) => ((v % supercell[axis]) + supercell[axis]) % supercell[axis]);
    return (w[0] * supercell[1] + w[1]) * supercell[2] + w[2];
};

// The tensor's six components per q, in this order: xx, yy, zz, xy, xz, yz.
const TENSOR_PAIRS = [[0, 0], [1, 1], [2, 2], [0, 1], [0, 2], [1, 2]];

/** The 3×3 power tensor at linear index q of elementSpectrum's `tensor`. */
export const tensorAt = (tensor, q) => {
    const [xx, yy, zz, xy, xz, yz] = tensor.subarray(6 * q, 6 * q + 6);
    return [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]];
};

/**
 * The weights of the independent-cell null of one element's power (F1):
 * the eigenvalues, largest first, of W^½ (N/(N − 1) Σ_c y(c) y(c)ᵀ) W^½, with
 * y(c) the displacements of the element's m moving sites in cell c (a
 * missing copy is a zero) and W = diag(1 / (m n_s tr C_s)); their sum is the
 * mean power over q ≠ 0, N/(N − 1). Weights below NULL_WEIGHT_RTOL of the
 * largest are left out; empty when no site moves. Mirrors independent_cell_null.
 */
export const independentCellNull = (parsed, name) => {
    const supercell = parsed.supercell.map(Number);
    const nQ = supercell[0] * supercell[1] * supercell[2];
    const linear = gridLinear(supercell);
    const moving = parsed.sites.filter((site) => site.element === name && siteVariance(site.displacements) > 0);
    if (!moving.length) return [];
    const width = 3 * moving.length;
    // z: the cells' scaled displacements, nQ rows of width columns.
    const z = new Float64Array(nQ * width);
    moving.forEach((site, s) => {
        const u = site.displacements;
        const scale = Math.sqrt(1 / (moving.length * u.length * siteVariance(u)));
        site.cells.forEach((cell, row) => {
            const base = linear(cell) * width + 3 * s;
            for (let k = 0; k < 3; k += 1) z[base + k] = u[row][k] * scale;
        });
    });
    const factor = nQ / (nQ - 1);
    const covariance = new Float64Array(width * width);
    for (let c = 0; c < nQ; c += 1) {
        const offset = c * width;
        for (let i = 0; i < width; i += 1) {
            const zi = z[offset + i];
            if (zi === 0) continue;
            for (let j = i; j < width; j += 1) covariance[i * width + j] += zi * z[offset + j];
        }
    }
    for (let i = 0; i < width; i += 1) {
        for (let j = i; j < width; j += 1) {
            covariance[i * width + j] *= factor;
            covariance[j * width + i] = covariance[i * width + j];
        }
    }
    const values = symmetricEigenvalues(covariance, width);
    return values.filter((value) => value > values[0] * NULL_WEIGHT_RTOL);
};

/**
 * One element's displacement spectrum over the supercell's cell grid, from a
 * parsed configuration: { power (Float64Array, the power averaged over the
 * element's sites that move, C order), tensor (Float64Array, 6 components
 * per q: the real part of their cross-spectral tensor normalised the same
 * way; tensorAt reads one), missing, used }. Mirrors element_spectrum.
 */
export const elementSpectrum = (parsed, name) => {
    const supercell = parsed.supercell.map(Number);
    const nQ = supercell[0] * supercell[1] * supercell[2];
    const linear = gridLinear(supercell);
    const power = new Float64Array(nQ);
    const tensor = new Float64Array(6 * nQ);
    let missing = 0;
    let used = 0;
    parsed.sites.forEach((site) => {
        if (site.element !== name) return;
        const u = site.displacements;
        const n = u.length;
        const variance = siteVariance(u);
        if (!(variance > 0)) return;
        missing += nQ - n;
        const transforms = [0, 1, 2].map((component) => {
            const re = new Float64Array(nQ);
            const im = new Float64Array(nQ);
            site.cells.forEach((cell, row) => { re[linear(cell)] = u[row][component]; });
            for (let axis = 0; axis < 3; axis += 1) dftAxis(re, im, supercell, axis);
            return { re, im };
        });
        transforms.forEach(({ re, im }) => {
            for (let q = 0; q < nQ; q += 1) power[q] += (re[q] * re[q] + im[q] * im[q]) / (n * variance);
        });
        // Re(U_j U_k*) = re_j re_k + im_j im_k.
        for (let q = 0; q < nQ; q += 1) {
            TENSOR_PAIRS.forEach(([j, k], slot) => {
                tensor[6 * q + slot] += (transforms[j].re[q] * transforms[k].re[q] + transforms[j].im[q] * transforms[k].im[q]) / (n * variance);
            });
        }
        used += 1;
    });
    if (used) {
        for (let q = 0; q < nQ; q += 1) power[q] /= used;
        for (let i = 0; i < tensor.length; i += 1) tensor[i] /= used;
    }
    return { power, tensor, missing, used };
};

/**
 * Per-element displacement power at every wavevector of the supercell grid,
 * from a parsed configuration (siteDisplacementsFromRmc6f). Mirrors
 * rmc_toolkits.modulations.displacement_modulations; camelCase payload. The
 * significance family is every element that moves, whichever `element`
 * names (F3).
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
    const every = [...new Set(labels)].sort();
    let elements;
    if (element != null) {
        const symbol = capitalize(String(element).trim());
        if (!labels.includes(symbol)) {
            const species = [...new Set(sites.flatMap((site) => Object.keys(site.elementCounts || { [site.element]: 1 })))].sort();
            throw new Error(`no site holds ${symbol}; elements: ${species.join(', ')}`);
        }
        elements = [symbol];
    } else {
        elements = every;
    }
    // The family (F3): every element with a site that moves, as a call
    // without `element` tests them.
    const familyElements = every.filter((name) => sites.some((site) => site.element === name && siteVariance(site.displacements) > 0));
    const unitCell = parsed.latticeVectors.map((row, axis) => row.map((value) => value / supercell[axis]));
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
    const linear = gridLinear(supercell);
    const unravel = (q) => [Math.floor(q / (supercell[1] * supercell[2])), Math.floor(q / supercell[2]) % supercell[1], q % supercell[2]];
    // The compared modes: one per ±q pair (the smaller grid index), q ≠ 0, no long wave.
    const compared = [];
    for (let q = 1; q < nQ; q += 1) {
        const index = unravel(q);
        const partner = linear(index.map((v) => -v));
        if (q <= partner && !isLongWave(index.map((v, axis) => folded(v, supercell[axis])))) compared.push(q);
    }
    // The diffuse test's modes: the compared modes that are not zone-boundary points.
    const diffuseModes = compared.filter((q) => !unravel(q).every((v, axis) => v === 0 || 2 * v === supercell[axis]));
    const diffuseQ = diffuseModes.map((q) => unravel(q).map((v, axis) => folded(v, supercell[axis])));
    const spectra = new Map();
    for (const name of elements) {
        const { power, tensor, missing, used } = elementSpectrum(parsed, name);
        if (!used) {
            result.elements.push({ element: name, sites: 0, verdict: untestedVerdict() });
            continue;
        }
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
                rank: reference.length ? higher + 1 : null,
                polarisation: modePolarisation(tensorAt(tensor, q), unitCell)
            };
        };
        const zoneBoundary = boundary.map((flags) => {
            const entry = describe(linear(flags.map((flag, axis) => (flag ? supercell[axis] / 2 : 0))));
            entry.q = flags.map((flag) => 0.5 * flag);
            return entry;
        });
        const order = Array.from({ length: nQ }, (_, q) => q).sort((a, b) => (power[b] - power[a]) || (a - b));
        const strongest = [];
        const strongestKeys = [];
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
            strongestKeys.push(key);
            if (strongest.length === top) break;
        }
        const weights = independentCellNull(parsed, name);
        const weightSum = weights.reduce((sum, value) => sum + value, 0);
        const weightSquares = weights.reduce((sum, value) => sum + value * value, 0);
        result.elements.push({
            element: name,
            sites: used,
            missingCopies: missing,
            medianPower: median,
            comparedModes: reference.length,
            // The null's effective number of independent components, (Σλ)²/Σλ².
            nullDof: (weightSum * weightSum) / weightSquares,
            zoneBoundary,
            strongest
        });
        spectra.set(name, { weights, power, strongest: strongestKeys });
    }
    decide(result, familyElements, boundary.length, compared.length, spectra, { modes: diffuseModes, q: diffuseQ });
    return result;
};

// --- the verdict (modulations.py _decide; keep the two in sync) ------------------

// Ranks 1..K of `values` in ascending order; a run within PERCENTILE_RTOL of
// its smallest shares the run's mean rank (modulations.py _average_ranks).
const averageRanks = (values) => {
    const order = Array.from(values, (_, i) => i).sort((i, j) => (values[i] - values[j]) || (i - j));
    const ranks = new Float64Array(values.length);
    let start = 0;
    while (start < values.length) {
        const first = values[order[start]];
        let end = start + 1;
        while (end < values.length && values[order[end]] <= first + PERCENTILE_RTOL * Math.abs(first)) end += 1;
        for (let k = start; k < end; k += 1) ranks[order[k]] = 0.5 * (start + 1 + end);
        start = end;
    }
    return ranks;
};

const dot = (a, b) => {
    let sum = 0;
    for (let i = 0; i < a.length; i += 1) sum += a[i] * b[i];
    return sum;
};

// Modified Gram–Schmidt over `columns` in order: { basis, kept }, a column
// left out when its part orthogonal to the earlier ones is within
// DIFFUSE_RANK_RTOL of its norm (modulations.py _orthonormal_columns).
const orthonormalColumns = (columns) => {
    const basis = [];
    const kept = [];
    columns.forEach((column, index) => {
        const vector = Float64Array.from(column);
        const norm = Math.sqrt(dot(vector, vector));
        basis.forEach((unit) => {
            const along = dot(unit, vector);
            for (let i = 0; i < vector.length; i += 1) vector[i] -= along * unit[i];
        });
        const rest = Math.sqrt(dot(vector, vector));
        if (norm > 0 && rest > DIFFUSE_RANK_RTOL * norm) {
            basis.push(vector.map((value) => value / rest));
            kept.push(index);
        }
    });
    return { basis, kept };
};

// A small dense system by Gaussian elimination with partial pivoting (modulations.py _solve_dense).
const solveDense = (matrix, rhs) => {
    const a = matrix.map((row) => [...row]);
    const b = [...rhs];
    const n = b.length;
    for (let column = 0; column < n; column += 1) {
        let pivot = column;
        for (let row = column + 1; row < n; row += 1) if (Math.abs(a[row][column]) > Math.abs(a[pivot][column])) pivot = row;
        if (pivot !== column) {
            [a[column], a[pivot]] = [a[pivot], a[column]];
            [b[column], b[pivot]] = [b[pivot], b[column]];
        }
        for (let row = column + 1; row < n; row += 1) {
            const factor = a[row][column] / a[column][column];
            for (let k = column; k < n; k += 1) a[row][k] -= factor * a[column][k];
            b[row] -= factor * b[column];
        }
    }
    const x = new Array(n).fill(0);
    for (let row = n - 1; row >= 0; row -= 1) {
        let sum = 0;
        for (let k = row + 1; k < n; k += 1) sum += a[row][k] * x[k];
        x[row] = (b[row] - sum) / a[row][row];
    }
    return x;
};

/**
 * The diffuse test (plan item 4.4): the normal scores of the ranks of the
 * power over the modes compared (one per ±q pair, no long wave, no
 * zone-boundary point, none found ordered) regressed on the neighbour-cell
 * lags' templates cos(2π q·R) (with the constant; a lag the grid does not
 * resolve is left out): χ² = (K − 1) |P a|² / |a|², χ² with r degrees of
 * freedom under independent cells, Šidák over `family` elements. `power`
 * holds those modes' powers, `q` their fractional wavevectors. Returns
 * { modes, df, chiSquare, sigma, correlations (Γ per DIFFUSE_LAGS from the
 * least-squares fit of the power, β_R / (2 · mean power)), lag, correlation,
 * phase }. Mirrors diffuse_order.
 */
export const diffuseOrder = (power, q, family) => {
    const modes = power.length;
    const out = {
        modes, df: null, chiSquare: null, sigma: null, correlations: DIFFUSE_LAGS.map(() => null), lag: null, correlation: null, phase: null
    };
    if (!modes) return out;
    const columns = [Array(modes).fill(1), ...DIFFUSE_LAGS.map((lag) => q.map((v) => Math.cos(2 * Math.PI * (v[0] * lag[0] + v[1] * lag[1] + v[2] * lag[2]))))];
    const { basis, kept } = orthonormalColumns(columns);
    const lags = kept.filter((index) => index > 0).map((index) => index - 1);
    if (!lags.length || kept[0] !== 0 || modes < DIFFUSE_MODES_PER_TERM * (lags.length + 1)) return out;
    out.df = lags.length;
    const design = kept.map((index) => columns[index]);
    const gram = design.map((a) => design.map((b) => dot(a, b)));
    const beta = solveDense(gram, design.map((a) => dot(a, power)));
    const mean = power.reduce((sum, value) => sum + value, 0) / modes;
    if (mean > 0) lags.forEach((lag, position) => { out.correlations[lag] = beta[position + 1] / (2 * mean); });
    const ranks = averageRanks(power);
    const raw = Array.from(ranks, (rank) => normalQuantile((rank - 0.5) / modes));
    const centre = raw.reduce((sum, value) => sum + value, 0) / modes;
    const scores = raw.map((value) => value - centre);
    const total = dot(scores, scores);
    if (!(total > 0) || !(mean > 0)) return out;
    let projection = 0;
    basis.slice(1).forEach((unit) => { projection += dot(unit, scores) ** 2; });
    const chiSquare = ((modes - 1) * projection) / total;
    out.chiSquare = chiSquare;
    out.sigma = reported(sidakSigma(regularizedGamma(lags.length / 2, chiSquare / 2).upper, family));
    // The strongest correlation: the largest |Γ|, the first lag on a tie.
    const magnitude = lags.map((lag) => Math.abs(out.correlations[lag]));
    const largest = Math.max(...magnitude);
    const strongest = lags[magnitude.findIndex((value) => value >= largest * (1 - PERCENTILE_RTOL))];
    out.lag = [...DIFFUSE_LAGS[strongest]];
    out.correlation = out.correlations[strongest];
    out.phase = out.correlation > 0 ? 'in_phase' : 'antiphase';
    return out;
};

// diffuseOrder's finding: answer, because and one Evidence record (the
// cluster engine's form, evidenceClears); the answer is the record's clears.
// Mirrors _diffuse_finding.
const diffuseFinding = (diffuse) => {
    if (diffuse.df == null) return { answer: null, because: 'too_few_modes', evidence: [] };
    const magnitude = diffuse.correlation != null ? Math.abs(diffuse.correlation) : null;
    const record = {
        statistic: 'neighbourCellCorrelation',
        against: 'independentCells',
        lag: diffuse.lag,
        correlation: diffuse.correlation,
        phase: diffuse.phase,
        chiSquare: diffuse.chiSquare,
        df: diffuse.df,
        modes: diffuse.modes,
        observed: magnitude,
        null: 0,
        z: diffuse.sigma,
        change: magnitude,
        changeKind: 'difference',
        bar: { z: DIFFUSE_SIGMA, change: DIFFUSE_CORRELATION },
        side: 'above'
    };
    record.clears = evidenceClears(record);
    if (diffuse.chiSquare == null || magnitude == null) return { answer: null, because: 'undefined_statistic', evidence: [record] };
    let because;
    if (record.clears) because = diffuse.phase === 'in_phase' ? 'correlated_in_phase' : 'correlated_antiphase';
    else if (diffuse.sigma == null || diffuse.sigma < DIFFUSE_SIGMA) because = 'below_sigma';
    else because = 'below_correlation';
    return { answer: record.clears, because, evidence: [record] };
};

// One point against its bar and the element's independent-cell null
// (`weights`). The power is ratio × median, as the tool formed it: σ far
// below the background is ill-conditioned, and this keeps it bit for bit the
// value the tool reported.
const testPoint = (point, test, entry, family, weights) => {
    const { ratio } = point;
    const power = ratio != null ? ratio * entry.medianPower : null;
    const isDefined = power != null && power > 0 && weights.length > 0 && (test === 'zone' || family >= 1);
    let sigma = null;
    let ordered;
    let because;
    if (!isDefined) {
        ordered = null;
        because = 'undefined_statistic';
    } else {
        sigma = test === 'zone' ? zoneBoundarySigma(power, weights, family) : generalModeSigma(power, weights, family);
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
    return { cellDoubling, generalOrder, longWave, diffuseOrder: diffuseFinding(entry.diffuse) };
};

// An element with no displacement to analyse (zero spread everywhere).
const untestedVerdict = () => ({
    ...Object.fromEntries(['cellDoubling', 'generalOrder', 'longWave'].map((name) => [name, modulationFinding(null, 'no_displacements')])),
    diffuseOrder: { answer: null, because: 'no_displacements', evidence: [] }
});

// Every tested point's σ, ordered, because and bar, and each element's
// verdict, over the family: the `points` zone-boundary points and `modes`
// compared modes of each of `familyElements` — every element the
// configuration has with displacements, whichever the call reports (F3) —
// and one diffuse test per element. `spectra` holds each tested element's
// null weights, power and strongest modes' grid indices; `diffuse` the
// diffuse test's modes and their wavevectors.
const decide = (result, familyElements, points, modes, spectra, diffuse) => {
    const tested = result.elements.filter((entry) => entry.sites);
    const family = {
        zoneBoundaryPoints: points * familyElements.length,
        comparedModes: modes * familyElements.length,
        diffuse: familyElements.length,
        elements: [...familyElements]
    };
    tested.forEach((entry) => {
        const spectrum = spectra.get(entry.element);
        entry.zoneBoundary.forEach((point) => {
            point.kind = 'zone_boundary';
            testPoint(point, 'zone', entry, family.zoneBoundaryPoints, spectrum.weights);
        });
        entry.strongest.forEach((mode) => {
            if (mode.kind === 'general') testPoint(mode, 'general', entry, family.comparedModes, spectrum.weights);
            else Object.assign(mode, { sigma: null, ordered: null, because: mode.kind, bar: null });
        });
        // The diffuse test asks what the modes found ordered leave.
        const found = new Set(spectrum.strongest.filter((_, i) => entry.strongest[i].kind === 'general' && entry.strongest[i].ordered === true));
        const keep = diffuse.modes.map((q) => !found.has(q));
        entry.diffuse = diffuseOrder(
            diffuse.modes.filter((_, i) => keep[i]).map((q) => spectrum.power[q]),
            diffuse.q.filter((_, i) => keep[i]),
            family.diffuse
        );
        entry.verdict = elementVerdict(entry);
    });
    result.family = family;
    result.rules = { ...MODULATION_RULES };
};
