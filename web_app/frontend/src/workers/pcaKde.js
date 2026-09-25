// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// PCA + KDE engine for RMC displacement clouds (browser / static mode).
//
// PCA and KDE are standard statistical tools; the specific analysis here (per-site
// PCA of RMC displacement clouds -> thermal ellipsoid + 3D KDE isosurface with
// wall projections) follows Maksim Eremenko's PCA_KDE utilities:
// https://github.com/MaximEremenko/Utilities/tree/main/RMCProfileUtilities/PCA_KDE
// We followed that approach and reimplemented it independently -- not a port of
// his KDE.js (which evaluates a full multivariate Gaussian KDE).
//
// This is the JavaScript port of `rmc_toolkits/pca_kde.py`; the Python module is
// the source of truth and documents the physics and the separability argument.
// The short version: RMCProfile tags every atom with its reference site and box
// copy, so subtracting a site's mean position leaves one displacement cloud per
// crystallographic site. Its covariance is the anisotropic displacement tensor
// (the thermal ellipsoid) and a Gaussian KDE of the cloud is the smooth density
// the ellipsoid approximates.
//
// The KDE is evaluated on a grid aligned with the cloud's principal axes. With
// scipy's covariance-scaled bandwidth H = factor^2 * C, both C and H are
// diagonal in that PCA frame, so the 3D Gaussian factorizes into three 1D
// Gaussians and the volume is a tensor product of three (grid x N) kernel
// matrices. That turns N * grid^3 exponentials into 3 * N * grid, which is what
// makes a 48^3 volume interactive in a browser worker; the contraction itself
// is the same estimator scipy.stats.gaussian_kde defines, to round-off.

import { parseAtomLine } from '../rmc6f.js';

// --- Linear algebra on 3x3 symmetric matrices --------------------------------

const EIGENVALUE_FLOOR_RATIO = 1e-8;
const DEGENERATE_RATIO = 1e-6;
// Absolute floor on the largest displacement variance (A^2): below (1e-4 A)^2 a
// site has no displacement at all (an *AVERAGE.rmc6f or ideal configuration). Same
// constant and rule as ZERO_SPREAD_VARIANCE in pca_kde.py.
export const ZERO_SPREAD_VARIANCE = 1e-8;

// --- Chi-square(3) quantile ------------------------------------------------------
// The squared Mahalanobis radius of a 3D Gaussian is chi-square with 3 degrees of
// freedom, whose CDF has the closed form F(x) = erf(sqrt(x/2)) - sqrt(2x/pi) e^(-x/2)
// = P(3/2, x/2), the regularised lower incomplete gamma function. It is evaluated
// here without cancellation -- by its positive series below t = a + 1 and by the
// continued fraction of the complement Q above (Numerical Recipes gser/gcf) -- and
// inverted by safeguarded Newton, so k(p) = sqrt(F^-1(p)) equals the server's
// sqrt(scipy.stats.chi2.ppf(p, 3)) to ~1e-15 (pinned at 1e-10 by the tests).
const GAMMA_A = 1.5;
const LOG_GAMMA_A = Math.log(Math.sqrt(Math.PI) / 2);   // ln Gamma(3/2)

// { lower: F(x), upper: 1 - F(x) }, each computed directly (never as 1 - tiny).
const chiSquare3Tails = (x) => {
    if (!(x > 0)) return { lower: 0, upper: 1 };
    const t = x / 2;
    const prefactor = Math.exp(-t + GAMMA_A * Math.log(t) - LOG_GAMMA_A);
    if (t < GAMMA_A + 1) {
        // P(a, t) = e^-t t^a / Gamma(a) * sum_k t^k / (a (a+1) ... (a+k)), all terms positive.
        let term = 1 / GAMMA_A;
        let sum = term;
        let ap = GAMMA_A;
        for (let k = 0; k < 1000; k += 1) {
            ap += 1;
            term *= t / ap;
            sum += term;
            if (term <= sum * 1e-17) break;
        }
        const lower = prefactor * sum;
        return { lower, upper: 1 - lower };
    }
    // Q(a, t) by the modified-Lentz continued fraction.
    const TINY = 1e-300;
    let b = t + 1 - GAMMA_A;
    let c = 1 / TINY;
    let d = 1 / b;
    let h = d;
    for (let i = 1; i < 1000; i += 1) {
        const an = -i * (i - GAMMA_A);
        b += 2;
        d = an * d + b;
        if (Math.abs(d) < TINY) d = TINY;
        c = b + an / c;
        if (Math.abs(c) < TINY) c = TINY;
        d = 1 / d;
        const delta = d * c;
        h *= delta;
        if (Math.abs(delta - 1) <= 1e-16) break;
    }
    const upper = prefactor * h;
    return { lower: 1 - upper, upper };
};

/** Exact chi-square(3) quantile F^-1(p), 0 < p < 1 (scipy.stats.chi2.ppf(p, 3)). */
export const chiSquare3Quantile = (probability) => {
    const p = Number(probability);
    if (!(p > 0 && p < 1)) throw new Error('probability must lie strictly between 0 and 1');
    const q = 1 - p;   // exact for p >= 1/2 (Sterbenz); used only there
    // Residual F(x) - p from whichever tail keeps full relative precision.
    const residual = (x) => {
        const { lower, upper } = chiSquare3Tails(x);
        return p <= 0.5 ? lower - p : q - upper;
    };
    let lo = 0;
    let hi = 1;
    while (residual(hi) < 0 && hi < 1e4) { lo = hi; hi *= 2; }
    let x = 0.5 * (lo + hi);
    for (let iter = 0; iter < 300; iter += 1) {
        const r = residual(x);
        if (r === 0) break;
        if (r > 0) hi = x; else lo = x;
        const density = Math.sqrt(x / (2 * Math.PI)) * Math.exp(-x / 2);
        let next = x - r / density;
        // Safeguard: fall back to bisection whenever Newton leaves the bracket.
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        const step = Math.abs(next - x);
        x = next;
        if (step <= 4e-16 * x || hi - lo <= 4e-16 * hi) break;
    }
    return x;
};

// Ellipsoid scale factor k such that k*sigma encloses `probability`:
// k = sqrt(chi2_3^-1(p)), 1.5381722 at the crystallographic 50% convention.
export const probabilityScale = (probability) => Math.sqrt(chiSquare3Quantile(probability));

// Symmetric 3x3 eigendecomposition by cyclic Jacobi rotation. Robust for the
// near-degenerate clouds (flat or linear disorder) that trip analytic formulas,
// and three iterations of a 3x3 sweep are negligible next to the KDE itself.
// The stopping test is RELATIVE to the matrix's Frobenius norm, so a matrix of
// any scale (1e-28 A^2 round-off included) is rotated to the same precision
// instead of being returned undiagonalised.
const jacobiEigenSymmetric = (matrix) => {
    const a = matrix.map((row) => row.slice());
    const v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let norm = 0;
    for (let i = 0; i < 3; i += 1) for (let j = 0; j < 3; j += 1) norm += a[i][j] * a[i][j];
    const tolerance = 1e-15 * Math.sqrt(norm);
    for (let sweep = 0; sweep < 50; sweep += 1) {
        const off = Math.abs(a[0][1]) + Math.abs(a[0][2]) + Math.abs(a[1][2]);
        if (off <= tolerance) break;
        for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
            if (Math.abs(a[p][q]) < 1e-300) continue;
            const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
            const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
            const c = 1 / Math.sqrt(t * t + 1);
            const s = t * c;
            for (let k = 0; k < 3; k += 1) {
                const akp = a[k][p];
                const akq = a[k][q];
                a[k][p] = c * akp - s * akq;
                a[k][q] = s * akp + c * akq;
            }
            for (let k = 0; k < 3; k += 1) {
                const apk = a[p][k];
                const aqk = a[q][k];
                a[p][k] = c * apk - s * aqk;
                a[q][k] = s * apk + c * aqk;
            }
            for (let k = 0; k < 3; k += 1) {
                const vkp = v[k][p];
                const vkq = v[k][q];
                v[k][p] = c * vkp - s * vkq;
                v[k][q] = s * vkp + c * vkq;
            }
        }
    }
    const values = [a[0][0], a[1][1], a[2][2]];
    const vectors = [0, 1, 2].map((col) => [v[0][col], v[1][col], v[2][col]]);
    return { values, vectors };
};

// Descending eigenvalues with sign- and handedness-canonicalized axes (rows), so
// results are reproducible: largest-magnitude component of each axis positive,
// then the third axis flipped if needed to keep a right-handed frame. Mirrors
// `_canonical_axes` / `_eigen_decomposition` in the Python engine.
export const eigenDecomposition = (covariance) => {
    const { values, vectors } = jacobiEigenSymmetric(covariance);
    const order = [0, 1, 2].sort((i, j) => values[j] - values[i]);
    const eigenvalues = order.map((i) => Math.max(values[i], 0));
    const axes = order.map((i) => {
        const axis = vectors[i].slice();
        let lead = 0;
        for (let k = 1; k < 3; k += 1) {
            if (Math.abs(axis[k]) > Math.abs(axis[lead])) lead = k;
        }
        const sign = axis[lead] < 0 ? -1 : 1;
        return axis.map((value) => value * sign);
    });
    const det = axes[0][0] * (axes[1][1] * axes[2][2] - axes[1][2] * axes[2][1])
        - axes[0][1] * (axes[1][0] * axes[2][2] - axes[1][2] * axes[2][0])
        + axes[0][2] * (axes[1][0] * axes[2][1] - axes[1][1] * axes[2][0]);
    if (det < 0) axes[2] = axes[2].map((value) => -value);
    return { eigenvalues, axes };
};

const covariance3 = (points) => {
    const n = points.length;
    const mean = [0, 0, 0];
    points.forEach((point) => {
        mean[0] += point[0];
        mean[1] += point[1];
        mean[2] += point[2];
    });
    mean[0] /= n;
    mean[1] /= n;
    mean[2] /= n;
    const cov = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    points.forEach((point) => {
        const dx = point[0] - mean[0];
        const dy = point[1] - mean[1];
        const dz = point[2] - mean[2];
        cov[0][0] += dx * dx; cov[0][1] += dx * dy; cov[0][2] += dx * dz;
        cov[1][1] += dy * dy; cov[1][2] += dy * dz; cov[2][2] += dz * dz;
    });
    const denom = Math.max(n - 1, 1);
    cov[0][0] /= denom; cov[0][1] /= denom; cov[0][2] /= denom;
    cov[1][1] /= denom; cov[1][2] /= denom; cov[2][2] /= denom;
    cov[1][0] = cov[0][1]; cov[2][0] = cov[0][2]; cov[2][1] = cov[1][2];
    return { mean, cov };
};

// Standard errors of the eigenvalue gap an axis must clear to count as resolved
// (AXIS_RESOLUTION_SIGMAS in pca_kde.py): a truly degenerate pair passes in < 1%.
const AXIS_RESOLUTION_SIGMAS = 3;

// Kurtosis of a centered cloud in its PCA frame -- the port of `_shape_statistics`:
//  - excessKurtosis[a] = m4/m2^2 - 3 along each principal axis (population
//    moments), null where undefined: no spread (lambda_1 < ZERO_SPREAD_VARIANCE)
//    or a collapsed axis (lambda_a < DEGENERATE_RATIO * lambda_1), where m2 is
//    round-off and the ratio 0/0. One rule, no floor constant.
//  - nonGaussianity: Mardia's b2 = mean[(u^T S^-1 u)^2] over the d defined axes,
//    normalised 3 (b2 - d(d+2)) / (d(d+2)) = (b2 - 15)/5 for d = 3. Affine
//    invariant (so noise in a near-isotropic site's frame cannot move it) and,
//    for any elliptical distribution, the excess kurtosis along every direction.
//  - axisResolved[a]: the eigenvalue gap to each neighbour exceeds
//    AXIS_RESOLUTION_SIGMAS standard errors, SE(lambda_a) = sqrt((m4 - m2^2)/n);
//    otherwise the axis direction (its kappa, its crystal orientation) is noise.
const shapeStatistics = (points, mean, axes, eigenvalues) => {
    const n = points.length;
    const m2 = [0, 0, 0];
    const m4 = [0, 0, 0];
    const project = (point, out) => {
        const dx = point[0] - mean[0];
        const dy = point[1] - mean[1];
        const dz = point[2] - mean[2];
        for (let a = 0; a < 3; a += 1) out[a] = dx * axes[a][0] + dy * axes[a][1] + dz * axes[a][2];
        return out;
    };
    const q = [0, 0, 0];
    points.forEach((point) => {
        project(point, q);
        for (let a = 0; a < 3; a += 1) {
            const p2 = q[a] * q[a];
            m2[a] += p2;
            m4[a] += p2 * p2;
        }
    });
    const count = Math.max(n, 1);
    for (let a = 0; a < 3; a += 1) { m2[a] /= count; m4[a] /= count; }

    const hasSpread = eigenvalues[0] >= ZERO_SPREAD_VARIANCE;
    const defined = [0, 1, 2].map((a) => hasSpread && eigenvalues[a] >= DEGENERATE_RATIO * eigenvalues[0]);
    const excessKurtosis = [0, 1, 2].map((a) => (defined[a] ? m4[a] / (m2[a] * m2[a]) - 3 : null));

    const inverse = [0, 1, 2].map((a) => (defined[a] ? 1 / m2[a] : 0));
    let b2 = 0;
    points.forEach((point) => {
        project(point, q);
        const r2 = q[0] * q[0] * inverse[0] + q[1] * q[1] * inverse[1] + q[2] * q[2] * inverse[2];
        b2 += r2 * r2;
    });
    b2 /= count;
    const d = defined.filter(Boolean).length;
    const reference = d * (d + 2);
    const nonGaussianity = d > 0 ? (3 * (b2 - reference)) / reference : null;

    const error = [0, 1, 2].map((a) => Math.sqrt(Math.max(m4[a] - m2[a] * m2[a], 0) / count));
    const gap = [0, 1].map((a) => (m2[a] - m2[a + 1]) > AXIS_RESOLUTION_SIGMAS * Math.hypot(error[a], error[a + 1]));
    const axisResolved = [gap[0], gap[0] && gap[1], gap[1]].map((value) => value && hasSpread);
    return { excessKurtosis, nonGaussianity, axisResolved };
};

// --- Bandwidth and sampling ---------------------------------------------------

const bandwidthFactor = (method, count, dimensions) => {
    if (typeof method === 'number') {
        if (!(Number.isFinite(method) && method > 0)) throw new Error('numeric bandwidth must be a positive finite number');
        return method;
    }
    const name = String(method).toLowerCase();
    if (name === 'scott') return count ** (-1 / (dimensions + 4));
    if (name === 'silverman') return (count * (dimensions + 2) / 4) ** (-1 / (dimensions + 4));
    throw new Error("bw must be 'scott', 'silverman', or a positive number");
};

const mulberry32 = (seed) => {
    let value = seed >>> 0;
    return () => {
        value += 0x6d2b79f5;
        let mixed = value;
        mixed = Math.imul(mixed ^ (mixed >>> 15), mixed | 1);
        mixed ^= mixed + Math.imul(mixed ^ (mixed >>> 7), mixed | 61);
        return ((mixed ^ (mixed >>> 14)) >>> 0) / 4294967296;
    };
};

const subsample = (points, limit, seed) => {
    if (points.length <= limit) return points;
    const random = mulberry32(seed);
    const indices = Array.from({ length: points.length }, (_, index) => index);
    for (let index = 0; index < limit; index += 1) {
        const swap = index + Math.floor(random() * (indices.length - index));
        [indices[index], indices[swap]] = [indices[swap], indices[index]];
    }
    return indices.slice(0, limit).sort((i, j) => i - j).map((index) => points[index]);
};

// Same cap as pca_kde.py. It binds for pooled clouds and for a single site in a
// box of >= 28 cells per edge; this draw (mulberry32 Fisher-Yates) differs from
// numpy's, so above the cap the two runtimes agree only statistically.
export const MAX_PCA_FIT_POINTS = 20000;

// One 1D Gaussian kernel matrix, laid out row-major as grid x N in a Float64Array
// (grid rows of N samples); `kernel[i * n + m]`.
const kernelMatrix = (axisCoords, projected, bandwidth) => {
    const grid = axisCoords.length;
    const n = projected.length;
    const kernel = new Float64Array(grid * n);
    for (let i = 0; i < grid; i += 1) {
        const gx = axisCoords[i];
        const base = i * n;
        for (let m = 0; m < n; m += 1) {
            const scaled = (gx - projected[m]) / bandwidth;
            kernel[base + m] = Math.exp(-0.5 * scaled * scaled);
        }
    }
    return kernel;
};

const isoLevels = (density, cellVolume, probabilities) => {
    const flat = Array.from(density).sort((a, b) => b - a);
    let running = 0;
    const cumulative = new Float64Array(flat.length);
    for (let i = 0; i < flat.length; i += 1) {
        running += flat[i] * cellVolume;
        cumulative[i] = running;
    }
    const mass = flat.length ? cumulative[cumulative.length - 1] : 0;

    const massLevels = [];
    if (mass > 0) {
        probabilities.forEach((p) => {
            const target = p * mass;
            let lo = 0;
            let hi = cumulative.length - 1;
            let idx = cumulative.length - 1;
            while (lo <= hi) {
                const mid = (lo + hi) >> 1;
                if (cumulative[mid] >= target) { idx = mid; hi = mid - 1; } else { lo = mid + 1; }
            }
            massLevels.push({ p, level: flat[idx] });
        });
    }

    let vmin = Infinity;
    let vmax = -Infinity;
    for (let i = 0; i < density.length; i += 1) {
        if (density[i] < vmin) vmin = density[i];
        if (density[i] > vmax) vmax = density[i];
    }
    if (!Number.isFinite(vmin)) { vmin = 0; vmax = 0; }
    const densityLevels = probabilities.map((p) => ({ p, level: vmin + p * (vmax - vmin) }));
    return { massLevels, densityLevels, mass, vmin, vmax };
};

// 2D KDE of the cloud projected onto a principal plane. The marginal of the
// separable 3D volume keeps the 3D bandwidth sub-block, so it is a single
// (grid x N) . (N x grid) product of two kernel matrices -- the honest shadow of
// the displayed volume, not an independently re-optimized 2D estimate.
const projection = (kernels, bandwidths, axisCoords, first, second, count) => {
    const kf = kernels[first];
    const ks = kernels[second];
    const gridF = axisCoords[first].length;
    const gridS = axisCoords[second].length;
    const n = count;
    const norm = 1 / (count * 2 * Math.PI * bandwidths[first] * bandwidths[second]);
    const density = [];
    let vmax = 0;
    for (let i = 0; i < gridF; i += 1) {
        const row = new Array(gridS);
        const baseF = i * n;
        for (let j = 0; j < gridS; j += 1) {
            const baseS = j * n;
            let sum = 0;
            for (let m = 0; m < n; m += 1) sum += kf[baseF + m] * ks[baseS + m];
            const value = sum * norm;
            row[j] = value;
            if (value > vmax) vmax = value;
        }
        density.push(row);
    }
    return {
        density,
        extent: [axisCoords[first][0], axisCoords[first][gridF - 1],
            axisCoords[second][0], axisCoords[second][gridS - 1]],
        axes: [first, second],
        bandwidth: [bandwidths[first], bandwidths[second]],
        vmax
    };
};

/**
 * PCA statistics and a separable 3D Gaussian KDE volume for one displacement
 * cloud. `points` is an array of [x, y, z] Cartesian (Angstrom) offsets. The
 * returned `density` is a flat Float64Array in C order over (PC1, PC2, PC3):
 * index = (i * grid + j) * grid + k. Mirrors `pca_kde_volume` in Python.
 */
export const pcaKdeVolume = (points, options = {}) => {
    const {
        bw = 'scott',
        bwScale = 1,
        grid: gridOption = 48,
        extent = 3,
        cubicBox = false,
        probability = 0.5,
        probabilities = Array.from({ length: 101 }, (_, i) => i / 100),
        projections = true,
        maxFitPoints = MAX_PCA_FIT_POINTS,
        rngSeed = 0
    } = options;

    if (!Array.isArray(points) || points.length < 4) {
        throw new Error('a 3D KDE needs at least four points');
    }
    if (!points.every((point) => Number.isFinite(point[0]) && Number.isFinite(point[1]) && Number.isFinite(point[2]))) {
        throw new Error('displacement cloud contains non-finite coordinates');
    }
    // NaN fails every `> 0` test but Infinity passes it, so each parameter is
    // checked for finiteness too (as pca_kde_volume does).
    if (!Number.isFinite(Number(gridOption))) throw new Error('grid must be a finite number');
    const grid = Math.max(8, Math.min(Math.round(gridOption), 128));
    if (!(Number.isFinite(bwScale) && bwScale > 0)) throw new Error('bwScale must be a positive finite number');
    if (!(Number.isFinite(extent) && extent > 0)) throw new Error('extent must be a positive finite number');

    const total = points.length;
    const fit = subsample(points, maxFitPoints, rngSeed);
    const count = fit.length;

    const { mean, cov } = covariance3(fit);
    const decomposition = eigenDecomposition(cov);
    const axes = decomposition.axes;
    let eigenvalues = decomposition.eigenvalues;

    // No spread at all (an average/ideal configuration): nothing to estimate.
    const largest = eigenvalues[0];
    if (!(largest >= ZERO_SPREAD_VARIANCE)) {
        throw new Error('displacement cloud has zero spread (RMS below 1e-4 A on every axis)');
    }
    // A flat direction would make the bandwidth singular; floor it, and say so.
    const ratio = eigenvalues[2] / largest;
    const rawEigenvalues = eigenvalues;
    eigenvalues = eigenvalues.map((value) => Math.max(value, largest * EIGENVALUE_FLOOR_RATIO));

    const factor = bandwidthFactor(bw, count, 3) * bwScale;
    const sigma = eigenvalues.map(Math.sqrt);
    const bandwidths = sigma.map((value) => factor * value);

    // The volume (and mass, iso levels, walls) is always sampled on the per-axis
    // box, whose nodes resolve every axis's kernel; cubicBox only sizes the display
    // box (boxHalfWidths), exactly as in pca_kde_volume.
    const broaden = Math.sqrt(1 + factor * factor);
    const halfWidths = sigma.map((value) => extent * value * broaden);
    const maxHalf = Math.max(...halfWidths);
    const boxHalfWidths = cubicBox ? [maxHalf, maxHalf, maxHalf] : halfWidths.slice();

    // Project the centered cloud onto the principal axes: projected[axis][m].
    const projected = [new Float64Array(count), new Float64Array(count), new Float64Array(count)];
    for (let m = 0; m < count; m += 1) {
        const dx = fit[m][0] - mean[0];
        const dy = fit[m][1] - mean[1];
        const dz = fit[m][2] - mean[2];
        for (let axis = 0; axis < 3; axis += 1) {
            projected[axis][m] = dx * axes[axis][0] + dy * axes[axis][1] + dz * axes[axis][2];
        }
    }

    const axisCoords = halfWidths.map((half) => {
        const line = new Float64Array(grid);
        const step = grid > 1 ? (2 * half) / (grid - 1) : 0;
        for (let i = 0; i < grid; i += 1) line[i] = -half + i * step;
        return line;
    });
    const kernels = [0, 1, 2].map((axis) => kernelMatrix(axisCoords[axis], projected[axis], bandwidths[axis]));

    // density[i,j,k] = norm * sum_m Kx[i,m] Ky[j,m] Kz[k,m], assembled one PC3
    // slice at a time so no N x grid^3 temporary is ever built. The inner loop
    // reuses a (grid x N) buffer of Kx[i,m]*Kz[k,m] shared across all j.
    const norm = 1 / (count * (2 * Math.PI) ** 1.5 * bandwidths[0] * bandwidths[1] * bandwidths[2]);
    const kx = kernels[0];
    const ky = kernels[1];
    const kz = kernels[2];
    const density = new Float64Array(grid * grid * grid);
    const xzRow = new Float64Array(count);
    for (let k = 0; k < grid; k += 1) {
        const kzBase = k * count;
        for (let i = 0; i < grid; i += 1) {
            const kxBase = i * count;
            for (let m = 0; m < count; m += 1) xzRow[m] = kx[kxBase + m] * kz[kzBase + m];
            const outBase = (i * grid) * grid + k;
            for (let j = 0; j < grid; j += 1) {
                const kyBase = j * count;
                let sum = 0;
                for (let m = 0; m < count; m += 1) sum += xzRow[m] * ky[kyBase + m];
                density[outBase + j * grid] = sum * norm;
            }
        }
    }

    const cellVolume = axisCoords.reduce((product, coords) => product * (coords[1] - coords[0]), 1);
    const { massLevels, densityLevels, mass, vmin, vmax } = isoLevels(density, cellVolume, probabilities);

    const { excessKurtosis, nonGaussianity, axisResolved } = shapeStatistics(fit, mean, axes, rawEigenvalues);
    const scale = probabilityScale(probability);
    const result = {
        count: total,
        fitCount: count,
        mean,
        covariance: cov,
        eigenvalues,
        axes,
        rms: sigma,
        semiAxes: sigma.map((value) => scale * value),
        probability,
        uIso: (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]) / 3,
        bIso: 8 * Math.PI * Math.PI * ((eigenvalues[0] + eigenvalues[1] + eigenvalues[2]) / 3),
        anisotropy: sigma[0] / sigma[2],
        excessKurtosis,
        axisResolved,
        nonGaussianity,
        degenerate: ratio < DEGENERATE_RATIO,
        zeroSpread: false,
        bw,
        bwScale,
        factor,
        bandwidth: bandwidths,
        grid,
        extent,
        cubicBox,
        halfWidths,
        boxHalfWidths,
        axisCoords: axisCoords.map((coords) => Array.from(coords)),
        cellVolume,
        density,
        vmin,
        vmax,
        mass,
        massLevels,
        densityLevels,
        browserPcaKde: true
    };

    if (projections) {
        result.projections = {
            pc12: projection(kernels, bandwidths, axisCoords, 0, 1, count),
            pc13: projection(kernels, bandwidths, axisCoords, 0, 2, count),
            pc23: projection(kernels, bandwidths, axisCoords, 1, 2, count)
        };
    }
    return result;
};

// --- Site extraction ----------------------------------------------------------

const readCellVectors = (text) => {
    const lines = text.split(/\r?\n/);
    let latticeVectors = null;
    let supercell = null;
    lines.forEach((line, index) => {
        const parts = line.trim().split(/\s+/).filter(Boolean);
        if (!parts.length) return;
        if (parts[0] === 'Supercell') supercell = parts.slice(-3).map(Number);
        if (parts[0] === 'Lattice') {
            latticeVectors = [lines[index + 1], lines[index + 2], lines[index + 3]]
                .map((row) => row.trim().split(/\s+/).map(Number));
        }
    });
    if (!latticeVectors || !supercell) throw new Error('Missing lattice or supercell metadata');
    return { latticeVectors, supercell };
};

/**
 * Parse an `.rmc6f` file into per-site Cartesian displacement clouds. Each
 * atom's offset from its own box copy is `coords - cellIndices / supercell`,
 * folded over the supercell boundary; subtracting the site mean leaves the
 * displacement about the average structure, then the supercell lattice maps
 * fractional offsets to Cartesian Angstrom. Mirrors `load_site_displacements`.
 */
// Default fold-and-cluster distance (Å) for reconstructing sites from an old file
// that carries no reference-site or cell columns. Chosen below typical bond lengths
// but well above thermal spread, so genuine sites separate; exposed as a UI knob.
export const DEFAULT_CLUSTER_THRESHOLD = 1.5;

// Element symbols as Python's str.capitalize() normalises them in iter_rmc6f_atoms
// ('SE' -> 'Se'), so both engines label and pool the same species.
const capitalizeElement = (symbol) => {
    const text = String(symbol);
    return text.charAt(0).toUpperCase() + text.slice(1).toLowerCase();
};

// Code-point string order (Python's default), never locale-dependent.
const byCodePoint = (a, b) => (a < b ? -1 : a > b ? 1 : 0);

// A site's composition and label: species counts in name order, and the majority
// species (ties to the alphabetically first) -- the rule of `_site_compositions`.
const siteComposition = (atomElements) => {
    const tally = new Map();
    atomElements.forEach((element) => tally.set(element, (tally.get(element) || 0) + 1));
    const names = [...tally.keys()].sort(byCodePoint);
    const elementCounts = Object.fromEntries(names.map((name) => [name, tally.get(name)]));
    const element = names.reduce((best, name) => (best === null || tally.get(name) > tally.get(best) ? name : best), null) ?? '';
    return { element, elementCounts, mixed: names.length > 1 };
};

// Turn one site's fractional offsets (each atom's position within its own cell, in
// supercell-fractional units) into the site record: mean-centered Cartesian
// displacements (Å) and the site's within-unit-cell fractional position. Offsets
// must already be unwrapped (no boundary split) so the plain mean is the true mean.
// `atomElements` is each row's own species (a reference number can carry several).
const buildSite = (referenceNumber, atomElements, offsets, latticeVectors, supercell, copiesPerCell = null) => {
    const n = offsets.length;
    const mean = [0, 0, 0];
    offsets.forEach((offset) => { mean[0] += offset[0]; mean[1] += offset[1]; mean[2] += offset[2]; });
    mean[0] /= n; mean[1] /= n; mean[2] /= n;
    // Centered fractional offset mapped to Cartesian Angstrom via the box.
    const displacements = offsets.map((offset) => {
        const df = [offset[0] - mean[0], offset[1] - mean[1], offset[2] - mean[2]];
        return [
            df[0] * latticeVectors[0][0] + df[1] * latticeVectors[1][0] + df[2] * latticeVectors[2][0],
            df[0] * latticeVectors[0][1] + df[1] * latticeVectors[1][1] + df[2] * latticeVectors[2][1],
            df[0] * latticeVectors[0][2] + df[1] * latticeVectors[1][2] + df[2] * latticeVectors[2][2]
        ];
    });
    const siteFractional = mean.map((value, i) => {
        const frac = (value * supercell[i]) % 1;
        return (frac + 1) % 1;
    });
    const { element, elementCounts, mixed } = siteComposition(atomElements);
    return {
        referenceNumber, element, elementCounts, mixed, atomElements,
        count: n, displacements, siteFractional, copiesPerCell
    };
};

// numpy's np.round: round half to even, so a tie folds the same way in both engines.
const roundHalfEven = (x) => {
    const r = Math.round(x);
    return Math.abs(x % 1) === 0.5 && r % 2 !== 0 ? r - 1 : r;
};

// Circular mean of period-1 values: arg(sum exp(2 pi i v)) / 2 pi. Mirrors
// `_circular_site_centres`; used only as the reference an unwrap is taken about.
const circularMean = (values) => {
    let sin = 0;
    let cos = 0;
    values.forEach((value) => { sin += Math.sin(2 * Math.PI * value); cos += Math.cos(2 * Math.PI * value); });
    return Math.atan2(sin, cos) / (2 * Math.PI);
};

// Current path: RMCProfile tags every atom with its reference site and box copy, so
// grouping by reference number and subtracting the cell origin gives each cloud.
// The offset `coords - cellIndices / supercell` is the site position plus its
// displacement, known only modulo one supercell period; it is unwrapped about the
// site's own circular mean (o -= round(o - centre)), never about zero -- a fold
// about zero tears a site at x = 1/2 in a one-cell-thick box (or near x = 1 in a
// two-cell one) into halves a box edge apart. Mirrors `load_site_displacements`.
const sitesByReferenceNumber = (atoms, latticeVectors, supercell) => {
    const clouds = new Map();   // referenceNumber -> { elements: [...], offsets: [[dfx,dfy,dfz], ...] }
    atoms.forEach(({ element, referenceNumber, coords, cellIndices }) => {
        const offset = coords.map((value, i) => value - cellIndices[i] / supercell[i]);
        let cloud = clouds.get(referenceNumber);
        if (!cloud) { cloud = { elements: [], offsets: [] }; clouds.set(referenceNumber, cloud); }
        cloud.elements.push(element);
        cloud.offsets.push(offset);
    });
    const referenceNumbers = [...clouds.keys()].sort((a, b) => a - b);
    const sites = referenceNumbers.map((referenceNumber) => {
        const { elements, offsets } = clouds.get(referenceNumber);
        const centre = [0, 1, 2].map((i) => circularMean(offsets.map((offset) => offset[i])));
        const unwrapped = offsets.map((offset) => offset.map((value, i) => value - roundHalfEven(value - centre[i])));
        return buildSite(referenceNumber, elements, unwrapped, latticeVectors, supercell);
    });
    return { referenceNumbers, sites };
};

// Periodic single-linkage clustering of unit-cell fractional points by minimum-image
// Cartesian distance. A uniform grid with bins at least `thresholdA` wide bounds each
// point's neighbour search to its own and adjacent bins (wrapped). All copies of a
// site fold into the same few bins, so the candidate pairs still grow as copies^2;
// what keeps that affordable is testing a pair's distance (a 27-image loop) only
// when the two points are not yet in the same cluster -- a cheap union-find lookup
// that is true for almost every pair once a site has linked up, and that cannot
// change the result (the union would be a no-op). Returns arrays of point indices,
// one per cluster.
const clusterPeriodic = (points, unitVec, thresholdA) => {
    const n = points.length;
    const parent = Array.from({ length: n }, (_, i) => i);
    const find = (x) => { let r = x; while (parent[r] !== r) r = parent[r]; while (parent[x] !== r) { const next = parent[x]; parent[x] = r; x = next; } return r; };
    const union = (a, b) => { const ra = find(a); const rb = find(b); if (ra !== rb) parent[ra] = rb; };

    // Grid over the (fractional) unit cell. Binning by each axis's PERPENDICULAR
    // width (cell volume / opposite-face area), not its vector length, keeps a
    // fractional step of 1/bins spanning >= thresholdA even for an oblique cell, so
    // two points within thresholdA always share a bin or an adjacent one. For an
    // orthogonal cell the perpendicular width equals the edge length.
    const cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    const norm = (u) => Math.hypot(u[0], u[1], u[2]);
    const volume = Math.abs(dot(unitVec[0], cross(unitVec[1], unitVec[2]))) || 1;
    const bins = [0, 1, 2].map((i) => {
        const perpWidth = volume / (norm(cross(unitVec[(i + 1) % 3], unitVec[(i + 2) % 3])) || 1);
        return Math.max(1, Math.floor(perpWidth / Math.max(thresholdA, 1e-9)));
    });
    const wrap = (value) => ((value % 1) + 1) % 1;
    const binOf = (uf) => uf.map((v, i) => Math.min(bins[i] - 1, Math.floor(wrap(v) * bins[i])));
    const keyOf = (b) => `${b[0]},${b[1]},${b[2]}`;

    const buckets = new Map();
    points.forEach((uf, idx) => {
        const key = keyOf(binOf(uf));
        const bucket = buckets.get(key);
        if (bucket) bucket.push(idx); else buckets.set(key, [idx]);
    });

    const thr2 = thresholdA * thresholdA;
    // Is the minimum-image distance under the FULL unit-cell metric below the
    // threshold? Per-axis reduction gives the primary image. For an orthogonal cell
    // (edges mutually perpendicular) |sum d_i a_i|^2 = sum d_i^2 |a_i|^2, so that
    // primary image IS the closest copy and one evaluation decides. An oblique
    // cell's closest copy can be diagonal, which per-axis rounding alone would miss,
    // so the 27 neighbouring images are searched -- stopping at the first within
    // the threshold, since only the yes/no answer is used.
    const orthogonal = [[0, 1], [0, 2], [1, 2]].every(([i, j]) => (
        Math.abs(dot(unitVec[i], unitVec[j])) <= 1e-12 * norm(unitVec[i]) * norm(unitVec[j])
    ));
    const within = (a, b) => {
        let f0 = a[0] - b[0]; f0 -= Math.round(f0);
        let f1 = a[1] - b[1]; f1 -= Math.round(f1);
        let f2 = a[2] - b[2]; f2 -= Math.round(f2);
        const lo = orthogonal ? 0 : -1;
        const hi = orthogonal ? 0 : 1;
        for (let ix = lo; ix <= hi; ix += 1) {
            for (let iy = lo; iy <= hi; iy += 1) {
                for (let iz = lo; iz <= hi; iz += 1) {
                    const d0 = f0 + ix; const d1 = f1 + iy; const d2 = f2 + iz;
                    const x = d0 * unitVec[0][0] + d1 * unitVec[1][0] + d2 * unitVec[2][0];
                    const y = d0 * unitVec[0][1] + d1 * unitVec[1][1] + d2 * unitVec[2][1];
                    const z = d0 * unitVec[0][2] + d1 * unitVec[1][2] + d2 * unitVec[2][2];
                    if (x * x + y * y + z * z < thr2) return true;
                }
            }
        }
        return false;
    };
    points.forEach((uf, idx) => {
        const b = binOf(uf);
        const seen = new Set();
        for (let dx = -1; dx <= 1; dx += 1) {
            for (let dy = -1; dy <= 1; dy += 1) {
                for (let dz = -1; dz <= 1; dz += 1) {
                    const nb = [
                        ((b[0] + dx) % bins[0] + bins[0]) % bins[0],
                        ((b[1] + dy) % bins[1] + bins[1]) % bins[1],
                        ((b[2] + dz) % bins[2] + bins[2]) % bins[2]
                    ];
                    const key = keyOf(nb);
                    if (seen.has(key)) continue;   // few bins -> neighbours alias; scan each once
                    seen.add(key);
                    const bucket = buckets.get(key);
                    if (!bucket) continue;
                    bucket.forEach((j) => {
                        if (j > idx && find(idx) !== find(j) && within(uf, points[j])) union(idx, j);
                    });
                }
            }
        }
    });

    const groups = new Map();
    for (let i = 0; i < n; i += 1) {
        const root = find(i);
        const group = groups.get(root);
        if (group) group.push(i); else groups.set(root, [i]);
    }
    return [...groups.values()];
};

// Fallback path for old files without site/cell columns: fold every atom into a
// single unit cell, cluster per element, and treat each cluster as one site. The
// expected copy count is the supercell product (one image per cell), so a cluster of
// that size is a clean crystallographic site while a multiple flags close/merged or
// orientationally-disordered atoms (e.g. a rotor shell) — surfaced as count/copies.
const sitesByClustering = (atoms, latticeVectors, supercell, thresholdA) => {
    const copiesPerCell = Math.max(1, Math.round(supercell[0] * supercell[1] * supercell[2]));
    // Unit-cell vectors (supercell vectors / counts); distances use the full metric.
    const unitVec = latticeVectors.map((row, i) => row.map((value) => value / Math.max(supercell[i], 1)));
    const folded = atoms.map((atom) => ({
        element: atom.element,
        uf: atom.coords.map((value, i) => { const f = (value * supercell[i]) % 1; return (f + 1) % 1; })
    }));

    const elements = [...new Set(folded.map((atom) => atom.element))].sort();
    const clusters = [];   // { element, members: [uf, ...], centroid: [x,y,z] }
    elements.forEach((element) => {
        const pts = folded.filter((atom) => atom.element === element).map((atom) => atom.uf);
        clusterPeriodic(pts, unitVec, thresholdA).forEach((indices) => {
            const members = indices.map((i) => pts[i]);
            // Unwrap into the frame centred on the cluster's CIRCULAR mean per axis, not
            // an arbitrary member: a wide cluster (e.g. an orientationally-disordered
            // rotor shell spanning more than half a cell edge) would be split by a
            // member-relative min-image and its displacements corrupted, whereas every
            // member lies within half a cell of the circular centre. The circular mean
            // also handles a compact cluster that straddles a cell boundary.
            const TAU = 2 * Math.PI;
            const centre = [0, 1, 2].map((i) => {
                let cos = 0; let sin = 0;
                members.forEach((uf) => { cos += Math.cos(TAU * uf[i]); sin += Math.sin(TAU * uf[i]); });
                const m = Math.atan2(sin, cos) / TAU;
                return m - Math.floor(m);
            });
            const unwrapped = members.map((uf) => uf.map((v, i) => { let d = v - centre[i]; d -= Math.round(d); return centre[i] + d; }));
            const centroid = [0, 1, 2].map((i) => unwrapped.reduce((s, uf) => s + uf[i], 0) / unwrapped.length);
            clusters.push({ element, unwrapped, centroid });
        });
    });

    // Stable ordering: element, then folded centroid (x, y, z) — deterministic across
    // runs so a site keeps its synthetic reference number when the knob is unchanged.
    clusters.sort((a, b) => (
        a.element.localeCompare(b.element)
        || (((a.centroid[0] % 1 + 1) % 1) - ((b.centroid[0] % 1 + 1) % 1))
        || (((a.centroid[1] % 1 + 1) % 1) - ((b.centroid[1] % 1 + 1) % 1))
        || (((a.centroid[2] % 1 + 1) % 1) - ((b.centroid[2] % 1 + 1) % 1))
    ));

    const sites = clusters.map((cluster, index) => {
        const offsets = cluster.unwrapped.map((uf) => uf.map((v, i) => v / supercell[i]));
        const atomElements = offsets.map(() => cluster.element);
        return buildSite(index + 1, atomElements, offsets, latticeVectors, supercell, copiesPerCell);
    });
    return { referenceNumbers: sites.map((site) => site.referenceNumber), sites };
};

/**
 * Parse an `.rmc6f` file into per-site Cartesian displacement clouds. Files that
 * carry the reference-site and cell columns are grouped by reference number; older
 * files that carry only coordinates are reconstructed by folding into one unit cell
 * and clustering (see `sitesByClustering`), and the result is flagged `reconstructed`.
 */
export const siteDisplacementsFromRmc6f = (text, { clusterThreshold = DEFAULT_CLUSTER_THRESHOLD } = {}) => {
    const { latticeVectors, supercell } = readCellVectors(text);
    const atoms = [];
    let inAtoms = false;
    text.split(/\r?\n/).forEach((line) => {
        const parts = line.trim().split(/\s+/).filter(Boolean);
        if (!parts.length) return;
        if (parts[0] === 'Atoms:') { inAtoms = true; return; }
        if (!inAtoms) return;
        const atom = parseAtomLine(parts);
        if (atom) atoms.push({ ...atom, element: capitalizeElement(atom.element) });
    });
    if (atoms.length === 0) throw new Error('No atoms found in structure');

    // Choose the path by majority so a single malformed line can't flip a normal,
    // site-tagged file onto the reconstruction path: if most atoms carry reference
    // and cell columns, group by reference number (dropping any untagged strays);
    // otherwise reconstruct sites by folding every atom into one cell and clustering.
    const tagged = atoms.filter((atom) => atom.referenceNumber !== null && atom.cellIndices !== null);
    const useReferenceNumbers = tagged.length > atoms.length / 2;
    const { referenceNumbers, sites } = useReferenceNumbers
        ? sitesByReferenceNumber(tagged, latticeVectors, supercell)
        : sitesByClustering(atoms, latticeVectors, supercell, clusterThreshold);

    // Per-atom element + supercell-fraction coordinate, in file order. The
    // site extraction above discards absolute positions, but the triplets
    // (bond-angle) request needs every atom's identity and position; keeping
    // them here rides the same parse cache with no second pass over the text.
    const atomList = {
        elements: atoms.map((atom) => atom.element),
        fractional: atoms.map((atom) => atom.coords)
    };

    return { referenceNumbers, sites, latticeVectors, supercell, reconstructed: !useReferenceNumbers, atomList };
};

/** Anisotropic displacement tensor + ellipsoid for every site, in one pass. */
export const siteEllipsoids = (sites, probability = 0.5) => {
    const scale = probabilityScale(probability);
    return sites.map((site) => {
        const { mean, cov } = covariance3(site.displacements);
        const { eigenvalues, axes } = eigenDecomposition(cov);
        const uEq = (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]) / 3;
        const largest = Math.max(eigenvalues[0], 1e-30);
        const zeroSpread = eigenvalues[0] < ZERO_SPREAD_VARIANCE;
        const { excessKurtosis, nonGaussianity, axisResolved } = shapeStatistics(site.displacements, mean, axes, eigenvalues);
        return {
            referenceNumber: site.referenceNumber,
            element: site.element,
            mixed: Boolean(site.mixed),
            elementCounts: site.elementCounts ?? { [site.element]: site.count },
            count: site.count,
            copiesPerCell: site.copiesPerCell ?? null,
            siteFractional: site.siteFractional,
            covariance: cov,
            eigenvalues,
            // A zero-spread site's eigenvectors are round-off: no axes.
            axes: zeroSpread ? null : axes,
            rms: eigenvalues.map(Math.sqrt),
            semiAxes: eigenvalues.map((value) => scale * Math.sqrt(value)),
            probability,
            uIso: uEq,
            bIso: 8 * Math.PI * Math.PI * uEq,
            rmsIso: Math.sqrt(Math.max(uEq, 0)),
            anisotropy: zeroSpread ? null : Math.sqrt(largest / Math.max(eigenvalues[2], 1e-30)),
            excessKurtosis,
            axisResolved,
            nonGaussianity,
            degenerate: zeroSpread || eigenvalues[2] / largest < DEGENERATE_RATIO,
            zeroSpread
        };
    });
};

/**
 * One site's cloud, one element's pooled atoms, or every atom -- the port of
 * `displacement_cloud`. Element pooling selects each atom by its OWN species
 * (`atomElements`), so a mixed-occupancy site contributes only the matching
 * atoms. Returns `{ cloud, site }`, `site` being the record for a reference number.
 */
export const displacementCloud = (parsed, { referenceNumber = null, element = null } = {}) => {
    if (referenceNumber !== null) {
        const site = parsed.sites.find((entry) => entry.referenceNumber === referenceNumber);
        if (!site) throw new Error(`Unknown reference number ${referenceNumber}`);
        return { cloud: site.displacements, site };
    }
    if (element !== null && element !== '' && element !== 'all') {
        const wanted = String(element).toLowerCase();
        const cloud = [];
        parsed.sites.forEach((site) => {
            const own = site.atomElements ?? site.displacements.map(() => site.element);
            site.displacements.forEach((row, i) => { if (own[i].toLowerCase() === wanted) cloud.push(row); });
        });
        if (!cloud.length) {
            throw new Error(`Unknown element ${element}; available: ${siteSpecies(parsed.sites).join(', ')}`);
        }
        return { cloud, site: null };
    }
    const cloud = [];
    parsed.sites.forEach((site) => { site.displacements.forEach((row) => cloud.push(row)); });
    return { cloud, site: null };
};

/** Every species present across the sites (minority species of mixed sites included), sorted. */
export const siteSpecies = (sites) => {
    const names = new Set();
    sites.forEach((site) => Object.keys(site.elementCounts ?? { [site.element]: 1 }).forEach((name) => names.add(name)));
    return [...names].sort(byCodePoint);
};

/** One site's (or one element's pooled) cloud, then `pcaKdeVolume`. */
export const sitePcaKde = (parsed, { referenceNumber = null, element = null, ...options } = {}) => {
    const { cloud, site } = displacementCloud(parsed, { referenceNumber, element });
    const result = pcaKdeVolume(cloud, options);
    if (site) {
        result.referenceNumber = site.referenceNumber;
        result.element = site.element;
        result.elementCounts = site.elementCounts;
        result.mixed = Boolean(site.mixed);
        result.siteFractional = site.siteFractional;
    } else if (element) {
        result.element = String(element);
    }
    return result;
};
