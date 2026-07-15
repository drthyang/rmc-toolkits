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

// Chi-square(3) inverse-CDF sampled on a fine grid, so the ellipsoid scale for a
// given enclosed probability needs no special-function library. k = sqrt(q).
const CHI2_3_TABLE = [
    [0.10, 0.5843744], [0.20, 1.0051740], [0.30, 1.4236522], [0.40, 1.8691684],
    [0.50, 2.3659739], [0.6827, 3.5058779], [0.70, 3.6648530], [0.80, 4.6415889],
    [0.90, 6.2513886], [0.95, 7.8147279], [0.99, 11.3448667], [0.9973, 14.1560750]
];

export const probabilityScale = (probability) => {
    const p = Number(probability);
    if (!(p > 0 && p < 1)) throw new Error('probability must lie strictly between 0 and 1');
    const table = CHI2_3_TABLE;
    if (p <= table[0][0]) return Math.sqrt(table[0][1]);
    if (p >= table[table.length - 1][0]) return Math.sqrt(table[table.length - 1][1]);
    for (let i = 1; i < table.length; i += 1) {
        if (p <= table[i][0]) {
            const [p0, q0] = table[i - 1];
            const [p1, q1] = table[i];
            const t = (p - p0) / (p1 - p0);
            return Math.sqrt(q0 + t * (q1 - q0));
        }
    }
    return Math.sqrt(table[table.length - 1][1]);
};

// Symmetric 3x3 eigendecomposition by cyclic Jacobi rotation. Robust for the
// near-degenerate clouds (flat or linear disorder) that trip analytic formulas,
// and three iterations of a 3x3 sweep are negligible next to the KDE itself.
const jacobiEigenSymmetric = (matrix) => {
    const a = matrix.map((row) => row.slice());
    const v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for (let sweep = 0; sweep < 50; sweep += 1) {
        const off = Math.abs(a[0][1]) + Math.abs(a[0][2]) + Math.abs(a[1][2]);
        if (off < 1e-18) break;
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

// Per-axis excess kurtosis of a centered cloud in its PCA frame (0 = Gaussian).
// Positive means a peaked, fat-tailed distribution whose covariance ellipsoid is
// wider than the KDE isosurface -- the anharmonicity signal the view reveals.
const excessKurtosisPca = (points, mean, axes) => {
    const n = points.length;
    const m2 = [0, 0, 0];
    const m4 = [0, 0, 0];
    points.forEach((point) => {
        const dx = point[0] - mean[0];
        const dy = point[1] - mean[1];
        const dz = point[2] - mean[2];
        for (let a = 0; a < 3; a += 1) {
            const p = dx * axes[a][0] + dy * axes[a][1] + dz * axes[a][2];
            const p2 = p * p;
            m2[a] += p2;
            m4[a] += p2 * p2;
        }
    });
    return [0, 1, 2].map((a) => (m4[a] / n) / Math.max((m2[a] / n) ** 2, 1e-30) - 3);
};

// --- Bandwidth and sampling ---------------------------------------------------

const bandwidthFactor = (method, count, dimensions) => {
    if (typeof method === 'number') {
        if (!(method > 0)) throw new Error('numeric bandwidth must be positive');
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
    const grid = Math.max(8, Math.min(Math.round(gridOption), 128));
    if (!(bwScale > 0)) throw new Error('bwScale must be positive');
    if (!(extent > 0)) throw new Error('extent must be positive');

    const total = points.length;
    const fit = subsample(points, maxFitPoints, rngSeed);
    const count = fit.length;

    const { mean, cov } = covariance3(fit);
    const decomposition = eigenDecomposition(cov);
    const axes = decomposition.axes;
    let eigenvalues = decomposition.eigenvalues;

    const largest = eigenvalues[0];
    if (!(largest > 0)) throw new Error('displacement cloud has zero spread');
    const ratio = eigenvalues[2] / largest;
    eigenvalues = eigenvalues.map((value) => Math.max(value, largest * EIGENVALUE_FLOOR_RATIO));

    const factor = bandwidthFactor(bw, count, 3) * bwScale;
    const sigma = eigenvalues.map(Math.sqrt);
    const bandwidths = sigma.map((value) => factor * value);

    const broaden = Math.sqrt(1 + factor * factor);
    let halfWidths = sigma.map((value) => extent * value * broaden);
    if (cubicBox) {
        const maxHalf = Math.max(...halfWidths);
        halfWidths = [maxHalf, maxHalf, maxHalf];
    }

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

    const excessKurtosis = excessKurtosisPca(fit, mean, axes);
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
        nonGaussianity: (excessKurtosis[0] + excessKurtosis[1] + excessKurtosis[2]) / 3,
        degenerate: ratio < DEGENERATE_RATIO,
        bw,
        bwScale,
        factor,
        bandwidth: bandwidths,
        grid,
        extent,
        cubicBox,
        halfWidths,
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
export const siteDisplacementsFromRmc6f = (text) => {
    const { latticeVectors, supercell } = readCellVectors(text);
    const clouds = new Map();   // referenceNumber -> { element, offsets: [[dfx,dfy,dfz], ...] }
    let inAtoms = false;
    text.split(/\r?\n/).forEach((line) => {
        const parts = line.trim().split(/\s+/).filter(Boolean);
        if (!parts.length) return;
        if (parts[0] === 'Atoms:') { inAtoms = true; return; }
        if (!inAtoms) return;
        const atom = parseAtomLine(parts);
        if (!atom) return;
        const { element, referenceNumber, coords, cellIndices } = atom;
        const offset = coords.map((value, i) => {
            let delta = value - cellIndices[i] / supercell[i];
            delta -= Math.round(delta);
            return delta;
        });
        let cloud = clouds.get(referenceNumber);
        if (!cloud) { cloud = { element, offsets: [] }; clouds.set(referenceNumber, cloud); }
        cloud.offsets.push(offset);
    });
    if (clouds.size === 0) throw new Error('No atoms found in structure');

    const referenceNumbers = [...clouds.keys()].sort((a, b) => a - b);
    const sites = referenceNumbers.map((referenceNumber) => {
        const { element, offsets } = clouds.get(referenceNumber);
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
        return { referenceNumber, element, count: n, displacements, siteFractional };
    });

    return { referenceNumbers, sites, latticeVectors, supercell };
};

/** Anisotropic displacement tensor + ellipsoid for every site, in one pass. */
export const siteEllipsoids = (sites, probability = 0.5) => {
    const scale = probabilityScale(probability);
    return sites.map((site) => {
        const { mean, cov } = covariance3(site.displacements);
        const { eigenvalues, axes } = eigenDecomposition(cov);
        const uEq = (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]) / 3;
        const largest = Math.max(eigenvalues[0], 1e-30);
        const excessKurtosis = excessKurtosisPca(site.displacements, mean, axes);
        return {
            referenceNumber: site.referenceNumber,
            element: site.element,
            count: site.count,
            siteFractional: site.siteFractional,
            covariance: cov,
            eigenvalues,
            axes,
            rms: eigenvalues.map(Math.sqrt),
            semiAxes: eigenvalues.map((value) => scale * Math.sqrt(value)),
            probability,
            uIso: uEq,
            bIso: 8 * Math.PI * Math.PI * uEq,
            rmsIso: Math.sqrt(Math.max(uEq, 0)),
            anisotropy: Math.sqrt(largest / Math.max(eigenvalues[2], 1e-30)),
            excessKurtosis,
            nonGaussianity: (excessKurtosis[0] + excessKurtosis[1] + excessKurtosis[2]) / 3,
            degenerate: eigenvalues[2] / largest < DEGENERATE_RATIO
        };
    });
};

/** One site's (or one element's pooled) cloud, then `pcaKdeVolume`. */
export const sitePcaKde = (parsed, { referenceNumber = null, element = null, ...options } = {}) => {
    let cloud = [];
    let tagged = null;
    if (referenceNumber !== null) {
        tagged = parsed.sites.find((site) => site.referenceNumber === referenceNumber);
        if (!tagged) throw new Error(`Unknown reference number ${referenceNumber}`);
        cloud = tagged.displacements;
    } else if (element !== null && element !== '' && element !== 'all') {
        const matches = parsed.sites.filter(
            (site) => site.element.toLowerCase() === String(element).toLowerCase()
        );
        if (!matches.length) throw new Error(`Unknown element ${element}`);
        matches.forEach((site) => { cloud = cloud.concat(site.displacements); });
    } else {
        parsed.sites.forEach((site) => { cloud = cloud.concat(site.displacements); });
    }

    const result = pcaKdeVolume(cloud, options);
    if (tagged) {
        result.referenceNumber = tagged.referenceNumber;
        result.element = tagged.element;
        result.siteFractional = tagged.siteFractional;
    } else if (element) {
        result.element = String(element);
    }
    return result;
};
