// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { computeDensityGpu, shouldUseGpu } from './gpuKde.js';
import { isInSlab } from './slabSelection.js';

const CUBE_CORNERS = [
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
    [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]
];

const CUBE_EDGES = [
    [0, 1], [0, 2], [1, 3], [2, 3],
    [4, 5], [4, 6], [5, 7], [6, 7],
    [0, 4], [1, 5], [2, 6], [3, 7]
];

const dot = (a, b) => a.reduce((sum, value, index) => sum + value * b[index], 0);
const add = (a, b) => a.map((value, index) => value + b[index]);
const subtract = (a, b) => a.map((value, index) => value - b[index]);
const scale = (vector, factor) => vector.map((value) => value * factor);
const cross = (a, b) => [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
];

const vectorLength = (vector) => Math.sqrt(dot(vector, vector));
const normalize = (vector, fallback = [0, 0, 1]) => {
    const length = vectorLength(vector);
    if (length <= 1e-9) return fallback;
    return vector.map((value) => value / length);
};

const makeFreePlaneBasis = (normal) => {
    const reference = Math.abs(normal[0]) < 0.85 ? [1, 0, 0] : [0, 1, 0];
    const u = normalize(subtract(reference, scale(normal, dot(reference, normal))), [0, 1, 0]);
    const v = normalize(cross(normal, u), [0, 0, 1]);
    return { u, v };
};

const planeSectionVertices = (normal, offset) => {
    const vertices = [];
    CUBE_EDGES.forEach(([start, end]) => {
        const p0 = CUBE_CORNERS[start];
        const p1 = CUBE_CORNERS[end];
        const d0 = dot(p0, normal) - offset;
        const d1 = dot(p1, normal) - offset;
        if (Math.abs(d0) <= 1e-9) vertices.push(p0);
        if (Math.abs(d1) <= 1e-9) vertices.push(p1);
        if (d0 * d1 < 0) {
            const t = d0 / (d0 - d1);
            vertices.push(add(p0, scale(subtract(p1, p0), t)));
        }
    });

    const unique = [];
    vertices.forEach((vertex) => {
        if (!unique.some((existing) => vectorLength(subtract(vertex, existing)) <= 1e-8)) {
            unique.push(vertex);
        }
    });
    if (unique.length < 3) return [];

    const { u, v } = makeFreePlaneBasis(normal);
    const center = scale(unique.reduce((acc, vertex) => add(acc, vertex), [0, 0, 0]), 1 / unique.length);
    return unique.sort((a, b) => {
        const aDelta = subtract(a, center);
        const bDelta = subtract(b, center);
        return Math.atan2(dot(aDelta, v), dot(aDelta, u)) - Math.atan2(dot(bDelta, v), dot(bDelta, u));
    });
};

// Order of preference among an atom's periodic images: fewest, then earliest
// shifts. |offset|^2 * 27 plus the offset's position in the (ox, oy, oz) loop
// below, so the unshifted atom ranks first and every rank is unique; the same
// numbers as _image_rank() in kde.py.
export const imageRank = (ox, oy, oz) => (ox * ox + oy * oy + oz * oz) * 27 + (ox + 1) * 9 + (oy + 1) * 3 + (oz + 1);

// Folding atoms into one unit cell drops their periodic neighbors, so a KDE
// evaluated near a cell face misses the density that should wrap around from
// the opposite face. Tile images from neighbor cells within `margin` of the
// unit cube; sourceIndex maps every kept point back to its source atom so the
// kernel can be normalized by unique atoms rather than by image duplicates,
// and imageRank (see above) picks each atom's representative for the bandwidth.
export const augmentPeriodicImages = (points, margin) => {
    const augmented = [];
    const sourceIndex = [];
    const imageRanks = [];
    points.forEach((point, index) => {
        for (let ox = -1; ox <= 1; ox += 1) {
            for (let oy = -1; oy <= 1; oy += 1) {
                for (let oz = -1; oz <= 1; oz += 1) {
                    const x = point.x + ox;
                    const y = point.y + oy;
                    const z = point.z + oz;
                    if (
                        x >= -margin && x <= 1 + margin
                        && y >= -margin && y <= 1 + margin
                        && z >= -margin && z <= 1 + margin
                    ) {
                        augmented.push({ x, y, z });
                        sourceIndex.push(index);
                        imageRanks.push(imageRank(ox, oy, oz));
                    }
                }
            }
        }
    });
    return { augmented, sourceIndex, imageRanks };
};

// Slab rows (every image in the slab, for the density sum) and one row per
// source atom (its lowest-ranked image in the slab, for the bandwidth), as
// _source_atom_rows() in kde.py; atoms come out in source-index order.
const makeSlab = ({ points, sourceIndex, imageRanks, normal, uVector, vVector, range, zCenter, thickness }) => {
    const depthSpan = range[1] - range[0] || 1;
    const slab = [];
    const atoms = new Map();
    points.forEach((point, index) => {
        const fraction = [point.x, point.y, point.z];
        const normalizedDepth = (dot(fraction, normal) - range[0]) / depthSpan;
        if (isInSlab(normalizedDepth, zCenter, thickness)) {
            const row = [dot(fraction, uVector), dot(fraction, vVector)];
            slab.push(row);
            const source = sourceIndex ? sourceIndex[index] : index;
            const rank = imageRanks ? imageRanks[index] : index;
            const current = atoms.get(source);
            if (!current || rank < current.rank) atoms.set(source, { rank, row });
        }
    });
    const atomRows = [...atoms.entries()].sort(([a], [b]) => a - b).map(([, entry]) => entry.row);
    return { slab, atoms: atomRows, sourceCount: atoms.size };
};

const randomUnit = (seed) => {
    let value = seed >>> 0;
    return () => {
        value += 0x6D2B79F5;
        let mixed = value;
        mixed = Math.imul(mixed ^ (mixed >>> 15), mixed | 1);
        mixed ^= mixed + Math.imul(mixed ^ (mixed >>> 7), mixed | 61);
        return ((mixed ^ (mixed >>> 14)) >>> 0) / 4294967296;
    };
};

const sampleWithoutReplacement = (items, limit, seed = 0) => {
    if (items.length <= limit) return items;
    const random = randomUnit(seed);
    const indices = Array.from({ length: items.length }, (_, index) => index);
    for (let index = 0; index < limit; index += 1) {
        const swapIndex = index + Math.floor(random() * (indices.length - index));
        [indices[index], indices[swapIndex]] = [indices[swapIndex], indices[index]];
    }
    return indices.slice(0, limit).map((index) => items[index]);
};

// Why a slab produced no density. rmc_toolkits/kde.py (KDE_MESSAGES) returns the
// same strings for the same conditions, checked in the same order, so the two
// runtimes either draw the same kernel or decline with the same reason.
export const KDE_MESSAGES = {
    empty: 'No atoms in this slab.',
    bandwidth: 'The bandwidth must be a positive finite number.',
    tooFew: 'Fewer than 5 slab rows: too few atoms for a 2D KDE.',
    fewUnique: 'The slab atoms occupy fewer than 3 distinct in-plane positions, '
        + 'so their covariance (the KDE bandwidth) is undefined.',
    collinear: 'The slab atoms are collinear in the slice plane, '
        + 'so their covariance (the KDE bandwidth) is singular.',
    singular: 'The slab covariance is singular to within round-off, '
        + 'so the KDE bandwidth is undefined.'
};

// Warnings attached to a drawn map, as { code, message }: the same codes,
// strings, thresholds and order as kde.py (KDE_WARNINGS, KERNEL_SUBGRID_RATIO,
// UNRESOLVED_MASS_LIMIT). `subgrid`: a Gaussian sampled at spacing h keeps its
// integral to ~1 % while sigma >= h/2. `unresolved`: the grid-summed linear
// density (nodes times the cell area; about 1 for a resolved map) is below the
// limit, so the map is the tails of kernels that fall between the nodes.
export const KERNEL_SUBGRID_RATIO = 0.5;
export const UNRESOLVED_MASS_LIMIT = 1e-6;
export const KDE_WARNINGS = {
    subgrid: 'The kernel is narrower than half a grid step along its minor axis, so the map is '
        + 'aliased: peak values, contours and the integrated density depend on the grid size. '
        + 'Raise the bandwidth or the grid.',
    unresolved: 'The kernel falls between the grid nodes: the grid holds less than a millionth of '
        + 'the slab\'s density, so the map is round-off and is neither contoured nor drawn. '
        + 'Raise the bandwidth or the grid.'
};

const kernelWarnings = (kernel, gridStep) => (
    kernel.sigmaMinor < KERNEL_SUBGRID_RATIO * gridStep
        ? [{ code: 'subgrid', message: KDE_WARNINGS.subgrid }]
        : []
);

// Decline a covariance with 1 - rho^2 <= this limit (rho = the in-plane
// correlation coefficient) as numerically singular: below it the Cholesky pivot
// is at the level of summation round-off, so its sign would depend on the
// summation order. Same constant and test as kde.py (COVARIANCE_CONDITION_LIMIT).
export const COVARIANCE_CONDITION_LIMIT = 1e-10;

// np.unique(points, axis=0).shape[0] >= minimum, without materializing the set
// beyond what the answer needs. String keys are exact for doubles (shortest
// round-trip form) and fold -0 onto 0, matching numpy's float comparison.
const hasDistinctPoints = (points, minimum) => {
    const seen = new Set();
    for (let index = 0; index < points.length; index += 1) {
        seen.add(`${points[index][0]},${points[index][1]}`);
        if (seen.size >= minimum) return true;
    }
    return false;
};

// np.linalg.matrix_rank(points - mean) >= 2 for an N x 2 point set, with numpy's
// default tolerance S.max() * max(N, 2) * eps. The singular values come from a
// twice-orthogonalized Gram-Schmidt QR of the two centered columns (as accurate
// as Householder QR) and the closed-form SVD of the 2 x 2 triangle, so exactly
// collinear points give sigma_min at round-off level, as numpy's SVD does.
export const hasTwoDimensionalSpread = (points) => {
    const n = points.length;
    if (n < 2) return false;
    let meanU = 0;
    let meanV = 0;
    for (let index = 0; index < n; index += 1) {
        meanU += points[index][0];
        meanV += points[index][1];
    }
    meanU /= n;
    meanV /= n;
    const x = new Float64Array(n);
    const y = new Float64Array(n);
    let r00Squared = 0;
    for (let index = 0; index < n; index += 1) {
        x[index] = points[index][0] - meanU;
        y[index] = points[index][1] - meanV;
        r00Squared += x[index] * x[index];
    }
    const r00 = Math.sqrt(r00Squared);
    if (!(r00 > 0)) return false;
    let r01 = 0;
    for (let pass = 0; pass < 2; pass += 1) {
        let projection = 0;
        for (let index = 0; index < n; index += 1) projection += (x[index] / r00) * y[index];
        for (let index = 0; index < n; index += 1) y[index] -= projection * (x[index] / r00);
        r01 += projection;
    }
    let r11Squared = 0;
    for (let index = 0; index < n; index += 1) r11Squared += y[index] * y[index];
    const r11 = Math.sqrt(r11Squared);
    const frobeniusSquared = r00 * r00 + r01 * r01 + r11 * r11;
    const determinant = Math.abs(r00 * r11);
    const sigmaMax = 0.5 * (Math.sqrt(frobeniusSquared + 2 * determinant)
        + Math.sqrt(Math.max(0, frobeniusSquared - 2 * determinant)));
    const sigmaMin = sigmaMax > 0 ? determinant / sigmaMax : 0;
    return sigmaMin > sigmaMax * Math.max(n, 2) * Number.EPSILON;
};

// Sample covariance with the n-1 divisor (numpy.cov / scipy's gaussian_kde).
export const covariance = (samples) => {
    const n = samples.length;
    const mean = samples.reduce((acc, sample) => [acc[0] + sample[0], acc[1] + sample[1]], [0, 0])
        .map((value) => value / Math.max(n, 1));
    let c00 = 0;
    let c01 = 0;
    let c11 = 0;
    samples.forEach(([x, y]) => {
        const dx = x - mean[0];
        const dy = y - mean[1];
        c00 += dx * dx;
        c01 += dx * dy;
        c11 += dy * dy;
    });
    const denom = Math.max(n - 1, 1);
    return { c00: c00 / denom, c01: c01 / denom, c11: c11 / denom };
};

// Lower Cholesky factor of a 2 x 2 covariance, or null when the covariance is
// not safely positive definite: _well_conditioned() in kde.py, then the
// non-positive (or NaN) pivot where LAPACK's potrf raises LinAlgError.
const cholesky2 = ({ c00, c01, c11 }) => {
    if (!(c00 > 0 && c11 > 0)) return null;
    if (!(1 - (c01 * c01) / (c00 * c11) > COVARIANCE_CONDITION_LIMIT)) return null;
    const l00 = Math.sqrt(c00);
    const l10 = c01 / l00;
    const pivot = c11 - l10 * l10;
    if (!(pivot > 0)) return null;
    return { l00, l10, l11: Math.sqrt(pivot) };
};

// Kernel matrix H and its principal sigmas: the same closed form as
// _kernel_summary() in kde.py (minor eigenvalue from det(H) / lambda_max).
const kernelSummary = (h00, h01, h11, rootDet) => {
    const lambdaMajor = 0.5 * (h00 + h11) + Math.hypot(0.5 * (h00 - h11), h01);
    const lambdaMinor = lambdaMajor > 0 ? (rootDet * rootDet) / lambdaMajor : 0;
    return {
        covariance: [[h00, h01], [h01, h11]],
        sigmaMinor: Math.sqrt(Math.max(lambdaMinor, 0)),
        sigmaMajor: Math.sqrt(Math.max(lambdaMajor, 0))
    };
};

// The Gaussian kernel of scipy.stats.gaussian_kde with a scalar bw_method f:
// H = f^2 C, evaluated exactly (no ridge, no fallback). `cov` is C: the n-1
// covariance of the slab's source atoms (one row per atom, see makeSlab), so the
// periodic images and the subsample leave the kernel alone. The density at a
// node is
//   normalizer * sum_i exp(-0.5 |W (p - p_i)|^2)
// over the summed rows p_i, with W = L^-1 the inverse of the lower Cholesky
// factor L = f chol(C) of H (the whitening SciPy applies through cho_cov) and
// normalizer = weight / (2 pi det(L)). `weight` is the per-row weight: 1 / (rows
// summed), times slab rows / unique source atoms (>= 1) when the rows include
// periodic images, so the map is normalized per source atom. Returns null when
// C is not safely positive definite, where kde.py declines.
export const makeKernel = (cov, factor, weight) => {
    const chol = cholesky2(cov);
    if (!chol) return null;
    const l00 = factor * chol.l00;
    const l10 = factor * chol.l10;
    const l11 = factor * chol.l11;
    const scaleFactor = factor * factor;
    return {
        w00: 1 / l00,
        w10: -l10 / (l00 * l11),
        w11: 1 / l11,
        normalizer: weight / (2 * Math.PI * l00 * l11),
        ...kernelSummary(cov.c00 * scaleFactor, cov.c01 * scaleFactor, cov.c11 * scaleFactor, l00 * l11)
    };
};

const extractContours = ({ density, grid, xMin, xMax, yMin, yMax, vmin, vmax, levels = 8 }) => {
    if (levels <= 0 || !(vmax > vmin)) return [];
    const xStep = (xMax - xMin) / Math.max(grid - 1, 1);
    const yStep = (yMax - yMin) / Math.max(grid - 1, 1);
    const contours = [];

    const interpolate = (level, a, b) => {
        const denom = b.value - a.value;
        const t = Math.abs(denom) <= 1e-12 ? 0.5 : (level - a.value) / denom;
        return [a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t];
    };

    for (let levelIndex = 1; levelIndex <= levels; levelIndex += 1) {
        const level = vmin + (levelIndex / (levels + 1)) * (vmax - vmin);
        const lines = [];
        for (let y = 0; y < grid - 1; y += 1) {
            for (let x = 0; x < grid - 1; x += 1) {
                const corners = [
                    { x: xMin + x * xStep, y: yMin + y * yStep, value: density[y][x] },
                    { x: xMin + (x + 1) * xStep, y: yMin + y * yStep, value: density[y][x + 1] },
                    { x: xMin + (x + 1) * xStep, y: yMin + (y + 1) * yStep, value: density[y + 1][x + 1] },
                    { x: xMin + x * xStep, y: yMin + (y + 1) * yStep, value: density[y + 1][x] }
                ];
                const edgePoints = [];
                [[0, 1], [1, 2], [2, 3], [3, 0]].forEach(([start, end]) => {
                    const a = corners[start];
                    const b = corners[end];
                    if ((a.value < level && b.value >= level) || (b.value < level && a.value >= level)) {
                        edgePoints.push(interpolate(level, a, b));
                    }
                });
                if (edgePoints.length === 2) {
                    lines.push(edgePoints);
                } else if (edgePoints.length === 4) {
                    lines.push([edgePoints[0], edgePoints[1]], [edgePoints[2], edgePoints[3]]);
                }
            }
        }
        if (lines.length) contours.push({ level, lines });
    }
    return contours;
};

// The density map is the worker's hot loop: for every grid cell, sum the
// kernel over every sample (O(grid^2 * samples)). This is the CPU
// implementation, used directly on devices without WebGPU and as the fallback
// whenever the GPU path is unavailable or errors. The WGSL shader in gpuKde.js
// evaluates the same expression from the same kernel fields (w00, w10, w11,
// normalizer); packKdeParams() is the one place that hands them over.
export const computeDensityCpu = ({ samples, kernel, grid, xMin, yMin, xStep, yStep }) => {
    const density = Array.from({ length: grid }, () => new Array(grid).fill(0));
    const { w00, w10, w11, normalizer } = kernel;
    for (let y = 0; y < grid; y += 1) {
        const gy = yMin + y * yStep;
        for (let x = 0; x < grid; x += 1) {
            const gx = xMin + x * xStep;
            let sum = 0;
            for (let index = 0; index < samples.length; index += 1) {
                const dx = gx - samples[index][0];
                const dy = gy - samples[index][1];
                const w0 = w00 * dx;
                const w1 = w10 * dx + w11 * dy;
                const exponent = -0.5 * (w0 * w0 + w1 * w1);
                if (exponent > -60) sum += Math.exp(exponent);
            }
            density[y][x] = sum * normalizer;
        }
    }
    return density;
};

const FIT_LIMIT = 6000;

export const computeKde = async (payload) => {
    const {
        points,
        normal,
        uVector,
        vVector,
        range,
        zCenter,
        thickness,
        bandwidth,
        gridSize,
        logScale
    } = payload;
    const projectedCorners = CUBE_CORNERS.map((corner) => [dot(corner, uVector), dot(corner, vVector)]);
    const xValues = projectedCorners.map(([x]) => x);
    const yValues = projectedCorners.map(([, y]) => y);
    const xMin = Math.min(...xValues);
    const xMax = Math.max(...xValues);
    const yMin = Math.min(...yValues);
    const yMax = Math.max(...yValues);
    // The bandwidth factor is used as given, like kde.py: no substitution and no
    // floor. A non-positive or non-finite value declines the slab below.
    const factor = Number(bandwidth);
    const validBandwidth = typeof bandwidth === 'number' && Number.isFinite(factor) && factor > 0;
    // Margin must cover the kernel reach (sigma scales with bandwidth times the
    // data spread, which is O(1) in fractional units) and the slab depth, so
    // both the in-plane density and the depth selection wrap correctly.
    const margin = Math.min(0.5, Math.max(0.1, 2 * (validBandwidth ? factor : 0), thickness));
    const { augmented, sourceIndex, imageRanks } = augmentPeriodicImages(points, margin);
    const { slab, atoms, sourceCount } = makeSlab({
        points: augmented,
        sourceIndex,
        imageRanks,
        normal,
        uVector,
        vVector,
        range,
        zCenter,
        thickness
    });
    const grid = Math.max(16, Math.min(Number(gridSize) || 120, 260));
    const xStep = (xMax - xMin) / Math.max(grid - 1, 1);
    const yStep = (yMax - yMin) / Math.max(grid - 1, 1);
    let density = null;
    let fitCount = 0;
    let kernelInfo = null;
    let warnings = [];
    let message = KDE_MESSAGES.empty;
    let backend = 'cpu';

    // Decline reasons, in the order kde_slice() checks them.
    if (slab.length === 0) {
        message = KDE_MESSAGES.empty;
    } else if (!validBandwidth) {
        message = KDE_MESSAGES.bandwidth;
    } else if (slab.length < 5) {
        message = KDE_MESSAGES.tooFew;
    } else if (!hasDistinctPoints(atoms, 3)) {
        message = KDE_MESSAGES.fewUnique;
    } else if (!hasTwoDimensionalSpread(atoms)) {
        message = KDE_MESSAGES.collinear;
    } else {
        const samples = sampleWithoutReplacement(slab, FIT_LIMIT, 0);
        const imageFactor = slab.length / Math.max(sourceCount, 1);
        const kernel = makeKernel(covariance(atoms), factor, imageFactor / samples.length);
        if (!kernel) {
            message = KDE_MESSAGES.singular;
        } else {
            message = null;
            fitCount = samples.length;
            kernelInfo = {
                covariance: kernel.covariance,
                sigmaMinor: kernel.sigmaMinor,
                sigmaMajor: kernel.sigmaMajor
            };
            warnings = kernelWarnings(kernel, Math.max(xStep, yStep));
            const args = { samples, kernel, grid, xMin, yMin, xStep, yStep };

            // Run the density map on the GPU when the workload is large enough to
            // amortize the setup cost. Any failure or unavailability falls back to
            // the CPU loop, which evaluates the same kernel in float64.
            let mapped = null;
            if (shouldUseGpu(grid, samples.length)) {
                try {
                    mapped = await computeDensityGpu(args);
                } catch {
                    mapped = null;
                }
            }
            density = mapped ?? computeDensityCpu(args);
            if (mapped) backend = 'gpu';
        }
    }

    // Degenerate slabs never reach the density solver; emit a flat grid so the
    // scaling and contour passes below always have a well-formed array.
    if (!density) {
        density = Array.from({ length: grid }, () => new Array(grid).fill(0));
    }

    // Contour a drawn map whose grid resolves the atoms, judged on the linear
    // values (as kde.py does): a log map whose peak is below 1 is still
    // contoured, a map of kernel tails between the nodes is not.
    let linearSum = 0;
    let vmin = Infinity;
    let vmax = -Infinity;
    for (let y = 0; y < grid; y += 1) {
        for (let x = 0; x < grid; x += 1) {
            linearSum += density[y][x];
            if (logScale) density[y][x] = Math.log10(density[y][x] + 1e-12);
            vmin = Math.min(vmin, density[y][x]);
            vmax = Math.max(vmax, density[y][x]);
        }
    }
    let resolved = false;
    if (kernelInfo) {
        resolved = linearSum * xStep * yStep >= UNRESOLVED_MASS_LIMIT;
        if (!resolved) warnings.push({ code: 'unresolved', message: KDE_WARNINGS.unresolved });
    }

    const depthSpan = range[1] - range[0] || 1;
    const centerDepth = range[0] + zCenter * depthSpan;
    const planeVertices = planeSectionVertices(normal, centerDepth);
    return {
        density,
        extent: [xMin, xMax, yMin, yMax],
        grid,
        z: zCenter,
        dz: thickness,
        bw: validBandwidth ? factor : null,
        log: logScale,
        slabCount: sourceCount,
        fitCount,
        kernel: kernelInfo,
        message,
        warnings,
        vmin: Number.isFinite(vmin) ? vmin : 0,
        vmax: Number.isFinite(vmax) ? vmax : 0,
        contours: resolved ? extractContours({ density, grid, xMin, xMax, yMin, yMax, vmin, vmax }) : [],
        center: zCenter,
        thickness,
        normal,
        uVector,
        vVector,
        planeVertices,
        planePolygon: planeVertices.map((vertex) => [dot(vertex, uVector), dot(vertex, vVector)]),
        browserKde: true,
        backend
    };
};

// Guarded so the module can be imported by tests outside a worker context.
if (typeof self !== 'undefined' && typeof self.postMessage === 'function') {
    self.onmessage = async (event) => {
        try {
            const result = await computeKde(event.data);
            self.postMessage({ id: event.data.id, result });
        } catch (error) {
            self.postMessage({ id: event.data.id, error: error.message || 'Browser KDE computation failed' });
        }
    };
}
