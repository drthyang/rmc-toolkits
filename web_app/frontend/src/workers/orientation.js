// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Orientation distribution of RMC displacement vectors, hex-binned on a sphere
// (browser / static mode).
//
// This is the JavaScript port of `rmc_toolkits/orientation.py`; the Python
// module is the source of truth and documents the physics. The short version:
// keep only the *direction* of each displacement `dr = r_i - r_avg`, bin those
// unit vectors into the cells of a Goldberg polyhedron (the dual of a geodesic
// icosahedron: 10*nu^2 + 2 cells, hexagons plus the 12 pentagons any hexagonal
// tiling of a sphere must carry), and divide each cell's count by its exact
// solid angle. `enhancement = 4*pi*density` is 1 everywhere for an isotropic
// site, so the map reads directly as "this direction is Nx more likely than
// chance"; `zScore` is the Poisson significance that keeps an over-binned map
// from being read as structure. The map is never antipodally folded -- a +u/-u
// imbalance (static off-centring, odd anharmonicity) is precisely the signal
// the ellipsoid cannot show, and `antipodalAsymmetry` quantifies it.
//
// Keep this file in sync with orientation.py -- same constants, same outputs.

import { eigenDecomposition } from './pcaKde.js';

export const MIN_FREQUENCY = 1;
export const MAX_FREQUENCY = 64;
export const NEGLIGIBLE_AMPLITUDE = 1e-9;
export const SMOOTHING_ALPHA = 0.5;
export const DEFAULT_TARGET_PER_CELL = 12;

const WEIGHTS = ['count', 'amplitude', 'amplitude2'];
const FRAMES = ['cartesian', 'pca'];

const PHI = (1 + Math.sqrt(5)) / 2;

// --- small vector helpers -----------------------------------------------------

const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
const cross = (u, v) => [
    u[1] * v[2] - u[2] * v[1],
    u[2] * v[0] - u[0] * v[2],
    u[0] * v[1] - u[1] * v[0]
];
const norm = (u) => Math.hypot(u[0], u[1], u[2]);
const normalize = (u) => {
    const length = Math.max(norm(u), 1e-300);
    return [u[0] / length, u[1] / length, u[2] / length];
};

// Inverse of a 3x3 matrix given as rows; used for the gnomonic face inversion.
const inverse3 = (m) => {
    const [[a, b, c], [d, e, f], [g, h, i]] = m;
    const A = e * i - f * h;
    const B = c * h - b * i;
    const C = b * f - c * e;
    const det = a * A + d * B + g * C;
    return [
        [A / det, B / det, C / det],
        [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
        [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det]
    ];
};

// --- icosahedron + geodesic subdivision ---------------------------------------

// 12 unit vertices and 20 outward, counter-clockwise faces, derived (not
// hard-coded) exactly as in the Python engine so cell indices match it.
const icosahedron = () => {
    const raw = [];
    for (const s1 of [1, -1]) {
        for (const s2 of [1, -1]) {
            raw.push([0, s1, s2 * PHI]);
            raw.push([s1, s2 * PHI, 0]);
            raw.push([s2 * PHI, 0, s1]);
        }
    }
    let edge = Infinity;
    for (let a = 0; a < 12; a += 1) {
        for (let b = a + 1; b < 12; b += 1) {
            const d = norm([raw[a][0] - raw[b][0], raw[a][1] - raw[b][1], raw[a][2] - raw[b][2]]);
            if (d > 1e-9 && d < edge) edge = d;
        }
    }
    const adjacent = (a, b) => {
        const d = norm([raw[a][0] - raw[b][0], raw[a][1] - raw[b][1], raw[a][2] - raw[b][2]]);
        return Math.abs(d - edge) < 1e-9 * Math.max(edge, 1);
    };
    const faces = [];
    for (let a = 0; a < 12; a += 1) {
        for (let b = a + 1; b < 12; b += 1) {
            if (!adjacent(a, b)) continue;
            for (let c = b + 1; c < 12; c += 1) {
                if (adjacent(a, c) && adjacent(b, c)) faces.push([a, b, c]);
            }
        }
    }
    if (faces.length !== 20) throw new Error(`icosahedron construction produced ${faces.length} faces`);

    const vertices = raw.map(normalize);
    const oriented = faces.map(([a, b, c]) => {
        const ab = [vertices[b][0] - vertices[a][0], vertices[b][1] - vertices[a][1], vertices[b][2] - vertices[a][2]];
        const ac = [vertices[c][0] - vertices[a][0], vertices[c][1] - vertices[a][1], vertices[c][2] - vertices[a][2]];
        const normal = cross(ab, ac);
        const centroid = [
            vertices[a][0] + vertices[b][0] + vertices[c][0],
            vertices[a][1] + vertices[b][1] + vertices[c][1],
            vertices[a][2] + vertices[b][2] + vertices[c][2]
        ];
        return dot(normal, centroid) < 0 ? [a, c, b] : [a, b, c];
    });
    return { vertices, faces: oriented };
};

// Orientation-independent key for a point on the icosahedron edge u--v.
const edgeKey = (u, v, position, nu) => (
    u < v ? `e:${u}:${v}:${position}` : `e:${v}:${u}:${nu - position}`
);

// Geodesic vertices, the per-face lattice map, and the CCW triangles. Shared
// points are merged by exact combinatorial key (corner / edge id), never by
// rounding coordinates. Mirrors `_subdivide`.
const subdivide = (vertices, faces, nu) => {
    const lattice = faces.map(() => {
        const rows = [];
        for (let i = 0; i <= nu; i += 1) rows.push(new Int32Array(nu + 1).fill(-1));
        return rows;
    });
    const indexOf = new Map();
    const points = [];

    faces.forEach(([ia, ib, ic], faceId) => {
        const A = vertices[ia];
        const B = vertices[ib];
        const C = vertices[ic];
        for (let i = 0; i <= nu; i += 1) {
            for (let j = 0; j <= nu - i; j += 1) {
                const k = nu - i - j;
                let key;
                if (i === nu) key = `v:${ib}`;
                else if (j === nu) key = `v:${ic}`;
                else if (k === nu) key = `v:${ia}`;
                else if (k === 0) key = edgeKey(ib, ic, j, nu);
                else if (i === 0) key = edgeKey(ia, ic, j, nu);
                else if (j === 0) key = edgeKey(ia, ib, i, nu);
                else key = `f:${faceId}:${i}:${j}`;
                let index = indexOf.get(key);
                if (index === undefined) {
                    const point = normalize([
                        (k * A[0] + i * B[0] + j * C[0]) / nu,
                        (k * A[1] + i * B[1] + j * C[1]) / nu,
                        (k * A[2] + i * B[2] + j * C[2]) / nu
                    ]);
                    index = points.length;
                    indexOf.set(key, index);
                    points.push(point);
                }
                lattice[faceId][i][j] = index;
            }
        }
    });

    const expected = 10 * nu * nu + 2;
    if (points.length !== expected) {
        throw new Error(`geodesic subdivision produced ${points.length} vertices, expected ${expected}`);
    }

    const triangles = [];
    for (let f = 0; f < faces.length; f += 1) {
        for (let i = 0; i < nu; i += 1) {
            for (let j = 0; j < nu - i; j += 1) {
                triangles.push([lattice[f][i][j], lattice[f][i + 1][j], lattice[f][i][j + 1]]);
            }
        }
    }
    for (let f = 0; f < faces.length; f += 1) {
        for (let i = 0; i < nu - 1; i += 1) {
            for (let j = 0; j < nu - 1 - i; j += 1) {
                triangles.push([lattice[f][i + 1][j], lattice[f][i + 1][j + 1], lattice[f][i][j + 1]]);
            }
        }
    }
    return { centers: points, lattice, triangles };
};

// Right-handed tangent frame (e1, e2, n) at a unit normal, so atan2 of the
// (e2, e1) components orders points counter-clockwise seen from outside.
const tangentBasis = (normal) => {
    const reference = [0, 0, 0];
    let least = 0;
    if (Math.abs(normal[1]) < Math.abs(normal[least])) least = 1;
    if (Math.abs(normal[2]) < Math.abs(normal[least])) least = 2;
    reference[least] = 1;
    const projection = dot(reference, normal);
    const e1 = normalize([
        reference[0] - normal[0] * projection,
        reference[1] - normal[1] * projection,
        reference[2] - normal[2] * projection
    ]);
    return { e1, e2: cross(normal, e1) };
};

// Sort each cell's incident items counter-clockwise about the cell centre.
// `entries` is an array of {owner, point, payload}; returns per-owner payload
// lists in CCW order. Mirrors `_angular_order`.
const angularOrder = (centers, entries) => {
    const byOwner = centers.map(() => []);
    entries.forEach(({ owner, point, payload }) => {
        const center = centers[owner];
        const { e1, e2 } = tangentBasis(center);
        const radial = dot(point, center);
        const local = [
            point[0] - center[0] * radial,
            point[1] - center[1] * radial,
            point[2] - center[2] * radial
        ];
        byOwner[owner].push({ angle: Math.atan2(dot(local, e2), dot(local, e1)), payload });
    });
    return byOwner.map((list) => list.sort((a, b) => a.angle - b.angle).map((item) => item.payload));
};

// Solid angle of a spherical polygon (unit vertices) by a signed triangle fan,
// in the van Oosterom-Strackee form. A repeated vertex contributes exactly 0.
const sphericalPolygonArea = (polygon) => {
    const a = polygon[0];
    let total = 0;
    for (let i = 1; i < polygon.length - 1; i += 1) {
        const b = polygon[i];
        const c = polygon[i + 1];
        const numerator = dot(a, cross(b, c));
        const denominator = 1 + dot(a, b) + dot(b, c) + dot(c, a);
        total += 2 * Math.atan2(numerator, denominator);
    }
    return total;
};

// Nearest cell centre for one unit direction: gnomonic seed (exact ray/cone
// face test + largest-remainder lattice rounding) then a greedy walk on the
// adjacency graph. Mirrors `_assign`.
const assignOne = (tiling, direction) => {
    const { faceInverse, lattice, frequency: nu, centers, neighbors } = tiling;
    let bestFace = 0;
    let bestMin = -Infinity;
    let bestLam = null;
    for (let f = 0; f < faceInverse.length; f += 1) {
        const m = faceInverse[f];
        const l0 = m[0][0] * direction[0] + m[0][1] * direction[1] + m[0][2] * direction[2];
        const l1 = m[1][0] * direction[0] + m[1][1] * direction[1] + m[1][2] * direction[2];
        const l2 = m[2][0] * direction[0] + m[2][1] * direction[1] + m[2][2] * direction[2];
        const least = Math.min(l0, l1, l2);
        if (least > bestMin) {
            bestMin = least;
            bestFace = f;
            bestLam = [l0, l1, l2];
        }
    }
    const sum = Math.max(bestLam[0] + bestLam[1] + bestLam[2], 1e-300);
    const raw = bestLam.map((value) => nu * Math.min(Math.max(value / sum, 0), 1));
    const floor = raw.map(Math.floor);
    let deficit = nu - floor[0] - floor[1] - floor[2];
    const order = [0, 1, 2].sort((a, b) => (raw[b] - floor[b]) - (raw[a] - floor[a]));
    const integral = floor.slice();
    for (let slot = 0; slot < 3 && deficit > 0; slot += 1, deficit -= 1) {
        integral[order[slot]] += 1;
    }

    let current = lattice[bestFace][integral[1]][integral[2]];
    if (current < 0) current = 0;

    for (let round = 0; round < 8; round += 1) {
        let best = dot(centers[current], direction);
        let next = current;
        const row = neighbors[current];
        for (let k = 0; k < row.length; k += 1) {
            const candidate = row[k];
            if (candidate < 0) continue;
            const value = dot(centers[candidate], direction);
            if (value > best) {
                best = value;
                next = candidate;
            }
        }
        if (next === current) break;
        current = next;
    }
    return current;
};

// --- tiling construction (cached) ---------------------------------------------

const tilingCache = new Map();

/**
 * Hex-and-12-pentagon (Goldberg) tiling of the unit sphere at geodesic
 * frequency `nu`: 10*nu^2 + 2 cells. Cached per frequency. Mirrors
 * `goldberg_tiling` -- same construction order, so cell indices agree with the
 * Python engine.
 */
export const goldbergTiling = (frequency = 8) => {
    const nu = Math.trunc(frequency);
    if (!(nu >= MIN_FREQUENCY && nu <= MAX_FREQUENCY)) {
        throw new Error(`frequency must lie in [${MIN_FREQUENCY}, ${MAX_FREQUENCY}]`);
    }
    const cached = tilingCache.get(nu);
    if (cached) return cached;

    const { vertices, faces } = icosahedron();
    const { centers, lattice, triangles } = subdivide(vertices, faces, nu);
    const cellCount = centers.length;

    // Dual: each triangle's circumcentre (the normalized CCW plane normal) is a
    // polygon vertex of its three corner cells.
    const circumcenters = triangles.map(([p0, p1, p2]) => {
        const u = [centers[p1][0] - centers[p0][0], centers[p1][1] - centers[p0][1], centers[p1][2] - centers[p0][2]];
        const v = [centers[p2][0] - centers[p0][0], centers[p2][1] - centers[p0][1], centers[p2][2] - centers[p0][2]];
        return normalize(cross(u, v));
    });

    const polygonEntries = [];
    triangles.forEach((triangle, t) => {
        triangle.forEach((owner) => {
            polygonEntries.push({ owner, point: circumcenters[t], payload: circumcenters[t] });
        });
    });
    const polygonsRagged = angularOrder(centers, polygonEntries);
    const sizes = polygonsRagged.map((polygon) => polygon.length);
    const pentagonCount = sizes.filter((size) => size === 5).length;
    if (pentagonCount !== 12 || sizes.some((size) => size !== 5 && size !== 6)) {
        throw new Error('dual tiling is not a Goldberg polyhedron (bad cell degrees)');
    }
    // Pad pentagons by repeating the last vertex (a zero-area fan step).
    const polygons = polygonsRagged.map((polygon) => (
        polygon.length === 5 ? [...polygon, polygon[4]] : polygon
    ));

    // Directed edges: each ordered pair occurs exactly once over CCW triangles.
    const edgeEntries = [];
    triangles.forEach(([a, b, c]) => {
        edgeEntries.push({ owner: a, point: centers[b], payload: b });
        edgeEntries.push({ owner: b, point: centers[c], payload: c });
        edgeEntries.push({ owner: c, point: centers[a], payload: a });
    });
    const neighborsRagged = angularOrder(centers, edgeEntries);
    if (neighborsRagged.some((row, cell) => row.length !== sizes[cell])) {
        throw new Error('cell adjacency disagrees with the dual polygon degrees');
    }
    const neighbors = neighborsRagged.map((row) => (
        row.length === 5 ? [...row, -1] : row
    ));

    const areas = polygons.map(sphericalPolygonArea);
    const total = areas.reduce((sum, area) => sum + area, 0);
    if (Math.abs(total - 4 * Math.PI) > 1e-9 * 4 * Math.PI) {
        throw new Error(`cell areas sum to ${total}, not 4*pi`);
    }

    // Columns of the face matrix are the corner vertices: u = M . lambda.
    const faceInverse = faces.map(([a, b, c]) => inverse3([
        [vertices[a][0], vertices[b][0], vertices[c][0]],
        [vertices[a][1], vertices[b][1], vertices[c][1]],
        [vertices[a][2], vertices[b][2], vertices[c][2]]
    ]));

    const tiling = {
        frequency: nu,
        centers,
        polygons,
        sizes,
        areas,
        neighbors,
        isPentagon: sizes.map((size) => size === 5),
        faceInverse,
        lattice,
        cellCount,
        antipode: null
    };

    // The tiling is centrosymmetric, so -c is exactly another cell centre. The
    // antipode map powers the antipodal-asymmetry readout (u vs -u imbalance).
    tiling.antipode = centers.map((center) => assignOne(tiling, [-center[0], -center[1], -center[2]]));
    const worst = tiling.antipode.reduce((acc, other, cell) => Math.max(
        acc,
        Math.abs(centers[other][0] + centers[cell][0]),
        Math.abs(centers[other][1] + centers[cell][1]),
        Math.abs(centers[other][2] + centers[cell][2])
    ), 0);
    if (worst > 1e-9) throw new Error('tiling is not centrosymmetric; the antipode map would be wrong');

    tilingCache.set(nu, tiling);
    return tiling;
};

/** Cell index for each direction ([x, y, z], need not be normalized). */
export const assignCells = (tiling, directions) => (
    directions.map((direction) => assignOne(tiling, normalize(direction)))
);

/**
 * Largest frequency whose cells still average `targetPerCell` points -- the
 * over-binning guard. Mirrors `recommended_frequency`.
 */
export const recommendedFrequency = (nPoints, { targetPerCell = DEFAULT_TARGET_PER_CELL, maxFrequency = 24 } = {}) => {
    if (!(nPoints > 0)) return MIN_FREQUENCY;
    const cells = Math.max(12, nPoints / Math.max(Math.trunc(targetPerCell), 1));
    const frequency = Math.round(Math.sqrt(Math.max(cells - 2, 10) / 10));
    return Math.min(Math.max(frequency, MIN_FREQUENCY), Math.min(maxFrequency, MAX_FREQUENCY));
};

// --- histogram ----------------------------------------------------------------

// Neighbour diffusion on the cell graph; exactly mass-conserving.
const smooth = (mass, neighbors, passes) => {
    let current = mass;
    for (let pass = 0; pass < passes; pass += 1) {
        const updated = current.map((value) => (1 - SMOOTHING_ALPHA) * value);
        current.forEach((value, cell) => {
            const row = neighbors[cell];
            let degree = 0;
            for (let k = 0; k < row.length; k += 1) if (row[k] >= 0) degree += 1;
            const share = SMOOTHING_ALPHA * value / Math.max(degree, 1);
            for (let k = 0; k < row.length; k += 1) {
                if (row[k] >= 0) updated[row[k]] += share;
            }
        });
        current = updated;
    }
    return current;
};

// numpy's default (linear-interpolation) quantile of an unsorted sample.
const quantile = (values, q) => {
    const sorted = [...values].sort((a, b) => a - b);
    const position = q * (sorted.length - 1);
    const low = Math.floor(position);
    const high = Math.min(low + 1, sorted.length - 1);
    return sorted[low] + (position - low) * (sorted[high] - sorted[low]);
};

const covariance3 = (points) => {
    const n = points.length;
    const mean = [0, 0, 0];
    points.forEach((point) => { mean[0] += point[0]; mean[1] += point[1]; mean[2] += point[2]; });
    mean[0] /= n; mean[1] /= n; mean[2] /= n;
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
    return cov;
};

/**
 * Solid-angle distribution of the directions of `vectors` (array of [x, y, z]
 * Cartesian Angstrom displacements; only directions are used). Mirrors
 * `orientation_histogram` in Python -- same options, same output keys:
 * `density` integrates to 1 over the sphere, `enhancement = 4*pi*density` is 1
 * for an isotropic cloud, `zScore` is the Poisson significance per cell.
 */
export const orientationHistogram = (vectors, options = {}) => {
    const {
        frequency = null,
        weight = 'count',
        minAmplitude = 0,
        minAmplitudeQuantile = 0,
        smoothing = 0,
        frame = 'cartesian',
        geometry = true,
        targetPerCell = DEFAULT_TARGET_PER_CELL
    } = options;

    if (!Array.isArray(vectors) || vectors.some((row) => !row || row.length !== 3)) {
        throw new Error('vectors must be an array of [x, y, z] rows');
    }
    if (!WEIGHTS.includes(weight)) throw new Error(`weight must be one of ${WEIGHTS.join(', ')}`);
    if (!FRAMES.includes(frame)) throw new Error(`frame must be one of ${FRAMES.join(', ')}`);
    if (!(minAmplitudeQuantile >= 0 && minAmplitudeQuantile < 1)) {
        throw new Error('minAmplitudeQuantile must lie in [0, 1)');
    }

    const totalPoints = vectors.length;
    const amplitude = vectors.map((row) => norm(row));

    // PCA rotation fitted on every vector, before any amplitude cut, so the
    // frame is the same one the ellipsoid view reports for this site.
    let pcaAxes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    if (totalPoints >= 4) {
        pcaAxes = eigenDecomposition(covariance3(vectors)).axes;
    }

    let threshold = Number(minAmplitude);
    if (minAmplitudeQuantile > 0 && amplitude.length) {
        threshold = Math.max(threshold, quantile(amplitude, minAmplitudeQuantile));
    }
    const cutoff = Math.max(threshold, NEGLIGIBLE_AMPLITUDE);
    const keptIndices = [];
    for (let i = 0; i < totalPoints; i += 1) {
        if (amplitude[i] > cutoff) keptIndices.push(i);
    }
    const used = keptIndices.length;
    if (used < 1) throw new Error('no displacement vectors survive the amplitude cutoff');

    let directions = keptIndices.map((i) => [
        vectors[i][0] / amplitude[i],
        vectors[i][1] / amplitude[i],
        vectors[i][2] / amplitude[i]
    ]);
    const keptAmplitude = keptIndices.map((i) => amplitude[i]);
    if (frame === 'pca') {
        directions = directions.map((u) => [
            dot(u, pcaAxes[0]), dot(u, pcaAxes[1]), dot(u, pcaAxes[2])
        ]);
    }

    const resolvedFrequency = frequency === null || frequency === undefined
        ? recommendedFrequency(used, { targetPerCell })
        : Math.trunc(frequency);
    const tiling = goldbergTiling(resolvedFrequency);
    const cellCount = tiling.cellCount;

    let weights;
    if (weight === 'amplitude') weights = keptAmplitude;
    else if (weight === 'amplitude2') weights = keptAmplitude.map((value) => value * value);
    else weights = new Array(used).fill(1);

    const counts = new Float64Array(cellCount);
    let mass = new Float64Array(cellCount);
    directions.forEach((direction, index) => {
        const cell = assignOne(tiling, direction);
        counts[cell] += 1;
        mass[cell] += weights[index];
    });

    if (smoothing > 0) {
        mass = Float64Array.from(smooth(Array.from(mass), tiling.neighbors, Math.trunc(smoothing)));
    }

    let totalMass = 0;
    for (let cell = 0; cell < cellCount; cell += 1) totalMass += mass[cell];
    if (!(totalMass > 0)) throw new Error('displacement weights sum to zero');

    const density = new Array(cellCount);
    const enhancement = new Array(cellCount);
    for (let cell = 0; cell < cellCount; cell += 1) {
        density[cell] = mass[cell] / (totalMass * tiling.areas[cell]);
        enhancement[cell] = density[cell] * 4 * Math.PI;
    }

    // Poisson significance from the *raw* counts (smoothing would correlate
    // cells and overstate the evidence).
    const expected = new Array(cellCount);
    const zScore = new Array(cellCount);
    for (let cell = 0; cell < cellCount; cell += 1) {
        expected[cell] = used * tiling.areas[cell] / (4 * Math.PI);
        zScore[cell] = (counts[cell] - expected[cell]) / Math.sqrt(Math.max(expected[cell], 1e-12));
    }

    // Antipodal (inversion) asymmetry: sum_pairs |n(u) - n(-u)| / N -- 0 for an
    // inversion-symmetric cloud, 1 for a fully one-sided one. The null level is
    // what pure Poisson noise alone would produce (|X - Y| of two iid
    // Poisson(m) cells averages ~2*sqrt(m/pi)); report both so a small
    // asymmetry is not over-read.
    let asymmetrySum = 0;
    for (let cell = 0; cell < cellCount; cell += 1) {
        asymmetrySum += Math.abs(counts[cell] - counts[tiling.antipode[cell]]);
    }
    const antipodalAsymmetry = 0.5 * asymmetrySum / used;
    const antipodalAsymmetryNull = Math.sqrt(cellCount / (Math.PI * used));

    // Orientation tensor T = <u u^T> (I/3 for a uniform sphere), from the
    // points, not the binned map -- resolution independent.
    const tensor = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    let weightSum = 0;
    directions.forEach((u, index) => {
        const w = weights[index];
        weightSum += w;
        for (let a = 0; a < 3; a += 1) {
            for (let b = 0; b < 3; b += 1) tensor[a][b] += w * u[a] * u[b];
        }
    });
    for (let a = 0; a < 3; a += 1) {
        for (let b = 0; b < 3; b += 1) tensor[a][b] /= weightSum;
    }
    const tensorDecomposition = eigenDecomposition(tensor);

    let peak = 0;
    for (let cell = 1; cell < cellCount; cell += 1) {
        if (enhancement[cell] > enhancement[peak]) peak = cell;
    }
    let vmin = Infinity;
    let vmax = -Infinity;
    let emptyCells = 0;
    let zSquares = 0;
    for (let cell = 0; cell < cellCount; cell += 1) {
        if (enhancement[cell] < vmin) vmin = enhancement[cell];
        if (enhancement[cell] > vmax) vmax = enhancement[cell];
        if (counts[cell] === 0) emptyCells += 1;
        zSquares += zScore[cell] * zScore[cell];
    }
    let amplitudeSum = 0;
    let amplitudeSquares = 0;
    keptAmplitude.forEach((value) => { amplitudeSum += value; amplitudeSquares += value * value; });

    const result = {
        frequency: tiling.frequency,
        cellCount,
        pentagonCount: 12,
        totalPoints,
        usedPoints: used,
        rejectedPoints: totalPoints - used,
        amplitudeCutoff: Math.max(threshold, 0),
        weight,
        smoothing: Math.trunc(smoothing),
        frame,
        pcaAxes,
        centers: tiling.centers,
        areas: tiling.areas,
        sizes: tiling.sizes,
        antipode: tiling.antipode,
        counts: Array.from(counts),
        mass: Array.from(mass),
        density,
        enhancement,
        expected,
        zScore,
        vmin,
        vmax,
        meanCount: used / cellCount,
        emptyFraction: emptyCells / cellCount,
        meanAmplitude: amplitudeSum / used,
        rmsAmplitude: Math.sqrt(amplitudeSquares / used),
        antipodalAsymmetry,
        antipodalAsymmetryNull,
        orientationTensor: tensor,
        orientationEigenvalues: tensorDecomposition.eigenvalues,
        orientationAxes: tensorDecomposition.axes,
        orientationAnisotropy: 3 * tensorDecomposition.eigenvalues[0] - 1,
        peakCell: peak,
        peakDirection: tiling.centers[peak],
        peakEnhancement: enhancement[peak],
        peakZScore: zScore[peak],
        significance: Math.sqrt(zSquares / cellCount),
        recommendedFrequency: recommendedFrequency(used, { targetPerCell }),
        browserOrientation: true
    };
    if (geometry) {
        result.polygons = tiling.polygons.map((polygon, cell) => polygon.slice(0, tiling.sizes[cell]));
        result.neighbors = tiling.neighbors;
    }
    return result;
};

/**
 * `orientationHistogram` for one site (or one element's pooled sites) of a
 * parsed `.rmc6f` (the object `siteDisplacementsFromRmc6f` returns). Mirrors
 * `site_orientation_histogram`.
 */
export const siteOrientationHistogram = (parsed, { referenceNumber = null, element = null, ...options } = {}) => {
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

    const result = orientationHistogram(cloud, options);
    if (tagged) {
        result.referenceNumber = tagged.referenceNumber;
        result.element = tagged.element;
        result.siteFractional = tagged.siteFractional;
    } else if (element) {
        result.element = String(element);
    }
    return result;
};
