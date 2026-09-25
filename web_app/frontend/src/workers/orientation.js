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
// chance"; `zScore` is each cell's local (uncorrected) z, and the calibrated
// readouts are the `...Significance` fields (one-sided normal deviates of
// exact tail probabilities). The map is never antipodally folded -- a +u/-u
// imbalance about the site mean (skewness: odd anharmonicity, or unequal
// occupation of opposite off-centre wells) is precisely the signal the
// ellipsoid cannot show, and `antipodalAsymmetry` quantifies it. Displacements
// are measured from the configuration's own site mean, so a coherent shift of
// every copy (ordered off-centring) is invisible here by construction.
//
// Keep this file in sync with orientation.py -- same constants, same outputs.

import { eigenDecomposition } from './pcaKde.js';

export const MIN_FREQUENCY = 1;
export const MAX_FREQUENCY = 64;
export const NEGLIGIBLE_AMPLITUDE = 1e-9;
export const SMOOTHING_ALPHA = 0.5;
export const DEFAULT_TARGET_PER_CELL = 12;
// Relative tolerance of the tied-peak rule (lowest index within it of the
// maximum wins). Mirrors PEAK_TIE_RTOL.
export const PEAK_TIE_RTOL = 1e-9;
// Tail probabilities are floored here before conversion to a normal deviate,
// so |z| <= 37.04 and no infinity reaches the payload. Mirrors
// SIGNIFICANCE_TAIL_FLOOR.
export const SIGNIFICANCE_TAIL_FLOOR = 1e-300;
// Antipodal-asymmetry flag: A > null mean + this many null SDs. Mirrors
// ASYMMETRY_FLAG_SIGMA.
export const ASYMMETRY_FLAG_SIGMA = 3;
// Isotropic expectation of 3*lambda1 - 1 is this / sqrt(N_eff) to leading
// order (9 / sqrt(10 pi)). Mirrors ISOTROPIC_ANISOTROPY_SCALE.
export const ISOTROPIC_ANISOTROPY_SCALE = 9 / Math.sqrt(10 * Math.PI);

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

// --- special functions (scipy stand-ins for the significance readouts) -------
//
// Python uses scipy.special.gammainc / gammaincc / ndtri; these ports agree
// with them to <= 1e-12 relative on the small tail (pinned against scipy
// values in orientationFixes.test.js), far inside the 1e-9 golden tolerance.

const LOG_SQRT_2PI = 0.5 * Math.log(2 * Math.PI);

// Stirling-series remainder lnGamma(x) - [(x - 0.5) ln x - x + ln sqrt(2 pi)],
// accurate to ~1e-16 for x >= 15.
const stirlingCorrection = (x) => {
    const r = 1 / x;
    const r2 = r * r;
    return r * (1 / 12 - r2 * (1 / 360 - r2 * (1 / 1260 - r2 * (1 / 1680 - r2 * (1 / 1188 - r2 * (691 / 360360))))));
};

export const logGamma = (x) => {
    if (!(x > 0)) throw new Error('logGamma needs x > 0');
    let shift = 0;
    let z = x;
    while (z < 15) { shift += Math.log(z); z += 1; }
    return (z - 0.5) * Math.log(z) - z + LOG_SQRT_2PI + stirlingCorrection(z) - shift;
};

// log(x^a e^-x / Gamma(a)), the prefactor of both incomplete-gamma branches,
// written for a >= 15 as -a*(u - log1p(u)) + ... (u = x/a - 1) so the large
// cancelling terms never meet in floating point.
const logGammaPrefactor = (a, x) => {
    if (a < 15) return a * Math.log(x) - x - logGamma(a);
    const u = (x - a) / a;
    let uMinusLog1p;
    if (Math.abs(u) < 0.25) {
        // u - log1p(u) = sum_{k>=2} (-1)^k u^k / k
        let term = u * u;
        let sum = 0;
        for (let k = 2; k < 200; k += 1) {
            const piece = term / k;
            sum += k % 2 === 0 ? piece : -piece;
            if (Math.abs(piece) < 1e-17 * Math.abs(sum)) break;
            term *= u;
        }
        uMinusLog1p = sum;
    } else {
        uMinusLog1p = u - Math.log1p(u);
    }
    return -a * uMinusLog1p + 0.5 * Math.log(a) - LOG_SQRT_2PI - stirlingCorrection(a);
};

const GAMMA_EPS = 1e-16;
const GAMMA_MAX_ITER = 100000;

// Regularized incomplete gamma: { lower: P(a, x), upper: Q(a, x) }, each side
// computed directly where it is the small one (series for x < a + 1, Lentz
// continued fraction otherwise) so both tails keep full relative accuracy.
export const regularizedGamma = (a, x) => {
    if (!(a > 0) || !(x >= 0)) throw new Error('regularizedGamma needs a > 0, x >= 0');
    if (x === 0) return { lower: 0, upper: 1 };
    if (x === Infinity) return { lower: 1, upper: 0 };
    const logPrefactor = logGammaPrefactor(a, x);
    if (x < a + 1) {
        let ap = a;
        let term = 1 / a;
        let sum = term;
        for (let n = 0; n < GAMMA_MAX_ITER; n += 1) {
            ap += 1;
            term *= x / ap;
            sum += term;
            if (Math.abs(term) < Math.abs(sum) * GAMMA_EPS) break;
        }
        const lower = sum * Math.exp(logPrefactor);
        return { lower, upper: 1 - lower };
    }
    const tiny = 1e-300;
    let b = x + 1 - a;
    let c = 1 / tiny;
    let d = 1 / b;
    let h = d;
    for (let i = 1; i < GAMMA_MAX_ITER; i += 1) {
        const an = -i * (i - a);
        b += 2;
        d = an * d + b;
        if (Math.abs(d) < tiny) d = tiny;
        c = b + an / c;
        if (Math.abs(c) < tiny) c = tiny;
        d = 1 / d;
        const delta = d * c;
        h *= delta;
        if (Math.abs(delta - 1) < GAMMA_EPS) break;
    }
    const upper = Math.exp(logPrefactor) * h;
    return { lower: 1 - upper, upper };
};

// Wichura's AS 241 (PPND16): the standard normal quantile to ~1e-16.
export const normalQuantile = (p) => {
    if (!(p > 0 && p < 1)) {
        if (p === 0) return -Infinity;
        if (p === 1) return Infinity;
        throw new Error('normalQuantile needs 0 <= p <= 1');
    }
    const q = p - 0.5;
    if (Math.abs(q) <= 0.425) {
        const r = 0.180625 - q * q;
        return q * (((((((2509.0809287301227 * r + 33430.57558358813) * r + 67265.7709270087) * r
            + 45921.95393154987) * r + 13731.69376550946) * r + 1971.5909503065513) * r
            + 133.14166789178438) * r + 3.3871328727963665)
            / (((((((5226.495278852546 * r + 28729.085735721943) * r + 39307.89580009271) * r
            + 21213.794301586597) * r + 5394.196021424751) * r + 687.1870074920579) * r
            + 42.31333070160091) * r + 1);
    }
    let r = q < 0 ? p : 1 - p;
    r = Math.sqrt(-Math.log(r));
    let value;
    if (r <= 5) {
        r -= 1.6;
        value = (((((((0.0007745450142783414 * r + 0.022723844989269184) * r + 0.2417807251774506) * r
            + 1.2704582524523684) * r + 3.6478483247632045) * r + 5.769497221460691) * r
            + 4.630337846156546) * r + 1.4234371107496835)
            / (((((((1.0507500716444169e-9 * r + 0.0005475938084995345) * r + 0.015198666563616457) * r
            + 0.14810397642748008) * r + 0.6897673349851) * r + 1.6763848301838038) * r
            + 2.053191626637759) * r + 1);
    } else {
        r -= 5;
        value = (((((((2.0103343992922881e-7 * r + 0.000027115555687434876) * r + 0.0012426609473880784) * r
            + 0.026532189526576124) * r + 0.29656057182850487) * r + 1.7848265399172913) * r
            + 5.463784911164114) * r + 6.657904643501103)
            / (((((((2.0442631033899397e-15 * r + 1.421511758316446e-7) * r + 0.000018463183175100548) * r
            + 0.0007868691311456133) * r + 0.014875361290850615) * r + 0.1369298809227358) * r
            + 0.599832206555888) * r + 1);
    }
    return q < 0 ? -value : value;
};

/**
 * One-sided normal deviate z with P(Z >= z) = upper; `lower` = 1 - upper
 * computed independently so both tails stay accurate. Mirrors `_normal_deviate`.
 */
export const normalDeviate = (upper, lower) => (
    upper <= 0.5
        ? -normalQuantile(Math.max(upper, SIGNIFICANCE_TAIL_FLOOR))
        : normalQuantile(Math.max(lower, SIGNIFICANCE_TAIL_FLOOR))
);

// Look-elsewhere-corrected Poisson significance of one cell's raw count:
// exact upper tail P(X >= n | e) = P(n, e), Sidak-corrected over `trials`
// cells in log space. Mirrors `_peak_significance`.
const peakSignificance = (count, expected, trials) => {
    let local = 1;
    let localLower = 0;
    if (count > 0) {
        const tails = regularizedGamma(count, expected);
        local = tails.lower;
        localLower = tails.upper;
    }
    let logLower;
    if (local < 0.5) logLower = trials * Math.log1p(-local);
    else if (localLower > 0) logLower = trials * Math.log(localLower);
    else logLower = -Infinity;
    const corrected = -Math.expm1(logLower);
    return { local, corrected, deviate: normalDeviate(corrected, Math.exp(logLower)) };
};

// c_k = C(2k, k) / 4^k for k = 0..kmax by the exact recurrence
// c_k = c_{k-1} (2k - 1) / (2k) -- the same sequential product as
// `_central_binomial`, so the two engines agree bit for bit.
const centralBinomial = (kmax) => {
    const table = new Float64Array(kmax + 1);
    table[0] = 1;
    for (let k = 1; k <= kmax; k += 1) table[k] = table[k - 1] * ((2 * k - 1) / (2 * k));
    return table;
};

// A = sum_pairs |n+ - n-| / N and its exact inversion-symmetric null,
// conditional on each antipodal pair's total T: X ~ Bin(T, 1/2), so
// D = |2X - T| has E[D] = T c_{floor(T/2)} and E[D^2] = T. Mirrors
// `_antipodal_asymmetry_test`.
const antipodalAsymmetryTest = (counts, antipode, used) => {
    let maxTotal = 0;
    for (let cell = 0; cell < counts.length; cell += 1) {
        if (cell < antipode[cell]) maxTotal = Math.max(maxTotal, counts[cell] + counts[antipode[cell]]);
    }
    const central = centralBinomial(Math.floor(maxTotal / 2));
    let observed = 0;
    let meanSum = 0;
    let varianceSum = 0;
    for (let cell = 0; cell < counts.length; cell += 1) {
        const other = antipode[cell];
        if (cell >= other) continue;
        const total = counts[cell] + counts[other];
        const pairMean = total * central[Math.floor(total / 2)];
        observed += Math.abs(counts[cell] - counts[other]);
        meanSum += pairMean;
        varianceSum += Math.max(total - pairMean * pairMean, 0);
    }
    const value = observed / used;
    const nullMean = meanSum / used;
    const nullSd = Math.sqrt(varianceSum) / used;
    if (!(nullSd > 0)) return { value, null: nullMean, nullSd, z: null, significant: false };
    const z = (value - nullMean) / nullSd;
    return { value, null: nullMean, nullSd, z, significant: z > ASYMMETRY_FLAG_SIGMA };
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
// `entries` is an array of {owner, point, payload, key}; returns per-owner
// payload lists in CCW order, each cycle rotated to start at its smallest
// integer `key` (unique per owner). The raw atan2 start is not reproducible
// across engines -- an item on the -e1 ray sits at +pi or -pi depending on a
// 1e-17 round-off sign -- so the combinatorial start is what keeps the
// exported tables identical. Mirrors `_angular_order`.
const angularOrder = (centers, entries) => {
    const byOwner = centers.map(() => []);
    entries.forEach(({ owner, point, payload, key }) => {
        const center = centers[owner];
        const { e1, e2 } = tangentBasis(center);
        const radial = dot(point, center);
        const local = [
            point[0] - center[0] * radial,
            point[1] - center[1] * radial,
            point[2] - center[2] * radial
        ];
        byOwner[owner].push({ angle: Math.atan2(dot(local, e2), dot(local, e1)), payload, key });
    });
    return byOwner.map((list) => {
        list.sort((a, b) => a.angle - b.angle);
        let start = 0;
        for (let i = 1; i < list.length; i += 1) if (list[i].key < list[start].key) start = i;
        return [...list.slice(start), ...list.slice(0, start)].map((item) => item.payload);
    });
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
    const nu = Math.trunc(Number(frequency));
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
            polygonEntries.push({ owner, point: circumcenters[t], payload: circumcenters[t], key: t });
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
        edgeEntries.push({ owner: a, point: centers[b], payload: b, key: b });
        edgeEntries.push({ owner: b, point: centers[c], payload: c, key: c });
        edgeEntries.push({ owner: c, point: centers[a], payload: a, key: a });
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

// True where the first non-zero component is negative: for any u != 0
// exactly one of u, -u is flagged (the zero vector is unflagged). Mirrors
// `_opposite_hemisphere`.
const oppositeHemisphere = (u) => (
    u[0] < 0 || (u[0] === 0 && (u[1] < 0 || (u[1] === 0 && u[2] < 0)))
);

// Exactly inversion-equivariant assignment of a unit direction:
// assign(-u) === antipode[assign(u)] for every u, including exact Voronoi
// ties (the crystal axes at odd nu, <111> when 3 does not divide nu), which a
// first-maximum tie-break would otherwise resolve non-antipodally and make a
// centrosymmetric cloud read as one-sided. Mirrors `assign_cells`.
const assignCell = (tiling, u) => (
    oppositeHemisphere(u)
        ? tiling.antipode[assignOne(tiling, [-u[0], -u[1], -u[2]])]
        : assignOne(tiling, u)
);

/** Cell index for each direction ([x, y, z], need not be normalized). */
export const assignCells = (tiling, directions) => (
    directions.map((direction) => assignCell(tiling, normalize(direction)))
);

// `targetPerCell` as the integer the resolution guard divides by. Mirrors
// `_validated_target`.
const validatedTarget = (targetPerCell) => {
    const value = Number(targetPerCell);
    if (!Number.isFinite(value) || value < 1) throw new Error('targetPerCell must be a finite number >= 1');
    return Math.trunc(value);
};

/**
 * Largest frequency whose cells still average `targetPerCell` points -- the
 * over-binning guard, a floor (the average occupancy never drops below the
 * target, except that fewer than 12 * targetPerCell points still get the
 * 12-cell dodecahedron). Mirrors `recommended_frequency`, including the exact
 * integer boundary test, so both engines pick the same tiling.
 */
export const recommendedFrequency = (nPoints, { targetPerCell = DEFAULT_TARGET_PER_CELL, maxFrequency = 24 } = {}) => {
    const target = validatedTarget(targetPerCell);
    if (!Number.isFinite(Number(maxFrequency)) || Math.trunc(Number(maxFrequency)) < MIN_FREQUENCY) {
        throw new Error(`maxFrequency must be a finite number >= ${MIN_FREQUENCY}`);
    }
    if (!Number.isFinite(Number(nPoints))) throw new Error('nPoints must be a finite number');
    const cap = Math.min(Math.trunc(Number(maxFrequency)), MAX_FREQUENCY);
    if (!(nPoints > 0)) return MIN_FREQUENCY;
    const seed = Math.floor(Math.sqrt(Math.max(nPoints / target - 2, 0) / 10));
    let frequency = Math.min(Math.max(seed, MIN_FREQUENCY), cap);
    while (frequency < cap && target * (10 * (frequency + 1) ** 2 + 2) <= nPoints) frequency += 1;
    while (frequency > MIN_FREQUENCY && target * (10 * frequency * frequency + 2) > nPoints) frequency -= 1;
    return frequency;
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
 * for an isotropic cloud, `zScore` is the local (uncorrected) z per cell and
 * `peakSignificance` the look-elsewhere-corrected peak test.
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
    if (!Number.isFinite(Number(minAmplitude))) throw new Error('minAmplitude must be a finite number');
    if (!(Number.isFinite(Number(smoothing)) && Number(smoothing) >= 0)) {
        throw new Error('smoothing must be a finite, non-negative number of passes');
    }
    if (frequency !== null && frequency !== undefined && !Number.isFinite(Number(frequency))) {
        throw new Error(`frequency must lie in [${MIN_FREQUENCY}, ${MAX_FREQUENCY}]`);
    }
    validatedTarget(targetPerCell);
    // A NaN/inf row would otherwise yield NaN PCA axes or a TypeError deep in
    // the cell assignment (and the Python engine used to fail differently on
    // it). Reject it by name instead -- same message as the Python engine.
    let badRows = 0;
    let firstBad = -1;
    vectors.forEach((row, index) => {
        if (!(Number.isFinite(row[0]) && Number.isFinite(row[1]) && Number.isFinite(row[2]))) {
            badRows += 1;
            if (firstBad < 0) firstBad = index;
        }
    });
    if (badRows > 0) {
        throw new Error(`vectors contain ${badRows} non-finite (NaN/inf) row(s); first at row ${firstBad}`);
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
    // Mean |dr| of the atoms that moved into each cell -- the radial-relief
    // quantity (independent of the color weighting, so shape and color carry
    // separate information). Smoothing applies to the numerator and denominator
    // sums, not the ratio, so a smoothed relief stays the mean of the same
    // smoothed population the color shows.
    let amplitudeSumCell = new Float64Array(cellCount);
    directions.forEach((direction, index) => {
        const cell = assignCell(tiling, direction);
        counts[cell] += 1;
        mass[cell] += weights[index];
        amplitudeSumCell[cell] += keptAmplitude[index];
    });
    let countField = Float64Array.from(counts);

    if (smoothing > 0) {
        const passes = Math.trunc(smoothing);
        mass = Float64Array.from(smooth(Array.from(mass), tiling.neighbors, passes));
        amplitudeSumCell = Float64Array.from(smooth(Array.from(amplitudeSumCell), tiling.neighbors, passes));
        countField = Float64Array.from(smooth(Array.from(countField), tiling.neighbors, passes));
    }

    const cellMeanAmplitude = new Array(cellCount);
    for (let cell = 0; cell < cellCount; cell += 1) {
        cellMeanAmplitude[cell] = countField[cell] > 1e-12
            ? amplitudeSumCell[cell] / Math.max(countField[cell], 1e-12)
            : 0;
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

    // Local z from the *raw* counts (smoothing would correlate cells and
    // overstate the evidence) -- a description, not a test; the peak test is
    // peakSignificance below.
    const expected = new Array(cellCount);
    const zScore = new Array(cellCount);
    for (let cell = 0; cell < cellCount; cell += 1) {
        expected[cell] = used * tiling.areas[cell] / (4 * Math.PI);
        zScore[cell] = (counts[cell] - expected[cell]) / Math.sqrt(Math.max(expected[cell], 1e-12));
    }

    // Antipodal (inversion) asymmetry: sum_pairs |n(u) - n(-u)| / N -- 0 for an
    // inversion-symmetric cloud, 1 for a fully one-sided one -- with its exact
    // null conditional on the observed pair totals (antipodalAsymmetryTest).
    const asymmetryTest = antipodalAsymmetryTest(counts, tiling.antipode, used);

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

    // Bingham's test of uniformity: S = 7.5 n_eff |T - tr(T)/3 I|_F^2 ~ chi^2_5,
    // n_eff = (sum w)^2 / sum w^2. Mirrors the Python engine.
    let weightSquares = 0;
    weights.forEach((w) => { weightSquares += w * w; });
    const effectivePoints = weightSum * weightSum / weightSquares;
    const traceThird = (tensor[0][0] + tensor[1][1] + tensor[2][2]) / 3;
    let deviatorSquares = 0;
    for (let a = 0; a < 3; a += 1) {
        for (let b = 0; b < 3; b += 1) {
            const value = tensor[a][b] - (a === b ? traceThird : 0);
            deviatorSquares += value * value;
        }
    }
    const binghamStatistic = 7.5 * effectivePoints * deviatorSquares;
    const binghamTails = regularizedGamma(2.5, binghamStatistic / 2);

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
    // Tie-tolerant argmax: the lowest index within PEAK_TIE_RTOL of the
    // maximum, so round-off in symmetry-equal cell areas cannot pick the peak.
    let peak = -1;
    let peakTieCount = 0;
    const tieFloor = vmax * (1 - PEAK_TIE_RTOL);
    for (let cell = 0; cell < cellCount; cell += 1) {
        if (enhancement[cell] >= tieFloor) {
            if (peak < 0) peak = cell;
            peakTieCount += 1;
        }
    }
    const peakTest = peakSignificance(counts[peak], expected[peak], cellCount);

    // Whole-map test: Pearson's X^2 = sum z^2 vs chi^2 with C - 1 degrees of
    // freedom, as a one-sided normal deviate. Mirrors the Python engine.
    const degreesOfFreedom = cellCount - 1;
    const mapTails = regularizedGamma(degreesOfFreedom / 2, zSquares / 2);
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
        cellMeanAmplitude,
        antipodalAsymmetry: asymmetryTest.value,
        antipodalAsymmetryNull: asymmetryTest.null,
        antipodalAsymmetryNullSd: asymmetryTest.nullSd,
        antipodalAsymmetryZ: asymmetryTest.z,
        antipodalAsymmetrySignificant: asymmetryTest.significant,
        orientationTensor: tensor,
        orientationEigenvalues: tensorDecomposition.eigenvalues,
        orientationAxes: tensorDecomposition.axes,
        orientationAnisotropy: 3 * tensorDecomposition.eigenvalues[0] - 1,
        orientationEffectivePoints: effectivePoints,
        orientationAnisotropyNull: ISOTROPIC_ANISOTROPY_SCALE / Math.sqrt(effectivePoints),
        orientationBinghamStatistic: binghamStatistic,
        orientationBinghamPValue: binghamTails.upper,
        orientationAnisotropySignificance: normalDeviate(binghamTails.upper, binghamTails.lower),
        peakCell: peak,
        peakTieCount,
        peakDirection: tiling.centers[peak],
        peakEnhancement: enhancement[peak],
        // Local, uncorrected Gaussian z -- kept for API compatibility; the
        // calibrated readout is peakSignificance (see _peak_significance).
        peakZScore: zScore[peak],
        peakCount: counts[peak],
        peakExpected: expected[peak],
        peakLocalPValue: peakTest.local,
        peakPValue: peakTest.corrected,
        peakSignificance: peakTest.deviate,
        // Legacy RMS of the local z (1 +/- 1/sqrt(2C) for noise) -- NOT a
        // sigma level; the calibrated readout is mapSignificance.
        significance: Math.sqrt(zSquares / cellCount),
        mapChiSquare: zSquares,
        mapDegreesOfFreedom: degreesOfFreedom,
        mapPValue: mapTails.upper,
        mapSignificance: normalDeviate(mapTails.upper, mapTails.lower),
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
    // '' and 'all' are the pooled-everything default, exactly as in Python.
    const pooledAll = element === null || element === undefined || element === '' || element === 'all';
    if (referenceNumber !== null) {
        tagged = parsed.sites.find((site) => site.referenceNumber === referenceNumber);
        if (!tagged) throw new Error(`Unknown reference number ${referenceNumber}`);
        cloud = tagged.displacements;
    } else if (!pooledAll) {
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
    } else if (!pooledAll) {
        result.element = String(element);
    }
    return result;
};
