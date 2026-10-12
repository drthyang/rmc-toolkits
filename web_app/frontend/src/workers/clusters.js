// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Cluster distortion engine — static-mode port of rmc_toolkits/clusters.py
// (the source of truth; keep the two in sync, the parity fixture pins them).
//
// A cluster is a finite group of one element's atoms bonded within a length
// window in the AVERAGE structure (the Ta4 tetrahedra of GaTa4Se8, a dimer, a
// trimer). Templates come from the average site positions, so a cluster in
// the configuration is always the same atoms however far its bonds stretch.
// Its edges are measured in the configuration and compared with a shuffle
// null: each member keeps its own site's displacement distribution but takes
// it from a random other cluster, which removes only the correlation between
// the atoms of one cluster. Each statistic comes back as observed, null mean,
// null sd over `shuffles` shuffles and z = (observed − mean) / sd. The
// shuffles are counter-based hash permutations, identical in both engines.

import { eigenDecomposition, siteSpecies } from './pcaKde.js';
import { WINDOW_TOL, countPatterns, perpendicularWidths } from './triplets.js';

export const DEFAULT_SHUFFLES = 24;
export const SHUFFLE_SEED = 0x5eed;

// The site's average position consistent with its atoms' cell indices
// (pcaKde.js `siteCellFractional`, unwrapped): a cluster pairs the atoms of
// given cells, so a site at fraction 0 must not wrap to 0.9999 and pull its
// partners from the next cell. clusters.py `_cell_positions`.
const cellPosition = (site) => site.siteCellFractional ?? site.siteFractional;
export const AXIS_CONE_DEG = 15;
// Opposite edge pairs of a tetrahedron, as indices into its sorted edge list
// (0,1) (0,2) (0,3) (1,2) (1,3) (2,3).
export const TETRAHEDRON_OPPOSITE = [[0, 5], [1, 4], [2, 3]];

/** The lowbias32 integer hash (C. Wellons), on unsigned 32-bit integers. */
export const mix32 = (value) => {
    let x = value >>> 0;
    x = Math.imul(x ^ (x >>> 16), 0x7feb352d) >>> 0;
    x = Math.imul(x ^ (x >>> 15), 0x846ca68b) >>> 0;
    return (x ^ (x >>> 16)) >>> 0;
};

/** A permutation of 0..n−1 from hashed keys, ties by index (shuffle_permutation). */
export const shufflePermutation = (n, shuffle, template, member, seed = SHUFFLE_SEED) => {
    let head = mix32((seed ^ mix32(shuffle + 1)) >>> 0);
    head = mix32((head ^ (template + 1)) >>> 0);
    head = mix32((head ^ (member + 1)) >>> 0);
    const keys = new Float64Array(n);
    for (let c = 0; c < n; c += 1) keys[c] = mix32((head ^ c) >>> 0);
    return Array.from({ length: n }, (_, c) => c).sort((a, b) => (keys[a] - keys[b]) || (a - b));
};

const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const length = (v) => Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
// Row vector × matrix (rows = basis vectors), as numpy's `frac @ unit`.
const toCartesian = (f, unit) => [0, 1, 2].map((k) => f[0] * unit[0][k] + f[1] * unit[1][k] + f[2] * unit[2][k]);
const mean = (values) => values.reduce((sum, value) => sum + value, 0) / values.length;
const median = (values) => {
    const sorted = [...values].sort((a, b) => a - b);
    const mid = sorted.length >> 1;
    return sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
};
const capitalize = (symbol) => (symbol ? symbol[0].toUpperCase() + symbol.slice(1).toLowerCase() : symbol);

/** Unit vectors of the <100>, <110> and <111> lattice directions (Cartesian). */
export const latticeFamilies = (unit) => {
    const [a, b, c] = unit;
    const neg = (v) => v.map((x) => -x);
    const groups = {
        100: [a, b, c],
        110: [add(a, b), sub(a, b), add(a, c), sub(a, c), add(b, c), sub(b, c)],
        111: [add(add(a, b), c), sub(add(a, b), c), add(sub(a, b), c), add(add(neg(a), b), c)]
    };
    return Object.fromEntries(Object.entries(groups).map(([name, vectors]) => [
        name, vectors.map((v) => v.map((x) => x / length(v)))
    ]));
};

const sameShift = (a, b) => a[0] === b[0] && a[1] === b[1] && a[2] === b[2];

// Finite clusters of the given sites under the window (mirrors _templates).
const buildTemplates = (sites, members, unit, [lo, hi]) => {
    const widths = perpendicularWidths(unit);
    const reach = widths.map((width) => Math.ceil(hi / width) + 1);
    const shifts = [];
    for (let i = -reach[0]; i <= reach[0]; i += 1) {
        for (let j = -reach[1]; j <= reach[1]; j += 1) {
            for (let k = -reach[2]; k <= reach[2]; k += 1) shifts.push([i, j, k]);
        }
    }
    const lower = Math.max(lo - WINDOW_TOL, 0);
    const upper = hi + WINDOW_TOL;
    const adjacency = new Map(members.map((s) => [s, []]));
    let extended = false;
    for (const s of members) {
        for (const t of members) {
            for (const shift of shifts) {
                const d = length(toCartesian(add(sub(cellPosition(sites[t]), cellPosition(sites[s])), shift), unit));
                if (!(d >= lower && d <= upper)) continue;
                if (s === t) {
                    if (shift.some((value) => value !== 0)) extended = true;
                    continue;
                }
                adjacency.get(s).push([t, shift]);
            }
        }
    }
    const templates = [];
    const seen = new Set();
    let isolated = 0;
    for (const anchor of members) {
        if (seen.has(anchor)) continue;
        if (!adjacency.get(anchor).length) {
            seen.add(anchor);
            isolated += 1;
            continue;
        }
        const offsets = new Map([[anchor, [0, 0, 0]]]);
        const queue = [anchor];
        while (queue.length) {
            const s = queue.shift();
            for (const [t, shift] of adjacency.get(s)) {
                const reached = add(offsets.get(s), shift);
                if (offsets.has(t)) {
                    if (!sameShift(offsets.get(t), reached)) extended = true;
                    continue;
                }
                offsets.set(t, reached);
                queue.push(t);
            }
        }
        const order = [...offsets.keys()].sort((a, b) => a - b);
        order.forEach((s) => seen.add(s));
        const index = new Map(order.map((s, i) => [s, i]));
        const edgeKeys = new Set();
        for (const s of order) {
            for (const [t, shift] of adjacency.get(s)) {
                if (index.has(t) && sameShift(sub(offsets.get(t), offsets.get(s)), shift)) {
                    const [i, j] = [index.get(s), index.get(t)].sort((a, b) => a - b);
                    edgeKeys.add(`${i},${j}`);
                }
            }
        }
        const edges = [...edgeKeys].map((key) => key.split(',').map(Number))
            .sort((a, b) => (a[0] - b[0]) || (a[1] - b[1]));
        templates.push({ sites: order, offsets: order.map((s) => offsets.get(s)), edges });
    }
    return { templates, extended, isolated };
};

const cellLinear = (cell, supercell) => {
    const wrapped = cell.map((value, axis) => ((value % supercell[axis]) + supercell[axis]) % supercell[axis]);
    return (wrapped[0] * supercell[1] + wrapped[1]) * supercell[2] + wrapped[2];
};

// Shape statistics of clusters given as member positions [cluster][member] (mirrors _statistics).
const clusterStatistics = (positions, edges, families) => {
    const nClusters = positions.length;
    const members = positions[0].length;
    const lengths = positions.map((cluster) => edges.map(([i, j]) => length(sub(cluster[j], cluster[i]))));
    const ladders = lengths.map((row) => [...row].sort((a, b) => a - b));
    const rowMeans = lengths.map(mean);
    let total = 0;
    lengths.forEach((row) => row.forEach((value) => { total += value; }));
    const nEdges = edges.length;
    const ladder = Array.from({ length: nEdges }, (_, k) => ladders.reduce((sum, row) => sum + row[k], 0) / nClusters);
    let shape = 0;
    lengths.forEach((row, c) => row.forEach((value) => { shape += (value - rowMeans[c]) ** 2; }));
    const grandMean = mean(rowMeans);
    const stats = {
        meanEdge: total / (nClusters * nEdges),
        ladder,
        ladderWidth: ladders.reduce((sum, row) => sum + (row[nEdges - 1] - row[0]), 0) / nClusters,
        shapeRms: Math.sqrt(shape / (nClusters * nEdges)),
        breathingSd: Math.sqrt(rowMeans.reduce((sum, value) => sum + (value - grandMean) ** 2, 0) / nClusters)
    };
    if (members === 4 && nEdges === 6) {
        let eSum = 0;
        let tSum = 0;
        let trig = 0;
        let valid = 0;
        lengths.forEach((row) => {
            const sums = TETRAHEDRON_OPPOSITE.map(([a, b]) => row[a] + row[b]);
            const diffs = TETRAHEDRON_OPPOSITE.map(([a, b]) => row[a] - row[b]);
            const sumMean = mean(sums);
            eSum += sums.reduce((acc, value) => acc + (value - sumMean) ** 2, 0);
            const norm2 = diffs.reduce((acc, value) => acc + value * value, 0);
            tSum += norm2;
            if (norm2 > 0) {
                trig += diffs.reduce((acc, value) => acc + value ** 4, 0) / (norm2 * norm2);
                valid += 1;
            }
        });
        stats.eRms = Math.sqrt(eSum / nClusters / 2);
        stats.t2Rms = Math.sqrt(tSum / nClusters / 3);
        stats.t2Trigonality = valid ? trig / valid : null;
    }
    if (members >= 2) {
        const cone = Math.cos(AXIS_CONE_DEG * Math.PI / 180);
        const axes = { longAxis: [], shortAxis: [] };
        positions.forEach((cluster) => {
            if (members === 2) {
                axes.longAxis.push(sub(cluster[1], cluster[0]));
                return;
            }
            const centre = [0, 1, 2].map((k) => cluster.reduce((sum, p) => sum + p[k], 0) / members);
            const gyration = [0, 1, 2].map((i) => [0, 1, 2].map((j) => (
                cluster.reduce((sum, p) => sum + (p[i] - centre[i]) * (p[j] - centre[j]), 0) / members
            )));
            const { axes: principal } = eigenDecomposition(gyration);
            axes.longAxis.push(principal[0]);
            axes.shortAxis.push(principal[2]);
        });
        Object.entries(axes).forEach(([name, vectors]) => {
            if (!vectors.length) return;
            stats[name] = Object.fromEntries(Object.entries(families).map(([family, directions]) => {
                let hits = 0;
                vectors.forEach((vector) => {
                    const norm = length(vector);
                    const best = Math.max(...directions.map((d) => Math.abs(dot(vector, d)) / norm));
                    if (best >= cone) hits += 1;
                });
                return [family, hits / vectors.length];
            }));
        });
    }
    return stats;
};

// Mean, sd (n − 1) and z of every statistic, elementwise (mirrors _null_summary).
const nullSummary = (observed, nulls) => {
    const combine = (values, obs) => {
        const kept = values.filter((value) => value != null);
        if (obs == null || kept.length < 2) return [null, null, null];
        const m = mean(kept);
        const sd = Math.sqrt(kept.reduce((sum, value) => sum + (value - m) ** 2, 0) / (kept.length - 1));
        return [m, sd, sd > 0 ? (obs - m) / sd : null];
    };
    const means = {};
    const sds = {};
    const zs = {};
    Object.entries(observed).forEach(([key, obs]) => {
        if (Array.isArray(obs)) {
            const parts = obs.map((value, k) => combine(nulls.map((n) => n[key][k]), value));
            [means[key], sds[key], zs[key]] = [0, 1, 2].map((i) => parts.map((part) => part[i]));
        } else if (obs && typeof obs === 'object') {
            const parts = Object.fromEntries(Object.keys(obs).map((k) => [k, combine(nulls.map((n) => n[key][k]), obs[k])]));
            [means[key], sds[key], zs[key]] = [0, 1, 2].map((i) => Object.fromEntries(
                Object.entries(parts).map(([k, part]) => [k, part[i]])
            ));
        } else {
            [means[key], sds[key], zs[key]] = combine(nulls.map((n) => n[key]), obs);
        }
    });
    return { means, sds, zs };
};

const invert3 = (m) => {
    const [[a, b, c], [d, e, f], [g, h, i]] = m;
    const A = e * i - f * h;
    const B = -(d * i - f * g);
    const C = d * h - e * g;
    const det = a * A + b * B + c * C;
    return [
        [A / det, -(b * i - c * h) / det, (b * f - c * e) / det],
        [B / det, (a * i - c * g) / det, -(a * f - c * d) / det],
        [C / det, -(a * h - b * g) / det, (a * e - b * d) / det]
    ];
};

/**
 * Cluster templates, observed shape statistics and their shuffle null, from a
 * parsed configuration (siteDisplacementsFromRmc6f). Mirrors
 * rmc_toolkits.clusters.cluster_distortion; camelCase payload.
 */
export const clusterDistortion = (parsed, element, [rmin, rmax], { shuffles = DEFAULT_SHUFFLES, seed = SHUFFLE_SEED } = {}) => {
    const lo = Number(rmin);
    const hi = Number(rmax);
    if (!Number.isFinite(lo) || !Number.isFinite(hi) || lo < 0 || hi <= lo) {
        throw new Error(`window must be 0 <= min < max, got [${rmin}, ${rmax}]`);
    }
    const sites = parsed.sites;
    if (parsed.reconstructed || sites.some((site) => !Array.isArray(site.cells))) {
        throw new Error('cluster analysis needs the reference-site and cell columns of the .rmc6f file');
    }
    if (!(shuffles >= 2)) throw new Error(`shuffles must be at least 2, got ${shuffles}`);
    const symbol = capitalize(String(element || '').trim());
    const members = sites.map((site, s) => (site.element === symbol ? s : -1)).filter((s) => s >= 0);
    if (!members.length) throw new Error(`no site holds ${symbol}; elements: ${siteSpecies(sites).join(', ')}`);

    const supercell = parsed.supercell.map(Number);
    const lattice = parsed.latticeVectors;
    const unit = lattice.map((row, i) => row.map((value) => value / supercell[i]));
    const families = latticeFamilies(unit);
    const nCells = supercell[0] * supercell[1] * supercell[2];
    // rowOf[s][cell] = the site's local atom index in that cell, −1 if none.
    const rowOf = new Map(members.map((s) => {
        const table = new Int32Array(nCells).fill(-1);
        sites[s].cells.forEach((cell, k) => { table[cellLinear(cell, supercell)] = k; });
        return [s, table];
    }));

    const { templates, extended, isolated } = buildTemplates(sites, members, unit, [lo, hi]);
    const result = {
        element: symbol,
        window: [lo, hi],
        shuffles,
        sites: members.map((s) => sites[s].referenceNumber),
        extended,
        isolatedSites: isolated,
        templates: [],
        kinds: [],
        parseWarning: parsed.parseWarning ?? null
    };
    if (extended || !templates.length) return result;

    const cellOf = (c) => [Math.floor(c / (supercell[1] * supercell[2])), Math.floor(c / supercell[2]) % supercell[1], c % supercell[2]];
    const perTemplate = templates.map((template, tIndex) => {
        const rowsAll = template.sites.map((s, m) => {
            const table = rowOf.get(s);
            return Array.from({ length: nCells }, (_, c) => table[cellLinear(add(cellOf(c), template.offsets[m]), supercell)]);
        });
        const complete = [];
        for (let c = 0; c < nCells; c += 1) if (rowsAll.every((rows) => rows[c] >= 0)) complete.push(c);
        const rows = rowsAll.map((all) => complete.map((c) => all[c]));
        const base = template.sites.map((s, m) => toCartesian(add(cellPosition(sites[s]), template.offsets[m]), unit));
        const idealEdges = template.edges.map(([i, j]) => length(sub(base[j], base[i])));
        const edgeSe = template.edges.map(([i, j]) => {
            const span = sub(base[j], base[i]);
            const direction = span.map((x) => x / length(span));
            let variance = 0;
            for (const s of [template.sites[i], template.sites[j]]) {
                const cloud = sites[s].displacements.map((u) => dot(u, direction));
                variance += mean(cloud.map((x) => x * x)) / cloud.length;
            }
            return Math.sqrt(variance);
        });
        result.templates.push({
            referenceNumbers: template.sites.map((s) => sites[s].referenceNumber),
            offsets: template.offsets,
            edges: template.edges,
            idealEdges,
            edgeSe,
            clusters: complete.length,
            incomplete: nCells - complete.length
        });
        return { tIndex, template, rows, base };
    });

    const kinds = new Map();
    perTemplate.forEach(({ template, rows }) => {
        const key = `${template.sites.length},${template.edges.length}`;
        kinds.set(key, (kinds.get(key) || 0) + rows[0].length);
    });
    result.kinds = [...kinds.entries()]
        .map(([key, clusters]) => {
            const [nMembers, nEdges] = key.split(',').map(Number);
            return { members: nMembers, edges: nEdges, clusters };
        })
        .sort((a, b) => b.clusters - a.clusters);
    let dominant = null;
    for (const [key, clusters] of kinds) {
        const nMembers = Number(key.split(',')[0]);
        if (!dominant || clusters > dominant.clusters || (clusters === dominant.clusters && nMembers < dominant.members)) {
            dominant = { key, clusters, members: nMembers };
        }
    }
    const chosen = perTemplate.filter(({ template }) => `${template.sites.length},${template.edges.length}` === dominant.key);
    const [kindMembers, kindEdges] = dominant.key.split(',').map(Number);
    result.kind = { members: kindMembers, edges: kindEdges, clusters: dominant.clusters, templates: chosen.length };

    const positionsFor = (shuffle) => {
        const out = [];
        chosen.forEach(({ tIndex, template, rows, base }) => {
            const n = rows[0].length;
            const orders = template.sites.map((_, m) => (
                shuffle == null ? null : shufflePermutation(n, shuffle, tIndex, m, seed)
            ));
            for (let c = 0; c < n; c += 1) {
                out.push(template.sites.map((s, m) => add(base[m], sites[s].displacements[rows[m][orders[m] ? orders[m][c] : c]])));
            }
        });
        return out;
    };
    const edges = chosen[0].template.edges;
    const observedPositions = positionsFor(null);
    const observed = clusterStatistics(observedPositions, edges, families);
    const nulls = Array.from({ length: shuffles }, (_, k) => clusterStatistics(positionsFor(k), edges, families));
    const { means, sds, zs } = nullSummary(observed, nulls);
    Object.assign(result, { observed, nullMean: means, nullSd: sds, z: zs });

    const chosenTemplates = chosen.map(({ tIndex }) => result.templates[tIndex]);
    const sortedIdeal = chosenTemplates.map((t) => [...t.idealEdges].sort((a, b) => a - b));
    const spreads = chosenTemplates.map((t) => Math.max(...t.idealEdges) - Math.min(...t.idealEdges));
    const errors = chosenTemplates.map((t) => median(t.edgeSe));
    result.ideal = {
        ladder: sortedIdeal[0].map((_, k) => mean(sortedIdeal.map((row) => row[k]))),
        spread: Math.max(...spreads),
        edgeSe: median(errors),
        spreadOverSe: Math.max(...spreads.map((spread, i) => spread / Math.max(errors[i], 1e-12)))
    };

    // Outliers: template edges outside the window; contacts in the window that
    // are not template bonds, over every atom of the element.
    const lower = Math.max(lo - WINDOW_TOL, 0);
    const upper = hi + WINDOW_TOL;
    const insideWindow = (value) => value >= lower && value <= upper;
    let inside = 0;
    observedPositions.forEach((cluster) => edges.forEach(([i, j]) => { if (insideWindow(length(sub(cluster[j], cluster[i])))) inside += 1; }));
    const inverse = invert3(lattice);
    const fractional = [];
    members.forEach((s) => {
        const site = sites[s];
        site.displacements.forEach((u, k) => {
            const shift = toCartesian(u, inverse);
            fractional.push([0, 1, 2].map((axis) => {
                const value = (cellPosition(site)[axis] + site.cells[k][axis]) / supercell[axis] + shift[axis];
                return ((value % 1) + 1) % 1;
            }));
        });
    });
    const rowsLocal = fractional.map((_, index) => index);
    const directed = countPatterns(fractional, rowsLocal, rowsLocal, lattice, [[lo, hi]]).reduce((sum, value) => sum + value, 0);
    let templateBondsInside = 0;
    perTemplate.forEach(({ template, rows, base }) => {
        for (let c = 0; c < rows[0].length; c += 1) {
            const cluster = template.sites.map((s, m) => add(base[m], sites[s].displacements[rows[m][c]]));
            template.edges.forEach(([i, j]) => { if (insideWindow(length(sub(cluster[j], cluster[i])))) templateBondsInside += 1; });
        }
    });
    const extra = Math.floor(directed / 2) - templateBondsInside;
    result.outliers = {
        brokenEdgeFraction: 1 - inside / (observedPositions.length * edges.length),
        extraContacts: extra,
        extraContactsPerAtom: fractional.length ? 2 * extra / fractional.length : 0
    };
    return result;
};
