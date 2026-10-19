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
// The engine also decides what its statistics show: `verdict`, always present
// (the verdict section below, clusters.py's), which the Copilot's
// analyze_clusters only renders.

import { eigenDecomposition, siteEllipsoids, siteSpecies } from './pcaKde.js';
import { WINDOW_TOL, countPatterns, perpendicularWidths } from './triplets.js';

export const DEFAULT_SHUFFLES = 24;
export const SHUFFLE_SEED = 0x5eed;

// The site's average position consistent with its atoms' cell indices
// (pcaKde.js `siteCellFractional`, unwrapped): a cluster pairs the atoms of
// given cells, so a site at fraction 0 must not wrap to 0.9999 and pull its
// partners from the next cell. clusters.py `_cell_positions`.
const cellPosition = (site) => site.siteCellFractional ?? site.siteFractional;
export const AXIS_CONE_DEG = 15;
// The size statistics: read against the independent-atom null only (the
// shape null keeps every cluster's own size).
export const SIZE_STATISTICS = ['meanEdge', 'breathingSd'];
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
// `edges[c]` holds cluster c's member index pairs: templates of one kind need
// not share their edges (two chains of three, one bonded 0-1-2 and one 1-0-2).
// Mirrors _statistics / _edge_vectors.
const clusterStatistics = (positions, edges, families) => {
    const nClusters = positions.length;
    const members = positions[0].length;
    const lengths = positions.map((cluster, c) => edges[c].map(([i, j]) => length(sub(cluster[j], cluster[i]))));
    const ladders = lengths.map((row) => [...row].sort((a, b) => a - b));
    const rowMeans = lengths.map(mean);
    let total = 0;
    lengths.forEach((row) => row.forEach((value) => { total += value; }));
    const nEdges = lengths[0].length;
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
    if (extended || !templates.length) {
        result.verdict = noClusterVerdict(extended ? 'extended_network' : 'no_clusters');
        return result;
    }

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
        return { tIndex, template, rows, rowsAll, base };
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
        // Most clusters first, ties to the smaller kind: the dominant kind leads.
        .sort((a, b) => (b.clusters - a.clusters) || (a.members - b.members));
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
    // Each cluster's breathing: its least-squares scale about the template
    // centroid, ε_c = Σ_m d_mc·r_m / Σ_m |r_m|² (r_m = base_m − centroid), and
    // the displacements with it taken out, d'_mc = d_mc − ε_c r_m. The fit also
    // takes out a share β_m = |r_m|² / Σ_k |r_k|² of each member's own radial
    // noise, so the residual's radial component is scaled back by
    // 1/√(1 − β_m) (clusters.py: exact for members of equal radial variance).
    const breathing = chosen.map(({ template, rows, base }) => {
        const centroid = [0, 1, 2].map((k) => mean(base.map((point) => point[k])));
        const radial = base.map((point) => sub(point, centroid));
        const lengths2 = radial.map((r) => dot(r, r));
        const norm2 = lengths2.reduce((sum, value) => sum + value, 0);
        const unit = radial.map((r, m) => (lengths2[m] > 0 ? r.map((x) => x / Math.sqrt(lengths2[m])) : [0, 0, 0]));
        const boost = lengths2.map((value) => 1 / Math.sqrt(1 - value / norm2) - 1);
        const n = rows[0].length;
        const scale = new Float64Array(n);
        const residual = [];
        for (let c = 0; c < n; c += 1) {
            const displaced = template.sites.map((s, m) => sites[s].displacements[rows[m][c]]);
            const eps = displaced.reduce((sum, d, m) => sum + dot(d, radial[m]), 0) / norm2;
            scale[c] = eps;
            residual.push(displaced.map((d, m) => {
                const left = sub(d, radial[m].map((x) => eps * x));
                const along = dot(left, unit[m]);
                return add(left, unit[m].map((x) => boost[m] * along * x));
            }));
        }
        return { radial, scale, residual };
    });
    // The shape null: every cluster keeps exactly its own breathing, and its
    // members' breathing-free displacements come from random other clusters
    // (the same permutations as the independent-atom null), with the
    // breathing their combination happens to carry taken out again.
    const shapePositionsFor = (shuffle) => {
        const out = [];
        chosen.forEach(({ tIndex, template, rows, base }, index) => {
            const { radial, scale, residual } = breathing[index];
            const norm2 = radial.reduce((sum, r) => sum + dot(r, r), 0);
            const n = rows[0].length;
            const orders = template.sites.map((_, m) => shufflePermutation(n, shuffle, tIndex, m, seed));
            for (let c = 0; c < n; c += 1) {
                const mixed = template.sites.map((_, m) => residual[orders[m][c]][m]);
                const stray = mixed.reduce((sum, d, m) => sum + dot(d, radial[m]), 0) / norm2;
                out.push(mixed.map((d, m) => add(add(base[m], d), radial[m].map((x) => (scale[c] - stray) * x))));
            }
        });
        return out;
    };
    // Each cluster is measured with its own template's edges.
    const edges = chosen.flatMap(({ template, rows }) => rows[0].map(() => template.edges));
    const observedPositions = positionsFor(null);
    const observed = clusterStatistics(observedPositions, edges, families);
    const nulls = Array.from({ length: shuffles }, (_, k) => clusterStatistics(positionsFor(k), edges, families));
    const { means, sds, zs } = nullSummary(observed, nulls);
    const shape = nullSummary(observed, Array.from({ length: shuffles }, (_, k) => clusterStatistics(shapePositionsFor(k), edges, families)));
    // The shape null holds each cluster's size by construction: its size
    // statistics are no test, and are left out.
    [shape.means, shape.sds, shape.zs].forEach((summary) => SIZE_STATISTICS.forEach((key) => { delete summary[key]; }));
    Object.assign(result, {
        observed, nullMean: means, nullSd: sds, z: zs,
        shapeNull: { mean: shape.means, sd: shape.sds, z: shape.zs }
    });

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
    observedPositions.forEach((cluster, c) => edges[c].forEach(([i, j]) => { if (insideWindow(length(sub(cluster[j], cluster[i])))) inside += 1; }));
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
    // Every template bond whose two atoms are there, in incomplete clusters
    // too: they are not contacts between clusters (clusters.py).
    perTemplate.forEach(({ template, rowsAll, base }) => {
        template.edges.forEach(([i, j]) => {
            const si = template.sites[i];
            const sj = template.sites[j];
            for (let c = 0; c < nCells; c += 1) {
                const ri = rowsAll[i][c];
                const rj = rowsAll[j][c];
                if (ri < 0 || rj < 0) continue;
                const span = sub(add(base[j], sites[sj].displacements[rj]), add(base[i], sites[si].displacements[ri]));
                if (insideWindow(length(span))) templateBondsInside += 1;
            }
        });
    });
    const extra = Math.floor(directed / 2) - templateBondsInside;
    result.outliers = {
        brokenEdgeFraction: 1 - inside / (observedPositions.length * kindEdges),
        extraContacts: extra,
        extraContactsPerAtom: fractional.length ? 2 * extra / fractional.length : 0
    };

    // The analysed kind's member sites, from the clouds this engine holds:
    // their site rows decide which are split (one rule, isSplitSite).
    const memberReferences = [...new Set(chosen.flatMap(({ template }) => template.sites.map((s) => sites[s].referenceNumber)))]
        .sort((a, b) => a - b);
    const memberRows = siteEllipsoids(memberReferences.map((reference) => sites.find((site) => site.referenceNumber === reference)));
    result.verdict = clusterVerdict(result, memberRows);
    return result;
};

// --- the verdict (clusters.py's verdict section; keep the two in sync) -----------
//
// Each question is a Finding, { answer, because, evidence }: answer true/false,
// or null where the question does not apply or its statistic is undefined
// (never "not significant"); because, a code from CLUSTER_BECAUSE; evidence,
// every statistic the answer consulted, those that clear their bar first, then
// the rest in a fixed order. An Evidence record is
//
//   { statistic, against, observed, null, z, change, changeKind, bar, clears, side?, strict?, site? }
//
// against: 'shapeNull' | 'independentNull' | 'samplingError' | 'siteKurtosis' |
// 'window'; change = observed / null − 1 ('ratio') or observed − null
// ('difference'); bar maps a field of the record to its threshold, side the
// way it is cleared ('above' ≥ bar, 'below' ≤ −bar, absent |value| ≥ bar;
// strict: > instead of ≥), and clears = evidenceClears(record). Decisions use
// the unrounded statistics; every number is finite or null. The rules were the
// Copilot's summarizeClusters' until 2026-10-19, ported rule for rule.

/** A cluster statistic departs from its null only when both its z and its effect size clear these. */
export const CLUSTER_Z = 3;
/** Relative change (observed / null − 1) of a width or an rms. */
export const CLUSTER_EFFECT = 0.05;
/** Change (observed − null) of the tetrahedral T2 trigonality (range 1/3 to 1). */
export const TRIGONALITY_EFFECT = 0.02;
/** Observed / null fraction of clusters with a gyration axis near one lattice family (one-sided). */
export const AXIS_ENHANCEMENT = 1.1;
/** The average cluster is statically distorted when its edges spread over at least this many sampling errors. */
export const STATIC_SPREAD_SE = 5;
/** Template bonds outside the window, or contacts linking clusters per atom, above this: the window. */
export const CLUSTER_OUTLIER_FRACTION = 0.05;
/** Clusters of other kinds (sizes) above this share of all: the window may cut through or merge clusters. */
export const MIXED_SIZE_SHARE = 0.1;
/**
 * A site is split (two positions) when its PC1 is resolved and its PC1 excess
 * kurtosis is at most −0.5 — two equal wells ±δ of thermal width σ separate at
 * δ = σ (κ = −2δ⁴/(δ² + σ²)² = −0.5) — and at least SPLIT_SITE_SE sampling
 * errors (√(24/n)) below zero. The one rule for the app: list_sites and the
 * split member sites (clusters.py's SPLIT_SITE_KURTOSIS, SPLIT_SITE_SE).
 */
export const SPLIT_SITE_KURTOSIS = -0.5;
export const SPLIT_SITE_SE = 3;

export const CLUSTER_RULES = Object.freeze({
    z: CLUSTER_Z,
    effect: CLUSTER_EFFECT,
    trigonalityEffect: TRIGONALITY_EFFECT,
    axisEnhancement: AXIS_ENHANCEMENT,
    staticSpreadSe: STATIC_SPREAD_SE,
    outlierFraction: CLUSTER_OUTLIER_FRACTION,
    mixedSizeShare: MIXED_SIZE_SHARE,
    splitKurtosis: SPLIT_SITE_KURTOSIS,
    splitSe: SPLIT_SITE_SE
});

const NO_CLUSTERS = ['extended_network', 'no_clusters'];
/** Every `because` each finding can give (closed sets, pinned by the parity fixture). */
export const CLUSTER_BECAUSE = Object.freeze({
    distorted: ['edge_spread', 'axis_family', 'trigonality', 'less_distorted', 'within_null', 'undefined_statistic', ...NO_CLUSTERS],
    lessDistorted: ['below_null', 'above_null', 'within_null', 'single_edge', 'undefined_statistic', ...NO_CLUSTERS],
    breathing: ['above_null', 'below_null', 'within_null', 'undefined_statistic', ...NO_CLUSTERS],
    rigid: ['less_distorted', 'stiffer_bonds', 'split_sites', 'pattern', 'within_null', 'undefined_statistic', ...NO_CLUSTERS],
    staticDistortion: ['above_sampling_error', 'within_sampling_error', ...NO_CLUSTERS],
    splitMemberSites: ['half_or_more_split', 'fewer_than_half_split', 'none_split', ...NO_CLUSTERS],
    trigonalLean: ['leans_trigonal', 'leans_one_edge', 'within_null', 'not_tetrahedra', 'undefined_statistic', ...NO_CLUSTERS],
    favouredAxes: ['enhanced', 'within_null', 'undefined_statistic', ...NO_CLUSTERS],
    'favouredAxes.families': ['enhanced', 'within_null', 'undefined_statistic'],
    'window.extendedNetwork': ['extended_network', 'finite_clusters', 'no_clusters'],
    'window.mixedSizes': ['above_bar', 'within_bar', ...NO_CLUSTERS],
    'window.bondsOutsideWindow': ['above_bar', 'within_bar', ...NO_CLUSTERS],
    'window.contactsBetweenClusters': ['above_bar', 'within_bar', ...NO_CLUSTERS]
});
export const EVIDENCE_AGAINST = Object.freeze(['shapeNull', 'independentNull', 'samplingError', 'siteKurtosis', 'window']);
export const CHANGE_KINDS = Object.freeze(['ratio', 'difference']);

/**
 * Whether an Evidence record clears its bar, from the record alone: every
 * field the bar names defined and past its threshold the record's way.
 */
export const evidenceClears = (record) => {
    const entries = Object.entries(record?.bar || {});
    if (!entries.length) return false;
    const strict = Boolean(record.strict);
    return entries.every(([field, threshold]) => {
        const value = record[field];
        if (value == null) return false;
        if (record.side === 'above') return strict ? value > threshold : value >= threshold;
        if (record.side === 'below') return strict ? value < -threshold : value <= -threshold;
        return strict ? Math.abs(value) > threshold : Math.abs(value) >= threshold;
    });
};

const finiteOrNull = (value) => (Number.isFinite(value) ? value : null);

// One Evidence record; its change from observed and null by its kind.
const evidence = (statistic, against, observedValue, nullValue, zValue, changeKind, bar, side = null, extra = {}) => {
    const observedNumber = finiteOrNull(observedValue);
    const nullNumber = finiteOrNull(nullValue);
    let change = null;
    if (changeKind != null && observedNumber != null && nullNumber != null) {
        if (changeKind === 'ratio') change = nullNumber !== 0 ? finiteOrNull(observedNumber / nullNumber - 1) : null;
        else change = finiteOrNull(observedNumber - nullNumber);
    }
    const record = {
        statistic, against, ...extra,
        observed: observedNumber, null: nullNumber, z: finiteOrNull(zValue), change, changeKind, bar: { ...bar }
    };
    if (side != null) record.side = side;
    record.clears = evidenceClears(record);
    return record;
};

/**
 * The split rule on one site row (siteEllipsoids, the 'sites' request):
 * observed = the PC1 excess kurtosis, z = it in sampling errors of a
 * Gaussian's, κ/√(24/n) — null unless PC1 is resolved (an unresolved axis's
 * direction, and so its kurtosis, is noise). Clears at z ≤ −SPLIT_SITE_SE and
 * κ ≤ SPLIT_SITE_KURTOSIS. Mirrors split_site_evidence.
 */
export const splitSiteEvidence = (row) => {
    const kappa = finiteOrNull(row?.excessKurtosis?.[0]);
    const count = Number(row?.count) || 0;
    const resolved = Boolean(row?.axisResolved?.[0]);
    const z = kappa != null && resolved && count > 0 ? kappa / Math.sqrt(24 / count) : null;
    return evidence('pc1ExcessKurtosis', 'siteKurtosis', kappa, 0, z, 'difference',
        { z: SPLIT_SITE_SE, change: -SPLIT_SITE_KURTOSIS }, 'below', { site: row.referenceNumber });
};

/** Whether one site row is split (two positions): its splitSiteEvidence clears. The app's one rule (list_sites, the clusters). */
export const isSplitSite = (row) => splitSiteEvidence(row).clears;

const defined = (record) => record != null && record.z != null && record.change != null;

// A Finding; its evidence with the records that clear first, the rest after, each in the given order.
const finding = (answer, because, records = [], extra = {}) => ({
    answer, because, ...extra,
    evidence: [...records.filter((record) => record.clears), ...records.filter((record) => !record.clears)]
});

/** An extended network or no clusters: every shape finding null, `because` says which. */
const noClusterVerdict = (because) => {
    const empty = (extra = {}) => finding(null, because, [], extra);
    return {
        rules: { ...CLUSTER_RULES },
        distorted: empty(),
        lessDistorted: empty(),
        breathing: empty(),
        rigid: empty(),
        staticDistortion: empty(),
        splitMemberSites: empty({ sites: [], memberSites: [] }),
        trigonalLean: empty(),
        favouredAxes: empty({ families: [] }),
        window: {
            extendedNetwork: finding(because === 'extended_network', because),
            mixedSizes: empty(),
            bondsOutsideWindow: empty(),
            contactsBetweenClusters: empty()
        }
    };
};

/**
 * The findings of a clusterDistortion result that has clusters (mirrors
 * cluster_verdict). `memberRows` are the analysed kind's member sites' rows
 * of siteEllipsoids, in reference order. Shape statistics against the shape
 * null, the size statistic (breathing sd) against independent atoms.
 */
export const clusterVerdict = (result, memberRows) => {
    const { observed, kind } = result;
    const nulls = { shapeNull: result.shapeNull, independentNull: { mean: result.nullMean, z: result.z } };
    const nullTest = (path, against, effect, side, changeKind = 'ratio') => {
        let values = [observed, nulls[against].mean, nulls[against].z];
        path.forEach((part) => {
            values = values.map((value) => (value && typeof value === 'object' ? value[part] : null));
        });
        return evidence(path.join('.'), against, ...values, changeKind, { z: CLUSTER_Z, change: effect }, side);
    };
    const widthUp = nullTest(['ladderWidth'], 'shapeNull', CLUSTER_EFFECT, 'above');
    const widthDown = nullTest(['ladderWidth'], 'shapeNull', CLUSTER_EFFECT, 'below');
    const breathUp = nullTest(['breathingSd'], 'independentNull', CLUSTER_EFFECT, 'above');
    const breathDown = nullTest(['breathingSd'], 'independentNull', CLUSTER_EFFECT, 'below');
    const tetrahedra = finiteOrNull(observed.eRms) != null;
    const trig = tetrahedra ? nullTest(['t2Trigonality'], 'shapeNull', TRIGONALITY_EFFECT, null, 'difference') : null;
    // obs / null ≥ AXIS_ENHANCEMENT exactly when obs / null − 1 ≥ AXIS_ENHANCEMENT − 1
    // (the subtraction is exact for ratios in [0.5, 2]).
    const familyTests = [];
    [['longAxis', 'long'], ['shortAxis', 'short']].forEach(([axisKey, axis]) => {
        if (!observed[axisKey] || typeof observed[axisKey] !== 'object') return;
        Object.keys(observed[axisKey]).forEach((family) => {
            familyTests.push({ axis, family, record: nullTest([axisKey, family], 'shapeNull', AXIS_ENHANCEMENT - 1, 'above') });
        });
    });
    const familyRecords = familyTests.map(({ record }) => record);
    const favoured = familyRecords.filter((record) => record.clears);
    const trigonal = trig != null && trig.clears;
    const pattern = favoured.length > 0 || trigonal;

    // Distorted: the edge spread above the shape null, or a pattern whatever the spread.
    const consulted = [widthUp, ...(trig != null ? [trig] : []), ...familyRecords];
    let distorted;
    if (widthUp.clears) distorted = [true, 'edge_spread'];
    else if (favoured.length) distorted = [true, 'axis_family'];
    else if (trigonal) distorted = [true, 'trigonality'];
    else if (!consulted.some(defined)) distorted = [null, 'undefined_statistic'];
    else if (widthDown.clears) distorted = [false, 'less_distorted'];
    else distorted = [false, 'within_null'];

    let less;
    if (kind.edges === 1) less = [null, 'single_edge']; // a dimer: its one edge has no spread
    else if (!defined(widthDown)) less = [null, 'undefined_statistic'];
    else if (widthDown.clears) less = [true, 'below_null'];
    else if (widthUp.clears) less = [false, 'above_null'];
    else less = [false, 'within_null'];

    let breathing;
    if (!defined(breathUp)) breathing = [null, 'undefined_statistic'];
    else if (breathUp.clears) breathing = [true, 'above_null'];
    else if (breathDown.clears) breathing = [false, 'below_null'];
    else breathing = [false, 'within_null'];

    // Split member sites: at least half of the analysed kind's member sites.
    const splitRecords = memberRows.map(splitSiteEvidence);
    const splitSites = splitRecords.filter((record) => record.clears).map((record) => record.site);
    const memberSites = memberRows.map((row) => row.referenceNumber);
    let split;
    if (splitSites.length && 2 * splitSites.length >= memberSites.length) split = [true, 'half_or_more_split'];
    else if (splitSites.length) split = [false, 'fewer_than_half_split'];
    else split = [false, 'none_split'];

    // Rigid: less distorted (and no pattern) or a size stiffer than independent
    // atoms — never with split sites (the null already holds the amplitude).
    const isSplit = split[0] === true;
    let rigid;
    if (!defined(widthDown) && !defined(breathDown)) rigid = [null, 'undefined_statistic'];
    else if (isSplit && (widthDown.clears || breathDown.clears)) rigid = [false, 'split_sites'];
    else if (!isSplit && widthDown.clears && !pattern) rigid = [true, 'less_distorted'];
    else if (!isSplit && breathDown.clears) rigid = [true, 'stiffer_bonds'];
    else if (widthDown.clears && pattern) rigid = [false, 'pattern'];
    else rigid = [false, 'within_null'];

    // The average cluster: its edge spread (null 0) in sampling errors, the
    // largest spread / error over the templates.
    const { ideal } = result;
    const staticRecord = evidence('idealSpread', 'samplingError', ideal.spread, 0, ideal.spreadOverSe, null,
        { z: STATIC_SPREAD_SE }, 'above');
    const staticAnswer = staticRecord.clears ? [true, 'above_sampling_error'] : [false, 'within_sampling_error'];

    let trigonalLean;
    if (trig == null) trigonalLean = [null, 'not_tetrahedra'];
    else if (!defined(trig)) trigonalLean = [null, 'undefined_statistic'];
    else if (trig.clears) trigonalLean = [true, trig.change < 0 ? 'leans_trigonal' : 'leans_one_edge'];
    else trigonalLean = [false, 'within_null'];

    const families = familyTests.map(({ axis, family, record }) => {
        let answer;
        if (record.clears) answer = [true, 'enhanced'];
        else if (defined(record)) answer = [false, 'within_null'];
        else answer = [null, 'undefined_statistic'];
        return { axis, family, answer: answer[0], because: answer[1] };
    });
    let axes;
    if (favoured.length) axes = [true, 'enhanced'];
    else if (familyRecords.some(defined)) axes = [false, 'within_null'];
    else axes = [null, 'undefined_statistic'];

    // The window: other cluster sizes, bonds outside it, contacts between clusters.
    const total = result.kinds.reduce((sum, entry) => sum + entry.clusters, 0);
    const share = total ? result.kinds.slice(1).reduce((sum, entry) => sum + entry.clusters, 0) / total : null;
    const windowFinding = (statistic, value, bar) => {
        const record = evidence(statistic, 'window', value, null, null, null, { observed: bar }, 'above', { strict: true });
        return finding(record.clears, record.clears ? 'above_bar' : 'within_bar', [record]);
    };
    const { outliers } = result;
    return {
        rules: { ...CLUSTER_RULES },
        distorted: finding(...distorted, consulted),
        lessDistorted: finding(...less, [widthDown]),
        breathing: finding(...breathing, [breathUp]),
        rigid: finding(...rigid, [widthDown, breathDown]),
        staticDistortion: finding(...staticAnswer, [staticRecord]),
        splitMemberSites: finding(...split, splitRecords, { sites: splitSites, memberSites }),
        trigonalLean: finding(...trigonalLean, trig != null ? [trig] : []),
        favouredAxes: finding(...axes, familyRecords, { families }),
        window: {
            extendedNetwork: finding(false, 'finite_clusters'),
            mixedSizes: windowFinding('otherKindsShare', share, MIXED_SIZE_SHARE),
            bondsOutsideWindow: windowFinding('brokenEdgeFraction', outliers.brokenEdgeFraction, CLUSTER_OUTLIER_FRACTION),
            contactsBetweenClusters: windowFinding('extraContactsPerAtom', outliers.extraContactsPerAtom, CLUSTER_OUTLIER_FRACTION)
        }
    };
};
