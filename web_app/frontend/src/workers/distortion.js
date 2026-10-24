// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Polyhedral distortion — static-mode port of rmc_toolkits/distortion.py (the
// source of truth; keep the two in sync, the parity fixture
// distortion_fixture.json pins them: tests/generate_distortion_fixture.py).
//
// The ligand atoms within an inclusive window [rMin, rMax] of a centre atom are
// the vertices of its polyhedron, found by the bond-angle engine's neighbour
// search (triplets.js neighborBonds: the linked-cell periodic search with
// explicit image shifts, the same WINDOW_TOL). The polyhedra of the most common
// coordination n among the centres with a ligand (the lowest count on a tie)
// are measured, the rest counted. Per polyhedron:
//
//   DI        = (1/n) Σ |l_i − l̄| / l̄                       Baur 1974, any n
//   ⟨λ⟩       = (1/n) Σ (l_i / l₀)²                          Robinson et al. 1971,
//               l₀ of the regular polyhedron of the same volume   n = 6 or 4
//   σ²_θ      = Σ (θ_i − θ₀)² / (m − 1) over the m cis angles    n = 6 or 4
//
// An octahedron's three trans pairs are the perfect matching (one of 15) with
// the largest angle sum (ties within TRANS_TIE_DEG to the first), its 12 cis
// angles the rest, θ₀ = 90°, V = |d₁·(d₂×d₃)|/6 over its trans diagonals; a
// tetrahedron's 6 angles are all cis, θ₀ = arccos(−1/3), V = |e₁·(e₂×e₃)|/6.
// The same indices are read on the average structure: the basis of
// structureFromRmc6f (one site per reference number at its `mean`, in the unit
// cell), a site counting as its majority species (ties to the code-point first).
// Each statistic is summarised as { count, mean, sd, p10, median, p90 }.

import { coordinationSummaryOf, countPatterns, neighborBonds, WINDOW_TOL } from './triplets.js';

// Work budget of one app request (mirrors APP_MAX_BONDS in distortion.py; the
// fixture pins the two equal): the bonds the window holds, counted before any
// is stored. The engine is unrestricted unless the caller passes `maxBonds`.
export const APP_MAX_BONDS = 1_000_000;

export const OCTAHEDRAL_ANGLE_DEG = 90;
const toDegrees = 180 / Math.PI;
export const TETRAHEDRAL_ANGLE_DEG = Math.acos(-1 / 3) * toDegrees;

// l₀³ per unit volume of the regular polyhedron.
const OCTAHEDRON_RADIUS_CUBED_PER_VOLUME = 3 / 4;
const TETRAHEDRON_RADIUS_CUBED_PER_VOLUME = 9 * Math.sqrt(3) / 8;

export const TRANS_TIE_DEG = 1e-9;
export const QUANTILES = Object.freeze([0.1, 0.5, 0.9]);

// coordination → [shape, ideal cis angle].
const IDEAL_POLYHEDRA = new Map([[6, ['octahedron', OCTAHEDRAL_ANGLE_DEG]], [4, ['tetrahedron', TETRAHEDRAL_ANGLE_DEG]]]);

const pairsOf = (n) => {
    const pairs = [];
    for (let i = 0; i < n; i += 1) for (let j = i + 1; j < n; j += 1) pairs.push([i, j]);
    return pairs;
};

// Vertex pairs (i < j), in itertools.combinations' order.
export const PAIRS = Object.freeze({ 4: pairsOf(4), 6: pairsOf(6) });

// Every perfect matching: the first item pairs with each later one in turn,
// the rest recursively (distortion.py _perfect_matchings, the same order).
const perfectMatchings = (items) => {
    if (!items.length) return [[]];
    const [first, ...rest] = items;
    const matchings = [];
    rest.forEach((partner, k) => {
        for (const tail of perfectMatchings([...rest.slice(0, k), ...rest.slice(k + 1)])) {
            matchings.push([[first, partner], ...tail]);
        }
    });
    return matchings;
};

// The 15 ways to choose an octahedron's trans pairs, as indices into PAIRS[6].
export const OCTAHEDRON_MATCHINGS = Object.freeze(perfectMatchings([0, 1, 2, 3, 4, 5]).map((matching) => (
    matching.map(([i, j]) => PAIRS[6].findIndex(([a, b]) => a === i && b === j))
)));

const capitalize = (symbol) => {
    const text = String(symbol).trim();
    return text.charAt(0).toUpperCase() + text.slice(1).toLowerCase();
};

const byCodePoint = (a, b) => (a < b ? -1 : a > b ? 1 : 0);

/** The species holding most of a site's copies, ties to the code-point first symbol. */
export const siteSpecies = (elementCounts) => Object.entries(elementCounts).reduce((best, entry) => (
    entry[1] > best[1] || (entry[1] === best[1] && entry[0] < best[0]) ? entry : best
))[0];

const quantile = (ordered, q) => {
    const position = (ordered.length - 1) * q;
    const lower = Math.floor(position);
    const upper = Math.min(lower + 1, ordered.length - 1);
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower);
};

/** { count, mean, sd, p10, median, p90 } of the finite values; null when there is none. */
export const summarize = (input) => {
    const values = input.filter((value) => Number.isFinite(value));
    if (!values.length) return null;
    const ordered = Float64Array.from(values).sort();
    let sum = 0;
    for (const value of values) sum += value;
    const mean = sum / values.length;
    let squares = 0;
    for (const value of values) squares += (value - mean) * (value - mean);
    return {
        count: values.length,
        mean,
        sd: values.length > 1 ? Math.sqrt(squares / (values.length - 1)) : null,
        p10: quantile(ordered, QUANTILES[0]),
        median: quantile(ordered, QUANTILES[1]),
        p90: quantile(ordered, QUANTILES[2])
    };
};

// a · (b × c), in the Python engine's evaluation order.
const triple = (a, b, c) => {
    const crossX = b[1] * c[2] - b[2] * c[1];
    const crossY = b[2] * c[0] - b[0] * c[2];
    const crossZ = b[0] * c[1] - b[1] * c[0];
    return (a[0] * crossX + a[1] * crossY) + a[2] * crossZ;
};

const minus = (a, b) => [a.vx - b.vx, a.vy - b.vy, a.vz - b.vz];

/**
 * DI, ⟨λ⟩, σ²_θ and the mean trans angle of one polyhedron: `bonds` its
 * centre-to-ligand bonds ({ vx, vy, vz, length }, triplets.js neighborBonds).
 * ⟨λ⟩ (NaN for a flat polyhedron) and σ²_θ only for 4 or 6 bonds, the trans
 * angle only for 6 (null otherwise). Port of distortion.py polyhedron_indices.
 */
export const polyhedronIndices = (bonds) => {
    const n = bonds.length;
    let sum = 0;
    for (const bond of bonds) sum += bond.length;
    const meanLength = sum / n;
    let spread = 0;
    for (const bond of bonds) spread += Math.abs(bond.length - meanLength);
    const result = { distortionIndex: spread / n / meanLength, quadraticElongation: null, angleVariance: null, transAngle: null };
    const ideal = IDEAL_POLYHEDRA.get(n);
    if (!ideal) return result;
    const pairs = PAIRS[n];
    // The bond-angle engine's evaluation order (bitwise the Python engine's).
    const angles = pairs.map(([i, j]) => {
        const one = bonds[i];
        const two = bonds[j];
        const cosine = (one.vx * two.vx + one.vy * two.vy + one.vz * two.vz) / (one.length * two.length);
        return Math.acos(Math.min(1, Math.max(-1, cosine))) * toDegrees;
    });
    const cis = angles.map(() => true);
    let radiusCubed;
    if (n === 6) {
        const sums = OCTAHEDRON_MATCHINGS.map(([a, b, c]) => (angles[a] + angles[b]) + angles[c]);
        const best = Math.max(...sums);
        const trans = OCTAHEDRON_MATCHINGS[sums.findIndex((value) => value >= best - TRANS_TIE_DEG)];
        trans.forEach((k) => { cis[k] = false; });
        result.transAngle = ((angles[trans[0]] + angles[trans[1]]) + angles[trans[2]]) / 3;
        const [d1, d2, d3] = trans.map((k) => minus(bonds[pairs[k][0]], bonds[pairs[k][1]]));
        radiusCubed = OCTAHEDRON_RADIUS_CUBED_PER_VOLUME * (Math.abs(triple(d1, d2, d3)) / 6);
    } else {
        const [e1, e2, e3] = [1, 2, 3].map((k) => minus(bonds[k], bonds[0]));
        radiusCubed = TETRAHEDRON_RADIUS_CUBED_PER_VOLUME * (Math.abs(triple(e1, e2, e3)) / 6);
    }
    let squares = 0;
    let cisCount = 0;
    angles.forEach((angle, k) => {
        if (!cis[k]) return;
        const deviation = angle - ideal[1];
        squares += deviation * deviation;
        cisCount += 1;
    });
    result.angleVariance = squares / (cisCount - 1);
    const radius = Math.cbrt(radiusCubed);
    let elongation = 0;
    for (const bond of bonds) {
        const ratio = bond.length / radius;
        elongation += ratio * ratio;
    }
    elongation /= n;
    result.quadraticElongation = radius > 0 && Number.isFinite(elongation) ? elongation : NaN;
    return result;
};

const INDEX_KEYS = ['distortionIndex', 'quadraticElongation', 'angleVariance', 'transAngle'];

// The summaries of polyhedra of one coordination (`polyhedra`: each its bonds).
const block = (polyhedra) => {
    if (!polyhedra.length) {
        return { polyhedra: 0, bondLength: null, distortionIndex: null, quadraticElongation: null, angleVariance: null, transAngle: null };
    }
    const indices = polyhedra.map(polyhedronIndices);
    const out = { polyhedra: polyhedra.length, bondLength: summarize(polyhedra.flatMap((bonds) => bonds.map((bond) => bond.length))) };
    INDEX_KEYS.forEach((key) => {
        out[key] = indices[0][key] === null ? null : summarize(indices.map((entry) => entry[key]));
    });
    return out;
};

const finiteOrNull = (value) => (value !== null && Number.isFinite(value) ? value : null);

const fold = (row) => row.map((value) => value - Math.floor(value));

// The average structure's polyhedra (distortion.py _average_structure).
const averageStructure = ({ basis, latticeVectors, supercell }, centre, ligand, window, n) => {
    const empty = { ...block([]), reason: null, sites: [] };
    if (!Array.isArray(basis) || !basis.length) return { ...empty, reason: 'no_basis' };
    const labels = basis.map((row) => siteSpecies(row.elementCounts));
    const centreSites = [];
    const ligandSites = [];
    labels.forEach((label, index) => {
        if (label === centre) centreSites.push(index);
        if (label === ligand) ligandSites.push(index);
    });
    if (!centreSites.length) return { ...empty, reason: 'no_centre_sites' };
    // The unit cell's vectors, as structureFromRmc6f and structure_basis take them.
    const unit = latticeVectors.map((row, i) => row.map((value) => value / Math.max(supercell[i], 1)));
    const wrapped = basis.map((row) => fold(row.mean));
    const perCentre = neighborBonds(wrapped, centreSites, ligandSites, unit, [window]);
    const sites = centreSites.map((site, local) => {
        const row = basis[site];
        const bonds = perCentre[local];
        const record = {
            referenceNumber: row.referenceNumber,
            element: labels[site],
            mixed: Object.keys(row.elementCounts).length > 1,
            copies: row.count,
            coordination: bonds.length,
            bondLengths: bonds.map((bond) => bond.length).sort((a, b) => a - b),
            distortionIndex: null,
            quadraticElongation: null,
            angleVariance: null,
            transAngle: null
        };
        if (bonds.length) {
            const indices = polyhedronIndices(bonds);
            INDEX_KEYS.forEach((key) => { record[key] = finiteOrNull(indices[key]); });
        }
        return record;
    });
    const matching = perCentre.filter((bonds) => bonds.length === n);
    return { ...block(matching), reason: matching.length ? null : 'coordination_differs', sites };
};

const validateWindow = (rMin, rMax) => {
    const lo = Number(rMin);
    const hi = Number(rMax);
    if (rMin === null || rMax === null || rMin === undefined || rMax === undefined
        || !Number.isFinite(lo) || !Number.isFinite(hi) || lo < 0 || hi <= lo) {
        throw new Error(`window must be 0 <= r_min < r_max, got [${rMin}, ${rMax}]`);
    }
    return [lo, hi];
};

/**
 * Baur's DI, Robinson's ⟨λ⟩ and σ²_θ of the `centre`-centred `ligand`
 * polyhedra in [rMin, rMax] Å, for the configuration and its average
 * structure. `config`: { fractional, elements } of every atom, the supercell's
 * `latticeVectors` and `supercell`, and the average-structure `basis`
 * (structureFromRmc6f's; [] for a coords-only file). `maxBonds` (the worker
 * passes APP_MAX_BONDS) refuses a window holding more bonds, before any is
 * stored. Mirrors rmc_toolkits.distortion.polyhedral_distortion.
 */
export const polyhedralDistortion = (config, { centre, ligand, rMin, rMax, maxBonds = null }) => {
    const [lo, hi] = validateWindow(rMin, rMax);
    const centreSymbol = capitalize(centre);
    const ligandSymbol = capitalize(ligand);
    const { fractional, elements, latticeVectors } = config;
    if (fractional.some((row) => row.some((value) => !Number.isFinite(value)))) {
        throw new Error('fractional coordinates contain non-finite values');
    }
    const symbols = elements.map(capitalize);
    if (symbols.length !== fractional.length) throw new Error(`${symbols.length} elements for ${fractional.length} coordinates`);
    const available = [...new Set(symbols)].sort(byCodePoint).join(', ');
    for (const [role, symbol] of [['centre', centreSymbol], ['ligand', ligandSymbol]]) {
        if (!symbols.includes(symbol)) {
            throw new Error(`unknown element '${symbol}' for the ${role}: the configuration holds ${available}`);
        }
    }
    // Fold every coordinate into [0, 1); the image bookkeeping restores the
    // true relative geometry (as in the bond-angle engine).
    const wrapped = fractional.map(fold);
    const centreRows = [];
    const ligandRows = [];
    symbols.forEach((symbol, row) => {
        if (symbol === centreSymbol) centreRows.push(row);
        if (symbol === ligandSymbol) ligandRows.push(row);
    });
    if (maxBonds != null) {
        const patterns = countPatterns(wrapped, centreRows, ligandRows, latticeVectors, [[lo, hi]]);
        let counted = 0;
        for (let centreIndex = 0; centreIndex < centreRows.length; centreIndex += 1) counted += patterns[centreIndex * 2 + 1];
        if (counted > maxBonds) {
            throw new Error(`${centreSymbol}-${ligandSymbol} within ${lo}-${hi} A holds ${counted.toLocaleString('en-US')} bonds, `
                + `over the limit of ${Number(maxBonds).toLocaleString('en-US')} for one request; narrow the window`);
        }
    }
    const perCentre = neighborBonds(wrapped, centreRows, ligandRows, latticeVectors, [[lo, hi]]);
    const maxCount = perCentre.reduce((acc, bonds) => Math.max(acc, bonds.length), 0);
    // coordination[k] = centre atoms with exactly k ligands in the window.
    const coordination = new Array(maxCount + 1).fill(0);
    for (const bonds of perCentre) coordination[bonds.length] += 1;
    let n = 0;
    for (let count = 1; count < coordination.length; count += 1) {
        if (coordination[count] > (n ? coordination[n] : 0)) n = count;
    }
    if (!n) throw new Error(`no polyhedra: no ${centreSymbol} atom has a ${ligandSymbol} neighbour within [${lo}, ${hi}] A`);
    if (n < 2) {
        throw new Error(`no polyhedra: the most common ${centreSymbol} coordination within [${lo}, ${hi}] A is one ${ligandSymbol} `
            + '(a polyhedron needs at least 2 ligands)');
    }
    const chosen = perCentre.filter((bonds) => bonds.length === n);
    const [shape, ideal] = IDEAL_POLYHEDRA.get(n) ?? [null, null];
    return {
        centre: centreSymbol,
        ligand: ligandSymbol,
        window: [lo, hi],
        centres: centreRows.length,
        coordination,
        coordinationSummary: coordinationSummaryOf(coordination),
        polyhedron: {
            coordination: n,
            shape,
            idealAngle: ideal,
            cisAngles: shape === null ? null : PAIRS[n].length - (n === 6 ? 3 : 0)
        },
        skipped: centreRows.length - chosen.length,
        configuration: block(chosen),
        average: averageStructure(config, centreSymbol, ligandSymbol, [lo, hi], n),
        rules: { windowTol: WINDOW_TOL, transTieDeg: TRANS_TIE_DEG, tetrahedralAngle: TETRAHEDRAL_ANGLE_DEG, quantiles: [...QUANTILES] }
    };
};
