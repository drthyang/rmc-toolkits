// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Synthetic RMC configurations with planted local structure and a known truth:
// the cases of the Copilot validation (docs/algorithms/copilot-validation.md).
// Each case is a supercell written as an .rmc6f (reference site and cell index
// on every atom line, coordinates wrapped into [0, 1) as RMCProfile writes
// them) plus its partial g(r), so a case loads as a run folder in the
// Workbench and in the harness. Everything is deterministic: a seeded PRNG,
// fixed parameters, no data files.
//
// Every kind has a DEV case (used to improve tools and skills) and a TEST
// case (the reported numbers only): the same physics with other elements,
// directions, magnitudes or structures. Never tune on TEST: the split test
// pins the TEST definitions by hash and keeps their ids out of the tools and
// skills.

// --- deterministic randomness -----------------------------------------------------

export const mulberry32 = (seed) => {
    let state = seed >>> 0;
    return () => {
        state = (state + 0x6D2B79F5) >>> 0;
        let t = state;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
};

const gaussian = (rng) => {
    let u = 0;
    while (u === 0) u = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rng());
};
const gauss3 = (rng, sigma) => [gaussian(rng) * sigma, gaussian(rng) * sigma, gaussian(rng) * sigma];
const add = (a, b) => a.map((value, i) => value + b[i]);
const scale = (a, s) => a.map((value) => value * s);
const unit = (v) => {
    const n = Math.hypot(...v);
    return v.map((value) => value / n);
};

// --- the .rmc6f writer and the partial g(r) ---------------------------------------

/**
 * One configuration. `cell` is the unit cell's edge lengths (Å, orthogonal),
 * `sites` its reference sites `{ el, frac }`, and `shift(siteIndex, cellState,
 * rng)` an atom's displacement from its site in Å; `cellState(cellIndex, rng)`
 * draws what the atoms of one cell share (a cluster's breathing, a JT axis).
 */
export const buildConfiguration = ({ cell, supercell, sites, shift, cellState = () => null, seed }) => {
    const rng = mulberry32(seed);
    const box = cell.map((length, i) => length * supercell[i]);
    const lines = [];
    const atoms = [];
    let count = 0;
    for (let ix = 0; ix < supercell[0]; ix += 1) {
        for (let iy = 0; iy < supercell[1]; iy += 1) {
            for (let iz = 0; iz < supercell[2]; iz += 1) {
                const n = [ix, iy, iz];
                const state = cellState(n, rng);
                sites.forEach((site, index) => {
                    const d = shift(index, state, rng, n);
                    const frac = [0, 1, 2].map((k) => {
                        const value = ((site.frac[k] + n[k] + d[k] / cell[k]) / supercell[k]) % 1;
                        return value < 0 ? value + 1 : value;
                    });
                    count += 1;
                    atoms.push({ el: site.el, cart: frac.map((value, k) => value * box[k]) });
                    lines.push(`${String(count).padStart(6)} ${site.el.padStart(2)} [1] ${frac.map((value) => value.toFixed(12)).join(' ')} `
                        + `${String(index + 1).padStart(4)} ${ix} ${iy} ${iz}`);
                });
            }
        }
    }
    const header = [
        '(Version 6f format configuration file)',
        `Number of atoms: ${count}`,
        `Supercell dimensions: ${supercell.join(' ')}`,
        'Lattice vectors (Ang):',
        `${box[0].toFixed(6)} 0.000000 0.000000`,
        `0.000000 ${box[1].toFixed(6)} 0.000000`,
        `0.000000 0.000000 ${box[2].toFixed(6)}`,
        'Atoms:'
    ];
    return { text: `${[...header, ...lines].join('\n')}\n`, atoms, box };
};

/**
 * The partial g(r) of every element pair (labels "A-B", elements in site
 * order), as RMCProfile's *_PDFpartials.csv: g → 1 at large r, B neighbours
 * per A = 4π ρ_B ∫ r² g dr. Each distance is a Gaussian of `kernel` Å (a
 * finite-resolution histogram), minimum image in the orthogonal box.
 */
export const partialsCsv = (atoms, box, elements, { rMax = 8, dr = 0.01, kernel = 0.015 } = {}) => {
    const reach = Math.min(rMax, Math.min(...box) / 2);
    const bins = Math.floor(reach / dr);
    const index = new Map(elements.map((el, i) => [el, i]));
    const pairs = [];
    elements.forEach((_, i) => elements.forEach((__, j) => { if (j >= i) pairs.push([i, j]); }));
    const pairIndex = (i, j) => pairs.findIndex(([a, b]) => (a === Math.min(i, j) && b === Math.max(i, j)));
    const lookup = elements.map((_, i) => elements.map((__, j) => pairIndex(i, j)));
    const counts = pairs.map(() => new Float64Array(bins));
    const spread = Math.ceil((3 * kernel) / dr);
    const weights = [];
    for (let k = -spread; k <= spread; k += 1) weights.push(Math.exp(-0.5 * ((k * dr) / kernel) ** 2));
    const total = weights.reduce((sum, w) => sum + w, 0);
    const kinds = atoms.map((atom) => index.get(atom.el));
    for (let a = 0; a < atoms.length; a += 1) {
        const pa = atoms[a].cart;
        for (let b = a + 1; b < atoms.length; b += 1) {
            const pb = atoms[b].cart;
            let r2 = 0;
            for (let k = 0; k < 3; k += 1) {
                let d = pb[k] - pa[k];
                d -= box[k] * Math.round(d / box[k]);
                r2 += d * d;
            }
            if (r2 >= reach * reach) continue;
            const centre = Math.sqrt(r2) / dr;
            const bin = Math.floor(centre);
            const target = counts[lookup[kinds[a]][kinds[b]]];
            for (let k = -spread; k <= spread; k += 1) {
                const j = bin + k;
                if (j >= 0 && j < bins) target[j] += weights[k + spread] / total;
            }
        }
    }
    const volume = box[0] * box[1] * box[2];
    const number = elements.map((el) => atoms.filter((atom) => atom.el === el).length);
    const header = `   r (Ang),${pairs.map(([i, j]) => `     ${elements[i]}-${elements[j]}`).join(',')}`;
    const rows = [];
    for (let k = 0; k < bins; k += 1) {
        const r = (k + 0.5) * dr;
        const shell = 4 * Math.PI * r * r * dr;
        const values = pairs.map(([i, j], p) => {
            const ordered = i === j ? 2 * counts[p][k] : counts[p][k];
            return (volume * ordered) / (number[i] * number[j] * shell);
        });
        rows.push(`   ${r.toFixed(5)} ,${values.map((value) => `         ${value.toFixed(15)} `).join(',')}, `);
    }
    return `${[header, ...rows].join('\n')}\n`;
};

// --- structures ----------------------------------------------------------------------

// ABX3 perovskite, A at the origin (a site at fractional 0 in every
// perovskite case), B at the body centre, X on the face centres. Thermal
// motion: A and B isotropic, X flattened across its B–X bond (U⊥ > U∥, as
// in real perovskites): anisotropic but Gaussian, which is not order.
const perovskiteSites = ([A, B, X]) => [
    { el: A, frac: [0, 0, 0] },
    { el: B, frac: [0.5, 0.5, 0.5] },
    { el: X, frac: [0.5, 0.5, 0], bond: 2 },
    { el: X, frac: [0.5, 0, 0.5], bond: 1 },
    { el: X, frac: [0, 0.5, 0.5], bond: 0 }
];
const THERMAL = { A: 0.08, B: 0.06, X_PAR: 0.05, X_PERP: 0.09 };
const perovskiteThermal = (sites, index, rng) => {
    const site = sites[index];
    if (index === 0) return gauss3(rng, THERMAL.A);
    if (index === 1) return gauss3(rng, THERMAL.B);
    const d = gauss3(rng, THERMAL.X_PERP);
    d[site.bond] = gaussian(rng) * THERMAL.X_PAR;
    return d;
};

const perovskiteSymmetry = { space_group: 'Pm-3m', chiral: false, centrosymmetric: true, polar: false };

// --- the case kinds ------------------------------------------------------------------

const KINDS = {
    // Two-site disorder: each B atom sits ±delta off its site along `axis`,
    // the sign random per atom. The average site stays put (cubic).
    offcentre: (p) => {
        const sites = perovskiteSites(p.elements);
        const u = unit(p.axis);
        const config = buildConfiguration({
            cell: [p.a, p.a, p.a],
            supercell: p.supercell,
            sites,
            seed: p.seed,
            shift: (index, _, rng) => {
                const d = perovskiteThermal(sites, index, rng);
                return index === 1 ? add(d, scale(u, rng() < 0.5 ? p.delta : -p.delta)) : d;
            }
        });
        return {
            config,
            truth: {
                summary: `${p.elements[1]} off-centred ±${p.delta} Å along [${p.axis.join(' ')}] at random (two-site disorder); the average structure stays cubic`,
                features: [{
                    kind: 'off_centre', element: p.elements[1], accept: ['off_centre'],
                    direction_uvw: p.axis, direction_mode: 'axis', magnitude_A: p.delta
                }],
                tolerated: [],
                symmetry: { ...perovskiteSymmetry, primary: false },
                traps: [`the ${p.elements[1]} ellipsoid alone cannot tell a split site from large thermal motion: read the kurtosis and the direction map`]
            }
        };
    },

    // Ordered antiphase displacements of B: +delta·u where q·n is even, −delta·u
    // where it is odd. Each B site's folded cloud is split ±delta exactly as
    // two-site disorder is; only the correlation between cells differs. A
    // stronger box-scale wave of A (the decoy) is the strongest mode.
    zone_boundary: (p) => {
        const sites = perovskiteSites(p.elements);
        const u = unit(p.axis);
        const decoyAxis = unit(p.decoy.axis);
        const config = buildConfiguration({
            cell: [p.a, p.a, p.a],
            supercell: p.supercell,
            sites,
            seed: p.seed,
            shift: (index, _, rng, n) => {
                const d = perovskiteThermal(sites, index, rng);
                if (index === 1) {
                    const phase = p.q.reduce((sum, q, k) => sum + q * n[k], 0);
                    return add(d, scale(u, p.delta * Math.cos(2 * Math.PI * phase)));
                }
                if (index === 0) {
                    const phase = p.decoy.q.reduce((sum, q, k) => sum + q * n[k], 0);
                    return add(d, scale(decoyAxis, p.decoy.amp * Math.cos(2 * Math.PI * phase)));
                }
                return d;
            }
        });
        return {
            config,
            truth: {
                summary: `${p.elements[1]} displaced ±${p.delta} Å along [${p.axis.join(' ')}] in an ordered pattern with q = (${p.q.join(', ')}) `
                    + `(a superstructure the folded average cannot show); ${p.elements[0]} carries a stronger box-scale wave at q = (${p.decoy.q.join(', ')}) as a decoy`,
                features: [{
                    kind: 'modulation', element: p.elements[1], accept: ['modulation'],
                    q: p.q, direction_uvw: p.axis, direction_mode: 'axis'
                }],
                tolerated: [],
                symmetry: { ...perovskiteSymmetry, primary: false },
                traps: [
                    `${p.elements[1]}'s folded cloud is split ±${p.delta} Å like two-site disorder: only the cell-to-cell correlation (analyze_modulations) shows it is ordered`,
                    `${p.elements[0]}'s strongest mode is a long wave at the box scale, not a crystallographic superstructure`
                ]
            }
        };
    },

    // Finite M clusters (a triangle in the (111) plane or a tetrahedron)
    // centred on the body centre, an anion X at the origin. Each cluster is
    // scaled about its centre by (1 + ε), ε ~ N(0, breathe): the cluster size
    // fluctuates collectively while its shape stays regular.
    breathing: (p) => {
        const directions = p.shape === 'triangle'
            ? [[2, -1, -1], [-1, 2, -1], [-1, -1, 2]].map(unit)
            : [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]].map(unit);
        // Centre-to-vertex distance for the edge: triangle e/√3, tetrahedron e·√(3/8).
        const radius = p.shape === 'triangle' ? p.edge / Math.sqrt(3) : p.edge * Math.sqrt(3 / 8);
        const sites = [
            { el: p.anion, frac: [0, 0, 0] },
            ...directions.map((d) => ({ el: p.element, frac: d.map((value) => 0.5 + (radius * value) / p.a), dir: d }))
        ];
        const config = buildConfiguration({
            cell: [p.a, p.a, p.a],
            supercell: p.supercell,
            sites,
            seed: p.seed,
            cellState: (_, rng) => gaussian(rng) * p.breathe,
            shift: (index, epsilon, rng) => {
                const d = gauss3(rng, p.sigma);
                return index === 0 ? d : add(d, scale(sites[index].dir, radius * epsilon));
            }
        });
        const name = p.shape === 'triangle' ? `${p.element}3 triangles` : `${p.element}4 tetrahedra`;
        return {
            config,
            truth: {
                summary: `${name} (edge ${p.edge} Å) that breathe: each cluster scaled by 1 + ε, sd(ε) = ${p.breathe} `
                    + `(edge sd ${(p.edge * p.breathe).toFixed(3)} Å from breathing); their shape stays regular`,
                features: [{ kind: 'cluster_breathing', element: p.element, accept: ['cluster_breathing'] }],
                tolerated: [],
                // The average group of these cells is not pinned by construction.
                symmetry: null,
                traps: ['a size fluctuation is not a shape distortion: the edge spread (shape) matches independent atoms, only the mean edge (breathing) does not']
            }
        };
    },

    // Isolated BX6 octahedra (B at the body centre, A at the origin). Each
    // octahedron is Jahn–Teller elongated, 2 long + 4 short bonds, its long
    // axis along x, y or z at random: the average octahedron is regular at
    // (long + 2 short)/3, the local B–X shell is split 4 + 2.
    jahn_teller: (p) => {
        const [A, B, X] = p.elements;
        const mean = (p.long + 2 * p.short) / 3;
        const axes = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]];
        const sites = [
            { el: A, frac: [0, 0, 0] },
            { el: B, frac: [0.5, 0.5, 0.5] },
            ...axes.map((v) => ({ el: X, frac: v.map((value) => 0.5 + (mean * value) / p.a), axis: v }))
        ];
        const config = buildConfiguration({
            cell: [p.a, p.a, p.a],
            supercell: p.supercell,
            sites,
            seed: p.seed,
            cellState: (_, rng) => Math.floor(rng() * 3),
            shift: (index, longAxis, rng) => {
                const d = gauss3(rng, p.sigma);
                if (index < 2) return d;
                const v = sites[index].axis;
                const bond = v[longAxis] !== 0 ? p.long : p.short;
                return add(d, scale(v, bond - mean));
            }
        });
        return {
            config,
            truth: {
                summary: `${B}${X}6 octahedra Jahn–Teller elongated: 4 × ${p.short} + 2 × ${p.long} Å (Δ = ${(p.long - p.short).toFixed(2)} Å), `
                    + `long axis along x, y or z at random; the average octahedron is regular at ${mean.toFixed(3)} Å`,
                features: [{
                    kind: 'jahn_teller', element: B, accept: ['jahn_teller'],
                    magnitude_A: Number((p.long - p.short).toFixed(3)), direction_uvw: [0, 0, 1], direction_mode: 'cubic'
                }],
                // True but partial descriptions of the same distortion.
                tolerated: [{ kind: 'off_centre', element: X }, { kind: 'cluster_distortion', element: X }],
                symmetry: { ...perovskiteSymmetry, primary: false },
                traps: [
                    `the ${B}–${X} first shell has two peaks: a window around the first peak holds about 4 ${X}, the whole shell 6`,
                    'the average octahedron is regular: the distortion is local only'
                ]
            }
        };
    },

    // Rigid M2 dimers with the first site at fractional 0 along `axis`:
    // about half its atoms wrap to ≈ 1 in the file. Both atoms of a dimer
    // share one random shift (sd sigmaCommon) plus a little of their own.
    rigid_dimer: (p) => {
        const first = [0.5, 0.5, 0.5];
        first[p.axis] = 0;
        const second = [...first];
        second[p.axis] = p.bond / p.a;
        const sites = [{ el: p.element, frac: first }, { el: p.element, frac: second }];
        const config = buildConfiguration({
            cell: [p.a, p.a, p.a],
            supercell: p.supercell,
            sites,
            seed: p.seed,
            cellState: (_, rng) => gauss3(rng, p.sigmaCommon),
            shift: (_, common, rng) => add(common, gauss3(rng, p.sigmaOwn))
        });
        const axisName = 'xyz'[p.axis];
        return {
            config,
            truth: {
                summary: `${p.element}2 dimers (bond ${p.bond} Å along ${axisName}) that move as rigid units: both atoms share one shift `
                    + `(${p.sigmaCommon} Å) plus ${p.sigmaOwn} Å of their own; the first site sits at ${axisName} = 0`,
                features: [{ kind: 'rigid_cluster', element: p.element, accept: ['rigid_cluster'] }],
                tolerated: [],
                symmetry: { space_group: 'P4/mmm', chiral: false, centrosymmetric: true, polar: false, primary: false },
                traps: [`the site at ${axisName} = 0: about half its atoms are written near ${axisName} = 1 (the cell wrap)`]
            }
        };
    },

    // Thermal motion only (anisotropic Gaussian X), in a composition whose
    // literature has ordered off-centring: the answer is "no order".
    thermal_null: (p) => {
        const sites = perovskiteSites(p.elements);
        const config = buildConfiguration({
            cell: [p.a, p.a, p.a],
            supercell: p.supercell,
            sites,
            seed: p.seed,
            shift: (index, _, rng) => perovskiteThermal(sites, index, rng)
        });
        return {
            config,
            truth: {
                summary: 'thermal motion only (independent Gaussian displacements; X flattened across its bond): no local order',
                features: [],
                tolerated: [],
                symmetry: { ...perovskiteSymmetry, primary: false },
                traps: [
                    `${p.elements.join('')} is known in the literature for off-centred ${p.elements[1]}: this configuration has none`,
                    'the X ellipsoids are anisotropic but Gaussian: anisotropic thermal motion is not order'
                ]
            }
        };
    },

    // A coherent shift of every B atom by delta along `axis`: the average
    // structure itself loses symmetry (P4mm along [001], R3m along [111]).
    // Smaller than the Detected SG card's 0.2 Å default tolerance.
    symmetry_lowering: (p) => {
        const sites = perovskiteSites(p.elements);
        const u = unit(p.axis);
        const config = buildConfiguration({
            cell: [p.a, p.a, p.a],
            supercell: p.supercell,
            sites,
            seed: p.seed,
            shift: (index, _, rng) => {
                const d = perovskiteThermal(sites, index, rng);
                return index === 1 ? add(d, scale(u, p.delta)) : d;
            }
        });
        return {
            config,
            truth: {
                summary: `every ${p.elements[1]} shifted by ${p.delta} Å along [${p.axis.join(' ')}] (coherent): the average structure is ${p.group}, polar`,
                features: [{
                    kind: 'symmetry_lowering', element: p.elements[1], accept: ['symmetry_lowering', 'off_centre'],
                    direction_uvw: p.axis, direction_mode: 'axis', magnitude_A: p.delta
                }],
                tolerated: [],
                symmetry: { space_group: p.group, chiral: false, centrosymmetric: false, polar: true, primary: true },
                traps: [`at the Detected SG card's default 0.2 Å the shift (${p.delta} Å) is inside the tolerance: the card says Pm-3m`]
            }
        };
    },

    // Thermal motion only, on a structure given by its sites (the symmetry
    // cases: an achiral non-centrosymmetric group, or a chiral one).
    sites_only: (p) => {
        const config = buildConfiguration({
            cell: [p.a, p.a, p.a],
            supercell: p.supercell,
            sites: p.sites,
            seed: p.seed,
            shift: (_, __, rng) => gauss3(rng, p.sigma)
        });
        return {
            config,
            truth: {
                summary: `${p.structure}, thermal motion only: ${p.group}, ${p.facts.chiral ? 'chiral (a Sohncke group)' : 'achiral'}, `
                    + `${p.facts.centrosymmetric ? 'centrosymmetric' : 'non-centrosymmetric'}, ${p.facts.polar ? 'polar' : 'non-polar'}`,
                features: [],
                tolerated: [],
                symmetry: { space_group: p.group, ...p.facts, primary: true },
                traps: p.traps
            }
        };
    }
};

// Zincblende (F-43m): M at the fcc points, X at (¼,¼,¼) + fcc.
const FCC = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]];
const zincblende = (M, X) => [
    ...FCC.map((t) => ({ el: M, frac: t })),
    ...FCC.map((t) => ({ el: X, frac: t.map((value) => value + 0.25) }))
];
// B20 (P2₁3): (x,x,x), (½+x,½−x,−x), (−x,½+x,½−x), (½−x,−x,½+x).
const b20 = (x) => [[x, x, x], [0.5 + x, 0.5 - x, -x], [-x, 0.5 + x, 0.5 - x], [0.5 - x, -x, 0.5 + x]]
    .map((v) => v.map((value) => ((value % 1) + 1) % 1));
// Sulvanite (P-43m): V at the origin, Cu on the edge centres, S at (x,x,x) and its -43m images.
const sulvanite = (x) => [
    { el: 'V', frac: [0, 0, 0] },
    { el: 'Cu', frac: [0.5, 0, 0] }, { el: 'Cu', frac: [0, 0.5, 0] }, { el: 'Cu', frac: [0, 0, 0.5] },
    ...[[x, x, x], [-x, -x, x], [-x, x, -x], [x, -x, -x]].map((v) => ({ el: 'S', frac: v.map((value) => ((value % 1) + 1) % 1) }))
];

const ACHIRAL_TRAP = 'no inversion centre does not make a structure chiral: a group with mirrors or rotoinversions is achiral';
const CHIRAL_TRAP = 'a chiral (Sohncke) group: no mirrors, no inversion, no rotoinversions';

// --- the cases -----------------------------------------------------------------------

/** The question every case is asked (the same for every arm and model). */
export const CASE_QUESTION = 'Analyse the local structure of this RMC configuration. Is there order or distortion beyond '
    + 'the average structure (off-centred atoms, correlated clusters, a superstructure, distorted polyhedra)? '
    + 'What is the average structure\'s space group, and is it chiral, centrosymmetric, polar? Base every claim '
    + 'on the Workbench\'s analyses, and say so plainly if there is no such order.';

export const DEV_CASES = [
    { id: 'dev-offcentre', kind: 'offcentre', params: { elements: ['K', 'Nb', 'O'], a: 4.0, supercell: [8, 8, 8], axis: [0, 0, 1], delta: 0.2, seed: 101 } },
    {
        id: 'dev-zone-boundary', kind: 'zone_boundary',
        params: { elements: ['Sr', 'Ti', 'O'], a: 3.9, supercell: [8, 8, 8], q: [0, 0, 0.5], axis: [0, 0, 1], delta: 0.1, decoy: { axis: [0, 1, 0], q: [0.125, 0, 0], amp: 0.12 }, seed: 102 }
    },
    { id: 'dev-breathing', kind: 'breathing', params: { shape: 'triangle', element: 'Mo', anion: 'O', a: 7.0, edge: 2.8, breathe: 0.04, sigma: 0.05, supercell: [8, 8, 8], seed: 103 } },
    { id: 'dev-jahn-teller', kind: 'jahn_teller', params: { elements: ['K', 'Cu', 'F'], a: 8.0, short: 1.93, long: 2.25, sigma: 0.05, supercell: [6, 6, 6], seed: 104 } },
    { id: 'dev-cell-wrap', kind: 'rigid_dimer', params: { element: 'V', a: 5.0, axis: 0, bond: 2.25, sigmaCommon: 0.08, sigmaOwn: 0.01, supercell: [8, 8, 8], seed: 105 } },
    { id: 'dev-thermal-null', kind: 'thermal_null', params: { elements: ['Ba', 'Ti', 'O'], a: 4.0, supercell: [8, 8, 8], seed: 106 } },
    { id: 'dev-symmetry-lowering', kind: 'symmetry_lowering', params: { elements: ['K', 'Ta', 'O'], a: 4.0, supercell: [8, 8, 8], axis: [0, 0, 1], delta: 0.07, group: 'P4mm', seed: 107 } },
    {
        id: 'dev-chirality-trap', kind: 'sites_only',
        params: {
            structure: 'zincblende CuBr', a: 5.69, supercell: [6, 6, 6], sigma: 0.07, sites: zincblende('Cu', 'Br'), group: 'F-43m',
            facts: { chiral: false, centrosymmetric: false, polar: false }, traps: [ACHIRAL_TRAP], seed: 108
        }
    },
    {
        id: 'dev-chiral', kind: 'sites_only',
        params: {
            structure: 'B20 FeSi', a: 4.49, supercell: [6, 6, 6], sigma: 0.06,
            sites: [...b20(0.1358).map((frac) => ({ el: 'Fe', frac })), ...b20(0.844).map((frac) => ({ el: 'Si', frac }))],
            group: 'P2_13', facts: { chiral: true, centrosymmetric: false, polar: false }, traps: [CHIRAL_TRAP], seed: 109
        }
    }
];

export const TEST_CASES = [
    { id: 'test-offcentre', kind: 'offcentre', params: { elements: ['Sr', 'Ge', 'O'], a: 3.95, supercell: [8, 8, 8], axis: [1, 1, 0], delta: 0.16, seed: 201 } },
    {
        id: 'test-zone-boundary', kind: 'zone_boundary',
        params: { elements: ['Ca', 'Zr', 'O'], a: 4.1, supercell: [8, 8, 8], q: [0.5, 0.5, 0.5], axis: [1, 0, 0], delta: 0.09, decoy: { axis: [0, 0, 1], q: [0, 0.125, 0], amp: 0.12 }, seed: 202 }
    },
    { id: 'test-breathing', kind: 'breathing', params: { shape: 'tetrahedron', element: 'Nb', anion: 'Se', a: 7.4, edge: 2.95, breathe: 0.035, sigma: 0.05, supercell: [8, 8, 8], seed: 203 } },
    { id: 'test-jahn-teller', kind: 'jahn_teller', params: { elements: ['Rb', 'Mn', 'Cl'], a: 9.5, short: 2.48, long: 2.85, sigma: 0.06, supercell: [6, 6, 6], seed: 204 } },
    { id: 'test-cell-wrap', kind: 'rigid_dimer', params: { element: 'Ti', a: 5.2, axis: 2, bond: 2.29, sigmaCommon: 0.07, sigmaOwn: 0.012, supercell: [8, 8, 8], seed: 205 } },
    { id: 'test-thermal-null', kind: 'thermal_null', params: { elements: ['Pb', 'Ti', 'O'], a: 3.97, supercell: [8, 8, 8], seed: 206 } },
    { id: 'test-symmetry-lowering', kind: 'symmetry_lowering', params: { elements: ['Ca', 'Ti', 'O'], a: 3.85, supercell: [8, 8, 8], axis: [1, 1, 1], delta: 0.06, group: 'R3m', seed: 207 } },
    {
        id: 'test-chirality-trap', kind: 'sites_only',
        params: {
            structure: 'sulvanite Cu3VS4', a: 5.39, supercell: [6, 6, 6], sigma: 0.07, sites: sulvanite(0.235), group: 'P-43m',
            facts: { chiral: false, centrosymmetric: false, polar: false }, traps: [ACHIRAL_TRAP], seed: 208
        }
    },
    {
        id: 'test-chiral', kind: 'sites_only',
        params: {
            structure: 'B20 MnSi', a: 4.558, supercell: [6, 6, 6], sigma: 0.06,
            sites: [...b20(0.138).map((frac) => ({ el: 'Mn', frac })), ...b20(0.845).map((frac) => ({ el: 'Si', frac }))],
            group: 'P2_13', facts: { chiral: true, centrosymmetric: false, polar: false }, traps: [CHIRAL_TRAP], seed: 209
        }
    }
];

export const ALL_CASES = [...DEV_CASES, ...TEST_CASES];
export const splitOf = (id) => (TEST_CASES.some((item) => item.id === id) ? 'test' : (DEV_CASES.some((item) => item.id === id) ? 'dev' : null));
export const casesOf = (split) => {
    if (split === 'dev') return DEV_CASES;
    if (split === 'test') return TEST_CASES;
    throw new Error(`unknown split '${split}' (dev or test)`);
};

const elementsInOrder = (sites) => [...new Set(sites.map((site) => site.el))];

/**
 * A case built: the .rmc6f text, its partial g(r) CSV, the truth, and the
 * file names a run folder uses. Built lazily (a perovskite case is ~2,500
 * atoms; its g(r) takes about a second).
 */
export const buildCase = (item, { partials = true } = {}) => {
    const make = KINDS[item.kind];
    if (!make) throw new Error(`unknown case kind '${item.kind}'`);
    const { config, truth } = make(item.params);
    const elements = elementsInOrder(config.atoms);
    return {
        id: item.id,
        split: splitOf(item.id),
        kind: item.kind,
        params: item.params,
        question: CASE_QUESTION,
        truth,
        elements,
        rmc6f: config.text,
        rmc6fName: `${item.id}.rmc6f`,
        partialsName: `${item.id}_PDFpartials.csv`,
        partials: partials ? partialsCsv(config.atoms, config.box, elements) : null
    };
};
