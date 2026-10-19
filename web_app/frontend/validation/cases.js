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
// skills. A dev-only case (`devOnly`, saying why) has no test twin:
// dev-fit-unconverged, the one case with fit files (a χ² log, an F(Q) and its
// G(r)), for the fit tools and the rmc-fit-assessment skill, with a kind of
// its own; and the judgment subset (`judgment: true` and a `situation`, plan
// item 1.6 part ii), whose cases live in judgmentCases.js and build with the
// kinds below (options of existing kinds, and elongated_clusters). Its
// held-out matches, at the end of TEST_CASES, were written by a fresh author
// from judgment-spec.md alone, never from these cases or their results.

import { JUDGMENT_DEV_CASES } from './judgmentCases';

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
 * `header` lines go right after the format line (RMCProfile's run counters,
 * "Number of moves generated: …"); none by default.
 */
export const buildConfiguration = ({ cell, supercell, sites, shift, cellState = () => null, seed, header: extraHeader = [] }) => {
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
        ...extraHeader,
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
 * The partial g(r) of every element pair (elements in site order): g → 1 at
 * large r, B neighbours per A = 4π ρ_B ∫ r² g dr. Each distance is a Gaussian
 * of `kernel` Å (a finite-resolution histogram), minimum image in the
 * orthogonal box. Returns `{ r, pairs: [[i, j]], g: [[…] per pair], elements,
 * number, volume }` (bin centres r = (k + ½)·dr); formatPartials writes it.
 */
export const partialGr = (atoms, box, elements, { rMax = 8, dr = 0.01, kernel = 0.015 } = {}) => {
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
    const r = [];
    const g = pairs.map(() => []);
    for (let k = 0; k < bins; k += 1) {
        const at = (k + 0.5) * dr;
        const shell = 4 * Math.PI * at * at * dr;
        r.push(at);
        pairs.forEach(([i, j], p) => {
            const ordered = i === j ? 2 * counts[p][k] : counts[p][k];
            g[p].push((volume * ordered) / (number[i] * number[j] * shell));
        });
    }
    return { r, pairs, g, elements, number, volume };
};

/** A partialGr as RMCProfile's *_PDFpartials.csv (labels "A-B"). */
export const formatPartials = ({ r, pairs, g, elements }) => {
    const header = `   r (Ang),${pairs.map(([i, j]) => `     ${elements[i]}-${elements[j]}`).join(',')}`;
    const rows = r.map((at, k) => `   ${at.toFixed(5)} ,${g.map((values) => `         ${values[k].toFixed(15)} `).join(',')}, `);
    return `${[header, ...rows].join('\n')}\n`;
};

/** The partial g(r) of every element pair, as RMCProfile's *_PDFpartials.csv (partialGr, written). */
export const partialsCsv = (atoms, box, elements, options = {}) => formatPartials(partialGr(atoms, box, elements, options));

// --- fit outputs (RMCProfile's X-ray file family) --------------------------------
//
// What an RMCProfile run that fits one X-ray measurement writes beside its
// configuration, in the formats of the committed Demo run
// (web_app/frontend/public/demo): <stem>_FQ1.csv (Q, F(Q)_RMC, F(Q)_Expt),
// <stem>_FT_XFQ1.csv (its Fourier transform: r(A), X_ray-calc,
// X_ray_exp_renorm) and <stem>-00.log (the χ² history: clock time, moves
// accepted and generated, one χ² column per fit term, the last one the X-ray
// real-space term). The "calculated" columns are the configuration's own
// functions; the "experimental" ones are made from them with a planted
// difference, so the truth about the fit is known.

// The X-ray G(r) of a partialGr, Faber–Ziman normalised as RMCProfile's xPDF
// (f ≈ Z at Q → 0, constant): G = Σ_{i≤j} (2 − δ_ij) c_i c_j f_i f_j (g_ij − 1) / (Σ c f)²,
// so G → −1 below the closest approach and → 0 at large r.
const xrayGr = ({ g, pairs, elements, number }, f) => {
    const total = number.reduce((sum, n) => sum + n, 0);
    const c = number.map((n) => n / total);
    const mean = elements.reduce((sum, el, i) => sum + c[i] * f[el], 0);
    const weights = pairs.map(([i, j]) => ((i === j ? 1 : 2) * c[i] * c[j] * f[elements[i]] * f[elements[j]]) / (mean * mean));
    return g[0].map((_, k) => pairs.reduce((sum, _pair, p) => sum + weights[p] * (g[p][k] - 1), 0));
};

// F(Q) = ρ0 ∫ 4π r² G(r) M(r) sin(Qr)/(Qr) dr over the partials' range, M the
// Lorch window (no truncation ripples): the reciprocal-space function whose
// transform RMCProfile writes as <stem>_FT_XFQ1.csv.
const xrayFq = (r, G, rho0, q) => {
    const dr = r[1] - r[0];
    const reach = r[r.length - 1] + dr / 2;
    const kernel = r.map((at, k) => {
        const x = (Math.PI * at) / reach;
        return 4 * Math.PI * at * at * G[k] * (Math.sin(x) / x) * dr;
    });
    return q.map((Q) => rho0 * r.reduce((sum, at, k) => sum + (kernel[k] * Math.sin(Q * at)) / (Q * at), 0));
};

// RMCProfile's .log number format, Fortran E10.3: 0.161E+00.
const fortranE = (value) => {
    if (value === 0) return '0.000E+00';
    let exponent = Math.floor(Math.log10(Math.abs(value))) + 1;
    let mantissa = Math.round((Math.abs(value) / 10 ** exponent) * 1000);
    if (mantissa >= 1000) {
        mantissa = 100;
        exponent += 1;
    }
    const sign = exponent < 0 ? '-' : '+';
    return `${value < 0 ? '-' : ''}0.${String(mantissa).padStart(3, '0')}E${sign}${String(Math.abs(exponent)).padStart(2, '0')}`;
};

// A χ² that is still falling: a fast early drop, then a slow tail that has
// not levelled off (floor + fast·e^(−t/τfast) + slow·e^(−t/τslow), t the log row).
const chiLaw = ({ floor, fast, tauFast, slow, tauSlow }) => (t) => floor + fast * Math.exp(-t / tauFast) + slow * Math.exp(-t / tauSlow);

// The run's move counters, one entry per log row: moves generated at a steady
// rate, accepted ever more rarely as χ² falls; `total` is what the .rmc6f
// header says at the end of the run (no randomness: the header does not
// depend on the partial g(r)).
const runMoves = (fit) => {
    let accepted = 0;
    let generated = 0;
    const rows = [];
    for (let t = 0; t < fit.rows; t += 1) {
        generated += fit.movesPerRow;
        accepted += Math.round(fit.movesPerRow * (fit.accept.floor + fit.accept.start * Math.exp(-t / fit.accept.tau)));
        rows.push({ accepted, generated });
    }
    return { rows, total: { generated, tried: Math.round(generated * fit.triedShare), accepted, seconds: fit.rows * fit.minutes * 60 } };
};

/**
 * The run's fit files from its partial g(r) and the plan `fit` (the case's
 * params): the X-ray G(r) and F(Q) of the configuration, the "measured" ones
 * with the planted missing feature and noise, and a χ² log whose last term is
 * still falling. Returns `{ fq, gr, log }` (file texts).
 */
const xrayFitFiles = (partial, { f, fit }) => {
    const rng = mulberry32(fit.seed);
    const calcG = xrayGr(partial, f);
    const { centre, height, sd } = fit.missed;
    const missed = partial.r.map((at) => height * Math.exp(-0.5 * ((at - centre) / sd) ** 2));
    const rho0 = partial.number.reduce((sum, n) => sum + n, 0) / partial.volume;

    // G(r) from rMin, as RMCProfile writes it above the closest approach.
    const grRows = partial.r.map((at, k) => [at, calcG[k], calcG[k] + missed[k] + gaussian(rng) * fit.grNoise])
        .filter(([at]) => at >= fit.rMin);
    const gr = ['r(A) ,   X_ray-calc  ,   X_ray_exp_renorm',
        ...grRows.map((row) => row.map((value) => `   ${value.toFixed(14)}       `).join(','))].join('\n');

    const q = Array.from({ length: Math.round((fit.qMax - fit.qMin) / fit.dq) + 1 }, (_, k) => fit.qMin + k * fit.dq);
    const calcF = xrayFq(partial.r, calcG, rho0, q);
    const dataF = xrayFq(partial.r, calcG.map((value, k) => value + missed[k]), rho0, q).map((value) => value + gaussian(rng) * fit.fqNoise);
    const fq = ['       Q,      F(Q)_RMC,  F(Q)_Expt',
        ...q.map((Q, k) => [Q, calcF[k], dataF[k]].map((value) => `${value.toFixed(7).padStart(15)} `).join(',').trimEnd())].join('\n');

    // The log: one row per `minutes` of clock time (HHMMSS.sss, as
    // RMCProfile writes it), the cumulative moves, one χ² column per fit term.
    const laws = fit.chi.map(chiLaw);
    const rows = runMoves(fit).rows.map(({ accepted, generated }, t) => {
        const seconds = (fit.startS + t * fit.minutes * 60) % 86400;
        const clock = [Math.floor(seconds / 3600), Math.floor((seconds % 3600) / 60), seconds % 60].map((v) => String(v).padStart(2, '0')).join('');
        const chi = laws.map((law) => fortranE(law(t) * Math.exp(gaussian(rng) * fit.chiNoise)));
        return `${clock}.000  ${String(accepted).padStart(10)}  ${String(generated).padStart(10)}  ${chi.join('  ')}`;
    });
    const log = [
        `Time       moves_acc moves_gen  ${fit.chiColumns.join('  ')}`,
        `h/m/s/.th    WEIGHT PARAMETERS  ${fit.chiColumns.map(() => '0.100E+01').join('  ')}`,
        ...rows
    ].join('\n');
    return { fq: `${fq}\n`, gr: `${gr}\n`, log: `${log}\n` };
};

// The .rmc6f header's run counters (RMCProfile writes them after the format line).
const movesHeader = (moves) => [
    `Number of moves generated:           ${moves.generated}`,
    `Number of moves tried:               ${moves.tried}`,
    `Number of moves accepted:            ${moves.accepted}`,
    `Accumulated time (s) in running loop: ${moves.seconds.toFixed(2)}`
];

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

// The weak zone-boundary order of the judgment subset (KINDS.zone_boundary
// with `weak: true`, plan item 1.6 ii): ordered antiphase B displacements
// whose amplitude is below half the B thermal rms — scoring.js's 0.03 Å,
// below which a static shift cannot be told from the thermal width — and no
// decoy. The truth names q only: no reading of one site resolves the axis.
const weakZoneBoundary = (p) => {
    if (p.decoy != null) throw new Error('zone_boundary: a weak order has no decoy (decoy: null)');
    if (!(p.delta < THERMAL.B / 2)) throw new Error(`zone_boundary: a weak order is below ${THERMAL.B / 2} Å, got ${p.delta}`);
    const sites = perovskiteSites(p.elements);
    const u = unit(p.axis);
    const config = buildConfiguration({
        cell: [p.a, p.a, p.a],
        supercell: p.supercell,
        sites,
        seed: p.seed,
        shift: (index, _, rng, n) => {
            const d = perovskiteThermal(sites, index, rng);
            if (index !== 1) return d;
            const phase = p.q.reduce((sum, q, k) => sum + q * n[k], 0);
            return add(d, scale(u, p.delta * Math.cos(2 * Math.PI * phase)));
        }
    });
    const B = p.elements[1];
    const q = `(${p.q.join(', ')})`;
    return {
        config,
        truth: {
            summary: `${B} displaced ±${p.delta} Å along [${p.axis.join(' ')}] in an ordered pattern with q = ${q}, far below its thermal `
                + `width (${THERMAL.B} Å rms): a weak superstructure at the edge of significance`,
            features: [{ kind: 'modulation', element: B, accept: ['modulation'], q: p.q }],
            tolerated: [],
            symmetry: { ...perovskiteSymmetry, primary: false },
            traps: [
                `the order is real but weak: analyze_modulations puts ${B} at q = ${q} near its 3σ bar, so a fixed 3σ rule is a coin flip; `
                    + 'an honest answer gives the significance and says how confident it is',
                `${B}'s folded cloud looks thermal: no ellipsoid, kurtosis or direction map shows an amplitude this far below the thermal width, `
                    + 'and nothing resolves its axis'
            ]
        }
    };
};

// --- the case kinds ------------------------------------------------------------------

// The eight <111> directions (±1, ±1, ±1): an axis of this family names them all.
const isBodyDiagonal = (axis) => Array.isArray(axis) && axis.length === 3 && axis.every((value) => Math.abs(value) === 1);

const KINDS = {
    // Two-site disorder: each B atom sits ±delta off its site along `axis`,
    // the sign random per atom. The average site stays put (cubic).
    //
    // `wells: 8` (the judgment subset, plan item 1.6 ii): eight-site disorder
    // instead, each B atom delta off its site along one of the eight <111>
    // directions at random (±delta/√3 on each axis, the three signs
    // independent; `axis` is a <111> direction, the family's name). The eight
    // wells average to an isotropic ellipsoid and orientation tensor, so the
    // split shows only in the kurtosis and the direction map, and the B–X
    // shell splits 3 + 3. Without `wells` (or 2) every random draw is the
    // two-site case's, so those cases build exactly as before.
    offcentre: (p) => {
        const sites = perovskiteSites(p.elements);
        const u = unit(p.axis);
        const eight = p.wells === 8;
        if (p.wells != null && p.wells !== 2 && !eight) throw new Error(`offcentre: wells is 2 or 8, got ${p.wells}`);
        if (eight && !isBodyDiagonal(p.axis)) throw new Error('offcentre: eight wells need a <111> axis');
        const component = p.delta / Math.sqrt(3);
        const config = buildConfiguration({
            cell: [p.a, p.a, p.a],
            supercell: p.supercell,
            sites,
            seed: p.seed,
            shift: (index, _, rng) => {
                const d = perovskiteThermal(sites, index, rng);
                if (index === 1 && eight) return d.map((value) => value + (rng() < 0.5 ? component : -component));
                return index === 1 ? add(d, scale(u, rng() < 0.5 ? p.delta : -p.delta)) : d;
            }
        });
        const [, B, X] = p.elements;
        if (eight) {
            return {
                config,
                truth: {
                    summary: `${B} off-centred ${p.delta} Å along one of the eight <111> directions at random (eight-site disorder); `
                        + 'the average structure stays cubic',
                    features: [{
                        kind: 'off_centre', element: B, accept: ['off_centre'],
                        direction_uvw: p.axis, direction_mode: 'cubic', magnitude_A: p.delta
                    }],
                    tolerated: [],
                    symmetry: { ...perovskiteSymmetry, primary: false },
                    traps: [
                        `the eight wells average to an isotropic ${B} ellipsoid and an isotropic orientation tensor: no principal axis is resolved, `
                            + 'so the ellipsoid shows no split; the Mardia kurtosis is strongly negative, and the direction map, anisotropic at '
                            + 'its default tiling, resolves the eight <111> lobes only on a finer one',
                        `the ${B}–${X} first shell is split 3 + 3 by the off-centring itself (three short, three long bonds), not by a Jahn–Teller distortion`
                    ]
                }
            };
        }
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
    //
    // `weak: true` (the judgment subset, plan item 1.6 ii): the same order with
    // an amplitude far below the B thermal width and no decoy (`decoy: null`):
    // the folded cloud looks thermal, no single-site reading resolves the
    // displacement axis (so the truth names only q), and the modulation
    // engine's calibrated significance sits near its 3σ bar. A case without
    // `weak` keeps its decoy and builds exactly as before.
    zone_boundary: (p) => {
        if (p.weak) return weakZoneBoundary(p);
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

    // The judgment subset (plan item 1.6 ii), the one kind of its own: finite
    // M6 octahedra (a vertex on each cubic axis about the body centre), an
    // anion X at the origin, as the breathing kind places its clusters. Each
    // octahedron is elongated along x, y or z at random — its two apical
    // atoms pushed out by `elongation` Å — so the local octahedron has 4 short
    // (equatorial, `edge`) and 8 long edges, and the average one is regular,
    // its edge between the two. The M–M shell is split across a trough, and
    // it closes at CN 4, an octahedron: the Workbench joins the two peaks and
    // asks when a window cuts them, while the first peak alone holds no edge
    // of the average octahedron. No other kind plants a cluster distortion.
    elongated_clusters: (p) => {
        const axes = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]];
        // Equatorial radius (short edge = r0·√2); the average radius r0 + elongation/3.
        const r0 = p.edge / Math.SQRT2;
        const mean = r0 + p.elongation / 3;
        const sites = [
            { el: p.anion, frac: [0, 0, 0] },
            ...axes.map((v) => ({ el: p.element, frac: v.map((value) => 0.5 + (mean * value) / p.a), axis: v }))
        ];
        const config = buildConfiguration({
            cell: [p.a, p.a, p.a],
            supercell: p.supercell,
            sites,
            seed: p.seed,
            cellState: (_, rng) => Math.floor(rng() * 3),
            shift: (index, longAxis, rng) => {
                const d = gauss3(rng, p.sigma);
                if (index === 0) return d;
                const v = sites[index].axis;
                return add(d, scale(v, (v[longAxis] !== 0 ? r0 + p.elongation : r0) - mean));
            }
        });
        const long = Math.hypot(r0, r0 + p.elongation);
        const M = p.element;
        return {
            config,
            truth: {
                summary: `${M}6 octahedra each elongated along x, y or z at random: locally 4 short + 8 long edges (${p.edge} / ${long.toFixed(2)} Å), `
                    + `the average octahedron regular (edge ${(mean * Math.SQRT2).toFixed(2)} Å)`,
                features: [{ kind: 'cluster_distortion', element: M, accept: ['cluster_distortion'] }],
                // True but partial descriptions of the same distortion: each
                // vertex sits on one of two positions along its own axis, and
                // each M's own M–M shell is split.
                tolerated: [{ kind: 'off_centre', element: M }, { kind: 'jahn_teller', element: M }],
                // A regular octahedron at the body centre, the anion at the origin: m-3m.
                symmetry: { space_group: 'Pm-3m', chiral: false, centrosymmetric: true, polar: false, primary: false },
                traps: [
                    `the ${M}–${M} first shell is split across a trough (the short and the long edges) and closes at 4 ${M} per ${M}, an octahedron: `
                        + 'the first peak alone holds no edge of the average octahedron, so no clusters there; which window is the user\'s choice, '
                        + 'and the Workbench asks when a window cuts the shell',
                    `the ${M} split sites are the elongation seen one atom at a time, not independent off-centring`
                ]
            }
        };
    },

    // Isolated BX6 octahedra (B at the body centre, A at the origin). Each
    // octahedron is Jahn–Teller elongated, 2 long + 4 short bonds, its long
    // axis along x, y or z at random: the average octahedron is regular at
    // (long + 2 short)/3, the local B–X shell is split 4 + 2.
    //
    // `offCentre: { axis, delta }` (the judgment subset, plan item 1.6 ii):
    // each B also sits ±delta off the octahedron's centre along `axis`, the
    // sign random per atom — two features of one polyhedron at once. An
    // off-centred B splits its own shell 1 + 4 + 1 about the mean distance;
    // the Jahn–Teller elongation still puts 4 X below it and 2 above. Without
    // `offCentre` every random draw is as before.
    jahn_teller: (p) => {
        const [A, B, X] = p.elements;
        const mean = (p.long + 2 * p.short) / 3;
        const axes = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]];
        const sites = [
            { el: A, frac: [0, 0, 0] },
            { el: B, frac: [0.5, 0.5, 0.5] },
            ...axes.map((v) => ({ el: X, frac: v.map((value) => 0.5 + (mean * value) / p.a), axis: v }))
        ];
        const off = p.offCentre ? unit(p.offCentre.axis) : null;
        const config = buildConfiguration({
            cell: [p.a, p.a, p.a],
            supercell: p.supercell,
            sites,
            seed: p.seed,
            cellState: (_, rng) => Math.floor(rng() * 3),
            shift: (index, longAxis, rng) => {
                const d = gauss3(rng, p.sigma);
                if (index === 1 && off) return add(d, scale(off, rng() < 0.5 ? p.offCentre.delta : -p.offCentre.delta));
                if (index < 2) return d;
                const v = sites[index].axis;
                const bond = v[longAxis] !== 0 ? p.long : p.short;
                return add(d, scale(v, bond - mean));
            }
        });
        const jt = {
            kind: 'jahn_teller', element: B, accept: ['jahn_teller'],
            magnitude_A: Number((p.long - p.short).toFixed(3)), direction_uvw: [0, 0, 1], direction_mode: 'cubic'
        };
        const summary = `${B}${X}6 octahedra Jahn–Teller elongated: 4 × ${p.short} + 2 × ${p.long} Å (Δ = ${(p.long - p.short).toFixed(2)} Å), `
            + `long axis along x, y or z at random; the average octahedron is regular at ${mean.toFixed(3)} Å`;
        const traps = [
            `the ${B}–${X} first shell has two peaks: a window around the first peak holds about 4 ${X}, the whole shell 6`,
            'the average octahedron is regular: the distortion is local only'
        ];
        return {
            config,
            truth: off ? {
                summary: `${summary}; every ${B} also sits ±${p.offCentre.delta} Å off the centre along [${p.offCentre.axis.join(' ')}] at random (two-site disorder)`,
                features: [jt, {
                    kind: 'off_centre', element: B, accept: ['off_centre'],
                    direction_uvw: p.offCentre.axis, direction_mode: 'axis', magnitude_A: p.offCentre.delta
                }],
                tolerated: [{ kind: 'off_centre', element: X }, { kind: 'cluster_distortion', element: X }],
                symmetry: { ...perovskiteSymmetry, primary: false },
                traps: [
                    ...traps,
                    `two distortions of one octahedron: an off-centred ${B} splits its own shell 1 + 4 + 1 about the mean distance, which cannot put `
                        + `4 ${X} below the mean and 2 above; reading the split shell as the off-centring's alone misses the Jahn–Teller elongation`
                ]
            } : {
                summary,
                features: [jt],
                // True but partial descriptions of the same distortion.
                tolerated: [{ kind: 'off_centre', element: X }, { kind: 'cluster_distortion', element: X }],
                symmetry: { ...perovskiteSymmetry, primary: false },
                traps
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

    // A run that has not finished, with its fit files (dev only). The
    // configuration is rock salt M X (Fm-3m) with isotropic thermal motion
    // only: the structural answer is "no order". The truth is about the FIT:
    // the χ² history still falls (a fast drop, then a slow tail that has not
    // levelled off: the watchdog says improving, though the plot's end looks
    // flat), and the "measured" X-ray G(r) holds a peak the configuration does
    // not reproduce (a Gaussian of `missed.height` at `missed.centre` Å, data
    // above the model); in F(Q), the same measurement, that one peak is a slow
    // oscillation strongest at low Q.
    unconverged_fit: (p) => {
        const [M, X] = p.elements;
        const sites = [...FCC.map((t) => ({ el: M, frac: t })), ...FCC.map((t) => ({ el: X, frac: t.map((value) => (value + 0.5) % 1) }))];
        const moves = runMoves(p.fit).total;
        const config = buildConfiguration({
            cell: [p.a, p.a, p.a],
            supercell: p.supercell,
            sites,
            seed: p.seed,
            header: movesHeader(moves),
            shift: (index, _, rng) => gauss3(rng, p.sigma[sites[index].el])
        });
        const { centre, height, sd } = p.fit.missed;
        const last = p.fit.chiColumns.at(-1);
        return {
            config,
            truth: {
                summary: `rock-salt ${M}${X}, thermal motion only (independent isotropic Gaussian displacements): no local order; Fm-3m. `
                    + `The run has not converged: the ${last} χ² is still falling at the end of its history (a slow tail after a fast drop). `
                    + `The measured X-ray G(r) has a peak at ${centre} Å (height ${height}, sd ${sd} Å, data above the model) that the configuration `
                    + 'does not reproduce; F(Q) is the same measurement in reciprocal space',
                features: [],
                tolerated: [],
                symmetry: { space_group: 'Fm-3m', chiral: false, centrosymmetric: true, polar: false, primary: false },
                // Not scored (the verdict has no fit fields): what the fit
                // tools and the answer checks should find, for the record.
                fit: {
                    convergence: 'improving',
                    chi_column: last,
                    misfit: { dataset: 'xpdf', r_A: centre, height, sd_A: sd, sign: 'data above the model' },
                    same_measurement: ['xray_sq', 'xpdf'],
                    accepted_moves_per_atom: Number((moves.accepted / (sites.length * p.supercell.reduce((n, k) => n * k, 1))).toFixed(2))
                },
                traps: [
                    `the end of the χ² plot looks flat, but ${last} is still falling: the watchdog and assess_convergence say improving, so the run has not converged`,
                    `the G(r) misfit is one feature at ${centre} Å that the configuration lacks, not noise spread over r`,
                    `F(Q) and its Fourier transform G(r) are one measurement: the F(Q) misfit sits at low Q because that one missing peak is a slow `
                        + `oscillation in Q (period 2π/${centre} ≈ ${((2 * Math.PI) / centre).toFixed(1)} Å⁻¹, fading as 1/Q), not a scale or box-size error`,
                    'a misfit is not order in the configuration: every analysis of the configuration shows thermal motion only'
                ]
            },
            fit: (partial) => xrayFitFiles(partial, p)
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
    },
    // Dev only (plan item 4.6): no test case has fit data. Scored on its
    // structure as every case is (no order, Fm-3m); its fit truth is for the
    // fit tools, the answer checks and the record.
    {
        id: 'dev-fit-unconverged', kind: 'unconverged_fit',
        devOnly: 'fit data: exercises assess_convergence, compare_fit and the rmc-fit-assessment skill (plan item 4.6)',
        params: {
            elements: ['Mg', 'O'], a: 4.212, supercell: [6, 6, 6], sigma: { Mg: 0.062, O: 0.066 }, seed: 110,
            // X-ray scattering factors at Q → 0 (f ≈ Z).
            f: { Mg: 12, O: 8 },
            fit: {
                seed: 111,
                // G(r) written from rMin Å; the peak the configuration lacks.
                rMin: 1.0, grNoise: 0.004, missed: { centre: 2.55, height: 0.15, sd: 0.08 },
                qMin: 0.5, qMax: 25, dq: 0.02, fqNoise: 0.002,
                // 400 log rows 10 minutes apart from 08:00; moves.
                rows: 400, minutes: 10, startS: 8 * 3600, movesPerRow: 540, triedShare: 0.97,
                accept: { floor: 0.05, start: 0.25, tau: 40 },
                // χ² per fit term (the last is the history the Workbench charts).
                chiColumns: ['F(Q)_1', 'X_ray_(R)1'],
                chi: [
                    { floor: 0.02, fast: 0.5, tauFast: 30, slow: 0.03, tauSlow: 500 },
                    { floor: 0.01, fast: 1.0, tauFast: 25, slow: 0.02, tauSlow: 600 }
                ],
                chiNoise: 0.002
            }
        }
    },
    // The judgment subset (plan item 1.6 ii), dev only: judgmentCases.js.
    ...JUDGMENT_DEV_CASES
];

/**
 * The judgment subset's situations (plan item 1.6 ii): where a fixed rule
 * that joins the engines' readings is ambiguous or wrong. Each judgment case
 * names one; validation/judgment-spec.md describes them for the author of the
 * held-out matches.
 */
export const JUDGMENT_SITUATIONS = Object.freeze({
    two_features: 'two features at once, where a rule that lets one finding explain a reading misses the other',
    marginal: 'a signal at 2.5–4σ, where a fixed 3σ rule is a coin flip and an honest answer says how confident it is',
    two_readings: 'evidence that two tools read differently, where the right answer needs both readings',
    user_choice: 'a choice only the user can make, where the Workbench asks, the simulated user skips, and the answer must state its default and its consequence'
});

export const TEST_CASES = [
    { id: 'test-offcentre', kind: 'offcentre', params: { elements: ['Sr', 'Ge', 'O'], a: 3.95, supercell: [8, 8, 8], axis: [1, 1, 0], delta: 0.16, seed: 201 } },
    {
        id: 'test-zone-boundary', kind: 'zone_boundary',
        params: { elements: ['Ca', 'Zr', 'O'], a: 4.1, supercell: [8, 8, 8], q: [0.5, 0.5, 0.5], axis: [1, 0, 0], delta: 0.09, decoy: { axis: [0, 0, 1], q: [0, 0.125, 0], amp: 0.12 }, seed: 202 }
    },
    { id: 'test-breathing', kind: 'breathing', params: { shape: 'tetrahedron', element: 'Fe', anion: 'Cl', a: 6.90, edge: 2.75, breathe: 0.035, sigma: 0.05, supercell: [8, 8, 8], seed: 732 } },
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
    },
    // The judgment subset's held-out matches (plan item 1.6 ii), one per dev
    // judgment case, written blind from judgment-spec.md alone (protocol
    // 2026-10-19.11). Each seed is the first of its own sequence that its
    // fixed rule accepts, from the case's engine readings only
    // (docs/algorithms/copilot-validation.md).
    {
        id: 'test-two-features', kind: 'jahn_teller', judgment: true, situation: 'two_features',
        params: {
            elements: ['Rb', 'Cr', 'Br'], a: 10.6, short: 2.53, long: 2.97, sigma: 0.055, supercell: [7, 7, 7],
            offCentre: { axis: [0, 1, 0], delta: 0.11 }, seed: 211
        }
    },
    {
        id: 'test-marginal-above', kind: 'zone_boundary', judgment: true, situation: 'marginal',
        params: { elements: ['Ba', 'Hf', 'O'], a: 4.17, supercell: [10, 10, 10], q: [0.5, 0, 0.5], axis: [1, 1, 0], delta: 0.0128, weak: true, decoy: null, seed: 227 }
    },
    {
        id: 'test-marginal-below', kind: 'zone_boundary', judgment: true, situation: 'marginal',
        params: { elements: ['Ba', 'Hf', 'O'], a: 4.17, supercell: [10, 10, 10], q: [0.5, 0, 0.5], axis: [1, 1, 0], delta: 0.0128, weak: true, decoy: null, seed: 223 }
    },
    {
        id: 'test-two-readings', kind: 'offcentre', judgment: true, situation: 'two_readings',
        params: { elements: ['Cs', 'Ge', 'Cl'], a: 5.43, supercell: [8, 8, 8], axis: [1, 1, 1], delta: 0.23, wells: 8, seed: 231 }
    },
    {
        id: 'test-user-choice', kind: 'elongated_clusters', judgment: true, situation: 'user_choice',
        params: { element: 'Re', anion: 'Te', a: 8.4, edge: 2.62, elongation: 0.36, sigma: 0.035, supercell: [7, 7, 7], seed: 241 }
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
 * A case built: the .rmc6f text, its partial g(r) CSV, its fit files when the
 * case has fit data (dev-fit-unconverged), and the truth. Built lazily (a
 * perovskite case is ~2,500 atoms; its g(r) takes about a second). Its run
 * folder's names are runFolder.js's RUN_NAME and RUN_FILES, the same for
 * every case: no name is built from the case id, which names its kind. A
 * judgment case's truth also says so (`judgment: { situation }`), so its
 * records (which keep the truth; the model never sees it) say which
 * situation they test, and the report reads the subset from them.
 */
export const buildCase = (item, { partials = true } = {}) => {
    const make = KINDS[item.kind];
    if (!make) throw new Error(`unknown case kind '${item.kind}'`);
    if (item.judgment && !JUDGMENT_SITUATIONS[item.situation]) throw new Error(`${item.id}: unknown judgment situation '${item.situation}'`);
    const { config, truth: planted, fit } = make(item.params);
    const truth = item.judgment ? { ...planted, judgment: { situation: item.situation } } : planted;
    const elements = elementsInOrder(config.atoms);
    const partial = partials ? partialGr(config.atoms, config.box, elements) : null;
    return {
        id: item.id,
        split: splitOf(item.id),
        kind: item.kind,
        params: item.params,
        question: CASE_QUESTION,
        truth,
        elements,
        rmc6f: config.text,
        partials: partial ? formatPartials(partial) : null,
        // A case with fit data: its RMCProfile fit files ({ fq, gr, log }),
        // made from the partial g(r), so only when the partials are built.
        fitFiles: partial && fit ? fit(partial) : null
    };
};
