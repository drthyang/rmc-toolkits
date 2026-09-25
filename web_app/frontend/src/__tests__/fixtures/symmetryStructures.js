// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/__tests__/fixtures/symmetryStructures.js
//
// Well-known crystal structures and the operations that re-describe them (another
// cell, another axis order, another origin, RMC-like noise, shuffled site order), for
// the symmetry-finder tests. Every structure is written in its ITA standard setting;
// the helpers produce the variants a real .rmc6f can arrive in.

import { readFileSync } from 'node:fs';

import { structureFromRmc6f } from '../../browserData.js';

export const wrap = (v) => { const y = v - Math.floor(v); return y > 1 - 1e-9 ? 0 : y; };
export const at = (el, x, y, z) => ({ el, frac: [wrap(x), wrap(y), wrap(z)] });

export const cubicCell = (a) => [[a, 0, 0], [0, a, 0], [0, 0, a]];
export const tetragonalCell = (a, c) => [[a, 0, 0], [0, a, 0], [0, 0, c]];
export const hexagonalCell = (a, c) => [[a, 0, 0], [-a / 2, (a * Math.sqrt(3)) / 2, 0], [0, 0, c]];

const FCC = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]];
const BCC = [[0, 0, 0], [0.5, 0.5, 0.5]];

/** Place `motif` ([el, x, y, z] rows) at every centring translation. */
export function centred(translations, motif) {
    const basis = [];
    for (const [cx, cy, cz] of translations) {
        for (const [el, x, y, z] of motif) basis.push(at(el, x + cx, y + cy, z + cz));
    }
    return basis;
}

export const STRUCTURES = {
    rocksalt: () => ({ A: cubicCell(5.64), basis: centred(FCC, [['Na', 0, 0, 0], ['Cl', 0.5, 0, 0]]), symbol: 'Fm-3m', number: 225 }),
    cscl: () => ({ A: cubicCell(4.12), basis: [at('Cs', 0, 0, 0), at('Cl', 0.5, 0.5, 0.5)], symbol: 'Pm-3m', number: 221 }),
    perovskite: () => ({
        A: cubicCell(3.905),
        basis: [at('Sr', 0, 0, 0), at('Ti', 0.5, 0.5, 0.5), at('O', 0.5, 0.5, 0), at('O', 0.5, 0, 0.5), at('O', 0, 0.5, 0.5)],
        symbol: 'Pm-3m', number: 221,
    }),
    diamond: () => ({ A: cubicCell(3.567), basis: centred(FCC, [['C', 0, 0, 0], ['C', 0.25, 0.25, 0.25]]), symbol: 'Fd-3m', number: 227 }),
    zincblende: () => ({ A: cubicCell(5.41), basis: centred(FCC, [['Zn', 0, 0, 0], ['S', 0.25, 0.25, 0.25]]), symbol: 'F-43m', number: 216 }),
    bccIron: () => ({ A: cubicCell(2.87), basis: centred(BCC, [['Fe', 0, 0, 0]]), symbol: 'Im-3m', number: 229 }),
    rutile: () => ({
        A: tetragonalCell(4.594, 2.959),
        basis: [
            at('Ti', 0, 0, 0), at('Ti', 0.5, 0.5, 0.5),
            at('O', 0.305, 0.305, 0), at('O', -0.305, -0.305, 0), at('O', 0.805, 0.195, 0.5), at('O', 0.195, 0.805, 0.5),
        ],
        symbol: 'P4_2/mnm', number: 136,
    }),
    hcp: () => ({ A: hexagonalCell(3.21, 5.21), basis: [at('Mg', 1 / 3, 2 / 3, 0.25), at('Mg', 2 / 3, 1 / 3, 0.75)], symbol: 'P6_3/mmc', number: 194 }),
    wurtzite: () => ({
        A: hexagonalCell(3.25, 5.21),
        basis: [at('Zn', 1 / 3, 2 / 3, 0), at('Zn', 2 / 3, 1 / 3, 0.5), at('O', 1 / 3, 2 / 3, 0.382), at('O', 2 / 3, 1 / 3, 0.882)],
        symbol: 'P6_3mc', number: 186,
    }),
    // Bi-type (A7) in hexagonal axes, obverse R centring.
    bismuth: () => ({
        A: hexagonalCell(4.546, 11.862),
        basis: centred([[0, 0, 0], [2 / 3, 1 / 3, 1 / 3], [1 / 3, 2 / 3, 2 / 3]], [['Bi', 0, 0, 0.2341], ['Bi', 0, 0, -0.2341]]),
        symbol: 'R-3m', number: 166,
    }),
    // GdFeO3-type Pnma perovskite (CaTiO3 values).
    pnmaPerovskite: () => {
        const ops = [
            [(x, y, z) => [x, y, z]], [(x, y, z) => [-x + 0.5, -y, z + 0.5]], [(x, y, z) => [-x, y + 0.5, -z]],
            [(x, y, z) => [x + 0.5, -y + 0.5, -z + 0.5]], [(x, y, z) => [-x, -y, -z]], [(x, y, z) => [x + 0.5, y, -z + 0.5]],
            [(x, y, z) => [x, -y + 0.5, z]], [(x, y, z) => [-x + 0.5, y + 0.5, z + 0.5]],
        ].map((o) => o[0]);
        return {
            A: [[5.44, 0, 0], [0, 7.64, 0], [0, 0, 5.38]],
            basis: orbits(ops, [['Ca', 0.0065, 0.25, -0.0135], ['Ti', 0, 0, 0.5], ['O', 0.571, 0.25, 0.0086], ['O', 0.2891, 0.0373, 0.7109]]),
            symbol: 'Pnma', number: 62,
        };
    },
};

/** Expand asymmetric-unit rows through coordinate-triplet functions, dropping duplicates. */
export function orbits(ops, rows) {
    const basis = [];
    for (const [el, x, y, z] of rows) {
        for (const op of ops) {
            const q = op(x, y, z).map(wrap);
            const dup = basis.some((s) => s.el === el && s.frac.every((v, i) => Math.abs((((v - q[i]) % 1) + 1.5) % 1 - 0.5) < 1e-6));
            if (!dup) basis.push({ el, frac: q });
        }
    }
    return basis;
}

const matMul = (A, B) => A.map((r) => [0, 1, 2].map((j) => r[0] * B[0][j] + r[1] * B[1][j] + r[2] * B[2][j]));
const det3 = (m) => m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
export function inv3(m) {
    const d = det3(m);
    const c = (i, j) => {
        const r = [0, 1, 2].filter((x) => x !== i);
        const s = [0, 1, 2].filter((x) => x !== j);
        return ((i + j) % 2 ? -1 : 1) * (m[r[0]][s[0]] * m[r[1]][s[1]] - m[r[0]][s[1]] * m[r[1]][s[0]]);
    };
    return [0, 1, 2].map((i) => [0, 1, 2].map((j) => c(j, i) / d));
}

/**
 * The same crystal described on a new cell whose basis vectors are the ROWS of M in the
 * old fractional basis (a'_i = Σ_j M_ij a_j). M may be a permutation, a supercell (det > 1)
 * or a primitive/other cell (det < 1); an origin shift (old fractional) is applied first.
 */
export function redescribe({ A, basis }, M, originShift = [0, 0, 0]) {
    const Ap = matMul(M, A);
    const Minv = inv3(M);
    const out = [];
    const span = 3;
    for (const s of basis) {
        for (let i = -span; i <= span; i++) for (let j = -span; j <= span; j++) for (let k = -span; k <= span; k++) {
            const x = [s.frac[0] + i - originShift[0], s.frac[1] + j - originShift[1], s.frac[2] + k - originShift[2]];
            // x = Mᵀ x'  ⇒  x' = x·M⁻¹ (row vector)
            const xp = [0, 1, 2].map((c) => x[0] * Minv[0][c] + x[1] * Minv[1][c] + x[2] * Minv[2][c]);
            if (xp.some((v) => v < -1e-7 || v >= 1 - 1e-7)) continue;
            const w = xp.map(wrap);
            const dup = out.some((o) => o.el === s.el && o.frac.every((v, c) => Math.abs((((v - w[c]) % 1) + 1.5) % 1 - 0.5) < 1e-6));
            if (!dup) out.push({ el: s.el, frac: w });
        }
    }
    return { A: Ap, basis: out };
}

// Proper (det +1) axis relabellings: new (a, b, c) = rows in the old basis.
export const AXIS_SETTINGS = {
    abc: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    bca: [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
    cab: [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
    'ba-c': [[0, 1, 0], [1, 0, 0], [0, 0, -1]],
    'a-cb': [[1, 0, 0], [0, 0, -1], [0, 1, 0]],
    '-cba': [[0, 0, -1], [0, 1, 0], [1, 0, 0]],
};

/** Seeded Gaussian noise (Å, isotropic) on every site, as an RMC average structure has. */
export function withNoise({ A, basis }, sigma, seed = 1) {
    let s = seed >>> 0;
    const rand = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
    const gauss = () => { let u = 0; while (u === 0) u = rand(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rand()); };
    const Ainv = inv3(A);
    return {
        A,
        basis: basis.map((site) => {
            const d = [gauss() * sigma, gauss() * sigma, gauss() * sigma];
            const df = [0, 1, 2].map((j) => d[0] * Ainv[0][j] + d[1] * Ainv[1][j] + d[2] * Ainv[2][j]);
            return { ...site, frac: site.frac.map((v, i) => wrap(v + df[i])) };
        }),
    };
}

/** Deterministic shuffle of the site order. */
export function shuffled(basis, seed = 1) {
    let s = seed >>> 0;
    const rand = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
    const out = basis.slice();
    for (let i = out.length - 1; i > 0; i--) { const j = Math.floor(rand() * (i + 1)); [out[i], out[j]] = [out[j], out[i]]; }
    return out;
}

/**
 * Products of an operation set that are missing from it (R equal, t equal modulo the
 * lattice within `ttol` in every fractional component). Empty for a group.
 */
export function closureDefects(ops, ttol = 0.02) {
    const has = (R, t) => ops.some((o) => o.R.every((r, i) => r.every((v, j) => v === R[i][j]))
        && o.t.every((v, i) => Math.abs((((v - t[i]) % 1) + 1.5) % 1 - 0.5) < ttol));
    const missing = [];
    for (const a of ops) {
        for (const b of ops) {
            const R = a.R.map((row) => [0, 1, 2].map((j) => row[0] * b.R[0][j] + row[1] * b.R[1][j] + row[2] * b.R[2][j]));
            const t = [0, 1, 2].map((i) => wrap(a.R[i][0] * b.t[0] + a.R[i][1] * b.t[1] + a.R[i][2] * b.t[2] + a.t[i]));
            if (!has(R, t)) missing.push({ R, t });
        }
    }
    return missing;
}

/** The bundled demo run (GaTa4Se8, F-43m, 52-site average basis). */
export function demoStructure() {
    const text = readFileSync(new URL('../../../public/demo/GTS_250K.rmc6f', import.meta.url), 'utf8');
    return structureFromRmc6f({ text });
}
