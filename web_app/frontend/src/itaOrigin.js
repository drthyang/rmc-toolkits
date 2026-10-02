// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/itaOrigin.js
//
// The origin shift that turns a space group's exact operations, in the standard cell it was
// named in, into International Tables' own (itaOperations.js): the CIF export's last step
// (averageStructure.js), so that the H–M symbol alone describes the written structure and
// every orbit can be given its Wyckoff letter.
//
// Shifting coordinates by p (x' = x + p) changes each translation to t + (I − R)·p. The
// shifted group equals ITA's when each ITA generator {R_g | τ_g} is in it, i.e. when
//   (I − R_g)·p ≡ τ_g − t_g + c_g   (mod ℤ³)      for some centring vector c_g,
// with t_g our translation for the rotation R_g (both groups share the centring, checked
// first). The stacked integer matrix M = [I − R_g] is brought to diagonal form U·M·V = D by
// unimodular row and column operations, so with q = V⁻¹·p the congruences decouple:
// D_ii·q_i ≡ (U·rhs)_i, solved by q_i = ((U·rhs)_i + m)/D_ii for m = 0 … D_ii − 1, a free
// q_i where D_ii = 0 (a polar direction), and no solution when a zero row of D meets a
// non-integer (U·rhs)_i. That lists every solution — the ITA origins equivalent under the
// group's normalizer included — and each candidate is checked against the full group.
// Of the valid shifts, the one that moves the structure least (Cartesian, counted from the
// .rmc6f origin) is returned; nothing moves along a polar axis.

import { itaGenerators, itaOperations } from './itaOperations.js';
import { parseCoordinateForm } from './wyckoff.js';

const cyc = (x) => x - Math.round(x);
const rotKey = (R) => R.flat().join(',');
const mulV = (M, v) => M.map((row) => row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
const sameMod1 = (u, v) => Math.abs(cyc(u[0] - v[0])) < 1e-6 && Math.abs(cyc(u[1] - v[1])) < 1e-6 && Math.abs(cyc(u[2] - v[2])) < 1e-6;

/**
 * Diagonal form of an integer m×n matrix by unimodular row (U, m×m) and column (V, n×n)
 * operations: U·M·V = D with D zero off the diagonal and D_ii ≥ 0. (The divisibility chain
 * of the Smith form is not needed to solve congruences.)
 */
export function diagonalForm(M) {
  const m = M.length;
  const n = M[0].length;
  const D = M.map((row) => row.slice());
  const U = Array.from({ length: m }, (_, i) => Array.from({ length: m }, (_, j) => (i === j ? 1 : 0)));
  const V = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
  const swapRows = (X, a, b) => { [X[a], X[b]] = [X[b], X[a]]; };
  const swapCols = (X, a, b) => { for (const row of X) [row[a], row[b]] = [row[b], row[a]]; };
  const addRow = (X, to, from, k) => { for (let j = 0; j < X[to].length; j++) X[to][j] += k * X[from][j]; };
  const addCol = (X, to, from, k) => { for (const row of X) row[to] += k * row[from]; };
  let rank = 0;
  for (let t = 0; t < Math.min(m, n); t++) {
    for (;;) {
      let pi = -1;
      let pj = -1;
      for (let i = t; i < m; i++) {
        for (let j = t; j < n; j++) {
          if (D[i][j] && (pi < 0 || Math.abs(D[i][j]) < Math.abs(D[pi][pj]))) { pi = i; pj = j; }
        }
      }
      if (pi < 0) return { D, U, V, rank };
      swapRows(D, t, pi); swapRows(U, t, pi);
      swapCols(D, t, pj); swapCols(V, t, pj);
      let clean = true;
      for (let i = t + 1; i < m; i++) {
        const k = Math.round(D[i][t] / D[t][t]);
        if (k) { addRow(D, i, t, -k); addRow(U, i, t, -k); }
        if (D[i][t]) clean = false;
      }
      for (let j = t + 1; j < n; j++) {
        const k = Math.round(D[t][j] / D[t][t]);
        if (k) { addCol(D, j, t, -k); addCol(V, j, t, -k); }
        if (D[t][j]) clean = false;
      }
      if (clean) break;
    }
    if (D[t][t] < 0) { D[t] = D[t].map((v) => (v ? -v : 0)); U[t] = U[t].map((v) => (v ? -v : 0)); }
    rank = t + 1;
  }
  return { D, U, V, rank };
}

// Cartesian length of a fractional vector in the cell with rows A.
const cartLength = (A, d) => Math.hypot(
  d[0] * A[0][0] + d[1] * A[1][0] + d[2] * A[2][0],
  d[0] * A[0][1] + d[1] * A[1][1] + d[2] * A[2][1],
  d[0] * A[0][2] + d[1] * A[1][2] + d[2] * A[2][2],
);

// The shortest representative of v modulo the lattice and centring, then with its component
// along the directions every rotation fixes (a polar axis) removed in the cell's metric.
function shortest(v, centring, fixedBasis, A) {
  let best = null;
  for (const c of centring) {
    const w = v.map((x, i) => cyc(x - c[i]));
    const length = cartLength(A, w);
    if (!best || length < best.length - 1e-12) best = { w, length };
  }
  let w = best.w;
  if (fixedBasis.length) {
    // w − B·(BᵀGB)⁻¹·BᵀG·w, G = A·Aᵀ: the G-orthogonal projection off the fixed subspace.
    const G = [0, 1, 2].map((i) => [0, 1, 2].map((j) => A[i][0] * A[j][0] + A[i][1] * A[j][1] + A[i][2] * A[j][2]));
    const g = (u, x) => u.reduce((s, ui, i) => s + ui * (G[i][0] * x[0] + G[i][1] * x[1] + G[i][2] * x[2]), 0);
    const k = fixedBasis.length;
    const BtGB = fixedBasis.map((u) => fixedBasis.map((x) => g(u, x)));
    const rhs = fixedBasis.map((u) => g(u, w));
    const coef = k === 1 ? [rhs[0] / BtGB[0][0]] : solveSmall(BtGB, rhs);
    w = w.map((x, i) => x - fixedBasis.reduce((s, u, j) => s + coef[j] * u[i], 0));
  }
  return w;
}

// k×k linear solve (k ≤ 3) by Gaussian elimination.
function solveSmall(M, b) {
  const k = b.length;
  const a = M.map((row, i) => [...row, b[i]]);
  for (let c = 0; c < k; c++) {
    let p = c;
    for (let r = c + 1; r < k; r++) if (Math.abs(a[r][c]) > Math.abs(a[p][c])) p = r;
    [a[c], a[p]] = [a[p], a[c]];
    for (let r = 0; r < k; r++) {
      if (r === c) continue;
      const f = a[r][c] / a[c][c];
      for (let j = c; j <= k; j++) a[r][j] -= f * a[c][j];
    }
  }
  return a.map((row, i) => row[k] / row[i]);
}

// Independent directions fixed by every rotation: columns of the Reynolds projector
// P = (1/|P|)·ΣR (a polar axis or plane), as a basis.
function fixedDirections(rotations) {
  const P = [0, 1, 2].map((i) => [0, 1, 2].map((j) => rotations.reduce((s, R) => s + R[i][j], 0) / rotations.length));
  const basis = [];
  for (let j = 0; j < 3; j++) {
    const v = [P[0][j], P[1][j], P[2][j]];
    if (Math.hypot(...v) < 1e-9) continue;
    const trial = [...basis, v];
    // Keep v if it is independent of the basis so far (Gram determinant > 0).
    const gram = trial.map((u) => trial.map((x) => u[0] * x[0] + u[1] * x[1] + u[2] * x[2]));
    const det = trial.length === 1 ? gram[0][0]
      : trial.length === 2 ? gram[0][0] * gram[1][1] - gram[0][1] * gram[1][0]
        : gram[0][0] * (gram[1][1] * gram[2][2] - gram[1][2] * gram[2][1])
          - gram[0][1] * (gram[1][0] * gram[2][2] - gram[1][2] * gram[2][0])
          + gram[0][2] * (gram[1][0] * gram[2][1] - gram[1][1] * gram[2][0]);
    if (det > 1e-9) basis.push(v);
  }
  return basis;
}

/**
 * The shift onto ITA's standard origin for group `number`.
 * @param {{R:number[][], t:number[]}[]} ops  the exact operations in the group's standard
 *   (naming) cell, centring included
 * @param {number} number  ITA number
 * @param {number[][]} A  the standard cell's rows (Å), for "nearest"
 * @param {number[]} [offset]  the shift already applied (standard-cell fractions), so that
 *   the returned one moves the structure least counted from the .rmc6f origin
 * @returns {{ shift:number[], total:number[], totalA:number }|null}  shift p (x' = x + p);
 *   total = offset + p, reduced, and its Cartesian length. null when the rotations or the
 *   centring differ from ITA's, or no shift matches.
 */
export function itaOriginShift(ops, number, A, offset = [0, 0, 0]) {
  const spec = itaGenerators(number);
  const ita = itaOperations(number);
  if (!spec || !ita || !ops.length) return null;

  const ours = new Map();          // rotation → our translations
  for (const { R, t } of ops) {
    const k = rotKey(R);
    if (!ours.has(k)) ours.set(k, []);
    ours.get(k).push(t);
  }
  const theirs = new Map();
  for (const { R, t } of ita) {
    const k = rotKey(R);
    if (!theirs.has(k)) theirs.set(k, []);
    theirs.get(k).push(t);
  }
  if (ours.size !== theirs.size || [...ours.keys()].some((k) => !theirs.has(k))) return null;
  const centring = theirs.get(rotKey([[1, 0, 0], [0, 1, 0], [0, 0, 1]]));
  const ourPure = ours.get(rotKey([[1, 0, 0], [0, 1, 0], [0, 0, 1]])) ?? [];
  if (ourPure.length !== centring.length || !ourPure.every((t) => centring.some((c) => sameMod1(t, c)))) return null;
  for (const k of ours.keys()) if (ours.get(k).length !== theirs.get(k).length) return null;

  const fixedBasis = fixedDirections([...ours.keys()].map((k) => k.split(',').map(Number)).map((f) => [f.slice(0, 3), f.slice(3, 6), f.slice(6, 9)]));
  const valid = (p) => ops.every(({ R, t }) => {
    const moved = [0, 1, 2].map((i) => t[i] + p[i] - (R[i][0] * p[0] + R[i][1] * p[1] + R[i][2] * p[2]));
    return theirs.get(rotKey(R)).some((tau) => sameMod1(moved, tau));
  });

  // ITA's generators: with the shared centring they generate its group, so matching them
  // matches the group (every candidate is still checked against all operations).
  const genOps = spec.generators.map((text) => {
    const { R, t } = parseCoordinateForm(text);
    return { R, tau: t };
  });
  if (!genOps.length) {
    // P1: every shift matches; keep the structure where it is.
    const total = shortest(offset, centring, fixedBasis, A);
    return { shift: total.map((v, i) => v - offset[i]), total, totalA: cartLength(A, total) };
  }

  const M = [];
  for (const { R } of genOps) for (let i = 0; i < 3; i++) M.push([0, 1, 2].map((j) => (i === j ? 1 : 0) - R[i][j]));
  const { D, U, V, rank } = diagonalForm(M);

  const candidates = [];
  const choose = (g, picks) => {
    if (g === genOps.length) {
      const rhs = [];
      genOps.forEach(({ R, tau }, k) => {
        const t = ours.get(rotKey(R))[0];
        for (let i = 0; i < 3; i++) rhs.push(tau[i] - t[i] + picks[k][i]);
      });
      const Ur = U.map((row) => row.reduce((s, u, j) => s + u * rhs[j], 0));
      for (let i = rank; i < Ur.length; i++) if (Math.abs(cyc(Ur[i])) > 1e-6) return;
      const ranges = [0, 1, 2].map((i) => (i < rank ? D[i][i] : 1));
      for (let a = 0; a < ranges[0]; a++) for (let b = 0; b < ranges[1]; b++) for (let c = 0; c < ranges[2]; c++) {
        const m = [a, b, c];
        const q = [0, 1, 2].map((i) => (i < rank ? (Ur[i] + m[i]) / D[i][i] : 0));
        candidates.push(mulV(V, q));
      }
      return;
    }
    for (const c of centring) choose(g + 1, [...picks, c]);
  };
  choose(0, []);

  let best = null;
  for (const p of candidates) {
    if (!valid(p)) continue;
    const total = shortest(p.map((v, i) => v + offset[i]), centring, fixedBasis, A);
    const shift = total.map((v, i) => v - offset[i]);
    if (!valid(shift)) continue;
    const totalA = cartLength(A, total);
    if (!best || totalA < best.totalA - 1e-9) best = { shift, total, totalA };
  }
  return best;
}
