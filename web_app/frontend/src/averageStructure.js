// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/averageStructure.js
//
// The structure behind the Detected SG card's "Download CIF": the RMC model folded into one
// unit cell and averaged over every copy of each site (browserData.structureFromRmc6f), then
// averaged over the orbits of the space group picked on the tolerance ladder, and written in
// the standard cell that group was named in. Coordinates are measured averages, not the
// idealized values of a structure type: a free coordinate keeps its measured value, and only
// what the group fixes (a special position, the site-symmetry form of U) is exact.
//
//   1. Exact operations (exactGroup). The finder's translations are least-squares estimates
//      (symmetry.js refineOperation), so the product of two operations matches a third only
//      within the residual. Each product defines an integer 2-cocycle
//      n(g,h) = round(t_g + R_g·t_h − t_gh); s_g = (1/|G|)·Σ_h n(g,h) solves
//      s_g + R_g·s_h − s_gh = n(g,h) exactly, and the detected translations differ from s by
//      a coboundary (R_g − I)·o plus noise, with o = −(1/|G|)·Σ_g (t_g − s_g).
//      t̂_g = s_g + (R_g − I)·o is the exact group nearest the detected one, on the data's
//      own origin.
//   2. Origin (chooseOrigin). t̂ is exact but in general irrational. The origin moves by the
//      smallest δ (Cartesian) that puts every translation of the output cell on the 1/48
//      grid, so the operations print as x+1/2, -y+1/4, … The candidates are the data's own
//      origin with its translations snapped (δ of the order of the noise for a box built on
//      a standard origin) and points on the symmetry elements and their intersections
//      (where International Tables puts origins), for a box on any other origin. For a group
//      named in a standard cell the structure then moves on to ITA's OWN origin
//      (itaOrigin.js, against itaOperations.js): the equivalent one nearest the .rmc6f
//      origin, none along a polar axis. The CIF then lists ITA's operations and every orbit's
//      Wyckoff letter is read exactly from the table. The total shift is reported in the CIF.
//   3. Orbit averaging (averageOrbit). Every operation g carries the orbit representative to
//      the nearest orbit member; that member's mean position and covariance are carried back
//      by g⁻¹ and averaged, weighted by the member's atom count. The terms are permuted by the
//      representative's stabilizer, so the average lies exactly on its special position and
//      U has exactly its site-symmetry form (an explicit stabilizer projection guards the
//      round-off). U is the pooled second moment of every atom of the orbit about the
//      symmetrized position: within-site spread plus the scatter of the member means, which
//      is how a distortion the picked group averages away shows up.
//   4. Output (symmetryAveragedStructure). Positions, U and the operations are carried into
//      the standard cell of the group (setting.Q from spaceGroupSymbol.js), the metric is
//      averaged over the point group (Rᵀ·G·R), and each orbit is written once, at the image
//      that fits its tabulated Wyckoff form when there is one.
//
// Fractional coordinates are COLUMN vectors, as in symmetry.js: x' = R·x + t; cell rows
// A = [a, b, c] (Å). Covariances are fractional (cell fractions²) until the CIF's U_ij.

import { analyseSymmetry } from './symmetryModel.js';
import { applySetting, POINT_GROUP_SYSTEM } from './spaceGroupSymbol.js';
import { inv3 } from './symmetry.js';
import { wyckoffPositions, fitsForm } from './wyckoff.js';
import { ORIGIN_CHOICE_2, itaGenerators, itaOperations } from './itaOperations.js';
import { itaOriginShift } from './itaOrigin.js';

const I3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const ZERO3 = () => [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
const mul = (X, Y) => X.map((row) => [0, 1, 2].map((j) => row[0] * Y[0][j] + row[1] * Y[1][j] + row[2] * Y[2][j]));
const mulV = (M, v) => M.map((row) => row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
const transpose = (M) => [0, 1, 2].map((i) => [0, 1, 2].map((j) => M[j][i]));
const add = (u, v) => [u[0] + v[0], u[1] + v[1], u[2] + v[2]];
const sub = (u, v) => [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
const sub3 = (X, Y) => X.map((row, i) => row.map((v, j) => v - Y[i][j]));
const cyc = (x) => x - Math.round(x);
const rotKey = (R) => R.flat().join(',');
const intInverse = (R) => inv3(R).map((row) => row.map((v) => Math.round(v)));
const isIdentity = (R) => rotKey(R) === '1,0,0,0,1,0,0,0,1';

/** x mod 1 in [0, 1), with values within 1e-9 of an integer set to 0. */
export const wrapTidy = (x) => {
  const y = x - Math.floor(x);
  return y < 1e-9 || y > 1 - 1e-9 ? 0 : y;
};

/** Cartesian length (Å) of a fractional vector d in the cell with rows A. */
const cartLength = (A, d) => Math.hypot(
  d[0] * A[0][0] + d[1] * A[1][0] + d[2] * A[2][0],
  d[0] * A[0][1] + d[1] * A[1][1] + d[2] * A[2][1],
  d[0] * A[0][2] + d[1] * A[1][2] + d[2] * A[2][2],
);
const cycDist = (A, u, v) => cartLength(A, [cyc(u[0] - v[0]), cyc(u[1] - v[1]), cyc(u[2] - v[2])]);

/** Lattice rows of the cell with basis Q (columns, in the fractional basis of A). */
export const cellRows = (A, Q) => [0, 1, 2].map((i) => [0, 1, 2].map((k) => Q[0][i] * A[0][k] + Q[1][i] * A[1][k] + Q[2][i] * A[2][k]));

/* ── 1. exact operations ─────────────────────────────────────────────────── */

/**
 * The exact space group nearest a detected one. `ops` is a closed group {R, t} in the
 * fractional basis of A (rows, Å), t the refined translations (any lift).
 *
 * @returns {{ ops:{R,t}[], s:number[][], origin:number[], lifts:number[][], defect:number }|null}
 *   ops    — exact operations, same order (t unwrapped: s + (R − I)·origin)
 *   s      — the exact translations at the cocycle's own origin (denominators divide |G|)
 *   origin — o, the data's origin relative to s's
 *   lifts  — the detected translations the cocycle was computed from
 *   defect — the largest |t_g + R_g·t_h − t_gh − n(g,h)|, cell fractions: how far the
 *            detected set is from closing exactly
 * null when a product's rotation is missing or the products miss by a quarter of a cell or
 * more (not a group).
 */
export function exactGroup(ops, A) {
  const n = ops.length;
  if (!n) return null;
  const byR = new Map();
  ops.forEach((o, k) => {
    const key = rotKey(o.R);
    if (!byR.has(key)) byR.set(key, []);
    byR.get(key).push(k);
  });
  const s = ops.map(() => [0, 0, 0]);
  let defect = 0;
  for (let i = 0; i < n; i++) {
    const a = ops[i];
    for (let j = 0; j < n; j++) {
      const b = ops[j];
      const same = byR.get(rotKey(mul(a.R, b.R)));
      if (!same) return null;
      const c = add(mulV(a.R, b.t), a.t);
      let k = same[0];
      let best = Infinity;
      for (const m of same) {
        const d = cycDist(A, c, ops[m].t);
        if (d < best) { best = d; k = m; }
      }
      for (let q = 0; q < 3; q++) {
        const v = c[q] - ops[k].t[q];
        const r = Math.round(v);
        defect = Math.max(defect, Math.abs(v - r));
        s[i][q] += r;
      }
    }
  }
  if (!(defect < 0.25)) return null;
  for (const v of s) for (let q = 0; q < 3; q++) v[q] /= n;
  const sum = [0, 0, 0];
  ops.forEach((o, g) => { for (let q = 0; q < 3; q++) sum[q] += o.t[q] - s[g][q]; });
  const origin = sum.map((v) => -v / n);
  return {
    ops: ops.map((o, g) => ({ R: o.R, t: add(s[g], sub(mulV(o.R, origin), origin)) })),
    s,
    origin,
    lifts: ops.map((o) => o.t.slice()),
    defect,
  };
}

/* ── 2. origin ───────────────────────────────────────────────────────────── */

// Denominators of the translations ITA prints (all divide 48).
const NICE_DENOMINATORS = [1, 2, 3, 4, 6, 8, 12, 16, 24, 48];

/** Smallest crystallographic denominator d (divides 48) with x·d an integer, else null. */
export function niceDenominator(x) {
  for (const d of NICE_DENOMINATORS) if (Math.abs(cyc(x * d)) < 1e-6) return d;
  return null;
}

// Order n of a rotation (R^n = I, n ≤ 6) and the intrinsic (screw / glide) part of {R | t}:
// w = (1/n)·Σ_{k<n} R^k·t, the translation of {R | t}^n shared out over its n steps.
function intrinsicPart(R, t) {
  let power = R;
  let sum = t.slice();
  let n = 1;
  while (!isIdentity(power) && n < 6) {
    const step = mulV(power, t);
    sum = add(sum, step);
    power = mul(power, R);
    n += 1;
  }
  return sum.map((v) => v / n);
}

// Least-squares point f with M_k·f = l_k for every (M, l) in `rows` — a point on a symmetry
// element, or where two meet. A small ridge keeps a rank-deficient system (a rotation axis,
// a mirror plane) well posed and picks the solution nearest the current origin.
function leastSquaresPoint(rows) {
  const N = ZERO3();
  const b = [0, 0, 0];
  for (const { M, l } of rows) {
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) N[i][j] += M[0][i] * M[0][j] + M[1][i] * M[1][j] + M[2][i] * M[2][j];
      b[i] += M[0][i] * l[0] + M[1][i] * l[1] + M[2][i] * l[2];
    }
  }
  for (let i = 0; i < 3; i++) N[i][i] += 1e-10;
  return mulV(inv3(N), b);
}

/**
 * The origin shift δ (given-cell fractions; x' = x + δ) that puts every translation of the
 * output cell's operations on the 1/48 grid — the smallest such shift (Cartesian, through
 * A) among the candidates:
 *   • none (the exact group on the data's own origin),
 *   • the data's origin with every detected translation snapped to the grid — the answer
 *     for a box built on a standard origin, δ of the order of the noise,
 *   • a point on each symmetry element, and where two elements of different rotations
 *     meet — where International Tables puts its origins, for a box on any origin.
 * `exact` is exactGroup's result; Qinv maps given fractions to output fractions.
 *
 * @returns {{ shift:number[], shiftA:number, nice:boolean }}  nice false (and no shift)
 *   when no candidate is on the grid: the operations then print with decimals.
 */
export function chooseOrigin(exact, A, Qinv) {
  const { ops } = exact;
  const onGrid = (delta) => ops.every(({ R, t }) => mulV(Qinv, add(t, sub(delta, mulV(R, delta))))
    .every((x) => niceDenominator(x) !== null));

  const candidates = [[0, 0, 0]];
  const snapped = exact.lifts.map((t) => t.map((v) => Math.round(v * 48) / 48));
  const mean = [0, 1, 2].map((q) => snapped.reduce((sum, t) => sum + t[q], 0) / snapped.length);
  candidates.push(add(mean, exact.origin));
  // The origin moved onto an element (x' = x − f) leaves that operation its intrinsic part.
  const elements = [];
  const seen = new Set();
  for (const { R, t } of ops) {
    if (isIdentity(R)) continue;
    const element = { M: sub3(I3, R), l: sub(t, intrinsicPart(R, t)) };
    candidates.push(leastSquaresPoint([element]).map((v) => -v));
    if (!seen.has(rotKey(R))) { seen.add(rotKey(R)); elements.push(element); }
  }
  for (let i = 0; i < elements.length; i++) {
    for (let j = i + 1; j < elements.length; j++) candidates.push(leastSquaresPoint([elements[i], elements[j]]).map((v) => -v));
  }

  const best = candidates
    .map((delta) => ({ delta, length: cartLength(A, delta) }))
    .sort((p, q) => p.length - q.length)
    .find(({ delta }) => onGrid(delta));
  if (!best) return { shift: [0, 0, 0], shiftA: 0, nice: false };
  return { shift: best.delta, shiftA: best.length, nice: true };
}

/* ── 3. orbit averaging ──────────────────────────────────────────────────── */

/**
 * One orbit's symmetry average, in the given cell. `ops` are exact operations; `sites` the
 * orbit members [{ pos, cov, weight }] — pos the mean position (fractional, any lift), cov
 * its fractional covariance, weight its atom count; member 0 is the representative. `tol`
 * (Å) is the tolerance the group was picked at: an operation that maps the averaged
 * position within it of itself belongs to its site symmetry.
 *
 * @returns {{ x:number[], V:number[][], maxShift:number, rmsShift:number, matched:number }}
 *   x, V — the representative's averaged position (near member 0) and pooled covariance;
 *   maxShift / rmsShift — how far the members' means lie from their symmetrized positions
 *   (Å); matched — members some operation reached (the rest do not contribute).
 */
export function averageOrbit(ops, sites, A, tol) {
  const xr = sites[0].pos;
  const terms = ops.map((g) => {
    const p = add(mulV(g.R, xr), g.t);
    let member = 0;
    let bestLength = Infinity;
    let bestPoint = null;
    sites.forEach((site, s) => {
      const d = [cyc(site.pos[0] - p[0]), cyc(site.pos[1] - p[1]), cyc(site.pos[2] - p[2])];
      const length = cartLength(A, d);
      if (length < bestLength) { bestLength = length; member = s; bestPoint = add(p, d); }
    });
    const Ri = intInverse(g.R);
    return { member, Ri, back: mulV(Ri, sub(bestPoint, g.t)) };
  });

  let weight = 0;
  const sum = [0, 0, 0];
  for (const { member, back } of terms) {
    const w = sites[member].weight;
    weight += w;
    for (let q = 0; q < 3; q++) sum[q] += w * back[q];
  }
  let x = sum.map((v) => v / weight);

  let V = ZERO3();
  let maxShift = 0;
  let sumSq = 0;
  for (const { member, Ri, back } of terms) {
    const w = sites[member].weight;
    const C = mul(mul(Ri, sites[member].cov), transpose(Ri));
    const d = sub(back, x);
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) V[i][j] += w * (C[i][j] + d[i] * d[j]);
    const shift = cartLength(A, d);
    maxShift = Math.max(maxShift, shift);
    sumSq += w * shift * shift;
  }
  V = V.map((row) => row.map((v) => v / weight));

  // The stabilizer projection: exact invariance whatever the round-off.
  const stabilizer = ops.filter((h) => cycDist(A, add(mulV(h.R, x), h.t), x) <= tol);
  if (stabilizer.length > 1) {
    const px = [0, 0, 0];
    let pV = ZERO3();
    for (const h of stabilizer) {
      const p = add(mulV(h.R, x), h.t);
      for (let q = 0; q < 3; q++) px[q] += p[q] + Math.round(x[q] - p[q]);
      const C = mul(mul(h.R, V), transpose(h.R));
      pV = pV.map((row, i) => row.map((v, j) => v + C[i][j]));
    }
    x = px.map((v) => v / stabilizer.length);
    V = pV.map((row) => row.map((v) => v / stabilizer.length));
  }
  return {
    x,
    V,
    maxShift,
    rmsShift: Math.sqrt(sumSq / weight),
    matched: new Set(terms.map((t) => t.member)).size,
  };
}

/* ── 4. output ───────────────────────────────────────────────────────────── */

/** Cell parameters (Å, degrees) and volume (Å³) from a metric tensor G = A·Aᵀ. */
export function cellParameters(G) {
  const a = Math.sqrt(G[0][0]);
  const b = Math.sqrt(G[1][1]);
  const c = Math.sqrt(G[2][2]);
  const angle = (x) => (Math.acos(Math.max(-1, Math.min(1, x))) * 180) / Math.PI;
  const det = G[0][0] * (G[1][1] * G[2][2] - G[1][2] * G[2][1])
    - G[0][1] * (G[1][0] * G[2][2] - G[1][2] * G[2][0])
    + G[0][2] * (G[1][0] * G[2][1] - G[1][1] * G[2][0]);
  return {
    a, b, c,
    alpha: angle(G[1][2] / (b * c)),
    beta: angle(G[0][2] / (a * c)),
    gamma: angle(G[0][1] / (a * b)),
    volume: Math.sqrt(Math.max(det, 0)),
  };
}

// The operations of the output cell as coset representatives (one per rotation, identity
// first, its translation the smallest of its coset) × centring translations.
function orderedOperations(ops, translations) {
  const reps = new Map();
  for (const { R, t } of ops) {
    const key = rotKey(R);
    // The lexicographically smallest translation of the coset (the identity gets 0).
    const best = translations
      .map((tau) => t.map((v, i) => wrapTidy(v + tau[i])))
      .reduce((m, u) => (u[0] - m[0] || u[1] - m[1] || u[2] - m[2]) < 0 ? u : m);
    if (!reps.has(key)) reps.set(key, { R, t: best });
  }
  const det = (R) => R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1])
    - R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0])
    + R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]);
  const sorted = [...reps.values()].sort((p, q) => (isIdentity(q.R) - isIdentity(p.R))
    || (det(q.R) - det(p.R)) || (rotKey(p.R) < rotKey(q.R) ? -1 : rotKey(p.R) > rotKey(q.R) ? 1 : 0));
  const centring = translations.map((tau) => tau.map(wrapTidy))
    .sort((u, v) => u[0] - v[0] || u[1] - v[1] || u[2] - v[2]);
  const out = [];
  for (const tau of centring) for (const { R, t } of sorted) out.push({ R, t: t.map((v, i) => wrapTidy(v + tau[i])) });
  return { ops: out, rotations: sorted.map((o) => o.R) };
}

// Distinct images of x under the operations (mod 1), each with the rotation that made it.
function orbitImages(ops, x) {
  const images = [];
  for (const { R, t } of ops) {
    const p = add(mulV(R, x), t).map(wrapTidy);
    if (!images.some((q) => Math.abs(cyc(q.x[0] - p[0])) < 1e-6 && Math.abs(cyc(q.x[1] - p[1])) < 1e-6
      && Math.abs(cyc(q.x[2] - p[2])) < 1e-6)) images.push({ x: p, R });
  }
  return images;
}

// The representative written to the CIF: among the images that fit the orbit's tabulated
// Wyckoff form (when it has a letter), the lexicographically smallest.
function pickRepresentative(images, spaceGroupNumber, letter) {
  const row = letter && spaceGroupNumber ? wyckoffPositions(spaceGroupNumber).find((r) => r.letter === letter) : null;
  const fitting = row ? images.filter((im) => fitsForm(row.form, im.x, 1e-5)) : [];
  const pool = fitting.length ? fitting : images;
  const key = (v) => v.map((c) => Math.round(c * 1e7));
  return pool.reduce((best, im) => {
    const a = key(im.x), b = key(best.x);
    return (a[0] - b[0] || a[1] - b[1] || a[2] - b[2]) < 0 ? im : best;
  });
}

const gcd = (a, b) => (b ? gcd(b, a % b) : a);

// On ITA's origin a Wyckoff position is read straight from the table: the one row of the
// orbit's multiplicity whose coordinate form some image fits exactly (two positions never
// share a point, so at most one fits).
function exactWyckoffLetter(number, images) {
  const rows = wyckoffPositions(number)
    .filter((row) => row.multiplicity === images.length && images.some((im) => fitsForm(row.form, im.x, 1e-6)));
  return rows.length === 1 ? rows[0].letter : null;
}

/**
 * The symmetry-averaged structure of `structure` (browserData.structureFromRmc6f) in the
 * space group the finder reports at `tol` (Å) — what the Detected SG card shows when that
 * tolerance is picked. Throws an Error with a user-facing message when there is nothing to
 * export (no basis, a structure the finder skips, no operation, a set that is not a group).
 *
 * @returns {{
 *   spaceGroup: { label, symbol, number, pointGroup, system, centring, standard },
 *   cell: { a, b, c, alpha, beta, gamma, volume },
 *   operations: { R:number[][], t:number[] }[],
 *   sites: { element, elements:{element, occupancy}[], multiplicity, wyckoff, x:number[],
 *            U:number[][], Ueq:number }[],
 *   formula: { counts: Record<string, number>, Z:number },
 *   provenance: { tol, maxResidual, nSpace, Q, ratio, originShift, originShiftA, niceOrigin,
 *                 maxShiftA, rmsShiftA, unmatched, sites, atoms, supercell, source }
 * }}
 *   U is in Å² on the CIF's convention (U_ij = ⟨Δx_i Δx_j⟩ / (a*_i a*_j)), Ueq = tr(U·G)/3.
 */
export function symmetryAveragedStructure(structure, tol) {
  const analysis = analyseSymmetry(structure, tol);
  if (!analysis) throw new Error('This structure has no average-structure basis to export.');
  if (analysis.skipped) throw new Error(analysis.reason);
  const { A, sg, found, positions } = analysis;
  if (!sg.ops.length || !found.length) throw new Error('No symmetry operation fits at this tolerance.');
  const exact = exactGroup(sg.ops, A);
  if (!exact) throw new Error(`The operations of ${sg.spaceGroup} do not close into an exact group.`);

  // The output cell: the standard cell the group was named in, else the .rmc6f unit cell.
  const standard = !!sg.setting;
  const Q = standard ? sg.setting.Q : I3;
  const Qinv = standard ? sg.setting.Qinv : I3;
  const Aout = cellRows(A, Q);
  const origin = chooseOrigin(exact, A, Qinv);
  const delta = origin.shift;

  // The exact operations on the shifted origin (x' = x + δ): t' = t̂ + (I − R)·δ.
  const ops = exact.ops.map(({ R, t }) => ({ R, t: add(t, sub(delta, mulV(R, delta))) }));
  const pure = [[0, 0, 0]];
  for (const { R, t } of ops) {
    if (!isIdentity(R)) continue;
    const w = t.map(wrapTidy);
    if (!pure.some((u) => u.every((v, i) => Math.abs(cyc(v - w[i])) < 1e-6))) pure.push(w);
  }
  const setting = applySetting(ops.map(({ R, t }) => ({ R, t: t.map(wrapTidy) })), pure, Q);
  if (!setting) throw new Error(`The operations of ${sg.spaceGroup} cannot be carried into its standard cell.`);
  // ITA's own origin (itaOrigin.js): of the equivalent ones, the one nearest the .rmc6f
  // origin. The CIF then lists ITA's operations themselves. When the group has no standard
  // cell, or its operations cannot be matched to ITA's, the 1/48-grid origin stays.
  const number = standard ? sg.spaceGroupNumber : null;
  const ita = number ? itaOriginShift(setting.ops, number, Aout, mulV(Qinv, delta)) : null;
  const shiftOut = ita ? ita.shift : [0, 0, 0];
  const { ops: outOps, rotations } = ita
    ? orderedOperations(itaOperations(number), itaGenerators(number).centring)
    : orderedOperations(setting.ops, setting.translations);

  // Metric averaged over the point group: a cell that holds the group within the tolerance
  // but not exactly (a strained lattice) is written with the group's own metric.
  const G0 = mul(Aout, transpose(Aout));
  let G = ZERO3();
  for (const R of rotations) {
    const C = mul(mul(transpose(R), G0), R);
    G = G.map((row, i) => row.map((v, j) => v + C[i][j] / rotations.length));
  }
  const Gstar = inv3(G);
  const aStar = [0, 1, 2].map((i) => Math.sqrt(Gstar[i][i]));

  const nCells = (structure.supercell || [1, 1, 1]).reduce((p, v) => p * Math.max(v, 1), 1);
  let maxShift = 0;
  let shiftSq = 0;
  let shiftWeight = 0;
  let unmatched = 0;
  const sites = found.map((orbit, i) => {
    const members = orbit.index.map((k) => {
      const b = structure.basis[k];
      return {
        pos: add(b.mean ?? b.frac, delta),
        cov: b.covFrac ?? ZERO3(),
        weight: b.count ?? 1,
        elementCounts: b.elementCounts ?? { [b.el]: b.count ?? 1 },
      };
    });
    const avg = averageOrbit(ops, members, A, tol);
    const weight = members.reduce((w, m) => w + m.weight, 0);
    maxShift = Math.max(maxShift, avg.maxShift);
    shiftSq += avg.rmsShift ** 2 * weight;
    shiftWeight += weight;
    unmatched += members.length - avg.matched;

    const xOut = add(mulV(Qinv, avg.x), shiftOut);
    const VOut = mul(mul(Qinv, avg.V), transpose(Qinv));
    const images = orbitImages(outOps, xOut);
    // A letter goes with its multiplicity: dropped if the averaged orbit no longer has it.
    const letter = ita
      ? exactWyckoffLetter(number, images)
      : (positions[i].letter && positions[i].multiplicity === images.length ? positions[i].letter : null);
    const rep = pickRepresentative(images, number, letter);
    const V = mul(mul(rep.R, VOut), transpose(rep.R));
    const U = V.map((row, r) => row.map((v, c) => v / (aStar[r] * aStar[c])));
    const Ueq = (V[0][0] * G[0][0] + V[1][1] * G[1][1] + V[2][2] * G[2][2]
      + 2 * (V[0][1] * G[0][1] + V[0][2] * G[0][2] + V[1][2] * G[1][2])) / 3;

    const counts = {};
    for (const m of members) for (const [el, c] of Object.entries(m.elementCounts)) counts[el] = (counts[el] || 0) + c;
    const elements = Object.entries(counts)
      .filter(([, c]) => c > 0)
      .sort(([p, cp], [q, cq]) => cq - cp || (p < q ? -1 : 1))
      .map(([element, c]) => ({ element, occupancy: c / (members.length * nCells) }));
    return {
      element: elements[0]?.element ?? orbit.element,
      elements,
      multiplicity: images.length,
      wyckoff: letter,
      x: rep.x,
      U,
      Ueq,
    };
  });

  // Write order: by element, then special positions first, then coordinates.
  sites.sort((p, q) => (p.element < q.element ? -1 : p.element > q.element ? 1 : 0)
    || p.multiplicity - q.multiplicity
    || p.x[0] - q.x[0] || p.x[1] - q.x[1] || p.x[2] - q.x[2]);

  const cellCounts = {};
  for (const site of sites) {
    for (const { element, occupancy } of site.elements) {
      cellCounts[element] = (cellCounts[element] || 0) + site.multiplicity * occupancy;
    }
  }
  // Z: the common factor of the cell's atom counts, or, with partial occupancies, of the
  // site multiplicities (Na0.94K0.06Cl in rocksalt is still Z = 4).
  const integral = Object.values(cellCounts).every((v) => Math.abs(v - Math.round(v)) < 1e-3);
  const Z = (integral
    ? Object.values(cellCounts).map((v) => Math.round(v))
    : sites.map((site) => site.multiplicity)).reduce(gcd, 0) || 1;

  return {
    spaceGroup: {
      label: sg.spaceGroup,
      symbol: standard ? sg.spaceGroup : null,
      number: standard ? sg.spaceGroupNumber : null,
      pointGroup: sg.pointGroup,
      system: POINT_GROUP_SYSTEM[sg.pointGroup] ?? null,
      centring: setting.letter,
      standard,
    },
    cell: cellParameters(G),
    operations: outOps,
    sites,
    formula: { counts: Object.fromEntries(Object.entries(cellCounts).map(([el, v]) => [el, v / Z])), Z },
    provenance: {
      tol,
      maxResidual: sg.maxResidual,
      nSpace: sg.nSpace,
      Q,
      ratio: setting.ratio,
      originShift: add(delta, mulV(Q, shiftOut)),
      originShiftA: ita ? ita.totalA : origin.shiftA,
      niceOrigin: origin.nice || !!ita,
      itaOrigin: !!ita,
      originChoice: ita && ORIGIN_CHOICE_2.has(number) ? 2 : null,
      maxShiftA: maxShift,
      rmsShiftA: Math.sqrt(shiftSq / Math.max(shiftWeight, 1e-300)),
      unmatched,
      sites: structure.basis.length,
      atoms: structure.totalAtoms ?? null,
      supercell: structure.supercell ?? null,
      source: structure.source ?? null,
    },
  };
}
