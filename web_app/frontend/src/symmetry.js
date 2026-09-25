// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/symmetry.js
//
// Pure, fully-offline crystal symmetry finder (a lightweight FINDSYM-like tool),
// ported from the RMC-phonon-dynamics app (web/src/math/symmetry.js). Given the
// conventional cell A_conv (rows = lattice vectors, Å) + a basis and a tolerance,
// it returns the space-group operations {R|t} that map the basis onto itself,
// plus the residual of that fit, so the symmetry can be traced as a function of
// tolerance on a disordered RMC average structure — no external dependency, no
// WASM, runs client-side in the static dashboard.
//
// Method (spglib-lite, bounded):
//   1. Point operations = integer matrices R (|det R| = 1; entries in {-1,0,1} in a reduced
//      basis of the lattice, which finds all of them on any cell) whose Cartesian lattice
//      strain is within the tolerance (latticeStrain).
//   2. For each R, candidate translations are seeded from atom images, t = x_b − R·x_a0,
//      refined by least squares over every site, and kept when {R|t} maps every atom
//      onto a same-element atom within `tol` (cartesian Å); the residual is the worst site.
//   3. Along the tolerance ladder only CLOSED groups are reported (groupsByThreshold).
//   4. A group is named by spaceGroupSymbol.js in a standard setting it searches for
//      (axis orders, cells built from the symmetry elements, supercell reduction), reading
//      the screw/glide part of each t, and checked against the 230-group table of
//      spaceGroupTable.js. What cannot be named reliably is reported as its crystal class.
//
// What this still does NOT do, and FINDSYM does: no origin shift, no search in the
// primitive cell of the crystal's own translation lattice (rotations that do not map the
// GIVEN cell's lattice onto itself — the cubic 3-folds of a 2×2×1 supercell — are never
// tested; see spaceGroupHM's lower-bound flag), and no idealized structure output.
//
// Fractional coords are COLUMN vectors here: x' = R·x + t. Lattice rows: A = [a1,a2,a3].

import { hmSymbolInStandardSetting, centeringOfOps, latticeFullyTested, reduceBasis } from './spaceGroupSymbol.js';
import { spaceGroupNumber, canonicalSymbol, pointGroupOfSymbol } from './spaceGroupTable.js';

/** Determinant of a 3×3 matrix (rows). */
export function det3(m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

const wrap01 = (x) => x - Math.floor(x);
const nearestInt = (x) => x - Math.round(x);
const meanEdge = (A) => (Math.hypot(A[0][0], A[0][1], A[0][2]) + Math.hypot(A[1][0], A[1][1], A[1][2]) + Math.hypot(A[2][0], A[2][1], A[2][2])) / 3;

/** Metric tensor G = A·Aᵀ (A rows = lattice vectors, cartesian). */
export function metricTensor(A) {
  const G = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++)
    G[i][j] = A[i][0] * A[j][0] + A[i][1] * A[j][1] + A[i][2] * A[j][2];
  return G;
}

// Rᵀ · G · R for a 3×3 integer R and symmetric G.
function conjugate(R, G) {
  // M = Rᵀ G R.  (Rᵀ)_{ij} = R_{ji}
  const GR = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++)
    GR[i][j] = G[i][0] * R[0][j] + G[i][1] * R[1][j] + G[i][2] * R[2][j];
  const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++)
    M[i][j] = R[0][i] * GR[0][j] + R[1][i] * GR[1][j] + R[2][i] * GR[2][j];
  return M;
}

/** Inverse of a 3×3 matrix (NaN/Infinity entries for a singular one). */
export function inv3(m) {
  const d = det3(m);
  const c = (i, j) => {
    const r0 = i === 0 ? 1 : 0, r1 = i === 2 ? 1 : 2, s0 = j === 0 ? 1 : 0, s1 = j === 2 ? 1 : 2;
    return ((i + j) % 2 ? -1 : 1) * (m[r0][s0] * m[r1][s1] - m[r0][s1] * m[r1][s0]);
  };
  return [[c(0, 0) / d, c(1, 0) / d, c(2, 0) / d], [c(0, 1) / d, c(1, 1) / d, c(2, 1) / d], [c(0, 2) / d, c(1, 2) / d, c(2, 2) / d]];
}

/**
 * Cartesian lattice strain of the point operation R (Å): how far the cell-edge vectors
 * are displaced when R is used as if it were an isometry of the lattice.
 *
 * With x' = R·x acting on fractional columns, the Cartesian map is M = Aᵀ·R·A⁻ᵀ and its
 * Green strain is E = ½(MᵀM − I) = ½·A⁻¹·D·A⁻ᵀ with D = RᵀGR − G. Edge a_i = Aᵀe_i is
 * displaced by |E·a_i| = ½·|A⁻¹·D·e_i|; the strain is the largest of the three. It is 0
 * for an exact lattice symmetry and is on the same Å scale as an atomic-position residual,
 * so a strained cell is tested against the same tolerance as the atoms.
 */
export function latticeStrain(A, R, G = metricTensor(A), Ainv = inv3(A)) {
  const M = conjugate(R, G);
  let worst = 0;
  for (let i = 0; i < 3; i++) {
    const d0 = M[0][i] - G[0][i], d1 = M[1][i] - G[1][i], d2 = M[2][i] - G[2][i];
    const v = 0.5 * Math.hypot(
      Ainv[0][0] * d0 + Ainv[0][1] * d1 + Ainv[0][2] * d2,
      Ainv[1][0] * d0 + Ainv[1][1] * d1 + Ainv[1][2] * d2,
      Ainv[2][0] * d0 + Ainv[2][1] * d1 + Ainv[2][2] * d2,
    );
    if (!(v <= worst)) worst = v;          // NaN propagates, so a broken lattice rejects all
  }
  return worst;
}

// reduceBasis() lives in spaceGroupSymbol.js (the naming step needs it too); re-exported here.
export { reduceBasis };

// Every lattice rotation R (|det R| = 1) whose lattice strain is ≤ tol (Å), with that
// strain, as integer matrices in the GIVEN basis. They are enumerated as {-1,0,1}
// matrices in a reduced basis of the lattice (reduceBasis), where that range is
// complete, and carried back, so an oblique cell loses none of them. The strain is a
// floor on the operation's residual, so a strained cell shows up in the ladder at the
// tolerance that absorbs the strain.
function latticeCandidates(A, tol) {
  const M = reduceBasis(A);
  const Ar = [0, 1, 2].map((i) => [0, 1, 2].map((k) => M[i][0] * A[0][k] + M[i][1] * A[1][k] + M[i][2] * A[2][k]));
  const Mt = [[M[0][0], M[1][0], M[2][0]], [M[0][1], M[1][1], M[2][1]], [M[0][2], M[1][2], M[2][2]]];
  const Mi = inv3(M).map((row) => row.map((x) => Math.round(x)));
  const MiT = [[Mi[0][0], Mi[1][0], Mi[2][0]], [Mi[0][1], Mi[1][1], Mi[2][1]], [Mi[0][2], Mi[1][2], Mi[2][2]]];
  const mul = (X, Y) => X.map((row) => [0, 1, 2].map((j) => row[0] * Y[0][j] + row[1] * Y[1][j] + row[2] * Y[2][j]));
  const G = metricTensor(Ar);
  const Ainv = inv3(Ar);
  const out = [];
  const v = [-1, 0, 1];
  const R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  // Iterate all 3^9 sign/zero patterns; keep unimodular ones that preserve the metric.
  for (let code = 0; code < 19683; code++) {
    let c = code;
    for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) { R[a][b] = v[c % 3]; c = (c / 3) | 0; }
    const d = det3(R);
    if (d !== 1 && d !== -1) continue;
    const strain = latticeStrain(Ar, R, G, Ainv);
    // x_reduced = M⁻ᵀ·x, so R in the given basis is Mᵀ·R·M⁻ᵀ.
    if (strain <= tol) out.push({ R: mul(mul(Mt, R), MiT).map((row) => row.map((x) => Math.round(x))), strain });
  }
  return out;
}

/**
 * Lattice point operations: integer R (|det| = 1; entries in {-1,0,1} in a reduced
 * basis) that preserve the metric to within a Cartesian lattice strain of `tol` Å (see
 * latticeStrain), in the basis of A.
 */
export function latticePointOps(A, tol = 0.01) {
  return latticeCandidates(A, tol).map(c => c.R);
}

// Apply R (integer) to a fractional column vector.
function applyR(R, x) {
  return [
    R[0][0] * x[0] + R[0][1] * x[1] + R[0][2] * x[2],
    R[1][0] * x[0] + R[1][1] * x[1] + R[1][2] * x[2],
    R[2][0] * x[0] + R[2][1] * x[1] + R[2][2] * x[2],
  ];
}

// Cartesian distance between two fractional points (minimal image), rows A.
function cartDist(fa, fb, A) {
  const d = [nearestInt(fa[0] - fb[0]), nearestInt(fa[1] - fb[1]), nearestInt(fa[2] - fb[2])];
  const c = [
    d[0] * A[0][0] + d[1] * A[1][0] + d[2] * A[2][0],
    d[0] * A[0][1] + d[1] * A[1][1] + d[2] * A[2][1],
    d[0] * A[0][2] + d[1] * A[1][2] + d[2] * A[2][2],
  ];
  return Math.hypot(c[0], c[1], c[2]);
}

// Minimum-image fractional offset (component-wise nearest integer) and its Cartesian length.
function offset(from, to, A) {
  const d = [nearestInt(to[0] - from[0]), nearestInt(to[1] - from[1]), nearestInt(to[2] - from[2])];
  const c0 = d[0] * A[0][0] + d[1] * A[1][0] + d[2] * A[2][0];
  const c1 = d[0] * A[0][1] + d[1] * A[1][1] + d[2] * A[2][1];
  const c2 = d[0] * A[0][2] + d[1] * A[1][2] + d[2] * A[2][2];
  return { d, dist: Math.hypot(c0, c1, c2) };
}

/**
 * Same-element sites binned on a fractional grid, so the partner search only visits the
 * bins an image within `radius` Å can reach. A component-wise minimum-image distance
 * d ≤ radius bounds every fractional component by radius·|b_i| (b_i = column i of A⁻¹,
 * the reciprocal vector), so with bins at least that wide the bin of the image and its
 * ±1 neighbours (cyclically) hold every candidate. Fewer than 3 bins along an axis → one.
 * `near(el, p)` returns the candidate sites for an image p.
 */
function partnerIndex(byEl, A, radius) {
  const Ainv = inv3(A);
  const n = [0, 1, 2].map((i) => {
    const k = Math.floor(1 / (radius * Math.hypot(Ainv[0][i], Ainv[1][i], Ainv[2][i])));
    return k >= 3 ? Math.min(k, 48) : 1;
  });
  const bins = new Map();
  const cellOf = (p, i) => Math.min(n[i] - 1, Math.floor(wrap01(p[i]) * n[i]));
  for (const [el, sites] of byEl) {
    const cells = Array.from({ length: n[0] * n[1] * n[2] }, () => []);
    for (const s of sites) cells[cellOf(s.frac, 0) + n[0] * (cellOf(s.frac, 1) + n[1] * cellOf(s.frac, 2))].push(s);
    bins.set(el, cells);
  }
  const span = n.map((k) => (k === 1 ? [0] : [-1, 0, 1]));
  const out = [];
  return {
    near(el, p) {
      const cells = bins.get(el);
      if (!cells) return [];
      if (cells.length === 1) return cells[0];
      out.length = 0;
      const c0 = cellOf(p, 0), c1 = cellOf(p, 1), c2 = cellOf(p, 2);
      for (const d2 of span[2]) for (const d1 of span[1]) for (const d0 of span[0]) {
        const i0 = (c0 + d0 + n[0]) % n[0], i1 = (c1 + d1 + n[1]) % n[1], i2 = (c2 + d2 + n[2]) % n[2];
        for (const s of cells[i0 + n[0] * (i1 + n[1] * i2)]) out.push(s);
      }
      return out;
    },
  };
}

// Map every site through {R|t} and pair its image with the nearest same-element site
// (candidates from `index.near`, see partnerIndex). Returns null as soon as one image has
// no partner within `radius` Å; otherwise the mean fractional offset (partner − image)
// over all sites and the worst Cartesian distance.
function matchImages(R, t, basis, index, A, radius) {
  const mean = [0, 0, 0];
  let worst = 0;
  for (const s of basis) {
    const img = applyR(R, s.frac);
    img[0] += t[0]; img[1] += t[1]; img[2] += t[2];
    let best = Infinity;
    let bestD = null;
    for (const o of index.near(s.el, img)) {
      const { d, dist } = offset(img, o.frac, A);
      if (dist < best) { best = dist; bestD = d; }
    }
    if (!(best <= radius)) return null;
    mean[0] += bestD[0]; mean[1] += bestD[1]; mean[2] += bestD[2];
    if (best > worst) worst = best;
  }
  const n = basis.length;
  return { shift: [mean[0] / n, mean[1] / n, mean[2] / n], worst };
}

// Least-squares translation of {R|t0} and its residual. The seed t0 comes from ONE atom
// pair and carries both atoms' displacement; the translation that minimises the summed
// squared Cartesian mismatch over all matched pairs is t0 + mean(partner − image) (the
// metric is common to every pair). Re-pair and re-centre until the shift vanishes (at
// most four passes); the residual is the worst site's distance at the final translation. The seed is paired within 2·tol, since it can be off by the noise of the
// two atoms that defined it, but the refined operation must fit within tol. `!(x <= tol)`
// also rejects NaN (an unparseable lattice).
function refineOperation(R, t0, basis, index, A, tol) {
  let t = t0.slice();
  let m = matchImages(R, t, basis, index, A, 2 * tol);
  for (let pass = 0; pass < 4 && m; pass++) {
    const moved = Math.abs(m.shift[0]) + Math.abs(m.shift[1]) + Math.abs(m.shift[2]);
    t = [wrap01(t[0] + m.shift[0]), wrap01(t[1] + m.shift[1]), wrap01(t[2] + m.shift[2])];
    m = matchImages(R, t, basis, index, A, 2 * tol);
    if (moved < 1e-12) break;
  }
  if (!m || !(m.worst <= tol)) return null;
  return { t, residual: m.worst };
}

/**
 * Space-group operations of (A, basis) within a cartesian tolerance `tol` (Å).
 * basis: [{ el, frac:[x,y,z] }]. Returns operations + summary counts + residual.
 *
 * `latticeTol` (Å, default `tol`) bounds the Cartesian lattice strain of a point operation
 * (latticeStrain); an operation's residual is the larger of its lattice strain and its
 * worst atomic mismatch, so a strained cell is judged on the same scale as the atoms.
 *
 * @returns {{ ops:{R,t,residual}[], nSpace, nPoint, order, maxResidual }}
 *   nPoint : distinct rotation parts present in the space group (its point group).
 *   nSpace : total {R|t} (= point-group order × #centering-type cosets for the cell).
 */
export function findSpaceGroupOps(A, basis, tol = 0.1, latticeTol = tol) {
  if (!basis || basis.length === 0) return { ...UNDETERMINED, ops: [], order: 0, maxResidual: Number.NaN };
  const { ops, maxResidual } = detectOperations(A, basis, tol, latticeTol);
  // The raw set need not be closed (Step 11): classify the largest closed group in it.
  const walk = groupsByThreshold(ops, A, Infinity);
  if (!walk.length) return { ...UNDETERMINED, ops, order: ops.length, maxResidual: Number.NaN };
  const group = walk[walk.length - 1].members.map(k => ops[k]);
  return { ops, order: ops.length, maxResidual, ...classifyOperations(group, tol / meanEdge(A), { closed: true, A }) };
}

/**
 * The result when nothing can be analysed: no basis, or no operation — not even the
 * identity — survives (a lattice with a non-finite entry, or a singular one, rejects every
 * candidate in Step 7). Never "P1 No. 1": that would name a structure never seen.
 */
export const UNDETERMINED = Object.freeze({
  centering: null, pointGroup: '—', spaceGroup: 'undetermined', spaceGroupNumber: null,
  nSpace: 0, nPoint: 0, nTrans: 0, setting: null,
});

// Every candidate operation of (A, basis) within `tol`, with its residual (Steps 7–9).
function detectOperations(A, basis, tol, latticeTol) {
  const pointOps = latticeCandidates(A, Math.min(latticeTol, tol));
  const byEl = new Map();
  for (const s of basis) { if (!byEl.has(s.el)) byEl.set(s.el, []); byEl.get(s.el).push(s); }

  const ops = [];
  let maxResidual = 0;
  // Use the rarest element for candidate translations (fewest partners → fastest).
  let refEl = basis[0].el;
  for (const [el, arr] of byEl) if (arr.length < byEl.get(refEl).length) refEl = el;
  const refAtom = byEl.get(refEl)[0];
  const index = partnerIndex(byEl, A, 2 * tol);

  for (const { R, strain } of pointOps) {
    const Ra0 = applyR(R, refAtom.frac);
    const tSeen = [];
    for (const cand of byEl.get(refEl)) {
      const seed = [wrap01(cand.frac[0] - Ra0[0]), wrap01(cand.frac[1] - Ra0[1]), wrap01(cand.frac[2] - Ra0[2])];
      if (tSeen.some(u => cartDist(u, seed, A) < tol)) continue;   // an accepted op already covers this seed
      const op = refineOperation(R, seed, basis, index, A, tol);
      if (!op) continue;
      if (tSeen.some(u => cartDist(u, op.t, A) < tol)) continue;   // same op reached from another seed
      tSeen.push(op.t);
      const residual = Math.max(op.residual, strain);
      ops.push({ R, t: op.t, residual });
      if (residual > maxResidual) maxResidual = residual;
    }
  }
  return { ops, maxResidual };
}

/* ── group closure ───────────────────────────────────────────────────────────
 * The operations kept at a residual threshold are only a CANDIDATE set: with noise
 * every operation of the true group has its own residual, so a threshold keeps an
 * arbitrary subset of it, and a subset of the right size is not a group. Only a set
 * closed under composition (modulo lattice translations) is ever classified. */

const rotKey = (R) => R.flat().join(',');

/** {Ra|ta}·{Rb|tb} = {Ra·Rb | Ra·tb + ta}, translation reduced mod 1. */
function composeOps(a, b) {
  const R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++)
    R[i][j] = a.R[i][0] * b.R[0][j] + a.R[i][1] * b.R[1][j] + a.R[i][2] * b.R[2][j];
  const t = applyR(a.R, b.t);
  return { R, t: [wrap01(t[0] + a.t[0]), wrap01(t[1] + a.t[1]), wrap01(t[2] + a.t[2])] };
}

// Slack (Å) on the product match, for round-off in exact structures.
const PRODUCT_SLACK = 1e-5;
// Above this many candidate operations the elimination pass (quadratic per threshold)
// is skipped and the groups come from growth alone.
const ELIMINATION_MAX_OPS = 256;

// Canonical sort key of an operation, so ties in residual are broken the same way
// whatever order the operations were found in.
const opKey = (o) => `${rotKey(o.R)}|${o.t.map(v => Math.round(wrap01(v) * 1e4) % 1e4).join(',')}`;

/**
 * Lazily evaluated multiplication table of the detected operations: product(i, j) is the
 * index of the operation equal to ops[i]·ops[j], or −1. Each operation maps every site to
 * within its residual ρ of a same-element site, so the product maps every site to within
 * ρ_i + ρ_j and differs from a same-rotation operation k by at most ρ_i + ρ_j + ρ_k when
 * both send a site to the same partner. Detected translations of one rotation are at
 * least the pass tolerance apart (Step 8 dedup), so the nearest same-rotation operation
 * is the only candidate.
 */
function productTable(ops, A) {
  const n = ops.length;
  const byR = new Map();
  ops.forEach((o, k) => { const key = rotKey(o.R); if (!byR.has(key)) byR.set(key, []); byR.get(key).push(k); });
  const memo = new Map();
  const product = (i, j) => {
    const code = i * n + j;
    const hit = memo.get(code);
    if (hit !== undefined) return hit;
    const p = composeOps(ops[i], ops[j]);
    let best = -1, bestD = Infinity;
    for (const k of byR.get(rotKey(p.R)) || []) { const d = cartDist(p.t, ops[k].t, A); if (d < bestD) { bestD = d; best = k; } }
    const out = best >= 0 && bestD <= ops[i].residual + ops[j].residual + ops[best].residual + PRODUCT_SLACK ? best : -1;
    memo.set(code, out);
    return out;
  };
  return { n, product };
}

/**
 * A closed set of operations grown by generators. `extend(x, alive)` replaces the set by
 * the group generated by it and x, if every element of that group is alive; otherwise it
 * leaves the set unchanged and returns the operation that blocked it (−1 when a product
 * does not exist among the detected operations at all, so x can never join).
 *
 * Closure by generators: starting from the current set (closed under the old generators),
 * add H·x, then right-multiply every new element by every generator. A finite set that
 * holds the identity and is closed under right multiplication by the generators is the
 * group they generate, so no all-pairs check is needed.
 */
function growingGroup(product, n, identity) {
  const inH = new Uint8Array(n);
  const members = [identity];
  const gens = [];
  inH[identity] = 1;
  const mark = new Uint8Array(n);
  return {
    members, gens, inH,
    extend(x, alive) {
      if (inH[x]) return { ok: true };
      const added = [x];
      mark[x] = 1;
      const all = [...gens, x];
      let blocker = null;
      const push = (p) => {
        if (p < 0 || !alive[p]) { blocker = p; return false; }
        if (!inH[p] && !mark[p]) { mark[p] = 1; added.push(p); }
        return true;
      };
      let ok = true;
      for (let i = 0; ok && i < members.length; i++) ok = push(product(members[i], x));
      for (let i = 0; ok && i < added.length; i++) {
        for (let g = 0; ok && g < all.length; g++) ok = push(product(added[i], all[g]));
      }
      for (const k of added) mark[k] = 0;
      if (!ok) return { ok: false, blocker };
      for (const k of added) { inH[k] = 1; members.push(k); }
      gens.push(x);
      return { ok: true };
    },
  };
}

/**
 * Largest closed subset of the operations flagged in `alive`, by elimination: while some
 * product of two kept operations is missing from the kept set, drop the operation with
 * the WORST residual among those taking part in such a product. The identity never takes
 * part in one (e·x = x). Quadratic in the candidate count, so used only for small sets
 * (ELIMINATION_MAX_OPS); it finds groups that growth from a smaller group cannot reach
 * because they do not contain it.
 */
function closeUnder(product, n, ops, alive, keys) {
  const live = Uint8Array.from(alive);
  const idx = [];
  for (let k = 0; k < n; k++) if (live[k]) idx.push(k);
  const bad = new Int32Array(n);
  const into = new Map();                 // target k → pairs [i, j] whose product it is
  const defective = (i, j) => { const k = product(i, j); return k < 0 || !live[k]; };
  for (const i of idx) {
    for (const j of idx) {
      const k = product(i, j);
      if (k < 0 || !live[k]) { bad[i]++; bad[j]++; } else {
        if (!into.has(k)) into.set(k, []);
        into.get(k).push(i, j);
      }
    }
  }
  for (;;) {
    let worst = -1;
    for (const k of idx) {
      if (!live[k] || bad[k] <= 0) continue;
      if (worst < 0 || ops[k].residual > ops[worst].residual
        || (ops[k].residual === ops[worst].residual && keys[k] > keys[worst])) worst = k;
    }
    if (worst < 0) break;
    // Pairs involving `worst` disappear with it; pairs whose product it was become defective.
    for (const j of idx) {
      if (!live[j] || j === worst) continue;
      if (defective(worst, j)) bad[j]--;
      if (defective(j, worst)) bad[j]--;
    }
    live[worst] = 0;
    bad[worst] = 0;
    const pairs = into.get(worst) || [];
    for (let p = 0; p < pairs.length; p += 2) {
      const i = pairs[p], j = pairs[p + 1];
      if (live[i] && live[j]) { bad[i]++; bad[j]++; }
    }
  }
  return live;
}

/**
 * Walk the distinct residual thresholds ≤ `limit`, tight → loose, and return the closed
 * group holding at each. The group only ever grows: at each threshold the newly admitted
 * operations (and any that were waiting on them) are offered, in order of residual, to
 * the current group, and kept when the group they generate with it is entirely admitted.
 * For small candidate sets the largest closed subset found by elimination (closeUnder)
 * replaces the grown group when it is larger. A group closed at a tighter threshold is
 * still closed at a looser one, so the operation count never falls.
 *
 * @returns {{ r:number, members:number[] }[]}  members = indices into `ops`, same array
 *   object for consecutive thresholds holding the same group.
 */
function groupsByThreshold(ops, A, limit) {
  const n = ops.length;
  const order = ops.map((_, k) => k).filter(k => ops[k].residual <= limit + 1e-9);
  if (!order.length) return [];
  const keys = ops.map(opKey);
  order.sort((a, b) => ops[a].residual - ops[b].residual || (keys[a] < keys[b] ? -1 : keys[a] > keys[b] ? 1 : 0));
  // The identity: R = I with the translation nearest 0.
  let identity = -1;
  for (const k of order) if (isIdentityR(ops[k].R) && (identity < 0 || cartDist(ops[k].t, [0, 0, 0], A) < cartDist(ops[identity].t, [0, 0, 0], A))) identity = k;
  if (identity < 0) return [];
  const { product } = productTable(ops, A);
  const alive = new Uint8Array(n);
  const rank = new Int32Array(n);
  order.forEach((k, i) => { rank[k] = i; });

  let group = growingGroup(product, n, identity);
  let pending = [];                        // alive, not in the group, not blocked
  const waiting = new Map();               // blocker → ops waiting for it to be admitted
  const out = [];
  let members = group.members.slice();

  // Offer the pending operations to the group, best residual first. One that fails
  // waits for the operation that blocked it; the block stays valid while the group only
  // grows (the failing product involves elements that stay in it), and a product that
  // does not exist at all blocks for good.
  const offer = () => {
    pending.sort((a, b) => rank[a] - rank[b]);
    for (const x of pending) {
      if (group.inH[x]) continue;
      const res = group.extend(x, alive);
      if (!res.ok && res.blocker >= 0) {
        if (!waiting.has(res.blocker)) waiting.set(res.blocker, []);
        waiting.get(res.blocker).push(x);
      }
    }
    pending = [];
  };

  let i = 0;
  while (i < order.length) {
    const r = ops[order[i]].residual;
    while (i < order.length && ops[order[i]].residual <= r + 1e-9) {
      const k = order[i++];
      alive[k] = 1;
      if (!group.inH[k]) pending.push(k);
      const released = waiting.get(k);
      if (released) { pending.push(...released); waiting.delete(k); }
    }
    offer();
    if (n <= ELIMINATION_MAX_OPS) {
      const live = closeUnder(product, n, ops, alive, keys);
      let size = 0;
      for (let k = 0; k < n; k++) size += live[k];
      if (size > group.members.length) {
        // Switch to the larger group: rebuild its generators, and re-offer everything else.
        group = growingGroup(product, n, identity);
        for (const k of order) if (live[k] && !group.inH[k]) group.extend(k, live);
        waiting.clear();
        pending = order.filter(k => alive[k] && !group.inH[k]);
        offer();
      }
    }
    if (group.members.length !== members.length) members = group.members.slice();
    out.push({ r, members });
  }
  return out;
}

/**
 * Whether an operation set is closed under composition modulo lattice translations,
 * each product matched to a same-rotation member within `ttol` in every fractional
 * component. For callers holding an arbitrary set; the ladder and headline build their
 * sets closed (groupsByThreshold) and skip this.
 */
export function isClosedSet(ops, ttol = 1e-3) {
  const byR = new Map();
  for (const o of ops) { const key = rotKey(o.R); if (!byR.has(key)) byR.set(key, []); byR.get(key).push(o); }
  const near = (u, v) => Math.abs(nearestInt(u[0] - v[0])) <= ttol && Math.abs(nearestInt(u[1] - v[1])) <= ttol
    && Math.abs(nearestInt(u[2] - v[2])) <= ttol;
  for (const a of ops) {
    for (const b of ops) {
      const p = composeOps(a, b);
      const same = byR.get(rotKey(p.R));
      if (!same || !same.some(o => near(o.t, p.t))) return false;
    }
  }
  return true;
}

const isIdentityR = (R) => R[0][0] === 1 && R[1][1] === 1 && R[2][2] === 1
  && R[0][1] === 0 && R[0][2] === 0 && R[1][0] === 0 && R[1][2] === 0 && R[2][0] === 0 && R[2][1] === 0;

/**
 * Classify a set of operations {R,t} into point group + centering → H–M symbol.
 * Only a group is named: pass `{ closed: true }` when the set was built closed
 * (groupsByThreshold); otherwise closure is checked here, products matched within
 * 3·tolFrac per fractional component. `A` (the cell's lattice rows, Å) lets the naming
 * search cells other than the given one (spaceGroupSymbol.js → derivedBases).
 *
 * `centering` is the centering of the GIVEN cell read from all its pure translations,
 * or null when they are not exactly a Bravais centering (a supercell of the true cell).
 * `setting` is the standard cell the symbol belongs to (null when there is none, or when
 * Wyckoff positions cannot be placed — see spaceGroupHM).
 */
export function classifyOperations(ops, tolFrac = 0.02, { closed = false, A = null } = {}) {
  const rotMap = new Map();
  const transSeen = new Set();          // distinct pure (identity-rotation) translations
  for (const { R, t } of ops) {
    const key = R.flat().join(',');
    if (!rotMap.has(key)) rotMap.set(key, R);
    // Folded mod 1 before keying: 0.9997 and 0 are the same translation.
    if (isIdentityR(R)) transSeen.add(t.map(x => Math.round(wrap01(x) * 1000) % 1000).join(','));
  }
  const centering = centeringOfOps(ops);
  const pointGroup = pointGroupOf([...rotMap.values()]);
  const base = { centering, pointGroup, nSpace: ops.length, nPoint: rotMap.size, nTrans: transSeen.size };
  const group = isValidGroup(base) && (closed || isClosedSet(ops, Math.max(3 * tolFrac, 1e-6)));
  // Naming reads every operation's screw/glide part, so it is the expensive step and
  // meaningless for a set that is not a group.
  // Lattice symmetries that hold at least as well as the group's own operations.
  const groupResidual = ops.reduce((m, o) => (o.residual > m ? o.residual : m), 0);
  const sg = group ? spaceGroupHM(centering, pointGroup, ops, A, Math.max(groupResidual, 1e-3)) : { symbol: 'not a group', number: null, setting: null };
  return { ...base, spaceGroup: sg.symbol, spaceGroupNumber: sg.number, setting: sg.setting ?? null };
}

const POINT_GROUP_ORDER = {
  '1': 1, '-1': 2, '2': 2, 'm': 2, '2/m': 4, '222': 4, 'mm2': 4, 'mmm': 8,
  '4': 4, '-4': 4, '4/m': 8, '422': 8, '4mm': 8, '-42m': 8, '4/mmm': 16,
  '3': 3, '-3': 6, '32': 6, '3m': 6, '-3m': 12, '6': 6, '-6': 6, '6/m': 12,
  '622': 12, '6mm': 12, '-6m2': 12, '6/mmm': 24,
  '23': 12, 'm-3': 24, '432': 24, '-43m': 24, 'm-3m': 48,
};
// Cheap pre-filter for a group (closure itself is checked separately — isClosedSet, or
// by construction in groupsByThreshold): the op count equals point-group order × the
// actual number of pure translations, as it must for any closed set (a tiled supercell
// has extra translations, so they are counted rather than assumed). There is no
// centering allow-list: a subgroup keeps its parent cell's centering (R3m or I2/a in an
// F-cubic cell), which is only a matter of which cell names it (spaceGroupHM).
function isValidGroup(cls) {
  const pg = POINT_GROUP_ORDER[cls.pointGroup] || 0;
  return pg !== 0 && cls.nSpace === pg * (cls.nTrans || 1);
}

/**
 * Symmetry-vs-tolerance ladder in ONE detection pass. Detect at the loosest
 * tolerance (all candidate ops with their residuals), then threshold: an op is a
 * candidate at tolerance t iff its residual ≤ t, and the rung is the largest CLOSED
 * group among the candidates (groupsByThreshold). Distinct thresholds → the rungs,
 * each a group over a tolerance range (merged when the symbol repeats). Monotonic:
 * looser tol ⇒ equal or more operations. The first threshold is 0 (the identity).
 *
 * @returns {{from:number, to:number, spaceGroup:string, spaceGroupNumber:number|null,
 *            pointGroup:string, nSpace:number}[]} bricks, tight→loose.
 */
export function symmetryLadder(A, basis, tolMax = 1.5, latticeTol = tolMax) {
  if (!basis || !basis.length) return [];
  const { ops } = detectOperations(A, basis, tolMax, latticeTol);
  if (!ops.length) return [];
  const tolFrac = tolMax / meanEdge(A);
  const walk = groupsByThreshold(ops, A, tolMax);
  const bricks = [];
  let lastMembers = null, cls = null;
  for (let i = 0; i < walk.length; i++) {
    const { r, members } = walk[i];
    const to = i + 1 < walk.length ? walk[i + 1].r : tolMax;  // the last group holds for all looser tol
    if (members !== lastMembers) {
      cls = classifyOperations(members.map(k => ops[k]), tolFrac, { closed: true, A });
      lastMembers = members;
    }
    const last = bricks[bricks.length - 1];
    if (last && last.spaceGroup === cls.spaceGroup) { last.to = to; last.nSpace = cls.nSpace; continue; }
    bricks.push({ from: r, to, spaceGroup: cls.spaceGroup, spaceGroupNumber: cls.spaceGroupNumber, pointGroup: cls.pointGroup, nSpace: cls.nSpace });
  }
  return bricks;
}

/**
 * The space group holding at cartesian tolerance `tol` (Å), with its operations: the
 * same closed group the ladder shows at `tol` (groupsByThreshold over one detection
 * pass at `tol`). `maxResidual` is the worst residual of the operations returned.
 *
 * @returns {{ centering, pointGroup, spaceGroup, spaceGroupNumber, nSpace, nPoint,
 *             maxResidual, ops:{R,t,residual}[] }}
 */
export function spaceGroupAtTolerance(A, basis, tol = 0.2, latticeTol = Math.max(tol, 1e-3)) {
  const empty = { ...UNDETERMINED, maxResidual: Number.NaN, ops: [] };
  if (!basis || !basis.length) return empty;
  const { ops: all } = detectOperations(A, basis, Math.max(tol, 1e-3), latticeTol);
  const walk = groupsByThreshold(all, A, tol);
  if (!walk.length) return empty;
  const ops = walk[walk.length - 1].members.map(k => all[k]);
  const tolFrac = Math.max(tol, 1e-6) / meanEdge(A);
  const maxResidual = ops.reduce((m, o) => Math.max(m, o.residual), 0);
  return { ...classifyOperations(ops, tolFrac, { closed: true, A }), maxResidual, ops };
}

/* ── Space-group identification ──────────────────────────────────────────────
 * Classify the detected operations into a point group + centering → Hermann–
 * Mauguin symbol. Rotation TYPE is basis-independent (det & trace are similarity
 * invariants): proper (det +1) trace 3,2,1,0,−1 → 1,6,4,3,2-fold; improper
 * (det −1) trace −3,1,−1,0,−2 → inversion, mirror m, −4, −3, −6. The point group
 * is then fixed by (crystal system, order, has-inversion, which proper folds).
 * The full H–M symbol comes from spaceGroupSymbol.js, which reads each operation's
 * screw/glide component so non-symmorphic groups are named as themselves (Pnma,
 * I4/mcm, Fd-3m) rather than as their symmorphic parent. */

export function classifyRotation(R) {
  const d = det3(R), t = R[0][0] + R[1][1] + R[2][2];
  if (d === 1) return t === 3 ? '1' : t === 2 ? '6' : t === 1 ? '4' : t === 0 ? '3' : '2';
  return t === -3 ? '-1' : t === 1 ? 'm' : t === -1 ? '-4' : t === 0 ? '-3' : '-6';
}

// The distinct rotation parts of a space group → its point-group H–M symbol.
// The crystal CLASS is derived from the rotation content itself (not the lattice
// metric), so a structure whose symmetry is a proper subgroup of its lattice — the
// generic case along the tolerance ladder — is classified correctly.
export function pointGroupOf(rotations) {
  const h = { '1': 0, '2': 0, '3': 0, '4': 0, '6': 0, '-1': 0, 'm': 0, '-3': 0, '-4': 0, '-6': 0 };
  for (const R of rotations) h[classifyRotation(R)]++;
  const order = rotations.length, inv = h['-1'] > 0, nm = h.m;
  if (h['3'] >= 8) {                                     // cubic: 4 three-fold axes
    if (order === 12) return '23';
    if (order === 48) return 'm-3m';
    return inv ? 'm-3' : (h['4'] > 0 ? '432' : '-43m');
  }
  if (h['6'] > 0 || h['-6'] > 0) {                       // hexagonal
    if (order === 24) return '6/mmm';
    if (order === 12) return inv ? '6/m' : (h['6'] > 0 ? (nm > 0 ? '6mm' : '622') : '-6m2');
    return h['6'] > 0 ? '6' : '-6';
  }
  if (h['4'] > 0 || h['-4'] > 0) {                       // tetragonal
    if (order === 16) return '4/mmm';
    if (order === 8) return inv ? '4/m' : (h['4'] > 0 ? (nm > 0 ? '4mm' : '422') : '-42m');
    return h['4'] > 0 ? '4' : '-4';
  }
  if (h['3'] > 0 || h['-3'] > 0) {                       // trigonal
    if (order === 12) return '-3m';
    if (order === 6) return inv ? '-3' : (nm > 0 ? '3m' : '32');
    return '3';
  }
  if (h['2'] > 0 || nm > 0) {                            // ortho / mono
    if (order === 8) return 'mmm';
    if (order === 4) return inv ? '2/m' : (h['2'] >= 3 ? '222' : 'mm2');
    return h['2'] > 0 ? '2' : 'm';
  }
  return inv ? '-1' : '1';                               // triclinic
}

/** Label shown when a group cannot be named reliably: its crystal class, no number. */
export const classLabel = (pointGroup) => `${pointGroup} class`;
/** Label for a verified group that may be a subgroup of the true one (see spaceGroupHM). */
export const lowerBoundLabel = (symbol) => `≥ ${symbol}`;

/**
 * Hermann–Mauguin symbol + ITA number for a detected group.
 *
 * With the operations in hand this reads their screw and glide components and names the
 * actual space group (Pnma, I4/mcm, Fd-3m), accepting only a tabulated symbol of the
 * detected crystal class and centering (hmSymbolInStandardSetting). Triclinic groups
 * are P1 or P-1 whatever cell describes them. Anything that cannot be named that way is
 * reported as its crystal class with no number — never as the symmorphic `centering +
 * point group` string, which spells a real but different group for most classes.
 *
 * Whatever the branch, when the operations' full translation lattice has rotations that
 * are not integer matrices of the given cell (a supercell of the true cell: the cubic
 * 3-folds of a 2×2×1 cell), those were never tested, and the result — symbol, P1/P-1 or
 * crystal class — is only a lower bound: '≥ …', no number, no setting. `latticeTol` (Å)
 * is the strain up to which a rotation of that lattice counts (the group's worst residual).
 *
 * @returns {{ symbol:string, number:number|null, standard:boolean, setting:object|null }}
 */
export function spaceGroupHM(centering, pointGroup, ops, A = null, latticeTol = 1e-3) {
  const holohedry = (As, tol) => latticeCandidates(As, tol).map(c => c.R);
  if (pointGroup === '1' || pointGroup === '-1') {
    // P1 / P-1 whatever cell describes them; Wyckoff positions only in a primitive cell.
    const symbol = pointGroup === '1' ? 'P1' : 'P-1';
    if (A && ops?.length && !latticeFullyTested(ops, A, holohedry, latticeTol)) {
      return { symbol: lowerBoundLabel(symbol), number: null, standard: false, setting: null };
    }
    const primitive = centering === 'P';
    const I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    return { symbol, number: spaceGroupNumber(symbol), standard: true, setting: primitive ? { Q: I, Qinv: I, translations: [[0, 0, 0]], ratio: 1 } : null };
  }
  if (!ops || !ops.length) return { symbol: classLabel(pointGroup), number: null, standard: false, setting: null };
  const found = hmSymbolInStandardSetting(ops, centering, pointGroup, pointGroupOfSymbol, { A, holohedry, latticeTol });
  if (!found.symbol) {
    const untested = A && !latticeFullyTested(ops, A, holohedry, latticeTol);
    return { symbol: untested ? lowerBoundLabel(classLabel(pointGroup)) : classLabel(pointGroup), number: null, standard: false, setting: null };
  }
  const symbol = canonicalSymbol(found.symbol) ?? found.symbol;
  // A group named in a cell whose lattice symmetries the given cell could not all test is
  // only a lower bound: shown as such, with no number and no Wyckoff letters.
  if (found.complete === false) return { symbol: lowerBoundLabel(symbol), number: null, standard: false, setting: null };
  return { symbol, number: spaceGroupNumber(found.symbol), standard: true, setting: found.setting };
}

/**
 * Partition the basis into symmetry orbits under the operations `ops`: two sites
 * are in the same orbit if some {R|t} maps one onto the other (same element,
 * within `tol` Å). Orbits = the Wyckoff structure of the (average) crystal — the
 * sets of sites the detected symmetry says are equivalent.
 *
 * Each orbit also carries its multiplicity (size), site-symmetry point group (the
 * stabilizer of a representative — robust, no tables), and a representative frac.
 *
 * @returns {{ index:number[], element:string, size:number, site:string, rep:number[] }[]}
 *   orbits, largest first.
 */
export function siteOrbits(A, basis, ops, tol = 0.1) {
  const n = basis.length;
  const parent = Array.from({ length: n }, (_, i) => i);
  const find = (x) => { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
  const union = (a, b) => { const ra = find(a), rb = find(b); if (ra !== rb) parent[ra] = rb; };

  const byEl = new Map();
  basis.forEach((s, i) => { if (!byEl.has(s.el)) byEl.set(s.el, []); byEl.get(s.el).push({ frac: s.frac, i }); });
  const index = partnerIndex(byEl, A, tol);

  for (const { R, t } of ops) {
    for (let i = 0; i < n; i++) {
      const img = applyR(R, basis[i].frac);
      img[0] = wrap01(img[0] + t[0]); img[1] = wrap01(img[1] + t[1]); img[2] = wrap01(img[2] + t[2]);
      let best = -1, bestD = tol;
      for (const o of index.near(basis[i].el, img)) { const d = cartDist(img, o.frac, A); if (d < bestD) { bestD = d; best = o.i; } }
      if (best >= 0) union(i, best);
    }
  }
  const groups = new Map();
  for (let i = 0; i < n; i++) { const r = find(i); if (!groups.has(r)) groups.set(r, []); groups.get(r).push(i); }
  return [...groups.values()]
    .map(index => {
      const rep = basis[index[0]];
      // Site symmetry = rotation parts of the ops that FIX the representative point.
      const stab = []; const seen = new Set();
      for (const { R, t } of ops) {
        const img = applyR(R, rep.frac);
        img[0] = wrap01(img[0] + t[0]); img[1] = wrap01(img[1] + t[1]); img[2] = wrap01(img[2] + t[2]);
        if (cartDist(img, rep.frac, A) < tol) { const k = R.flat().join(','); if (!seen.has(k)) { seen.add(k); stab.push(R); } }
      }
      return { index, element: rep.el, size: index.length, site: pointGroupOf(stab), rep: rep.frac };
    })
    .sort((a, b) => b.size - a.size);
}

// Wyckoff letters live in wyckoff.js, which matches an orbit's multiplicity and site
// symmetry against the tabulated positions of the detected group.
