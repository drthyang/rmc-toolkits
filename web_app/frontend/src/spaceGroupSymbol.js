// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/spaceGroupSymbol.js
//
// Hermann–Mauguin space-group symbol construction from a detected operation set.
//
// symmetry.js finds the operations {R|t} of a structure, but named the result by
// concatenating the centering letter with the point group — correct only for the
// 73 symmorphic groups. A structure in Pnma, I4/mcm or Fd-3m has exactly the right
// operations detected yet came back labelled Pmmm, I4/mmm, Fm-3m, because the
// translation parts were dropped at naming time. This module keeps them.
//
// Each operation is split into an INTRINSIC translation (the screw/glide component)
// and a location part. The intrinsic part is the projection of t onto the invariant
// subspace of R:
//
//     t_intrinsic = (1/n) · Σ_{k=0..n-1} Rᵏ·t ,     n = order of R (Rⁿ = I)
//
// — the axis for a rotation, the plane for a reflection. It is independent of where
// the origin sits, which is what makes it a reliable element label: a non-zero part
// along the axis makes a rotation a SCREW (n_m), a non-zero part lying in the mirror
// plane makes a reflection a GLIDE (a, b, c, n, d, e). Since t is only defined modulo
// a lattice vector, t_intrinsic is well defined modulo a lattice vector too (the sum
// telescopes to a full lattice translation), so it is always reduced mod 1.
//
// The symbol is then assembled positionally: every crystal system has an ordered list
// of symmetry DIRECTIONS (blickrichtungen) and each position of the H–M symbol names
// the element found along that direction. Axis direction for a rotation; plane NORMAL
// for a reflection.
//
// Fractional coords are COLUMN vectors, matching symmetry.js: x' = R·x + t.

/* ── small integer/matrix helpers (kept local so this module has no imports) ─── */

const IDENTITY = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

export function det3i(m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

function matMul(A, B) {
  const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++)
    M[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
  return M;
}

function matVec(R, v) {
  return [
    R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
    R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
    R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2],
  ];
}

const sameMat = (A, B) => A.every((row, i) => row.every((v, j) => v === B[i][j]));
const trace = (R) => R[0][0] + R[1][1] + R[2][2];
const wrap1 = (x) => { const y = x - Math.floor(x); return y > 1 - 1e-6 ? 0 : y; };

/** Order n of an integer rotation matrix: the smallest n ≥ 1 with Rⁿ = I (n ≤ 6). */
export function rotationOrder(R) {
  let M = R;
  for (let n = 1; n <= 6; n++) {
    if (sameMat(M, IDENTITY)) return n;
    M = matMul(M, R);
  }
  return 0;                                     // not a crystallographic rotation
}

/* ── intrinsic (screw / glide) translation ──────────────────────────────────── */

/**
 * Intrinsic translation of {R|t}: (1/n)·Σ Rᵏ·t, reduced mod 1. Zero for a pure
 * rotation/reflection, (m/n) along the axis for a screw, the glide vector for a glide.
 */
export function intrinsicTranslation(R, t) {
  const n = rotationOrder(R);
  if (!n) return [0, 0, 0];
  const sum = [0, 0, 0];
  let M = IDENTITY;
  for (let k = 0; k < n; k++) {
    const v = matVec(M, t);
    sum[0] += v[0]; sum[1] += v[1]; sum[2] += v[2];
    M = matMul(M, R);
  }
  return sum.map((v) => wrap1(v / n));
}

/* ── characteristic direction ───────────────────────────────────────────────── */

const gcd2 = (a, b) => { a = Math.abs(a); b = Math.abs(b); while (b) { [a, b] = [b, a % b]; } return a; };

/** Reduce an integer direction by its gcd and fix the sign (first non-zero positive). */
export function normalizeDirection(d) {
  const g = gcd2(gcd2(d[0], d[1]), d[2]) || 1;
  const r = d.map((v) => v / g);
  const lead = r.find((v) => v !== 0) || 1;
  // `|| 0` keeps a negated zero from surfacing as -0.
  return lead < 0 ? r.map((v) => -v || 0) : r;
}

// Candidate integer directions, small components — covers every crystallographic
// blickrichtung including the hexagonal <2 1 0> family.
const CANDIDATE_DIRECTIONS = (() => {
  const out = [];
  const seen = new Set();
  for (let x = -2; x <= 2; x++) for (let y = -2; y <= 2; y++) for (let z = -2; z <= 2; z++) {
    if (!x && !y && !z) continue;
    const d = normalizeDirection([x, y, z]);
    const k = d.join(',');
    if (!seen.has(k)) { seen.add(k); out.push(d); }
  }
  return out;
})();

/**
 * The direction that labels an operation's position in the H–M symbol: the rotation
 * axis for a proper rotation (R·d = d), the plane normal for a reflection and the
 * axis for a rotoinversion (R·d = −d). Null for the identity and for inversion,
 * which have no single direction.
 */
export function characteristicDirection(R) {
  if (sameMat(R, IDENTITY)) return null;
  if (R[0][0] === -1 && R[1][1] === -1 && R[2][2] === -1
    && !R[0][1] && !R[0][2] && !R[1][0] && !R[1][2] && !R[2][0] && !R[2][1]) return null;
  const want = det3i(R) === 1 ? 1 : -1;
  for (const d of CANDIDATE_DIRECTIONS) {
    const v = matVec(R, d);
    if (v[0] === want * d[0] && v[1] === want * d[1] && v[2] === want * d[2]) return d;
  }
  return null;
}

/* ── element classification ─────────────────────────────────────────────────── */

const QUARTERS = [0, 0.25, 0.5, 0.75];
const snapQuarter = (x) => {
  const y = wrap1(x);
  let best = 0, bd = Infinity;
  for (const q of QUARTERS) { const d = Math.min(Math.abs(y - q), Math.abs(y - q - 1)); if (d < bd) { bd = d; best = q; } }
  return best;
};

/**
 * Glide letter for a reflection whose intrinsic translation is `ti` (fractional).
 * m (none) · a/b/c (half a lattice vector) · n (diagonal) · d (quarter, diamond).
 */
export function glideLetter(ti) {
  const q = ti.map(snapQuarter);
  if (q.some((v) => v === 0.25 || v === 0.75)) return 'd';
  const halves = q.map((v, i) => (v === 0.5 ? i : -1)).filter((i) => i >= 0);
  if (halves.length === 0) return 'm';
  if (halves.length === 1) return ['a', 'b', 'c'][halves[0]];
  return 'n';
}

/**
 * Screw index m of a proper rotation of order n about `axis`: the intrinsic
 * translation is (m/n) of the shortest lattice vector along the axis. Returned in
 * 0…n−1, measured for a right-handed rotation about +axis so that enantiomorphic
 * pairs (4₁/4₃, 3₁/3₂, 6₁/6₅, 6₂/6₄) come out distinct.
 */
export function screwIndex(R, t, axis, n) {
  const ti = intrinsicTranslation(R, t);
  const dd = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
  const frac = (ti[0] * axis[0] + ti[1] * axis[1] + ti[2] * axis[2]) / dd;
  let m = Math.round(frac * n) % n;
  if (m < 0) m += n;
  if (m === 0 || n < 3) return m;
  // Handedness: for n ≥ 3 the rotation by +2π/n and its inverse carry m and n−m.
  // Report the index of the +2π/n (right-handed about +axis) member.
  return rotationSense(R, axis) >= 0 ? m : (n - m) % n;
}

/**
 * Sign of the rotation sense of a proper rotation about `axis`, from the handedness
 * of {axis, v, R·v}. Basis-independent up to the handedness of the lattice basis.
 */
export function rotationSense(R, axis) {
  for (const v of CANDIDATE_DIRECTIONS) {
    // pick a v not parallel to the axis
    const cx = axis[1] * v[2] - axis[2] * v[1];
    const cy = axis[2] * v[0] - axis[0] * v[2];
    const cz = axis[0] * v[1] - axis[1] * v[0];
    if (!cx && !cy && !cz) continue;
    const w = matVec(R, v);
    const d = det3i([axis, v, w]);
    if (d !== 0) return Math.sign(d);
  }
  return 1;
}

/**
 * Classify one operation {R|t} into a symmetry element.
 * @returns {{kind, order, direction, label}|null} null for the identity, a pure
 *   lattice/centering translation and for inversion (none occupy a symbol position).
 */
export function classifyElement(R, t) {
  const dir = characteristicDirection(R);
  if (!dir) return null;
  const n = rotationOrder(R);
  const proper = det3i(R) === 1;
  if (proper) {
    const m = screwIndex(R, t, dir, n);
    return { kind: m ? 'screw' : 'rotation', order: n, direction: dir, label: m ? `${n}_${m}` : `${n}` };
  }
  if (trace(R) === 1) {                         // reflection: n = 2, invariant plane
    const letter = glideLetter(intrinsicTranslation(R, t));
    return { kind: letter === 'm' ? 'mirror' : 'glide', order: 2, direction: dir, label: letter };
  }
  // rotoinversion −4 (n = 4), −3 and −6 (n = 6)
  const bar = n === 4 ? 4 : (trace(R) === 0 ? 3 : 6);
  return { kind: 'rotoinversion', order: bar, direction: dir, label: `-${bar}` };
}

/* ── symbol assembly ────────────────────────────────────────────────────────── */

const dirKey = (d) => normalizeDirection(d).join(',');

// Ordered symmetry directions per crystal system. Each position is the family of
// equivalent directions whose elements share one slot of the H–M symbol.
const SYSTEM_DIRECTIONS = {
  triclinic: [],
  monoclinic: [[[0, 1, 0]]],                                  // re-aimed at the real unique axis
  orthorhombic: [[[1, 0, 0]], [[0, 1, 0]], [[0, 0, 1]]],
  tetragonal: [[[0, 0, 1]], [[1, 0, 0], [0, 1, 0]], [[1, 1, 0], [1, -1, 0]]],
  trigonal: [[[0, 0, 1]], [[1, 0, 0], [0, 1, 0], [1, 1, 0]], [[1, -1, 0], [1, 2, 0], [2, 1, 0]]],
  hexagonal: [[[0, 0, 1]], [[1, 0, 0], [0, 1, 0], [1, 1, 0]], [[1, -1, 0], [1, 2, 0], [2, 1, 0]]],
  cubic: [
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[1, 1, 1], [1, -1, -1], [1, -1, 1], [1, 1, -1]],
    [[1, 1, 0], [1, -1, 0], [0, 1, 1], [0, 1, -1], [1, 0, 1], [1, 0, -1]],
  ],
};

export const POINT_GROUP_SYSTEM = {
  '1': 'triclinic', '-1': 'triclinic',
  '2': 'monoclinic', 'm': 'monoclinic', '2/m': 'monoclinic',
  '222': 'orthorhombic', 'mm2': 'orthorhombic', 'mmm': 'orthorhombic',
  '4': 'tetragonal', '-4': 'tetragonal', '4/m': 'tetragonal', '422': 'tetragonal',
  '4mm': 'tetragonal', '-42m': 'tetragonal', '4/mmm': 'tetragonal',
  '3': 'trigonal', '-3': 'trigonal', '32': 'trigonal', '3m': 'trigonal', '-3m': 'trigonal',
  '6': 'hexagonal', '-6': 'hexagonal', '6/m': 'hexagonal', '622': 'hexagonal',
  '6mm': 'hexagonal', '-6m2': 'hexagonal', '6/mmm': 'hexagonal',
  '23': 'cubic', 'm-3': 'cubic', '432': 'cubic', '-43m': 'cubic', 'm-3m': 'cubic',
};

const PLANE_PRIORITY = ['m', 'e', 'a', 'b', 'c', 'n', 'd'];

/**
 * Plane letters for one symbol position, best first (H–M priority m > e > a > b > c
 * > n > d). More than one is genuinely possible and the priority order is not always
 * the answer: a centred lattice puts several glides in the same plane — I-centring
 * gives I4/mcm both b- and c-glides perpendicular to [100], and ITA writes c. So the
 * options are returned ranked and the caller picks the one the space-group table
 * recognises, rather than committing here.
 */
function planeOptions(letters) {
  if (!letters.length) return [];
  const set = new Set(letters);
  // Two different axial glides sharing a plane are the "e" double glide of ITA 5th
  // ed. (Cmce, Aem2, …) — offered as a candidate, not forced.
  if (['a', 'b', 'c'].filter((l) => set.has(l)).length >= 2) set.add('e');
  return PLANE_PRIORITY.filter((l) => set.has(l));
}

/** Best axis label for one symbol position: highest order, rotation before screw. */
function pickAxis(elements) {
  const proper = elements.filter((e) => e.kind === 'rotation' || e.kind === 'screw');
  const bars = elements.filter((e) => e.kind === 'rotoinversion');
  const maxOrder = proper.reduce((a, e) => Math.max(a, e.order), 0);
  const barOrder = bars.reduce((a, e) => Math.max(a, e.order), 0);
  // −3 is always written in place of 3; −4 / −6 only when no proper rotation of
  // that order is present.
  if (barOrder === 3 && maxOrder <= 3) return { label: '-3', order: 3, isBar: true };
  if (maxOrder >= 4 || (maxOrder === 3 && barOrder !== 6)) {
    const same = proper.filter((e) => e.order === maxOrder);
    const pure = same.find((e) => e.kind === 'rotation');
    return { label: (pure || same[0]).label, order: maxOrder, isBar: false };
  }
  if (barOrder >= 4) return { label: `-${barOrder}`, order: barOrder, isBar: true };
  if (maxOrder >= 2) {
    const same = proper.filter((e) => e.order === maxOrder);
    const pure = same.find((e) => e.kind === 'rotation');
    return { label: (pure || same[0]).label, order: maxOrder, isBar: false };
  }
  return null;
}

/** Axis labels for one position, best first. */
function axisOptions(elements) {
  const pick = pickAxis(elements);
  if (!pick) return [];
  const out = [pick];
  // 4₁/4₃, 3₁/3₂, 6₁/6₅ and 6₂/6₄ are enantiomorphic pairs told apart only by the
  // handedness of the rotation. screwIndex() resolves that, but a left-handed input
  // cell would flip it, so the partner is offered as a fallback candidate.
  const m = /^(\d)_(\d)$/.exec(pick.label);
  if (m) {
    const n = Number(m[1]);
    const k = Number(m[2]);
    const alt = (n - k) % n;
    if (alt && alt !== k) out.push({ ...pick, label: `${n}_${alt}` });
  }
  return out;
}

/**
 * Candidate Hermann–Mauguin symbols for a detected operation set, best first, in the
 * setting the operations are given in.
 *
 * Several positions admit more than one defensible letter (see planeOptions), so this
 * returns a ranked list rather than one answer; `hmSymbolInStandardSetting` picks the
 * first that the space-group table recognises. The first entry is always the one the
 * plain H–M priority rules give.
 *
 * @param {{R:number[][], t:number[]}[]} ops   all operations of the group
 * @param {string} centering  Bravais centering letter (P, A, B, C, I, F, R)
 * @param {string} pointGroup one of the 32 H–M point-group symbols
 * @returns {string[]} ranked candidate symbols
 */
export function hmSymbolCandidates(ops, centering, pointGroup, limit = 96) {
  const system = POINT_GROUP_SYSTEM[pointGroup];
  if (!system || system === 'triclinic') return [`${centering}${pointGroup}`];

  const elements = [];
  for (const { R, t } of ops) {
    const e = classifyElement(R, t);
    if (e) elements.push(e);
  }

  let families = SYSTEM_DIRECTIONS[system];
  if (system === 'monoclinic') {
    // Aim the single position at the unique axis actually present.
    const e = elements.find((x) => x.order === 2);
    families = [[e ? e.direction : [0, 1, 0]]];
  }

  const positions = families.map((family) => {
    const keys = new Set(family.map(dirKey));
    const here = elements.filter((e) => keys.has(dirKey(e.direction)));
    return {
      axes: axisOptions(here),
      planes: planeOptions(here.filter((e) => e.kind === 'mirror' || e.kind === 'glide').map((e) => e.label)),
    };
  });

  // Rendering rules for the SHORT symbol:
  //  · a position with an axis of order > 2 AND a perpendicular plane is written
  //    "axis/plane" (4/m, 6₃/m — and in monoclinic, where the single position always
  //    carries the axis, 2₁/c);
  //  · otherwise the plane wins over a 2-fold axis (Pnma, not P2₁/n2₁/m2₁/a);
  //  · cubic never writes the combined form — Pm-3m, not P4/m-32/m;
  //  · a rotoinversion outranks the plane it implies. −6 IS 3/m, so the mirror
  //    perpendicular to c comes free with it and P-6m2 keeps the −6 rather than
  //    reporting that mirror.
  const render = (axis, plane) => {
    if (!axis && !plane) return '1';
    if (!axis) return plane;
    if (!plane) return axis.label;
    if (system === 'cubic') return plane;
    if (axis.isBar) return axis.label;
    if (system === 'monoclinic') return `${axis.label}/${plane}`;
    if (axis.order > 2) return `${axis.label}/${plane}`;
    return plane;
  };

  const perPosition = positions.map((p) => {
    const out = [];
    for (const a of (p.axes.length ? p.axes : [null])) {
      for (const pl of (p.planes.length ? p.planes : [null])) {
        const s = render(a, pl);
        if (!out.includes(s)) out.push(s);
      }
    }
    return out;
  });

  // Ranked cross-product: earlier positions vary slowest, so the first combination is
  // the all-best one.
  let combos = [[]];
  for (const opts of perPosition) {
    const next = [];
    for (const c of combos) for (const o of opts) if (next.length < limit) next.push([...c, o]);
    combos = next;
  }

  const out = [];
  const push = (s) => { if (s && !out.includes(s)) out.push(s); };
  for (const parts of combos) {
    // A trailing "1" placeholder is dropped (R3m, Pa-3) but an embedded one is
    // significant — P3m1 and P31m are different groups — so offer both spellings.
    const trimmed = [...parts];
    while (trimmed.length && trimmed[trimmed.length - 1] === '1') trimmed.pop();
    push(centering + trimmed.join(''));
    push(centering + parts.join(''));
  }
  return out;
}

/**
 * Best-guess Hermann–Mauguin symbol plus the full ranked candidate list.
 * @returns {{ symbol:string, candidates:string[] }}
 */
export function hmSymbol(ops, centering, pointGroup) {
  const candidates = hmSymbolCandidates(ops, centering, pointGroup);
  return { symbol: candidates[0] ?? `${centering}${pointGroup}`, candidates };
}

/**
 * Whether every symmetry element of `ops` lies along one of the crystal system's
 * symmetry directions.
 *
 * False means the operations are in an axis setting the H–M positions cannot describe,
 * and any symbol built from them is meaningless rather than merely non-standard. The
 * case that matters in practice: a low-symmetry subgroup detected inside a cubic parent
 * cell — a rung part-way up the tolerance ladder — keeps the parent's axes, so its
 * 3-fold runs along ⟨111⟩ where the hexagonal-axes positions expect [001], and the
 * unplaceable elements would silently drop out of the symbol.
 */
export function coversAllElements(ops, pointGroup) {
  const system = POINT_GROUP_SYSTEM[pointGroup];
  if (!system) return false;
  if (system === 'triclinic') return true;

  const elements = [];
  for (const { R, t } of ops) {
    const e = classifyElement(R, t);
    if (e) elements.push(e);
  }
  let families = SYSTEM_DIRECTIONS[system];
  if (system === 'monoclinic') {
    const unique = elements.find((x) => x.order === 2);
    families = [[unique ? unique.direction : [0, 1, 0]]];
  }
  const keys = new Set(families.flat().map(dirKey));
  return elements.every((e) => keys.has(dirKey(e.direction)));
}

/* ── setting search ─────────────────────────────────────────────────────────── */

// The 6 axis permutations. A monoclinic or orthorhombic cell coming out of an RMC
// box need not be in the standard setting (Pbnm rather than Pnma), and the symbol
// is only meaningful once the axes are ordered the standard way.
const PERMUTATIONS = [
  [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
  [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
  [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
  [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
  [[1, 0, 0], [0, 0, 1], [0, 1, 0]],
  [[0, 0, 1], [0, 1, 0], [1, 0, 0]],
];

function inverseIntMatrix(P) {
  const d = det3i(P);
  const c = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
    const m = [];
    for (let a = 0; a < 3; a++) { if (a === i) continue; const row = []; for (let b = 0; b < 3; b++) { if (b === j) continue; row.push(P[a][b]); } m.push(row); }
    c[j][i] = ((i + j) % 2 ? -1 : 1) * (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / d;
  }
  return c;
}

/** Re-express the operations in a permuted axis setting: R' = P⁻¹RP, t' = P⁻¹t. */
export function transformOps(ops, P) {
  const Pi = inverseIntMatrix(P);
  return ops.map(({ R, t }) => ({
    R: matMul(matMul(Pi, R), P).map((row) => row.map((v) => Math.round(v))),
    t: matVec(Pi, t).map(wrap1),
  }));
}

/**
 * H–M symbol for an operation set, searching axis settings until one is recognised.
 *
 * Monoclinic and orthorhombic groups are only named in a standard setting, and an
 * RMC cell is under no obligation to be in one. `isKnown` decides which spellings
 * count as standard (injected so this module stays table-free); when nothing is
 * recognised the identity-setting symbol is returned with `standard: false`.
 *
 * @returns {{ symbol:string, standard:boolean, permutation:number }}
 */
export function hmSymbolInStandardSetting(ops, centering, pointGroup, isKnown) {
  const accept = (s) => (typeof isKnown === 'function' ? isKnown(s) : false);
  const direct = hmSymbolCandidates(ops, centering, pointGroup);
  for (const cand of direct) if (accept(cand)) return { symbol: cand, standard: true, permutation: 0 };

  const system = POINT_GROUP_SYSTEM[pointGroup];
  if (system === 'monoclinic' || system === 'orthorhombic') {
    for (let i = 1; i < PERMUTATIONS.length; i++) {
      const moved = transformOps(ops, PERMUTATIONS[i]);
      for (const cand of hmSymbolCandidates(moved, centering, pointGroup)) {
        if (accept(cand)) return { symbol: cand, standard: true, permutation: i };
      }
    }
  }
  // Nothing standard matched. A symbol assembled from elements the positions could not
  // place is not a non-standard spelling, it is nonsense ("P2m"), so fall back to the
  // crystal class — which is at least true, if less specific.
  if (!coversAllElements(ops, pointGroup)) {
    return { symbol: `${centering}${pointGroup}`, standard: false, permutation: 0, placed: false };
  }
  return { symbol: direct[0] ?? `${centering}${pointGroup}`, standard: false, permutation: 0, placed: true };
}
