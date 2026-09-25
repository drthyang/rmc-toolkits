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

// What each symmetry-direction family of a crystal system may carry, as element labels
// reduced to their rotation part: '2' (2-fold axis, rotation or screw), '3', '4', '6'
// (rotations or screws of that order), '-3', '-4', '-6' (rotoinversions), 'm' (mirror or
// glide, by its normal). A conventional cell puts every element in a family that allows
// it — the 4-fold of a tetragonal group along [001], the 3-folds of a cubic group along
// <111>. Covering a family direction is not enough: on a primitive cubic-F cell the
// 4-folds lie along primitive <111>, a direction the cubic symbol reserves for 3-folds.
const FAMILY_ALLOWS = {
  monoclinic: [['2', 'm']],
  orthorhombic: [['2', 'm'], ['2', 'm'], ['2', 'm']],
  tetragonal: [['2', '4', '-4', 'm'], ['2', 'm'], ['2', 'm']],
  trigonal: [['3', '-3'], ['2', 'm'], ['2', 'm']],
  hexagonal: [['2', '3', '6', '-3', '-6', 'm'], ['2', 'm'], ['2', 'm']],
  cubic: [['2', '4', '-4', 'm'], ['3', '-3'], ['2', 'm']],
};

const elementType = (e) => (e.kind === 'mirror' || e.kind === 'glide' ? 'm'
  : e.kind === 'rotoinversion' ? `-${e.order}` : `${e.order}`);

/**
 * Whether the operations are in a conventional setting of their crystal system: every
 * symmetry element lies along a direction family of the system that may carry an element
 * of its type (FAMILY_ALLOWS), with a monoclinic unique axis along b. Stronger than
 * coversAllElements, which only asks that the direction belong to SOME family. Only a set
 * that passes can be named positionally.
 */
export function elementsFitSetting(ops, pointGroup) {
  const system = POINT_GROUP_SYSTEM[pointGroup];
  if (!system) return false;
  if (system === 'triclinic') return true;
  const elements = [];
  for (const { R, t } of ops) {
    const e = classifyElement(R, t);
    if (e) elements.push(e);
  }
  // The standard monoclinic setting has its unique axis along b. A cell with it along a
  // or c spells the same short symbol (P2_1/c) but is not the setting the Wyckoff tables
  // describe, so it is not accepted here; the axis orders tried next move it to b.
  const families = SYSTEM_DIRECTIONS[system];
  const familyKeys = families.map((family) => new Set(family.map(dirKey)));
  const allows = FAMILY_ALLOWS[system];
  return elements.every((e) => {
    const key = dirKey(e.direction);
    const type = elementType(e);
    return familyKeys.some((keys, i) => keys.has(key) && allows[i].includes(type));
  });
}

/* ── setting search ─────────────────────────────────────────────────────────
 * An RMC cell is under no obligation to be a conventional cell: the axes may be in
 * another order (Pbnm rather than Pnma), the principal axis along a or b, the cell a
 * primitive or other centred cell choice, or a supercell of the true cell. A symbol is
 * only meaningful in a standard setting, so the operations are re-expressed in a set of
 * candidate cells and named in the first one that is conventional. */

// The 6 axis orders as PROPER basis changes (det +1; an odd permutation flips one
// axis), so a screw's handedness (4_1 vs 4_3) survives the change. Columns = the new
// basis vectors in the old fractional basis.
const PERMUTATIONS = [
  [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
  [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
  [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
  [[0, 1, 0], [1, 0, 0], [0, 0, -1]],
  [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],
  [[0, 0, 1], [0, -1, 0], [1, 0, 0]],
];

/** Determinant and inverse of a real 3×3 matrix. */
function det3f(m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}
function inv3f(m) {
  const d = det3f(m);
  const c = (i, j) => {
    const r0 = i === 0 ? 1 : 0, r1 = i === 2 ? 1 : 2, s0 = j === 0 ? 1 : 0, s1 = j === 2 ? 1 : 2;
    return ((i + j) % 2 ? -1 : 1) * (m[r0][s0] * m[r1][s1] - m[r0][s1] * m[r1][s0]);
  };
  return [[c(0, 0) / d, c(1, 0) / d, c(2, 0) / d], [c(0, 1) / d, c(1, 1) / d, c(2, 1) / d], [c(0, 2) / d, c(1, 2) / d, c(2, 2) / d]];
}

/**
 * Integer matrix M (rows = new basis vectors in the given basis, det +1) such that M·A is
 * a reduced basis of the same lattice: repeatedly shorten each vector by integer
 * multiples of the others (pairwise size reduction, then b_k ± b_i ± b_j), shortest
 * first, until nothing shortens. In a reduced basis every lattice rotation has entries in
 * {-1, 0, 1}; in an oblique one (a cubic lattice on a 45° cell) some need ±2.
 */
export function reduceBasis(A) {
  const M = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  const vec = (m) => [0, 1, 2].map((k) => m[0] * A[0][k] + m[1] * A[1][k] + m[2] * A[2][k]);
  const len2 = (m) => { const v = vec(m); return v[0] * v[0] + v[1] * v[1] + v[2] * v[2]; };
  const shorter = (cand, cur) => len2(cand) < len2(cur) * (1 - 1e-9);
  for (let iter = 0; iter < 200; iter++) {
    let changed = false;
    M.sort((p, q) => len2(p) - len2(q));
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        if (i === j) continue;
        const u = vec(M[i]), v = vec(M[j]);
        const mu = Math.round((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / len2(M[i]));
        if (!mu) continue;
        const cand = M[j].map((x, k) => x - mu * M[i][k]);
        if (shorter(cand, M[j])) { M[j] = cand; changed = true; }
      }
    }
    for (let k = 0; k < 3; k++) {
      const [i, j] = [0, 1, 2].filter((x) => x !== k);
      for (const si of [-1, 1]) for (const sj of [-1, 1]) {
        const cand = M[k].map((x, c) => x + si * M[i][c] + sj * M[j][c]);
        if (shorter(cand, M[k])) { M[k] = cand; changed = true; }
      }
    }
    if (!changed) break;
  }
  if (det3i(M) < 0) M[2] = M[2].map((x) => -x);
  return M;
}

/**
 * The lattice vectors n + τ of ℤ³ + T (given fractional basis; τ ∈ T, the snapped pure
 * translations), shortest first, with n over a ±3 box of a REDUCED basis of ℤ³
 * (reduceBasis) and τ centred in that basis. Every short vector of the lattice is in it
 * whatever the shape of the given cell — a ±2 box of an oblique given cell misses some
 * (the conventional axes of a cubic crystal on a sheared cell) — and ±3 covers the
 * c axis of an R lattice on a reduced rhombohedral basis (3·r − 2a − b).
 * @returns {{v:number[], len:number}[]}  len = Cartesian length (Å) through A
 */
function shortLatticeVectors(translations, A, range = 3) {
  const M = reduceBasis(A);
  const Mi = inv3f(M).map((row) => row.map((x) => Math.round(x)));
  const out = [];
  for (const tau of translations) {
    // τ in reduced coordinates (x_r = M⁻ᵀ·x), centred on the origin.
    const tr = [0, 1, 2].map((i) => Mi[0][i] * tau[0] + Mi[1][i] * tau[1] + Mi[2][i] * tau[2]);
    const tc = tr.map((x) => x - Math.round(x));
    for (let a = -range; a <= range; a++) for (let b = -range; b <= range; b++) for (let c = -range; c <= range; c++) {
      const r = [a + tc[0], b + tc[1], c + tc[2]];
      // back to the given basis: x = Mᵀ·x_r
      const v = [0, 1, 2].map((i) => M[0][i] * r[0] + M[1][i] * r[1] + M[2][i] * r[2]);
      if (Math.abs(v[0]) + Math.abs(v[1]) + Math.abs(v[2]) < 1e-9) continue;
      const x = [0, 1, 2].map((k) => v[0] * A[0][k] + v[1] * A[1][k] + v[2] * A[2][k]);
      out.push({ v, len: Math.hypot(x[0], x[1], x[2]) });
    }
  }
  return out.sort((p, q) => p.len - q.len);
}

/** Re-express the operations in another basis Q (columns = new basis vectors): R' = Q⁻¹RQ, t' = Q⁻¹t. */
export function transformOps(ops, Q) {
  const Qi = inv3f(Q);
  return ops.map(({ R, t }) => ({
    R: matMul(matMul(Qi, R), Q).map((row) => row.map((v) => Math.round(v))),
    t: matVec(Qi, t).map(wrap1),
  }));
}

const nearInt = (x, tol = 1e-6) => Math.abs(x - Math.round(x)) <= tol;
const cyc = (x) => x - Math.round(x);                    // signed distance to the nearest integer
const nearMod1 = (u, v, tol) => Math.abs(cyc(u[0] - v[0])) <= tol && Math.abs(cyc(u[1] - v[1])) <= tol
  && Math.abs(cyc(u[2] - v[2])) <= tol;

/**
 * The group's pure translations (operations with R = I), snapped to exact fractions, with
 * the zero translation first; null when they do not snap. They form a finite group of
 * order n (zero included), so each has an order d dividing n — d·τ is a lattice vector —
 * and is snapped to round(d·τ)/d with the SMALLEST such d for which every component lies
 * within `tol` of the 1/d grid. d is capped at 1/(2·tol) (25), where neighbouring grid
 * points are 2·tol apart and a snap would be a guess. Reading the grid from the group order
 * covers every centering and any supercell fraction up to 1/25 (fifths, sevenths, ...)
 * with the same loose tolerance a noisy structure needs: a refined translation is a mean
 * over all sites, but at a loose τ with few sites it can sit 0.005 from its fraction.
 */
export function pureTranslations(ops, tol = 0.02) {
  const raw = [];
  for (const { R, t } of ops) {
    if (!sameMat(R, IDENTITY)) continue;
    const w = t.map(wrap1);
    if (!raw.some((u) => nearMod1(u, w, 1e-3))) raw.push(w);
  }
  if (!raw.some((u) => nearMod1(u, [0, 0, 0], tol))) raw.push([0, 0, 0]);
  const n = raw.length;
  const maxD = Math.floor(1 / (2 * tol) + 1e-9);
  const orders = [];
  for (let d = 1; d <= Math.min(n, maxD); d++) if (n % d === 0) orders.push(d);
  const out = [];
  for (const y of raw) {
    let snapped = null;
    for (const d of orders) {
      const k = y.map((x) => Math.round(x * d));
      if (y.every((x, i) => Math.abs(x - k[i] / d) <= tol)) { snapped = k.map((ki) => ((ki % d) + d) % d / d); break; }
    }
    if (!snapped) return null;
    if (!out.some((u) => nearMod1(u, snapped, 1e-9))) out.push(snapped);
  }
  if (out.length !== n) return null;              // two translations snapped together
  out.sort((u, v) => u[0] - v[0] || u[1] - v[1] || u[2] - v[2]);   // lexicographic: zero first
  return out;
}

// Bravais centering vectors of a conventional cell. R is recognised in its reverse
// setting only so that it can be refused: the Wyckoff tables use the obverse one.
const BRAVAIS = [
  ['P', []],
  ['A', [[0, 0.5, 0.5]]], ['B', [[0.5, 0, 0.5]]], ['C', [[0.5, 0.5, 0]]],
  ['I', [[0.5, 0.5, 0.5]]],
  ['F', [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]],
  ['R', [[2 / 3, 1 / 3, 1 / 3], [1 / 3, 2 / 3, 2 / 3]]],
  ['R(reverse)', [[1 / 3, 2 / 3, 1 / 3], [2 / 3, 1 / 3, 2 / 3]]],
];

/**
 * Bravais centering letter of a translation set (fractional, mod 1, with the zero
 * translation), or null when the set is not EXACTLY a Bravais centering: the finer
 * translation lattice of a supercell of the true cell (a perovskite in a 2×2×2 cell has
 * all eight (i/2, j/2, k/2)) contains the F vectors but is not F-centred.
 */
export function bravaisCentering(translations, tol = 1e-3) {
  const distinct = [];
  for (const t of translations) if (!distinct.some((u) => nearMod1(u, t, tol))) distinct.push(t);
  for (const [letter, vecs] of BRAVAIS) {
    if (distinct.length !== vecs.length + 1) continue;
    if (!distinct.some((u) => nearMod1(u, [0, 0, 0], tol))) continue;
    if (vecs.every((v) => distinct.some((u) => nearMod1(u, v, tol)))) return letter;
  }
  return null;
}

/**
 * Centering letter of the cell the operations are written in, read from ALL their pure
 * translations: 'P', 'A', 'B', 'C', 'I', 'F' or 'R' (either R setting), or null when the
 * translations are not exactly a Bravais centering (a supercell of the true cell).
 */
export function centeringOfOps(ops) {
  const translations = pureTranslations(ops);
  if (!translations) return null;
  const letter = bravaisCentering(translations);
  return letter === 'R(reverse)' ? 'R' : letter;
}

// Centerings a standard (ITA) setting of each crystal system uses.
const STANDARD_CENTERING = {
  monoclinic: ['P', 'C'], orthorhombic: ['P', 'A', 'C', 'I', 'F'], tetragonal: ['P', 'I'],
  trigonal: ['P', 'R'], hexagonal: ['P'], cubic: ['P', 'I', 'F'],
};

/**
 * The group re-expressed in the cell with basis Q (columns, old fractional coordinates).
 * Null unless Q is a right-handed basis of lattice vectors in which every rotation is an
 * integer matrix. `translations` are the group's pure translations (snapped, with 0).
 *
 * @returns {{ ops, translations, letter, Q, Qinv, ratio }|null}  ratio = new/old volume
 */
export function applySetting(ops, translations, Q) {
  const d = det3f(Q);
  if (!(d > 1e-9)) return null;
  for (let c = 0; c < 3; c++) {
    const v = [Q[0][c], Q[1][c], Q[2][c]];
    if (!translations.some((tau) => nearInt(v[0] - tau[0]) && nearInt(v[1] - tau[1]) && nearInt(v[2] - tau[2]))) return null;
  }
  const Qi = inv3f(Q);
  const reps = new Map();                   // one operation per rotation: the group is closed
  for (const op of ops) {
    const key = op.R.flat().join(',');
    if (reps.has(key)) continue;
    const Rn = matMul(matMul(Qi, op.R), Q);
    if (!Rn.every((row) => row.every((v) => nearInt(v)))) return null;
    reps.set(key, { R: Rn.map((row) => row.map((v) => Math.round(v))), t: matVec(Qi, op.t) });
  }
  const want = translations.length * d;
  if (!nearInt(want, 1e-6) || Math.round(want) < 1) return null;
  const count = Math.round(want);
  // The new cell's translations are ℤ³ + T modulo the new cell's lattice Qℤ³: the finite
  // group generated, mod 1, by Q⁻¹ of the given basis vectors and of T (closed by adding
  // generators until nothing new appears). No search box, so an oblique Q loses none.
  const clean = (v) => v.map((x) => { const y = wrap1(x); return nearInt(y, 1e-6) ? 0 : y; });
  const gens = [[1, 0, 0], [0, 1, 0], [0, 0, 1], ...translations].map((g) => clean(matVec(Qi, g)));
  const Tn = [[0, 0, 0]];
  for (let k = 0; k < Tn.length; k++) {
    for (const g of gens) {
      const v = clean([Tn[k][0] + g[0], Tn[k][1] + g[1], Tn[k][2] + g[2]]);
      if (!Tn.some((u) => nearMod1(u, v, 1e-6))) {
        Tn.push(v);
        if (Tn.length > count) return null;
      }
    }
  }
  if (Tn.length !== count) return null;
  const newOps = [];
  for (const { R, t } of reps.values()) {
    for (const tau of Tn) newOps.push({ R, t: [wrap1(t[0] + tau[0]), wrap1(t[1] + tau[1]), wrap1(t[2] + tau[2])] });
  }
  return { ops: newOps, translations: Tn, letter: bravaisCentering(Tn), Q, Qinv: Qi, ratio: d };
}

// Rotation type from determinant and trace (similarity invariants).
const isProperOrder = (R, n) => det3i(R) === 1 && trace(R) === { 2: -1, 3: 0, 4: 1, 6: 2 }[n];
const isMirror = (R) => det3i(R) === -1 && trace(R) === 1;
const isBar4 = (R) => det3i(R) === -1 && trace(R) === -1;

/**
 * Candidate conventional cells built from the symmetry elements themselves, for a group
 * whose own cell is not conventional: the principal axis along a or b, a centred or
 * primitive cell choice, or a supercell of the true cell (the lattice vectors used are
 * those of the full translation lattice, so a supercell reduces to the true cell).
 *
 *   cubic        a, b, c along the three 4-fold (-4 for -43m, 2-fold for 23/m-3) axes
 *   tetragonal   c along the 4 / -4 axis; a the shortest lattice vector ⊥ c; b = 4·a
 *   trig./hex.   c along the 3-fold; a = ± the shortest lattice vector ⊥ c; b = 3·a
 *   orthorh.     the three 2-fold axes / mirror normals, in all six orders
 *   monoclinic   b along the unique axis; (a, c) every unimodular pair from the two
 *                shortest lattice vectors ⊥ b (covers the cell choices and P2_1/n)
 *
 * Each basis vector is the SHORTEST lattice vector along its direction. Returned as
 * column matrices Q; the caller keeps those that are proper lattice bases.
 */
function derivedBases(ops, translations, pointGroup, A) {
  const system = POINT_GROUP_SYSTEM[pointGroup];
  const rotations = [];
  const seen = new Set();
  for (const { R } of ops) { const k = R.flat().join(','); if (!seen.has(k)) { seen.add(k); rotations.push(R); } }
  // Short lattice vectors of ℤ³ + T, shortest first (Cartesian length through A).
  const lattice = shortLatticeVectors(translations, A);
  const same = (u, v) => Math.abs(u[0] - v[0]) < 1e-9 && Math.abs(u[1] - v[1]) < 1e-9 && Math.abs(u[2] - v[2]) < 1e-9;
  const shortest = (pred) => lattice.find(({ v }) => pred(v))?.v ?? null;
  const fixedBy = (R, sign) => (v) => same(matVec(R, v), v.map((x) => sign * x));
  const parallel = (u, v) => Math.abs(u[1] * v[2] - u[2] * v[1]) < 1e-9 && Math.abs(u[2] * v[0] - u[0] * v[2]) < 1e-9
    && Math.abs(u[0] * v[1] - u[1] * v[0]) < 1e-9;
  const columns = (a, b, c) => [[a[0], b[0], c[0]], [a[1], b[1], c[1]], [a[2], b[2], c[2]]];
  const neg = (v) => v.map((x) => -x);
  const rightHanded = (a, b, c) => (det3f(columns(a, b, c)) > 0 ? columns(a, b, c) : columns(a, b, neg(c)));
  // Distinct axis lines of the operations matching `pick`, each as its shortest lattice
  // vector (R·v = sign·v: +1 along a rotation axis, −1 along a rotoinversion axis or a
  // mirror normal).
  const axes = (pick, sign) => {
    const out = [];
    for (const R of rotations) {
      if (!pick(R)) continue;
      const v = shortest(fixedBy(R, sign));
      if (v && !out.some((u) => parallel(u, v))) out.push(v);
    }
    return out;
  };
  const out = [];
  if (system === 'cubic') {
    const lines = pointGroup === '-43m' ? axes(isBar4, -1)
      : (pointGroup === '432' || pointGroup === 'm-3m') ? axes((R) => isProperOrder(R, 4), 1)
        : axes((R) => isProperOrder(R, 2), 1);
    if (lines.length === 3) out.push(rightHanded(...lines));
  } else if (system === 'tetragonal') {
    const R4 = rotations.find((R) => isProperOrder(R, 4)) || rotations.find(isBar4);
    if (R4) {
      const c = shortest(fixedBy(R4, det3i(R4)));
      const R2 = matMul(R4, R4);
      const a = shortest(fixedBy(R2, -1));
      if (a && c) {
        for (const a0 of [a, matVec(R4, a).map((x, i) => x + a[i])]) {
          const b0 = matVec(R4, a0);
          out.push(det3f(columns(a0, b0, c)) > 0 ? columns(a0, b0, c) : columns(a0, neg(b0), c));
        }
      }
    }
  } else if (system === 'trigonal' || system === 'hexagonal') {
    const R3 = rotations.find((R) => isProperOrder(R, 3));
    if (R3) {
      const c = shortest(fixedBy(R3, 1));
      const inPlane = (v) => { const w = matVec(R3, v); const u = matVec(R3, w); return same([v[0] + w[0] + u[0], v[1] + w[1] + u[1], v[2] + w[2] + u[2]], [0, 0, 0]); };
      const a = shortest(inPlane);
      if (a && c) {
        for (const a0 of [a, neg(a)]) {
          const b1 = matVec(R3, a0);
          const b2 = matVec(R3, b1);
          out.push(det3f(columns(a0, b1, c)) > 0 ? columns(a0, b1, c) : columns(a0, b2, c));
        }
      }
    }
  } else if (system === 'orthorhombic') {
    const lines = [...axes((R) => isProperOrder(R, 2), 1), ...axes(isMirror, -1)]
      .filter((v, i, all) => all.findIndex((u) => parallel(u, v)) === i);
    if (lines.length === 3) {
      for (const [i, j, k] of [[0, 1, 2], [1, 2, 0], [2, 0, 1], [1, 0, 2], [0, 2, 1], [2, 1, 0]]) out.push(rightHanded(lines[i], lines[j], lines[k]));
    }
  } else if (system === 'monoclinic') {
    const U = rotations.find((R) => isProperOrder(R, 2)) || rotations.find(isMirror);
    if (U) {
      const proper = det3i(U) === 1;
      const b = shortest(fixedBy(U, proper ? 1 : -1));
      const plane = fixedBy(U, proper ? -1 : 1);
      const v1 = shortest(plane);
      const v2 = v1 ? shortest((v) => plane(v) && !parallel(v, v1)) : null;
      if (b && v1 && v2) {
        const combos = [];
        for (let p = -1; p <= 1; p++) for (let q = -1; q <= 1; q++) if (p || q) combos.push([p, q]);
        for (const [p1, q1] of combos) {
          for (const [p2, q2] of combos) {
            if (Math.abs(p1 * q2 - q1 * p2) !== 1) continue;
            const a = [0, 1, 2].map((i) => p1 * v1[i] + q1 * v2[i]);
            const c = [0, 1, 2].map((i) => p2 * v1[i] + q2 * v2[i]);
            out.push(det3f(columns(a, b, c)) > 0 ? columns(a, b, c) : columns(a, neg(b), c));
          }
        }
      }
    }
  }
  return out;
}

// Pairs whose members have the same element types along the same directions: the I
// centering turns each 2-fold into a 2_1 half a cell away and back. [symmorphic, other]
const LOCATION_PAIRS = [['I222', 'I2_12_12_1'], ['I23', 'I2_13']];

/**
 * Whether the 2-fold rotations along a, b and c (standard setting, diagonal R) have a
 * common fixed point — true for I222 and I23, whose 2-folds meet at the origin, false for
 * I2_12_12_1 and I2_13, whose 2-folds never meet. A pure (non-screw) 2-fold along axis i,
 * {R|t} with t_i ≡ 0, fixes the line 2p_j ≡ t_j (j ≠ i); lines along a and b meet the one
 * along c at a common point iff their translations agree mod 1 on the shared components.
 * Independent of the origin.
 */
export function twoFoldsMeet(ops, tol = 0.05) {
  const pure = [0, 1, 2].map((i) => ops.filter(({ R, t }) => R.every((row, r) => row.every((v, c) => v === (r === c ? (r === i ? 1 : -1) : 0)))
    && Math.abs(cyc(t[i])) <= tol));
  const agree = (u, v) => Math.abs(cyc(u - v)) <= tol;
  for (const ox of pure[0]) {
    for (const oy of pure[1]) {
      if (!agree(ox.t[2], oy.t[2])) continue;
      for (const oz of pure[2]) if (agree(ox.t[1], oz.t[1]) && agree(oy.t[0], oz.t[0])) return true;
    }
  }
  return false;
}

// The member of a location-degenerate pair the operations actually form.
function resolveLocationPair(symbol, ops) {
  const pair = LOCATION_PAIRS.find((p) => p.includes(symbol));
  if (!pair) return symbol;
  return twoFoldsMeet(ops) ? pair[0] : pair[1];
}

/**
 * H–M symbol for a closed group, searching cells until one is conventional.
 *
 * The group is re-expressed (applySetting) in: its own cell, the five other axis orders,
 * and cells built from its symmetry elements (derivedBases). In each, the centering letter
 * is read from the FULL translation lattice of that cell (bravaisCentering) and must be
 * one a standard setting of the crystal system uses. A candidate symbol is accepted only
 * if it is a tabulated symbol OF THE DETECTED CLASS (`classOf` returns the crystal class
 * of a tabulated symbol, or null — injected so this module stays table-free), starts with
 * that letter, and every element lies where its type belongs (elementsFitSetting).
 * Positional assembly on a cell that is not conventional can spell another group's symbol
 * ("Pmmm" for rocksalt on its primitive cell), which is why the class, centering and
 * placement are checked, not only the spelling.
 *
 * When nothing is accepted the symbol is null: the caller reports the crystal class,
 * never a symbol that was not verified. `setting` is the cell the symbol belongs to
 * ({ Q, Qinv, translations, ratio, ops }), for placing Wyckoff positions.
 *
 * @param {{R,t}[]} ops  every operation of the group (pure translations included)
 * @param {string|null} centering  unused (kept for call compatibility; the letter is
 *   read from the translations in each candidate cell)
 * @param {string} pointGroup
 * @param {(symbol:string)=>string|null} classOf
 * @param {{A?:number[][], holohedry?:function, latticeTol?:number}} [options]  A = the
 *   cell's lattice rows (Å), used to pick the shortest lattice vectors of derived cells
 *   (without it only the axis orders are tried); holohedry(A', tol) = the lattice rotations
 *   of a cell within a strain of tol Å, used to tell whether the translation lattice has
 *   symmetries the given cell could not test (`complete`); latticeTol = that strain
 *   tolerance (the group's worst residual)
 * @returns {{ symbol:string|null, standard:boolean, setting:object|null, placed:boolean,
 *   complete?:boolean }}
 */
export function hmSymbolInStandardSetting(ops, centering, pointGroup, classOf, { A = null, holohedry = null, latticeTol = 1e-3 } = {}) {
  const system = POINT_GROUP_SYSTEM[pointGroup];
  const none = { symbol: null, standard: false, setting: null, placed: coversAllElements(ops, pointGroup) };
  if (!system || system === 'triclinic' || typeof classOf !== 'function') return none;
  const translations = pureTranslations(ops);
  if (!translations) return none;

  const allowed = STANDARD_CENTERING[system];
  const tryBasis = (Q) => {
    const setting = applySetting(ops, translations, Q);
    if (!setting || !allowed.includes(setting.letter)) return null;
    if (!elementsFitSetting(setting.ops, pointGroup)) return null;
    for (const cand of hmSymbolCandidates(setting.ops, setting.letter, pointGroup)) {
      if (classOf(cand) === pointGroup && cand.startsWith(setting.letter)) {
        return { symbol: resolveLocationPair(cand, setting.ops), setting };
      }
    }
    return null;
  };
  const bases = [...PERMUTATIONS];
  if (A) bases.push(...derivedBases(ops, translations, pointGroup, A));
  for (const Q of bases) {
    const found = tryBasis(Q);
    if (found) {
      const complete = !A || typeof holohedry !== 'function' || allLatticeOpsTested(translations, A, holohedry, latticeTol);
      return { symbol: found.symbol, standard: true, setting: found.setting, placed: true, complete };
    }
  }
  return none;
}

/**
 * A primitive basis (columns, in the given cell's fractional basis) of the lattice
 * ℤ³ + T, from its successive minima — in three dimensions the shortest vector, the
 * shortest one not parallel to it and the shortest one not coplanar with both always form
 * a basis. Null if they do not have the primitive volume 1/|T| (never expected).
 */
function primitiveBasis(translations, A) {
  const vectors = shortLatticeVectors(translations, A);
  const cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
  const v1 = vectors[0]?.v;
  if (!v1) return null;
  const v2 = vectors.find(({ v }) => cross(v1, v).some((x) => Math.abs(x) > 1e-9))?.v;
  if (!v2) return null;
  const n = cross(v1, v2);
  const found = vectors.find(({ v }) => Math.abs(n[0] * v[0] + n[1] * v[1] + n[2] * v[2]) > 1e-9);
  if (!found) return null;
  let v3 = found.v;
  const vol = n[0] * v3[0] + n[1] * v3[1] + n[2] * v3[2];
  if (Math.abs(Math.abs(vol) - 1 / translations.length) > 1e-6) return null;
  if (vol < 0) v3 = v3.map((x) => -x);
  return [[v1[0], v2[0], v3[0]], [v1[1], v2[1], v3[1]], [v1[2], v2[2], v3[2]]];
}

/**
 * Whether every symmetry of the crystal's translation lattice could have been tested in
 * the cell the operations were found in. The finder only tries rotations that map the
 * GIVEN cell's lattice ℤ³ onto itself (integer matrices there). When that cell is a
 * supercell of the translation lattice ℤ³ + T — a perovskite in a √2×√2×2 or 2×2×1 cell,
 * rocksalt in a 1×1×2 cell — rotations of the full lattice that do not preserve ℤ³ (the
 * cubic 3-folds there) were never tried, and the group found is only a lower bound.
 * `holohedry(A', tol)` lists the lattice rotations of a cell with lattice rows A' (strain
 * ≤ tol Å); they are taken on a primitive basis of ℤ³ + T and carried to the given cell.
 */
function allLatticeOpsTested(translations, A, holohedry, tol) {
  const P = primitiveBasis(translations, A);
  if (!P) return false;
  const Pi = inv3f(P);
  const Ap = [0, 1, 2].map((i) => [0, 1, 2].map((k) => P[0][i] * A[0][k] + P[1][i] * A[1][k] + P[2][i] * A[2][k]));
  for (const R of holohedry(Ap, tol)) {
    const back = matMul(matMul(P, R), Pi);
    if (!back.every((row) => row.every((v) => nearInt(v)))) return false;
  }
  return true;
}

/**
 * Whether the lattice symmetries of the operations' full translation lattice could all be
 * tested in the given cell (see allLatticeOpsTested) — for callers naming a group without
 * the setting search (P1 / P-1). False when a pure translation cannot be snapped.
 */
export function latticeFullyTested(ops, A, holohedry, tol = 1e-3) {
  const translations = pureTranslations(ops);
  return !!translations && allLatticeOpsTested(translations, A, holohedry, tol);
}
