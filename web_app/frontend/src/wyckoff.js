// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/wyckoff.js
//
// Wyckoff letters for the symmetry orbits found by symmetry.js.
//
// An orbit already knows its multiplicity and its site-symmetry point group, both derived
// from the operations rather than looked up. Those two together usually name the Wyckoff
// position outright. When they don't — Pm-3m has 3c and 3d, both multiplicity 3 with site
// symmetry 4/mmm — the tie is broken by the position's COORDINATE FORM: 3c is (0,½,½) and
// 3d is (½,0,0), so testing which form the orbit's atoms actually fit separates them.
//
// Two wrinkles the form test has to survive. The representative of a detected orbit is
// whichever atom happened to come first, not the one International Tables prints, so every
// member of the orbit is tried and any match counts. And a form like (x,x,z) constrains
// its coordinates only modulo a lattice vector, so the fit is attempted over integer shifts
// of the point as well.
//
// A letter is returned only when exactly one position fits. The app never shifts the origin
// to match the standard setting, so a structure described on a different origin genuinely
// cannot be assigned letters — reporting none is the honest answer, and the caller falls
// back to showing multiplicity and site symmetry.

import { WYCKOFF_DATA } from './wyckoffTable.js';

const parsedCache = new Map();

const fraction = (s) => (s.includes('/') ? Number(s.split('/')[0]) / Number(s.split('/')[1]) : Number(s));

/** Parse a coordinate form ("x,1/4,z") into { R, t }, so the point is R·(x,y,z) + t. */
export function parseCoordinateForm(text) {
  const R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  const t = [0, 0, 0];
  text.replace(/\s+/g, '').split(',').forEach((part, i) => {
    for (const token of part.match(/[+-]?[^+-]+/g) || []) {
      const sign = token.startsWith('-') ? -1 : 1;
      const body = token.replace(/^[+-]/, '');
      const axis = 'xyz'.indexOf(body[body.length - 1]);
      if (axis >= 0) {
        const coefficient = body.slice(0, -1);
        R[i][axis] = sign * (coefficient === '' ? 1 : fraction(coefficient));
      } else {
        t[i] += sign * fraction(body);
      }
    }
  });
  return { R, t };
}

/**
 * Wyckoff positions of a space group: [{ letter, multiplicity, site, coordinate, form }].
 * Empty when the group is not tabulated. Parsed on first use and cached.
 */
export function wyckoffPositions(spaceGroupNumber) {
  if (parsedCache.has(spaceGroupNumber)) return parsedCache.get(spaceGroupNumber);
  const packed = WYCKOFF_DATA[spaceGroupNumber];
  const rows = !packed ? [] : packed.split(';').map((row) => {
    // letter : multiplicity : site symmetry : coordinate  (the coordinate holds the commas)
    const [letter, multiplicity, site, coordinate] = row.split(':');
    return { letter, multiplicity: Number(multiplicity), site, coordinate, form: parseCoordinateForm(coordinate) };
  });
  parsedCache.set(spaceGroupNumber, rows);
  return rows;
}

/* ── fitting a point to a coordinate form ───────────────────────────────────── */

// Solve (RᵀR + λI)·v = Rᵀd for the free parameters. R is often rank-deficient — a fully
// fixed position such as (0,0,0) has R = 0 — so a small ridge keeps it well posed; the
// residual check below is what actually decides the fit.
function solveRidge(R, d) {
  const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  const b = [0, 0, 0];
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      let s = 0;
      for (let k = 0; k < 3; k++) s += R[k][i] * R[k][j];
      A[i][j] = s + (i === j ? 1e-9 : 0);
    }
    let s = 0;
    for (let k = 0; k < 3; k++) s += R[k][i] * d[k];
    b[i] = s;
  }
  const det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
    - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
    + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
  if (Math.abs(det) < 1e-18) return [0, 0, 0];
  const col = (c, v) => A.map((row, i) => row.map((x, j) => (j === c ? v[i] : x)));
  const d3 = (m) => m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  return [d3(col(0, b)) / det, d3(col(1, b)) / det, d3(col(2, b)) / det];
}

const SHIFTS = (() => {
  const out = [];
  for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) for (let c = -1; c <= 1; c++) out.push([a, b, c]);
  return out;
})();

/** Whether `point` (fractional) can be written in the coordinate form, within `tol`. */
export function fitsForm(form, point, tol = 0.05) {
  for (const shift of SHIFTS) {
    const p = [point[0] + shift[0], point[1] + shift[1], point[2] + shift[2]];
    const d = [p[0] - form.t[0], p[1] - form.t[1], p[2] - form.t[2]];
    const v = solveRidge(form.R, d);
    let ok = true;
    for (let i = 0; i < 3 && ok; i++) {
      const fitted = form.R[i][0] * v[0] + form.R[i][1] * v[1] + form.R[i][2] * v[2] + form.t[i];
      if (Math.abs(fitted - p[i]) > tol) ok = false;
    }
    if (ok) return true;
  }
  return false;
}

/**
 * Wyckoff letter for each orbit, or null where it cannot be pinned down.
 *
 * @param {number|null} spaceGroupNumber
 * @param {{size:number, site:string, index:number[]}[]} orbits  from siteOrbits()
 * @param {{frac:number[]}[]} basis  the sites the orbit indices point into
 * @param {number} tolFrac  position tolerance in cell fractions
 * @returns {(string|null)[]} one letter per orbit, in the same order
 */
export function assignWyckoffLetters(spaceGroupNumber, orbits, basis, tolFrac = 0.05) {
  const rows = spaceGroupNumber ? wyckoffPositions(spaceGroupNumber) : [];
  if (!rows.length) return orbits.map(() => null);
  return orbits.map((orbit) => {
    const candidates = rows.filter((r) => r.multiplicity === orbit.size && r.site === orbit.site);
    if (candidates.length === 1) return candidates[0].letter;
    if (candidates.length === 0) return null;
    const members = orbit.index.map((i) => basis[i]?.frac).filter(Boolean);
    const fits = candidates.filter((c) => members.some((p) => fitsForm(c.form, p, tolFrac)));
    return fits.length === 1 ? fits[0].letter : null;
  });
}

/** Number of space groups with a tabulated Wyckoff table. */
export const tabulatedGroupCount = () => Object.keys(WYCKOFF_DATA).length;
