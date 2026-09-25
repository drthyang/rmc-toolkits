// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/symmetryModel.js
//
// Glue between the parsed RMC structure (browserData.structureFromRmc6f) and the
// pure symmetry finder (symmetry.js): builds the conventional cell + basis and
// returns a space-group description + tolerance ladder for the UI.

import { spaceGroupAtTolerance, symmetryLadder, siteOrbits } from './symmetry.js';
import { assignWyckoffLetters } from './wyckoff.js';

/**
 * Largest average-structure basis the finder analyses. It runs synchronously on the main
 * thread (ModelSummary's useMemo), and a box with one reference site per atom — a glass,
 * or an imported configuration declared as a 1×1×1 supercell — is not a unit-cell
 * configuration anyway; above this the card says so instead of freezing the page.
 */
export const MAX_SYMMETRY_SITES = 2000;

const tooLarge = (structure) => structure.basis.length > MAX_SYMMETRY_SITES;

/** Conventional unit cell A_conv (rows, Å) = supercell lattice / supercell dims. */
export function conventionalCell(structure) {
  const { latticeVectors, supercell } = structure;
  return latticeVectors.map((row, i) => row.map((v) => v / Math.max(supercell[i], 1)));
}

/** Mean cell edge (Å), for converting a cartesian tolerance into cell fractions. */
const meanEdge = (A) => A.reduce((sum, row) => sum + Math.hypot(row[0], row[1], row[2]), 0) / 3;

const wrap01 = (x) => { const y = x - Math.floor(x); return y > 1 - 1e-9 ? 0 : y; };

/**
 * Wyckoff letters for the orbits, read in the standard cell the group was named in
 * (`setting`: Q = new basis vectors as columns in the given cell's fractional basis,
 * Qinv, the new cell's pure translations, and ratio = new/old cell volume). The tables
 * describe the standard cell, so each orbit's members are carried into it (x' = Q⁻¹x,
 * plus every translation of that cell, so any tabulated representative can be matched)
 * and its multiplicity scaled to that cell. With no setting — an unnamed group, a lower
 * bound, a triclinic group in a centred cell — every letter is withheld.
 */
function lettersInSetting(sg, found, basis, A, tol) {
  const { setting } = sg;
  if (!sg.spaceGroupNumber || !setting) return found.map(() => null);
  const { Q, Qinv, translations, ratio } = setting;
  const As = [0, 1, 2].map((i) => [0, 1, 2].map((k) => Q[0][i] * A[0][k] + Q[1][i] * A[1][k] + Q[2][i] * A[2][k]));
  const moved = [];
  const orbits = found.map((o) => {
    const size = o.size * ratio;
    const index = [];
    for (const i of o.index) {
      const x = basis[i].frac;
      const p = [0, 1, 2].map((r) => Qinv[r][0] * x[0] + Qinv[r][1] * x[1] + Qinv[r][2] * x[2]);
      for (const tau of translations) {
        index.push(moved.length);
        moved.push({ frac: p.map((v, r) => wrap01(v + tau[r])) });
      }
    }
    // A fractional multiplicity means the orbit does not fit the cell: no letter.
    return { size: Math.abs(size - Math.round(size)) < 1e-6 ? Math.round(size) : -1, site: o.site, index };
  });
  return assignWyckoffLetters(sg.spaceGroupNumber, orbits, moved, tol / meanEdge(As));
}

/**
 * Space group of a parsed structure at a cartesian tolerance `tol` (Å):
 *   { spaceGroup, spaceGroupNumber, pointGroup, centering, nSpace, nPoint,
 *     maxResidual, orbits:[{ element, size, site, rep, wyckoff }] }
 * Returns null when the structure has no basis (not yet loaded / no reference sites).
 * `size` is the orbit's multiplicity in the GIVEN cell; `wyckoff` is the letter in the
 * standard cell the group is named in (the multiplicity there can differ).
 */
export function describeSymmetry(structure, tol = 0.2) {
  if (!structure?.basis?.length || !structure?.latticeVectors) return null;
  if (tooLarge(structure)) {
    const n = structure.basis.length;
    // Rendered as-is by the card: spaceGroup is the headline, pointGroup its subtitle.
    return {
      skipped: true,
      reason: `The average structure has ${n} reference sites; symmetry detection runs in the browser `
        + `and is limited to ${MAX_SYMMETRY_SITES} sites (a box with one reference site per atom is not `
        + 'a unit-cell configuration).',
      spaceGroup: 'not analysed',
      spaceGroupNumber: null,
      pointGroup: `${n} sites > ${MAX_SYMMETRY_SITES} limit`,
      centering: null,
      nSpace: '—',
      nPoint: 0,
      maxResidual: Number.NaN,   // nothing was fitted — never 0, which would read as an exact fit
      orbits: [],
    };
  }
  const A = conventionalCell(structure);
  const sg = spaceGroupAtTolerance(A, structure.basis, tol);   // a closed group, or 'undetermined'
  // No operation at all (a broken lattice): no orbits either — not one orbit per site.
  const found = sg.ops.length ? siteOrbits(A, structure.basis, sg.ops, tol) : [];
  const letters = lettersInSetting(sg, found, structure.basis, A, tol);
  const orbits = found.map((o, i) => ({
    element: o.element,
    size: o.size,
    site: o.site,
    rep: o.rep,
    wyckoff: letters[i],
    // Indices into structure.basis, so callers can aggregate per-site data
    // (e.g. rms displacements) over each orbit's member sites.
    members: o.index,
  }));
  return {
    spaceGroup: sg.spaceGroup,
    spaceGroupNumber: sg.spaceGroupNumber,
    pointGroup: sg.pointGroup,
    centering: sg.centering,
    nSpace: sg.nSpace,
    nPoint: sg.nPoint,
    maxResidual: sg.maxResidual,
    orbits,
  };
}

/** Symmetry-vs-tolerance ladder (bricks tight→loose) for the structure. */
export function toleranceLadder(structure, tolMax = 1.0) {
  if (!structure?.basis?.length || !structure?.latticeVectors || tooLarge(structure)) return [];
  return symmetryLadder(conventionalCell(structure), structure.basis, tolMax);
}

/** Wyckoff label for an orbit: multiplicity + letter, or multiplicity + site symmetry. */
export function orbitLabel(orbit) {
  return `${orbit.size}${orbit.wyckoff || ` (${orbit.site})`}`;
}
