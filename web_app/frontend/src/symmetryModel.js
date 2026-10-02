// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/symmetryModel.js
//
// Glue between the parsed RMC structure (browserData.structureFromRmc6f) and the
// pure symmetry finder (symmetry.js): builds the conventional cell + basis and
// returns a space-group description + tolerance ladder for the UI.

import { spaceGroupAtTolerance, symmetryLadder, siteOrbits, operationEstimate } from './symmetry.js';
import { assignWyckoffLetters } from './wyckoff.js';

/**
 * Largest average-structure basis the finder analyses. It runs synchronously on the main
 * thread (ModelSummary's useMemo), and a box with one reference site per atom — a glass,
 * or an imported configuration declared as a 1×1×1 supercell — is not a unit-cell
 * configuration anyway; above this the card says so instead of freezing the page.
 */
export const MAX_SYMMETRY_SITES = 2000;

const tooLarge = (structure) => structure.basis.length > MAX_SYMMETRY_SITES;

/**
 * Most candidate operations (lattice rotations × pure translations, operationEstimate) the
 * page analyses. A correctly declared cell has at most 48 × 4 = 192; this allows twice
 * that (a 2×2×2 supercell of a primitive cubic cell). A crystalline box declared as a
 * 1×1×1 supercell has one pure translation per repeat unit, and the candidates then grow
 * with the square of the box: a 4×4×4 rocksalt box (12 288) took seconds, a noisy 3×3×3
 * one (5184) minutes, on the main thread.
 */
export const MAX_SYMMETRY_OPS = 384;

/**
 * The loosest tolerance the card's ladder explores (ModelSummary calls toleranceLadder
 * with it). The operation budget is judged there for the headline too, so the card
 * either analyses a structure at every tolerance or explains why it does not.
 */
const LADDER_TOL_MAX = 1.0;

/** A 'not analysed' description, rendered as-is by the card (headline + subtitle). */
const notAnalysed = (subtitle, reason) => ({
  skipped: true,
  reason,
  spaceGroup: 'not analysed',
  spaceGroupNumber: null,
  pointGroup: subtitle,
  centering: null,
  nSpace: '—',
  nPoint: 0,
  maxResidual: Number.NaN,   // nothing was fitted — never 0, which would read as an exact fit
  orbits: [],
});

// The operation budget at the ladder's loosest tolerance (or `tol`, if looser).
const overBudget = (A, basis, tol) => operationEstimate(A, basis, Math.max(tol, LADDER_TOL_MAX, 1e-3), MAX_SYMMETRY_OPS);

/** Conventional unit cell A_conv (rows, Å) = supercell lattice / supercell dims. */
export function conventionalCell(structure) {
  const { latticeVectors, supercell } = structure;
  return latticeVectors.map((row, i) => row.map((v) => v / Math.max(supercell[i], 1)));
}

/** Mean cell edge (Å), for converting a cartesian tolerance into cell fractions. */
const meanEdge = (A) => A.reduce((sum, row) => sum + Math.hypot(row[0], row[1], row[2]), 0) / 3;

const wrap01 = (x) => { const y = x - Math.floor(x); return y > 1 - 1e-9 ? 0 : y; };

/**
 * Wyckoff positions of the orbits, read in the standard cell the group was named in
 * (`setting`: Q = new basis vectors as columns in the given cell's fractional basis,
 * Qinv, the new cell's pure translations, and ratio = new/old cell volume). The tables
 * describe the standard cell, so each orbit's members are carried into it (x' = Q⁻¹x,
 * plus every translation of that cell, so any tabulated representative can be matched)
 * and its multiplicity scaled to that cell. Returns one { letter, multiplicity } per
 * orbit — both of the standard cell, since a label pairs them — or nulls where no letter
 * is assigned. With no setting — an unnamed group, a lower bound, a triclinic group in a
 * centred cell — every letter is withheld.
 */
function lettersInSetting(sg, found, basis, A, tol) {
  const { setting } = sg;
  const none = { letter: null, multiplicity: null };
  if (!sg.spaceGroupNumber || !setting) return found.map(() => none);
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
  const letters = assignWyckoffLetters(sg.spaceGroupNumber, orbits, moved, tol / meanEdge(As));
  return letters.map((letter, i) => (letter ? { letter, multiplicity: orbits[i].size } : none));
}

/**
 * The finder's full result at a cartesian tolerance `tol` (Å), for callers that need more
 * than the card's summary (the symmetry-averaged CIF export, averageStructure.js):
 *   { A, sg, found, positions }
 * A = the conventional cell (rows, Å); sg = spaceGroupAtTolerance's closed group with its
 * operations {R, t, residual} in A's fractional basis and the standard `setting` it was
 * named in (null when there is none); found = siteOrbits (index lists into
 * structure.basis); positions = one { letter, multiplicity } per orbit in that setting.
 * Returns null without a basis, and the card's 'not analysed' description (skipped: true)
 * when the structure is too large or over the operation budget.
 */
export function analyseSymmetry(structure, tol = 0.2) {
  if (!structure?.basis?.length || !structure?.latticeVectors) return null;
  if (tooLarge(structure)) {
    const n = structure.basis.length;
    return notAnalysed(`${n} sites > ${MAX_SYMMETRY_SITES} limit`,
      `The average structure has ${n} reference sites; symmetry detection runs in the browser `
      + `and is limited to ${MAX_SYMMETRY_SITES} sites (a box with one reference site per atom is not `
      + 'a unit-cell configuration).');
  }
  const A = conventionalCell(structure);
  const budget = overBudget(A, structure.basis, tol);
  if (budget.exceeds) {
    const t = budget.translations;
    return notAnalysed(`≥ ${t} translations per cell`,
      `At least ${t} pure translations map the average structure onto itself within `
      + `${Math.max(tol, LADDER_TOL_MAX)} Å, so the cell is a supercell of the structure's repeat unit `
      + '(check the .rmc6f "Supercell dimensions"). With its lattice rotations that is more than '
      + `${MAX_SYMMETRY_OPS} candidate operations, too many to analyse on the page.`);
  }
  // Detected in the ladder's own pass (pairing and near-duplicate radii of LADDER_TOL_MAX),
  // then walked to tol: the headline is then exactly the ladder's group at tol.
  const sg = spaceGroupAtTolerance(A, structure.basis, tol, Math.max(tol, LADDER_TOL_MAX));   // a closed group, or 'undetermined'
  // No operation at all (a broken lattice): no orbits either — not one orbit per site.
  const found = sg.ops.length ? siteOrbits(A, structure.basis, sg.ops, tol) : [];
  const positions = lettersInSetting(sg, found, structure.basis, A, tol);
  return { A, sg, found, positions };
}

/**
 * Space group of a parsed structure at a cartesian tolerance `tol` (Å):
 *   { spaceGroup, spaceGroupNumber, pointGroup, centering, nSpace, nPoint,
 *     maxResidual, orbits:[{ element, size, site, rep, wyckoff, wyckoffMultiplicity, members }] }
 * Returns null when the structure has no basis (not yet loaded / no reference sites).
 * `size` is the orbit's multiplicity in the GIVEN cell; `wyckoff` is the letter in the
 * standard cell the group is named in and `wyckoffMultiplicity` the multiplicity there
 * (they can differ from `size`: label a position with the pair, as orbitLabel does).
 */
export function describeSymmetry(structure, tol = 0.2) {
  const analysis = analyseSymmetry(structure, tol);
  if (!analysis || analysis.skipped) return analysis;
  const { sg, found, positions } = analysis;
  const orbits = found.map((o, i) => ({
    element: o.element,
    size: o.size,
    site: o.site,
    rep: o.rep,
    wyckoff: positions[i].letter,
    // The multiplicity that goes with the letter: the orbit's size in the standard cell the
    // letter is read in, which differs from `size` when that is not the given cell (4 Ga
    // in an F-cubic cell are R3m's 3a). null when there is no letter.
    wyckoffMultiplicity: positions[i].multiplicity,
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
  const A = conventionalCell(structure);
  if (overBudget(A, structure.basis, tolMax).exceeds) return [];
  return symmetryLadder(A, structure.basis, tolMax);
}

/**
 * Wyckoff label for an orbit: the standard cell's multiplicity + letter ('3a'), or, with no
 * letter, the given cell's multiplicity + site symmetry ('4 (3m)').
 */
export function orbitLabel(orbit) {
  if (orbit.wyckoff) return `${orbit.wyckoffMultiplicity ?? orbit.size}${orbit.wyckoff}`;
  return `${orbit.size} (${orbit.site})`;
}
