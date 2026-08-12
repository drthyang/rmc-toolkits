// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/symmetryModel.js
//
// Glue between the parsed RMC structure (browserData.structureFromRmc6f) and the
// pure symmetry finder (symmetry.js): builds the conventional cell + basis and
// returns a space-group description + tolerance ladder for the UI.

import { spaceGroupAtTolerance, symmetryLadder, siteOrbits } from './symmetry.js';
import { assignWyckoffLetters } from './wyckoff.js';

/** Conventional unit cell A_conv (rows, Å) = supercell lattice / supercell dims. */
export function conventionalCell(structure) {
  const { latticeVectors, supercell } = structure;
  return latticeVectors.map((row, i) => row.map((v) => v / Math.max(supercell[i], 1)));
}

/** Mean cell edge (Å), for converting a cartesian tolerance into cell fractions. */
const meanEdge = (A) => A.reduce((sum, row) => sum + Math.hypot(row[0], row[1], row[2]), 0) / 3;

/**
 * Space group of a parsed structure at a cartesian tolerance `tol` (Å):
 *   { spaceGroup, spaceGroupNumber, pointGroup, centering, nSpace, nPoint,
 *     maxResidual, orbits:[{ element, size, site, rep, wyckoff }] }
 * Returns null when the structure has no basis (not yet loaded / no reference sites).
 */
export function describeSymmetry(structure, tol = 0.2) {
  if (!structure?.basis?.length || !structure?.latticeVectors) return null;
  const A = conventionalCell(structure);
  const sg = spaceGroupAtTolerance(A, structure.basis, tol);   // always a valid group
  const found = siteOrbits(A, structure.basis, sg.ops, tol);
  const letters = assignWyckoffLetters(sg.spaceGroupNumber, found, structure.basis, tol / meanEdge(A));
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
  if (!structure?.basis?.length || !structure?.latticeVectors) return [];
  return symmetryLadder(conventionalCell(structure), structure.basis, tolMax);
}

/** Wyckoff label for an orbit: multiplicity + letter, or multiplicity + site symmetry. */
export function orbitLabel(orbit) {
  return `${orbit.size}${orbit.wyckoff || ` (${orbit.site})`}`;
}
