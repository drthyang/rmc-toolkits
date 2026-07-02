// web_app/frontend/src/symmetryModel.js
//
// Glue between the parsed RMC structure (browserData.structureFromRmc6f) and the
// pure symmetry finder (symmetry.js): builds the conventional cell + basis and
// returns a space-group description + tolerance ladder for the UI.

import { spaceGroupAtTolerance, symmetryLadder, siteOrbits, wyckoffLetter } from './symmetry';

/** Conventional unit cell A_conv (rows, Å) = supercell lattice / supercell dims. */
export function conventionalCell(structure) {
  const { latticeVectors, supercell } = structure;
  return latticeVectors.map((row, i) => row.map((v) => v / Math.max(supercell[i], 1)));
}

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
  const orbits = siteOrbits(A, structure.basis, sg.ops, tol).map((o) => ({
    element: o.element,
    size: o.size,
    site: o.site,
    rep: o.rep,
    wyckoff: wyckoffLetter(sg.spaceGroupNumber, sg.centering, o.size, o.site, o.rep),
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
