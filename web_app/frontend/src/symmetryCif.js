// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/symmetryCif.js
//
// A tolerance-ladder brick of the Detected SG card → its CIF file: the glue between the
// card (ModelSummary.jsx), the symmetry-averaged structure (averageStructure.js) and the
// CIF text (cifWriter.js).

import { symmetryAveragedStructure } from './averageStructure.js';
import { writeCif } from './cifWriter.js';
import { sanitizeFilename } from './figureExport.js';

/**
 * The tolerance (Å) a brick's CIF is computed at: its lower edge, the tightest tolerance at
 * which the group holds (plus 1e-6 Å, as the finder's orbit test is strict). Orbits and site
 * symmetries are then exactly what the group needs — a looser value, such as the brick's
 * middle that the card selects, could merge sites that are merely close (a split site).
 */
export const brickTolerance = (brick) => brick.from + 1e-6;

/**
 * The CIF of `structure` (browserData.structureFromRmc6f) in the group of `brick`
 * (symmetryModel.toleranceLadder). Throws with a user-facing message when there is nothing
 * to export.
 * @returns {{ filename:string, text:string }}
 */
export function brickCif(structure, brick, { date = new Date().toISOString().slice(0, 10) } = {}) {
  const model = symmetryAveragedStructure(structure, brickTolerance(brick));
  const stem = (structure.source?.split('/').pop() || 'structure').replace(/\.rmc6f$/i, '');
  const name = sanitizeFilename(`${stem}_${brick.spaceGroup}`);
  const text = writeCif(model, { dataName: name, date, tolRange: [brick.from, brick.to] });
  return { filename: `${name}.cif`, text };
}
