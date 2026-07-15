// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Robust field extraction for one RMCProfile `.rmc6f` atom line, tolerant of the
// current format and of older files that drop columns:
//   full:    id  element  [type]  x  y  z  ref  cellx  celly  cellz   (10 fields)
//   full:    id  element          x  y  z  ref  cellx  celly  cellz   ( 9 fields)
//   coords:  id  element  [label]  x  y  z                            ( 5-6 fields)
// In the full forms the reference number and the three cell indices are always the
// last four fields and the fractional coordinates the three before them, so indexing
// from the END tolerates any number of label columns between element and coordinates.
// The oldest files carry neither a site (reference) column nor per-atom cell indices,
// so the three fractional coordinates are simply the final fields; `referenceNumber`
// and `cellIndices` come back null and the site/cell membership is reconstructed
// downstream by folding every atom into one unit cell and clustering.
//
// `parts` is the whitespace-split, non-empty tokens of a line. Returns null for a
// non-atom / malformed line (too few fields, or non-numeric coordinates) so callers
// can skip it.
export const parseAtomLine = (parts) => {
    const n = parts.length;
    // Minimum: id, element, x, y, z.
    if (n < 5) return null;
    const element = parts[1];
    // Coords-only short form ( id element [label] x y z, 5-6 fields): no site or cell
    // columns, so the fractional coordinates are the final three fields and a leading
    // `[label]` is ignored. Restricting this to short lines keeps a truncated full
    // line (7-8 fields) from being misread as coordinates.
    if (n <= 6) {
        const coords = [Number(parts[n - 3]), Number(parts[n - 2]), Number(parts[n - 1])];
        if (!coords.every(Number.isFinite)) return null;
        return { element, coords, referenceNumber: null, cellIndices: null };
    }
    // Full format ( >= 9 fields): ref + 3 cell indices are the last four fields and the
    // fractional coordinates the three before them. 7-8 field lines are malformed.
    if (n < 9) return null;
    const coords = [Number(parts[n - 7]), Number(parts[n - 6]), Number(parts[n - 5])];
    const referenceNumber = Number(parts[n - 4]);
    const cellIndices = [Number(parts[n - 3]), Number(parts[n - 2]), Number(parts[n - 1])];
    if (!Number.isFinite(referenceNumber)) return null;
    if (!coords.every(Number.isFinite) || !cellIndices.every(Number.isFinite)) return null;
    return { element, coords, referenceNumber, cellIndices };
};
