// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Whether a structure read looks like a configuration RMCProfile is still
// writing — used by Live Data to keep the previous complete model summary.
//
// A file caught mid-write is short of its header's `Number of atoms:` by the
// lines not yet written (and its last line may be cut, which reads as one
// unparsed line). A COMPLETE file whose atom has blown up to NaN/Inf/****
// has every line present: its non-finite lines are skipped by the parser but
// they are not missing, so such a read is shown with its parse warning rather
// than held back as "still being written" (review: parsers, Dashboard guard).
export const isIncompleteStructure = (structure) => {
    const report = structure?.parseReport;
    if (!report || !Number.isFinite(report.declaredAtoms)) return false;
    const present = (report.parsedAtoms || 0) + (report.coordsOnlyAtoms || 0) + (report.nonFiniteLines || 0);
    return present < report.declaredAtoms;
};
