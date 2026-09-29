// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// How a stat rail lays out stats that may not fit on one line (StatRail
// measures, ui.css draws). `runs` are the natural widths of the runs of stats
// between band starts — one run when the rail marks no bands, at most three —
// and `available` the width of the line.
//
// - Everything fits: { mode: null } — the plain one-line row.
// - No bands: { mode: 'wrapped' } — an aligned grid.
// - Bands: { mode: 'banded', breakA, breakB } — the runs packed into rows
//   greedily, a row broken before the first / second band start only when
//   that band does not fit beside the previous one.
export const planStatRows = (runs, available) => {
    const needed = runs.reduce((sum, width) => sum + width, 0);
    if (needed <= available) return { mode: null, breakA: false, breakB: false };
    if (runs.length < 2) return { mode: 'wrapped', breakA: false, breakB: false };
    const breaks = [false, false];
    let row = runs[0];
    runs.slice(1, 3).forEach((width, index) => {
        if (row + width <= available) {
            row += width;
        } else {
            breaks[index] = true;
            row = width;
        }
    });
    return { mode: 'banded', breakA: breaks[0], breakB: breaks[1] };
};
