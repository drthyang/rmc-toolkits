// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Series colors for InteractivePlot, in the order they are handed out. Its own
// module rather than an export from the component: a component file that also
// exports constants breaks Fast Refresh (react-refresh/only-export-components).
//
// Callers can use it to color a guide line to match the data series it belongs
// to — guide series never consume a palette slot, so the Nth data series in a
// plot is PLOT_PALETTE[N].
export const PLOT_PALETTE = [
    '#1f6fd6',
    '#e8590c',
    '#099268',
    '#d6336c',
    '#6741d9',
    '#66a80f',
    '#0c8599',
    '#e67700'
];

// Neutral stroke for a guide that is not tied to a particular series.
export const GUIDE_STROKE = '#98a2b3';

// Bond-role colors on the Bond Geometry page: the A–B bond and, when the two
// bonds differ, the B–C bond. The same pair marks the 3D bond sticks, the
// split window guides and window chips, and the chips' bond dashes — and they
// are the first two plot colors, so a partial g(r) curve (series N is
// PLOT_PALETTE[N]) wears the color of the bond it brackets.
export const BOND_COLORS = { ab: PLOT_PALETTE[0], bc: PLOT_PALETTE[1] };
