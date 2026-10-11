// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Which edge an InfoBadge popover aligns to, so the open popover stays inside
// the window: the requested edge when it fits, else the other one when that
// fits, else the requested one. Phones show every popover as a sheet across
// the screen (ui.css), so there the requested edge is kept.

const SHEET_QUERY = '(max-width: 760px)';
// The popover's widest box (ui.css max-width) and the gap it keeps from the
// window's edge, in rem.
const POPOVER_MAX_REM = 20;
const EDGE_GAP_REM = 0.75;

/** Pure: the edge for a trigger spanning [left, right] in a window this wide. */
export const popoverEdge = ({ align, left, right, viewportWidth, popoverWidth, gap }) => {
    const fitsStart = left + popoverWidth <= viewportWidth - gap;
    const fitsEnd = right - popoverWidth >= gap;
    if (align === 'end') return fitsEnd || !fitsStart ? 'end' : 'start';
    return fitsStart || !fitsEnd ? 'start' : 'end';
};

/** The edge for a badge element as it sits in the window now. */
export const measurePopoverEdge = (element, align) => {
    if (typeof window === 'undefined' || !element) return align;
    if (window.matchMedia?.(SHEET_QUERY).matches) return align;
    const rect = element.getBoundingClientRect();
    const viewportWidth = document.documentElement.clientWidth || window.innerWidth;
    const rem = parseFloat(window.getComputedStyle(document.documentElement).fontSize) || 15;
    return popoverEdge({
        align,
        left: rect.left,
        right: rect.right,
        viewportWidth,
        popoverWidth: Math.min(POPOVER_MAX_REM * rem, viewportWidth),
        gap: EDGE_GAP_REM * rem
    });
};
