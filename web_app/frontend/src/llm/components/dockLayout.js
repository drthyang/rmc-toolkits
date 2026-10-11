// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Where the Copilot sits, and how big. Pure helpers for CopilotChat.
//
// - Docked (a window wide enough for the narrowest panel beside pages of
//   PAGE_MIN px — laptops, monitors, tablets in landscape — and taller than a
//   phone in landscape): a side panel beside the pages, resizable from its
//   left edge. The pages keep at least PAGE_MIN px, so they never fall below
//   the tablet layout; the host gives the panel its own grid column, so it
//   never covers a card. A landscape tablet is short: there a sheet would
//   leave a card and the conversation too little height each.
// - Sheet (tablets in portrait, phones, short windows): a bottom panel in
//   its own grid row, resizable from its top edge, as a fraction of the
//   window height.
// - Collapsed: a slim rail (docked) or bar (sheet) with an unread mark.
//
// The open state, the docked width and the sheet height persist in
// localStorage under DOCK_STORAGE_KEY.

export const DOCK_STORAGE_KEY = 'rmc-llm-dock-v1';
export const DOCK_MIN = 320;
export const DOCK_MAX = 760;
export const PAGE_MIN = 760;
export const DOCK_MIN_WINDOW = PAGE_MIN + DOCK_MIN;
export const DOCK_QUERY = `(min-width: ${DOCK_MIN_WINDOW}px) and (min-height: 541px)`;
export const DOCK_STEP = 24;
export const SHEET_MIN = 0.28;
export const SHEET_MAX = 0.72;
export const SHEET_DEFAULT = 0.46;
export const SHEET_STEP = 0.05;
// While a card is in focus on a sheet layout, the page keeps at least this
// share of the window above the sheet.
export const SHEET_FOCUS_MAX = 0.5;
// Air kept under a focused card, above the sheet (px).
const FOCUS_GAP = 16;

const finite = (value) => typeof value === 'number' && Number.isFinite(value);

/** The default docked width: 28 % of the window, 23–30rem (rem = root font size). */
export const defaultDockWidth = (viewportWidth, rem = 15) => {
    const scale = finite(rem) && rem > 0 ? rem : 15;
    const width = finite(viewportWidth) ? viewportWidth * 0.28 : 26 * scale;
    return Math.round(Math.min(30 * scale, Math.max(23 * scale, width)));
};

/** The widest dock that leaves the pages PAGE_MIN px (never below DOCK_MIN). */
export const maxDockWidth = (viewportWidth) => (
    finite(viewportWidth) ? Math.max(DOCK_MIN, Math.min(DOCK_MAX, viewportWidth - PAGE_MIN)) : DOCK_MAX
);

/** A docked width within [DOCK_MIN, maxDockWidth(viewport)], whole px. */
export const clampDockWidth = (width, viewportWidth) => {
    const upper = maxDockWidth(viewportWidth);
    const value = finite(width) ? width : defaultDockWidth(viewportWidth);
    return Math.round(Math.min(upper, Math.max(DOCK_MIN, value)));
};

/** A sheet height (fraction of the window) within [SHEET_MIN, SHEET_MAX]. */
export const clampSheet = (fraction) => {
    const value = finite(fraction) ? fraction : SHEET_DEFAULT;
    return Math.round(Math.min(SHEET_MAX, Math.max(SHEET_MIN, value)) * 1000) / 1000;
};

/** The persisted state, validated: `{ open, width, sheet }` (null width = default). */
export const readDockState = (storage, key = DOCK_STORAGE_KEY) => {
    const fallback = { open: false, width: null, sheet: SHEET_DEFAULT };
    try {
        const raw = storage?.getItem(key);
        if (!raw) return fallback;
        const parsed = JSON.parse(raw);
        if (!parsed || typeof parsed !== 'object') return fallback;
        return {
            open: parsed.open === true,
            width: finite(parsed.width) && parsed.width > 0 ? Math.round(parsed.width) : null,
            sheet: finite(parsed.sheet) ? clampSheet(parsed.sheet) : SHEET_DEFAULT
        };
    } catch {
        return fallback;
    }
};

export const writeDockState = (storage, patch, key = DOCK_STORAGE_KEY) => {
    try {
        const next = { ...readDockState(storage, key), ...patch };
        storage?.setItem(key, JSON.stringify(next));
        return next;
    } catch {
        // Private mode or a full quota: the layout still works, unremembered.
        return null;
    }
};

/**
 * A keyboard step on the resize handle: docked ←/→ widen/narrow the panel
 * (it grows to the left), the sheet ↑/↓ raise/lower it; Home/End go to the
 * limits. Returns the new value, or null for a key it does not handle.
 */
export const keyResize = (mode, key, value, viewportWidth) => {
    if (mode === 'docked') {
        const upper = maxDockWidth(viewportWidth);
        if (key === 'ArrowLeft') return clampDockWidth(value + DOCK_STEP, viewportWidth);
        if (key === 'ArrowRight') return clampDockWidth(value - DOCK_STEP, viewportWidth);
        if (key === 'Home') return DOCK_MIN;
        if (key === 'End') return upper;
        return null;
    }
    if (key === 'ArrowUp') return clampSheet(value + SHEET_STEP);
    if (key === 'ArrowDown') return clampSheet(value - SHEET_STEP);
    if (key === 'Home') return SHEET_MIN;
    if (key === 'End') return SHEET_MAX;
    return null;
};

/**
 * The sheet height that lets a focused card show whole above the sheet: the
 * card's height plus the space above its scroll area (`topOffset`: the
 * pages' header, or the pinned strip on a phone) must fit in the window above
 * the sheet. Never raises the sheet, never below SHEET_MIN (a card taller
 * than that shows its top), and at most SHEET_FOCUS_MAX while a card is in
 * focus.
 */
export const sheetForCard = ({ viewportHeight, topOffset = 0, cardHeight, current }) => {
    const now = clampSheet(current);
    if (!finite(viewportHeight) || viewportHeight <= 0 || !finite(cardHeight)) return Math.min(now, SHEET_FOCUS_MAX);
    const room = (viewportHeight - Math.max(0, topOffset) - cardHeight - FOCUS_GAP) / viewportHeight;
    return clampSheet(Math.min(now, SHEET_FOCUS_MAX, room));
};

/** The docked width a drag reaches: the panel grows as the pointer moves left. */
export const dragDockWidth = (startWidth, startX, x, viewportWidth) => (
    clampDockWidth(startWidth + (startX - x), viewportWidth)
);

/** The sheet height a drag reaches: the panel's top follows the pointer. */
export const dragSheet = (y, viewportHeight) => (
    finite(viewportHeight) && viewportHeight > 0 ? clampSheet((viewportHeight - y) / viewportHeight) : SHEET_DEFAULT
);
