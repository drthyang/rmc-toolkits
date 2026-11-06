// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Where a pointer readout sits beside the point it reads (CSS px inside its
// stage): right of and below the point, or left of it past the stage's middle,
// so a readout near the right edge (a tap on a phone) stays on the card
// instead of running off it.
export const TOOLTIP_OFFSET = { x: 14, y: 12 };

export const tooltipPlace = ({ x, y, w }) => (
    Number.isFinite(w) && x > w / 2
        ? { right: w - x + TOOLTIP_OFFSET.x, top: y + TOOLTIP_OFFSET.y }
        : { left: x + TOOLTIP_OFFSET.x, top: y + TOOLTIP_OFFSET.y }
);
