// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Axis range for a set of data values, padded by 5% so the extremes are not
// drawn on the frame. Non-finite values are skipped; a series with no finite
// values at all falls back to [0, 1].
//
// One pass instead of Math.min(...values): spreading a series onto the argument
// stack throws RangeError past ~10^5 points, and RMCProfile outputs get that long.
export const niceDomain = (values) => {
    let min = Infinity;
    let max = -Infinity;
    for (const value of values) {
        if (!Number.isFinite(value)) continue;
        if (value < min) min = value;
        if (value > max) max = value;
    }
    if (min > max) return [0, 1];
    if (min === max) {
        min -= 1;
        max += 1;
    }
    const pad = (max - min) * 0.05;
    return [min - pad, max + pad];
};
