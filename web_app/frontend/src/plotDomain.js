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

// A chart payload from GET /api/plot/data is usable only as an object with a
// `series` array. When a body is not valid JSON (the old Flask provider wrote
// bare NaN tokens), axios hands back the raw text as `response.data` — a
// string — and rendering it threw on the missing axis labels. Returns the
// message to show instead of a chart, or null for a usable payload.
export const plotPayloadError = (data) => {
    if (typeof data === 'string') return 'The server sent a plot payload that is not valid JSON';
    if (!data || typeof data !== 'object' || Array.isArray(data)) return 'The server sent no plot data';
    if (!Array.isArray(data.series)) return 'The plot payload has no series';
    return null;
};

// Index of the point nearest `target` along x among the points whose x AND y
// are finite, or -1 when there is none. Series may carry null/NaN gaps (a
// masked region arrives from Flask as JSON null); `null - target` would
// otherwise coerce to 0 and snap the hover onto a gap.
export const nearestFiniteIndex = (xs, ys, target) => {
    let best = -1;
    let bestDistance = Infinity;
    for (let index = 0; index < xs.length; index += 1) {
        const x = xs[index];
        const y = ys[index];
        if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
        const distance = Math.abs(x - target);
        if (distance < bestDistance) {
            bestDistance = distance;
            best = index;
        }
    }
    return best;
};
