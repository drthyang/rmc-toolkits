// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The AI Copilot's handle on the Workbench's 3D views: crystal-direction math
// (pure, no Three.js) and a registry the Three.js panels join while mounted.
//
// Conventions (the same as the panels' own axis snaps):
// - "Looking down [u v w]" puts the camera on the +[u v w] side of the target,
//   so that direction points out of the screen, toward the viewer.
// - [u v w] is a direct-lattice direction: u·a + v·b + w·c in the panels'
//   Cartesian frame (the unit-cell rows), so it is honest in an oblique cell.
// - Screen-up is the first of c, b, a that is not within ~18° of the line of
//   sight, with its component along the view removed (down a or b: c up; down
//   c: b up), as OrientationView's and PcaKdePage's snaps do.
//
// A panel registers a controller (useCopilotViewer.js) under its card id:
//   {
//     ready(): bool                     drawn, with data
//     cell(): [a, b, c] | null          unit-cell rows (Å, Cartesian)
//     describe(): { direction, up, … }  the current view, Cartesian
//     lookAlong(direction, up)          Cartesian unit vectors
//     reset()
//     centerOn?(ref), highlightSites?(refs), clearHighlight?(),
//     selectSite?(ref), sites?(): [{ ref, element }], pcaAxes?(): rows | null
//   }

import { invert3, matVec3, transpose3 } from './pcaCrystalFrame';

const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const norm = (v) => Math.hypot(v[0], v[1], v[2]);
const scale = (v, s) => [v[0] * s, v[1] * s, v[2] * s];
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const normalize = (v) => {
    const length = norm(v);
    return length > 1e-12 ? scale(v, 1 / length) : null;
};
const finiteVector = (v) => Array.isArray(v) && v.length === 3 && v.every(Number.isFinite);

export const angleDeg = (a, b) => {
    const na = normalize(a);
    const nb = normalize(b);
    if (!na || !nb) return null;
    return Math.acos(Math.max(-1, Math.min(1, dot(na, nb)))) * 180 / Math.PI;
};

const AXES = { a: 0, b: 1, c: 2 };
const CARTESIAN = { x: [1, 0, 0], y: [0, 1, 0], z: [0, 0, 1] };

/**
 * A model's spelling of a view direction → one of
 *   { uvw: [u, v, w] }   "[111]", "[1 -1 0]", "1,1,0", "[1-10]", "[1̄10]", "<111>", "[0.71 0.71 1]"
 *   { axis: 0|1|2 }      "a", "b", "c" (also "a axis", "along c")
 *   { pc: 0|1|2 }        "PC1", "PC2", "PC3"
 *   { cartesian: [..] }  "x", "y", "z"
 * or null. Without separators each digit is one index ("[110]" = 1 1 0),
 * a minus or a combining overline negating its digit.
 */
export const parseDirection = (raw) => {
    if (Array.isArray(raw)) return raw.length === 3 && raw.every((v) => Number.isFinite(Number(v))) ? { uvw: raw.map(Number) } : null;
    if (typeof raw !== 'string') return null;
    let text = raw.trim().toLowerCase().replace(/^(along|down|the)\s+/, '').replace(/[\s-]*(axis|direction)$/, '').trim();
    if (text in AXES) return { axis: AXES[text] };
    if (text in CARTESIAN) return { cartesian: CARTESIAN[text] };
    const pc = /^pc\s*([123])$/.exec(text);
    if (pc) return { pc: Number(pc[1]) - 1 };
    text = text.replace(/^[[<(]\s*|\s*[\]>)]$/g, '').replace(/[−–]/g, '-');
    if (!text) return null;
    // Overlined digits (1̄ = 1 + U+0304/U+0305) are negative.
    text = text.replace(/(\d)[̄̅]/g, '-$1');
    let parts;
    if (/[\s,;]/.test(text.trim())) {
        parts = text.split(/[\s,;]+/).filter(Boolean).map(Number);
    } else {
        const tokens = text.match(/-?\d(?:\.\d+)?/g);
        if (!tokens || tokens.join('') !== text.replace(/\+/g, '')) return null;
        parts = tokens.map(Number);
    }
    if (parts.length !== 3 || !parts.every(Number.isFinite) || parts.every((value) => value === 0)) return null;
    return { uvw: parts };
};

/** "[1 1 1]", "[1 -1 0]", "[0.71 0.71 1]" — a direction as the model reads it. */
export const directionLabel = (uvw) => (finiteVector(uvw)
    ? `[${uvw.map((value) => {
        const rounded = Math.abs(value - Math.round(value)) < 1e-9 ? Math.round(value) : Number(value.toFixed(2));
        return Object.is(rounded, -0) ? 0 : rounded;
    }).join(' ')}]`
    : null);

/** u·a + v·b + w·c (Cartesian, Å), or null. */
export const uvwToCartesian = (uvw, cell) => {
    if (!finiteVector(uvw) || !Array.isArray(cell) || cell.length !== 3) return null;
    return [0, 1, 2].map((axis) => uvw[0] * cell[0][axis] + uvw[1] * cell[1][axis] + uvw[2] * cell[2][axis]);
};

/** The direct-lattice components [u v w] of a Cartesian direction (M⁻ᵀ·d). */
export const cartesianToUvw = (vector, cell) => {
    if (!finiteVector(vector) || !Array.isArray(cell)) return null;
    const inverse = invert3(cell);
    return inverse ? matVec3(transpose3(inverse), vector) : null;
};

/**
 * The simplest integer [u v w] (largest index ≤ maxIndex) within toleranceDeg
 * of a Cartesian direction: `{ uvw, offDeg }`, else the direction's own
 * components scaled so the largest is ±1 (`{ uvw, offDeg: 0, exact: false }`).
 */
export const nearestLowIndex = (vector, cell, { maxIndex = 4, toleranceDeg = 2 } = {}) => {
    const raw = cartesianToUvw(vector, cell);
    if (!raw) return null;
    let best = null;
    for (let u = -maxIndex; u <= maxIndex; u += 1) {
        for (let v = -maxIndex; v <= maxIndex; v += 1) {
            for (let w = -maxIndex; w <= maxIndex; w += 1) {
                if (!u && !v && !w) continue;
                const candidate = uvwToCartesian([u, v, w], cell);
                const off = angleDeg(candidate, vector);
                if (off == null || off > toleranceDeg) continue;
                const size = Math.abs(u) + Math.abs(v) + Math.abs(w);
                if (!best || size < best.size - 1e-9 || (Math.abs(size - best.size) < 1e-9 && off < best.offDeg)) {
                    best = { uvw: [u, v, w], offDeg: off, size };
                }
            }
        }
    }
    if (best) return { uvw: best.uvw, offDeg: Number(best.offDeg.toFixed(2)), exact: true };
    const largest = Math.max(...raw.map(Math.abs)) || 1;
    return { uvw: raw.map((value) => Number((value / largest).toFixed(2)) + 0), offDeg: 0, exact: false };
};

// Screen-up within this angle of the line of sight is too steep to use.
const UP_LIMIT_COS = Math.cos((18 * Math.PI) / 180);

/**
 * Screen-up for a view along `direction` (Cartesian): the first of c, b, a
 * (or z, y, x without a cell) not within ~18° of the line of sight, with its
 * component along the view removed. Unit vector.
 */
export const viewUp = (direction, cell = null) => {
    const view = normalize(direction);
    if (!view) return [0, 1, 0];
    const candidates = cell ? [cell[2], cell[1], cell[0]] : [[0, 0, 1], [0, 1, 0], [1, 0, 0]];
    for (const candidate of [...candidates, [0, 1, 0], [1, 0, 0]]) {
        const unit = normalize(candidate);
        if (!unit || Math.abs(dot(unit, view)) > UP_LIMIT_COS) continue;
        const projected = normalize(sub(unit, scale(view, dot(unit, view))));
        if (projected) return projected;
    }
    return [0, 1, 0];
};

/** The Cartesian line of sight a parsed direction names, for a cell (and PC axes). */
export const directionVector = (parsed, cell, pcaAxes = null) => {
    if (!parsed) return null;
    if (parsed.uvw) return cell ? normalize(uvwToCartesian(parsed.uvw, cell)) : null;
    if (parsed.axis != null) return cell ? normalize(cell[parsed.axis]) : null;
    if (parsed.pc != null) return pcaAxes ? normalize(pcaAxes[parsed.pc]) : null;
    if (parsed.cartesian) return normalize(parsed.cartesian);
    return null;
};

/** Camera placement for a view: target + direction·distance, with `up`. */
export const posePlan = ({ target, direction, distance, up }) => {
    const view = normalize(direction);
    if (!view || !finiteVector(target) || !Number.isFinite(distance) || distance <= 0) return null;
    return {
        target: [...target],
        position: [0, 1, 2].map((axis) => target[axis] + view[axis] * distance),
        direction: view,
        up: normalize(up) || viewUp(view)
    };
};

/** A view as the model reads it: the low-index direction and screen-up. */
export const describeDirection = (direction, up, cell) => {
    if (!cell) return { view_along_cartesian: (normalize(direction) || []).map((v) => Number(v.toFixed(3)) + 0) };
    const along = nearestLowIndex(direction, cell);
    const upward = up ? nearestLowIndex(up, cell, { toleranceDeg: 5 }) : null;
    const described = {
        view_along: along ? `${along.exact ? '' : '≈'}${directionLabel(along.uvw)}` : null
    };
    if (along?.exact && along.offDeg > 0.05) described.off_by_deg = along.offDeg;
    if (upward) described.screen_up = `${upward.exact ? '' : '≈'}${directionLabel(upward.uvw)}`;
    return described;
};

// --- the registry ---------------------------------------------------------------

const viewers = new Map();

/** Register a mounted 3D view under its card id; returns the unregister function. */
export const registerViewer = (id, controller) => {
    viewers.set(id, controller);
    return () => {
        if (viewers.get(id) === controller) viewers.delete(id);
    };
};

export const getViewer = (id) => viewers.get(id) || null;

export const registeredViewers = () => [...viewers.keys()];

const abortError = () => new DOMException('Aborted', 'AbortError');

/**
 * The controller for `id` once it is registered (its page mounted) and, with
 * `ready`, once it has drawn its data. Resolves null after `timeoutMs`.
 */
export const waitForViewer = (id, { timeoutMs = 5000, signal = null, ready = true, pollMs = 80 } = {}) => new Promise((resolve, reject) => {
    const started = Date.now();
    const poll = () => {
        if (signal?.aborted) return reject(abortError());
        const viewer = getViewer(id);
        if (viewer && (!ready || viewer.ready())) return resolve(viewer);
        if (Date.now() - started >= timeoutMs) return resolve(null);
        setTimeout(poll, pollMs);
    };
    poll();
});
