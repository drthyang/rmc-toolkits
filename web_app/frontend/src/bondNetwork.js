// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The bonds the Bond Geometry page's folded cell draws (FoldedCellPanel), as
// plain numbers — no Three.js, so it is unit-tested. A bond set is the
// computed triplet's A–B (or B–C) window: from every average site of the
// first element in the unit cell to every site of the second within the
// window, periodic images one cell either way included (a bond may reach an
// image just outside the box: that is the real coordination). Same-element
// in-cell pairs count once.
//
// The search is a linked-cell one over the images (bins at least the window's
// maximum wide), so a large cell — a P1 box where every atom is its own
// site — stays linear in the sites instead of quadratic. What is drawn is
// capped (MAX_DRAWN_BONDS): above it the bonds nearest the cell centre are
// kept, a readable cluster rather than a solid block, and the panel says how
// many of how many it shows.

export const MAX_DRAWN_BONDS = 1500;
// The cell is split into at most this many bins along each axis, whatever
// the window: a tiny window in a huge cell would otherwise ask for millions.
const MAX_BINS_PER_AXIS = 64;

const IMAGE_SHIFTS = [];
for (let mx = -1; mx <= 1; mx += 1) {
    for (let my = -1; my <= 1; my += 1) {
        for (let mz = -1; mz <= 1; mz += 1) IMAGE_SHIFTS.push([mx, my, mz]);
    }
}

const toCartesian = (unit, fx, fy, fz, out, offset) => {
    out[offset] = fx * unit[0][0] + fy * unit[1][0] + fz * unit[2][0];
    out[offset + 1] = fx * unit[0][1] + fy * unit[1][1] + fz * unit[2][1];
    out[offset + 2] = fx * unit[0][2] + fy * unit[1][2] + fz * unit[2][2];
};

/**
 * The bonds of one set between the average sites of one unit cell.
 *
 * @param {Array<{referenceNumber, element, siteFractional}>} sites
 * @param {number[][]} unit   unit-cell vectors (rows, Cartesian Å)
 * @param {{elements: [string, string], window: [number, number]}} set
 * @returns {Float64Array} segments, six numbers per bond: from xyz, to xyz
 */
export const findSiteBonds = (sites, unit, set) => {
    const [elementA, elementB] = set.elements;
    const [rmin, rmax] = set.window;
    if (!(rmax > 0) || !(rmax >= rmin)) return new Float64Array(0);
    const sameElement = elementA === elementB;
    const from = sites.filter((site) => site.element === elementA);
    const to = sites.filter((site) => site.element === elementB);
    if (!from.length || !to.length) return new Float64Array(0);

    // Every image of every B site, in Cartesian.
    const count = to.length * IMAGE_SHIFTS.length;
    const points = new Float64Array(count * 3);
    const owner = new Int32Array(count);
    const inCell = new Uint8Array(count);
    let n = 0;
    to.forEach((site, index) => {
        const [fx, fy, fz] = site.siteFractional;
        IMAGE_SHIFTS.forEach(([mx, my, mz]) => {
            toCartesian(unit, fx + mx, fy + my, fz + mz, points, n * 3);
            owner[n] = index;
            inCell[n] = mx === 0 && my === 0 && mz === 0 ? 1 : 0;
            n += 1;
        });
    });

    // Bin them (CSR): bins at least rmax wide, so a bond is in a
    // neighbouring bin.
    const low = [Infinity, Infinity, Infinity];
    const high = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < count; i += 1) {
        for (let axis = 0; axis < 3; axis += 1) {
            const value = points[i * 3 + axis];
            if (value < low[axis]) low[axis] = value;
            if (value > high[axis]) high[axis] = value;
        }
    }
    const extent = Math.max(high[0] - low[0], high[1] - low[1], high[2] - low[2], 1e-9);
    const size = Math.max(rmax, extent / MAX_BINS_PER_AXIS);
    const dims = [0, 1, 2].map((axis) => Math.floor((high[axis] - low[axis]) / size) + 1);
    const binOf = (x, y, z) => {
        const ix = Math.floor((x - low[0]) / size);
        const iy = Math.floor((y - low[1]) / size);
        const iz = Math.floor((z - low[2]) / size);
        return [ix, iy, iz];
    };
    const flat = (ix, iy, iz) => (ix * dims[1] + iy) * dims[2] + iz;
    const starts = new Int32Array(dims[0] * dims[1] * dims[2] + 1);
    const binIndex = new Int32Array(count);
    for (let i = 0; i < count; i += 1) {
        const [ix, iy, iz] = binOf(points[i * 3], points[i * 3 + 1], points[i * 3 + 2]);
        binIndex[i] = flat(ix, iy, iz);
        starts[binIndex[i] + 1] += 1;
    }
    for (let b = 1; b < starts.length; b += 1) starts[b] += starts[b - 1];
    const fill = starts.slice(0, -1);
    const order = new Int32Array(count);
    for (let i = 0; i < count; i += 1) order[fill[binIndex[i]]++] = i;

    const out = [];
    const origin = new Float64Array(3);
    from.forEach((site) => {
        const [fx, fy, fz] = site.siteFractional;
        toCartesian(unit, fx, fy, fz, origin, 0);
        const [cx, cy, cz] = binOf(origin[0], origin[1], origin[2]);
        for (let ix = Math.max(0, cx - 1); ix <= Math.min(dims[0] - 1, cx + 1); ix += 1) {
            for (let iy = Math.max(0, cy - 1); iy <= Math.min(dims[1] - 1, cy + 1); iy += 1) {
                for (let iz = Math.max(0, cz - 1); iz <= Math.min(dims[2] - 1, cz + 1); iz += 1) {
                    const bin = flat(ix, iy, iz);
                    for (let k = starts[bin]; k < starts[bin + 1]; k += 1) {
                        const j = order[k];
                        const other = to[owner[j]];
                        if (inCell[j]) {
                            if (other.referenceNumber === site.referenceNumber) continue;
                            // Same-element in-cell pairs would double-draw.
                            if (sameElement && other.referenceNumber < site.referenceNumber) continue;
                        }
                        const dx = points[j * 3] - origin[0];
                        const dy = points[j * 3 + 1] - origin[1];
                        const dz = points[j * 3 + 2] - origin[2];
                        const length = Math.hypot(dx, dy, dz);
                        if (length < rmin || length > rmax || length === 0) continue;
                        out.push(origin[0], origin[1], origin[2], points[j * 3], points[j * 3 + 1], points[j * 3 + 2]);
                    }
                }
            }
        }
    });
    return Float64Array.from(out);
};

/**
 * Every set's bonds, capped at `maxBonds` in all: above the cap the bonds
 * whose midpoints lie nearest `center` are kept (ties by order, so the
 * choice is deterministic).
 *
 * @returns {{ sets: Array<{color, elements, window, segments: Float64Array}>,
 *             shown: number, total: number, limited: boolean }}
 */
export const buildBondNetwork = (sites, unit, bondSets, { center = [0, 0, 0], maxBonds = MAX_DRAWN_BONDS } = {}) => {
    const found = (bondSets ?? []).map((set) => ({ ...set, segments: findSiteBonds(sites ?? [], unit, set) }));
    const total = found.reduce((sum, set) => sum + set.segments.length / 6, 0);
    if (total <= maxBonds) return { sets: found, shown: total, total, limited: false };

    const ranked = [];
    found.forEach((set, setIndex) => {
        const s = set.segments;
        for (let b = 0; b < s.length / 6; b += 1) {
            const mx = (s[b * 6] + s[b * 6 + 3]) / 2 - center[0];
            const my = (s[b * 6 + 1] + s[b * 6 + 4]) / 2 - center[1];
            const mz = (s[b * 6 + 2] + s[b * 6 + 5]) / 2 - center[2];
            ranked.push({ setIndex, b, d: mx * mx + my * my + mz * mz });
        }
    });
    ranked.sort((x, y) => x.d - y.d);
    const keep = found.map(() => []);
    ranked.slice(0, maxBonds).forEach(({ setIndex, b }) => keep[setIndex].push(b));
    const sets = found.map((set, setIndex) => {
        const picked = keep[setIndex].sort((x, y) => x - y);
        const segments = new Float64Array(picked.length * 6);
        picked.forEach((b, i) => segments.set(set.segments.subarray(b * 6, b * 6 + 6), i * 6));
        return { ...set, segments };
    });
    return { sets, shown: maxBonds, total, limited: true };
};

/**
 * The drawn radius of a bond set's tubes (Å): a fraction of the bond length,
 * with a floor tied to the cell so a large cell's bonds stay visible when the
 * whole cell is in view.
 */
export const bondRadius = (window, cellSize) => {
    const mid = (window[0] + window[1]) / 2;
    return Math.max(0.035 * mid, 0.0025 * cellSize);
};

/**
 * The cloud's opacity for an element: full (0.55) with no bonds drawn or for
 * an element of a drawn bond, faded otherwise, so the selected triplet's
 * atoms and bonds stand out.
 */
export const CLOUD_OPACITY = 0.55;
export const CLOUD_OPACITY_FADED = 0.15;
export const cloudOpacity = (element, bondSets) => (
    !bondSets?.length || bondSets.some((set) => set.elements.includes(element))
        ? CLOUD_OPACITY
        : CLOUD_OPACITY_FADED
);
