// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Pure geometry/color helpers for the displacement-orientation sphere view:
// turn one orientation histogram (workers/orientation.js output) into the flat
// typed arrays a Three.js BufferGeometry consumes. Kept free of Three.js and
// DOM so the mapping is unit-testable (orientationSphere.test.js); the
// OrientationView component owns the scene and feeds these arrays in.

import { getLut, sampleColormap } from './colormaps';

// Clamp a 0..1 colormap coordinate.
const clamp01 = (t) => (t < 0 ? 0 : t > 1 ? 1 : t);

/**
 * Color coordinate (0..1) for one cell value. The scale runs 0 -> vmax so the
 * isotropic reference (enhancement = 1) keeps a stable position on the bar
 * whenever vmax is comparable, rather than being stretched per-frame.
 */
export const colorCoordinate = (value, vmax) => (vmax > 0 ? clamp01(value / vmax) : 0);

// Stable identity for a polygon vertex. Vertices are triangle circumcenters
// shared by (up to) three cells, but that identity is lost through the worker /
// JSON boundary, so key on rounded coordinates: bit-identical sources agree
// exactly, and 1e-9 rounding cannot merge distinct circumcenters (the closest
// pair at the maximum frequency is ~1e-2 apart).
const vertexKey = (v) => `${v[0].toFixed(9)},${v[1].toFixed(9)},${v[2].toFixed(9)}`;

/**
 * Per-cell radial scale factors for the amplitude relief: `1 + relief *
 * (cellMean/globalMean − 1)`, clamped, so a cell whose movers travel the
 * average distance stays on the unit sphere, farther-moving directions bulge
 * out and shorter ones dent in. Empty cells (mean 0) sit at the neutral 1.
 */
export const reliefFactors = (cellMeanAmplitude, globalMean, relief, clampLow = 0.3, clampHigh = 2.2) => (
    cellMeanAmplitude.map((mean) => {
        if (!(mean > 0) || !(globalMean > 0)) return 1;
        const factor = 1 + relief * (mean / globalMean - 1);
        return factor < clampLow ? clampLow : factor > clampHigh ? clampHigh : factor;
    })
);

/**
 * Per-vertex radius map for a relief surface. Each polygon vertex is shared by
 * up to three cells; averaging their relief factors at the shared vertex keeps
 * the surface continuous (crack-free) while every cell stays one flat color.
 * Returns a Map from vertexKey to radius.
 */
export const vertexRadii = (polygons, factors) => {
    const sums = new Map();
    polygons.forEach((polygon, cell) => {
        polygon.forEach((vertex) => {
            const key = vertexKey(vertex);
            const entry = sums.get(key);
            if (entry) { entry.total += factors[cell]; entry.count += 1; } else {
                sums.set(key, { total: factors[cell], count: 1 });
            }
        });
    });
    const radii = new Map();
    sums.forEach((entry, key) => { radii.set(key, entry.total / entry.count); });
    return radii;
};

/**
 * Triangle-fan mesh for every tiling cell, flat-colored by its value.
 *
 * `polygons` is the histogram's ragged cell-polygon list (5 or 6 unit vectors
 * per cell, CCW from outside); `values` one number per cell (enhancement).
 * Vertices are duplicated per cell so each cell is a single flat color -- the
 * crisp hex-bin look -- and `triangleCell` maps every triangle (in order) back
 * to its cell for raycast picking. `radii` (a vertexRadii Map) optionally
 * displaces each vertex radially for the amplitude relief.
 */
export const buildCellMesh = (polygons, values, colormap, vmax, radii = null) => {
    let triangles = 0;
    polygons.forEach((polygon) => { triangles += polygon.length - 2; });
    const positions = new Float32Array(triangles * 9);
    const colors = new Float32Array(triangles * 9);
    const triangleCell = new Int32Array(triangles);

    let t = 0;
    polygons.forEach((polygon, cell) => {
        const rgb = sampleColormap(colormap, colorCoordinate(values[cell], vmax));
        const r = rgb[0] / 255;
        const g = rgb[1] / 255;
        const b = rgb[2] / 255;
        for (let i = 1; i < polygon.length - 1; i += 1) {
            const base = t * 9;
            const fan = [polygon[0], polygon[i], polygon[i + 1]];
            for (let v = 0; v < 3; v += 1) {
                const scale = radii ? radii.get(vertexKey(fan[v])) ?? 1 : 1;
                positions[base + 3 * v] = fan[v][0] * scale;
                positions[base + 3 * v + 1] = fan[v][1] * scale;
                positions[base + 3 * v + 2] = fan[v][2] * scale;
                colors[base + 3 * v] = r;
                colors[base + 3 * v + 1] = g;
                colors[base + 3 * v + 2] = b;
            }
            triangleCell[t] = cell;
            t += 1;
        }
    });
    return { positions, colors, triangleCell, triangles };
};

/**
 * Line-segment positions tracing every cell boundary, radially inflated so the
 * outlines sit just above the mesh and never z-fight it. Interior edges are
 * shared by two cells and drawn twice -- harmless for a line overlay and far
 * simpler than deduplicating. `radii` matches buildCellMesh so the borders
 * follow the relief surface.
 */
export const buildCellOutline = (polygons, inflate = 1.002, radii = null) => {
    let edges = 0;
    polygons.forEach((polygon) => { edges += polygon.length; });
    const positions = new Float32Array(edges * 6);
    let e = 0;
    polygons.forEach((polygon) => {
        for (let i = 0; i < polygon.length; i += 1) {
            const a = polygon[i];
            const b = polygon[(i + 1) % polygon.length];
            const scaleA = (radii ? radii.get(vertexKey(a)) ?? 1 : 1) * inflate;
            const scaleB = (radii ? radii.get(vertexKey(b)) ?? 1 : 1) * inflate;
            const base = e * 6;
            positions[base] = a[0] * scaleA;
            positions[base + 1] = a[1] * scaleA;
            positions[base + 2] = a[2] * scaleA;
            positions[base + 3] = b[0] * scaleB;
            positions[base + 4] = b[1] * scaleB;
            positions[base + 5] = b[2] * scaleB;
            e += 1;
        }
    });
    return positions;
};

/**
 * CSS linear-gradient string sampling a colormap left (0) to right (vmax),
 * for the legend colorbar under the sphere.
 */
export const colorbarGradient = (colormap, stops = 24) => {
    const lut = getLut(colormap);
    const last = lut.length / 3 - 1;
    const parts = [];
    for (let s = 0; s <= stops; s += 1) {
        const index = Math.round((s / stops) * last) * 3;
        const percent = ((s / stops) * 100).toFixed(1);
        parts.push(`rgb(${lut[index]}, ${lut[index + 1]}, ${lut[index + 2]}) ${percent}%`);
    }
    return `linear-gradient(to right, ${parts.join(', ')})`;
};

/** Compact “[0.12, −0.99, 0.05]” form of a unit direction for readouts. */
export const formatDirection = (direction, digits = 2) => (
    `[${direction.map((value) => value.toFixed(digits)).join(', ')}]`
);
