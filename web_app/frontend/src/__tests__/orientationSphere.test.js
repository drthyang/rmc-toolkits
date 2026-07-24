// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, it, expect } from 'vitest';

import {
    buildCellMesh,
    buildCellOutline,
    colorCoordinate,
    colorbarGradient,
    formatDirection,
    reliefFactors,
    vertexRadii
} from '../orientationSphere.js';
import { orientationHistogram, goldbergTiling } from '../workers/orientation.js';
import { sampleColormap } from '../colormaps.js';

// A real histogram (with geometry) drives the mesh tests, so the helper is
// exercised against the exact ragged-polygon shape the engine emits.
const makeResult = () => {
    const points = [];
    let seed = 7;
    const rand = () => {
        seed = (seed * 1103515245 + 12345) % 2147483648;
        return seed / 2147483648 - 0.5;
    };
    for (let i = 0; i < 3000; i += 1) points.push([rand(), rand(), rand()]);
    return orientationHistogram(points, { frequency: 3, geometry: true });
};

describe('buildCellMesh', () => {
    const result = makeResult();
    const mesh = buildCellMesh(result.polygons, result.enhancement, 'viridis', result.vmax);

    it('emits one fan triangle per polygon vertex minus two', () => {
        const expected = result.sizes.reduce((sum, size) => sum + (size - 2), 0);
        expect(mesh.triangles).toBe(expected);
        expect(mesh.positions).toHaveLength(expected * 9);
        expect(mesh.colors).toHaveLength(expected * 9);
        expect(mesh.triangleCell).toHaveLength(expected);
    });

    it('keeps every vertex on the unit sphere', () => {
        for (let v = 0; v < mesh.positions.length; v += 3) {
            const norm = Math.hypot(mesh.positions[v], mesh.positions[v + 1], mesh.positions[v + 2]);
            expect(norm).toBeCloseTo(1, 5);
        }
    });

    it('flat-colors each cell (all its vertices share one color)', () => {
        const byCell = new Map();
        for (let t = 0; t < mesh.triangles; t += 1) {
            const cell = mesh.triangleCell[t];
            for (let v = 0; v < 3; v += 1) {
                const base = t * 9 + v * 3;
                const key = `${mesh.colors[base]},${mesh.colors[base + 1]},${mesh.colors[base + 2]}`;
                if (!byCell.has(cell)) byCell.set(cell, key);
                expect(key).toBe(byCell.get(cell));
            }
        }
        expect(byCell.size).toBe(result.cellCount);
    });

    it('maps the peak cell to the top of the colormap', () => {
        const peakTriangle = mesh.triangleCell.indexOf(result.peakCell);
        const base = peakTriangle * 9;
        const top = sampleColormap('viridis', colorCoordinate(result.vmax, result.vmax));
        expect(mesh.colors[base]).toBeCloseTo(top[0] / 255, 5);
        expect(mesh.colors[base + 1]).toBeCloseTo(top[1] / 255, 5);
        expect(mesh.colors[base + 2]).toBeCloseTo(top[2] / 255, 5);
    });

    it('covers every cell in the triangle -> cell map', () => {
        expect(new Set(mesh.triangleCell).size).toBe(result.cellCount);
    });
});

describe('buildCellOutline', () => {
    it('traces every polygon edge, inflated off the surface', () => {
        const tiling = goldbergTiling(2);
        const polygons = tiling.polygons.map((polygon, cell) => polygon.slice(0, tiling.sizes[cell]));
        const outline = buildCellOutline(polygons, 1.01);
        const edges = tiling.sizes.reduce((sum, size) => sum + size, 0);
        expect(outline).toHaveLength(edges * 6);
        for (let v = 0; v < outline.length; v += 3) {
            const norm = Math.hypot(outline[v], outline[v + 1], outline[v + 2]);
            expect(norm).toBeCloseTo(1.01, 5);
        }
    });
});

describe('amplitude relief', () => {
    it('maps mean amplitudes to clamped radial factors with a neutral empty cell', () => {
        const factors = reliefFactors([0.2, 0.1, 0, 0.9], 0.2, 1);
        expect(factors[0]).toBeCloseTo(1, 12);      // at the site average
        expect(factors[1]).toBeCloseTo(0.5, 12);    // half the average travel
        expect(factors[2]).toBe(1);                 // empty cell stays neutral
        expect(factors[3]).toBe(2.2);               // clamped at the ceiling
        // relief = 0 flattens everything back onto the unit sphere.
        reliefFactors([0.2, 0.1, 0, 0.9], 0.2, 0).forEach((factor) => expect(factor).toBe(1));
    });

    it('produces a crack-free relief surface (shared vertices share one radius)', () => {
        const tiling = goldbergTiling(3);
        const polygons = tiling.polygons.map((polygon, cell) => polygon.slice(0, tiling.sizes[cell]));
        // Strongly varying per-cell factors, so a per-cell (non-averaged)
        // radius would guarantee mismatches at shared vertices.
        const factors = polygons.map((_, cell) => 0.5 + (cell % 7) * 0.2);
        const radii = vertexRadii(polygons, factors);
        const mesh = buildCellMesh(polygons, factors, 'viridis', 2, radii);
        // Group every emitted vertex by its (pre-scale) direction; all copies
        // across cells must land at the identical radius.
        const byDirection = new Map();
        for (let v = 0; v < mesh.positions.length; v += 3) {
            const norm = Math.hypot(mesh.positions[v], mesh.positions[v + 1], mesh.positions[v + 2]);
            const key = [mesh.positions[v] / norm, mesh.positions[v + 1] / norm, mesh.positions[v + 2] / norm]
                .map((value) => value.toFixed(6)).join(',');
            if (!byDirection.has(key)) byDirection.set(key, norm);
            expect(norm).toBeCloseTo(byDirection.get(key), 9);
        }
        // And the radii genuinely vary (this is a relief, not a sphere).
        const norms = [...byDirection.values()];
        expect(Math.max(...norms) - Math.min(...norms)).toBeGreaterThan(0.2);
    });

    it('keeps the outline on the relief surface', () => {
        const tiling = goldbergTiling(2);
        const polygons = tiling.polygons.map((polygon, cell) => polygon.slice(0, tiling.sizes[cell]));
        const factors = polygons.map(() => 1.5);
        const radii = vertexRadii(polygons, factors);
        const outline = buildCellOutline(polygons, 1.01, radii);
        for (let v = 0; v < outline.length; v += 3) {
            const norm = Math.hypot(outline[v], outline[v + 1], outline[v + 2]);
            expect(norm).toBeCloseTo(1.5 * 1.01, 5);
        }
    });
});

describe('colorCoordinate', () => {
    it('scales 0..vmax into 0..1 and clamps (contrast 1 = plain linear)', () => {
        expect(colorCoordinate(0, 2)).toBe(0);
        expect(colorCoordinate(1, 2)).toBe(0.5);
        expect(colorCoordinate(3, 2)).toBe(1);
        expect(colorCoordinate(1, 0)).toBe(0);
    });

    it('pivots the contrast gain about the isotropic level (enhancement 1)', () => {
        const vmax = 4;
        const pivot = 1 / vmax; // the isotropic enhancement's color coordinate
        // The isotropic level is the fixed point of the gain for any contrast.
        expect(colorCoordinate(1, vmax, 0.5)).toBeCloseTo(pivot, 12);
        expect(colorCoordinate(1, vmax, 3)).toBeCloseTo(pivot, 12);
        // >1 pushes a super-isotropic cell higher, a sub-isotropic cell lower.
        const superBase = 2 / vmax;
        expect(colorCoordinate(2, vmax, 2)).toBeCloseTo(pivot + 2 * (superBase - pivot), 12);
        expect(colorCoordinate(2, vmax, 2)).toBeGreaterThan(colorCoordinate(2, vmax, 1));
        expect(colorCoordinate(0.5, vmax, 2)).toBeLessThan(colorCoordinate(0.5, vmax, 1));
        // Clamped to [0, 1] even at high contrast.
        expect(colorCoordinate(vmax, vmax, 3)).toBe(1);
        expect(colorCoordinate(0, vmax, 3)).toBe(0);
    });
});

describe('colorbarGradient', () => {
    it('builds a CSS gradient with the requested stops', () => {
        const gradient = colorbarGradient('viridis', 4);
        expect(gradient.startsWith('linear-gradient(to right, ')).toBe(true);
        expect(gradient.match(/rgb\(/g)).toHaveLength(5);
        expect(gradient).toContain('100.0%');
    });

    it('tracks the cell color transfer under contrast (bar and sphere agree)', () => {
        const vmax = 5;
        const contrast = 2.4;
        // Each stop's color must equal what a cell of that enhancement gets.
        const gradient = colorbarGradient('viridis', 10, contrast, vmax);
        const colors = [...gradient.matchAll(/rgb\((\d+), (\d+), (\d+)\) ([\d.]+)%/g)];
        expect(colors).toHaveLength(11);
        colors.forEach((match) => {
            const percent = Number(match[4]);
            const enhancement = (percent / 100) * vmax;
            const expected = sampleColormap('viridis', colorCoordinate(enhancement, vmax, contrast));
            expect(Number(match[1])).toBe(expected[0]);
            expect(Number(match[2])).toBe(expected[1]);
            expect(Number(match[3])).toBe(expected[2]);
        });
    });
});

describe('formatDirection', () => {
    it('formats a compact bracketed triple', () => {
        expect(formatDirection([0.1234, -0.9876, 0.05])).toBe('[0.12, -0.99, 0.05]');
    });
});
