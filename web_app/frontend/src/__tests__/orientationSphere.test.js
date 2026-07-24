// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, it, expect } from 'vitest';

import {
    buildCellMesh,
    buildCellOutline,
    colorCoordinate,
    colorbarGradient,
    formatDirection
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

describe('colorCoordinate', () => {
    it('scales 0..vmax into 0..1 and clamps', () => {
        expect(colorCoordinate(0, 2)).toBe(0);
        expect(colorCoordinate(1, 2)).toBe(0.5);
        expect(colorCoordinate(3, 2)).toBe(1);
        expect(colorCoordinate(1, 0)).toBe(0);
    });
});

describe('colorbarGradient', () => {
    it('builds a CSS gradient with the requested stops', () => {
        const gradient = colorbarGradient('viridis', 4);
        expect(gradient.startsWith('linear-gradient(to right, ')).toBe(true);
        expect(gradient.match(/rgb\(/g)).toHaveLength(5);
        expect(gradient).toContain('100.0%');
    });
});

describe('formatDirection', () => {
    it('formats a compact bracketed triple', () => {
        expect(formatDirection([0.1234, -0.9876, 0.05])).toBe('[0.12, -0.99, 0.05]');
    });
});
