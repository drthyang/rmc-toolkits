// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// dispA (the AI context's mean_disp_A / max_disp_A) must be the site's rms
// displacement in Å through the full cell metric. It used to combine per-axis
// fractional spreads with edge LENGTHS, as if the axes were orthogonal: +10%
// for an isotropic cloud in a hexagonal cell, +23% for fcc in its rhombohedral
// primitive cell, and a different value for the same crystal in another setting.

import { describe, expect, it } from 'vitest';
import { structureFromRmc6f } from '../browserData';

const N = 4;
const SITE = [0.3, 0.6, 0.2];

// Deterministic Cartesian displacements (Å) around the site: ±δ on x, y, z plus
// two diagonal pairs, cycled over the N³ copies — mean exactly zero.
const DISPLACEMENTS = [
    [0.1, 0, 0], [-0.1, 0, 0], [0, 0.1, 0], [0, -0.1, 0], [0, 0, 0.1], [0, 0, -0.1],
    [0.07, 0.07, 0], [-0.07, -0.07, 0]
];

const invert3 = (m) => {
    const [[a, b, c], [d, e, f], [g, h, i]] = m;
    const A = e * i - f * h; const B = -(d * i - f * g); const C = d * h - e * g;
    const det = a * A + b * B + c * C;
    return [
        [A / det, -(b * i - c * h) / det, (b * f - c * e) / det],
        [B / det, (a * i - c * g) / det, -(a * f - c * d) / det],
        [C / det, -(a * h - b * g) / det, (a * e - b * d) / det]
    ];
};

// Unit-cell rows a, b, c (Å); box = N × each row.
const configuration = (unitCell) => {
    const inverse = invert3(unitCell);   // frac (row) = cart (row) · A⁻¹
    const lines = [];
    const used = [];
    let id = 0;
    for (let cx = 0; cx < N; cx += 1) for (let cy = 0; cy < N; cy += 1) for (let cz = 0; cz < N; cz += 1) {
        const d = DISPLACEMENTS[id % DISPLACEMENTS.length];
        used.push(d);
        const dFrac = [0, 1, 2].map((j) => d[0] * inverse[0][j] + d[1] * inverse[1][j] + d[2] * inverse[2][j]);
        const box = [0, 1, 2].map((j) => (SITE[j] + dFrac[j] + [cx, cy, cz][j]) / N);
        id += 1;
        lines.push(`${id} Ga [1] ${box.map((v) => v.toFixed(12)).join(' ')} 1 ${cx} ${cy} ${cz}`);
    }
    const text = [
        `Number of atoms: ${id}`,
        `Supercell dimensions: ${N} ${N} ${N}`,
        'Lattice vectors (Ang):',
        ...unitCell.map((row) => row.map((v) => (v * N).toFixed(10)).join(' ')),
        'Atoms:',
        ...lines
    ].join('\n');
    // Cartesian rms about the mean, computed directly from the displacements used.
    const mean = [0, 1, 2].map((j) => used.reduce((s, d) => s + d[j], 0) / used.length);
    const rms = Math.sqrt(used.reduce((s, d) => s + d.reduce((t, v, j) => t + (v - mean[j]) ** 2, 0), 0) / used.length);
    return { text, rms };
};

const HEXAGONAL = [[5, 0, 0], [-2.5, 2.5 * Math.sqrt(3), 0], [0, 0, 8]];
const RHOMBOHEDRAL = [[0, 2.5, 2.5], [2.5, 0, 2.5], [2.5, 2.5, 0]];   // fcc primitive, 60°
const CUBIC = [[5, 0, 0], [0, 5, 0], [0, 0, 5]];

describe('dispA through the cell metric', () => {
    it.each([['hexagonal', HEXAGONAL], ['rhombohedral (60°)', RHOMBOHEDRAL], ['cubic', CUBIC]])(
        '%s cell: equals the Cartesian rms displacement',
        (_, unitCell) => {
            const { text, rms } = configuration(unitCell);
            const [site] = structureFromRmc6f({ path: 'cell.rmc6f', text }).basis;
            expect(site.dispA).toBeCloseTo(rms, 9);
        }
    );
});
