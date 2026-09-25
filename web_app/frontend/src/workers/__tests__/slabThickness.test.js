// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// z/dz are fractions of the unit cube's projection range along the normal, not
// of a cell edge (except for the a/b/c presets) and never Angstrom. The page
// prints the real slab thickness, dz * (|h| + |k| + |l|) * d_hkl.

import { describe, expect, it } from 'vitest';
import { slabThicknessAngstrom } from '../slabSelection';

const unit = (h) => {
    const length = Math.hypot(...h);
    return h.map((value) => value / length);
};
const range = (n) => {
    const corners = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]];
    const depths = corners.map((corner) => corner[0] * n[0] + corner[1] * n[1] + corner[2] * n[2]);
    return [Math.min(...depths), Math.max(...depths)];
};

describe('slab thickness in Angstrom', () => {
    const cubic = [[10.395, 0, 0], [0, 10.395, 0], [0, 0, 10.395]];

    it('is dz times the cell edge for a preset normal', () => {
        const n = [0, 0, 1];
        expect(slabThicknessAngstrom(0.08, n, range(n), cubic)).toBeCloseTo(0.8316, 10);
    });

    it('is dz (|h|+|k|+|l|) d_hkl for a custom plane: 1.41x and 1.73x the naive value in a cubic cell', () => {
        const n110 = unit([1, 1, 0]);
        const n111 = unit([1, 1, 1]);
        // d_110 = a / sqrt(2), d_111 = a / sqrt(3)
        expect(slabThicknessAngstrom(0.08, n110, range(n110), cubic)).toBeCloseTo(0.08 * 2 * 10.395 / Math.sqrt(2), 10);
        expect(slabThicknessAngstrom(0.08, n111, range(n111), cubic)).toBeCloseTo(0.08 * 3 * 10.395 / Math.sqrt(3), 10);
    });

    it('uses the interplanar spacing of an oblique cell', () => {
        // Hexagonal a = 5, c = 4: d_100 = a sqrt(3) / 2.
        const hexagonal = [[5, 0, 0], [-2.5, 5 * Math.sqrt(3) / 2, 0], [0, 0, 4]];
        const n = [1, 0, 0];
        expect(slabThicknessAngstrom(0.1, n, range(n), hexagonal)).toBeCloseTo(0.1 * 5 * Math.sqrt(3) / 2, 10);
    });
});
