// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Log scale changes the contour levels, never whether there are contours: the
// worker, like kde.py, judges positivity on the linear density. Python twin:
// tests/test_kde_contours.py (same Kronecker point cloud, same slice).

import { describe, expect, it } from 'vitest';
import { computeKde } from '../localKdeWorker';

const kroneckerPoints = (count) => Array.from({ length: count }, (_, i) => ({
    x: ((i + 1) * 0.7548776662466927) % 1,
    y: ((i + 1) * 0.5698402909980532) % 1,
    z: ((i + 1) * 0.6180339887498949) % 1
}));

// StructurePage's frame for a custom (1 1 1) normal (makeSliceConfig).
const n = 1 / Math.sqrt(3);
const u = [2 / Math.sqrt(6), -1 / Math.sqrt(6), -1 / Math.sqrt(6)];
const v = [n * u[2] - n * u[1], n * u[0] - n * u[2], n * u[1] - n * u[0]];

const slice = (overrides) => computeKde({
    points: kroneckerPoints(8000),
    normal: [n, n, n],
    uVector: u,
    vVector: v,
    range: [0, 3 * n],
    zCenter: 0.5,
    thickness: 0.08,
    bandwidth: 0.1,
    gridSize: 48,
    logScale: false,
    ...overrides
});

describe('contours in log mode', () => {
    it('keeps all eight contours on an oblique disordered slice whose log peak is negative', async () => {
        const logMap = await slice({ logScale: true });
        const linearMap = await slice({ logScale: false });
        expect(logMap.vmax).toBeLessThan(0);
        expect(linearMap.contours).toHaveLength(8);
        expect(logMap.contours).toHaveLength(8);
    });

    it('draws no contours for a declined slab in either mode', async () => {
        for (const logScale of [false, true]) {
            const result = await slice({ logScale, bandwidth: 0, normal: [0, 0, 1], uVector: [1, 0, 0], vVector: [0, 1, 0], range: [0, 1] });
            expect(result.contours).toEqual([]);
        }
    });
});
