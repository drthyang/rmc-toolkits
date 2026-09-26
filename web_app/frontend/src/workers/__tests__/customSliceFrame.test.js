// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Atomic Density page sends the in-plane frame of a custom (hkl) slice to
// /api/kde/slice, so the Flask map is drawn in the frame the browser worker,
// the Slab In Cell panel and the panel aspect use. Before, the page sent only
// the normal and Flask picked its own axes: the default (1 1 0) map came back
// rotated 90 degrees. tests/test_kde_custom_frame.py checks the route side
// with the same literal frames and the same parity-golden cases.

import { describe, expect, it } from 'vitest';
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { freePlaneBasis, kdeSliceQuery } from '../slabSelection.js';

const here = dirname(fileURLToPath(import.meta.url));
const unit = (vector) => {
    const length = Math.hypot(...vector);
    return vector.map((value) => value / length);
};
const half = Math.sqrt(0.5);
const close = (actual, expected) => actual.forEach((value, i) => expect(value).toBeCloseTo(expected[i], 14));

describe('custom slice frame', () => {
    it('matches the literal frames of the default planes (as test_kde_custom_frame.py)', () => {
        const f110 = freePlaneBasis(unit([1, 1, 0]));
        close(f110.u, [half, -half, 0]);
        close(f110.v, [0, 0, -1]);
        const f101 = freePlaneBasis(unit([1, 0, 1]));
        close(f101.u, [half, 0, -half]);
        close(f101.v, [0, 1, 0]);
    });

    it('matches the frames pinned in the browser-parity golden', () => {
        const fixture = JSON.parse(readFileSync(join(here, '..', '..', '__tests__', 'fixtures', 'kde_parity_fixture.json'), 'utf8'));
        const custom = fixture.cases.filter((c) => c.name.includes('(1'));
        expect(custom.map((c) => c.name)).toEqual(expect.arrayContaining(['demo Ta (110) bw=0.06', 'demo Se (101) bw=0.05']));
        for (const testCase of custom) {
            const { u, v } = freePlaneBasis(testCase.normal);
            close(u, testCase.u);
            close(v, testCase.v);
        }
    });

    it('sends the frame with a custom slice, and only the name with a preset', () => {
        const normal = unit([1, 1, 0]);
        const { u, v } = freePlaneBasis(normal);
        expect(kdeSliceQuery('custom', { normal, u, v })).toEqual({
            orientation: 'custom',
            nx: normal[0], ny: normal[1], nz: normal[2],
            ux: u[0], uy: u[1], uz: u[2],
            vx: v[0], vy: v[1], vz: v[2]
        });
        expect(kdeSliceQuery('c', { normal: [0, 0, 1], u: [1, 0, 0], v: [0, 1, 0] }))
            .toEqual({ orientation: 'c', nx: 0, ny: 0, nz: 1 });
    });

    it('StructurePage builds its Flask KDE query and its frame from these helpers', () => {
        const source = readFileSync(join(here, '..', '..', 'components', 'StructurePage.jsx'), 'utf8');
        expect(source).toMatch(/\.\.\.kdeSliceQuery\(sliceDirection, sliceConfig\)/);
        expect(source).toMatch(/const makeFreePlaneBasis = freePlaneBasis;/);
        expect(source).not.toMatch(/nx: sliceConfig\.normal\[0\]/);
    });
});
