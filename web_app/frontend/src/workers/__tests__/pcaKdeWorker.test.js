// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, it, expect } from 'vitest';

import { handlePcaMessage } from '../pcaKdeWorker.js';

// Minimal .rmc6f text with one reference site per element, n^3 box copies each.
const buildRmc6f = (elements, { supercell = [4, 4, 4], edge = 8, seed = 1 } = {}) => {
    let state = seed >>> 0;
    const rand = () => {
        state = (state * 1103515245 + 12345) & 0x7fffffff;
        return state / 0x7fffffff - 0.5;
    };
    const [sx, sy, sz] = supercell;
    const lines = [
        `Supercell dimensions ${sx} ${sy} ${sz}`,
        'Lattice vectors (Ang):',
        `${edge * sx} 0 0`,
        `0 ${edge * sy} 0`,
        `0 0 ${edge * sz}`,
        'Atoms:'
    ];
    let atom = 0;
    elements.forEach((element, ref) => {
        for (let ix = 0; ix < sx; ix += 1) {
            for (let iy = 0; iy < sy; iy += 1) {
                for (let iz = 0; iz < sz; iz += 1) {
                    atom += 1;
                    const base = [ix / sx, iy / sy, iz / sz];
                    const c = base.map((v) => v + 0.02 * rand());
                    lines.push(`${atom} ${element} [1] ${c[0].toFixed(8)} ${c[1].toFixed(8)} `
                        + `${c[2].toFixed(8)} ${ref + 1} ${ix} ${iy} ${iz}`);
                }
            }
        }
    });
    return lines.join('\n');
};

describe('pcaKdeWorker cache is content-addressed', () => {
    it('re-parses when the .rmc6f text changes (does not return the previous dataset)', async () => {
        const datasetA = buildRmc6f(['Ga', 'Se'], { seed: 1 });
        const datasetB = buildRmc6f(['Nb', 'Ta', 'O'], { seed: 2 });

        const a1 = await handlePcaMessage({ kind: 'sites' }, async () => datasetA);
        expect(a1.elements).toEqual(['Ga', 'Se']);

        // Loading a different dataset must reflect the new model, not the cached one.
        const b = await handlePcaMessage({ kind: 'sites' }, async () => datasetB);
        expect(b.elements).toEqual(['Nb', 'O', 'Ta']);
        expect(b.sites).toHaveLength(3);

        // And switching back returns the first dataset again (not a stale mix).
        const a2 = await handlePcaMessage({ kind: 'sites' }, async () => datasetA);
        expect(a2.elements).toEqual(['Ga', 'Se']);
        expect(a2.sites).toHaveLength(2);
    });

    it('serves a repeated identical request from cache (same result)', async () => {
        const dataset = buildRmc6f(['Ga', 'Se', 'Se'], { seed: 5 });
        const first = await handlePcaMessage({ kind: 'sites' }, async () => dataset);
        const second = await handlePcaMessage({ kind: 'sites' }, async () => dataset);
        expect(second.totalAtoms).toBe(first.totalAtoms);
        expect(second.referenceNumbers).toEqual(first.referenceNumbers);
    });

    it('kde requests also follow the current dataset', async () => {
        const datasetA = buildRmc6f(['Ga'], { seed: 7 });
        const datasetB = buildRmc6f(['Nb'], { seed: 8 });
        const kdeA = await handlePcaMessage({ kind: 'kde', referenceNumber: 1, grid: 12, projections: false }, async () => datasetA);
        expect(kdeA.element).toBe('Ga');
        const kdeB = await handlePcaMessage({ kind: 'kde', referenceNumber: 1, grid: 12, projections: false }, async () => datasetB);
        expect(kdeB.element).toBe('Nb');
    });
});
