// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The browser worker's 'triplets' boundary (pcaKdeWorker.handlePcaMessage):
// the same caps and validation as the Flask /api/triplets route.

import { describe, expect, it } from 'vitest';

import { handlePcaMessage } from '../pcaKdeWorker.js';
import { APP_MAX_ANGLES } from '../triplets.js';

// Simple-cubic .rmc6f of one element: n^3 atoms at spacing `spacing` Å.
const cubicRmc6f = (element, n, spacing) => {
    const edge = n * spacing;
    const lines = [
        `Supercell dimensions: ${n} ${n} ${n}`,
        'Lattice vectors (Ang):',
        `${edge} 0 0`,
        `0 ${edge} 0`,
        `0 0 ${edge}`,
        'Atoms:'
    ];
    let atom = 0;
    for (let ix = 0; ix < n; ix += 1) {
        for (let iy = 0; iy < n; iy += 1) {
            for (let iz = 0; iz < n; iz += 1) {
                atom += 1;
                // A tiny deterministic jitter keeps shells off exact distances.
                const jitter = ((atom * 7919) % 101) / 1e5;
                lines.push(`${atom} ${element} [1] ${((ix + jitter) / n).toFixed(10)} `
                    + `${((iy + jitter) / n).toFixed(10)} ${((iz + jitter) / n).toFixed(10)} 1 ${ix} ${iy} ${iz}`);
            }
        }
    }
    return lines.join('\n');
};

const request = (text, params) =>
    handlePcaMessage({ kind: 'triplets', ...params }, async () => text);

describe('triplets worker boundary: work budget', () => {
    it('refuses a spec over the shared angle budget before forming any angle', async () => {
        // 1000 atoms, ~400 bonds each inside 9.5 Å: ~8e7 angles, well under
        // the 15 Å rmax cap but over the budget (triplets.physics.1 et al.).
        const text = cubicRmc6f('Se', 10, 2.08);
        const started = Date.now();
        await expect(request(text, {
            end1: 'Se', apex: 'Se', end2: 'Se', r12Min: 0.5, r12Max: 9.5
        })).rejects.toThrow(/angles, over the limit of 50,000,000/);
        // Refused from the bond lists alone: no 8e7-angle pass happened.
        expect(Date.now() - started).toBeLessThan(4000);
        expect(APP_MAX_ANGLES).toBe(50_000_000);
    });

    it('computes a spec within the budget', async () => {
        const text = cubicRmc6f('Se', 6, 2.5);
        const result = await request(text, {
            end1: 'Se', apex: 'Se', end2: 'Se', r12Min: 2, r12Max: 2.6
        });
        expect(result.apexCount).toBe(216);
        expect(result.angleCount).toBe(216 * 15);
    });
});
