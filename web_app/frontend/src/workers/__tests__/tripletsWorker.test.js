// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The browser worker's 'triplets' boundary (pcaKdeWorker.handlePcaMessage):
// the same caps and validation as the Flask /api/triplets route.

import { describe, expect, it } from 'vitest';

import { handlePcaMessage } from '../pcaKdeWorker.js';
import { APP_MAX_ANGLES, tripletRequestFromInputs } from '../triplets.js';

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

describe('triplets worker boundary: missing window bounds', () => {
    // triplets.physics.4/22, numerics.14, parity.27: Number('') and
    // Number(null) are 0, so a cleared minimum silently became rmin = 0.
    const text = cubicRmc6f('Se', 6, 2.5);
    const base = { end1: 'Se', apex: 'Se', end2: 'Se', r12Min: 2, r12Max: 2.6 };

    for (const blank of ['', '   ', null]) {
        it(`rejects r12Min = ${JSON.stringify(blank)} instead of using 0`, async () => {
            await expect(request(text, { ...base, r12Min: blank }))
                .rejects.toThrow(/r12Min\/r12Max are required together; missing r12Min/);
        });
    }

    it('rejects a blank r12Max', async () => {
        await expect(request(text, { ...base, r12Max: '' })).rejects.toThrow(/missing r12Max/);
    });

    it('treats a blank B–C bound like the Flask route: required together', async () => {
        await expect(request(text, { ...base, r23Min: '', r23Max: 2.6 }))
            .rejects.toThrow(/r23Min\/r23Max are required together/);
        await expect(request(text, { ...base, r23Min: '  ', r23Max: 2.6 }))
            .rejects.toThrow(/r23Min\/r23Max are required together/);
    });

    it('rejects an empty binWidth rather than binning at 0', async () => {
        await expect(request(text, { ...base, binWidth: '' })).rejects.toThrow(/binWidth/);
    });

    it('the engine itself still rejects a blank bound if handed one', async () => {
        const { bondAngleSummary } = await import('../triplets.js');
        expect(() => bondAngleSummary([[0.5, 0.5, 0.5]], ['Se'], [[5, 0, 0], [0, 5, 0], [0, 0, 5]], {
            triplet: ['Se', 'Se', 'Se'], bond12: ['  ', 2]
        })).toThrow(/bond12 bounds must be finite/);
    });

    it('still accepts numeric strings (an HTTP-style request)', async () => {
        const result = await request(text, { ...base, r12Min: '2', r12Max: '2.6', binWidth: '1' });
        expect(result.bond12).toEqual([2, 2.6]);
        expect(result.angleCount).toBe(216 * 15);
    });
});

describe('page inputs → triplets request', () => {
    const inputs = {
        end1: 'Se', apex: 'Nb', end2: 'Se',
        r12Min: '2.2', r12Max: '2.9', split23: false, r23Min: '', r23Max: '', binWidth: '1.0'
    };

    it('builds numeric params, B–C only when split', () => {
        expect(tripletRequestFromInputs(inputs)).toEqual({
            end1: 'Se', apex: 'Nb', end2: 'Se', r12Min: 2.2, r12Max: 2.9, binWidth: 1
        });
        expect(tripletRequestFromInputs({ ...inputs, split23: true, r23Min: '2.5', r23Max: '3.1' }))
            .toMatchObject({ r23Min: 2.5, r23Max: 3.1 });
    });

    it('a cleared field is an error naming it, never 0', () => {
        expect(() => tripletRequestFromInputs({ ...inputs, r12Min: '' }))
            .toThrow(/A–B window minimum is empty/);
        expect(() => tripletRequestFromInputs({ ...inputs, r12Max: ' ' }))
            .toThrow(/A–B window maximum is empty/);
        expect(() => tripletRequestFromInputs({ ...inputs, split23: true, r23Min: '', r23Max: '3' }))
            .toThrow(/B–C window minimum is empty/);
        expect(() => tripletRequestFromInputs({ ...inputs, binWidth: '' }))
            .toThrow(/Bin width is empty/);
        expect(() => tripletRequestFromInputs({ ...inputs, r12Min: 'abc' }))
            .toThrow(/A–B window minimum is not a number/);
    });

    it('ignores the B–C fields while the split is off', () => {
        expect(tripletRequestFromInputs({ ...inputs, r23Min: '', r23Max: 'x' }).r23Min).toBeUndefined();
    });
});
