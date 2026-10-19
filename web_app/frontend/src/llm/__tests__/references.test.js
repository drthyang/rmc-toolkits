// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Numbers by reference (copilot/references.js), on their own: the store the
// verification fills an answer from, and turnState.js the earlier text
// segments. In a turn: verification.test.js and agentLoop.test.js.

import { describe, expect, it } from 'vitest';
import { createReferenceStore, fillReferences, formatReferenceValue } from '../copilot/references';

describe('numbers by reference', () => {
    const store = createReferenceStore({ context: { symmetry: { noise_floor_A: 0.026 } } });
    store.put('r1', { peak: { significance_sigma: 17.2, direction_uvw: [0, 0, -1] }, sites: [{ shift_A: 0.065, element: 'Ta' }] });
    store.put('r2', { local: { long_axis_within_15deg: { '<100>': { z: 12 } } } });

    it('fills each reference from the result it names, in prose and in JSON', () => {
        const out = fillReferences('Peak {{r1.peak.significance_sigma}}σ along {{ r1.peak.direction_uvw }}; {{r1.sites[0].element}} shifted '
            + '{{r1.sites.0.shift_A}} Å, floor {{ctx.symmetry.noise_floor_A}} Å; <100> z = {{r2.local.long_axis_within_15deg.<100>.z}}. '
            + 'VERDICT {"magnitude_A": {{r1.sites.0.shift_A}}, "direction_uvw": {{r1.peak.direction_uvw}}}', store);
        expect(out.text).toBe('Peak 17.2σ along [0, 0, -1]; Ta shifted 0.065 Å, floor 0.026 Å; <100> z = 12. '
            + 'VERDICT {"magnitude_A": 0.065, "direction_uvw": [0, 0, -1]}');
        expect(out).toMatchObject({ written: 8, filled: 8, unknown: [] });
        expect(JSON.parse(out.text.slice(out.text.indexOf('{')))).toEqual({ magnitude_A: 0.065, direction_uvw: [0, 0, -1] });
    });

    it('leaves a reference that names nothing as written, and anything else in braces alone', () => {
        const out = fillReferences('{{r9.x}} and {{r1.peak}} and {{r1.nope}}; $x^{{2}}$ stays.', store);
        expect(out.text).toBe('{{r9.x}} and {{r1.peak}} and {{r1.nope}}; $x^{{2}}$ stays.');
        expect(out).toMatchObject({ written: 3, filled: 0, unknown: ['r9.x', 'r1.peak', 'r1.nope'] });
        expect(fillReferences('No braces.', store)).toEqual({ text: 'No braces.', written: 0, filled: 0, unknown: [] });
        expect(formatReferenceValue(Number.NaN)).toBeNull();
        expect(formatReferenceValue({ a: 1 })).toBeNull();
    });

    it('keeps results under the refs the caller names; without a context, {{ctx…}} names nothing', () => {
        expect(store.size()).toBe(2);
        expect(store).not.toHaveProperty('add');
        expect(fillReferences('{{ctx.run}}', createReferenceStore()).unknown).toEqual(['ctx.run']);
    });
});
