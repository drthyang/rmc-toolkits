// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Numbers by reference (copilot/references.js), on their own: the store the
// verification fills an answer from, and turnState.js the earlier text
// segments. In a turn: verification.test.js and agentLoop.test.js.

import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
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

// The form a model writes from the root of the tool message it read,
// {"ref": "r5", "checks": […], "data": {…}}: {{r5.data.<path>}}. A model in
// the dock wrote this exact sentence on 2026-11-01 (polyhedral_distortion on
// the Demo's TaSe6 octahedra as r5), and gemma4:26b most of its references in
// the validation's first live runs; both were left unfilled.
describe('a reference written from the tool message\'s root', () => {
    const HERE = dirname(fileURLToPath(import.meta.url));
    // The Demo's polyhedral_distortion result exactly as the model reads it
    // (the model-visible snapshot of its tool message).
    const shown = JSON.parse(readFileSync(join(HERE, '../../__tests__/modelVisible/__snapshots__/tool-results/polyhedral_distortion.json'), 'utf8'));
    const USER_TEXT = 'has σ²θ = {{r5.data.average_structure.angle_variance_deg2.mean}} deg² and a trans angle of '
        + '{{r5.data.average_structure.trans_angle_deg.mean}}°,';

    const turnStore = () => {
        const store = createReferenceStore();
        ['r1', 'r2', 'r3', 'r4'].forEach((ref) => store.put(ref, { other: 1 }));
        store.put('r5', shown.content.data);
        return store;
    };

    it('fills the user\'s sentence with the engine\'s values', () => {
        expect(shown.tool).toBe('polyhedral_distortion');
        const out = fillReferences(USER_TEXT, turnStore());
        expect(out.text).toBe('has σ²θ = 111 deg² and a trans angle of 159.6°,');
        expect(out).toMatchObject({ written: 2, filled: 2, unknown: [] });
    });

    it('names the same value as the data path, in prose and in JSON', () => {
        const store = turnStore();
        expect(store.resolve('r5.data.configuration.trans_angle_deg.mean')).toBe(store.resolve('r5.configuration.trans_angle_deg.mean'));
        expect(store.resolve('r5.data.average_structure.bond_ladder_A')).toBe('[2.523, 2.525, 2.527, 2.74, 2.744, 2.747]');
        expect(store.resolve('r5.data.coordination.counts[6]')).toBe('15932');
        expect(store.resolve('r5.data')).toBeNull();
    });

    it('reads a result\'s own data key first, if a result ever had one', () => {
        const store = createReferenceStore();
        store.put('r1', { data: { x: 1 }, x: 2 });
        expect(store.resolve('r1.data.x')).toBe('1');
        expect(store.resolve('r1.x')).toBe('2');
    });
});

describe('explain: what a reference that names no value reaches', () => {
    const store = createReferenceStore({ context: { symmetry: { space_group: 'F-43m' } } });
    store.put('r1', { configuration: { distortion_index: { mean: 0.0469, sd: 0.017 } }, sites: [{ ref: 1 }, { ref: 2 }] });

    it('says why, and what is there, never what was meant', () => {
        expect(store.explain('r1.configuration.distortion_index.mean')).toBeNull();
        expect(store.explain('r1.data.configuration.distortion_index.mean')).toBeNull();
        expect(store.explain('r4.peak')).toEqual({ reason: 'no_result', at: '', keys: ['r1'] });
        expect(store.explain('r1.data.configuration.distortion_idx.mean'))
            .toEqual({ reason: 'no_key', at: 'r1.data.configuration', keys: ['distortion_index'] });
        expect(store.explain('r1.configuration.distortion_index'))
            .toEqual({ reason: 'not_a_value', at: 'r1.configuration.distortion_index', keys: ['mean', 'sd'] });
        expect(store.explain('r1.sites')).toEqual({ reason: 'not_a_value', at: 'r1.sites', keys: ['0…1'] });
        expect(store.explain('r1.sites.2.ref')).toEqual({ reason: 'no_key', at: 'r1.sites', keys: ['0…1'] });
        expect(store.explain('ctx.symmetry.point_group')).toEqual({ reason: 'no_key', at: 'ctx.symmetry', keys: ['space_group'] });
        expect(store.explain('r1.sites["<100>"]').reason).toBe('not_a_reference');
    });
});
