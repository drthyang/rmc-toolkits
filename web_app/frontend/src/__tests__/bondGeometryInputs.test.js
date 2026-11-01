// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Bond Geometry page's inputs as the Copilot's change tool reads and sets
// them (bondGeometryInputs.js): a change sets the controls it names and keeps
// the others, is refused in plain words when it contradicts itself or the
// page could not compute what it leaves, and reads as one line plus a line
// per control it changes.

import { describe, expect, it } from 'vitest';
import {
    BOND_GEOMETRY_LIMITS, DEFAULT_INPUTS, changedInputs, describeChange, inputsLine, inputsProblem, normalizeInputs,
    requestOf, sameInputs, sameRequest, toolInputs
} from '../bondGeometryInputs';

const ELEMENTS = ['Ga', 'Se', 'Ta'];
const PAGE = Object.freeze({ ...DEFAULT_INPUTS, end1: 'Se', apex: 'Ta', end2: 'Se' });
const change = (args, current = PAGE) => changedInputs(current, args, { elements: ELEMENTS });

describe('the page\'s inputs', () => {
    it('the boxes\' ranges are the app\'s caps: windows 0–15 Å, bins 0.1–45°', () => {
        expect(BOND_GEOMETRY_LIMITS).toMatchObject({ windowMinA: 0, windowMaxA: 15, binMinDeg: 0.1, binMaxDeg: 45 });
        expect(DEFAULT_INPUTS).toEqual({
            end1: '', apex: '', end2: '', r12Min: '2.00', r12Max: '3.00', split23: false, r23Min: '2.00', r23Max: '3.00', binWidth: '1.0'
        });
    });

    it('normalise to strings and one switch; compare box for box; make the page\'s request', () => {
        expect(normalizeInputs({ end1: 'Se', r12Min: 2.2, split23: 1 })).toMatchObject({ end1: 'Se', r12Min: '2.2', split23: true, binWidth: '1.0' });
        expect(sameInputs(PAGE, { ...PAGE })).toBe(true);
        expect(sameInputs(PAGE, { ...PAGE, r12Max: '3.0' })).toBe(false);
        expect(requestOf(PAGE)).toEqual({ end1: 'Se', apex: 'Ta', end2: 'Se', r12Min: 2, r12Max: 3, binWidth: 1 });
        expect(requestOf({ ...PAGE, r12Max: '' })).toBeNull();
        expect(sameRequest(requestOf(PAGE), requestOf({ ...PAGE, r12Max: '3.0' }))).toBe(true);
        expect(sameRequest(null, requestOf(PAGE))).toBe(false);
    });
});

describe('a change on the page (changedInputs)', () => {
    it('sets the controls it names and keeps the rest as the page has them', () => {
        expect(change({ r12_min: 2.21, r12_max: 2.65 })).toEqual({ ...PAGE, r12Min: '2.21', r12Max: '2.65' });
        expect(change({ end1: 'Ga', bin_width: 2 })).toEqual({ ...PAGE, end1: 'Ga', binWidth: '2' });
    });

    it('a B–C window turns Distinct B–C on; distinct_bc false turns it off and keeps the B–C boxes', () => {
        const split = change({ r23_min: 2.4, r23_max: 2.8 });
        expect(split).toMatchObject({ split23: true, r23Min: '2.4', r23Max: '2.8' });
        expect(change({ distinct_bc: false }, split)).toEqual({ ...split, split23: false });
        // Already split: distinct_bc true alone keeps its B–C window.
        expect(change({ distinct_bc: true }, split)).toEqual(split);
    });

    it('refuses arguments that contradict each other, in plain words', () => {
        expect(() => change({ r23_min: 2.4 })).toThrow('give r23_min and r23_max together, or neither');
        expect(() => change({ r23_min: 2.4, r23_max: 2.8, distinct_bc: false })).toThrow(/leave distinct_bc out, or set it true/);
        expect(() => change({ distinct_bc: true })).toThrow('give r23_min and r23_max to turn the distinct B–C window on');
    });

    it('refuses a page it could not compute, naming the box and the argument that sets it', () => {
        expect(() => change({ r12_max: 1.5 })).toThrow('the page cannot compute that: the A–B window would be 2.00–1.50 Å: its minimum must be below its maximum');
        expect(() => change({ r12_min: 3, r12_max: 3 })).toThrow(/its minimum must be below its maximum/);
        expect(() => change({ end1: 'Nb' })).toThrow('the page cannot compute that: A is Nb, which is not an element of this run (Ga, Se, Ta): give end1');
        expect(() => change({ bin_width: 0.05 })).toThrow(/0.05° bin width leaves the page's range, 0.1–45°/);
        expect(() => change({ r12_max: 16 })).toThrow(/leaves the page's range, 0–15 Å/);
        // What the user left in a box the change does not set.
        expect(() => change({ end1: 'Ga' }, { ...PAGE, r12Max: '' })).toThrow('the page cannot compute that: the A–B window maximum box is empty: give r12_max');
        expect(() => change({ end1: 'Ga' }, { ...PAGE, binWidth: 'x' })).toThrow('the bin width box holds "x", not a number: give bin_width');
        // A hidden B–C box is not read while the switch is off.
        expect(change({ end1: 'Ga' }, { ...PAGE, r23Max: '' })).toMatchObject({ end1: 'Ga', r23Max: '' });
        expect(inputsProblem({ ...PAGE, end1: '' })).toBe('A names no element yet: give end1');
        expect(inputsProblem(PAGE, { elements: ELEMENTS })).toBeNull();
    });
});

describe('how a change reads (describeChange, inputsLine, toolInputs)', () => {
    it('one line, the page after the change first; a detail line per control it changes', () => {
        const after = change({ end1: 'Ga', r12_min: 2.21, r12_max: 2.65, r23_min: 2.4, r23_max: 2.8, bin_width: 2 });
        expect(describeChange(PAGE, after)).toEqual({
            line: 'Ga–Ta–Se · A–B 2.21–2.65 Å, B–C 2.40–2.80 Å · 2° bins',
            details: [
                'Triplet: Se–Ta–Se → Ga–Ta–Se',
                'A–B window: 2.00–3.00 → 2.21–2.65 Å',
                'Distinct B–C window: off → on',
                'B–C window: 2.00–3.00 → 2.40–2.80 Å',
                'Bin width: 1° → 2°'
            ]
        });
        expect(describeChange(PAGE, PAGE)).toEqual({ line: 'Se–Ta–Se · 2.00–3.00 Å · 1° bins', details: ['No input changes: the page computes these inputs.'] });
        expect(inputsLine({ ...PAGE, r12Min: '2.205' })).toBe('Se–Ta–Se · 2.205–3.00 Å · 1° bins');
    });

    it('reads the boxes back as the tool\'s arguments, numbers as numbers', () => {
        expect(toolInputs(PAGE)).toEqual({ end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2, r12_max: 3, distinct_bc: false, bin_width: 1 });
        expect(toolInputs({ ...PAGE, split23: true, r23Min: '2.4', r23Max: '' })).toMatchObject({ distinct_bc: true, r23_min: 2.4, r23_max: null });
    });
});
