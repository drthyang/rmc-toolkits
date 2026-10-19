// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { checkShape, declaredKeys, describeShape, oneOf } from '../copilot/dataShape';

const SHAPE = {
    note: 'string',
    'element?': 'string',
    sites: [{ ref: 'number', shift_A: 'number|null', direction_uvw: oneOf(['number', 'number', 'number'], 'null'), 'zero_spread?': 'boolean' }],
    cn: { '*_per_*': 'number' },
    window_A: ['number', 'number']
};

const good = () => ({
    note: 'n',
    sites: [{ ref: 1, shift_A: 0.07, direction_uvw: [0, 0, 1] }, { ref: 2, shift_A: null, direction_uvw: null, zero_spread: true }],
    cn: { Se_per_Ta: 6, Ta_per_Se: 1.5 },
    window_A: [2.2, 3.2]
});

describe('checkShape', () => {
    it('accepts a value that fits, optional keys absent or present', () => {
        expect(checkShape(good(), SHAPE)).toEqual([]);
        expect(checkShape({ ...good(), element: 'Ta' }, SHAPE)).toEqual([]);
        expect(checkShape({ ...good(), cn: {} }, SHAPE)).toEqual([]);
    });

    it('names every mismatch by the path a reference writes', () => {
        const value = good();
        value.sites[1].shift_A = '0.2';
        value.sites[0].direction_uvw = [0, 1];
        delete value.note;
        value.extra = 1;
        value.cn.coordination = 4;
        value.window_A = [2.2, 'x'];
        expect(checkShape(value, SHAPE)).toEqual([
            'note: missing (a required key)',
            'sites.0.direction_uvw: expected oneOf([number, number, number], null), found an array of 2',
            'sites.1.shift_A: expected number|null, found "0.2" (string)',
            'cn.coordination: not declared in the shape',
            'window_A.1: expected number, found "x" (string)',
            'extra: not declared in the shape'
        ]);
    });

    it('reads values as JSON carries them: undefined keys are absent, non-finite numbers null', () => {
        expect(checkShape({ ...good(), element: undefined }, SHAPE)).toEqual([]);
        expect(checkShape({ a: NaN }, { a: 'number|null' })).toEqual([]);
        expect(checkShape({ a: Infinity }, { a: 'number' })).toEqual(['a: expected number, found Infinity (null in JSON)']);
        expect(checkShape({ a: undefined }, { a: 'number' })).toEqual(['a: missing (a required key)']);
        expect(checkShape([1, undefined], ['number|null'])).toEqual([]);
    });

    it('reports a tuple of the wrong length and a wrong container', () => {
        expect(checkShape({ ...good(), window_A: [1, 2, 3] }, SHAPE)).toEqual([
            'window_A: expected a tuple of 2 [number, number], found an array of 3'
        ]);
        expect(checkShape({ ...good(), sites: {} }, SHAPE)).toEqual(['sites: expected [{ ref, shift_A, direction_uvw, zero_spread? }], found an object']);
        expect(checkShape('text', SHAPE)).toEqual(['(data): expected { note, element?, sites, cn, window_A }, found "text" (string)']);
    });

    it('a oneOf reports the closest of its forms', () => {
        const shape = { by_range: [oneOf({ r_A: 'number', rw: 'number' }, { Q_A: 'number', rw: 'number' })] };
        expect(checkShape({ by_range: [{ r_A: 1, rw: 0.1 }, { Q_A: 2, rw: 0.2 }] }, shape)).toEqual([]);
        expect(checkShape({ by_range: [{ Q_A: 2, rw: 'x' }] }, shape)).toEqual([
            'by_range.0.rw: expected number, found "x" (string) (in the closest of 2 forms)'
        ]);
    });

    it('collects the declared keys it checked, as declaredKeys names them', () => {
        const seen = new Set();
        checkShape(good(), SHAPE, { seen });
        expect([...seen].sort()).toEqual([
            'cn', 'cn.*_per_*', 'note', 'sites', 'sites[].direction_uvw', 'sites[].ref', 'sites[].shift_A', 'sites[].zero_spread', 'window_A'
        ]);
        expect(declaredKeys(SHAPE).sort()).toEqual([...seen, 'element'].sort());
    });

    it('refuses notation it does not know', () => {
        expect(() => checkShape(1, 'nubmer')).toThrow(/unknown type "nubmer"/);
        expect(() => checkShape([], [])).toThrow(/declares no item shape/);
        expect(() => checkShape({}, { 'a?': 'number', a: 'string' })).toThrow(/declared twice/);
        expect(() => checkShape({}, { 'x_*?': 'number' })).toThrow(/optional already/);
        expect(() => declaredKeys({ a: 3 })).toThrow(/invalid data shape/);
        expect(() => oneOf('number')).toThrow(/two or more/);
    });

    it('describes a shape on one line', () => {
        expect(describeShape(SHAPE.sites)).toBe('[{ ref, shift_A, direction_uvw, zero_spread? }]');
        expect(describeShape(oneOf('string', ['number', 'number']))).toBe('oneOf(string, [number, number])');
    });
});
