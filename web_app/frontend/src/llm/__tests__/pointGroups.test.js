// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// What a point group implies (context/pointGroups.js), and what the
// group–subgroup check reads from a symbol (the local-symmetry interview,
// 2026-10-24): each point group's order, the subgroups of each (the 32 groups'
// subgroup diagram), and the point group of a space-group symbol however it
// is written — checked against the app's own table of the 230 groups.

import { describe, expect, it } from 'vitest';
import { isPointSubgroup, pointGroupFacts, pointGroupOfSymbol, pointGroupOrder } from '../context/pointGroups';
import { SPACE_GROUPS } from '../../spaceGroupTable';

const ALL = ['1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm', '4', '-4', '4/m', '422', '4mm', '-42m', '4/mmm', '3', '-3', '32', '3m', '-3m',
    '6', '-6', '6/m', '622', '6mm', '-6m2', '6/mmm', '23', 'm-3', '432', '-43m', 'm-3m'];

describe('point-group orders and subgroups', () => {
    it('knows the order of each of the 32 point groups', () => {
        expect(ALL.every((group) => pointGroupFacts(group))).toBe(true);
        expect(Object.fromEntries(ALL.map((group) => [group, pointGroupOrder(group)]))).toMatchObject({
            1: 1, '-1': 2, mmm: 8, '4/mmm': 16, '-3m': 12, '6/mmm': 24, '-43m': 24, 'm-3m': 48
        });
        expect(pointGroupOrder('not a group')).toBeNull();
    });

    it('reads a subgroup in any orientation, the group itself included, and every subgroup\'s order divides its parent\'s', () => {
        expect(isPointSubgroup('-42m', '-43m')).toBe(true);
        expect(isPointSubgroup('3m', '-43m')).toBe(true);
        expect(isPointSubgroup('mm2', '-6m2')).toBe(true);
        expect(isPointSubgroup('-3m', 'm-3m')).toBe(true);
        expect(isPointSubgroup('4/mmm', 'm-3m')).toBe(true);
        expect(isPointSubgroup('-43m', '-43m')).toBe(true);
        // -43m has no inversion centre and no four-fold rotation: 4/mmm, 4mm, mmm and m-3 are none of its subgroups.
        for (const group of ['4/mmm', '4mm', '4', 'mmm', 'm-3', '-1']) expect(isPointSubgroup(group, '-43m'), group).toBe(false);
        expect(isPointSubgroup('m-3m', '4mm')).toBe(false);
        expect(isPointSubgroup('6', 'm-3m')).toBe(false);
        expect(isPointSubgroup('x', 'm-3m')).toBeNull();
        for (const parent of ALL) {
            for (const sub of ALL) {
                if (isPointSubgroup(sub, parent)) expect(pointGroupOrder(parent) % pointGroupOrder(sub), `${sub} in ${parent}`).toBe(0);
            }
        }
    });
});

describe('the point group of a space-group symbol', () => {
    it('reads every one of the 230 groups\' symbols; with their screws unmarked, rightly or not at all', () => {
        const ambiguous = [];
        for (const [number, symbol, pointGroup, legacy] of SPACE_GROUPS) {
            expect({ number, symbol, read: pointGroupOfSymbol(symbol) }).toEqual({ number, symbol, read: pointGroup });
            if (legacy) expect(pointGroupOfSymbol(legacy), legacy).toBe(pointGroup);
            // Without the screws' marks two readings may disagree (P4322: 4₃22 or 4 3₂ 2): then no group, never a wrong one.
            const unmarked = pointGroupOfSymbol(symbol.replace(/_/g, ''));
            if (unmarked === null) ambiguous.push(symbol);
            else expect({ number, symbol, unmarked }).toEqual({ number, symbol, unmarked: pointGroup });
        }
        // P4322 (4₃22 or 4 3₂ 2) and P32 (3₂, or 3 and 2).
        expect(ambiguous).toEqual(['P4_322', 'P3_2']);
    });

    it('reads a symbol however it is written, and nothing that is not one', () => {
        expect(pointGroupOfSymbol('P-42₁m')).toBe('-42m');
        expect(pointGroupOfSymbol('P−4 2_1 m')).toBe('-42m');
        expect(pointGroupOfSymbol('$P\\bar{4}2_1m$')).toBe('-42m');
        expect(pointGroupOfSymbol('P4₂/mnm')).toBe('4/mmm');
        expect(pointGroupOfSymbol('F-43m (No. 216)')).toBe('-43m');
        expect(pointGroupOfSymbol('P 1 2/m 1')).toBe('2/m');
        expect(pointGroupOfSymbol('P3m1')).toBe('3m');
        expect(pointGroupOfSymbol('P-62m')).toBe('-6m2');
        for (const text of ['F-43', 'the', 'X4/mmm', 'P7', '', null]) expect(pointGroupOfSymbol(text), String(text)).toBeNull();
    });
});
