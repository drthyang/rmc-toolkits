// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The ITA standard-setting operations the CIF export moves a structure onto
// (itaOperations.js) are spglib's, group for group: the same full operation set, centring
// included, in origin choice 2 and on hexagonal R axes
// (fixture: tests/generate_ita_operations_fixture.py). Wyckoff rows are tied to the same
// operations by wyckoff.test.js, so letters and operations agree on the origin.

import { describe, expect, it } from 'vitest';

import { CENTRING, ORIGIN_CHOICE_2, itaGenerators, itaOperations } from '../itaOperations.js';
import { parseCoordinateForm } from '../wyckoff.js';
import { spaceGroupNumber } from '../spaceGroupTable.js';
import { SPACE_GROUP_FIXTURES } from './fixtures/spaceGroups.js';
import FIXTURE from './fixtures/ita_operations_fixture.json';

const key = ({ R, t }) => `${R.flat().join(',')}|${t.map((v) => Math.round((((v % 1) + 1) % 1) * 48) % 48).join(',')}`;

describe('ITA standard operations', () => {
    it('are spglib\'s for all 230 groups', () => {
        expect(Object.keys(FIXTURE.groups)).toHaveLength(230);
        for (let number = 1; number <= 230; number++) {
            const expected = new Set(FIXTURE.groups[number].operations.map((text) => key(parseCoordinateForm(text))));
            const actual = new Set(itaOperations(number).map(key));
            expect(actual.size, `No. ${number}`).toBe(expected.size);
            for (const k of actual) expect(expected.has(k), `No. ${number}: ${k}`).toBe(true);
        }
    });

    it('use origin choice 2 exactly where ITA gives two origins, and hexagonal R axes', () => {
        const two = Object.entries(FIXTURE.groups).filter(([, g]) => g.choice === '2').map(([n]) => Number(n));
        expect(new Set(two)).toEqual(ORIGIN_CHOICE_2);
        for (const [n, g] of Object.entries(FIXTURE.groups)) {
            if (itaGenerators(Number(n)).letter === 'R') expect(g.choice).toBe('H');
        }
    });

    it('match the symbols\' lattice letters and the test fixtures', () => {
        for (const f of SPACE_GROUP_FIXTURES) {
            const spec = itaGenerators(f.number);
            expect(spaceGroupNumber(f.symbol)).toBe(f.number);
            expect(spec.letter).toBe(f.symbol[0]);
            expect(spec.centring).toBe(CENTRING[spec.letter]);
            expect(itaOperations(f.number)).toHaveLength(f.multiplicity);
        }
        expect(itaGenerators(0)).toBeNull();
        expect(itaOperations(231)).toBeNull();
    });
});
