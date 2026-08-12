// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { readFileSync } from 'node:fs';

import { describe, it, expect } from 'vitest';

import { structureFromRmc6f } from '../browserData.js';
import { describeSymmetry } from '../symmetryModel.js';
import {
    parseCoordinateForm,
    fitsForm,
    wyckoffPositions,
    assignWyckoffLetters,
    tabulatedGroupCount,
} from '../wyckoff.js';
import { WYCKOFF_DATA } from '../wyckoffTable.js';
import { SPACE_GROUPS } from '../spaceGroupTable.js';
import { spaceGroupAtTolerance, siteOrbits } from '../symmetry.js';
import { SPACE_GROUP_FIXTURES, closeGroup, structureFor } from './fixtures/spaceGroups.js';

const cubic = (a) => [[a, 0, 0], [0, a, 0], [0, 0, a]];
const at = (el, x, y, z) => ({ el, frac: [x, y, z] });

describe('coordinate forms', () => {
    it('parses fixed, free and coupled coordinates', () => {
        expect(parseCoordinateForm('0,0,0')).toEqual({ R: [[0, 0, 0], [0, 0, 0], [0, 0, 0]], t: [0, 0, 0] });
        expect(parseCoordinateForm('x,1/4,z')).toEqual({
            R: [[1, 0, 0], [0, 0, 0], [0, 0, 1]], t: [0, 0.25, 0],
        });
        expect(parseCoordinateForm('x,2x,1/4')).toEqual({
            R: [[1, 0, 0], [2, 0, 0], [0, 0, 0]], t: [0, 0, 0.25],
        });
    });

    it('accepts a point of the right form and rejects one that is not', () => {
        const form = parseCoordinateForm('x,1/4,z');
        expect(fitsForm(form, [0.3, 0.25, 0.7])).toBe(true);
        expect(fitsForm(form, [0.3, 0.4, 0.7])).toBe(false);
    });

    it('honours coupling between coordinates', () => {
        const form = parseCoordinateForm('x,x,z');
        expect(fitsForm(form, [0.3, 0.3, 0.7])).toBe(true);
        expect(fitsForm(form, [0.3, 0.5, 0.7])).toBe(false);
    });

    it('fits across the cell boundary', () => {
        // 0.98 and 0.02 are the same coordinate modulo a lattice vector.
        expect(fitsForm(parseCoordinateForm('x,x,z'), [0.98, 0.02, 0.5])).toBe(true);
    });

    it('pins a fixed position', () => {
        const form = parseCoordinateForm('1/2,0,0');
        expect(fitsForm(form, [0.5, 0, 0])).toBe(true);
        expect(fitsForm(form, [0, 0.5, 0.5])).toBe(false);
    });
});

describe('wyckoff table', () => {
    it('covers substantially more than the four groups it used to', () => {
        expect(tabulatedGroupCount()).toBeGreaterThan(200);
    });

    it('only holds real space-group numbers', () => {
        const known = new Set(SPACE_GROUPS.map((g) => g[0]));
        for (const key of Object.keys(WYCKOFF_DATA)) expect(known.has(Number(key))).toBe(true);
    });

    it('parses every row and keeps letters unique within a group', () => {
        for (const key of Object.keys(WYCKOFF_DATA)) {
            const rows = wyckoffPositions(Number(key));
            expect(rows.length, `group ${key}`).toBeGreaterThan(0);
            expect(new Set(rows.map((r) => r.letter)).size, `group ${key}`).toBe(rows.length);
            for (const row of rows) {
                expect(Number.isFinite(row.multiplicity), `group ${key} ${row.letter}`).toBe(true);
                expect(row.site.length, `group ${key} ${row.letter}`).toBeGreaterThan(0);
            }
        }
    });

    it('tops out at the general multiplicity of each group', () => {
        for (const fixture of SPACE_GROUP_FIXTURES) {
            const rows = wyckoffPositions(fixture.number);
            if (!rows.length) continue;
            const largest = Math.max(...rows.map((r) => r.multiplicity));
            expect(largest, fixture.symbol).toBe(fixture.multiplicity);
        }
    });

    // Re-deriving 1724 positions across 230 groups is real work — well past the 5 s default.
    it('every tabulated position reproduces its multiplicity under the group', { timeout: 60000 }, () => {
        // The check the table was built with, re-run here so the committed data cannot
        // drift away from the operations it claims to describe.
        const wrap = (v) => { const y = ((v % 1) + 1) % 1; return y > 1 - 1e-9 ? 0 : y; };
        const near = (a, b) => Math.abs((((a - b) % 1) + 1.5) % 1 - 0.5) < 1e-6;
        const generic = [0.137, 0.213, 0.061];
        for (const fixture of SPACE_GROUP_FIXTURES) {
            const rows = wyckoffPositions(fixture.number);
            if (!rows.length) continue;
            const ops = closeGroup(fixture);
            for (const row of rows) {
                const { R, t } = row.form;
                const rep = [0, 1, 2].map((i) => wrap(R[i][0] * generic[0] + R[i][1] * generic[1] + R[i][2] * generic[2] + t[i]));
                const orbit = [];
                for (const op of ops) {
                    const q = [0, 1, 2].map((i) => wrap(op.R[i][0] * rep[0] + op.R[i][1] * rep[1] + op.R[i][2] * rep[2] + op.t[i]));
                    if (!orbit.some((p) => p.every((v, i) => near(v, q[i])))) orbit.push(q);
                }
                expect(orbit.length, `${fixture.symbol} ${row.multiplicity}${row.letter} (${row.coordinate})`).toBe(row.multiplicity);
            }
        }
    });
});

describe('the bundled demo model', () => {
    // The whole pipeline on a real RMC configuration: 52 000 atoms folded to a 52-site
    // basis, the space group found from it, and each orbit given a Wyckoff letter.
    it('reads GaTa4Se8 as F-43m and labels every orbit', { timeout: 60000 }, () => {
        const text = readFileSync(new URL('../../public/demo/GTS_250K.rmc6f', import.meta.url), 'utf8');
        const structure = structureFromRmc6f({ text });
        expect(structure.basis).toHaveLength(52);

        const symmetry = describeSymmetry(structure, 0.2);
        expect(symmetry).toMatchObject({ spaceGroup: 'F-43m', spaceGroupNumber: 216, nSpace: 96 });

        const labelled = symmetry.orbits.map((o) => `${o.element}:${o.size}${o.wyckoff ?? '?'}`).sort();
        expect(labelled).toEqual(['Ga:4c', 'Se:16e', 'Se:16e', 'Ta:16e']);
    });
});

describe('assigning letters to detected orbits', () => {
    it('separates Pm-3m 3c from 3d, which multiplicity and site symmetry cannot', () => {
        const basis = [
            at('Sr', 0, 0, 0), at('Ti', 0.5, 0.5, 0.5),
            at('O', 0.5, 0.5, 0), at('O', 0.5, 0, 0.5), at('O', 0, 0.5, 0.5),
        ];
        const A = cubic(3.905);
        const found = spaceGroupAtTolerance(A, basis, 0.05);
        const orbits = siteOrbits(A, basis, found.ops, 0.05);
        const letters = assignWyckoffLetters(found.spaceGroupNumber, orbits, basis, 0.02);
        const byElement = Object.fromEntries(orbits.map((o, i) => [o.element, letters[i]]));
        // Both 3c and 3d are multiplicity 3 with site symmetry 4/mmm; only the coordinate
        // form tells them apart, and only one orbit member is written as ITA writes it.
        expect(byElement).toEqual({ Sr: 'a', Ti: 'b', O: 'c' });
    });

    it('labels the rocksalt orbits', () => {
        const fcc = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]];
        const basis = [];
        for (const [x, y, z] of fcc) {
            basis.push(at('Na', x, y, z));
            basis.push(at('Cl', (x + 0.5) % 1, y, z));
        }
        const A = cubic(5.64);
        const found = spaceGroupAtTolerance(A, basis, 0.05);
        const orbits = siteOrbits(A, basis, found.ops, 0.05);
        const letters = assignWyckoffLetters(found.spaceGroupNumber, orbits, basis, 0.02);
        expect(found.spaceGroupNumber).toBe(225);
        expect(letters.filter(Boolean).sort()).toEqual(['a', 'b']);
    });

    it('returns no letter for an untabulated group rather than guessing', () => {
        expect(assignWyckoffLetters(null, [{ size: 1, site: '1', index: [0] }], [at('X', 0, 0, 0)]))
            .toEqual([null]);
    });

    it('assigns the general position of every fixture group', { timeout: 60000 }, () => {
        // A generic orbit is the general position by construction, so its letter must be
        // the last one in the group's table.
        for (const fixture of SPACE_GROUP_FIXTURES) {
            const rows = wyckoffPositions(fixture.number);
            if (!rows.length) continue;
            const { A, basis } = structureFor(fixture);
            const found = spaceGroupAtTolerance(A, basis, 0.05);
            if (found.spaceGroupNumber !== fixture.number) continue;
            const orbits = siteOrbits(A, basis, found.ops, 0.05);
            const letters = assignWyckoffLetters(found.spaceGroupNumber, orbits, basis, 0.02);
            const generalIndex = orbits.findIndex((o) => o.size === fixture.multiplicity && o.site === '1');
            if (generalIndex < 0) continue;
            const expected = rows.reduce((a, r) => (r.multiplicity > a.multiplicity ? r : a), rows[0]);
            expect(letters[generalIndex], `${fixture.symbol} general position`).toBe(expected.letter);
        }
    });
});
