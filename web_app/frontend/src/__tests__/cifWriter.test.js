// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The CIF text of the symmetry-averaged structure (cifWriter.js): CIF 1.1 syntax (ASCII,
// loops whose rows have one value per tag), the spelling CIF readers expect for symbols
// and operations, and a round trip — reading the file back regenerates the model's cell.

import { describe, expect, it } from 'vitest';

import { cifSpaceGroupName, fractionText, rationalText, symopText, writeCif } from '../cifWriter.js';
import { symmetryAveragedStructure } from '../averageStructure.js';
import { toleranceLadder } from '../symmetryModel.js';
import { parseCoordinateForm } from '../wyckoff.js';
import { demoStructure } from './fixtures/symmetryStructures.js';

// Minimal CIF reader for what writeCif emits: data items and loops of whitespace-separated
// values, single-quoted when they hold blanks.
function readCif(text) {
    const tokens = (line) => line.match(/'[^']*'|\S+/g).map((t) => t.replace(/^'|'$/g, ''));
    const items = {};
    const loops = [];
    const lines = text.split('\n');
    for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();
        if (!line || line.startsWith('#')) continue;
        if (line === 'loop_') {
            const tags = [];
            while (lines[i + 1]?.trim().startsWith('_')) tags.push(lines[++i].trim());
            const rows = [];
            while (lines[i + 1]?.trim() && !lines[i + 1].trim().startsWith('_') && lines[i + 1].trim() !== 'loop_') rows.push(tokens(lines[++i]));
            loops.push({ tags, rows });
        } else if (line.startsWith('_')) {
            const [tag, ...rest] = tokens(line);
            items[tag] = rest.join(' ');
        }
    }
    return { items, loops };
}

describe('CIF spelling', () => {
    it.each([
        ['P2_1/c', 'P 21/c'], ['Fd-3m', 'F d -3 m'], ['I4_1/amd', 'I 41/a m d'], ['P6_3mc', 'P 63 m c'],
        ['P-42_1m', 'P -4 21 m'], ['R-3m', 'R -3 m'], ['P2_12_12_1', 'P 21 21 21'], ['Pna2_1', 'P n a 21'],
        ['P3_121', 'P 31 2 1'], ['P1', 'P 1'], ['P-1', 'P -1'], ['F-43m', 'F -4 3 m'], ['Ia-3d', 'I a -3 d'],
    ])('%s → %s', (symbol, expected) => expect(cifSpaceGroupName(symbol)).toBe(expected));

    it('writes translations as reduced fractions', () => {
        expect(fractionText(0)).toBe('');
        expect(fractionText(1 - 1e-9)).toBe('');
        expect(fractionText(0.5)).toBe('1/2');
        expect(fractionText(0.75)).toBe('3/4');
        expect(fractionText(2 / 3)).toBe('2/3');
        expect(fractionText(5 / 48)).toBe('5/48');
        expect(fractionText(-0.25)).toBe('3/4');
        expect(fractionText(0.1234567)).toBe('0.123457');
        expect(rationalText(1.5)).toBe('3/2');
        expect(rationalText(-2)).toBe('-2');
    });

    it('writes operations as coordinate triplets', () => {
        expect(symopText([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0])).toBe('x, y, z');
        expect(symopText([[0, -1, 0], [1, -1, 0], [0, 0, 1]], [0, 0, 0.5])).toBe('-y, x-y, z+1/2');
        expect(symopText([[-1, 0, 0], [0, 1, 0], [0, 0, -1]], [0.5, 0.5, 0.25])).toBe('-x+1/2, y+1/2, -z+1/4');
    });
});

describe('writeCif origin choice', () => {
    it('spells origin choice 2 and hexagonal R axes in the symbol', () => {
        const base = symmetryAveragedStructure(demoStructure(), 0.032);
        const withSetting = (symbol, number, centring, originChoice) => writeCif({
            ...base,
            spaceGroup: { ...base.spaceGroup, symbol, number, centring },
            provenance: { ...base.provenance, itaOrigin: true, originChoice },
        });
        expect(withSetting('Fd-3m', 227, 'F', 2)).toMatch(/_space_group_name_H-M_alt\s+'F d -3 m :2'/);
        expect(withSetting('Fd-3m', 227, 'F', 2)).toMatch(/the ITA standard origin \(origin choice 2\)/);
        expect(withSetting('R-3m', 166, 'R', null)).toMatch(/_space_group_name_H-M_alt\s+'R -3 m :H'/);
        expect(withSetting('F-43m', 216, 'F', null)).toMatch(/_space_group_name_H-M_alt\s+'F -4 3 m'\n/);
    });
});

describe('writeCif', () => {
    const structure = demoStructure();
    const ladder = toleranceLadder(structure, 1.0);
    const exported = ladder.map((brick) => ({
        brick,
        model: symmetryAveragedStructure(structure, brick.from + 1e-6),
    }));

    it('is ASCII CIF 1.1 with one value per tag in every loop row', () => {
        for (const { brick, model } of exported) {
            const text = writeCif(model, { dataName: `GTS_250K_${brick.spaceGroup}`, date: '2026-10-02', tolRange: [brick.from, brick.to] });
            expect(/^[\x20-\x7e\n]*$/.test(text), brick.spaceGroup).toBe(true);
            expect(text).not.toMatch(/NaN|undefined|Infinity/);
            const { items, loops } = readCif(text);
            expect(items._cell_length_a).toMatch(/^\d+\.\d{5}$/);
            for (const { tags, rows } of loops) for (const row of rows) expect(row).toHaveLength(tags.length);
        }
    });

    it('reads back into the same cell: operations, sites and U', () => {
        const { model, brick } = exported[exported.length - 1];
        expect(brick.spaceGroup).toBe('F-43m');
        const text = writeCif(model, { dataName: 'GTS_250K_F-43m', date: '2026-10-02' });
        const { items, loops } = readCif(text);
        expect(items._space_group_IT_number).toBe('216');
        expect(items['_space_group_name_H-M_alt']).toBe('F -4 3 m');
        expect(items._chemical_formula_sum).toBe('Ga Se8 Ta4');
        expect(items._cell_formula_units_Z).toBe('4');

        const symops = loops.find((l) => l.tags.includes('_space_group_symop_operation_xyz'));
        expect(symops.rows).toHaveLength(96);
        symops.rows.forEach(([, xyz], i) => {
            const { R, t } = parseCoordinateForm(xyz);
            expect(R).toEqual(model.operations[i].R);
            t.forEach((v, k) => expect(Math.abs(v - model.operations[i].t[k])).toBeLessThan(1e-12));
        });

        const atoms = loops.find((l) => l.tags.includes('_atom_site_fract_x'));
        const col = (tag) => atoms.tags.indexOf(`_atom_site_${tag}`);
        expect(atoms.rows.map((r) => `${r[col('label')]} ${r[col('symmetry_multiplicity')]}${r[col('Wyckoff_symbol')]}`))
            .toEqual(['Ga1 4c', 'Se1 16e', 'Se2 16e', 'Ta1 16e']);
        atoms.rows.forEach((row, i) => {
            ['fract_x', 'fract_y', 'fract_z'].forEach((tag, k) => expect(Number(row[col(tag)])).toBeCloseTo(model.sites[i].x[k], 6));
            expect(row[col('adp_type')]).toBe('Uani');
        });

        const aniso = loops.find((l) => l.tags.includes('_atom_site_aniso_U_11'));
        aniso.rows.forEach((row, i) => {
            const U = model.sites[i].U;
            [U[0][0], U[1][1], U[2][2], U[0][1], U[0][2], U[1][2]].forEach((v, k) => expect(Number(row[k + 1])).toBeCloseTo(v, 5));
        });
    });

    it('records the provenance a reader needs in the header', () => {
        const { model, brick } = exported[exported.length - 1];
        const text = writeCif(model, { tolRange: [brick.from, brick.to] });
        expect(text).toMatch(/# Source: .*\(52000 atoms, 52 reference sites, supercell 10x10x10\)\./);
        expect(text).toMatch(/holds from 0\.031 to 1\.000 A on the tolerance ladder; orbits taken at 0\.031 A\./);
        expect(text).toMatch(/# Origin: x\(CIF\) = inv\(Q\) \(x\(\.rmc6f\) \+ d\)/);
        expect(text).toMatch(/the ITA standard origin, of the equivalent ones the nearest the \.rmc6f origin/);
        expect(text).toMatch(/The symmetry operations below are those of the ITA standard setting/);
    });

    it('writes a standard cell other than the .rmc6f one with its transformation', () => {
        const { model } = exported.find(({ brick }) => brick.spaceGroup === 'Cmm2');
        const text = writeCif(model);
        expect(text).toMatch(/# Cell: the standard cell the group is named in: a' = [-+abc]+, b' = [-+abc]+, c' = [-+abc]+/);
        expect(readCif(text).items._cell_length_a).toBe('14.65917');
    });

    it('omits the symbol and number when the group has no standard setting', () => {
        const { model } = exported[0];
        const unnamed = {
            ...model,
            spaceGroup: { ...model.spaceGroup, label: 'mm2 class', symbol: null, number: null, standard: false },
            provenance: { ...model.provenance, itaOrigin: false, originChoice: null },
        };
        const text = writeCif(unnamed);
        expect(text).not.toMatch(/_space_group_name_H-M_alt|_space_group_IT_number/);
        expect(text).toMatch(/no standard setting was found for this group/);
        expect(text).toMatch(/The symmetry operations below are authoritative; the origin is not the ITA standard one/);
        expect(readCif(text).loops.find((l) => l.tags.includes('_space_group_symop_operation_xyz')).rows.length).toBeGreaterThan(0);
    });
});
