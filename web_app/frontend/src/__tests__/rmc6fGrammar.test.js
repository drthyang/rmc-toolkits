// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The .rmc6f atom-line grammar, pinned on real-file variants — the same variants
// and expectations as tests/test_parsers_rmc6f_grammar.py, so the browser and the
// Python parser are held to one grammar. The base is the GaNb4Se8 5 K run in
// data/5K_try1 when present locally, else the committed demo run
// public/demo/GTS_250K.rmc6f, cut to its first ATOMS atoms with the header's
// `Number of atoms` adjusted. Expected values come from the base itself.

import { existsSync, readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { structureFromRmc6f } from '../browserData';
import { classifyAtomLine, parseFortranNumber, parseRmc6fAtoms, rmc6fParseWarning } from '../rmc6f';

const REAL_5K = new URL('../../../../data/5K_try1/GaNb4Se8_5K.rmc6f', import.meta.url);
const DEMO = new URL('../../public/demo/GTS_250K.rmc6f', import.meta.url);
const ATOMS = 300;

const base = () => {
    const text = readFileSync(existsSync(REAL_5K) ? REAL_5K : DEMO, 'utf8');
    const lines = text.split('\n');
    const marker = lines.findIndex((line) => line.trim() === 'Atoms:');
    const header = lines.slice(0, marker + 1)
        .map((line) => line.replace(/(Number of atoms:\s*)\d+/, `$1${ATOMS}`));
    return { header, atoms: lines.slice(marker + 1, marker + 1 + ATOMS) };
};

const { header: HEADER, atoms: ATOM_LINES } = base();
const FIRST = ATOM_LINES[0].trim().split(/\s+/);
const FIRST_COORDS = FIRST.slice(3, 6).map(Number);
const FIRST_REF = Number(FIRST[6]);
const SUPERCELL = HEADER.find((line) => line.startsWith('Supercell')).trim().split(/\s+/).slice(-3).map(Number);

const withAtoms = (atomLines, header = HEADER, newline = '\n') => [...header, ...atomLines].join(newline) + newline;
const editTokens = (edit) => ATOM_LINES.map((line) => edit(line.trim().split(/\s+/)).join('   '));
const fortranD = (value) => {
    const [mantissa, exponent] = Number(value).toExponential(14).split('e');
    const sign = exponent.startsWith('-') ? '-' : '+';
    return `${mantissa}D${sign}${exponent.replace(/^[+-]/, '').padStart(2, '0')}`;
};
const eNotation = (value) => Number(value).toExponential(15).toUpperCase();

const VARIANTS = {
    crlf: withAtoms(ATOM_LINES, HEADER, '\r\n'),
    cr_only: withAtoms(ATOM_LINES, HEADER, '\r'),
    tabs: withAtoms(ATOM_LINES.map((line) => line.trim().split(/\s+/).join('\t'))),
    bom: `\uFEFF${withAtoms(ATOM_LINES)}`,
    no_label: withAtoms(editTokens((t) => [...t.slice(0, 2), ...t.slice(3)])),
    split_label: withAtoms(editTokens((t) => [...t.slice(0, 2), '[', t[2].slice(1), ...t.slice(3)])),
    e_notation: withAtoms(editTokens((t) => [...t.slice(0, 3), ...t.slice(3, 6).map(eNotation), ...t.slice(6)])),
    fortran_d: withAtoms(editTokens((t) => [...t.slice(0, 3), ...t.slice(3, 6).map(fortranD), ...t.slice(6)])),
    trailing_blank_lines: withAtoms([...ATOM_LINES, '', '   ', '']),
    marker_space: withAtoms(ATOM_LINES, [...HEADER.slice(0, -1), 'Atoms :']),
    marker_lower: withAtoms(ATOM_LINES, [...HEADER.slice(0, -1), 'atoms:']),
    marker_suffix: withAtoms(ATOM_LINES, [...HEADER.slice(0, -1), 'Atoms (fractional coordinates):']),
    upper_element: withAtoms(editTokens((t) => [t[0], t[1].toUpperCase(), ...t.slice(2)])),
    extra_numeric_column: withAtoms(ATOM_LINES.map((line) => `${line}   2.500000`)),
    trailing_moment_token: withAtoms(ATOM_LINES.map((line) => `${line}   M:  2.500000`)),
    label_without_reference: withAtoms(editTokens((t) => [...t.slice(0, 6), ...t.slice(7)])),
    coords_only: withAtoms(editTokens((t) => t.slice(0, 6))),
};

const CLEAN = [
    'crlf', 'cr_only', 'tabs', 'bom', 'no_label', 'split_label', 'e_notation',
    'fortran_d', 'trailing_blank_lines', 'marker_space', 'marker_lower',
    'marker_suffix', 'upper_element',
];
const ALL_INVALID = ['extra_numeric_column', 'trailing_moment_token', 'label_without_reference'];

const structureOf = (name, text) => structureFromRmc6f({ path: `run/${name}.rmc6f`, text });

describe('.rmc6f layout variants', () => {
    it.each(CLEAN)('%s parses every atom identically', (name) => {
        const { atoms, report } = parseRmc6fAtoms(VARIANTS[name]);
        expect(report.hasAtomsSection).toBe(true);
        expect(report.declaredAtoms).toBe(ATOMS);
        expect(report.parsedAtoms).toBe(ATOMS);
        expect([report.invalidLines, report.nonFiniteLines]).toEqual([0, 0]);
        expect(rmc6fParseWarning(report)).toBeNull();
        expect(atoms).toHaveLength(ATOMS);
        atoms[0].coords.forEach((value, axis) => expect(Math.abs(value - FIRST_COORDS[axis])).toBeLessThan(1e-12));
        expect(atoms[0].referenceNumber).toBe(FIRST_REF);
        expect(atoms[0].element).toBe(FIRST[1].charAt(0).toUpperCase() + FIRST[1].slice(1).toLowerCase());

        const structure = structureOf(name, VARIANTS[name]);
        expect(structure.totalAtoms).toBe(ATOMS);
        expect(structure.parseWarning).toBeNull();
        expect(structure.basis.length).toBeGreaterThan(0);
    });

    it.each(ALL_INVALID)('%s is reported, never column-shifted or silently empty', (name) => {
        // An extra trailing field used to shift every column (y,z → x,y; a cell index
        // became the reference number) while every value stayed finite.
        const { atoms, report } = parseRmc6fAtoms(VARIANTS[name]);
        expect(atoms).toEqual([]);
        expect(report.invalidLines).toBe(ATOMS);
        expect(report.atomLines).toBe(ATOMS);
        const warning = rmc6fParseWarning(report);
        expect(warning).toContain(`${ATOMS} of ${ATOMS} atom lines unparsed`);
        expect(warning).toContain(`parsed 0 of ${ATOMS} atoms declared`);
        expect(() => structureOf(name, VARIANTS[name])).toThrow(/no atoms could be parsed.*unparsed/);
    });

    it('coords-only lines are atoms without reference or cell columns', () => {
        const { atoms, report } = parseRmc6fAtoms(VARIANTS.coords_only);
        expect(report.coordsOnlyAtoms).toBe(ATOMS);
        expect(rmc6fParseWarning(report)).toBeNull();
        expect(atoms[0].referenceNumber).toBeNull();
        expect(atoms[0].cellIndices).toBeNull();
        const structure = structureOf('coords_only', VARIANTS.coords_only);
        expect(structure.totalAtoms).toBe(ATOMS);
        expect(structure.atomIndices).toEqual({});
    });

    it('a truncated file (Live Data mid-write) reports the shortfall', () => {
        const text = withAtoms(ATOM_LINES);
        const cut = text.slice(0, Math.floor(text.length * 0.6));
        const structure = structureOf('truncated', cut);
        expect(structure.totalAtoms).toBeLessThan(ATOMS);
        expect(structure.parseWarning).toContain(`parsed ${structure.totalAtoms} of ${ATOMS} atoms declared in the header`);
    });

    it('non-finite coordinate lines are skipped and counted', () => {
        const lines = [...ATOM_LINES];
        const nan = lines[5].trim().split(/\s+/);
        nan[3] = 'NaN';
        lines[5] = nan.join('   ');
        const overflow = lines[9].trim().split(/\s+/);
        overflow[4] = '*********';
        lines[9] = overflow.join('   ');
        const { atoms, report } = parseRmc6fAtoms(withAtoms(lines));
        expect(atoms).toHaveLength(ATOMS - 2);
        expect(report.nonFiniteLines).toBe(2);
        expect(rmc6fParseWarning(report)).toContain('2 atom lines skipped for non-finite coordinates');
    });

    it('validates reference numbers and cell indices', () => {
        const lines = [...ATOM_LINES];
        const badCell = lines[1].trim().split(/\s+/);
        badCell[7] = String(SUPERCELL[0]);   // outside [0, N_x)
        const badRef = lines[2].trim().split(/\s+/);
        badRef[6] = '0';                     // reference numbers start at 1
        const badInt = lines[3].trim().split(/\s+/);
        badInt[8] = '1.5';                   // cell indices are integers
        lines.splice(1, 3, badCell.join('   '), badRef.join('   '), badInt.join('   '));
        const { atoms, report } = parseRmc6fAtoms(withAtoms(lines));
        expect(atoms).toHaveLength(ATOMS - 3);
        expect(report.invalidLines).toBe(3);
    });

    it('names the file when the lattice metadata is missing', () => {
        expect(() => structureOf('headless', ATOM_LINES.join('\n'))).toThrow('run/headless.rmc6f is missing lattice or supercell metadata');
    });
});

describe('parseFortranNumber', () => {
    it('accepts plain, E and D exponents', () => {
        expect(parseFortranNumber('0.117D-03')).toBe(0.117e-3);
        expect(parseFortranNumber('0.117E-03')).toBe(0.117e-3);
        expect(parseFortranNumber('-.5')).toBe(-0.5);
        expect(parseFortranNumber('12')).toBe(12);
    });

    it('maps non-finite tokens to NaN and rejects non-numbers', () => {
        ['NaN', 'nan', 'Inf', '-Infinity', '********'].forEach((token) => expect(parseFortranNumber(token)).toBeNaN());
        ['abc', '1_000', '0x10', '1.2.3', '', 'M:'].forEach((token) => expect(parseFortranNumber(token)).toBeNull());
    });
});

describe('classifyAtomLine', () => {
    it('keeps the anchored layouts and rejects the rest', () => {
        const split = (line) => line.trim().split(/\s+/);
        expect(classifyAtomLine(split('1 Ga [1] 0.1 0.2 0.3 4 0 0 0')).kind).toBe('atom');
        expect(classifyAtomLine(split('1 Ga 0.1 0.2 0.3 4 0 0 0')).kind).toBe('atom');
        expect(classifyAtomLine(split('1 Ga Ga1 0.1 0.2 0.3 4 0 0 0')).kind).toBe('atom');
        expect(classifyAtomLine(split('1 Ga [1] 0.1 0.2 0.3')).kind).toBe('coords');
        expect(classifyAtomLine(split('1 Ga [1] NaN 0.2 0.3 4 0 0 0')).kind).toBe('nonFinite');
        expect(classifyAtomLine(split('1 Ga [1] 0.1 0.2 0.3 4 0 0 0 2.5')).kind).toBe('invalid');
        expect(classifyAtomLine(split('1 Ga [1 0.1 0.2 0.3 4 0 0 0')).kind).toBe('invalid');
        expect(classifyAtomLine(split('1 Ga [1] 0.1 0.2 0.3 4 0 0 10'), [10, 10, 10]).kind).toBe('invalid');
        expect(classifyAtomLine(split('Magnetic moments: 1 2 3')).kind).toBe('invalid');
    });
});
