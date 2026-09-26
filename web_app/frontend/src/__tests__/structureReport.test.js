// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Live Data guard of the Dashboard's model summary: hold a read that is
// short of its header count, show one whose only problem is a non-finite atom.

import { describe, expect, it } from 'vitest';
import { parseRmc6fAtoms } from '../rmc6f';
import { isIncompleteStructure } from '../structureReport';

const header = (count) => [
    '(Version 6f format configuration file)',
    `Number of atoms:  ${count}`,
    'Supercell dimensions:  1 1 1',
    'Lattice vectors (Ang):',
    '    4.0 0.0 0.0',
    '    0.0 4.0 0.0',
    '    0.0 0.0 4.0',
    'Atoms:',
];
const line = (id, x) => `  ${id}  Nb  [1]  ${x}  0.25  0.25  1  0  0  0`;
const structureOf = (lines) => ({ parseReport: parseRmc6fAtoms(lines.join('\n')).report });

describe('isIncompleteStructure', () => {
    it('keeps a complete, clean read', () => {
        expect(isIncompleteStructure(structureOf([...header(3), line(1, 0.1), line(2, 0.2), line(3, 0.3)]))).toBe(false);
    });

    it('holds a read with atom lines still missing', () => {
        expect(isIncompleteStructure(structureOf([...header(3), line(1, 0.1), line(2, 0.2)]))).toBe(true);
    });

    it('holds a read whose last line was cut mid-write', () => {
        expect(isIncompleteStructure(structureOf([...header(3), line(1, 0.1), line(2, 0.2), '  3  Nb  [1]  0.3']))).toBe(true);
    });

    it('shows a complete read whose atom blew up to NaN, with its warning', () => {
        for (const bad of ['NaN', '****', 'Infinity']) {
            const structure = structureOf([...header(3), line(1, 0.1), line(2, bad), line(3, 0.3)]);
            expect(structure.parseReport.nonFiniteLines).toBe(1);
            expect(isIncompleteStructure(structure)).toBe(false);
        }
    });

    it('ignores a structure without a header count', () => {
        expect(isIncompleteStructure({})).toBe(false);
        expect(isIncompleteStructure({ parseReport: { declaredAtoms: null, parsedAtoms: 0 } })).toBe(false);
    });
});
