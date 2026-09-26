// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// `.rmc6f` header validation, shared with tests/test_parsers_rmc6f_header.py
// through one fixture: the same values and the same error text in both runtimes.
// A zero supercell, a NaN / collinear / overflowing lattice and Fortran D
// exponents used to give nonsense (NaN folds, singular clouds, all-zero
// results) or a bare NaN instead of a clear error.

import { describe, expect, it } from 'vitest';
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { readRmc6fCellVectors } from '../rmc6f';
import { structureFromRmc6f } from '../browserData';
import { siteDisplacementsFromRmc6f } from '../workers/pcaKde';

const here = dirname(fileURLToPath(import.meta.url));
const { cases } = JSON.parse(readFileSync(join(here, 'fixtures', 'rmc6f_header_cases.json'), 'utf8'));

describe('readRmc6fCellVectors (shared header cases)', () => {
    it('covers the fixture', () => {
        expect(cases.length).toBeGreaterThanOrEqual(15);
    });

    for (const testCase of cases) {
        it(testCase.name, () => {
            const name = `${testCase.name}.rmc6f`;
            if (testCase.error) {
                expect(() => readRmc6fCellVectors(testCase.text, name))
                    .toThrow(testCase.error.replace('{name}', name));
                // Exact text, not just a substring.
                try {
                    readRmc6fCellVectors(testCase.text, name);
                } catch (error) {
                    expect(error.message).toBe(testCase.error.replace('{name}', name));
                }
            } else {
                const { latticeVectors, supercell } = readRmc6fCellVectors(testCase.text, name);
                latticeVectors.forEach((row, i) => row.forEach((value, j) => {
                    expect(value).toBeCloseTo(testCase.latticeVectors[i][j], 12);
                }));
                expect(supercell).toEqual(testCase.supercell);
            }
        });
    }

    it('reads CRLF files the same way', () => {
        const testCase = cases.find((entry) => entry.name === 'fortran_d_exponents');
        const { supercell } = readRmc6fCellVectors(testCase.text.replace(/\n/g, '\r\n'), 'crlf.rmc6f');
        expect(supercell).toEqual([2, 2, 2]);
    });
});

describe('every browser consumer rejects a bad header', () => {
    const text = (name) => cases.find((entry) => entry.name === name).text;

    it('structureFromRmc6f (dashboard / KDE / 3D)', () => {
        expect(() => structureFromRmc6f({ name: 'run.rmc6f', text: text('zero_supercell') }))
            .toThrow('run.rmc6f: supercell dimensions must be three positive integers');
    });

    it('siteDisplacementsFromRmc6f (PCA, orientation and bond angles)', () => {
        expect(() => siteDisplacementsFromRmc6f(text('collinear_lattice'), { name: 'run.rmc6f' }))
            .toThrow('run.rmc6f: lattice vectors are singular (zero cell volume)');
    });
});
