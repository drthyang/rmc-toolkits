// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The browser's average-structure basis (browserData.js, structureFromRmc6f)
// against the server's (rmc_toolkits/structure_basis.py, the `basis` of Flask
// /api/structure) through the golden tests/generate_structure_basis_fixture.py
// writes; Python checks the committed file against its own engine
// (tests/test_structure_basis_fixture.py). The basis is what the Detected SG
// card, its CIF export and the Copilot's detect_symmetry and
// measure_symmetry_lowering read, so a server folder and the same files picked
// in the browser must give the same card, the same CIF and the same answers.
//
// Tolerances: the key sets, strings (el) and integers (referenceNumber, count,
// elementCounts) exactly; every other number to 1e-12 absolute plus 1e-12
// relative. The two runtimes take the same floating-point steps in the same
// order and differ only where their cos/sin/atan2 do (measured on the demo:
// 2.3e-16 in frac, 8e-17 in dispA, 1.3e-19 in covFrac). A changed convention
// (÷ (n − 1), a majority species, a fold about zero) moves values far more.

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { structureFromRmc6f } from '../browserData';
import { brickCif } from '../symmetryCif';
import { describeSymmetry, MAX_SYMMETRY_SITES, toleranceLadder } from '../symmetryModel';
import fixture from './fixtures/structure_basis_fixture.json';

const TOL = 1e-12;
const KEYS = ['count', 'covFrac', 'dispA', 'el', 'elementCounts', 'frac', 'mean', 'referenceNumber'];
const DEMO_PATH = 'Demo/GTS_250K.rmc6f';

const textOf = (source) => {
    if (source.demo) return readFileSync(new URL(`../../public/demo/${source.demo}`, import.meta.url), 'utf8');
    return source.synthetic ? fixture.synthetic : fixture.coordsOnly;
};

const expectNumbers = (actual, expected, path) => {
    if (Array.isArray(expected)) {
        expect(actual, `${path} length`).toHaveLength(expected.length);
        expected.forEach((value, index) => expectNumbers(actual[index], value, `${path}[${index}]`));
        return;
    }
    expect(typeof actual, path).toBe('number');
    expect(Math.abs(actual - expected), `${path}: ${actual} vs ${expected}`).toBeLessThanOrEqual(TOL + TOL * Math.abs(expected));
};

describe('average-structure basis parity (structure_basis_fixture.json)', () => {
    it('the server omits the basis exactly where the finder declines: one site cap', () => {
        // rmc_toolkits.structure_basis.MAX_BASIS_SITES, written by the generator
        // (and held to this file's JS constant by tests/test_structure_basis.py).
        expect(fixture.maxBasisSites).toBe(MAX_SYMMETRY_SITES);
    });

    it('covers the kinds of site the fixture promises', () => {
        const [demo, synthetic, coordsOnly] = fixture.cases.map((entry) => entry.basis);
        expect(demo).toHaveLength(52);
        expect(coordsOnly).toEqual([]);
        const mixed = synthetic.find((site) => site.referenceNumber === 1);
        expect(mixed.el).toBe('Ta');
        expect(mixed.elementCounts).toEqual({ Ta: 12, Nb: 20 });
        expect(synthetic.some((site) => site.count === 1)).toBe(true);
        expect(synthetic.some((site) => site.frac[0] < 0.01 && site.frac[2] > 0.99)).toBe(true);
    });

    for (const entry of fixture.cases) {
        it(entry.name, () => {
            const { basis } = structureFromRmc6f({ path: 'run.rmc6f', text: textOf(entry.source) });
            expect(basis).toHaveLength(entry.basis.length);
            basis.forEach((site, index) => {
                const expected = entry.basis[index];
                const path = `site ${expected.referenceNumber}`;
                expect(Object.keys(site).sort(), path).toEqual(KEYS);
                expect(site.el, `${path}.el`).toBe(expected.el);
                expect(site.referenceNumber, `${path}.referenceNumber`).toBe(expected.referenceNumber);
                expect(site.count, `${path}.count`).toBe(expected.count);
                expect(site.elementCounts, `${path}.elementCounts`).toEqual(expected.elementCounts);
                for (const key of ['frac', 'mean', 'dispA', 'covFrac']) expectNumbers(site[key], expected[key], `${path}.${key}`);
            });
        });
    }
});

describe('a server folder\'s basis gives the browser\'s card and CIF', () => {
    const text = textOf({ demo: 'GTS_250K.rmc6f' });
    const browser = structureFromRmc6f({ path: DEMO_PATH, text });
    // The same structure with the basis Flask's /api/structure returns.
    const server = { ...browser, basis: fixture.cases[0].basis };

    it('the same ladder and the same group, number and orbits at every rung', () => {
        const ladder = toleranceLadder(browser, 1.0);
        const rungs = (bricks) => bricks.map(({ spaceGroup, spaceGroupNumber, nSpace }) => ({ spaceGroup, spaceGroupNumber, nSpace }));
        expect(rungs(toleranceLadder(server, 1.0))).toEqual(rungs(ladder));
        expect(ladder.length).toBeGreaterThan(1);
        // Every orbit's members, site symmetry and Wyckoff label (its
        // representative's position is a basis frac, compared above).
        const orbits = (found) => found.orbits.map(({ rep: _rep, ...orbit }) => orbit);
        for (const brick of ladder) {
            const tolerance = (brick.from + brick.to) / 2;
            const fromServer = describeSymmetry(server, tolerance);
            const fromBrowser = describeSymmetry(browser, tolerance);
            expect(fromServer.spaceGroup).toBe(fromBrowser.spaceGroup);
            expect(fromServer.spaceGroupNumber).toBe(fromBrowser.spaceGroupNumber);
            expect(fromServer.nSpace).toBe(fromBrowser.nSpace);
            expect(fromBrowser.orbits.length).toBeGreaterThan(0);
            expect(orbits(fromServer)).toEqual(orbits(fromBrowser));
        }
    });

    it('the same CIF, byte for byte, for every brick of the ladder', () => {
        for (const brick of toleranceLadder(browser, 1.0)) {
            const fromServer = brickCif(server, brick, { date: '2026-10-19' });
            const fromBrowser = brickCif(browser, brick, { date: '2026-10-19' });
            expect(fromServer.filename).toBe(fromBrowser.filename);
            expect(fromServer.text).toBe(fromBrowser.text);
        }
    });
});
