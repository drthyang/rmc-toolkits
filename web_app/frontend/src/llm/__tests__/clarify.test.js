// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { ambiguousWindows, buildClarifier, readQuestion, resolveClarification } from '../copilot/clarify';
import { shellCatalogue } from '../context/shells';

// A cubic 10 Å cell with 100 A and 100 B atoms, and partial g(r)s built from
// Gaussian shells (normalised to 1 at large r).
const grid = Array.from({ length: 500 }, (_, i) => i * 0.02);
const gaussian = (r, centre, sigma, height) => height * Math.exp(-((r - centre) ** 2) / (2 * sigma ** 2));
const series = (label, shells, floor = 1) => ({
    label,
    x: grid,
    y: grid.map((r) => shells.reduce((sum, [centre, sigma, height]) => sum + gaussian(r, centre, sigma, height), 0)
        + (r > 4.5 ? floor : 0))
});
const structure = { elementCounts: { In: 100, Se: 100 }, latticeVectors: [[10, 0, 0], [0, 10, 0], [0, 0, 10]], supercell: [1, 1, 1] };
const partials = (list) => [{ name: 'p.csv', plotKind: 'pdf_partials', plotData: { series: list } }];

describe('readQuestion', () => {
    const elements = ['In', 'Se', 'Ta', 'S'];

    it('finds triplets and pairs of the run\'s elements, and single elements', () => {
        expect(readQuestion('What are the Se–Ta–Se angles?', elements)).toMatchObject({ angles: true, triplets: [['Se', 'Ta', 'Se']], pairs: [], elements: [] });
        expect(readQuestion('The Ta-Se coordination number', elements)).toMatchObject({ coordination: true, pairs: [['Ta', 'Se']] });
        expect(readQuestion('bond angles around Ta', elements).elements).toEqual(['Ta']);
        expect(readQuestion('Is TaSe6 an octahedron?', elements)).toMatchObject({ polyhedronCn: 6, elements: [] });
    });

    it('reads a word-like symbol only where an element belongs', () => {
        expect(readQuestion('In this run, what are the bond angles?', elements).elements).toEqual([]);
        expect(readQuestion('What are the bond angles around In?', elements).elements).toEqual(['In']);
        expect(readQuestion('Do the S atoms move?', elements).elements).toEqual(['S']);
    });

    it('notices windows, tolerances and site numbers the question names', () => {
        expect(readQuestion('angles in 2.3–2.7 Å', elements).namesAngstrom).toBe(true);
        expect(readQuestion('space group at 0.05 A', elements)).toMatchObject({ symmetry: true, namesAngstrom: true });
        expect(readQuestion('directions of site 3', elements).namesSite).toBe(true);
    });
});

describe('the shell catalog', () => {
    it('finds a split first shell and the coordination of each window', () => {
        const catalog = shellCatalogue({
            plotFiles: partials([series('In-Se', [[2.52, 0.06, 25], [2.74, 0.12, 12]])]),
            structure
        });
        const shell = catalog.shell('In', 'Se');
        expect(shell.peaksInShell).toHaveLength(2);
        const options = ambiguousWindows(shell);
        expect(options.map((option) => option.label)).toEqual(['Whole first shell', 'First peak only']);
        expect(options[0].cn - options[1].cn).toBeGreaterThan(0.15 * options[0].cn);
        // CN of Se around In and of In around Se: the same g, each partner's density.
        expect(catalog.shell('Se', 'In').cn).toBeCloseTo(shell.cn, 5);
    });

    it('asks where a shell runs into the next only when the cut changes the count', () => {
        const merged = shellCatalogue({
            plotFiles: partials([series('In-Se', [[2.6, 0.12, 20], [3.25, 0.2, 9]])]),
            structure
        });
        const options = ambiguousWindows(merged.shell('In', 'Se'));
        expect(options.map((option) => option.label)).toEqual(['Cut at the minimum', 'Cut 0.2 Å earlier', 'Cut 0.2 Å later']);
        const clean = shellCatalogue({ plotFiles: partials([series('In-Se', [[2.6, 0.05, 20]])]), structure });
        expect(ambiguousWindows(clean.shell('In', 'Se'))).toBeNull();
    });
});

describe('the clarifier reads the run\'s one shell reading (shells.js)', () => {
    it('its catalogue is the run\'s, and "First peak only" ends at the minimum between the two peaks', () => {
        const run = { plotFiles: partials([series('In-Se', [[2.52, 0.06, 25], [2.74, 0.12, 12]])]), structure };
        const clarifier = buildClarifier(run);
        expect(clarifier.catalog).toBe(shellCatalogue(run));
        const shell = clarifier.catalog.shell('In', 'Se');
        const [whole, first] = ambiguousWindows(shell);
        expect(whole.window).toEqual(shell.window);
        expect(first.window).toEqual(shell.firstPeak.window);
        expect(first.cn).toBe(shell.firstPeak.cn);
        expect(first.window[1]).toBeGreaterThan(shell.peaksInShell[0]);
        expect(first.window[1]).toBeLessThan(shell.peaksInShell[1]);
    });

    it('asks about a cut window on the CN margin, the same reading the catalogue gives', () => {
        const run = { plotFiles: partials([series('In-Se', [[2.52, 0.06, 25], [2.74, 0.12, 12]])]), structure };
        const clarifier = buildClarifier(run);
        const tool = { shellWindows: (args) => [{ centre: 'In', partner: 'Se', window: [args.lo, args.hi], apply: ([lo, hi]) => ({ ...args, lo, hi }) }] };
        const { firstPeak, window } = clarifier.catalog.shell('In', 'Se');
        expect(clarifier.catalog.window('In', 'Se', firstPeak.window).cnMargin.cut).toBe(true);
        expect(clarifier.beforeCall(tool, { lo: firstPeak.window[0], hi: firstPeak.window[1] }).kind).toBe('cut');
        expect(clarifier.catalog.window('In', 'Se', window).cnMargin.cut).toBe(false);
        expect(clarifier.beforeCall(tool, { lo: window[0], hi: window[1] })).toBeNull();
    });
});

describe('the Workbench\'s questions', () => {
    const run = {
        plotFiles: partials([series('In-Se', [[2.52, 0.06, 25], [2.74, 0.12, 12]]), series('In-In', [[2.9, 0.05, 10]])]),
        structure
    };
    const clarifier = buildClarifier(run);

    it('asks which neighbours for an element alone, nearest first, and settles by a named polyhedron', () => {
        const item = clarifier.beforeTurn('Compute the bond angles around In.');
        expect(item.kind).toBe('neighbours');
        expect(item.request.options.map((option) => option.label)).toEqual(['Se–In–Se', 'In–In–In']);
        expect(clarifier.beforeTurn('What is the coordination number of In?').request.options.at(-1).label).toBe('Every neighbour type');
    });

    it('asks which window for a split shell, unless the question settles it', () => {
        expect(clarifier.beforeTurn('Se–In–Se bond angles?').kind).toBe('window');
        expect(clarifier.beforeTurn('Se–In–Se bond angles in 2.4–2.9 Å?')).toBeNull();
    });

    it('asks which tolerance only when several groups hold above the noise floor', () => {
        const ladder = [
            { spaceGroup: 'P1', from: 0, to: 0.01 },
            { spaceGroup: 'R3m', from: 0.01, to: 0.08 },
            { spaceGroup: 'F-43m', from: 0.08, to: 1 }
        ];
        const item = buildClarifier({ ...run, symmetry: { ladder } }, { noiseFloor: 0.02 }).beforeTurn('What is the space group?');
        expect(item.request.options.map((option) => option.label)).toEqual(['F-43m', 'R3m']);
        expect(item.request.defaultIndex).toBe(0);
        expect(buildClarifier({ ...run, symmetry: { ladder } }, { noiseFloor: 0.07 }).beforeTurn('What is the space group?')).toBeNull();
        expect(buildClarifier({ ...run, symmetry: { ladder } }, { noiseFloor: 0.02 }).beforeTurn('What is the space group at 0.05 Å?')).toBeNull();
    });

    it('asks nothing without the run data it needs', () => {
        const bare = buildClarifier({});
        expect(bare.beforeTurn('Compute the bond angles around In.')).toBeNull();
        expect(bare.beforeCall({ shellWindows: () => [] }, {})).toBeNull();
    });

    it('turns a typed window into the arguments, and keeps the window when the text holds none', () => {
        const tool = { shellWindows: (args) => [{ centre: 'In', partner: 'Se', window: [args.lo, args.hi], apply: ([lo, hi]) => ({ ...args, lo, hi }) }] };
        const item = clarifier.beforeCall(tool, { lo: 2.4, hi: 2.62 });
        expect(item.kind).toBe('cut');
        expect(resolveClarification(item, { answer: 'use 2.3 to 3.0 Å', source: 'other' }, { args: { lo: 2.4, hi: 2.62 }, tool }).args)
            .toEqual({ lo: 2.3, hi: 3 });
        expect(resolveClarification(item, { answer: 'as you like', source: 'other' }, { args: { lo: 2.4, hi: 2.62 }, tool }))
            .toEqual({ args: { lo: 2.4, hi: 2.62 }, note: 'the user answered "as you like"; the window was kept as asked.' });
        // Skipped: the default (the whole shell) stands.
        expect(resolveClarification(item, { source: 'skipped' }, { args: { lo: 2.4, hi: 2.62 }, tool }).note)
            .toMatch(/^the user skipped, so the default stands: run in /);
        // A window covering the whole shell is not asked about.
        expect(clarifier.beforeCall(tool, { lo: 2.2, hi: 3.2 })).toBeNull();
    });
});

describe('a split shell whose bond lengths are well resolved (review L1)', () => {
    // Cu–F: 4 × 1.93 + 2 × 2.25 Å, σ 0.05 Å: the trough between the peaks
    // is deep, so the first minimum alone gave CN 4. Density 0.1 Å⁻³ (the
    // structure above): a narrow peak holds 4π ρ h σ √(2π) r0² atoms.
    const height = (cn, r0, sigma) => cn / (4 * Math.PI * 0.1 * sigma * Math.sqrt(2 * Math.PI) * r0 * r0);
    const jt = series('In-Se', [[1.93, 0.05, height(4, 1.93, 0.05)], [2.25, 0.05, height(2, 2.25, 0.05)]]);
    const run = { plotFiles: partials([jt]), structure };

    it('is one shell of CN 6, with the first peak alone (CN 4, to the trough) as the alternative', () => {
        const shell = shellCatalogue(run).shell('In', 'Se');
        expect(shell.cn).toBeCloseTo(6, 1);
        expect(shell.peaksInShell).toHaveLength(2);
        const options = ambiguousWindows(shell);
        expect(options.map((option) => option.label)).toEqual(['Whole first shell', 'First peak only']);
        expect(options[1].cn).toBeCloseTo(4, 1);
        expect(options[1].window[1]).toBeGreaterThan(2.0);
        expect(options[1].window[1]).toBeLessThan(2.2);
    });

    it('the Workbench asks before an analysis runs in the first peak alone, the whole octahedron the default', () => {
        const clarifier = buildClarifier(run);
        const tool = { shellWindows: (args) => [{ centre: 'In', partner: 'Se', window: [args.lo, args.hi], apply: ([lo, hi]) => ({ ...args, lo, hi }) }] };
        const item = clarifier.beforeCall(tool, { lo: 1.7, hi: 2.1 });
        expect(item.kind).toBe('cut');
        expect(item.request.question).toMatch(/holds 4 Se per In of the shell's 6/);
        expect(item.request.options[item.request.defaultIndex].label).toBe('Whole first shell');
        expect(clarifier.beforeTurn('Se–In–Se bond angles?').kind).toBe('window');
        expect(clarifier.beforeTurn('Is InSe6 a regular octahedron? Se–In–Se angles')).toBeNull();
    });
});
