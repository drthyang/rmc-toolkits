// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Browser parsers on the COMMITTED demo run (public/demo, GaTa4Se8 250 K), the
// counterpart of tests/test_parsers_demo_run.py. Expected values are read from
// the files themselves with plain string splitting, never from the parser.

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { detectPlotKind, plotDataFromText, structureFromRmc6f } from '../browserData';

const read = (name) => readFileSync(new URL(`../../public/demo/${name}`, import.meta.url), 'utf8');
const RMC6F = read('GTS_250K.rmc6f');
const headerValue = (key) => RMC6F.split('\n').find((line) => line.startsWith(key)).split(':')[1].trim().split(/\s+/);

const csvRows = (text) => text.split('\n').filter((line) => line.trim()).slice(1)
    .map((line) => line.split(',').map((cell) => cell.trim()).filter(Boolean).map(Number));
const conventionalR = (calc, expt) => Math.sqrt(
    calc.reduce((sum, value, index) => sum + (value - expt[index]) ** 2, 0) / expt.reduce((sum, value) => sum + value * value, 0)
);

describe('demo run in the browser parsers', () => {
    it('structure: header composition, supercell and sites', () => {
        const structure = structureFromRmc6f({ path: 'Demo/GTS_250K.rmc6f', text: RMC6F });
        const types = headerValue('Atom types present');
        const counts = headerValue('Number of each atom type').map(Number);
        expect(structure.totalAtoms).toBe(Number(headerValue('Number of atoms')[0]));
        types.forEach((element, index) => expect(structure.elementCounts[element]).toBe(counts[index]));
        expect(structure.supercell).toEqual(headerValue('Supercell dimensions').map(Number));
        expect(structure.parseWarning).toBeNull();
        const sites = new Set(RMC6F.split('\n').slice(RMC6F.split('\n').indexOf('Atoms:') + 1)
            .filter((line) => line.trim()).map((line) => line.trim().split(/\s+/)[6]));
        expect(structure.basis).toHaveLength(sites.size);
        expect(structure.moves.accepted).toBe(Number(headerValue('Number of moves accepted')[0]));
    });

    it.each(['GTS_250K_FQ1.csv', 'GTS_250K_FT_XFQ1.csv'])('%s: Rwp normalized by the experiment', (name) => {
        const text = read(name);
        const rows = csvRows(text);
        const plot = plotDataFromText({ name, plotKind: detectPlotKind(name), text });
        expect(plot.metrics.rwp).toBeCloseTo(conventionalR(rows.map((row) => row[1]), rows.map((row) => row[2])), 12);
        expect(plot.series[0].x).toHaveLength(rows.length);
    });
});
