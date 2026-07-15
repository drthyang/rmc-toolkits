// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import {
    buildLocalRun,
    detectPlotKind,
    fitTypeByFilename,
    parseRunSettings,
    plotDataFromText,
    plotMetadataFromFile,
    structureFromRmc6f,
} from '../browserData';

// Minimal .rmc6f fixture: 2×1×1 supercell of a 10 Å cubic cell, one Ga
// reference site whose two supercell copies sit at within-cell fraction
// 0.5 ± delta on x (exact 0.25 / 0.75 on y / z). Box coords are
// (withinFrac + cellIndex) / supercell.
const rmc6fText = (delta) => [
    '(Version 6f format configuration file)',
    'Number of atoms: 2',
    'Number of moves generated:           2689753',
    'Number of moves tried:               2685185',
    'Number of moves accepted:            586675',
    'Accumulated time (s) in running loop: 24958824.16',
    'Supercell dimensions: 2 1 1',
    'Lattice vectors (Ang):',
    ' 20.0 0.0 0.0',
    ' 0.0 10.0 0.0',
    ' 0.0 0.0 10.0',
    'Atoms:',
    ` 1 Ga [1] ${(0.5 + delta) / 2} 0.25 0.75 1 0 0 0`,
    ` 2 Ga [1] ${(0.5 - delta + 1) / 2} 0.25 0.75 1 1 0 0`
].join('\n');

describe('structureFromRmc6f site displacement (dispA)', () => {
    it('recovers the mean site and the rms displacement about it', () => {
        const structure = structureFromRmc6f({ path: 'test.rmc6f', text: rmc6fText(0.02) });
        expect(structure.basis).toHaveLength(1);
        const site = structure.basis[0];
        expect(site.el).toBe('Ga');
        expect(site.frac[0]).toBeCloseTo(0.5, 6);
        expect(site.frac[1]).toBeCloseTo(0.25, 6);
        expect(site.frac[2]).toBeCloseTo(0.75, 6);
        // Two points at ±0.02 in cell fraction on a 10 Å edge: circular std
        // √(−2 ln cos(2πδ))/2π · 10 ≈ 0.2003 Å; y and z contribute zero.
        expect(site.dispA).toBeCloseTo(0.2003, 3);
    });

    it('reports zero displacement for coincident copies', () => {
        const structure = structureFromRmc6f({ path: 'test.rmc6f', text: rmc6fText(0) });
        expect(structure.basis[0].dispA).toBe(0);
    });

    it('handles a boundary-wrapping site (mean at 0 ≡ 1)', () => {
        const text = [
            'Supercell dimensions: 2 1 1',
            'Lattice vectors (Ang):',
            ' 20.0 0.0 0.0',
            ' 0.0 10.0 0.0',
            ' 0.0 0.0 10.0',
            'Atoms:',
            ' 1 Ga [1] 0.495 0.25 0.75 1 0 0 0',   // within-cell 0.99
            ' 2 Ga [1] 0.505 0.25 0.75 1 1 0 0'    // within-cell 0.01
        ].join('\n');
        const site = structureFromRmc6f({ path: 'test.rmc6f', text }).basis[0];
        const wrapped = Math.min(site.frac[0], 1 - site.frac[0]);
        expect(wrapped).toBeLessThan(1e-6);
        expect(site.dispA).toBeCloseTo(0.1, 3);   // ±0.01 frac on a 10 Å edge
    });
});

describe('rmc6f moves metadata', () => {
    it('parses the header move counters', () => {
        const structure = structureFromRmc6f({ path: 'test.rmc6f', text: rmc6fText(0.02) });
        expect(structure.moves).toEqual({
            generated: 2689753,
            tried: 2685185,
            accepted: 586675,
            accumulatedTimeS: 24958824.16
        });
    });

    it('returns null moves when the header has no counters', () => {
        const text = rmc6fText(0).split('\n').filter((line) => !line.includes('moves') && !line.includes('Accumulated')).join('\n');
        expect(structureFromRmc6f({ path: 'test.rmc6f', text }).moves).toBeNull();
    });
});

describe('parseRunSettings (<stem>.dat run-control file)', () => {
    const datText = [
        'TITLE :: GaNb4Se8 5K',
        'MATERIAL :: GaNb4Se8',
        'PHASE :: cubic',
        'TEMPERATURE :: 5K',
        'INVESTIGATOR :: Tsung-Han',
        '',
        'NUMBER_DENSITY :: 0.046293 Angstrom^(-3)',
        'MINIMUM_DISTANCES :: 2.2 2.2 2.2 2.2 2.2 2.2 Angstrom',
        'MAXIMUM_MOVES :: 0.05 0.05 0.05 Angstrom',
        'TIME_LIMIT :: 10000000.00 MINUTES',
        'SAVE_PERIOD :: 10.00 MINUTES',
        'WEIGHT_OPTIMIZATION ::',
        '',
        'CURVATURE :: 1 2.0 2.3 1',
        '',
        'ATOMS :: Ga Nb Se',
        '',
        'FLAGS ::',
        '  > NO_MOVEOUT',
        '  > NO_SAVE_CONFIGURATIONS',
        '',
        'XRAY_RECIPROCAL_SPACE_DATA ::',
        '  > FILENAME :: scale_ft_rmc.fq',
        '  > DATA_TYPE :: F(Q)',
        '  > FIT_TYPE ::  F(Q)',
        '  > FITTED_SCALE',
        '',
        'END ::'
    ].join('\n');

    it('extracts the settings the assistant needs', () => {
        const settings = parseRunSettings(datText);
        expect(settings.title).toBe('GaNb4Se8 5K');
        expect(settings.temperature).toBe('5K');
        expect(settings.atoms).toEqual(['Ga', 'Nb', 'Se']);
        expect(settings.minimumDistancesA).toEqual([2.2, 2.2, 2.2, 2.2, 2.2, 2.2]);
        expect(settings.maximumMovesA).toEqual([0.05, 0.05, 0.05]);
        expect(settings.timeLimit).toBe('10000000.00 MINUTES');
        expect(settings.weightOptimization).toBe(true);
        expect(settings.flags).toEqual(['NO_MOVEOUT', 'NO_SAVE_CONFIGURATIONS']);
        expect(settings.datasets).toEqual([
            { block: 'XRAY_RECIPROCAL_SPACE_DATA', file: 'scale_ft_rmc.fq', fit_type: 'F(Q)' }
        ]);
    });

    it('returns null for non-control .dat files (e.g. numeric chi2.dat)', () => {
        expect(parseRunSettings('0.1 2.3 4.5\n0.2 2.1 4.4\n')).toBeNull();
        expect(parseRunSettings('')).toBeNull();
    });
});

describe('settings file selection (buildLocalRun)', () => {
    const makeFile = (name) => new File(['x'], name, { lastModified: 1 });
    const withPath = (name) => {
        const file = makeFile(name);
        Object.defineProperty(file, 'webkitRelativePath', { value: `run/${name}` });
        return file;
    };

    it('picks the .dat matching the structure stem, ignoring auxiliary .dat files', async () => {
        const run = await buildLocalRun([
            withPath('GaNb4Se8_5K.rmc6f'),
            withPath('GaNb4Se8_5K_FQ1.csv'),
            withPath('GaNb4Se8_5K.dat'),
            withPath('chi2.dat'),
            withPath('optimization.dat'),
            withPath('weights_update.dat')
        ]);
        expect(run.structureFile.name).toBe('GaNb4Se8_5K.rmc6f');
        expect(run.settingsFile.name).toBe('GaNb4Se8_5K.dat');
        expect(run.settingsFile.path).toBe('run/GaNb4Se8_5K.dat');
    });

    it('has no settings file when no stem-matched .dat exists', async () => {
        const run = await buildLocalRun([
            withPath('GaNb4Se8_5K.rmc6f'),
            withPath('GaNb4Se8_5K_FQ1.csv'),
            withPath('chi2.dat')
        ]);
        expect(run.settingsFile).toBeNull();
    });
});

describe('run-control fit-function labels on Dashboard plots', () => {
    it('detects any .gr / .sq / .fq file as a STOG plot, not just scale_ft.*', () => {
        expect(detectPlotKind('PMN_300k_rmc_v2.gr')).toBe('stog');
        expect(detectPlotKind('PMN_300k.sq')).toBe('stog');
        expect(detectPlotKind('data.fq')).toBe('stog');
        expect(detectPlotKind('scale_ft.gr')).toBe('stog');
        expect(detectPlotKind('something.csv')).not.toBe('stog');
    });

    it('parses FIT_TYPE per dataset and maps it by filename (fit type wins over data type)', () => {
        const dat = [
            'TITLE :: PMN',
            'NEUTRON_REAL_SPACE_DATA :: 3',
            '  > FILENAME :: PMN_300k_rmc_30_100_lOWTfurn_v2.gr',
            '  > DATA_TYPE :: G(r)',
            '  > FIT_TYPE :: D(r)',
            'XRAY_RECIPROCAL_SPACE_DATA ::',
            '  > FILENAME :: PMN_XFQ.fq',
            '  > DATA_TYPE :: F(Q)',
            '  > FIT_TYPE ::  F(Q)',
            'END ::',
        ].join('\n');
        const map = fitTypeByFilename(parseRunSettings(dat));
        expect(map.get('pmn_300k_rmc_30_100_lowtfurn_v2.gr')).toBe('D(r)');
        expect(map.get('pmn_xfq.fq')).toBe('F(Q)');
    });

    it('labels a STOG plot by its fit type when paired, else by extension', () => {
        const stogText = 'STOG header\n2\n0.00 1.0 1.1\n0.10 2.0 2.1\n';
        const withFit = { plotKind: 'stog', name: 'PMN_v2.gr', fitType: 'D(r)', text: stogText };
        const noFit = { plotKind: 'stog', name: 'PMN_v2.gr', text: stogText };
        expect(plotMetadataFromFile(withFit).title).toBe('D(r)');
        expect(plotDataFromText(withFit).yLabel).toBe('D(r)');
        expect(plotDataFromText(withFit).xLabel).toBe('r (Å)');
        expect(plotMetadataFromFile(noFit).title).toBe('G(r)');
        expect(plotDataFromText(noFit).yLabel).toBe('G(r)');
    });

    it('pairs the fit type from the stem-matched .dat onto the loaded .gr plot file (buildLocalRun)', async () => {
        const dat = [
            'TITLE :: PMN',
            'NEUTRON_REAL_SPACE_DATA :: 3',
            '  > FILENAME :: PMN_v2.gr',
            '  > FIT_TYPE :: D(r)',
            'END ::',
        ].join('\n');
        const withPath = (name, content) => {
            const file = new File([content], name, { lastModified: 1 });
            Object.defineProperty(file, 'webkitRelativePath', { value: `run/${name}` });
            return file;
        };
        const run = await buildLocalRun([
            withPath('PMN.rmc6f', 'x'),
            withPath('PMN.dat', dat),
            withPath('chi2.dat', '0.1 0.2\n0.3 0.4\n'),
            withPath('PMN_v2.gr', 'h1\nh2\n0 1\n0.1 2\n'),
        ]);
        const gr = run.files.find((file) => file.name === 'PMN_v2.gr');
        expect(gr).toBeTruthy();
        expect(gr.fitType).toBe('D(r)');
        expect(plotMetadataFromFile(gr).title).toBe('D(r)');
    });
});
