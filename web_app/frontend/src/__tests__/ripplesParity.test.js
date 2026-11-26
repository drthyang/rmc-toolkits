// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The browser's termination-ripple assessment (ripples.js through
// ripplesRun.js, on a run read as a picked folder is) against the server's
// (rmc_toolkits/ripples.py, what Flask's /api/ripples sends) through the
// golden tests/generate_ripples_fixture.py writes; Python checks the committed
// file against its own engine (tests/test_ripples_fixture.py). So a server
// folder and the same files picked in the browser give the Copilot's
// check_ripples the same fit-space risk, the same transform and the same
// r-resolution.
//
// Tolerances (goldenCompare.js, the rule of every golden): strings, booleans,
// nulls and integers exactly — so every risk, reason and source —; other
// numbers to its RTOL.

import { readFileSync } from 'node:fs';
import { beforeAll, describe, expect, it } from 'vitest';
import fixture from './fixtures/ripples_fixture.json';
import { expectGolden } from './goldenCompare';
import { makeRunFromEntries } from '../browserData';
import { qGrid, resolution, transformWeights } from '../ripples';
import { clearRipplesCache, ripplesForRun } from '../ripplesRun';

const entry = (folder, name, text) => ({ path: `${folder}/${name}`, file: new File([text], name.split('/').pop(), { lastModified: 0 }) });

const argsOf = (call) => ({ qmax: call.qmax ?? null, qmin: call.qmin ?? null, window: call.window ?? null });

const fitsOf = (run, names) => names.map((name) => run.allFiles.find((file) => file.name === name));

// One clean start: the results kept per file versions and arguments.
beforeAll(() => clearRipplesCache());

describe('the ripple assessment matches the server (ripples_fixture.json)', () => {
    it.each(fixture.demo.calls.map((call) => [JSON.stringify(call.args), call]))('the Demo, %s', async (name, call) => {
        const entries = fixture.demo.files.map((file) => entry(
            'Demo', file, readFileSync(new URL(`../../public/demo/${file}`, import.meta.url), 'utf8')
        ));
        const run = await makeRunFromEntries(entries);
        const result = await ripplesForRun(fitsOf(run, fixture.demo.fits), run.allFiles, argsOf(call.args));
        expectGolden(result, call.expected, `demo ${name}`);
    }, 60000);

    const calls = fixture.cases.flatMap((testCase) => testCase.calls.map((call) => [
        `${testCase.name} ${JSON.stringify(call.args)}`, testCase, call
    ]));
    it.each(calls)('%s', async (name, testCase, call) => {
        const entries = Object.entries(testCase.files).map(([file, key]) => entry(fixture.stem, file, fixture.texts[key]));
        const run = await makeRunFromEntries(entries);
        const result = await ripplesForRun(fitsOf(run, testCase.fits), run.allFiles, argsOf(call.args));
        expectGolden(result, call.expected, name);
    }, 60000);

    it('covers every risk, status and transform source the engine has', () => {
        const results = [...fixture.demo.calls, ...fixture.cases.flatMap((testCase) => testCase.calls)].map((call) => call.expected);
        const datasets = results.flatMap((result) => result.datasets);
        expect(new Set(datasets.map((d) => d.fit?.risk).filter(Boolean))).toEqual(new Set(['high', 'matched', 'unknown', 'none']));
        expect(new Set(datasets.map((d) => d.status))).toEqual(new Set(['assessed', 'unavailable', 'failed']));
        const sources = new Set(datasets.map((d) => d.transform?.source).filter(Boolean));
        for (const source of ['run-fq', 'q-file', 'stog', 'provenance', 'convolve', 'argument', 'xray-block']) {
            expect(sources).toContain(source);
        }
        expect(datasets.some((d) => d.transform?.input)).toBe(true);
        expect(datasets.some((d) => d.transform?.window === 'lorch')).toBe(true);
    });
});

describe('known answers in the port', () => {
    const dq = 0.01;

    it('a rectangle from 0 has the FWHM 3.79/Q_max, its main lobe at π/Q_max and the ripple period 2π/Q_max', () => {
        const q = qGrid(dq, dq, 2500);
        const res = resolution(q, transformWeights(q, dq, 'rectangle', 'none'));
        expect(res.fwhm * 25).toBeCloseTo(3.791, 2);
        expect(res.mainLobe).toBeCloseTo(Math.PI / 25, 3);
        expect(res.ripplePeriod).toBeCloseTo(2 * Math.PI / 25, 12);
    });
});
