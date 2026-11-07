// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The browser's termination-ripple assessment (ripples.js through
// ripplesRun.js, on a run read as a picked folder is) against the server's
// (rmc_toolkits/ripples.py, what Flask's /api/ripples sends) through the
// golden tests/generate_ripples_fixture.py writes; Python checks the committed
// file against its own engine (tests/test_ripples_fixture.py). So a server
// folder and the same files picked in the browser give the Copilot's
// check_ripples the same fit-space risk, the same transform and resolution,
// the same strong peaks and predicted ripples, the same low-r level and the
// same reading of every feature.
//
// Tolerances (goldenCompare.js, the rule of every golden): strings, booleans,
// nulls and integers exactly — so every class, risk, state and verdict —;
// other numbers to RIPPLES_RTOL = 1e-7 relative + 1e-12 absolute, not the
// usual 1e-9: the strong peaks are an iterative fit (Levenberg–Marquardt, then
// Gauss–Newton to the stationary point), and a flat minimum is located only to
// about √ε of its scale when the two engines sum in different orders (seen:
// 1.5e-9). The generator refuses a reading within 1e-6 of a threshold.

import { readFileSync } from 'node:fs';
import { beforeEach, describe, expect, it } from 'vitest';
import fixture from './fixtures/ripples_fixture.json';
import { expectGolden } from './goldenCompare';
import { makeRunFromEntries } from '../browserData';
import {
    analyseSeries, applyTransform, classifyFeature, gaussianSine, qGrid, resolution, ripplePattern, transformWeights
} from '../ripples';
import { clearRipplesCache, ripplesForRun } from '../ripplesRun';

const entry = (folder, name, text) => ({ path: `${folder}/${name}`, file: new File([text], name.split('/').pop(), { lastModified: 0 }) });

const argsOf = (call) => ({
    features: call.features ?? null, rRange: call.rRange ?? null, split: call.split ?? null, pair: call.pair ?? null,
    qmax: call.qmax ?? null, qmin: call.qmin ?? null, window: call.window ?? null
});

const RIPPLES_RTOL = 1e-7;

const fitsOf = (run, names) => names.map((name) => run.allFiles.find((file) => file.name === name));

beforeEach(() => clearRipplesCache());

describe('the ripple assessment matches the server (ripples_fixture.json)', () => {
    it.each(fixture.demo.calls.map((call) => [JSON.stringify(call.args), call]))('the Demo, %s', async (name, call) => {
        const entries = fixture.demo.files.map((file) => entry(
            'Demo', file, readFileSync(new URL(`../../public/demo/${file}`, import.meta.url), 'utf8')
        ));
        const run = await makeRunFromEntries(entries);
        const result = await ripplesForRun(fitsOf(run, fixture.demo.fits), run.allFiles, argsOf(call.args));
        expectGolden(result, call.expected, `demo ${name}`, { rtol: RIPPLES_RTOL });
    }, 60000);

    const calls = fixture.cases.flatMap((testCase) => testCase.calls.map((call) => [
        `${testCase.name} ${JSON.stringify(call.args)}`, testCase, call
    ]));
    it.each(calls)('%s', async (name, testCase, call) => {
        const entries = Object.entries(testCase.files).map(([file, key]) => entry(fixture.stem, file, fixture.texts[key]));
        const run = await makeRunFromEntries(entries);
        const result = await ripplesForRun(fitsOf(run, testCase.fits), run.allFiles, argsOf(call.args));
        expectGolden(result, call.expected, name, { rtol: RIPPLES_RTOL });
    }, 60000);

    it('covers every risk, every low-r state and every reading the engine has', () => {
        const results = [...fixture.demo.calls, ...fixture.cases.flatMap((testCase) => testCase.calls)].map((call) => call.expected);
        const datasets = results.flatMap((result) => result.datasets);
        expect(new Set(datasets.map((d) => d.fit?.risk).filter(Boolean))).toEqual(new Set(['high', 'matched', 'unknown', 'none']));
        expect(new Set(datasets.map((d) => d.status))).toEqual(new Set(['assessed', 'unavailable', 'failed']));
        const states = new Set(datasets.map((d) => d.lowR?.state).filter(Boolean));
        for (const state of ['measured', 'enforced']) expect(states).toContain(state);
        const classes = new Set(results.flatMap((result) => (result.features || []).flatMap((f) => f.datasets.map((d) => d.class))));
        for (const kind of ['ripple-consistent', 'larger-than-ripple', 'no-ripple-here']) expect(classes).toContain(kind);
        const sources = new Set(datasets.map((d) => d.transform?.source).filter(Boolean));
        for (const source of ['run-fq', 'q-file', 'stog', 'provenance', 'convolve', 'argument', 'xray-block']) {
            expect(sources).toContain(source);
        }
        expect(results.some((result) => result.pair)).toBe(true);
        expect(results.some((result) => (result.features || []).some((f) => f.across))).toBe(true);
    });
});

describe('known answers in the port', () => {
    const dq = 0.01;
    const peaks = [{ amplitude: 8, r: 2.5, sigma: 0.06 }, { amplitude: 6, r: 3.6, sigma: 0.08 }, { amplitude: 10, r: 4.4, sigma: 0.1 }];
    const r = Float64Array.from({ length: 1200 }, (_, i) => 0.01 * (i + 1));

    it('a rectangle from 0 has the FWHM 3.79/Q_max and its main lobe at π/Q_max', () => {
        const q = qGrid(dq, dq, 2500);
        const res = resolution(q, transformWeights(q, dq, 'rectangle', 'none'));
        expect(res.fwhm * 25).toBeCloseTo(3.791, 2);
        expect(res.mainLobe).toBeCloseTo(Math.PI / 25, 3);
    });

    it('recovers Gaussian peaks through the transform and predicts their ripples', () => {
        const q = qGrid(dq, dq, 2000);
        const weights = transformWeights(q, dq, 'rectangle', 'none');
        const spectrum = new Float64Array(q.length);
        for (const peak of peaks) gaussianSine(q, peak).forEach((value, k) => { spectrum[k] += value; });
        const y = applyTransform(q, weights, spectrum, r);
        const analysis = analyseSeries(r, y, { space: 'D', qMin: dq, qMax: 20, dq, rule: 'rectangle', window: 'none', rFirst: 2.3 });
        expect(analysis.peaks.map((peak) => Math.round(peak.r * 100) / 100)).toEqual([2.5, 3.6, 4.4]);
        analysis.peaks.forEach((peak, k) => {
            expect(Math.abs(peak.sigma - peaks[k].sigma)).toBeLessThan(0.002);
            const single = ripplePattern(q, weights, dq, 'rectangle', peaks[k], r);
            for (const ripple of peak.ripples) {
                const actual = single[Math.round(ripple.r / 0.01) - 1];
                expect(Math.abs(ripple.amplitude - actual)).toBeLessThan(0.02 * Math.abs(actual));
                expect(classifyFeature(analysis, ripple.r).class).toBe('ripple-consistent');
            }
        });
    }, 60000);
});
