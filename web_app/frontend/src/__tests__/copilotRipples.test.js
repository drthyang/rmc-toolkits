// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// check_ripples (src/copilotRipples.js) renders what the ripple engine decides:
// on the committed Demo end to end (the run module's in-process adapter: its
// ripple source over the Demo's files), and on the engine's results for the
// golden's synthetic runs (ripples_fixture.json: a neutron G(r) fitted in r
// against an untruncated model, an enforced low r, two totals of different
// Q_max, a partial g(r) at a predicted ripple, the failures) — the forms the
// Demo cannot reach.

import { describe, expect, it, vi } from 'vitest';
import fixture from './fixtures/ripples_fixture.json';
import { buildRippleTool, FEATURE_SENTENCES, RISK_SENTENCES, summarizeRipples } from '../copilotRipples';
import { checkShape } from '../llm/copilot/dataShape';
import { COPILOT_TOOL_EFFECTS } from '../llm/copilot/effects';
import { openRunFolder } from '../openRunFolder';
import { demoEntries } from './modelVisible/helpers';

const resultOf = (name, index = 0) => fixture.cases.find((testCase) => testCase.name.startsWith(name)).calls[index].expected;
const keysOf = (result) => new Map(result.datasets.map((dataset, k) => [dataset.file, `d${k + 1}`]));
const texts = (checks) => checks.map((check) => `${check.level}: ${check.text}`);

describe('check_ripples on the Demo, end to end', () => {
    let tool;

    it('is offered over the Demo\'s fit files, a compute tool', async () => {
        const run = await openRunFolder(demoEntries());
        tool = run.tools.find((entry) => entry.name === 'check_ripples');
        expect(tool).toBeTruthy();
        expect(COPILOT_TOOL_EFFECTS.check_ripples).toBe('compute');
        expect(tool.parameters.properties.dataset.enum).toEqual(['xray_sq', 'xpdf']);
        expect(tool.parameters.properties.pair.enum).toContain('Ta-Se');
        expect(tool.cards({})).toEqual(['plot:GTS_250K_FQ1.csv', 'plot:GTS_250K_FT_XFQ1.csv']);
    }, 60000);

    it('says the X-ray G(r) is compared through the same transform: no ripple is fitted', async () => {
        const { data, checks } = await tool.run({});
        expect(checkShape(data, tool.dataShape)).toEqual([]);
        const xpdf = data.datasets.find((row) => row.dataset === 'xpdf');
        expect(xpdf).toMatchObject({ risk: 'matched', fitted_in: 'r', Q_A: [0.8, 26.98], window: 'none' });
        expect(xpdf.resolution_A.fwhm).toBeCloseTo(3.79 / 26.98, 2);
        expect(xpdf.strong_peaks[0].r).toBeCloseTo(2.51, 1);
        expect(xpdf.low_r.state).toBe('measured');
        // The model, transformed alike, has only its own ripples below the first
        // distance: what the strong peaks predict there, within the approximation.
        expect(xpdf.low_r.model_rms / xpdf.low_r.predicted_rms).toBeGreaterThan(0.5);
        expect(xpdf.low_r.model_rms / xpdf.low_r.predicted_rms).toBeLessThan(2);
        expect(data.verdict).toMatch(/^No real-space fit asks the box to make ripples/);
        expect(checks[0].level).toBe('ok');
        expect(checks.some((check) => check.level === 'warn')).toBe(false);
    }, 60000);

    it('reads a feature at a strong peak as a resolved peak, and a split above the resolution', async () => {
        const { data } = await tool.run({ features: [2.51], split: [2.45, 2.6] });
        expect(data.features[0].in[0].class).toBe('resolved-peak');
        expect(data.split.below_resolution_in).toEqual([]);
        expect(data.next).toEqual([expect.stringMatching(/^One G\(r\) only/)]);
    }, 60000);

    it('refuses r_min without r_max before running anything', async () => {
        const ripples = vi.fn();
        const lone = buildRippleTool({ plotFiles: [{ name: 'a_FT_XFQ1.csv', path: 'a_FT_XFQ1.csv', plotKind: 'xpdf' }], ripples });
        const { checks } = await lone.run({ r_min: 2 });
        expect(checks).toEqual([{ level: 'fail', text: 'give r_min and r_max together.' }]);
        expect(ripples).not.toHaveBeenCalled();
    });

    it('is not offered without a ripple source or a fit file', () => {
        expect(buildRippleTool({ plotFiles: [{ name: 'a_FT_XFQ1.csv', plotKind: 'xpdf' }], ripples: null })).toBeNull();
        expect(buildRippleTool({ plotFiles: [{ name: 'a_PDFpartials.csv', plotKind: 'pdf_partials' }], ripples: vi.fn() })).toBeNull();
    });
});

describe('check_ripples renders the engine\'s readings (the golden\'s runs)', () => {
    it('a neutron G(r) fitted in r against the untruncated histogram: warned, with the fix', () => {
        const result = resultOf('a neutron G(r) from a stog input', 1);
        const { data, checks } = summarizeRipples(result, { keys: keysOf(result) });
        expect(data.datasets[0]).toMatchObject({ risk: 'high', fitted_in: 'r', window: 'lorch', Q_from: 'the stog input' });
        expect(texts(checks)).toContain(`warn: d1: ${RISK_SENTENCES['untruncated-model']}.`);
        expect(data.verdict).toMatch(/^The box is asked to fit termination ripples in d1/);
        expect(data.next[0]).toMatch(/^Fit the data in Q/);
        expect(data.split.below_resolution_in).toEqual(['d1']);
    });

    it('a low r set to the baseline (Auto StoG): the level is not read', () => {
        const result = resultOf('a neutron G(r) from Auto StoG');
        const { data, checks } = summarizeRipples(result, { keys: keysOf(result) });
        expect(data.datasets[0].low_r.state).toBe('enforced');
        expect(data.datasets[0].low_r.enforced_below_A).toBeCloseTo(2.0, 6);
        expect(texts(checks).some((text) => text.includes('set to the baseline'))).toBe(true);
    });

    it('no word on Q_max: not assessed, the dataset says to give qmax', () => {
        const result = resultOf('a neutron G(r) with no word on Q_max');
        const { data, checks } = summarizeRipples(result, { keys: keysOf(result) });
        expect(data.datasets[0]).toMatchObject({ status: 'unavailable', risk: 'high' });
        expect(data.datasets[0].why).toMatch(/give qmax/);
        expect(texts(checks).some((text) => text.startsWith('warn: d1: not assessed'))).toBe(true);
        expect(data.next.some((line) => line.startsWith('Give qmax'))).toBe(true);
    });

    it('two totals of different Q_max: each feature across them, and the box\'s partial g(r) at their ripples', () => {
        const result = resultOf('two totals of different Q_max');
        const { data } = summarizeRipples(result, { keys: keysOf(result) });
        const across = data.features.map((row) => row.across?.verdict ?? null);
        expect(across.every((verdict) => typeof verdict === 'string')).toBe(true);
        const pairResult = resultOf('two totals of different Q_max', 1);
        const pair = summarizeRipples(pairResult, { keys: keysOf(pairResult) });
        expect(pair.data.pair.pair).toBe('A-B');
        expect(pair.data.pair.features.some((row) => row.at_ripple.length > 0)).toBe(true);
        expect(pair.checks.some((check) => check.level === 'warn' && check.text.includes('the box may be making it'))).toBe(true);
        // One dataset kept: its rows only, nothing across.
        const only = summarizeRipples(result, { keys: keysOf(result), only: 'd2' });
        expect(only.data.datasets.map((row) => row.dataset)).toEqual(['d2']);
        expect(only.data.features.every((row) => row.in.every((reading) => reading.dataset === 'd2') && !row.across)).toBe(true);
    });

    it('the failures say why, and a fit file with no .dat is not called fitted in r', () => {
        const result = resultOf('the failures');
        const { data } = summarizeRipples(result, { keys: keysOf(result) });
        const why = Object.fromEntries(data.datasets.map((row) => [row.file, row.why ?? null]));
        expect(why['rmc_run_FT_XFQ1.csv']).toMatch(/does not reproduce/);
        expect(why['rmc_run_SQ1.csv']).toMatch(/could not read the run's files: rmc_run_SQ1\.csv/);
        expect(why['other_FT_XFQ1.csv']).toMatch(/no F\(Q\) file/);
        expect(data.datasets.find((row) => row.file === 'other_FT_XFQ1.csv').fitted_in).toBe('unknown');
    });

    it('has a sentence for every reading and every risk the engine gives', () => {
        const results = [...fixture.demo.calls, ...fixture.cases.flatMap((testCase) => testCase.calls)].map((call) => call.expected);
        for (const result of results) {
            for (const dataset of result.datasets) if (dataset.fit) expect(RISK_SENTENCES[dataset.fit.because], dataset.fit.because).toBeTruthy();
            for (const feature of result.features || []) {
                for (const reading of feature.datasets) expect(FEATURE_SENTENCES[reading.class], reading.class).toBeTruthy();
            }
        }
    });
});
