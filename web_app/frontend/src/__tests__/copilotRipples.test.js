// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// check_ripples (src/copilotRipples.js) renders what the ripple engine decides
// from the run's files and first principles: on the committed Demo end to end
// (the run module's in-process adapter: its ripple source over the Demo's
// files), and on the engine's results for the golden's synthetic runs
// (ripples_fixture.json: a neutron G(r) fitted in r against an untruncated
// model, Lorch-windowed, a Q_max the files do not give, an X-ray F(Q) beside
// two neutron G(r)s, the failures) — the forms the Demo cannot reach.

import { describe, expect, it, vi } from 'vitest';
import fixture from './fixtures/ripples_fixture.json';
import { buildRippleTool, RISK_SENTENCES, summarizeRipples } from '../copilotRipples';
import { checkShape } from '../llm/copilot/dataShape';
import { COPILOT_TOOL_EFFECTS } from '../llm/copilot/effects';
import { openRunFolder } from '../openRunFolder';
import { demoEntries } from './modelVisible/helpers';

const resultOf = (name, index = 0) => fixture.cases.find((testCase) => testCase.name.startsWith(name)).calls[index].expected;
const keysOf = (result) => new Map(result.datasets.map((dataset, k) => [dataset.file, `d${k + 1}`]));
const texts = (checks) => checks.map((check) => `${check.level}: ${check.text}`);

const FIT_IN_Q = 'Fit the data in Q (F(Q)/S(Q)) instead of r, or give RMCProfile the data\'s Q range with CONVOLVE :: Q_max on an unwindowed G(r), so the model is truncated as the data are.';
const GIVE_QMAX = 'Give qmax (and window, if the G(r) was transformed with Lorch\'s) for the G(r) whose Q_max the files do not give.';
const NOT_JUDGED = 'This check does not say whether a given peak or shoulder is a ripple. To test one: a ripple moves with Q_max and a real feature stays, so compare a G(r) transformed with a higher Q_max, a Lorch-windowed transform (which suppresses ripples) or the other radiation\'s data (another Q_max).';
const NEW_START = 'If a feature matters: one run with the model truncated as the data are (or fitted in Q), from a new start; if the feature does not return, it was a ripple.';

const resolutionLine = (key, { period, qmax, half, fwhm, lorch = false }) => `note: ${key}: r-resolution: ripples repeat every ${period} Å `
    + `(2π/Q_max, Q_max = ${qmax} Å⁻¹)${lorch ? ' (Lorch window: ripples suppressed, resolution coarser)' : ''}; a peak at r₀ `
    + `has its ripples between about r₀ ± (n + ½)·${half} Å and r₀ ± (n + 1)·${half} Å (n = 1, 2, …; π/Q_max), near the first `
    + 'for a sharp peak and toward the second for a broad one, their sign set by the peak\'s width; features '
    + `closer than the kernel's FWHM, ${fwhm} Å, are not resolved.`;

describe('check_ripples on the Demo, end to end', () => {
    let tool;

    it('is offered over the Demo\'s fit files, a compute tool, with dataset, qmax and window only', async () => {
        const run = await openRunFolder(demoEntries());
        tool = run.tools.find((entry) => entry.name === 'check_ripples');
        expect(tool).toBeTruthy();
        expect(COPILOT_TOOL_EFFECTS.check_ripples).toBe('compute');
        expect(Object.keys(tool.parameters.properties)).toEqual(['dataset', 'qmax', 'window']);
        expect(tool.parameters.properties.dataset.enum).toEqual(['xray_sq', 'xpdf']);
        expect(tool.cards({})).toEqual(['plot:GTS_250K_FQ1.csv', 'plot:GTS_250K_FT_XFQ1.csv']);
        expect(tool.label({ dataset: 'xpdf' })).toBe('Termination ripples xpdf');
    }, 60000);

    it('says the X-ray G(r) is compared through the same transform, and gives its r-resolution', async () => {
        const { data, checks } = await tool.run({});
        expect(checkShape(data, tool.dataShape)).toEqual([]);
        const xpdf = data.datasets.find((row) => row.dataset === 'xpdf');
        expect(xpdf).toMatchObject({ risk: 'matched', fitted_in: 'r', Q_A: [0.8, 26.98], window: 'none', Q_from: 'the run\'s F(Q) file' });
        expect(xpdf.resolution_A).toEqual({ fwhm: 0.138, nyquist: 0.116, ripple_period: 0.233 });
        expect(Object.keys(xpdf)).toEqual(['dataset', 'file', 'status', 'fitted_in', 'model', 'risk', 'Q_A', 'window', 'Q_from', 'resolution_A']);
        expect(data.verdict).toBe('No real-space fit asks the box to make ripples (each model is truncated as its data are).');
        expect(checks[0]).toEqual({ level: 'ok', text: data.verdict });
        expect(texts(checks)).toContain(resolutionLine('xpdf', { period: 0.233, qmax: 26.98, half: 0.116, fwhm: 0.138 }));
        // A fit in Q has no transform of its own to read: no r-resolution line.
        expect(texts(checks).some((text) => text.startsWith('note: xray_sq: r-resolution'))).toBe(false);
        expect(checks.some((check) => check.level === 'warn')).toBe(false);
        expect(data.next).toEqual([NOT_JUDGED]);
    }, 60000);

    it('passes only the transform arguments to the ripple source', async () => {
        const ripples = vi.fn(async () => ({ datasets: [] }));
        const lone = buildRippleTool({ plotFiles: [{ name: 'a_FT_XFQ1.csv', path: 'run/a_FT_XFQ1.csv', plotKind: 'xpdf' }], ripples });
        const { data, checks } = await lone.run({ qmax: 25, window: 'lorch' });
        expect(ripples).toHaveBeenCalledWith(['run/a_FT_XFQ1.csv'], { qmax: 25, qmin: null, window: 'lorch' }, { signal: undefined });
        expect(data).toEqual({ datasets: [], verdict: 'No fitted dataset to assess.' });
        expect(checks).toEqual([{ level: 'fail', text: 'No fitted dataset to assess.' }]);
    });

    it('is not offered without a ripple source or a fit file', () => {
        expect(buildRippleTool({ plotFiles: [{ name: 'a_FT_XFQ1.csv', plotKind: 'xpdf' }], ripples: null })).toBeNull();
        expect(buildRippleTool({ plotFiles: [{ name: 'a_PDFpartials.csv', plotKind: 'pdf_partials' }], ripples: vi.fn() })).toBeNull();
    });
});

describe('check_ripples renders the engine\'s readings (the golden\'s runs)', () => {
    it('a neutron G(r) fitted in r against the untruncated histogram, Lorch-windowed: warned, with the fix', () => {
        const result = resultOf('a neutron G(r) from a stog input');
        const { data, checks } = summarizeRipples(result, { keys: keysOf(result) });
        expect(data.datasets[0]).toMatchObject({ risk: 'high', fitted_in: 'r', window: 'lorch', Q_from: 'the stog input' });
        expect(texts(checks)).toEqual([
            'ok: The box is asked to fit termination ripples in d1 (fitted in r against an untruncated or differently truncated model).',
            `warn: d1: ${RISK_SENTENCES['untruncated-model']}.`,
            resolutionLine('d1', { period: 0.251, qmax: 25, half: 0.126, fwhm: 0.217, lorch: true })
        ]);
        expect(data.next).toEqual([FIT_IN_Q, NOT_JUDGED, NEW_START]);
    });

    it('no word on Q_max: not assessed, the dataset says to give qmax; given, a window is assumed', () => {
        const result = resultOf('a neutron G(r) with no word on Q_max');
        const { data, checks } = summarizeRipples(result, { keys: keysOf(result) });
        expect(data.datasets[0]).toMatchObject({ status: 'unavailable', risk: 'high' });
        expect(data.datasets[0].why).toMatch(/give qmax/);
        expect(texts(checks).some((text) => text.startsWith('warn: d1: not assessed'))).toBe(true);
        expect(texts(checks).some((text) => text.includes('r-resolution'))).toBe(false);
        expect(data.next).toEqual([FIT_IN_Q, GIVE_QMAX, NEW_START]);
        const given = summarizeRipples(resultOf('a neutron G(r) with no word on Q_max', 1), { keys: keysOf(result) });
        expect(given.data.datasets[0]).toMatchObject({ status: 'assessed', window: 'none (assumed)', Q_from: 'given (qmax)' });
        expect(given.data.next).toEqual([FIT_IN_Q, NOT_JUDGED, NEW_START]);
    });

    it('CONVOLVE with no word on the data\'s transform: whether the box fits ripples cannot be established', () => {
        const result = resultOf('a neutron G(r) with CONVOLVE');
        const { data, checks } = summarizeRipples(result, { keys: keysOf(result) });
        expect(data.verdict).toBe('Whether the box fits ripples in d1 cannot be established from the run\'s files.');
        expect(texts(checks)).toContain(`note: d1: ${RISK_SENTENCES['sinc-data-unverified']}.`);
        expect(data.next).toEqual([NOT_JUDGED]);
    });

    it('an X-ray F(Q) beside two neutron G(r)s: the riskiest first, one line each; one dataset kept on its own', () => {
        const result = resultOf('an X-ray F(Q) fitted in r beside');
        const { data, checks } = summarizeRipples(result, { keys: keysOf(result) });
        expect(data.datasets.map((row) => [row.dataset, row.risk])).toEqual([['d2', 'high'], ['d3', 'unknown'], ['d1', 'matched']]);
        expect(texts(checks).filter((text) => text.includes(': r-resolution:'))).toHaveLength(3);
        const only = summarizeRipples(result, { keys: keysOf(result), only: 'd1' });
        expect(only.data.datasets.map((row) => row.dataset)).toEqual(['d1']);
        expect(only.data.verdict).toMatch(/^No real-space fit asks the box/);
    });

    it('the failures say why, and a fit file with no .dat is not called fitted in r', () => {
        const result = resultOf('the failures');
        const { data } = summarizeRipples(result, { keys: keysOf(result) });
        const why = Object.fromEntries(data.datasets.map((row) => [row.file, row.why ?? null]));
        expect(why['rmc_run_FT_XFQ1.csv']).toMatch(/does not reproduce/);
        expect(why['rmc_run_SQ1.csv']).toMatch(/could not read the run's files: rmc_run_SQ1\.csv/);
        expect(why['other_FT_XFQ1.csv']).toMatch(/no F\(Q\) file/);
        expect(data.datasets.find((row) => row.file === 'other_FT_XFQ1.csv').fitted_in).toBe('unknown');
        expect(data.next).toBeUndefined();
    });

    it('has a sentence for every risk the engine gives', () => {
        const results = [...fixture.demo.calls, ...fixture.cases.flatMap((testCase) => testCase.calls)].map((call) => call.expected);
        for (const result of results) {
            for (const dataset of result.datasets) if (dataset.fit) expect(RISK_SENTENCES[dataset.fit.because], dataset.fit.because).toBeTruthy();
        }
    });
});
