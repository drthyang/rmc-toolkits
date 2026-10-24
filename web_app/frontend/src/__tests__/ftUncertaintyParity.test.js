// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The browser's G(r) error ribbon (ftUncertainty.js through ftUncertaintyRun.js,
// on a run read as a picked folder is) against the server's
// (rmc_toolkits/ft_uncertainty.py, the `uncertainty` of Flask's /api/plot/data)
// through the golden tests/generate_ft_uncertainty_fixture.py writes; Python
// checks the committed file against its own engine
// (tests/test_ft_uncertainty_fixture.py). So a server folder and the same files
// picked in the browser draw the same ribbon, say the same reason and quote the
// same reproduction residuals.
//
// Tolerances (goldenCompare.js, the rule of every golden): strings, booleans,
// nulls and integers exactly; other numbers to 1e-9 relative + 1e-12 absolute
// (the engines sum in different orders; a reproduction residual of 1e-14 is
// round-off, inside the absolute floor).

import { readFileSync } from 'node:fs';
import { beforeEach, describe, expect, it } from 'vitest';
import fixture from './fixtures/ft_uncertainty_fixture.json';
import pystog from './fixtures/pystog_golden.json';
import { expectGolden } from './goldenCompare';
import { makeRunFromEntries, parseDashboardPlotFile } from '../browserData';
import {
    boxConvolution, convolutionSigmaChange, correlationSummary, fileSigmaColumn, isSigmaLabel, rmcFt, rmcFtCovariance,
    rmcFtSigma, runDatasets, uniformStep, xrayDataset
} from '../ftUncertainty';
import { attachFtUncertainty, clearFtUncertaintyCache, ftUncertaintyForRun } from '../ftUncertaintyRun';

const entry = (folder, name, text) => ({ path: `${folder}/${name}`, file: new File([text], name, { lastModified: 0 }) });

const runOf = async (entries) => {
    const run = await makeRunFromEntries(entries);
    const transform = run.files.find((file) => file.plotKind === 'xpdf');
    return { run, transform };
};

const demoEntries = () => fixture.demo.files.map((name) => entry(
    'Demo', name, readFileSync(new URL(`../../public/demo/${name}`, import.meta.url), 'utf8')
));

beforeEach(() => clearFtUncertaintyCache());

describe('the G(r) error ribbon matches the server (ft_uncertainty_fixture.json)', () => {
    it('the demo run: no ribbon (its input file is not in the folder), the transform reproduced', async () => {
        const { run, transform } = await runOf(demoEntries());
        const payload = await ftUncertaintyForRun(transform, run.allFiles);
        expectGolden(payload, fixture.demo.expected, 'demo');
        expect(payload.reason).toBe('F(Q) carries no uncertainties: scale_ft_rmc.fq (named in GTS_250K.dat) is not in the run folder');
        expect(payload.reproduction.gr.residual).toBeLessThan(1e-13);
    });

    it.each(fixture.cases.map((testCase) => [testCase.name, testCase]))('%s', async (name, testCase) => {
        const entries = Object.entries(testCase.files).map(([file, key]) => entry(fixture.stem, file, fixture.texts[key]));
        const { run, transform } = await runOf(entries);
        const payload = await ftUncertaintyForRun(transform, run.allFiles);
        expectGolden(payload, testCase.expected, name);
    });

    it('the ribbon covers every outcome the server has', () => {
        expect(new Set(fixture.cases.map((testCase) => testCase.expected.status))).toEqual(new Set(['propagated', 'unavailable', 'failed']));
    });
});

describe('attachFtUncertainty: the parsed xPDF plot carries the payload, as Flask sends it', () => {
    it('adds the uncertainty to the transform\'s plotData and keeps the object while nothing changed', async () => {
        const { run } = await runOf(demoEntries());
        const parsed = await Promise.all(run.files.map((file) => (file.plotKind ? parseDashboardPlotFile(file) : file)));
        const first = await attachFtUncertainty(parsed, run.allFiles);
        const transform = first.find((file) => file.plotKind === 'xpdf');
        expectGolden(transform.plotData.uncertainty, fixture.demo.expected, 'attached');
        // Every other file is the same object; the transform's series are unchanged.
        first.forEach((file, index) => { if (file.plotKind !== 'xpdf') expect(file).toBe(parsed[index]); });
        expect(transform.plotData.series).toBe(parsed.find((file) => file.plotKind === 'xpdf').plotData.series);
        const again = await attachFtUncertainty(first, run.allFiles);
        expect(again.find((file) => file.plotKind === 'xpdf')).toBe(transform);
    });

    it('a file that cannot be read is a failure payload, never a throw', async () => {
        const broken = { name: 'x_FT_XFQ1.csv', path: 'run/x_FT_XFQ1.csv', plotKind: 'xpdf', plotData: { series: [] }, sourceFile: { text: () => Promise.reject(new Error('gone')) } };
        const [file] = await attachFtUncertainty([broken], [broken]);
        expect(file.plotData.uncertainty).toMatchObject({ status: 'failed', reason: 'could not compute the error ribbon: gone', sigma: null });
    });
});

describe('the engine port, against closed forms and the adapted pystog golden', () => {
    const q = Float64Array.from({ length: 1181 }, (_, k) => 0.4 + 0.02 * k);
    const r = Float64Array.from({ length: 200 }, (_, j) => 0.1 * (j + 1));
    const shells = [[2.5, 0.6, 0.006], [3.6, 0.45, 0.009], [4.4, 0.35, 0.011]];
    const f = Float64Array.from(q, (value) => {
        let s = 1;
        for (const [d, n, u] of shells) s += ((n * Math.sin(value * d)) / (value * d)) * Math.exp(-0.5 * u * value * value * d);
        return s - 1;
    });
    const dsq = Float64Array.from(q, (value) => 0.002 * (1 + 0.15 * value));

    it('the value is (π/2)·pystog plus the rectangle\'s endpoint halves', () => {
        const got = rmcFt(q, f, r);
        const scale = Math.max(...got.map(Math.abs));
        r.forEach((at, j) => {
            const ends = 0.01 * (q[0] * f[0] * Math.sin(q[0] * at) + q[q.length - 1] * f[f.length - 1] * Math.sin(q[q.length - 1] * at));
            expect(Math.abs(got[j] - (0.5 * Math.PI * pystog.plain.g[j] + ends))).toBeLessThan(1e-12 * scale);
        });
    });

    it('σ is pystog\'s once its endpoint variance weights are reconciled', () => {
        const got = rmcFtSigma(q, dsq, r);
        r.forEach((at, j) => {
            const e0 = (q[0] * dsq[0] * Math.sin(q[0] * at)) ** 2;
            const en = (q[q.length - 1] * dsq[dsq.length - 1] * Math.sin(q[q.length - 1] * at)) ** 2;
            const adapted = Math.sqrt((0.5 * Math.PI * pystog.plain.dg[j]) ** 2 + (0.02 ** 2) * (e0 + en) / 2);
            expect(Math.abs(got[j] / adapted - 1)).toBeLessThan(1e-12);
        });
    });

    it('the covariance diagonal is σ² and white integrand noise on the Nyquist grid is uncorrelated', () => {
        const big = 400;
        const dq = 0.06;
        const qq = Float64Array.from({ length: big + 1 }, (_, k) => dq * k);
        const sigma = Float64Array.from(qq, (value) => (value > 0 ? 0.02 / value : 1));
        const rr = Float64Array.from({ length: 30 }, (_, j) => ((j + 1) * Math.PI) / qq[big]);
        const cov = rmcFtCovariance(qq, sigma, rr);
        const expected = (dq * 0.02) ** 2 * (big / 2);
        const sd = rmcFtSigma(qq, sigma, rr);
        for (let i = 0; i < rr.length; i += 1) {
            expect(cov[i * rr.length + i]).toBe(sd[i] * sd[i]);
            expect(cov[i * rr.length + i] / expected).toBeCloseTo(1, 12);
            for (let j = i + 1; j < rr.length; j += 1) expect(Math.abs(cov[i * rr.length + j])).toBeLessThan(1e-12 * expected);
        }
    });

    it('σ const gives ≈ 5/9 of the Nyquist count of independent points', () => {
        const qq = Float64Array.from({ length: 2001 }, (_, k) => 0.01 * k);
        const rr = Float64Array.from({ length: 401 }, (_, i) => 10 + 0.05 * i);
        const summary = correlationSummary(qq, new Float64Array(qq.length).fill(0.01), rr);
        expect(summary.effectivePoints / summary.nyquistPoints).toBeGreaterThan((5 / 9) * 0.93);
        expect(summary.effectivePoints / summary.nyquistPoints).toBeLessThan((5 / 9) * 1.07);
    });

    it('the box convolution: the Toeplitz form is the dense kernel, and its σ change grows towards R', () => {
        const qq = Float64Array.from({ length: 300 }, (_, k) => 0.8 + 0.05 * k);
        const ff = Float64Array.from(qq, (value) => Math.sin(2.5 * value) / value);
        const radius = 25;
        const step = uniformStep(qq);
        const got = boxConvolution(qq, ff, radius);
        qq.forEach((Q, k) => {
            let sum = 0;
            qq.forEach((Qj, j) => {
                const d = Q - Qj;
                const first = d === 0 ? radius : Math.sin(d * radius) / d;
                sum += (first - Math.sin((Q + Qj) * radius) / (Q + Qj)) * step * Qj * ff[j];
            });
            expect(Math.abs(got[k] - sum / (Math.PI * Q))).toBeLessThan(1e-12);
        });
        const sigma = Float64Array.from(qq, (value) => 0.002 * (1 + 0.15 * value));
        const [inside, edge] = convolutionSigmaChange(qq, sigma, [3, 24], radius);
        expect(Math.abs(inside)).toBeLessThan(0.01);
        expect(Math.abs(edge)).toBeGreaterThan(0.01);
    });

    it('the run-control file and the σ header rule', () => {
        const settings = runDatasets('XRAY_RECIPROCAL_SPACE_DATA ::\n  > FILENAME :: a.fq\n  > CONVOLVE\n  > FITTED_SCALE\nEND ::\n');
        expect(xrayDataset(settings, 1)).toMatchObject({ file: 'a.fq', convolve: true, fittedScale: true, fittedOffset: false });
        expect(xrayDataset(settings, 2)).toBeNull();
        expect(fileSigmaColumn(['r(A)', 'X_ray-calc', 'X_ray_exp_renorm'])).toBeNull();
        expect(fileSigmaColumn(['r(A)', 'X_ray-calc', 'X_ray_exp_renorm', 'sigma'])).toBe(3);
        expect(fileSigmaColumn(['r', 'calc', 'exp', 'dG err'])).toBe(3);
        expect(isSigmaLabel('renorm')).toBe(false);
    });
});
