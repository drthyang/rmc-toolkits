// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { checkText, toolDescription } from '../copilot/method';
import { buildRunTools } from '../copilot/runTools';

const grid = (n, step) => Array.from({ length: n }, (_, i) => i * step);
const gaussian = (r, center, sigma, height) => height * Math.exp(-((r - center) ** 2) / (2 * sigma ** 2));

const partialsFile = () => {
    const x = grid(400, 0.02);
    return {
        plotKind: 'pdf_partials',
        plotData: {
            series: [
                { label: 'Ta-Se', x, y: x.map((r) => gaussian(r, 2.6, 0.06, 8) + gaussian(r, 4.4, 0.1, 3)) },
                { label: 'Ga-Se', x, y: x.map((r) => gaussian(r, 2.4, 0.05, 6)) },
                { label: 'total', x, y: x.map(() => 1) }
            ]
        }
    };
};

const toolMap = (run) => Object.fromEntries(buildRunTools(run).map((tool) => [tool.name, tool]));

describe('buildRunTools', () => {
    it('offers no tools for an empty run (read_run_data needs at least one section)', () => {
        expect(buildRunTools({})).toEqual([]);
        expect(buildRunTools(null)).toEqual([]);
    });

    it('offers read_run_data with only the sections the run has', async () => {
        const run = {
            plotFiles: [{ plotKind: 'neutron_sq', plotData: { title: 'S(Q)', metrics: { rwp: 0.04 } } }]
        };
        const { read_run_data: tool } = toolMap(run);
        expect(tool.parameters.properties.section.enum).toEqual(['datasets']);
        await expect(tool.run({ section: 'datasets' })).resolves.toEqual({
            data: { section: 'datasets', datasets: [{ kind: 'neutron_sq', title: 'S(Q)', rwp: 0.04 }] },
            checks: []
        });
    });

    it('classifies convergence and names the fit term', async () => {
        const y = Array.from({ length: 200 }, (_, i) => (i < 140 ? 3 - i * (2 / 140) : 1));
        const { assess_convergence: tool } = toolMap({
            rValueFile: { plotData: { chiColumn: 'X_ray_(R)1', series: [{ y }] } }
        });
        expect(tool.page).toBe('dashboard');
        const { data, checks } = await tool.run({});
        expect(data).toMatchObject({ column: 'X_ray_(R)1', status: 'converged', n_steps: 200, first: 3, last: 1 });
        expect(data.total_drop).toBeCloseTo(2, 6);
        expect(data.slope_per_step.last_10pct).toBe(0);
        expect(data.slope_per_step.last_50pct).toBeLessThan(0);
        expect(checks[0].text).toMatch(/X_ray_\(R\)1 term only/);
        expect(checks.some((check) => check.level === 'ok')).toBe(true);
    });

    it('warns on a diverging run with non-finite rows', async () => {
        const y = [2, 1.9, 1.8, null, 1.7, 1.75, 1.9, 2.4];
        const { assess_convergence: tool } = toolMap({ rValueFile: { plotData: { series: [{ y }] } } });
        const { data, checks } = await tool.run({});
        expect(data.status).toBe('diverging');
        expect(data.non_finite_steps).toBe(1);
        expect(checks.filter((check) => check.level === 'warn').map((check) => check.text)).toEqual([
            'Heuristic status: diverging.',
            '1 log rows are non-finite (NaN/Inf).'
        ]);
    });

    it('finds g(r) peaks, the first shell and a constraint-limited first peak', async () => {
        const run = {
            plotFiles: [partialsFile()],
            runSettings: { atoms: ['Ga', 'Se', 'Ta'], minimumDistancesA: [3.0, 2.38, 2.8, 3.0, 2.3, 3.0] }
        };
        const { find_gr_peaks: tool } = toolMap(run);
        // Only element-pair labels are offered.
        expect(tool.parameters.properties.pair.enum).toEqual(['Ta-Se', 'Ga-Se']);

        const tase = await tool.run({ pair: 'Ta-Se' });
        expect(tase.data.peaks_A.map((peak) => peak.r)).toEqual([2.6, 4.4]);
        // One resolved peak in the shell: no first-peak window (never r ± FWHM).
        expect(tase.data.first_peak_window_A).toBeUndefined();
        expect(tase.data.first_shell.first_peak_only).toBeUndefined();
        // The shell runs from the onset to the deep minimum before 4.4 Å.
        const [low, high] = tase.data.first_shell.window_A;
        expect(low).toBeLessThan(2.45);
        expect(high).toBeGreaterThan(2.8);
        expect(high).toBeLessThan(4.1);
        expect(tase.data.first_shell.peaks_in_shell_A).toEqual([2.6]);
        // No cell in this run: no coordination number, never a guessed one.
        expect(tase.data.first_shell.cn).toBeUndefined();
        // Se-Ta's constraint (2.3 Å, labelled Se-Ta by the .dat pair order) is matched to Ta-Se.
        expect(tase.data.min_distance_constraint_A).toBe(2.3);
        expect(tase.checks).toHaveLength(1);
        expect(tase.checks[0]).toMatchObject({ level: 'ok' });
        expect(tase.checks[0].text).toMatch(/^First peak at 2\.6 Å; first shell to [\d.]+ Å, mean 2\.6\d* Å\.$/);

        const gase = await tool.run({ pair: 'Ga-Se' });
        expect(gase.data.min_distance_constraint_A).toBe(2.38);
        expect(gase.checks[0]).toMatchObject({ level: 'warn' });
        expect(gase.checks[0].text).toMatch(/constraint-limited/);
    });

    it('integrates the coordination number, and reads the mean distance from r²g(r)', async () => {
        // A shell of exactly 6 B around each A at 2.6 Å (σ 0.12 Å): g = CN·N(r)/(4πρ r²),
        // so r²g integrates to 6 and has its mean at 2.6 Å, while g itself peaks below it.
        const L = 20;
        const counts = { Ta: 400, Se: 800 };
        const rhoSe = counts.Se / L ** 3;
        const x = grid(600, 0.01);
        const normal = (r) => Math.exp(-((r - 2.6) ** 2) / (2 * 0.12 ** 2)) / (0.12 * Math.sqrt(2 * Math.PI));
        const y = x.map((r) => (r > 0 ? 6 * normal(r) / (4 * Math.PI * rhoSe * r * r) : 0));
        const { find_gr_peaks: tool } = toolMap({
            structure: { latticeVectors: [[L, 0, 0], [0, L, 0], [0, 0, L]], elementCounts: counts },
            plotFiles: [{ plotKind: 'pdf_partials', plotData: { series: [{ label: 'Ta-Se', x, y }] } }]
        });
        const { data, checks } = await tool.run({ pair: 'Ta-Se' });
        expect(data.peaks_A[0].r).toBeLessThan(2.6);
        expect(data.first_shell.mean_d_A).toBeCloseTo(2.6, 2);
        expect(data.first_shell.sd_d_A).toBeCloseTo(0.12, 2);
        expect(data.first_shell.cn.Se_per_Ta).toBeCloseTo(6, 1);
        // Each Se has 6 · ρ(Ta)/ρ(Se) = 3 Ta.
        expect(data.first_shell.cn.Ta_per_Se).toBeCloseTo(3, 1);
        expect(data.first_shell.trough_g_fraction).toBeLessThan(0.02);
        expect(checks.map((check) => check.level)).toEqual(['ok']);
        // The CN of exactly the quoted window: it starts one grid step before
        // the onset (g at 1 % of its peak), so the ~0.2 % of the shell below
        // it is not counted (before 2026-10-19 the tool integrated from r = 0).
        expect(checks[0].text).toMatch(/CN 5\.99 Se per Ta, 3 Ta per Se/);
    });

    it('warns when the first shell is not separated from the next, with the cut sensitivity', async () => {
        // The demo's Ta–Ta: intracluster 3.0 Å, intercluster 4.3 Å, distances filling the gap.
        const L = 20;
        const x = grid(700, 0.01);
        const y = x.map((r) => gaussian(r, 3.0, 0.1, 5) + gaussian(r, 4.3, 0.1, 4) + (r > 3.2 && r < 4.1 ? 0.35 : 0));
        const { find_gr_peaks: tool } = toolMap({
            structure: { latticeVectors: [[L, 0, 0], [0, L, 0], [0, 0, L]], elementCounts: { Ta: 1000 } },
            plotFiles: [{ plotKind: 'pdf_partials', plotData: { series: [{ label: 'Ta-Ta', x, y }] } }]
        });
        const { data, checks } = await tool.run({ pair: 'Ta-Ta' });
        expect(data.first_shell.trough_g_fraction).toBeGreaterThan(0.02);
        const [earlier, later] = data.first_shell.cn_if_cut_0p2A_earlier_later;
        expect(earlier).toBeLessThan(data.first_shell.cn.Ta_per_Ta);
        expect(later).toBeGreaterThan(data.first_shell.cn.Ta_per_Ta);
        expect(checks[0].level).toBe('warn');
        expect(checks[0].text).toMatch(/not separated from the next/);
        expect(checks[0].text).toMatch(/depends on the cut/);
    });

    it('holds a split shell whole and says so', async () => {
        // Like the demo's Ta–Se: a sharp shell at 2.52 Å and a broad one at 2.74 Å, one TaSe6 octahedron.
        const x = grid(400, 0.01);
        const y = x.map((r) => gaussian(r, 2.52, 0.05, 8) + gaussian(r, 2.74, 0.19, 4));
        const { find_gr_peaks: tool } = toolMap({
            plotFiles: [{ plotKind: 'pdf_partials', plotData: { series: [{ label: 'Ta-Se', x, y }] } }]
        });
        const { data, checks } = await tool.run({ pair: 'Ta-Se' });
        const [first, second] = data.peaks_A;
        expect(data.first_shell.peaks_in_shell_A).toEqual([first.r, second.r]);
        // The shell window reaches past the second peak; the first peak alone
        // ends at the minimum between the two (approved 2026-10-19).
        expect(data.first_shell.window_A[1]).toBeGreaterThan(second.r + second.fwhm / 2);
        const alone = data.first_shell.first_peak_only;
        expect(alone.window_A[0]).toBe(data.first_shell.window_A[0]);
        expect(alone.window_A[1]).toBeGreaterThan(first.r);
        expect(alone.window_A[1]).toBeLessThan(second.r);
        expect(alone.trough_g_fraction).toBeGreaterThan(0.2);
        // No cell in this run: no CN for it.
        expect(alone.cn).toBeUndefined();
        expect(data.first_peak_window_A).toBeUndefined();
        expect(checks.map((check) => check.level)).toEqual(['note', 'ok']);
        expect(checks[0].text).toMatch(/first shell holds 2 peaks/);
    });

    it('gives the first peak alone its CN, beside the whole shell\'s, when the run has a cell', async () => {
        // The Demo's shape: 3 Se at 2.52 Å (σ 0.06) and 3 at 2.74 Å (σ 0.10) around each Ta.
        const L = 20;
        const counts = { Ta: 400, Se: 800 };
        const rhoSe = counts.Se / L ** 3;
        const x = grid(600, 0.01);
        const normal = (r, r0, s) => Math.exp(-((r - r0) ** 2) / (2 * s * s)) / (s * Math.sqrt(2 * Math.PI));
        const y = x.map((r) => (r > 0 ? (3 * normal(r, 2.52, 0.06) + 3 * normal(r, 2.74, 0.1)) / (4 * Math.PI * rhoSe * r * r) : 0));
        const { find_gr_peaks: tool } = toolMap({
            structure: { latticeVectors: [[L, 0, 0], [0, L, 0], [0, 0, L]], elementCounts: counts },
            plotFiles: [{ plotKind: 'pdf_partials', plotData: { series: [{ label: 'Ta-Se', x, y }] } }]
        });
        const { data, checks } = await tool.run({ pair: 'Ta-Se' });
        expect(data.first_shell.cn.Se_per_Ta).toBeCloseTo(6, 1);
        const alone = data.first_shell.first_peak_only;
        expect(alone.cn.Se_per_Ta).toBeGreaterThan(3);
        expect(alone.cn.Se_per_Ta).toBeLessThan(4);
        expect(checks[0].text).toBe(`The first shell holds 2 peaks (${data.first_shell.peaks_in_shell_A.join(', ')} Å): a split shell, e.g. two `
            + 'bond lengths of one distorted polyhedron. Use first_shell.window_A, not the first peak alone. The first peak alone '
            + `(first_peak_only, to ${alone.window_A[1]} Å) holds ${alone.cn.Se_per_Ta} Se per Ta of the shell's ${data.first_shell.cn.Se_per_Ta}.`);
    });

    it('joins a split shell across a deep trough when the two peaks make an octahedron (review L1)', async () => {
        // A Jahn–Teller CuF6: 4 × 1.93 + 2 × 2.25 Å, σ 0.05 Å; the old rule cut at the trough: CN 4.
        const L = 20;
        const counts = { Cu: 200, F: 1200 };
        const rhoF = counts.F / L ** 3;
        const x = grid(600, 0.005);
        const normal = (r, r0) => Math.exp(-((r - r0) ** 2) / (2 * 0.05 ** 2)) / (0.05 * Math.sqrt(2 * Math.PI));
        const y = x.map((r) => (r > 0 ? (4 * normal(r, 1.93) + 2 * normal(r, 2.25)) / (4 * Math.PI * rhoF * r * r) : 0));
        const { find_gr_peaks: tool } = toolMap({
            structure: { latticeVectors: [[L, 0, 0], [0, L, 0], [0, 0, L]], elementCounts: counts },
            plotFiles: [{ plotKind: 'pdf_partials', plotData: { series: [{ label: 'Cu-F', x, y }] } }]
        });
        const { data, checks } = await tool.run({ pair: 'Cu-F' });
        expect(data.first_shell.cn.F_per_Cu).toBeCloseTo(6, 1);
        expect(data.first_shell.window_A[1]).toBeGreaterThan(2.35);
        expect(data.first_shell.first_peak_only.cn.F_per_Cu).toBeCloseTo(4, 1);
        expect(data.first_shell.first_peak_only.window_A[1]).toBeLessThan(2.2);
        expect(checks.map((check) => check.level)).toEqual(['note', 'ok']);
        expect(checks[0].text).toMatch(/^One split shell across a deep trough: .* hold 4 F per Cu in the first peak alone .* and 6 together/);
    });

    it('flags a distribution cut at the closest-approach constraint even when its half-maximum clears it', async () => {
        // Ga–Se at 2.42 Å (FWHM ~0.2 Å) with a 2.2 Å constraint: g just above 2.2 Å is ~9 % of the peak.
        const x = grid(500, 0.01);
        const y = x.map((r) => (r < 2.2 ? 0 : gaussian(r, 2.42, 0.1, 8)));
        const { find_gr_peaks: tool } = toolMap({
            runSettings: { atoms: ['Ga', 'Se'], minimumDistancesA: [3.0, 2.2, 3.0] },
            plotFiles: [{ plotKind: 'pdf_partials', plotData: { series: [{ label: 'Ga-Se', x, y }] } }]
        });
        const { data, checks } = await tool.run({ pair: 'Ga-Se' });
        expect(data.peaks_A[0].r - data.peaks_A[0].fwhm / 2).toBeGreaterThan(2.25);
        expect(data.g_at_constraint_fraction).toBeGreaterThan(0.05);
        expect(checks.find((check) => check.level === 'warn').text).toMatch(/2\.2 Å closest-approach constraint cuts this distribution/);
    });

    it('no longer calls a distribution cut because its half maximum lies within 0.05 Å of the constraint (the FWHM clause, dropped)', async () => {
        // A sharp peak on a 0.005 Å grid, 0.07 Å above a 2.43 Å constraint: g
        // just above the constraint is ~1 % of the peak, so nothing piles up.
        const x = grid(2000, 0.005);
        const y = x.map((r) => gaussian(r, 2.5, 0.02, 8));
        const { find_gr_peaks: tool } = toolMap({
            runSettings: { atoms: ['Ga', 'Se'], minimumDistancesA: [3.0, 2.43, 3.0] },
            plotFiles: [{ plotKind: 'pdf_partials', plotData: { series: [{ label: 'Ga-Se', x, y }] } }]
        });
        const { data, checks } = await tool.run({ pair: 'Ga-Se' });
        expect(data.peaks_A[0].r - data.peaks_A[0].fwhm / 2).toBeLessThanOrEqual(2.43 + 0.05);
        expect(data.min_distance_constraint_A).toBe(2.43);
        expect(data.g_at_constraint_fraction).toBeUndefined();
        expect(checks.map((check) => check.level)).toEqual(['ok']);
    });

    it('compare_fit: locates the misfit by range and names the largest residual features', async () => {
        // A G(r) the model reproduces except a bump at 3.5 Å (model above the data).
        const x = grid(1000, 0.02).map((r) => r + 1);
        const exp = x.map((r) => Math.sin(3 * r) * Math.exp(-r / 15));
        const calc = exp.map((value, i) => value + gaussian(x[i], 3.5, 0.05, 0.3));
        const { compare_fit: tool } = toolMap({
            plotFiles: [{
                name: 'run_FT_XFQ1.csv',
                plotKind: 'xpdf',
                plotData: { series: [{ label: 'X_ray-calc', x, y: calc }, { label: 'X_ray_exp_renorm', x, y: exp }] }
            }]
        });
        expect(tool.parameters.properties.dataset.enum).toEqual(['xpdf']);
        const { data, checks } = await tool.run({ dataset: 'xpdf' });
        expect(data.axis).toBe('r');
        expect(data.by_range).toHaveLength(6);
        expect(data.by_range[0].share_of_misfit).toBeGreaterThan(0.95);
        expect(data.largest_residuals[0].r).toBeCloseTo(3.5, 1);
        expect(data.largest_residuals[0].model_minus_data).toBeGreaterThan(0);
        expect(checks).toHaveLength(1);
        expect(checks[0].level).toBe('ok');
        expect(checks[0].text).toMatch(/^Rw [\d.]+ over r = 1–20\.98 Å; the misfit is concentrated in r = 1–4\.33 Å/);
    });

    it('compare_fit: reads the roles from the labels, keys repeated kinds, and flags a poor weak-signal range', async () => {
        // Columns given (data, model): the labels decide. The high-Q third is weak and badly fitted.
        const x = grid(900, 0.03).map((q) => q + 0.5);
        const exp = x.map((q) => (q < 18 ? 2 * Math.sin(q) : 0.01 * Math.sin(q)));
        const calc = exp.map((value, i) => (x[i] < 18 ? value * 1.01 : value * 0.3));
        const file = (name) => ({
            name,
            plotKind: 'neutron_sq',
            plotData: { series: [{ label: 'S(Q)_Expt', x, y: exp }, { label: 'S(Q)_RMC', x, y: calc }] }
        });
        const { compare_fit: tool } = toolMap({ plotFiles: [file('a_SQ1.csv'), file('a_SQ2.csv')] });
        expect(tool.parameters.properties.dataset.enum).toEqual(['neutron_sq', 'neutron_sq_2']);
        const { data, checks } = await tool.run({ dataset: 'neutron_sq_2', n_ranges: 3 });
        expect(data.model_series).toBe('S(Q)_RMC');
        expect(data.data_series).toBe('S(Q)_Expt');
        expect(data.by_range.at(-1).rw).toBeGreaterThan(0.5);
        expect(checks.some((check) => /poorest at Q = /.test(check.text))).toBe(true);
    });

    it('reads the run\'s first shell whatever r_max is, and fails only when the run has none', async () => {
        // No caller picks the range (shells.js): below r_max there is no peak
        // to show, but the shell the context reads is the tool's too.
        const { find_gr_peaks: tool } = toolMap({ plotFiles: [partialsFile()] });
        const { data, checks } = await tool.run({ pair: 'Ta-Se', r_max: 1.5 });
        expect(data.peaks_A).toEqual([]);
        expect(data.first_shell.peaks_in_shell_A).toEqual([2.6]);
        expect(checks[0]).toEqual({
            level: 'note',
            text: 'No peak above 15 % of the maximum below 1.5 Å (r_max); the first shell, read from the peaks below 8 Å, lies beyond it.'
        });
        expect(checks[1].text).toMatch(/^First peak at 2\.6 Å; first shell to /);
        // No peak below 8 Å at all: no shell, and nothing to show below r_max either.
        const x = grid(600, 0.02);
        const { find_gr_peaks: far } = toolMap({
            plotFiles: [{ plotKind: 'pdf_partials', plotData: { series: [{ label: 'Ga-Ga', x, y: x.map((r) => gaussian(r, 9.5, 0.1, 5)) }] } }]
        });
        expect((await far.run({ pair: 'Ga-Ga' })).checks).toEqual([{ level: 'fail', text: 'No peak above 15 % of the maximum below 6 Å.' }]);
        const beyond = await far.run({ pair: 'Ga-Ga', r_max: 12 });
        expect(beyond.data.peaks_A.map((peak) => peak.r)).toEqual([9.5]);
        expect(beyond.data.first_shell).toBeUndefined();
        expect(beyond.checks[0]).toEqual({ level: 'note', text: 'No first shell is read: no peak below 8 Å, where the shells are read.' });
    });
});

describe('the cards each run tool reads', () => {
    it('names the host cards a call is about, from their kind and file', () => {
        const cards = [
            { id: 'model-info', kind: 'model' },
            { id: 'detected-sg', kind: 'symmetry' },
            { id: 'chi2-history', kind: 'r_value', file: 'run-01.log' },
            { id: 'plot:run_PDFpartials.csv', kind: 'pdf_partials', file: 'run_PDFpartials.csv' },
            { id: 'plot:run_SQ1.csv', kind: 'neutron_sq', file: 'run_SQ1.csv' }
        ];
        const x = grid(400, 0.05).map((q) => q + 0.5);
        const fit = {
            name: 'run_SQ1.csv',
            plotKind: 'neutron_sq',
            plotData: { series: [{ label: 'S(Q)_RMC', x, y: x.map(Math.sin) }, { label: 'S(Q)_Expt', x, y: x.map(Math.sin) }] }
        };
        const rValueFile = { name: 'run-01.log', plotData: { series: [{ y: [3, 2.5, 2.2, 2.1, 2.05, 2.0] }] } };
        const tools = Object.fromEntries(buildRunTools(
            { plotFiles: [{ ...partialsFile(), name: 'run_PDFpartials.csv' }, fit], rValueFile },
            { cards }
        ).map((tool) => [tool.name, tool]));
        expect(tools.assess_convergence.cards({})).toEqual(['chi2-history']);
        expect(tools.find_gr_peaks.cards({ pair: 'Ta-Se' })).toEqual(['plot:run_PDFpartials.csv']);
        expect(tools.compare_fit.cards({ dataset: 'neutron_sq' })).toEqual(['plot:run_SQ1.csv']);
        expect(tools.read_run_data.cards({ section: 'convergence' })).toEqual(['chi2-history']);
        expect(tools.read_run_data.cards({ section: 'datasets' })).toEqual([]);
        // Without a catalog, no cards; a file without a name matches none.
        expect(toolMap({ plotFiles: [partialsFile()] }).find_gr_peaks.cards({})).toEqual([]);
    });
});

// The method advice in the run tools (method.js, plan item 1.8): with method
// 'off' (the validation's arm C) a description says what the tool does and
// returns, and a check keeps its reading, without which window to use.
describe('the run tools\' method advice: arm C reads them without it', () => {
    const fitFile = () => {
        const x = grid(400, 0.05).map((q) => q + 0.5);
        return {
            name: 'run_SQ1.csv',
            plotKind: 'neutron_sq',
            plotData: { series: [{ label: 'S(Q)_RMC', x, y: x.map(Math.sin) }, { label: 'S(Q)_Expt', x, y: x.map(Math.sin) }] }
        };
    };

    it('find_gr_peaks and compare_fit: the advice comes back whole with method full, and is gone with off', () => {
        const tools = toolMap({ plotFiles: [partialsFile(), fitFile()] });
        const peaks = tools.find_gr_peaks;
        expect(toolDescription(peaks, 'full')).toContain('the window from onset to the first deep minimum (use first_shell.window_A for '
            + 'compute_bond_angles: it holds the whole shell, split shells included), mean distance');
        expect(toolDescription(peaks, 'off')).toBe(peaks.description);
        expect(peaks.description).toContain('the window from onset to the first deep minimum (it holds the whole shell, '
            + 'split shells included), mean distance');
        expect(peaks.description).not.toMatch(/compute_bond_angles|use first_shell/);
        const fit = tools.compare_fit;
        expect(toolDescription(fit, 'full')).toBe(`${fit.description} Use it for "where does the fit fail" and "is the misfit at short or long range"; `
            + 'run it once per dataset (a G(r) and its F(Q) are separate entries).');
        expect(toolDescription(fit, 'off')).toBe(fit.description);
        expect(fit.description).toMatch(/^Where the model misses the measured data, for one fitted dataset: .* size relative to the data's rms\)\.$/);
    });

    it('find_gr_peaks: a split shell\'s note keeps its reading; "use first_shell.window_A" is advice', async () => {
        const L = 20;
        const counts = { Ta: 400, Se: 800 };
        const rhoSe = counts.Se / L ** 3;
        const x = grid(600, 0.01);
        const normal = (r, r0, s) => Math.exp(-((r - r0) ** 2) / (2 * s * s)) / (s * Math.sqrt(2 * Math.PI));
        const y = x.map((r) => (r > 0 ? (3 * normal(r, 2.52, 0.06) + 3 * normal(r, 2.74, 0.1)) / (4 * Math.PI * rhoSe * r * r) : 0));
        const { find_gr_peaks: tool } = toolMap({
            structure: { latticeVectors: [[L, 0, 0], [0, L, 0], [0, 0, L]], elementCounts: counts },
            plotFiles: [{ plotKind: 'pdf_partials', plotData: { series: [{ label: 'Ta-Se', x, y }] } }]
        });
        const { checks } = await tool.run({ pair: 'Ta-Se' });
        expect(checks[0].method).toEqual([' Use first_shell.window_A, not the first peak alone.']);
        expect(checkText(checks[0], 'full')).toBe(checks[0].text);
        expect(checkText(checks[0], 'off')).toMatch(new RegExp('^The first shell holds 2 peaks \\([\\d., ]+ Å\\): a split shell, e\\.g\\. two '
            + 'bond lengths of one distorted polyhedron\\. The first peak alone \\(first_peak_only, to [\\d.]+ Å\\) holds [\\d.]+ Se per Ta '
            + 'of the shell\'s [\\d.]+\\.$'));
        expect(checkText(checks[1], 'off')).toBe(checks[1].text);
    });

    it('find_gr_peaks: a joined Jahn–Teller shell\'s note keeps "cuts the polyhedron"; when to use the first peak is advice', async () => {
        const L = 20;
        const counts = { Cu: 200, F: 1200 };
        const rhoF = counts.F / L ** 3;
        const x = grid(600, 0.005);
        const normal = (r, r0) => Math.exp(-((r - r0) ** 2) / (2 * 0.05 ** 2)) / (0.05 * Math.sqrt(2 * Math.PI));
        const y = x.map((r) => (r > 0 ? (4 * normal(r, 1.93) + 2 * normal(r, 2.25)) / (4 * Math.PI * rhoF * r * r) : 0));
        const { find_gr_peaks: tool } = toolMap({
            structure: { latticeVectors: [[L, 0, 0], [0, L, 0], [0, 0, L]], elementCounts: counts },
            plotFiles: [{ plotKind: 'pdf_partials', plotData: { series: [{ label: 'Cu-F', x, y }] } }]
        });
        const { checks } = await tool.run({ pair: 'Cu-F' });
        expect(checks[0].text).toMatch(/The first peak alone \(first_peak_only\) cuts the polyhedron; use it only if the user chooses it\.$/);
        expect(checks[0].method).toEqual(['; use it only if the user chooses it']);
        expect(checkText(checks[0], 'off')).toMatch(/^One split shell across a deep trough: .* \(e\.g\. Jahn–Teller\)\. The first peak alone \(first_peak_only\) cuts the polyhedron\.$/);
    });
});
