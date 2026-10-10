// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
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

    it('finds g(r) peaks, a first-peak window and a constraint-limited first peak', async () => {
        const run = {
            plotFiles: [partialsFile()],
            runSettings: { atoms: ['Ga', 'Se', 'Ta'], minimumDistancesA: [3.0, 2.38, 2.8, 3.0, 2.3, 3.0] }
        };
        const { find_gr_peaks: tool } = toolMap(run);
        // Only element-pair labels are offered.
        expect(tool.parameters.properties.pair.enum).toEqual(['Ta-Se', 'Ga-Se']);

        const tase = await tool.run({ pair: 'Ta-Se' });
        expect(tase.data.peaks_A.map((peak) => peak.r)).toEqual([2.6, 4.4]);
        expect(tase.data.first_peak_window_A[0]).toBeLessThan(2.6);
        expect(tase.data.first_peak_window_A[1]).toBeGreaterThan(2.6);
        // Se-Ta's constraint (2.3 Å, labelled Se-Ta by the .dat pair order) is matched to Ta-Se.
        expect(tase.data.min_distance_constraint_A).toBe(2.3);
        expect(tase.checks).toEqual([{ level: 'ok', text: 'First peak at 2.6 Å.' }]);

        const gase = await tool.run({ pair: 'Ga-Se' });
        expect(gase.data.min_distance_constraint_A).toBe(2.38);
        expect(gase.checks[0]).toMatchObject({ level: 'warn' });
        expect(gase.checks[0].text).toMatch(/constraint-limited/);
    });

    it('fails clearly when no peak is found below r_max', async () => {
        const { find_gr_peaks: tool } = toolMap({ plotFiles: [partialsFile()] });
        const { checks } = await tool.run({ pair: 'Ta-Se', r_max: 1.5 });
        expect(checks).toEqual([{ level: 'fail', text: 'No peak above 15 % of the maximum below 1.5 Å.' }]);
    });
});
