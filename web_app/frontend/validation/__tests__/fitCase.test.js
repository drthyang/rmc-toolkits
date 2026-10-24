// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The fit-data dev case (plan item 4.6, docs/AGENT_PLAN.md): dev-fit-unconverged's
// run folder carries RMCProfile's X-ray fit files besides its .rmc6f and
// partials, the run module opens them as it opens a user's run, and the fit
// tools, the context and the answer checks find the planted fit truth: a χ²
// history still falling, and a G(r) peak the configuration does not reproduce.
// The case is scored on its structure as every case is (no order, Fm-3m;
// arms.test.js pins A and A+); the fit truth is not scored. No model runs here.

import { beforeAll, describe, expect, it } from 'vitest';
import { buildRunTools, MISFIT_CONCENTRATION } from '../../src/llm/copilot/runTools';
import { toolResult } from '../../src/llm/copilot/tools';
import { createVerification } from '../../src/llm/copilot/verification';
import { buildRunContext } from '../../src/llm/context/runContext';
import { buildCase, DEV_CASES } from '../cases';
import { CONTEXT_SECTIONS, runLiveCase } from '../harness';
import { loadSkills } from '../run';
import { caseFiles, loadRun, RUN_FILES, RUN_NAME } from '../runFolder';
import { pinRecord } from './recordFiles';

const ITEM = DEV_CASES.find((item) => item.id === 'dev-fit-unconverged');
const MODEL = { name: 'scripted', label: 'scripted test model', model: 'scripted', baseUrl: 'none', temperature: 0 };

// A model that replies from a script: { calls: [[name, args], …] } or { text }.
const scripted = (replies) => {
    let index = 0;
    return async function* call() {
        const reply = replies[Math.min(index, replies.length - 1)];
        index += 1;
        if (reply.text) yield { content: reply.text };
        for (const [position, [name, args]] of (reply.calls || []).entries()) {
            yield { toolCall: { index: position, id: `call_${index}_${position}`, name, arguments: JSON.stringify(args) } };
        }
        yield { finishReason: reply.calls ? 'tool_calls' : 'stop' };
    };
};

describe('dev-fit-unconverged: a run with fit data (plan item 4.6)', () => {
    let built;
    let run;
    let context;
    let tools;

    beforeAll(async () => {
        built = buildCase(ITEM);
        run = await loadRun(built);
        context = buildRunContext(run.props);
        tools = Object.fromEntries(buildRunTools(run.props, { cards: [] }).map((tool) => [tool.name, tool]));
    }, 60000);

    it('is a dev-only case with a seeded configuration and seeded fit files, named as RMCProfile names them', () => {
        expect(ITEM.devOnly).toMatch(/assess_convergence, compare_fit/);
        const again = buildCase(ITEM);
        expect(again.rmc6f).toBe(built.rmc6f);
        expect(again.fitFiles).toEqual(built.fitFiles);
        expect(caseFiles(built).map((file) => file.name)).toEqual([
            `${RUN_NAME}.rmc6f`, `${RUN_NAME}_PDFpartials.csv`, `${RUN_NAME}_FQ1.csv`, `${RUN_NAME}_FT_XFQ1.csv`, `${RUN_NAME}-00.log`
        ]);
        // The fit files are made from the partial g(r); the .rmc6f (with its
        // run counters) does not depend on them.
        const lean = buildCase(ITEM, { partials: false });
        expect(lean.fitFiles).toBeNull();
        expect(lean.rmc6f).toBe(built.rmc6f);
        // No other dev case has fit data.
        DEV_CASES.filter((item) => item !== ITEM).forEach((item) => expect(buildCase(item, { partials: false }).fitFiles, item.id).toBeNull());
    });

    it('opens through the run module as a user\'s run: the kinds detected, the χ² history charted, the counters read', () => {
        expect(run.props.plotFiles.map((file) => [file.name, file.plotKind, file.parseError ?? null])).toEqual([
            [RUN_FILES.log, 'r_value', null],
            [RUN_FILES.fq, 'xray_sq', null],
            [RUN_FILES.gr, 'xpdf', null],
            [RUN_FILES.partials, 'pdf_partials', null]
        ]);
        expect(run.props.rValueFile.name).toBe(RUN_FILES.log);
        expect(run.props.rValueFile.plotData).toMatchObject({ chiColumn: 'X_ray_(R)1' });
        expect(run.props.rValueFile.plotData.series[0].y).toHaveLength(ITEM.params.fit.rows);
        expect(run.props.structure.moves).toMatchObject({ generated: 216000, tried: 209520, accumulatedTimeS: 240000 });
        expect(run.props.symmetry.toleranceA).toBe(0.2);
        expect(run.props.symmetry.spaceGroup).toMatch(/^Fm-3m/);
    });

    it('the context: convergence improving, both datasets with Rw, the FT as the same measurement, the sampling', () => {
        const sections = CONTEXT_SECTIONS.filter((key) => key in context);
        expect(sections).toEqual(['structure', 'configuration_optimization', 'symmetry', 'pair_correlations', 'datasets', 'partial_weights', 'convergence']);
        expect(context.convergence).toMatchObject({ watchdog_status: 'improving', column: 'X_ray_(R)1', n_steps: ITEM.params.fit.rows });
        // Still falling over the last tenth, though the plot's end looks flat.
        // Both χ² terms and their total still fall over the last tenth, each
        // by more than 3 × its noise (the rule of 2026-10-24).
        const { last_10pct: tenth } = context.convergence;
        expect(tenth.rows).toBe(40);
        expect(Object.keys(tenth.terms)).toEqual(['F(Q)_1', 'X_ray_(R)1']);
        [...Object.values(tenth.terms), tenth.total].forEach((entry) => {
            expect(entry.trend).toBe('falling');
            expect(entry.change).toBeLessThan(-0.02);
            expect(entry.change).toBeGreaterThan(-0.1);
            expect(-entry.change).toBeGreaterThan(3 * entry.noise_sd);
        });
        expect(context.convergence.last - context.convergence.first).toBeLessThan(-3);
        expect(context.datasets.map((entry) => entry.kind)).toEqual(['xray_sq', 'xpdf', 'pdf_partials']);
        expect(context.datasets[1]).toMatchObject({ same_measurement_as: 'XFQ1 (its Fourier transform)' });
        context.datasets.slice(0, 2).forEach((entry) => expect(entry.rwp).toBeGreaterThan(0));
        expect(Object.keys(context.partial_weights)).toEqual(['xray', 'note']);
        expect(context.configuration_optimization.accepted_moves_per_atom).toBe(built.truth.fit.accepted_moves_per_atom);
    });

    it('assess_convergence: improving, every χ² term and the total still falling over the last tenth, and its fact for the answer check', async () => {
        const result = await tools.assess_convergence.run({});
        expect(result.data).toMatchObject({ chart_column: 'X_ray_(R)1', status: built.truth.fit.convergence, window_rows: 40 });
        // Every term and the total still fall over the last tenth.
        expect([...result.data.terms.map((term) => [term.term, term.trend]), ['total', result.data.total.trend]])
            .toEqual([['F(Q)_1', 'falling'], ['X_ray_(R)1', 'falling'], ['total', 'falling']]);
        expect(result.facts).toEqual({ convergence: 'improving' });
    });

    it('compare_fit on G(r): the misfit is concentrated in the range holding the planted peak, its largest residual there, data above the model', async () => {
        const { centre, height, sd } = ITEM.params.fit.missed;
        const { data, checks } = await tools.compare_fit.run({ dataset: 'xpdf' });
        expect(data).toMatchObject({ kind: 'X-ray G(r)', axis: 'r', model_series: 'X_ray-calc', data_series: 'X_ray_exp_renorm' });
        const worst = data.by_range.reduce((best, entry) => (entry.share_of_misfit > best.share_of_misfit ? entry : best));
        expect(worst.r_A[0]).toBeLessThan(centre);
        expect(worst.r_A[1]).toBeGreaterThan(centre);
        expect(worst.share_of_misfit).toBeGreaterThanOrEqual(MISFIT_CONCENTRATION);
        const largest = data.largest_residuals.reduce((best, entry) => (Math.abs(entry.model_minus_data) > Math.abs(best.model_minus_data) ? entry : best));
        expect(Math.abs(largest.r - centre)).toBeLessThan(sd);
        expect(largest.model_minus_data).toBeLessThan(0);
        expect(Math.abs(largest.model_minus_data + height)).toBeLessThan(0.2 * height);
        expect(checks[0].text).toMatch(/the misfit is concentrated in r = /);
    });

    it('compare_fit on F(Q), the same measurement: the one missing peak is a misfit at low Q', async () => {
        const { data } = await tools.compare_fit.run({ dataset: 'xray_sq' });
        expect(data).toMatchObject({ kind: 'X-ray F(Q)', axis: 'Q', model_series: 'F(Q)_RMC', data_series: 'F(Q)_Expt' });
        const [lowest, ...rest] = data.by_range;
        expect(lowest.share_of_misfit).toBeGreaterThanOrEqual(MISFIT_CONCENTRATION);
        rest.forEach((entry) => expect(entry.share_of_misfit).toBeLessThan(lowest.share_of_misfit));
    });

    it('the answer check flags "converged", against the tool\'s status and against the context alone', async () => {
        const { facts } = await tools.assess_convergence.run({});
        // The checks through the verification: the run context given, the tool's facts recorded, no claim offered.
        const rules = (text, turnFacts) => {
            const verification = createVerification({ context, offered: [] });
            turnFacts.forEach((item) => verification.record({ name: 'assess_convergence', ...toolResult({ facts: item }) }, () => ({ data: null })));
            return verification.judge(text).done.audit.map((finding) => finding.rule);
        };
        expect(rules('The run has converged: the χ² history has plateaued.', [facts])).toEqual(['convergence']);
        expect(rules('The run has converged.', [])).toEqual(['convergence']);
        expect(rules('The run has not converged: χ² is still falling.', [facts])).toEqual([]);
    });

    it('in the harness: the fit tools and the skill are offered, the record states the fit sections, and "converged" is flagged, not scored', async () => {
        const answer = 'The run has converged and the fit is good. Thermal motion only: no order beyond the average structure.'
            + '\n\nVERDICT {"features": [], "space_group": "Fm-3m", "chiral": false, "centrosymmetric": true, "polar": false}';
        const record = await runLiveCase({
            built, arm: 'B', model: MODEL, skills: loadSkills(),
            // The fit tools, then the analyses a "no order" answer needs (the evidence gate).
            port: scripted([
                { calls: [['assess_convergence', {}], ['compare_fit', { dataset: 'xpdf' }]] },
                { calls: [['analyze_modulations', {}], ['analyze_displacement_directions', { site_ref: 1 }], ['analyze_displacement_directions', { site_ref: 5 }]] },
                { text: answer }
            ])
        });
        expect(record.error).toBeNull();
        expect(record.tools_offered).toEqual(expect.arrayContaining(['read_run_data', 'assess_convergence', 'compare_fit', 'find_gr_peaks']));
        expect(record.skills_offered).toContain('rmc-fit-assessment');
        expect(record.context.sections).toEqual(['structure', 'configuration_optimization', 'symmetry', 'pair_correlations', 'datasets', 'partial_weights', 'convergence']);
        expect(record.context.absent).toEqual(['run_settings', 'pca_displacements']);
        expect(record.steps.map((step) => step.name)).toEqual([
            'assess_convergence', 'compare_fit', 'analyze_modulations', 'analyze_displacement_directions', 'analyze_displacement_directions'
        ]);
        expect(record.steps.every((step) => step.status === 'ok')).toBe(true);
        expect(record.gates).toEqual([]);
        expect(record.final.audit.map((finding) => finding.rule)).toEqual(['convergence']);
        expect(record.final.outcome.status).toBe('not verified');
        // The score is the structure's: right, whatever the answer says of the fit.
        expect(record.score.score).toBe(1);
        expect(record.truth.fit).toMatchObject({ convergence: 'improving', chi_column: 'X_ray_(R)1', misfit: { dataset: 'xpdf', r_A: 2.55 } });
        await pinRecord('fitCase.in-the-harness', record);
    }, 60000);
});
