// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Every live tool's `data` against the shape its spec declares (`dataShape`,
// llm/copilot/dataShape.js). With numbers by reference the model writes
// {{r3.sites.0.shift_A}}, so the paths of a result are an interface: a key
// renamed, dropped or added must show here, and the PR that changes a tool's
// output changes its dataShape in the same diff.
//
// The tools run as the Copilot runs them (executeToolCall: arguments parsed
// and validated as the model's are), the engines in-thread (pcaKdeWorker's
// inProcessEngine), on the Demo run (as the dock gets it, modelVisible/
// helpers.js) and on every dev case of the validation (validation/runFolder.js,
// as the harness builds a run), with arguments derived from each run's own
// data: every read_run_data section, every partial g(r) pair, the whole first
// shell of each mixed pair (and a cut, a confirmed cut, two windows and an
// empty window), the clusters of each same-element shell, and so on. The fit
// tools (assess_convergence, compare_fit) run on the Demo and on the fit-data
// dev case (dev-fit-unconverged, plan item 4.6). The held-out test cases are
// never built here.
//
// Each result must fit its shape (no undeclared key, no required key
// missing, every value of its declared type), and every key a shape declares
// must be returned by some call, except the forms listed in NOT_REACHED that
// these runs cannot produce.

import { beforeAll, describe, expect, it } from 'vitest';
import { checkShape, declaredKeys } from '../llm/copilot/dataShape';
import { buildAskUserTool } from '../llm/copilot/askUser';
import { buildReadSkillTool } from '../llm/copilot/skills';
import { buildRunTools } from '../llm/copilot/runTools';
import { executeToolCall, toRequestTools } from '../llm/copilot/tools';
import { copilotSkills } from '../copilotSkills';
import { buildCopilotTools } from '../copilotTools';
import { inProcessEngine } from '../workers/pcaKdeWorker';
import { buildCase, DEV_CASES } from '../../validation/cases';
import { loadRun } from '../../validation/runFolder';
import { bundledSkills, loadDemoAssistantRun } from './modelVisible/helpers';

// Declared keys no call here returns, with the reason. Shrink it when a run
// reaches one; a key that turns up here unasked is stale or unexercised.
const NOT_REACHED = {
    assess_convergence: {
        non_finite_steps: 'the χ² histories here (the Demo\'s, dev-fit-unconverged\'s) are finite'
    },
    read_run_data: {
        'symmetry.note': 'only a structure the finder did not analyse',
        'symmetry.sites_omitted': 'no run has more than 12 orbits',
        'pca_displacements.sites[].degenerate': 'no degenerate site',
        'pca_displacements.sites[].zero_spread': 'no zero-spread site',
        'pca_displacements.sites[].mixed': 'no mixed-occupancy site',
        'pca_displacements.sites[].element_counts': 'no mixed-occupancy site',
        'pca_displacements.sites[].element_counts.*': 'no mixed-occupancy site',
        'convergence.non_finite_steps': 'the χ² histories here (the Demo\'s, dev-fit-unconverged\'s) are finite',
        'convergence.non_finite_note': 'the χ² histories here (the Demo\'s, dev-fit-unconverged\'s) are finite'
    },
    find_gr_peaks: {},
    // A G(r) with σ(r): no run here has an F(Q) with an error column (the
    // Demo's input file is not committed, dev-fit-unconverged has no .dat),
    // so each reads { sigma: null, why }; runTools.test.js runs the σ form.
    compare_fit: Object.fromEntries([
        'largest_residuals[].above_noise', 'largest_residuals[].relative_to_sigma', 'noise.above_noise_at',
        'noise.correlation', 'noise.correlation.effective_points', 'noise.correlation.neighbour', 'noise.correlation.points',
        'noise.from', 'noise.reproduction', 'noise.reproduction.residual', 'noise.reproduction.tolerance',
        'noise.rms_residual_over_sigma'
    ].map((key) => [key, 'no run here carries σ(r): no F(Q) with an error column'])),
    list_sites: {
        'sites[].zero_spread': 'no zero-spread site'
    },
    analyze_displacement_directions: {},
    compute_bond_angles: {},
    analyze_clusters: {
        'clusters.other_sizes': 'no window here gives clusters of two sizes',
        'clusters.other_sizes[].count': 'no window here gives clusters of two sizes',
        'clusters.other_sizes[].members': 'no window here gives clusters of two sizes'
    },
    analyze_modulations: {
        untestable_zone_boundary: 'every supercell axis here is even'
    },
    detect_symmetry: {},
    measure_symmetry_lowering: {},
    ask_user: {},
    read_skill: {}
};

// The tools on a run, as the harness offers them (no view tools), with the
// PCA table the PCA Ellipsoid page would hand the dock (read_run_data's
// pca_displacements section).
const runTools = async ({ props, request }) => {
    const pcaSites = await request('sites', {});
    const elements = Object.keys(props.structure.elementCounts || {}).sort();
    return [
        ...buildRunTools({ ...props, pcaSites }, { cards: [] }),
        ...buildCopilotTools({ request, structure: props.structure, elements, plotFiles: props.plotFiles, runSettings: props.runSettings })
    ];
};

const finiteWindow = (window) => Array.isArray(window) && window.every(Number.isFinite) && window[0] < window[1];

// Every call on one run, the arguments read from the run's own data.
const exercise = async (tools, structure) => {
    const byName = new Map(tools.map((tool) => [tool.name, tool]));
    const enumOf = (name, property) => byName.get(name)?.parameters?.properties?.[property]?.enum || [];
    const results = [];
    const call = async (name, args) => {
        const tool = byName.get(name);
        if (!tool) return null;
        const result = await executeToolCall(tool, JSON.stringify(args), { timeoutMs: 120000 });
        results.push({ tool, args, ...result });
        return result;
    };
    const elements = Object.keys(structure.elementCounts || {}).sort();
    const counts = structure.elementCounts || {};

    for (const section of enumOf('read_run_data', 'section')) await call('read_run_data', { section });
    await call('assess_convergence', {});

    const shells = [];
    for (const pair of enumOf('find_gr_peaks', 'pair')) {
        const result = await call('find_gr_peaks', { pair });
        const [a, b] = pair.split('-').map((element) => element.trim());
        if (finiteWindow(result?.data?.first_shell?.window_A)) shells.push({ a, b, window: result.data.first_shell.window_A });
    }
    const pairs = enumOf('find_gr_peaks', 'pair');
    if (pairs.length) await call('find_gr_peaks', { pair: pairs[0], r_max: 1, max_peaks: 1 });

    const datasets = enumOf('compare_fit', 'dataset');
    for (const dataset of datasets) await call('compare_fit', { dataset });
    if (datasets.length) await call('compare_fit', { dataset: datasets[0], n_ranges: 12 });

    const sites = await call('list_sites', {});
    if (elements.length) await call('list_sites', { element: elements[0] });
    const refs = (sites?.data?.sites || []).map((site) => site.ref);
    if (refs.length) {
        await call('analyze_displacement_directions', { site_ref: refs[0] });
        await call('analyze_displacement_directions', { site_ref: refs.at(-1), frequency: 4 });
    }
    if (elements.length) await call('analyze_displacement_directions', { element: elements[0] });
    await call('analyze_displacement_directions', { element: 'all' });

    // Bond angles around the less abundant element of each mixed pair.
    const mixed = shells.filter(({ a, b }) => a !== b).map(({ a, b, window }) => {
        const [apex, end] = (counts[a] || 0) <= (counts[b] || 0) ? [a, b] : [b, a];
        return { apex, end, window };
    });
    for (const { apex, end, window } of mixed) {
        await call('compute_bond_angles', { end1: end, apex, end2: end, r12_min: window[0], r12_max: window[1] });
    }
    if (mixed.length) {
        const { apex, end, window: [lo, hi] } = mixed[0];
        const triplet = { end1: end, apex, end2: end };
        const half = Number(((lo + hi) / 2).toFixed(3));
        // Two windows (an equal pair is read as one shared window).
        await call('compute_bond_angles', { ...triplet, r12_min: lo, r12_max: hi, r23_min: lo, r23_max: Number((hi + 0.05).toFixed(3)), bin_width: 5 });
        await call('compute_bond_angles', { ...triplet, r12_min: lo, r12_max: half });
        await call('compute_bond_angles', { ...triplet, r12_min: lo, r12_max: half, window_confirmed: true });
        await call('compute_bond_angles', { ...triplet, r12_min: 0.1, r12_max: 0.2 });
    }

    // Clusters of each element whose own first shell closes.
    for (const { a, b, window } of shells) {
        if (a === b && window[1] >= 0.5 && window[1] <= 8) await call('analyze_clusters', { element: a, r_min: window[0], r_max: window[1] });
    }
    if (elements.length) await call('analyze_clusters', { element: elements[0], r_min: 0.1, r_max: 0.6 });

    await call('analyze_modulations', {});
    if (elements.length) await call('analyze_modulations', { element: elements.at(-1) });

    for (const tolerance of [0.2, 0.001, 1]) await call('detect_symmetry', { tolerance_A: tolerance });
    await call('measure_symmetry_lowering', {});
    await call('measure_symmetry_lowering', { parent_tolerance_A: 0.2, tolerance_A: 0.001 });
    return results;
};

// ask_user and read_skill do not read the run: each reply form once, every
// skill and every reference it lists once.
const exerciseConversationTools = async () => {
    const results = [];
    const ask = buildAskUserTool();
    const question = {
        question: 'Which Ta–Se window?',
        options: [{ label: 'Whole first shell', description: '2.21–3.23 Å' }, { label: 'First peak only' }]
    };
    const replies = [
        [question, { answer: 'Whole first shell', source: 'option', index: 0 }],
        [question, { answer: 'First peak only', source: 'option', index: 1 }],
        [question, { answer: '2.2–3.0 Å', source: 'other', index: null }],
        [{ ...question, default: 'Whole first shell' }, { answer: null, source: 'skipped', index: null }],
        [question, { answer: null, source: 'skipped', index: null }]
    ];
    for (const [args, reply] of replies) {
        const result = await executeToolCall(ask, JSON.stringify(args), { context: { askUser: async () => reply } });
        results.push({ tool: ask, args, ...result });
    }
    const skills = bundledSkills();
    const read = buildReadSkillTool(skills, { load: copilotSkills.load });
    for (const skill of skills) {
        for (const args of [{ name: skill.name }, ...skill.references.map((reference) => ({ name: skill.name, reference: reference.name }))]) {
            const result = await executeToolCall(read, JSON.stringify(args), { timeoutMs: 60000 });
            results.push({ tool: read, args, ...result });
        }
    }
    return results;
};

const RUNS = ['demo', ...DEV_CASES.map((item) => item.id)];

describe('every live tool\'s data fits its declared dataShape', () => {
    const results = new Map();

    beforeAll(async () => {
        const demo = await loadDemoAssistantRun();
        const runs = [{
            name: 'demo',
            props: demo.run,
            request: inProcessEngine(demo.rmc6fText)
        }];
        // The dev cases only: the held-out test cases are never built here.
        for (const item of DEV_CASES) {
            // The harness's run (validation/runFolder.js → the run module).
            const run = await loadRun(buildCase(item));
            runs.push({ name: item.id, props: run.props, request: run.engine });
        }
        for (const run of runs) {
            const tools = await runTools(run);
            results.set(run.name, { tools, calls: await exercise(tools, run.props.structure) });
        }
        const conversation = await exerciseConversationTools();
        results.set('conversation', { tools: [...new Set(conversation.map((entry) => entry.tool))], calls: conversation });
    }, 120000);

    it.each([...RUNS, 'conversation'])('%s: every result fits its tool\'s dataShape', (name) => {
        const { calls } = results.get(name);
        expect(calls.length).toBeGreaterThan(0);
        const problems = calls.flatMap(({ tool, args, data, error }) => {
            const call = `${tool.name}(${JSON.stringify(args)})`;
            if (error) return [`${call} failed: ${error}`];
            return checkShape(data, tool.dataShape).map((problem) => `${call} ${problem}`);
        });
        expect(problems).toEqual([]);
    });

    it('every tool offered declares a dataShape and is called on some run', () => {
        const offered = new Set();
        const called = new Set();
        for (const { tools, calls } of results.values()) {
            tools.forEach((tool) => offered.add(tool));
            calls.forEach(({ tool }) => called.add(tool.name));
        }
        const names = [...new Set([...offered].map((tool) => tool.name))].sort();
        expect(names).toEqual(Object.keys(NOT_REACHED).sort());
        expect([...offered].filter((tool) => !tool.dataShape).map((tool) => tool.name)).toEqual([]);
        expect(names.filter((name) => !called.has(name))).toEqual([]);
    });

    it('every key a dataShape declares is returned by some call, except the forms NOT_REACHED lists', () => {
        const declared = new Map();
        const seen = new Map();
        for (const { tools, calls } of results.values()) {
            tools.forEach((tool) => {
                const keys = declared.get(tool.name) || new Set();
                declaredKeys(tool.dataShape).forEach((key) => keys.add(key));
                declared.set(tool.name, keys);
            });
            calls.forEach(({ tool, data, error }) => {
                if (error) return;
                const keys = seen.get(tool.name) || new Set();
                checkShape(data, tool.dataShape, { seen: keys });
                seen.set(tool.name, keys);
            });
        }
        const unseen = Object.fromEntries([...declared].map(([name, keys]) => [
            name, [...keys].filter((key) => !seen.get(name)?.has(key)).sort()
        ]));
        const expected = Object.fromEntries(Object.entries(NOT_REACHED).map(([name, keys]) => [name, Object.keys(keys).sort()]));
        expect(unseen).toEqual(expected);
    });

    it('dev-fit-unconverged (plan item 4.6): both fit tools run on its own fit files, every dataset', () => {
        const { tools, calls } = results.get('dev-fit-unconverged');
        const fitTools = ['assess_convergence', 'compare_fit'];
        expect(tools.filter((tool) => fitTools.includes(tool.name)).map((tool) => tool.name)).toEqual(fitTools);
        const fitCalls = calls.filter(({ tool }) => fitTools.includes(tool.name));
        expect(fitCalls.map(({ tool, args }) => `${tool.name}(${JSON.stringify(args)})`)).toEqual([
            'assess_convergence({})',
            'compare_fit({"dataset":"xray_sq"})',
            'compare_fit({"dataset":"xpdf"})',
            'compare_fit({"dataset":"xray_sq","n_ranges":12})'
        ]);
        fitCalls.forEach(({ error, data }) => {
            expect(error ?? null).toBeNull();
            expect(data).toBeTruthy();
        });
        // The other dev cases have no fit data, so no fit tool.
        DEV_CASES.filter((item) => item.id !== 'dev-fit-unconverged').forEach((item) => {
            expect(results.get(item.id).tools.some((tool) => fitTools.includes(tool.name)), item.id).toBe(false);
        });
    });

    it('the model never receives a dataShape: a request tool is its name, description and parameters', () => {
        const { tools } = results.get('demo');
        const sent = toRequestTools([...tools, ...results.get('conversation').tools]);
        sent.forEach((entry) => expect(Object.keys(entry.function).sort()).toEqual(['description', 'name', 'parameters']));
        expect(JSON.stringify(sent)).not.toMatch(/dataShape|number\|null/);
    });
});
