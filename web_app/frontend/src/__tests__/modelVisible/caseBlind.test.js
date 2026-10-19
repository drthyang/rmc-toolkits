// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// What the model sees, the harness's promise: nothing about the case except
// what the tools compute. A case id names the case's kind ("dev-zone-boundary"),
// so neither the id nor the kind may reach the model. The harness names every
// run and its files RUN_NAME (validation/runFolder.js); before protocol
// 2026-10-19.1 the run was named after the case id, and the context's `run`
// told the model its case kind. No model runs here.
//
//   1. The run folder as the harness loads it (validation/runFolder.js loadRun,
//      through the run module's in-process adapter src/openRunFolder.js): no
//      case id in the CopilotRun's names, paths or view.
//   2. The first request runLiveCase sends on every case, dev and test (the
//      sealed test cases where the loader may open them: CI's
//      RMC_VALIDATION_HOLDOUT=1; skipped locally), arms B and C, captured at
//      its model port: no case id, in any spelling, anywhere
//      in it; no case kind (nor a case id without its split), in any spelling,
//      in its run-specific parts — the context message, the question message
//      past the fixed case question, the request's other fields, and every
//      string of the tools that is not the same in every case. The fixed text
//      (the system prompt, the question with its verdict format, the tool
//      names and descriptions) names kinds of order on purpose: "jahn_teller"
//      is a verdict kind, measure_symmetry_lowering a tool.
//   3. Each dev case's tools, every one called once by a scripted model: no
//      case id in any request of the turn (tool messages, the Workbench's
//      questions, the gate's send-backs) or in any result the tools returned.
//      The dev cases only: the test cases' first requests are covered by 2,
//      and their tool results are not run here.
//
// Only each case's `id` and `kind` are read here: the test cases' parameters
// and truths stay unread (buildCase uses them to build the configuration; a
// failure names the case, the part and the token, never the content).

import { beforeAll, describe, expect, it } from 'vitest';
import { buildCase, DEV_CASES, TEST_CASE_IDS } from '../../../validation/cases';
import { holdoutEnabled, loadTestCases } from '../../../validation/holdoutLoader';
import { liveQuestion, runLiveCase } from '../../../validation/harness';
import { loadRun, RUN_NAME } from '../../../validation/runFolder';
import { loadSkills } from '../../../validation/run';
import { recordingModel } from './helpers';

const MODEL = { name: 'blind', label: 'recording fake model', model: 'blind', baseUrl: 'none', temperature: 0 };
const skills = loadSkills();

// A token in any spelling: its words joined by dashes, underscores, spaces or
// nothing, in any case ("dev-zone-boundary" ≡ "dev_zone_boundary" ≡ "Dev Zone Boundary").
const words = (token) => token.toLowerCase().split(/[^a-z0-9]+/).filter(Boolean);
const spelling = (token) => new RegExp(words(token).join('[-_\\s]*'), 'i');
const unique = (tokens) => [...new Map(tokens.map((token) => [words(token).join(' '), token])).values()];

// Every case the suite may build: the dev set, and the sealed test set where
// the loader may open it. The ids are public, so every id is looked for.
const HOLDOUT = holdoutEnabled();
const CASES = [...DEV_CASES, ...(HOLDOUT ? await loadTestCases() : [])];
const IDS = [...DEV_CASES.map((item) => item.id), ...TEST_CASE_IDS];
const ID_PATTERNS = IDS.map((id) => [id, spelling(id)]);
// The kinds, and each id without its split where that is more than one word
// ("cell-wrap"); a one-word stem is the kind itself or a word every context
// may use ("chiral", as in a point group's facts).
const KIND_TOKENS = unique([
    ...CASES.map((item) => item.kind),
    ...IDS.map((id) => id.replace(/^(dev|test)-/, '')).filter((stem) => words(stem).length > 1)
]);
const KIND_PATTERNS = KIND_TOKENS.map((token) => [token, spelling(token)]);

/** "<case>: <part> names '<token>'" for every pattern found in the text. */
const findings = (caseId, part, text, patterns) => patterns
    .filter(([, pattern]) => pattern.test(text))
    .map(([token]) => `${caseId}: ${part} names '${token}'`);

// Every string leaf of a value, keyed by its path.
const leaves = (value, path = '', out = new Map()) => {
    if (typeof value === 'string') out.set(path, value);
    else if (Array.isArray(value)) value.forEach((item, index) => leaves(item, `${path}[${index}]`, out));
    else if (value && typeof value === 'object') Object.entries(value).forEach(([key, item]) => leaves(item, `${path}.${key}`, out));
    return out;
};
// The tools' strings, keyed by tool name and path (the order sent may differ by case).
const toolLeaves = (tools) => new Map((tools || []).flatMap((tool) => (
    [...leaves(tool).entries()].map(([path, text]) => [`${tool.function?.name}${path}`, text])
)));

const firstRequest = async (built, arm) => {
    const model = recordingModel();
    const record = await runLiveCase({
        built, arm, model: MODEL, port: model.port, skills: arm === 'B' ? skills : null
    });
    expect(record.error).toBeNull();
    return model.requests[0];
};

describe('no case id or kind reaches the model (validation harness)', () => {
    // id → { built, B, C }: every case built once, its first request in each arm.
    const captured = new Map();

    beforeAll(async () => {
        for (const item of CASES) {
            const built = buildCase(item);
            captured.set(item.id, { built, B: await firstRequest(built, 'B'), C: await firstRequest(built, 'C') });
        }
    }, 120000);

    it.skipIf(!HOLDOUT)('captures every held-out case too (CI)', () => {
        expect(TEST_CASE_IDS.every((id) => captured.has(id))).toBe(true);
    });

    it('loads every case as a run named RUN_NAME, with no case id in a name or path', async () => {
        const found = [];
        for (const [id, { built }] of captured) {
            // The run as the Copilot is handed it (the run module's interface).
            const run = await loadRun(built, { ladder: false });
            expect(run.props.runName).toBe(RUN_NAME);
            const names = JSON.stringify({
                runName: run.props.runName,
                structure: run.props.structure.source,
                files: run.props.plotFiles.map((file) => [file.name, file.path]),
                view: run.view
            });
            found.push(...findings(id, 'the run folder', names, ID_PATTERNS));
        }
        expect(found).toEqual([]);
    });

    // Since plan item 4.6 the fit-data case (dev-fit-unconverged) offers the
    // fit tools, so every tool the sweep below calls is offered, and its
    // results checked, on some dev case.
    it('the dev cases between them offer every tool the sweep below calls', () => {
        const offered = (id) => captured.get(id).B.tools.map((tool) => tool.function.name);
        expect([...new Set(DEV_CASES.flatMap((item) => offered(item.id)))].sort()).toEqual(Object.keys(SWEEP).sort());
        expect(offered('dev-fit-unconverged')).toEqual(expect.arrayContaining(['assess_convergence', 'compare_fit']));
    });

    it.each(['B', 'C'])('arm %s: the first request names no case id, anywhere', (arm) => {
        const found = [...captured].flatMap(([id, requests]) => findings(id, 'the request', JSON.stringify(requests[arm]), ID_PATTERNS));
        expect(found).toEqual([]);
    });

    it.each(['B', 'C'])('arm %s: no case kind in the run-specific parts of the first request', (arm) => {
        const all = [...captured.values()];
        // The fixed text: the same in every case.
        const same = (pick) => new Set(all.map((entry) => pick(entry))).size === 1;
        const roles = all.map(({ [arm]: request }) => request.messages.map((message) => message.role).join(','));
        expect(new Set(roles)).toEqual(new Set(['system,user,assistant,user']));
        expect(same(({ [arm]: request }) => request.messages[0].content)).toBe(true);
        expect(same(({ [arm]: request }) => request.messages[2].content)).toBe(true);
        expect(same(({ built }) => liveQuestion(built))).toBe(true);
        const toolStrings = all.map(({ [arm]: request }) => toolLeaves(request.tools));
        const keys = new Set(toolStrings.flatMap((strings) => [...strings.keys()]));
        const fixed = new Set([...keys].filter((key) => (
            toolStrings.every((strings) => strings.has(key)) && new Set(toolStrings.map((strings) => strings.get(key))).size === 1
        )));
        expect(fixed.size).toBeGreaterThan(0);

        const found = [...captured].flatMap(([id, { built, [arm]: request }]) => {
            const { messages, tools, ...fields } = request;
            const question = liveQuestion(built);
            expect(messages[3].content.startsWith(question)).toBe(true);
            const runTools = [...toolLeaves(tools)].filter(([key]) => !fixed.has(key));
            return [
                ...findings(id, 'the context message', messages[1].content, KIND_PATTERNS),
                ...findings(id, 'the question message past the case question', messages[3].content.slice(question.length), KIND_PATTERNS),
                ...findings(id, 'the request fields', JSON.stringify(fields), KIND_PATTERNS),
                ...runTools.flatMap(([key, text]) => findings(id, `the tools (${key})`, text, KIND_PATTERNS))
            ];
        });
        expect(found).toEqual([]);
    });
});

// The calls for each tool the harness offers, from its definition: every
// section, pair, element and skill it lists, a site, a window of 1.5–3.2 Å
// (some calls fail or are refused, and their messages are checked as well).
// A tool missing here fails the sweep: add it, so its results are checked too.
const enumOf = (tool, name) => tool.parameters?.properties?.[name]?.enum || [];
const SWEEP = {
    read_run_data: (tool) => enumOf(tool, 'section').map((section) => ({ section })),
    find_gr_peaks: (tool) => enumOf(tool, 'pair').map((pair) => ({ pair })),
    assess_convergence: () => [{}],
    compare_fit: (tool) => enumOf(tool, 'dataset').map((dataset) => ({ dataset })),
    list_sites: (tool) => [{}, ...enumOf(tool, 'element').map((element) => ({ element }))],
    analyze_displacement_directions: () => [{ site_ref: 1 }, { element: 'all' }],
    compute_bond_angles: (tool) => enumOf(tool, 'apex').map((apex, index, elements) => {
        const end = elements.find((element) => element !== apex) || apex;
        return { end1: end, apex, end2: end, r12_min: 1.5, r12_max: 3.2 };
    }),
    analyze_clusters: (tool) => enumOf(tool, 'element').map((element) => ({ element, r_min: 2, r_max: 4 })),
    analyze_modulations: () => [{}],
    detect_symmetry: () => [{ tolerance_A: 0.05 }],
    measure_symmetry_lowering: () => [{}],
    ask_user: () => [{ question: 'Which window?', options: [{ label: 'The whole first shell' }, { label: 'The first peak only' }] }],
    read_skill: (tool) => enumOf(tool, 'name').map((name) => ({ name }))
};
const CALLS_PER_ROUND = 6;

describe('no case id in the tools\' results (validation harness, dev cases)', () => {
    it.each(DEV_CASES.map((item) => item.id))('%s: every tool once', async (id) => {
        const built = buildCase(DEV_CASES.find((item) => item.id === id));
        // The tools as runLiveCase offers them, read from its first request.
        const probe = recordingModel();
        await runLiveCase({ built, arm: 'B', model: MODEL, port: probe.port, skills });
        const offered = probe.requests[0].tools.map((tool) => tool.function);
        expect(offered.filter((tool) => !SWEEP[tool.name]).map((tool) => tool.name)).toEqual([]);
        const calls = offered.flatMap((tool) => SWEEP[tool.name](tool).map((args) => [tool.name, args]));
        const rounds = [];
        for (let start = 0; start < calls.length; start += CALLS_PER_ROUND) rounds.push({ calls: calls.slice(start, start + CALLS_PER_ROUND) });

        const model = recordingModel([...rounds, { text: 'Recorded.\n\nOutcome: every tool was called once.' }]);
        const record = await runLiveCase({
            built, arm: 'B', model: MODEL, port: model.port, skills,
            // Room for every call; a round whose calls all fail does not end the sweep.
            limits: { maxRounds: rounds.length + 1, maxCallsPerRound: CALLS_PER_ROUND, maxNoProgressRounds: Infinity }
        });
        expect(record.error).toBeNull();
        // Every call ran and its message reached the model.
        const modelSteps = record.steps.filter((step) => step.by === 'model');
        expect(modelSteps.map((step) => step.name)).toEqual(calls.map(([name]) => name));
        const toolMessages = model.requests.at(-1).messages.filter((message) => message.role === 'tool');
        expect(toolMessages.length).toBeGreaterThanOrEqual(calls.length);

        const found = [
            ...model.requests.flatMap((request, index) => findings(id, `request ${index + 1}`, JSON.stringify(request), ID_PATTERNS)),
            ...record.steps.flatMap((step) => findings(id, `step ${step.step} (${step.name})`,
                JSON.stringify({ data: step.data, checks: step.checks, error: step.error }), ID_PATTERNS))
        ];
        expect(found).toEqual([]);
    }, 60000);
});
