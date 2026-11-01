// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The user's report of 2026-11-01, as a regression: the dock showed
//   "has σ²θ = {{r5.data.average_structure.angle_variance_deg2.mean}} deg² and
//    a trans angle of {{r5.data.average_structure.trans_angle_deg.mean}}°,"
// raw. Here the same sentence comes out of a real turn on the Demo run — the
// dock's run tools, engine tools, ask_user, read_skill, the Workbench's own
// questions and the verification, the engines in-thread — with
// polyhedral_distortion (the TaSe6 octahedra in the whole first shell) as the
// fifth result, r5, whose data holds both paths. The model is scripted; no
// server is contacted. The verification fills both (llm copilot/references.js:
// the path written from the tool message's root), and a reference that names
// no value is sent back once with what its path reaches (verification/
// repair.js), never filled with a near miss.

import { beforeAll, describe, expect, it } from 'vitest';
import { runCopilotTurn } from '../llm/copilot/agentLoop';
import { buildAskUserTool } from '../llm/copilot/askUser';
import { buildClarifier } from '../llm/copilot/clarify';
import { buildReadSkillTool, skillIndex } from '../llm/copilot/skills';
import { buildRunTools } from '../llm/copilot/runTools';
import { createVerification } from '../llm/copilot/verification';
import { buildRunContext } from '../llm/context/runContext';
import { buildCopilotMessages } from '../llm/prompts/copilot';
import { copilotSkills } from '../copilotSkills';
import { buildCopilotCards } from '../copilotCards';
import { buildCopilotTools } from '../copilotTools';
import { inProcessEngine } from '../workers/pcaKdeWorker';
import { bundledSkills, loadDemoAssistantRun, recordingModel } from './modelVisible/helpers';

const WINDOW = { r12_min: 2.21, r12_max: 3.23 };
// Five results, in order: r1 … r5 (polyhedral_distortion).
const CALLS = [
    ['read_skill', { name: 'coordination-shells' }],
    ['find_gr_peaks', { pair: 'Ta-Se' }],
    ['list_sites', {}],
    ['compute_bond_angles', { end1: 'Se', apex: 'Ta', end2: 'Se', ...WINDOW }],
    ['polyhedral_distortion', { centre: 'Ta', ligand: 'Se', r_min: 2.21, r_max: 3.23 }]
];
// The user's sentence, verbatim, inside an answer.
const USER_SENTENCE = 'has σ²θ = {{r5.data.average_structure.angle_variance_deg2.mean}} deg² and a trans angle of '
    + '{{r5.data.average_structure.trans_angle_deg.mean}}°,';
const answer = (sentence) => `In 2.21–3.23 Å the average TaSe₆ octahedron ${sentence} against the configuration's `
    + '{{r5.configuration.angle_variance_deg2.mean}} deg² from one configuration.\n\nOutcome: the indices of both, as computed.';

describe('numbers by reference in a real turn on the Demo (the user\'s report, 2026-11-01)', () => {
    let setup;

    beforeAll(async () => {
        const { run, rmc6fText } = await loadDemoAssistantRun();
        const cards = buildCopilotCards({ plotFiles: run.plotFiles, rValueFile: run.rValueFile, symmetry: run.symmetry });
        const skills = bundledSkills();
        const tools = [
            ...buildRunTools(run, { cards }),
            ...buildCopilotTools({
                request: inProcessEngine(rmc6fText), structure: run.structure, elements: Object.keys(run.structure.elementCounts || {}).sort(),
                plotFiles: run.plotFiles, runSettings: run.runSettings
            }),
            buildAskUserTool(),
            buildReadSkillTool(skills, { load: copilotSkills.load })
        ];
        const context = buildRunContext(run);
        setup = { run, tools, context, index: skillIndex(skills), clarifier: buildClarifier(run, { noiseFloor: context.symmetry?.noise_floor_A ?? null }) };
    }, 120000);

    const turn = async (replies) => {
        const { tools, context, index, clarifier } = setup;
        const model = recordingModel([{ calls: CALLS }, ...replies]);
        const events = [];
        for await (const event of runCopilotTurn({
            model: model.port,
            messages: buildCopilotMessages({ context, userText: 'Compare the TaSe6 octahedra with the average structure.', askUser: true, skills: index }),
            tools,
            limits: { maxCallsPerRound: CALLS.length },
            toolContext: { askUser: async () => ({ answer: null, source: 'skipped', index: null }) },
            clarifier,
            verification: createVerification({ context, catalog: clarifier.catalog, offered: tools.map((tool) => tool.name) })
        })) events.push(event);
        return { events, done: events.at(-1), requests: model.requests };
    };

    it('fills the user\'s sentence with the engine\'s values: σ²θ 111 deg², trans angle 159.6°', async () => {
        const { events, done } = await turn([{ text: answer(USER_SENTENCE) }]);
        const steps = events.filter((event) => event.type === 'tool_end');
        expect(steps.map((step) => [step.name, step.ref])).toEqual(CALLS.map(([name], index) => [name, `r${index + 1}`]));
        const r5 = steps[4].data;
        expect(r5.average_structure).toMatchObject({ angle_variance_deg2: { mean: 111 }, trans_angle_deg: { mean: 159.6 } });
        expect(done.text).toBe('In 2.21–3.23 Å the average TaSe₆ octahedron has σ²θ = 111 deg² and a trans angle of 159.6°, '
            + 'against the configuration\'s 127 deg² from one configuration.\n\nOutcome: the indices of both, as computed.');
        expect(done.text).not.toContain('{{');
        expect(done.references).toEqual({ written: 3, filled: 3, unknown: [] });
        expect(done.gates).toEqual([]);
        expect(done.audit.filter((finding) => finding.rule === 'reference')).toEqual([]);
    });

    it('a path that names an object is sent back once with what it holds; the repaired answer is filled', async () => {
        const wrong = answer('has σ²θ = {{r5.data.average_structure.angle_variance_deg2}} deg² and a trans angle of '
            + '{{r5.data.average_structure.trans_angle_deg.mean}}°,');
        const { events, done, requests } = await turn([{ text: wrong }, { text: answer(USER_SENTENCE) }]);
        expect(events.find((event) => event.type === 'gate')).toMatchObject({
            id: 'reference', missing: ['r5.data.average_structure.angle_variance_deg2'], text: wrong
        });
        const nudge = requests.at(-1).messages.at(-1);
        expect(nudge.role).toBe('user');
        expect(nudge.content).toContain('- {{r5.data.average_structure.angle_variance_deg2}} names r5.data.average_structure.angle_variance_deg2, '
            + 'which holds mean, sd: name one value in it.');
        // The app did not pick the mean for the model.
        expect(nudge.content).not.toContain('111');
        expect(done.text).toContain('has σ²θ = 111 deg² and a trans angle of 159.6°,');
        expect(done.references).toEqual({ written: 3, filled: 3, unknown: [] });
        expect(done.gates).toEqual([{ round: 1, id: 'reference', missing: ['r5.data.average_structure.angle_variance_deg2'] }]);
    });
});
