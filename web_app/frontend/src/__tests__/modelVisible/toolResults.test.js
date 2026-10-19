// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// What the model sees, 4 of 4: one result of every run tool and engine tool
// on the Demo run, as the model reads it — the `tool` message runCopilotTurn
// (llm/copilot/agentLoop.js) appends after executeToolCall: the result fitted
// to its budget, its checks, its reference (`ref`) and any note. One file per
// call in __snapshots__/tool-results/.
//
// Each sample is one turn, set up as useCopilot sets one up (the dock's tools
// without the view tools, which need a screen; the Workbench's own questions;
// the evidence gate; the answer checks): a recording model calls the tool,
// and its second request holds the tool message. The engines run in-thread
// (pcaKdeWorker's handlePcaMessage, as src/openRunFolder.js runs them).
// Also: an ask_user answer, a read_skill result, and a window that cuts the
// Ta–Se shell, once as the app runs it (the Workbench asks, the user skips:
// the call runs in the whole shell, with a note) and once as the tool itself
// refuses it (no Workbench question: no shell catalog, or its questions used).
//
// A change to these files is a change to what the model sees: call it out in
// the PR.

import { beforeAll, describe, expect, it } from 'vitest';
import { runCopilotTurn } from '../../llm/copilot/agentLoop';
import { contextAuditFacts } from '../../llm/copilot/answerAudit';
import { buildAskUserTool } from '../../llm/copilot/askUser';
import { buildClarifier } from '../../llm/copilot/clarify';
import { buildEvidenceGate } from '../../llm/copilot/evidence';
import { buildReadSkillTool, skillIndex } from '../../llm/copilot/skills';
import { buildRunTools } from '../../llm/copilot/runTools';
import { buildRunContext } from '../../llm/context/runContext';
import { pointGroupFacts } from '../../llm/context/pointGroups';
import { buildCopilotMessages, shownContext } from '../../llm/prompts/copilot';
import { copilotSkills } from '../../copilotSkills';
import { buildCopilotCards } from '../../copilotCards';
import { buildCopilotTools } from '../../copilotTools';
import { handlePcaMessage } from '../../workers/pcaKdeWorker';
import { bundledSkills, jsonFile, loadDemoAssistantRun, recordingModel } from './helpers';

const QUESTION = 'Run the analysis.';
const CUT = { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.3, r12_max: 2.7 };

// [file, tool, arguments, options]; options.user answers ask_user (the
// question card's reply), options.workbench false runs without the
// Workbench's own questions.
const SAMPLES = [
    ['read_run_data', 'read_run_data', { section: 'convergence' }],
    ['assess_convergence', 'assess_convergence', {}],
    ['find_gr_peaks', 'find_gr_peaks', { pair: 'Ta-Se' }],
    ['compare_fit', 'compare_fit', { dataset: 'xray_sq' }],
    ['list_sites', 'list_sites', {}],
    ['analyze_displacement_directions', 'analyze_displacement_directions', { site_ref: 1 }],
    ['compute_bond_angles', 'compute_bond_angles', { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.21, r12_max: 3.23 }],
    ['compute_bond_angles.cut-window-workbench-asks', 'compute_bond_angles', CUT, { user: { answer: null, source: 'skipped', index: null } }],
    ['compute_bond_angles.cut-window-refused', 'compute_bond_angles', CUT, { workbench: false }],
    ['analyze_clusters', 'analyze_clusters', { element: 'Ta', r_min: 2.65, r_max: 3.33 }],
    ['analyze_modulations', 'analyze_modulations', {}],
    ['detect_symmetry', 'detect_symmetry', { tolerance_A: 0.2 }],
    ['measure_symmetry_lowering', 'measure_symmetry_lowering', {}],
    ['ask_user', 'ask_user', {
        question: 'Which Ta–Se window?',
        options: [
            { label: 'Whole first shell', description: '2.21–3.23 Å, both Ta–Se peaks' },
            { label: 'First peak only', description: '2.3–2.7 Å' }
        ],
        default: 'Whole first shell'
    }, { user: { answer: 'Whole first shell', source: 'option', index: 0 } }],
    ['read_skill', 'read_skill', { name: 'coordination-shells' }]
];

describe('tool results as the model reads them (the Demo run)', () => {
    let setup;

    beforeAll(async () => {
        const { run, rmc6fText } = await loadDemoAssistantRun();
        const request = (kind, params) => handlePcaMessage({ kind, ...params }, async () => rmc6fText);
        const cards = buildCopilotCards({ plotFiles: run.plotFiles, rValueFile: run.rValueFile, symmetry: run.symmetry });
        const skills = bundledSkills();
        const tools = [
            ...buildRunTools(run, { cards }),
            ...buildCopilotTools({ request, structure: run.structure, elements: Object.keys(run.structure.elementCounts || {}).sort() }),
            buildAskUserTool(),
            buildReadSkillTool(skills, { load: copilotSkills.load })
        ];
        const context = buildRunContext(run);
        setup = {
            run, tools, context,
            index: skillIndex(skills),
            clarifier: buildClarifier(run, { noiseFloor: context.symmetry?.noise_floor_A ?? null })
        };
    });

    it('offers the tools sampled here, and every run and engine tool is sampled', () => {
        const sampled = new Set(SAMPLES.map(([, tool]) => tool));
        expect(setup.tools.map((tool) => tool.name).filter((name) => !sampled.has(name))).toEqual([]);
    });

    it.each(SAMPLES)('%s', async (file, name, args, { user = null, workbench = true } = {}) => {
        const { run, tools, context, index, clarifier } = setup;
        const model = recordingModel([{ calls: [[name, args]] }, { text: 'Recorded.\n\nOutcome: a recorded tool result.' }]);
        const questions = [];
        const askUser = async (asked) => {
            questions.push(asked.question);
            return user || { answer: null, source: 'skipped', index: null };
        };
        const messages = buildCopilotMessages({
            context, history: [], userText: QUESTION, toolsEnabled: true, viewTools: false, askUser: true, skills: index
        });
        const events = [];
        for await (const event of runCopilotTurn({
            callModel: model.callModel,
            messages,
            tools,
            toolContext: { askUser },
            auditContext: contextAuditFacts(context, pointGroupFacts),
            clarifier: workbench ? clarifier : null,
            evidence: buildEvidenceGate({
                catalog: workbench ? clarifier.catalog : null,
                offered: tools.map((tool) => tool.name),
                elements: Object.keys(run.structure.elementCounts || {})
            }),
            referenceContext: shownContext(context)
        })) events.push(event);

        const second = model.requests[1];
        expect(second).toBeDefined();
        const message = second.messages.at(-1);
        expect(message).toMatchObject({ role: 'tool', tool_call_id: 'call_1_0' });
        await expect(jsonFile({
            tool: name,
            arguments: args,
            ...(questions.length ? { workbench_asked: questions, user_reply: user?.source ?? 'skipped' } : {}),
            content_chars: message.content.length,
            content: JSON.parse(message.content)
        })).toMatchFileSnapshot(`__snapshots__/tool-results/${file}.json`);
    });
});
