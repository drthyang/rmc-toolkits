// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// What the model sees, 2-3 of 4, the validation harness's path: the first
// request runLiveCase (validation/harness.js) sends on a dev case, captured at
// its model port (`port`: a recording port here, which records the turn's
// intent — messages, tools, tool choice — with its signal; validation/run.js
// gives providerModel, which adds the model's base URL, id and temperature
// and sends through streamCopilotChat), with the simulated user that skips
// every question. No model runs.
//
//   harness/first-request.arm-B.txt   arm B (skills): fields, tool names, messages
//   harness/first-request.arm-C.txt   arm C (no method), against arm B
//   harness/tools.arm-B.json          the `tools` field exactly as sent (arm B)
//   harness/tools.arm-C.json          the `tools` field exactly as sent (arm C)
//
// Arm C follows no method (the turn's method 'off', src/llm/copilot/method.js,
// plan item 1.8, protocol 2026-10-19.10): its first request is also checked
// for every skill name, the skills' step headings, the prompt's method
// sentences, EVIDENCE_RULES, the domain rules' inference parts and the tools'
// method advice.
//
// A change to these files is a change to what the model sees: call it out in
// the PR. The files' header line still names the capture point by its name
// before the turn module (runLiveCase's callModel): PORT below is kept as it
// was so the files stay byte-identical.

import { describe, expect, it } from 'vitest';
import { buildCase, DEV_CASES } from '../../../validation/cases';
import { runLiveCase } from '../../../validation/harness';
import { loadSkills } from '../../../validation/run';
import { loadRun } from '../../../validation/runFolder';
import { buildRunTools } from '../../llm/copilot/runTools';
import { parseSkillFiles, skillIndex } from '../../llm/copilot/skills';
import { EVIDENCE_RULES, METHOD_SENTENCES } from '../../llm/prompts/copilot';
import { GROUND_RULES } from '../../llm/prompts/system';
import { bundledSkills, jsonFile, promptVariants, recordingModel, renderRequest } from './helpers';

const CASE = 'dev-offcentre';
const MODEL = { name: 'snapshot', label: 'recording fake model', model: 'snapshot', baseUrl: 'none', temperature: 0.2 };
const variants = promptVariants(skillIndex(bundledSkills()));
const PORT = 'runLiveCase\'s callModel(request) (validation/harness.js), which validation/run.js binds to '
    + 'streamCopilotChat({ baseUrl, model, temperature, apiKey, session, ...request })';

const firstRequest = async (arm) => {
    const model = recordingModel();
    const record = await runLiveCase({
        built: buildCase(DEV_CASES.find((item) => item.id === CASE)), arm, model: MODEL,
        port: model.port, skills: arm === 'B' ? loadSkills() : null
    });
    expect(record.error).toBeNull();
    return { request: model.requests[0], record };
};

describe(`the harness's first request (${CASE})`, () => {
    it('arm B: the request and its tools', async () => {
        const { request, record } = await firstRequest('B');
        expect(request.tools.map((tool) => tool.function.name)).toEqual(record.tools_offered);
        await expect(jsonFile(request.tools)).toMatchFileSnapshot('__snapshots__/harness/tools.arm-B.json');
        await expect(renderRequest(request, {
            title: `First request of the validation harness, arm B (skills), dev case ${CASE}`,
            port: PORT, variants, toolsFile: 'harness/tools.arm-B.json'
        })).toMatchFileSnapshot('__snapshots__/harness/first-request.arm-B.txt');
    });

    it('arm C: the same request without the method', async () => {
        const [{ request: b }, { request: c, record }] = [await firstRequest('B'), await firstRequest('C')];
        expect(record.method).toBe('off');
        // Arm C's tools are arm B's without read_skill, with the same names and
        // parameters; a description loses its method advice and nothing else.
        const shared = b.tools.filter((tool) => tool.function.name !== 'read_skill');
        expect(c.tools.map((tool) => tool.function.name)).toEqual(shared.map((tool) => tool.function.name));
        expect(c.tools.map((tool) => tool.function.parameters)).toEqual(shared.map((tool) => tool.function.parameters));
        await expect(jsonFile(c.tools)).toMatchFileSnapshot('__snapshots__/harness/tools.arm-C.json');
        await expect(renderRequest(c, {
            title: `First request of the validation harness, arm C (no method), dev case ${CASE}`,
            port: PORT, variants, toolsFile: 'harness/tools.arm-C.json',
            base: b, baseName: 'first-request.arm-B.txt'
        })).toMatchFileSnapshot('__snapshots__/harness/first-request.arm-C.txt');
    });

    it('arm C names no skill or skill step and sends no method sentence, EVIDENCE_RULES or tool advice', async () => {
        const [{ request: b }, { request: c }] = [await firstRequest('B'), await firstRequest('C')];
        // Every text arm C's first request holds: the messages, and each tool's
        // description and parameter descriptions.
        const descriptions = (node) => (node && typeof node === 'object'
            ? Object.entries(node).flatMap(([key, value]) => (key === 'description' && typeof value === 'string' ? [value] : descriptions(value)))
            : []);
        const texts = (request) => [...request.messages.map((message) => message.content), ...descriptions(request.tools)];
        const sent = texts(c).join('\n');
        // Step headings that are plain domain terms the tool descriptions use
        // for what a tool computes, not a method stage: the local-order
        // draft's "Average structure" step (#71) is analyze_clusters',
        // analyze_modulations' and detect_symmetry's own words.
        const DOMAIN_TERMS = ['average structure'];
        const disk = loadSkills();
        const skills = parseSkillFiles(disk.files, { available: disk.references });
        // The skills' names, and each skill's step headings of more than one
        // word ("Find the shell", "Read the significances"); the one-word ones
        // ("Convergence", "Directions") are words any tool may need.
        const stages = skills.flatMap((skill) => [...(/^## Steps\s*\n([\s\S]*?)(?=^## |(?![\s\S]))/m.exec(skill.body)?.[1] || '')
            .matchAll(/^\d+\.\s+\*\*([^*]+?)[.:]?\*\*/gm)].map((match) => match[1].trim()).filter((stage) => stage.split(/\s+/).length > 1 && !DOMAIN_TERMS.includes(stage.toLowerCase())));
        expect(skills.map((skill) => skill.name)).toEqual(['coordination-shells', 'displacement-anisotropy', 'local-order', 'local-symmetry', 'rmc-fit-assessment']);
        expect(stages.length).toBeGreaterThan(15);
        for (const name of skills.map((skill) => skill.name)) expect(sent).not.toContain(name);
        for (const stage of stages) expect(sent.toLowerCase()).not.toContain(stage.toLowerCase());
        const system = c.messages[0].content;
        expect(system).not.toContain(EVIDENCE_RULES);
        expect(system).not.toContain('Skills —');
        for (const sentence of METHOD_SENTENCES) expect(sent).not.toContain(sentence);
        // The domain rules: every definition and conduct part, no inference part (system.js GROUND_RULES).
        for (const part of GROUND_RULES.flat()) {
            if (part.kind === 'inference') expect(sent).not.toContain(part.text.trim());
            else expect(system).toContain(part.text);
        }
        // Every piece of method advice the tools carry: in arm B's request,
        // in none of arm C's descriptions.
        const run = await loadRun(buildCase(DEV_CASES.find((item) => item.id === CASE)));
        const advice = [...buildRunTools(run.props, { cards: [] }), ...run.tools]
            .flatMap((tool) => [].concat(tool.method ?? []).map((part) => (typeof part === 'string' ? part : part.text).trim()));
        // find_gr_peaks, list_sites, compute_bond_angles, analyze_clusters (2),
        // analyze_modulations, detect_symmetry, measure_symmetry_lowering.
        expect(advice).toHaveLength(8);
        const armB = descriptions(b.tools).join('\n');
        for (const piece of advice) {
            expect(armB).toContain(piece);
            expect(sent).not.toContain(piece);
        }
        // B against C: the rest of the request is the same.
        expect(c.messages.slice(1)).toEqual(b.messages.slice(1));
    });
});

// The evidence gate's send-back (THY-18), as the model receives it: the first
// answer concludes "no order" with no analysis run, and the next request ends
// with the gate's message. Under arm C's method 'off' the gate checks the same
// rule but its words name no analysis and no tool.
//
//   harness/gate-send-back.arm-B.txt   the message, naming the analyses missing
//   harness/gate-send-back.arm-C.txt   the message, naming none
describe(`the evidence gate's send-back in the harness (${CASE})`, () => {
    const sendBack = async (arm) => {
        const model = recordingModel([
            { text: 'There is no local order beyond thermal motion.\n\nOutcome: thermal motion only.' },
            { text: 'Recorded.\n\nOutcome: a recorded request.' }
        ]);
        const record = await runLiveCase({
            built: buildCase(DEV_CASES.find((item) => item.id === CASE)), arm, model: MODEL,
            port: model.port, skills: arm === 'B' ? loadSkills() : null
        });
        expect(record.error).toBeNull();
        expect(record.requests.map((request) => request.why)).toEqual(['first', 'gate']);
        return { record, request: model.requests[1] };
    };
    const render = ({ record, request }, arm, label) => {
        const [gate] = record.gates;
        const message = request.messages.at(-1);
        return [
            `# The evidence gate's send-back in the validation harness, arm ${arm} (${label}), dev case ${CASE}`,
            '# The first answer concludes "no order beyond thermal motion" with no analysis run: the gate sends it back.',
            '',
            '## The gate, as recorded (the model never reads the record)',
            `id: ${gate.id}`,
            `missing: ${JSON.stringify(gate.missing)}`,
            '',
            `## The request it makes (why: gate), its last message · messages[${request.messages.length - 1}] · ${message.role}`,
            message.content,
            ''
        ].join('\n');
    };

    it('arm B names the analyses missing; arm C, the same gate, names none', async () => {
        const [b, c] = [await sendBack('B'), await sendBack('C')];
        // The same rule and the same missing analyses, for the record.
        expect(c.record.gates.map(({ id, missing }) => ({ id, missing }))).toEqual(b.record.gates.map(({ id, missing }) => ({ id, missing })));
        expect(c.record.gates[0].id).toBe('no-order');
        const tools = /\b(?:find_gr_peaks|list_sites|analyze_[a-z_]+|compute_bond_angles|detect_symmetry|measure_symmetry_lowering|compare_fit|read_run_data|assess_convergence|ask_user|read_skill)\b/;
        expect(b.request.messages.at(-1).content).toMatch(tools);
        expect(c.request.messages.at(-1).content).not.toMatch(tools);
        // Everything before the send-back is the arm's own first request and the withdrawn answer.
        expect(c.request.messages.at(-2)).toEqual(b.request.messages.at(-2));
        await expect(render(b, 'B', 'skills')).toMatchFileSnapshot('__snapshots__/harness/gate-send-back.arm-B.txt');
        await expect(render(c, 'C', 'no method')).toMatchFileSnapshot('__snapshots__/harness/gate-send-back.arm-C.txt');
    });
});
