// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The budget test of plan item 5.5 (docs/AGENT_PLAN.md, section E "Size
// against a 32k window"): does the worst realistic Copilot turn fit 32,768
// tokens, the window a user may load the local tier's smallest model with
// (qwen3.8-27b in LM Studio)? The validation's live runs no longer run at
// 32k: their run policy loads every LM Studio model with at least 48k
// (validation/runPolicy.js MIN_CONTEXT_TOKENS, the decision of 2026-10-19),
// which holds both worst turns below.
//
// Counted with the app's own estimator and constants, so the numbers are the
// context guard's (provider/contextGuard.js): estimatePromptTokens — messages
// and tool schemas at CHARS_PER_TOKEN (2.5), tool results and call arguments
// at JSON_CHARS_PER_TOKEN (1.4), MESSAGE_OVERHEAD_TOKENS a message — and
// planContext, which refuses a request whose estimate plus
// REPLY_RESERVE_TOKENS exceeds the model's window.
//
// The worst realistic turn of the dock is the request right after its
// analyses:
//   - the largest system prompt the dock sends (tools + view tools + ask_user
//     + the skills' index), with the dock's full tool set as sent (the JSON of
//     modelVisible/__snapshots__/dock/tools.json);
//   - the context block of the committed Demo run as the dock sends it
//     (modelVisible/__snapshots__/dock/first-request.txt), its acknowledgement
//     and the question;
//   - one round reading the two largest skills whole (read_skill has its own
//     result budget, SKILL_RESULT_BUDGET, and `keep`: a skill stays whole for
//     the question and does not count toward the turn's budget);
//   - one round of analyses at the agent loop's character budget
//     (DEFAULT_LIMITS: maxCallsPerRound results that together fill
//     turnCharBudget, none compacted, as they are the current round's);
//   - the reply reserve.
// The turn runs through the real agent loop (runCopilotTurn) with a scripted
// model, so every message is the loop's own: refs, checks, the skill text as
// read_skill returns it, the calls' echoed arguments. Only the analyses' data
// is filler, sized to the budget: the estimator counts characters, not what
// they say. Earlier questions (conversation memory) are not counted: this is
// a first question.
//
// The same is computed for the validation harness's arm B (skills, no view
// tools) on the dev case with the largest first request — never a held-out
// test case.
//
// TODAY THE WORST TURN DOES NOT FIT. Tool messages (the skills and the
// analyses) are counted at 1.4 characters a token, not section E's 2.5, so two
// skills and a full round of analyses outgrow the window. Plan item 5.5's other
// half, moving long skill material to on-demand references, is skill text: the
// user's decision (.claude/skills is the user's method). This test shrinks
// nothing. It pins today's numbers, so any change to the prompt, the tools,
// the Demo context, the skills or the loop's budgets shows here — update the
// pins in the same PR and quote the new headroom — and it keeps the real
// assertion as an expected failure: once the worst turn fits, `it.fails`
// reports it, and the test becomes a plain `it`.

import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { beforeAll, describe, expect, it } from 'vitest';
import { DEFAULT_LIMITS, runCopilotTurn } from '../llm/copilot/agentLoop';
import { buildReadSkillTool, READ_SKILL, SKILL_RESULT_BUDGET, skillIndex } from '../llm/copilot/skills';
import { buildRunContext } from '../llm/context/runContext';
import { buildCopilotMessages } from '../llm/prompts/copilot';
import {
    estimatePromptTokens, LMSTUDIO_REPLY_MIN_TOKENS, planContext, REPLY_RESERVE_TOKENS
} from '../llm/provider/contextGuard';
import { copilotSkills } from '../copilotSkills';
import { buildCase, DEV_CASES } from '../../validation/cases';
import { runLiveCase } from '../../validation/harness';
import { MODELS } from '../../validation/models';
import { loadSkills } from '../../validation/run';
import { MIN_CONTEXT_TOKENS } from '../../validation/runPolicy';
import { bundledSkills, loadDemoAssistantRun, promptVariants, recordingModel } from './modelVisible/helpers';

// A 32k window: qwen3.8-27b as a user may load it in LM Studio for the dock.
const WINDOW_TOKENS = 32768;

// --- today's numbers (tokens, estimatePromptTokens) --------------------------------
//
// Each term is its part estimated alone, so the terms sum to the whole within
// rounding. firstRequest: the turn's first request; prompt: the request after
// the analyses, estimated whole; needed: prompt + REPLY_RESERVE_TOKENS
// (planContext's `needed`, the Ollama rule); headroom: WINDOW_TOKENS − needed
// (negative: the overflow); lmStudioHeadroom: WINDOW_TOKENS − (prompt +
// LMSTUDIO_REPLY_MIN_TOKENS), the rule the guard applies on LM Studio, which
// serves qwen3.8-27b.

// The two largest skills as read_skill returns them (tool-message tokens).
const TWO_LARGEST_SKILLS_TODAY = Object.freeze([
    { name: 'coordination-shells', chars: 5650, tokens: 4044 },
    { name: 'local-symmetry', chars: 5189, tokens: 3715 }
]);

// Candidate #6 (2026-10-19) added 90 tokens: the Demo's pair_correlations
// read the one shell reading (shells.js; the context block +73) — peaks below
// 8 Å (Ga–Ga's 7.33 Å, Ga–Ta's 6.43 Å) and the first peak alone of the Ta–Se
// and Se–Se shells (first_peak_to_A, first_peak_cn) — and find_gr_peaks'
// description says its shell does not depend on r_max (the tools +16).
// Plan item 4.2 (2026-10-19) added 48: find_gr_peaks' description names a
// joined shell's split (its parts and long − short; the tools +48).
// Plan item 2.5 (2026-10-19) added 8: analyze_displacement_directions'
// description gives the peak direction ± its cell's radius (the tools +9).
// polyhedral_distortion (2026-10-24) added 698: its definition, 1,744
// characters with its method advice (the tools +698).
const DOCK_WORST_TURN_TODAY = Object.freeze({
    systemPrompt: 4762,
    toolJson: 7561,
    contextBlock: 2286,
    question: 70,
    skills: 7814,
    analyses: 8724,
    replyReserve: 8192,
    firstRequest: 14677,
    prompt: 31214,
    needed: 39406,
    headroom: -6638,
    lmStudioHeadroom: -2542
});
// By how much the dock's worst turn overflows the window today.
const DOCK_OVERFLOW_TOKENS_TODAY = 6638;

// Arm B's first request on each dev case is 9.7k–11.2k tokens; the largest,
// since plan item 4.6, is the fit-data case: its fit tools (assess_convergence,
// compare_fit) and its sampling, datasets, weights and convergence sections
// (before it, dev-symmetry-lowering: 10,391 tokens, 35,119 needed, 2,351 over).
const ARM_B_LARGEST_DEV_CASE = 'dev-fit-unconverged';
// find_gr_peaks' description (candidate #6) added 16 tokens to its tools, and
// 49 more with a joined shell's split (plan item 4.2); the directions tool's
// description 8 more with the cell radius (plan item 2.5; the request +9);
// polyhedral_distortion (2026-10-24) 693 more (the request +693).
const ARM_B_WORST_TURN_TODAY = Object.freeze({
    systemPrompt: 4422,
    toolJson: 5338,
    contextBlock: 1574,
    question: 610,
    skills: 7814,
    analyses: 8723,
    replyReserve: 8192,
    firstRequest: 11944,
    prompt: 28480,
    needed: 36672,
    headroom: -3904,
    lmStudioHeadroom: 192
});
// By how much arm B's worst turn (largest dev case) overflows the window today.
const ARM_B_OVERFLOW_TOKENS_TODAY = 3904;

// --- the turn ------------------------------------------------------------------------

const HERE = dirname(fileURLToPath(import.meta.url));
const snapshot = (path) => readFileSync(join(HERE, 'modelVisible/__snapshots__', path), 'utf8');

// The dock's question in the model-visible snapshot (dock.test.jsx).
const DOCK_QUESTION = 'Has this run converged?';

// The round of analyses: the engine tools whose results run longest on the
// Demo (the tool-results snapshots: clusters, modulations, the symmetry
// lowering) and the bond-angle call, whose arguments are the longest, with
// arguments as the model writes them. Elements come from the tool's own
// schema, so the calls validate on any run.
const ANALYSES = [
    ['analyze_clusters', (elements) => ({ element: elements.at(-1), r_min: 2.65, r_max: 3.33 })],
    ['compute_bond_angles', (elements) => ({ end1: elements[0], apex: elements.at(-1), end2: elements[0], r12_min: 2.21, r12_max: 3.23 })],
    ['analyze_modulations', () => ({})],
    ['measure_symmetry_lowering', () => ({})]
];

/** Data whose tool message ({ ref, data }, as the loop sends it) is exactly `chars` long. */
const fillerData = (chars, ref) => ({ filler: 'x'.repeat(chars - JSON.stringify({ ref, data: { filler: '' } }).length) });

/** Run a turn to its end; a failed step fails the test. */
const drain = async (turn) => {
    for await (const event of turn) {
        if (event.type === 'tool_end' && event.error) throw new Error(`${event.name}: ${event.error}`);
    }
};

/** Each skill's read_skill tool message (the loop's own), largest first. */
const skillMessages = async (skills, messages) => {
    const tools = [buildReadSkillTool(skills, { load: copilotSkills.load })];
    const out = [];
    for (const skill of skills) {
        const model = recordingModel([{ calls: [[READ_SKILL, { name: skill.name }]] }, { text: 'Recorded.' }]);
        await drain(runCopilotTurn({ model: model.port, messages, tools }));
        const message = model.requests[1].messages.at(-1);
        out.push({ skill, message, chars: message.content.length, tokens: estimatePromptTokens({ messages: [message] }) });
    }
    return out.sort((a, b) => b.tokens - a.tokens || a.skill.name.localeCompare(b.skill.name));
};

/**
 * The worst realistic turn on top of a first request ({ messages, tools } as
 * sent): the loop reads `readSkills`, then runs one round of analyses at the
 * turn's character budget; the request it sends next is estimated, whole and
 * by part.
 */
const worstTurn = async ({ first, skills, readSkills }) => {
    const { maxCallsPerRound, resultCharBudget, turnCharBudget } = DEFAULT_LIMITS;
    const perCall = turnCharBudget / maxCallsPerRound;
    // One round of maxCallsPerRound results reaches the turn's budget, each
    // within the per-result budget: the loop as it is today.
    expect(Number.isInteger(perCall) && perCall <= resultCharBudget).toBe(true);
    expect(ANALYSES).toHaveLength(maxCallsPerRound);

    const offered = first.tools.map((tool) => tool.function);
    const elements = offered.find((tool) => tool.name === 'compute_bond_angles').parameters.properties.end1.enum;
    const readSkill = buildReadSkillTool(skills, { load: copilotSkills.load });
    let ref = readSkills.length + 1;
    // Every tool as offered; read_skill is the app's own, the analyses answer
    // with filler at the per-call budget.
    const tools = offered.map((tool) => (tool.name === READ_SKILL ? readSkill : {
        name: tool.name,
        description: tool.description,
        parameters: tool.parameters,
        run: async () => {
            const data = fillerData(perCall, `r${ref}`);
            ref += 1;
            return { data, checks: [] };
        }
    }));
    const model = recordingModel([
        { calls: readSkills.map((name) => [READ_SKILL, { name }]) },
        { calls: ANALYSES.map(([name, args]) => [name, args(elements)]) },
        { text: 'Recorded.\n\nOutcome: a recorded turn.' }
    ]);
    await drain(runCopilotTurn({ model: model.port, messages: first.messages, tools }));
    const request = model.requests[2];
    // The tools went out exactly as in the first request.
    expect(request.tools).toEqual(first.tools);
    const opening = first.messages.length;
    const skillRound = request.messages.slice(opening, opening + 1 + readSkills.length);
    const analysisRound = request.messages.slice(opening + 1 + readSkills.length);
    expect(skillRound.map((message) => message.role)).toEqual(['assistant', ...readSkills.map(() => 'tool')]);
    expect(analysisRound.map((message) => message.role)).toEqual(['assistant', ...ANALYSES.map(() => 'tool')]);
    // The analyses fill the turn's budget exactly, uncompacted.
    expect(analysisRound.slice(1).map((message) => message.content.length)).toEqual(ANALYSES.map(() => perCall));

    const prompt = estimatePromptTokens({ messages: request.messages, tools: request.tools });
    const plan = planContext({ estimate: prompt, max: WINDOW_TOKENS });
    const [system, context, ...rest] = first.messages;
    return {
        skillRound,
        plan,
        numbers: {
            systemPrompt: estimatePromptTokens({ messages: [system] }),
            toolJson: estimatePromptTokens({ tools: request.tools }),
            contextBlock: estimatePromptTokens({ messages: [context] }),
            question: estimatePromptTokens({ messages: rest }),
            skills: estimatePromptTokens({ messages: skillRound }),
            analyses: estimatePromptTokens({ messages: analysisRound }),
            replyReserve: REPLY_RESERVE_TOKENS,
            firstRequest: estimatePromptTokens({ messages: first.messages, tools: first.tools }),
            prompt,
            needed: plan.needed,
            headroom: WINDOW_TOKENS - plan.needed,
            lmStudioHeadroom: WINDOW_TOKENS - (prompt + LMSTUDIO_REPLY_MIN_TOKENS)
        }
    };
};

describe('the worst realistic turn in a 32,768-token window (plan 5.5)', () => {
    let dock;
    let skillSizes;
    let dockTurn;
    let armB;
    let armBTurn;

    beforeAll(async () => {
        const skills = bundledSkills();
        const index = skillIndex(skills);
        const { run } = await loadDemoAssistantRun();
        dock = {
            messages: buildCopilotMessages({
                context: buildRunContext(run), history: [], userText: DOCK_QUESTION,
                toolsEnabled: true, viewTools: true, askUser: true, skills: index
            }),
            tools: JSON.parse(snapshot('dock/tools.json')),
            variants: promptVariants(index)
        };
        skillSizes = await skillMessages(skills, dock.messages);
        const readSkills = skillSizes.slice(0, 2).map((entry) => entry.skill.name);
        dockTurn = await worstTurn({ first: dock, skills, readSkills });

        // Arm B's first request on every dev case, as the harness sends it.
        const disk = loadSkills();
        const model = { name: 'budget', label: 'recording fake model', model: 'budget', baseUrl: 'none', temperature: 0.2 };
        armB = [];
        for (const item of DEV_CASES) {
            const recording = recordingModel();
            const record = await runLiveCase({ built: buildCase(item), arm: 'B', model, port: recording.port, skills: disk });
            expect(record.error).toBeNull();
            const [request] = recording.requests;
            armB.push({ id: item.id, request, tokens: estimatePromptTokens({ messages: request.messages, tools: request.tools }) });
        }
        armB.sort((a, b) => b.tokens - a.tokens || a.id.localeCompare(b.id));
        armBTurn = await worstTurn({ first: armB[0].request, skills, readSkills });
    }, 120000);

    it('measures what the dock sends: its largest prompt, its full tool set, the Demo context', () => {
        const [largest] = [...dock.variants.entries()].sort((a, b) => b[1].length - a[1].length);
        expect(largest[0]).toBe('tools-view-ask-skills');
        expect(dock.messages[0]).toEqual({ role: 'system', content: largest[1] });
        expect(snapshot('dock/first-request.txt')).toContain(dock.messages[1].content);
        expect(dock.tools.map((tool) => tool.function.name)).toEqual(expect.arrayContaining([READ_SKILL, 'ask_user', 'orient_3d_view']));
        expect(MODELS['qwen3.8-27b'].server).toBe('lmstudio');
    });

    it('reads the two largest skills whole and keeps them through the analyses', () => {
        expect(skillSizes.slice(0, 2).map(({ skill, chars, tokens }) => ({ name: skill.name, chars, tokens })))
            .toEqual(TWO_LARGEST_SKILLS_TODAY);
        for (const { skill, message } of skillSizes) {
            const { data } = JSON.parse(message.content);
            expect(message.content.length).toBeLessThanOrEqual(SKILL_RESULT_BUDGET);
            expect(data.truncated).toBeUndefined();
            expect(data.text.startsWith(skill.body.slice(0, 200))).toBe(true);
        }
        // Kept whole (same length, same text; only the ref differs) through
        // the round of analyses that followed.
        const kept = dockTurn.skillRound.slice(1);
        const read = skillSizes.slice(0, 2).map((entry) => entry.message);
        expect(kept.map((message) => message.content.length)).toEqual(read.map((message) => message.content.length));
        expect(kept.map((message) => JSON.parse(message.content).data)).toEqual(read.map((message) => JSON.parse(message.content).data));
    });

    it('pins the dock\'s worst turn today', () => {
        expect(dockTurn.numbers).toEqual(DOCK_WORST_TURN_TODAY);
        expect(dockTurn.plan).toEqual({ refuse: true, needed: DOCK_WORST_TURN_TODAY.needed, max: WINDOW_TOKENS });
        expect(-dockTurn.numbers.headroom).toBe(DOCK_OVERFLOW_TOKENS_TODAY);
    });

    // Fails today by DOCK_OVERFLOW_TOKENS_TODAY: 5.5's on-demand skill
    // references wait for the user's skill review. When this starts to pass,
    // `it.fails` reports it: make it a plain `it` and update the pins.
    it.fails('fits the dock\'s worst turn in 32,768 tokens (fails today: 5.5\'s on-demand references are the user\'s decision)', () => {
        expect(dockTurn.numbers.needed).toBeLessThanOrEqual(WINDOW_TOKENS);
        expect(dockTurn.plan.refuse).toBeUndefined();
    });

    it('pins arm B\'s worst turn on its largest dev case today', () => {
        expect(armB.every((entry) => entry.id.startsWith('dev-'))).toBe(true);
        expect(armB[0].id).toBe(ARM_B_LARGEST_DEV_CASE);
        expect(armB[0].request.tools.some((tool) => tool.function.name === READ_SKILL)).toBe(true);
        expect(armBTurn.numbers).toEqual(ARM_B_WORST_TURN_TODAY);
        expect(-armBTurn.numbers.headroom).toBe(ARM_B_OVERFLOW_TOKENS_TODAY);
    });

    it.fails('fits arm B\'s worst turn in 32,768 tokens (fails today: 5.5\'s on-demand references are the user\'s decision)', () => {
        expect(armBTurn.numbers.needed).toBeLessThanOrEqual(WINDOW_TOKENS);
        expect(armBTurn.plan.refuse).toBeUndefined();
    });

    // The run policy's minimum for an LM Studio model (48k) holds both worst
    // turns with the full reply reserve, under both of the guard's rules
    // (planContext's, and LM Studio's prompt + LMSTUDIO_REPLY_MIN_TOKENS):
    // the dock's with 9,746 tokens to spare, arm B's (dev-fit-unconverged since 4.6) with 12,480
    // (polyhedral_distortion took 698 and 693, 2026-10-24).
    it('fits both worst turns in the live runs\' LM Studio minimum (48k, validation/runPolicy.js)', () => {
        const window = MIN_CONTEXT_TOKENS.lmstudio;
        expect(window).toBe(49152);
        for (const turn of [dockTurn, armBTurn]) {
            const plan = planContext({ estimate: turn.numbers.prompt, max: window });
            expect(plan.refuse).toBeUndefined();
            expect(plan.needed).toBe(turn.numbers.needed);
            expect(turn.numbers.prompt + LMSTUDIO_REPLY_MIN_TOKENS).toBeLessThanOrEqual(window);
        }
        expect([window - dockTurn.numbers.needed, window - armBTurn.numbers.needed]).toEqual([9746, 12480]);
    });
});
