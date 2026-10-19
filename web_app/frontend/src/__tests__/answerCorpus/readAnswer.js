// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The answer corpus's runner (docs/AGENT_PLAN.md item 3.5; the corpus is the
// entries/ folder beside this file, one JSON file per real answer): the run
// an answer was given on (the Demo, or a dev case of the validation) as the
// Copilot is handed it, an entry's recorded steps run as the agent loop runs
// them, and ONE adapter from an answer to what the current checks read in it.
//
// `readAnswer` is the only place that knows the verification
// (llm/copilot/verification/), through its interface only: record the steps,
// then judge. It was rewritten when that one module replaced the evidence
// gate, the answer checks and the verdict (architecture candidate #4); the
// corpus pins stayed as they were.

import { copilotSkills } from '../../copilotSkills';
import { buildAskUserTool } from '../../llm/copilot/askUser';
import { buildClarifier } from '../../llm/copilot/clarify';
import { buildRunTools } from '../../llm/copilot/runTools';
import { buildReadSkillTool } from '../../llm/copilot/skills';
import { executeToolCall, stepStatus } from '../../llm/copilot/tools';
import { createVerification } from '../../llm/copilot/verification';
import { buildRunContext } from '../../llm/context/runContext';
import { buildCase, DEV_CASES } from '../../../validation/cases';
import { loadRun } from '../../../validation/runFolder';
import { bundledSkills, loadDemoAssistantRun } from '../modelVisible/helpers';

/** What an entry pins for the gate when its transcript does not show which tools ran. */
export const STEPS_UNKNOWN = 'steps unknown';

// A verification as the turn builds one (turn.js turnChecks: the run
// context, the clarifier's shell catalog, the offered tools), with `steps`
// recorded as the agent loop records them.
const verificationOf = (run, steps) => {
    const verification = createVerification({ context: run.context, catalog: run.catalog, offered: run.offered });
    steps.forEach((step) => verification.record(step, () => ({ data: null })));
    return verification;
};

/**
 * The adapter: an answer on `run` (openCorpusRun) → what the current checks
 * read in it, through the verification's interface.
 *   claims         the claims the answer makes (in a sentence about the run):
 *                  judged with nothing recorded, each is a missing-evidence
 *                  finding;
 *   gateSendBack   the claim the gate would send the answer back for, given
 *                  the steps that ran (`steps`, runSteps: judged while tools
 *                  can still run): { id, missing } or null; STEPS_UNKNOWN
 *                  when they are not known (`steps` null);
 *   auditFindings  the checks on the answer that fire (their rule names),
 *                  against the steps' facts and the run context's.
 */
export const readAnswer = (text, { run, steps = null }) => {
    const evidence = (finding) => finding.rule === 'evidence';
    const sendBack = steps ? verificationOf(run, steps).judge(text, { canSendBack: true }).sendBack : null;
    return {
        claims: verificationOf(run, []).judge(text).done.audit.filter(evidence).map((finding) => finding.id),
        gateSendBack: steps ? (sendBack ? { id: sendBack.id, missing: sendBack.missing } : null) : STEPS_UNKNOWN,
        auditFindings: verificationOf(run, steps || []).judge(text).done.audit
            .filter((finding) => !evidence(finding) && finding.rule !== 'reference')
            .map((finding) => finding.rule)
    };
};

/** The run an entry was answered on: 'Demo', or a dev case's id (validation/cases.js). */
export const runKeyOf = (entry) => (/^Demo\b/.test(entry.run) ? 'Demo' : entry.run);

// A run as the Copilot sees it (a CopilotRun): the Demo as App.jsx hands it
// to CopilotChat, or a dev case as the validation harness opens it (never a
// held-out case).
const openCopilotRun = async (key) => {
    if (key === 'Demo') return (await loadDemoAssistantRun()).copilot;
    const item = DEV_CASES.find((candidate) => candidate.id === key);
    if (!item) throw new Error(`the corpus reads answers on the Demo or a dev case, not "${key}"`);
    return loadRun(buildCase(item));
};

/**
 * The run (`key`: 'Demo' or a dev case's id) at the default view, with what
 * the dock builds for a turn (useCopilot.js, CopilotChat.jsx), as
 * `{ tools, offered, catalog, context }`: the run tools, the engine tools,
 * ask_user and read_skill (view tools offer no analysis and are left out)
 * and their names; the clarifier's shell catalog; and the run context (its
 * elements are the gate's).
 */
export const openCorpusRun = async (key = 'Demo') => {
    const copilot = await openCopilotRun(key);
    const run = copilot.props;
    const context = buildRunContext(run);
    const tools = [
        ...buildRunTools(run, { cards: [] }),
        ...copilot.tools,
        buildAskUserTool(),
        buildReadSkillTool(bundledSkills(), { load: copilotSkills.load })
    ].filter(Boolean);
    // As the turn builds it (turn.js turnChecks): the run's shell catalogue
    // with its run settings.
    const clarifier = buildClarifier({ plotFiles: run.plotFiles, structure: run.structure, runSettings: run.runSettings, symmetry: run.symmetry },
        { noiseFloor: context.symmetry?.noise_floor_A ?? null });
    return {
        tools,
        offered: tools.map((tool) => tool.name),
        catalog: clarifier.catalog,
        context
    };
};

/**
 * An entry's recorded steps, run as the agent loop runs a model's calls
 * (executeToolCall: arguments parsed and validated): `results` per step, and
 * `steps` as the loop records them with the verification — each with its
 * status, error and facts. A question, a skill read or a view tool is not run
 * here (it returns no facts) and is recorded as the loop records it: no
 * analysis.
 */
export const runSteps = async (tools, steps) => {
    const results = [];
    const recorded = [];
    for (const step of steps) {
        const tool = tools.find((candidate) => candidate.name === step.name);
        if (!tool) throw new Error(`the Demo offers no tool named ${step.name}`);
        const kind = { interactive: Boolean(tool.interactive), reading: Boolean(tool.keep), view: Boolean(tool.view) };
        if (kind.interactive || kind.reading || kind.view) {
            results.push({ step, result: null });
            recorded.push({ name: step.name, args: step.args || {}, status: 'ok', error: null, executed: true, ...kind });
            continue;
        }
        const result = await executeToolCall(tool, JSON.stringify(step.args || {}), { timeoutMs: 120_000 });
        results.push({ step, result });
        recorded.push({
            name: step.name, args: result.args, status: stepStatus(result), error: result.error || null, executed: true, ...kind, facts: result.facts
        });
    }
    return { results, steps: recorded };
};
