// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The answer corpus's runner (docs/AGENT_PLAN.md item 3.5; the corpus is the
// entries/ folder beside this file, one JSON file per real answer): the run
// an answer was given on (the Demo, or a dev case of the validation) as the
// Copilot is handed it, an entry's recorded steps run as the agent loop runs
// them, and ONE adapter from an answer to what the current checks read in it.
//
// `readAnswer` is the only place that knows today's verification modules:
// the evidence gate (llm/copilot/evidence.js) and the checks on the answer
// (llm/copilot/answerAudit.js), through their public functions only. When
// one verification module replaces evidence / answerAudit / verdict, this
// function is rewritten and the corpus pins stay as they are.

import { copilotSkills } from '../../copilotSkills';
import { auditAnswer, contextAuditFacts } from '../../llm/copilot/answerAudit';
import { buildAskUserTool } from '../../llm/copilot/askUser';
import { buildClarifier } from '../../llm/copilot/clarify';
import { buildEvidenceGate } from '../../llm/copilot/evidence';
import { buildRunTools } from '../../llm/copilot/runTools';
import { buildReadSkillTool } from '../../llm/copilot/skills';
import { executeToolCall } from '../../llm/copilot/tools';
import { pointGroupFacts } from '../../llm/context/pointGroups';
import { buildRunContext } from '../../llm/context/runContext';
import { buildCase, DEV_CASES } from '../../../validation/cases';
import { loadRun } from '../../../validation/runFolder';
import { bundledSkills, loadDemoAssistantRun } from '../modelVisible/helpers';

/** What an entry pins for the gate when its transcript does not show which tools ran. */
export const STEPS_UNKNOWN = 'steps unknown';

/**
 * The adapter: an answer on `run` (openCorpusRun) → what the current checks
 * read in it. The gate is built as the dock builds it for a turn
 * (useCopilot.js: the clarifier's shell catalog, the offered tools, the
 * context's elements), the audit's context facts likewise.
 *   claims         the gate's rules whose claim the answer makes (in a
 *                  sentence about the run), read as `unmet` reads them with
 *                  nothing run: every rule that applies is then unmet;
 *   gateSendBack   the rule the gate would send the answer back for, given
 *                  the steps that ran (`check`, the loop's call): { id,
 *                  missing } or null; STEPS_UNKNOWN when `ran` is null;
 *   auditFindings  the checks on the answer that fire (their rule names),
 *                  against the steps' facts and the run context's.
 */
export const readAnswer = (text, { run, ran = null, facts = [] }) => {
    const gate = buildEvidenceGate({ catalog: run.catalog, offered: run.offered, elements: run.elements });
    const sendBack = ran ? gate.check(text, { ran }) : null;
    return {
        claims: gate.unmet(text, []).map((rule) => rule.id),
        gateSendBack: ran ? (sendBack && { id: sendBack.id, missing: sendBack.missing }) : STEPS_UNKNOWN,
        auditFindings: auditAnswer(text, { facts, context: contextAuditFacts(run.context, pointGroupFacts) })
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
 * `{ tools, offered, catalog, elements, context }`: the run tools, the
 * engine tools, ask_user and read_skill (view tools offer no analysis and
 * are left out) and their names; the clarifier's shell catalog; the
 * context's elements; and the run context.
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
    // As the turn builds it (turn.js loopChecks): the run's shell catalogue
    // with its run settings.
    const clarifier = buildClarifier({ plotFiles: run.plotFiles, structure: run.structure, runSettings: run.runSettings, symmetry: run.symmetry },
        { noiseFloor: context.symmetry?.noise_floor_A ?? null });
    return {
        tools,
        offered: tools.map((tool) => tool.name),
        catalog: clarifier.catalog,
        elements: Object.keys(context.structure?.elements || {}),
        context
    };
};

/**
 * An entry's recorded steps, run as the agent loop runs a model's calls
 * (executeToolCall: arguments parsed and validated): `results` per step,
 * the tools' `facts`, and `ran` as the loop's ranCalls gives it to the gate
 * — executed without an error, and not a question, a skill read or a view
 * tool (those are not run here: they return no facts).
 */
export const runSteps = async (tools, steps) => {
    const results = [];
    for (const step of steps) {
        const tool = tools.find((candidate) => candidate.name === step.name);
        if (!tool) throw new Error(`the Demo offers no tool named ${step.name}`);
        if (tool.interactive || tool.keep || tool.view) {
            results.push({ step, analysis: false, result: null });
            continue;
        }
        results.push({ step, analysis: true, result: await executeToolCall(tool, JSON.stringify(step.args || {}), { timeoutMs: 120_000 }) });
    }
    const analyses = results.filter(({ analysis, result }) => analysis && !result.error);
    return {
        results,
        ran: analyses.map(({ step, result }) => ({ name: step.name, args: result.args })),
        facts: analyses.map(({ result }) => result.facts).filter(Boolean)
    };
};
