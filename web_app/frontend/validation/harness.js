// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The live arms: the Copilot itself on a case, exactly as the app assembles
// it (CopilotChat + useCopilot), minus the screen.
//
//   B — tools, skills (read_skill + the skills' index), guardrails (tool
//       refusals, the Workbench's window questions, the answer checks);
//   C — the same without the skills: the ablation.
//
// The view tools (highlight_cards, orient_3d_view) need a screen and are left
// out of both arms; every tool offered only computes. The simulated user
// skips every question (the model must use its stated default), and the run
// and its files have one neutral name for every case (runFolder.js RUN_NAME:
// a case id names its kind), so nothing about the planted truth reaches the
// model except through the tools. The record keeps the case id; the model
// never sees the record.
//
// The run is opened through the run module's in-process adapter
// (runFolder.js loadRun → src/openRunFolder.js), the app's own path from a
// folder's files to the Copilot's run props, tools and engine, with a view
// (HARNESS_VIEW: the app's defaults). The record keeps the view and the
// context sections the model was given, and that the question was one fresh
// question with no earlier turns; `pcaTable` adds the PCA table, as after a
// visit to the PCA Ellipsoid page, for a dev check of its effect (off by
// default).
//
// One run gives one record: every decision (round, tool, arguments, the
// reasoning and text the model wrote before it, the result it read), the
// final answer, its verdict and score, the answer checks against the
// contradictions the score found, where each reported number came from,
// and the safety counts. The record is the transcript's source.
//
// THY-18 (protocol 2026-10-17.2): the evidence gate is on (an answer whose
// analysis did not run goes back to the model; each send-back is recorded),
// numbers may be written by reference (the record keeps the answer as written
// and as filled; the VERDICT is parsed from the filled text), and the Outcome
// is the app's computed verdict, recorded beside the case score.

import { runCopilotTurn } from '../src/llm/copilot/agentLoop';
import { buildAskUserTool, ASK_USER } from '../src/llm/copilot/askUser';
import { buildClarifier } from '../src/llm/copilot/clarify';
import { buildEvidenceGate } from '../src/llm/copilot/evidence';
import { buildReadSkillTool, parseSkillFiles, READ_SKILL, skillIndex } from '../src/llm/copilot/skills';
import { buildRunTools } from '../src/llm/copilot/runTools';
import { contextAuditFacts } from '../src/llm/copilot/answerAudit';
import { buildRunContext } from '../src/llm/context/runContext';
import { pointGroupFacts } from '../src/llm/context/pointGroups';
import { buildCopilotMessages, shownContext } from '../src/llm/prompts/copilot';
import { HARNESS_VIEW, loadRun } from './runFolder';
import { parseVerdict, VERDICT_INSTRUCTIONS, VERDICT_MARK } from './verdict';
import { outcomeAgreement, scoreVerdict } from './scoring';
import { shownValues, traceNumbers } from './provenance';

// Bumped when the cases, question, scoring, harness, tools or skills change what
// a run means (2026-10-17.1: validation batch 1's tools; 2026-10-17.2:
// THY-18's evidence gate, numbers by reference and computed verdict;
// 2026-10-19.1: the run is named RUN_NAME, no longer after the case id, which
// told the model the case's kind).
export const PROTOCOL_VERSION = '2026-10-19.1';

export const ARMS = Object.freeze({
    B: { skills: true, label: 'Copilot with skills and guardrails' },
    C: { skills: false, label: 'Copilot without skills (ablation)' }
});

/** The simulated user: skips every question, so the model uses its own default. */
export const skippingUser = (log) => async (request) => {
    log.push({ question: request.question, options: request.options.map((option) => option.label), reply: 'skipped' });
    return { source: 'skipped' };
};

const clip = (text, limit) => (text && text.length > limit ? `${text.slice(0, limit)}… [${text.length - limit} more characters]` : text);

/** The question each live run is asked: the case question plus the verdict format. */
export const liveQuestion = (built) => `${built.question}\n\n${VERDICT_INSTRUCTIONS}`;

// The sections buildRunContext (src/llm/context/runContext.js) can give the
// model besides the run's name and Live Data flag, in its order. A run in the
// app can have them all (the Demo with the PCA table does: runFolder.test.js);
// a case has the structure, its symmetry and the partial g(r) only, unless the
// view adds the PCA table.
export const CONTEXT_SECTIONS = Object.freeze([
    'structure', 'configuration_optimization', 'run_settings', 'symmetry', 'pca_displacements',
    'pair_correlations', 'datasets', 'partial_weights', 'convergence'
]);

/**
 * What the model was told about the run, for the record: the context sections
 * it had and lacked, and that the question came with no earlier turns (the
 * app sends up to eight).
 */
export const contextRecord = (context) => ({
    sections: CONTEXT_SECTIONS.filter((key) => key in context),
    absent: CONTEXT_SECTIONS.filter((key) => !(key in context)),
    other_sections: Object.keys(context).filter((key) => key !== 'run' && key !== 'live_mode' && !CONTEXT_SECTIONS.includes(key)),
    earlier_turns: 0
});

/**
 * One live run. `callModel(request)` streams a model reply (contextGuard's
 * streamCopilotChat bound to a model, or a scripted fake in tests); `skills`
 * is `{ files, references, load }` as copilotSkills.js gives the app. `view`
 * (merged over HARNESS_VIEW) and `pcaTable` set what is "on screen" (runFolder.js
 * loadRun); the defaults change nothing.
 */
export const runLiveCase = async ({
    built, arm, model, runIndex = 0, callModel, skills = null, commit = null, limits = {}, signal = null, now = () => Date.now(),
    view = HARNESS_VIEW, pcaTable = false
}) => {
    if (!ARMS[arm]) throw new Error(`unknown arm '${arm}' (B or C)`);
    const startedAt = now();
    const run = await loadRun(built, view, { pcaTable });
    const { props } = run;
    const context = buildRunContext(props);
    const parsedSkills = ARMS[arm].skills && skills?.files ? parseSkillFiles(skills.files, { available: skills.references || null }) : [];
    const tools = [
        ...buildRunTools(props, { cards: [] }),
        ...run.tools,
        buildAskUserTool(),
        ARMS[arm].skills ? buildReadSkillTool(parsedSkills, { load: skills?.load }) : null
    ].filter(Boolean);
    if (tools.some((tool) => tool.view)) throw new Error('a view tool was offered: the live arms offer compute-only tools');
    const question = liveQuestion(built);
    const messages = buildCopilotMessages({
        context, history: [], userText: question, toolsEnabled: true, viewTools: false, askUser: true,
        skills: ARMS[arm].skills ? skillIndex(parsedSkills) : ''
    });
    const asked = [];
    const clarifier = buildClarifier({ plotFiles: props.plotFiles, structure: props.structure, symmetry: props.symmetry },
        { noiseFloor: context.symmetry?.noise_floor_A ?? null });

    const record = {
        id: `${built.id}__${arm}__${model.name}__r${runIndex + 1}`,
        protocol: PROTOCOL_VERSION,
        commit,
        case: built.id,
        split: built.split,
        arm,
        arm_label: ARMS[arm].label,
        model: { name: model.name, label: model.label, id: model.model, base_url: model.baseUrl, temperature: model.temperature },
        run_index: runIndex + 1,
        started_at: new Date(startedAt).toISOString(),
        question,
        tools_offered: tools.map((tool) => tool.name),
        skills_offered: parsedSkills.map((skill) => skill.name),
        // The view the run was opened with (the run module's, normalised) and
        // what the model was told: its context sections, one fresh question.
        view: run.view,
        context: contextRecord(context),
        steps: [],
        retracted: [],
        gates: [],
        rounds: 0,
        final: null,
        error: null
    };

    // The round's reasoning and text, the "why" of the calls that follow it.
    let round = { index: 0, reasoning: '', text: '' };
    const open = new Map();
    let gated = null;
    try {
        for await (const event of runCopilotTurn({
            callModel, messages, tools, limits, signal, toolsEnabled: true,
            toolContext: { askUser: skippingUser(asked) },
            auditContext: contextAuditFacts(context, pointGroupFacts),
            clarifier,
            evidence: buildEvidenceGate({ catalog: clarifier.catalog, offered: tools.map((tool) => tool.name), elements: run.elements }),
            referenceContext: shownContext(context)
        })) {
            if (event.type === 'round') round = { index: event.round, reasoning: '', text: '' };
            else if (event.type === 'reasoning') round.reasoning += event.delta;
            else if (event.type === 'text') round.text += event.delta;
            else if (event.type === 'gate') {
                // The evidence gate sent the answer back: its text is withdrawn next.
                gated = { round: round.index, id: event.id, what: event.what, missing: event.missing, notice: event.notice };
                record.gates.push(gated);
            } else if (event.type === 'retract') {
                // The loop withdrew this round's text (a looping reply, a call
                // written as JSON text, which it then runs, or an answer the
                // evidence gate sent back): kept for the reviewer, never shown
                // as the reason for a call.
                record.retracted.push({ round: round.index, text: round.text, ...(gated ? { by: `evidence gate: ${gated.id}` } : {}) });
                gated = null;
                round.text = '';
            }
            else if (event.type === 'tool_start') {
                open.set(event.id, {
                    step: record.steps.length + 1,
                    round: round.index + (event.app ? 0 : 1),
                    id: event.id,
                    name: event.name,
                    by: event.app ? 'workbench' : 'model',
                    args: event.args,
                    why: event.app ? 'the Workbench\'s own question (clarify.js)' : { reasoning: round.reasoning, text: round.text }
                });
            } else if (event.type === 'tool_end') {
                const step = open.get(event.id) || { step: record.steps.length + 1, id: event.id, name: event.name };
                open.delete(event.id);
                record.steps.push({
                    ...step, args: event.args ?? step.args, status: event.status, checks: event.checks || [],
                    data: event.data, error: event.error || null, ms: event.ms, cached: Boolean(event.cached)
                });
            } else if (event.type === 'tools_unsupported') {
                record.tools_unsupported = event.error;
            } else if (event.type === 'done') {
                record.rounds = event.rounds;
                record.final = {
                    text: event.text, raw_text: event.rawText ?? event.text, outcome: event.outcome, stopped_by: event.stoppedBy,
                    cut: Boolean(event.cut), audit: event.audit || [], references: event.references ?? null,
                    reasoning: round.reasoning, reason: event.reason ?? null
                };
            }
        }
    } catch (failure) {
        record.error = { name: failure?.name || 'Error', message: failure?.message || String(failure) };
        if (failure?.name === 'SpendCapError') {
            record.finished_at = new Date(now()).toISOString();
            throw Object.assign(failure, { record: finishRecord(record, built, context, asked, startedAt, now) });
        }
    }
    return finishRecord(record, built, context, asked, startedAt, now);
};

const finishRecord = (record, built, context, asked, startedAt, now) => {
    record.finished_at = new Date(now()).toISOString();
    record.wall_ms = now() - startedAt;
    record.questions = asked;
    const text = record.final?.text || '';
    const { verdict, error } = parseVerdict(text);
    record.verdict = verdict;
    record.verdict_error = error;
    record.score = scoreVerdict(built.truth, verdict);

    // The answer checks against the contradictions the score found.
    const audit = record.final?.audit || [];
    const rulesFound = new Set(audit.map((finding) => finding.rule));
    const covered = record.score.contradictions.filter((item) => item.rule !== 'none');
    record.audit_eval = {
        findings: audit,
        contradictions: record.score.contradictions,
        caught: covered.filter((item) => rulesFound.has(item.rule)).map((item) => item.what),
        missed: covered.filter((item) => !rulesFound.has(item.rule)).map((item) => item.what),
        not_covered: record.score.contradictions.filter((item) => item.rule === 'none').map((item) => item.what),
        // A finding with no matching contradiction in the verdict.
        unmatched_findings: audit.filter((finding) => !covered.some((item) => item.rule === finding.rule)).map((finding) => finding.rule)
    };

    // Where each reported number came from: what the model was shown. The
    // filled answer is traced (a number written by reference is the engine's,
    // so it always matches); the references themselves are counted apart.
    const steps = record.steps.filter((step) => step.by === 'model' || step.by === 'workbench');
    // The prose only: the VERDICT block restates numbers the score judges.
    const marks = [...text.matchAll(VERDICT_MARK)];
    const prose = marks.length ? text.slice(0, marks[0].index) : text;
    record.provenance = traceNumbers(prose, shownValues({ steps, context, question: record.question }));
    const references = record.final?.references || { written: 0, filled: 0, unknown: [] };
    record.provenance.by_reference = { written: references.written, filled: references.filled, unknown: references.unknown };

    // THY-18: what the gate sent back, the evidence still missing at the end,
    // and the computed verdict beside the case score.
    record.evidence = {
        gates: record.gates.map(({ round, id, missing }) => ({ round, id, missing })),
        missing: audit.filter((finding) => finding.rule === 'evidence').map((finding) => ({ id: finding.id, missing: finding.missing }))
    };
    record.outcome_eval = outcomeAgreement(record.final?.outcome?.status ?? null, record.score.score);

    // Safety: a window marked as the user's choice when the user chose none.
    const unsafe = [];
    record.steps.forEach((step, index) => {
        const askedBefore = record.steps.slice(0, index).some((earlier) => earlier.name === ASK_USER);
        if (step.name === 'compute_bond_angles' && step.args?.window_confirmed === true && !askedBefore) {
            unsafe.push({ step: step.step, what: 'window_confirmed: true with no question to the user before it' });
        }
    });
    record.safety = {
        unsafe_actions: unsafe,
        unknown_tool_calls: record.steps.filter((step) => step.error && /^unknown tool/.test(step.error)).length,
        failed_calls: record.steps.filter((step) => step.error).length,
        questions_asked: asked.length,
        skills_read: record.steps.filter((step) => step.name === READ_SKILL).map((step) => step.args?.name ?? step.args?.skill ?? null)
    };
    record.truth = built.truth;
    record.final_reasoning_excerpt = clip(record.final?.reasoning || '', 4000);
    return record;
};
