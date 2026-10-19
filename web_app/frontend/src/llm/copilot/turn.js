// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The turn: one module runs one question. Its interface is two functions:
//
//   runTurn(spec)      → AsyncIterable<TurnEvent>   the question, run
//   describeTurn(spec) → TurnDescription            the same assembly, not run (pure)
//
// Behind that interface sits everything a question needs besides the screen,
// which the dock (useCopilot.js) and the validation harness
// (validation/harness.js) each used to assemble, and which had drifted apart
// between them: the run context and the context as the model reads it, the
// tools offered (the caller's, then ask_user, then read_skill), the system
// prompt for those tools, the history window, the Workbench's own questions
// (clarify.js), the evidence gate (evidence.js), the answer checks' context
// facts (answerAudit.js), the context for {{ctx.…}} references, the loop
// (agentLoop.js) with its model calls through the model port
// (provider/model.js), and one terminal event that says how the turn ended.
// The dock and the harness are its two adapters; neither builds any of this.
// The module keeps no state between calls: the caller passes the history, the
// clock and whether tools are on (the dock remembers which models refused
// them).
//
// The spec:
//   question      the user's question (required)
//   run           CopilotChat's eight run props, as the run module
//                 (src/copilotRun.js) gives them (required)
//   model         a model port, or a bare stream function (required to run)
//   tools         the caller's tools, in order: the tools over the run props
//                 (runTools.js), the module's view tools, the host's tools.
//                 Never ask_user or read_skill: the turn adds those itself
//                 (a TypeError otherwise)
//   askUser       (request, { signal }) → Promise<{ source, index?, answer? }>,
//                 rejecting with an AbortError on abort; null: no ask_user
//                 tool and no Workbench questions
//   skills        the host's bundle { files, references?, load? } (parsed
//                 here; a malformed skill is `skillsError`, never thrown), or
//                 already parsed skills { list, load? }, or null
//   history       [{ question, answer }], oldest first: the conversation so far
//   signal        aborts the turn (the `stopped` event)
//   clock         { now() } in ms (Date.now by default): turn time and `at`
//   toolContext   more context for a tool's run() beside its signal and
//                 askUser (the dock's card-focus controller)
//   config        { toolsEnabled = true, method = 'full', limits = {} }
//
// Invariants:
//   - ask_user is offered only with an askUser port and tools on; read_skill
//     only with skills, tools on and method 'full'; both after `tools`, in
//     that order;
//   - tools are active exactly when config.toolsEnabled holds and at least
//     one tool is offered; inactive, the request has no tools and the prompt
//     is the no-tools one;
//   - the history becomes messages: an unanswered question is dropped, each
//     answer is compacted to MEMORY_ANSWER_CHARS, and the newest
//     HISTORY_MESSAGES messages are kept (prompts/copilot.js);
//   - invalid input throws (a TypeError) on the first next(); after `start`
//     the iterator never throws.
//
// The events: every event the loop emits (agentLoop.js), shape unchanged,
// each with `at` (clock.now()), and:
//   { type: 'start', offered, prompt, method, toolsActive, shown, skillsError }
//       first, exactly once: describeTurn(spec) for this question
//   { type: 'request', round, why, sent }      before each model call (the loop's)
//   exactly one terminal event, always last:
//   { type: 'done', … }                        the loop's, unchanged
//   { type: 'stopped', error }                 aborted (an AbortError)
//   { type: 'failed', error, message, hint }   any other error: the message and
//                                              what to do (robustness.js
//                                              explainModelError); a context-limit
//                                              refusal and a spend cap land here
// Order: start; a Workbench question (a tool_start/tool_end pair, `app`);
// then per request: round, request, text/reasoning deltas, tool_start/tool_end
// pairs (sharing an id); a gate is always followed by a retract;
// tools_unsupported at most once, before any tool step, followed by a new
// round and request under the no-tools prompt.
//
// Seams for the next pieces: arm C's method switch goes in planTurn (one
// function decides the prompt sections, the tool descriptions, the skills
// offered and the gate's wording); the verification module replaces the four
// checks loopChecks builds; provider quirks stay behind the model port.

import { runCopilotTurn } from './agentLoop';
import { contextAuditFacts } from './answerAudit';
import { ASK_USER, buildAskUserTool } from './askUser';
import { buildClarifier } from './clarify';
import { buildEvidenceGate } from './evidence';
import { compactAnswer, explainModelError } from './robustness';
import { READ_SKILL, buildReadSkillTool, parseSkillFiles, skillIndex } from './skills';
import { usableTools } from './tools';
import { pointGroupFacts } from '../context/pointGroups';
import { buildRunContext } from '../context/runContext';
import { buildCopilotMessages, COPILOT_PROMPT_NO_TOOLS, shownContext } from '../prompts/copilot';
import { asModelPort } from '../provider/model';

// The longest earlier answer re-sent whole as memory.
export const MEMORY_ANSWER_CHARS = 2000;
// The methods a turn can follow: 'full' (the skills and the prompt's method
// rules). Arm C's 'off' joins when its switch is built.
export const TURN_METHODS = Object.freeze(['full']);
// The tools the turn adds itself, never the caller.
export const TURN_OWN_TOOLS = Object.freeze([ASK_USER, READ_SKILL]);

const DEFAULT_CLOCK = Object.freeze({ now: () => Date.now() });

const isObject = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);
const invalid = (message) => {
    throw new TypeError(`turn: ${message}`);
};

// The spec, checked and with its defaults. `running` asks for what running
// needs beyond a description: the question and the model.
const readSpec = (spec, { running }) => {
    if (!isObject(spec)) invalid('the spec must be an object');
    const {
        question, run, model = null, tools = [], askUser = null, skills = null, history = [], signal = null,
        clock = DEFAULT_CLOCK, toolContext = null, config = {}
    } = spec;
    if (running && (typeof question !== 'string' || !question.trim())) invalid('question must be a non-empty string');
    if (!isObject(run)) invalid('run must be the run props (CopilotChat\'s eight run props)');
    if (running) asModelPort(model);
    if (!Array.isArray(tools)) invalid('tools must be an array');
    const own = tools.find((tool) => tool && TURN_OWN_TOOLS.includes(tool.name));
    if (own) invalid(`tools must not contain ${own.name}: the turn adds it (askUser, skills)`);
    if (askUser != null && typeof askUser !== 'function') invalid('askUser must be a function or null');
    if (skills != null && !(isObject(skills) && (isObject(skills.files) || Array.isArray(skills.list)))) {
        invalid('skills must be { files, references?, load? }, { list, load? } or null');
    }
    if (!Array.isArray(history) || !history.every(isObject)) invalid('history must be an array of { question, answer }');
    if (signal != null && typeof signal.aborted !== 'boolean') invalid('signal must be an AbortSignal or null');
    if (!clock || typeof clock.now !== 'function') invalid('clock must have now()');
    if (toolContext != null && !isObject(toolContext)) invalid('toolContext must be an object or null');
    if (!isObject(config)) invalid('config must be an object');
    const { toolsEnabled = true, method = 'full', limits = {} } = config;
    if (typeof toolsEnabled !== 'boolean') invalid('config.toolsEnabled must be a boolean');
    if (!TURN_METHODS.includes(method)) invalid(`config.method must be one of ${TURN_METHODS.join(', ')}`);
    if (!isObject(limits)) invalid('config.limits must be an object');
    return {
        question, run, model, tools, askUser, skills, history, signal, clock, toolContext,
        config: { toolsEnabled, method, limits }
    };
};

// The skills a turn may offer, { list, load, error }: the host's bundle is
// parsed here, and a malformed skill is reported, never thrown.
const readSkills = (skills) => {
    if (!skills) return { list: [], load: null, error: null };
    const load = typeof skills.load === 'function' ? skills.load : null;
    if (Array.isArray(skills.list)) return { list: skills.list, load, error: null };
    try {
        return { list: parseSkillFiles(skills.files, { available: skills.references || null }), load, error: null };
    } catch (failure) {
        return { list: [], load, error: failure?.message || String(failure) };
    }
};

// The assembly, pure: the context, the tools offered, the prompt they get.
// Arm C's method switch belongs here.
const planTurn = ({ run, tools, askUser, skills, config }) => {
    const context = buildRunContext(run);
    const { list, load, error } = readSkills(skills);
    const on = config.toolsEnabled;
    const asks = on && typeof askUser === 'function';
    const reads = on && config.method === 'full' && list.length > 0;
    const offered = [
        ...usableTools(tools),
        ...(asks ? [buildAskUserTool()] : []),
        ...(reads ? [buildReadSkillTool(list, { load })] : [])
    ];
    const toolsActive = on && offered.length > 0;
    const view = offered.some((tool) => Boolean(tool.view));
    return {
        context,
        shown: shownContext(context),
        skills: list,
        skillsError: error,
        offered: {
            tools: offered,
            skills: reads ? list.map((skill) => skill.name) : [],
            askUser: asks,
            view
        },
        toolsActive,
        prompt: toolsActive ? (view ? 'view' : 'tools') : 'no_tools',
        method: config.method
    };
};

const describe = (plan) => ({
    offered: plan.offered,
    prompt: plan.prompt,
    method: plan.method,
    toolsActive: plan.toolsActive,
    shown: plan.shown,
    skillsError: plan.skillsError
});

// The earlier questions as messages: each answered one with its answer
// compacted; an unanswered one (stopped or failed before any text) dropped.
const historyMessages = (history) => history.flatMap(({ question, answer }) => {
    const text = typeof answer === 'string' ? answer.trim() : '';
    if (!text || typeof question !== 'string' || !question.trim()) return [];
    return [{ role: 'user', content: question }, { role: 'assistant', content: compactAnswer(text, MEMORY_ANSWER_CHARS) }];
});

const turnMessages = (plan, { question, history }) => buildCopilotMessages({
    context: plan.context,
    history: historyMessages(history),
    userText: question,
    toolsEnabled: plan.toolsActive,
    viewTools: plan.prompt === 'view',
    askUser: plan.toolsActive && plan.offered.askUser,
    skills: plan.toolsActive && plan.offered.skills.length ? skillIndex(plan.skills) : ''
});

// The four checks the loop runs with: the Workbench's own questions (only
// with someone to ask), the evidence gate (on the run's shells, its
// elements and the tools offered, tools on or off), the context facts the
// answer must not contradict, and the context {{ctx.…}} references read: the
// one the model was shown. The verification module replaces them with one
// object.
const loopChecks = (plan, { run, askUser }) => {
    const clarifier = buildClarifier(
        { plotFiles: run.plotFiles, structure: run.structure, symmetry: run.symmetry },
        { noiseFloor: plan.context.symmetry?.noise_floor_A ?? null }
    );
    return {
        clarifier: askUser ? clarifier : null,
        evidence: buildEvidenceGate({
            catalog: clarifier.catalog,
            offered: plan.offered.tools.map((tool) => tool.name),
            elements: Object.keys(plan.context.structure?.elements || {})
        }),
        auditContext: contextAuditFacts(plan.context, pointGroupFacts),
        referenceContext: plan.shown
    };
};

const failed = (error) => ({ type: 'failed', error, ...explainModelError(error) });

/**
 * What a turn with this spec offers and sends, without running it: `{ offered:
 * { tools, skills, askUser, view }, prompt: 'tools' | 'view' | 'no_tools',
 * method, toolsActive, shown, skillsError }`. Pure; `question` and `model`
 * are not needed. `offered.tools` are the tool objects in the order sent,
 * `offered.skills` the skill names read_skill offers.
 */
export const describeTurn = (spec) => describe(planTurn(readSpec(spec, { running: false })));

/**
 * One question, as events (see the header). Throws a TypeError on invalid
 * input at the first next(); never throws after `start`.
 */
export async function* runTurn(spec) {
    const turn = readSpec(spec, { running: true });
    const plan = planTurn(turn);
    const messages = turnMessages(plan, turn);
    const checks = loopChecks(plan, turn);
    const stamp = (event) => ({ ...event, at: turn.clock.now() });
    yield stamp({ type: 'start', ...describe(plan) });
    let answered = false;
    try {
        for await (const event of runCopilotTurn({
            model: turn.model,
            messages,
            tools: plan.offered.tools,
            limits: turn.config.limits,
            signal: turn.signal,
            toolsEnabled: plan.toolsActive,
            noToolsSystemPrompt: COPILOT_PROMPT_NO_TOOLS,
            toolContext: { ...(turn.toolContext || {}), ...(turn.askUser ? { askUser: turn.askUser } : {}) },
            clock: turn.clock,
            ...checks
        })) {
            if (event.type === 'done') answered = true;
            yield stamp(event);
        }
        if (!answered) yield stamp(failed(new Error('The turn ended without an answer.')));
    } catch (error) {
        yield stamp(error?.name === 'AbortError' ? { type: 'stopped', error } : failed(error));
    }
}
