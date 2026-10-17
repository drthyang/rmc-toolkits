// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { toolBudgetNote } from '../prompts/copilot';
import { auditAnswer } from './answerAudit';
import { ASK_USER, MAX_QUESTIONS, answerResult } from './askUser';
import { resolveClarification } from './clarify';
import { parseOutcome } from './outcome';
import { createReferenceStore, fillReferences } from './references';
import { computeOutcome } from './verdict';
import { keyNumbers, repeatsItself, textToolCalls, withDeadlines } from './robustness';
import {
    accumulateToolCall, canonicalJson, finalizeToolCalls, parseToolArguments, validateToolArguments
} from './toolCalls';
import {
    executeToolCall, fitToBudget, stepStatus, toolCards, toolLabel, toRequestTools, usableTools
} from './tools';

// One Copilot turn: plan → act → observe → assess. The model is called with
// the tools; when it asks for some, each call is validated, run and answered
// with a `tool` message (result + deterministic checks), and the model is
// called again — until it answers without tools, or the round budget is used
// and a last request (tool_choice 'none') asks for the answer. An async
// generator of UI events; no React. `callModel` is injected (the page wraps
// provider/client.js streamChat), so tests script the model's replies.
//
// Events:
//   { type: 'round', round }                       a model request starts
//   { type: 'text', delta } / { type: 'reasoning', delta }
//   { type: 'tool_start', id, name, label, args, page, cards, interactive }
//   { type: 'tool_end', id, name, status, args, data, checks, error, ms, cached, pageParams }
//   { type: 'tools_unsupported', error }           retried without tools
//   { type: 'gate', id, what, missing, notice }    the evidence gate sent the answer back
//                                                  (a 'retract' withdraws its text)
//   { type: 'done', text, rawText, outcome, rounds, stoppedBy: 'answer' | 'max_rounds' | 'empty', cut, audit,
//     references, gates }
//                                                  text: the answer with its {{references}} filled in
//                                                  (rawText: as the model wrote it)
//                                                  cut: the last reply stopped at the length limit
//                                                  audit: findings on the answer (answer checks,
//                                                  unfilled references, missing evidence)
//                                                  outcome: computed from the checks (verdict.js)
// Errors (HTTP, network) and aborts are thrown to the caller.
//
// Budgets: at most `maxRounds` model requests with tool calls and
// `maxCallsPerRound` calls in each; one call may take `toolTimeoutMs`; after
// `maxTurnMs` no new round of tools starts. Each result is fitted to
// `resultCharBudget` characters, and the turn's tool messages together to
// `turnCharBudget` — older results are compacted to their checks first — so a
// small local model's context window keeps the run context and the rules.
// A tool may carry its own `resultBudget` (a skill is longer than an analysis
// result) and `keep: true` (its message is never compacted and does not count
// toward the turn budget: the agent follows a skill for the whole question).
// An interactive tool (ask_user) waits for the user: it has no time limit, the
// wait does not count toward `maxTurnMs`, and a turn asks at most
// `maxQuestions` questions, one per step.
//
// `clarifier` (clarify.js) lets the Workbench ask itself where the run data
// show that a choice changes the result: before the first request
// (`beforeTurn`: the answer joins the user's question for the model) and
// before an analysis runs in a window that cuts a shell (`beforeCall`: it runs
// in the window the user picks). These questions are ask_user steps
// (`app: true`) and count toward `maxQuestions`.
//
// Evidence, numbers and the verdict (THY-18). Every tool message carries a
// `ref` (r1, r2, …) the model quotes values by ({{r2.peak.significance_sigma}},
// references.js); the answer is shown, recorded and checked with them filled
// in from what the tools returned. `evidence` (evidence.js buildEvidenceGate)
// sends an answer that concludes something without its analysis back to the
// model, once per rule and at most `maxEvidenceGates` times a turn, when tools
// can still run; evidence still missing at the end is a finding. The Outcome
// is computed from the steps' checks and the findings (verdict.js); the
// model's Outcome line only explains it.
//
// On a local model (robustness.js): a model request that sends nothing for
// `stallMs`, or runs past `requestMaxMs`, is cancelled and explained; a reply
// with neither an answer nor a call (a reasoning-only reply), or one stopped
// because it kept writing the same passage (`repeatsItself`), gets one
// reminder, asked a little warmer (`antiRepeat` on the request: the host
// raises the temperature, as greedy decoding is what loops), before the turn
// ends on an explained empty reply; a call written
// as JSON text is run; and two rounds in a row that add nothing new (every
// call repeated or failed) end the analyses: the next request asks for the
// answer (`stoppedBy: 'no_progress'`), so repetition never spends the round
// budget.

export const DEFAULT_LIMITS = {
    maxRounds: 6,
    maxCallsPerRound: 4,
    toolTimeoutMs: 120000,
    maxTurnMs: 600000,
    resultCharBudget: 3000,
    turnCharBudget: 12000,
    maxQuestions: MAX_QUESTIONS,
    stallMs: 180000,
    requestMaxMs: 600000,
    maxNoProgressRounds: 2,
    maxEvidenceGates: 2
};

// A 4xx saying the model or server cannot take a `tools` field (Ollama:
// "<model> does not support tools"). Deliberately narrow: a 400 about a
// tool's schema is a real error and must be shown, not hidden by a fallback.
const TOOLS_UNSUPPORTED = /(does not|doesn't|do not) support (tools|tool calling|tool use|function calling)|(tools|tool calling|tool use|function calling) (is |are )?not supported/i;
export const isToolsUnsupportedError = (error) => (
    Number.isInteger(error?.status) && error.status >= 400 && error.status < 500
    && TOOLS_UNSUPPORTED.test(error.message || '')
);

// Assistant content that is only whitespace (models often stream a newline
// before a tool call) is sent as null: Anthropic rejects empty text blocks.
export const assistantContent = (text) => (typeof text === 'string' && text.trim() ? text : null);

const compactToolMessage = (content) => {
    try {
        const { ref, checks, error, data, compacted } = JSON.parse(content);
        // Already compacted: compacting again would lose its numbers.
        if (compacted) return content;
        const key = error ? null : keyNumbers(data);
        return JSON.stringify({
            // The ref still resolves: the app keeps the result it sent.
            ...(ref ? { ref } : {}),
            compacted: `this earlier result was shortened to fit the context; its checks${key ? ' and headline numbers' : ''} are kept`,
            ...(error ? { error } : {}),
            ...(checks ? { checks } : {}),
            ...(key ? { key } : {})
        });
    } catch {
        return JSON.stringify({ compacted: 'this earlier result was shortened to fit the context' });
    }
};

// What "Open <page>" should load on that page to show this step's analysis:
// the tool's own mapping of its validated arguments, or null.
const pageParamsFor = (tool, args) => {
    if (!tool || typeof tool.pageParams !== 'function' || !args || typeof args !== 'object') return null;
    try {
        return tool.pageParams(args) ?? null;
    } catch {
        return null;
    }
};

const EMPTY_FINAL = 'The analysis budget for this question was used before the model wrote an answer. '
    + 'The analyses it ran are listed above.';
const EMPTY_STALLED = 'The model kept repeating analyses that added nothing new, and wrote no answer. '
    + 'The analyses it ran are listed above.';

// When a reply holds neither an answer nor an analysis call.
const emptyNudge = (why) => ({
    role: 'user',
    content: `Your last reply ${why === 'loop'
        ? 'kept repeating the same passage, so it was stopped'
        : `had no answer and no analysis call${why === 'cut' ? ' (it reached the length limit while reasoning)' : ''}`}. `
        + 'Answer the question now from the context and the results above, briefly, or call one analysis if '
        + 'one is still missing. End with the Outcome line.'
});

const EMPTY_REASONS = {
    loop: 'the model kept repeating the same passage, even after a reminder',
    cut: 'the model spent its reply on reasoning and reached the length limit',
    empty: 'the model wrote no answer, even after a reminder'
};

// A note for the model goes into the last user message when that is the last
// message: some chat templates (Gemma's) refuse two user messages in a row.
const addNote = (messages, note) => {
    const last = messages[messages.length - 1];
    if (last?.role === 'user' && typeof last.content === 'string') {
        messages[messages.length - 1] = { ...last, content: `${last.content}\n\n${note.content}` };
    } else {
        messages.push(note);
    }
};

// How often (in new characters) a streaming reply is checked for a loop.
const LOOP_CHECK_CHARS = 400;

// When the last rounds added nothing new.
const noProgressNote = (rounds) => ({
    role: 'user',
    content: `Your last ${rounds} steps repeated earlier analyses or failed, so they added nothing new. `
        + 'Do not call them again. Answer now from the results above, say what is still missing, or ask the user '
        + 'with ask_user if a choice blocks you. End with the Outcome line.'
});

// The tool message the model reads: the result (fitted to the budget) and
// the checks, or the error. Checks are never trimmed: they are the verdict.
// Returns the message text and the fitted data, which the step shows as the
// result sent to the model.
const toolMessageContent = ({ data, checks, error, note }, budget, ref = null) => {
    const message = {};
    if (ref && !error) message.ref = ref;
    if (error) message.error = error;
    if (note) message.note = note;
    if (checks?.length) message.checks = checks;
    const head = JSON.stringify(message).length;
    if (!error) message.data = fitToBudget(data, Math.max(400, budget - head));
    return { content: JSON.stringify(message), data: error ? null : message.data };
};

export async function* runCopilotTurn({
    callModel,
    messages: initialMessages,
    tools = [],
    limits = {},
    signal,
    toolsEnabled = true,
    noToolsSystemPrompt = null,
    // Passed to every tool's run() beside its signal (the card-focus
    // controller for view tools).
    toolContext = null,
    // The run context's facts the answer must not contradict
    // (answerAudit.js contextAuditFacts).
    auditContext = null,
    // The Workbench's own questions (clarify.js buildClarifier).
    clarifier = null,
    // The evidence gate (evidence.js buildEvidenceGate), or null for none.
    evidence = null,
    // The run context, for {{ctx.…}} references.
    referenceContext = null
}) {
    const {
        maxRounds, maxCallsPerRound, toolTimeoutMs, maxTurnMs, resultCharBudget, turnCharBudget, maxQuestions,
        stallMs, requestMaxMs, maxNoProgressRounds, maxEvidenceGates
    } = { ...DEFAULT_LIMITS, ...limits };
    // The turn's referable results, its steps, and the evidence gate's state.
    const references = createReferenceStore({ context: referenceContext });
    const stepLog = [];
    const gatesFired = new Set();
    const gates = [];
    const ranCalls = () => stepLog.filter((step) => !step.error && !step.interactive && !step.reading)
        .map((step) => ({ name: step.name, args: step.args }));
    let noProgressRounds = 0;
    let stalled = false;
    let emptyNudges = 0;
    let warmer = false;
    let finalNoted = false;
    const startedAt = Date.now();
    // Time spent waiting for the user's answers, which the turn limit skips.
    let waitedMs = 0;
    let questions = 0;
    // The tools' facts this question, for the answer check.
    const facts = [];
    const audit = (answer) => auditAnswer(answer, { facts, context: auditContext });
    const available = toolsEnabled ? usableTools(tools) : [];
    const byName = new Map(available.map((tool) => [tool.name, tool]));
    const requestTools = available.length ? toRequestTools(available) : null;
    const messages = [...initialMessages];
    const cache = new Map();
    let useTools = Boolean(requestTools);
    let rounds = 0;
    let callCounter = 0;
    // Indices (into `messages`) of this turn's tool messages, oldest first.
    const toolMessages = [];

    // The user's question: the last user message as the turn starts.
    const questionIndex = messages.map((message) => message.role).lastIndexOf('user');
    // With the user's answers to the Workbench's questions appended: a window
    // the user chose there is the user's, as if written in the question.
    let questionText = questionIndex >= 0 && typeof messages[questionIndex].content === 'string' ? messages[questionIndex].content : '';
    const canAsk = () => Boolean(clarifier && useTools && typeof toolContext?.askUser === 'function' && questions < maxQuestions);
    let appQuestions = 0;
    // The window choices made in this turn (clarify.js), by (centre, partner, window).
    const clarifySession = { windows: new Map() };
    // One question the Workbench asks itself, shown as an ask_user step.
    async function* askFromApp(item) {
        appQuestions += 1;
        const id = `app_ask_${appQuestions}`;
        const args = {
            question: item.request.question,
            options: item.request.options,
            default: item.request.options[item.request.defaultIndex]?.label
        };
        yield { type: 'tool_start', id, name: ASK_USER, label: 'Question for you', args, page: null, cards: [], interactive: true, app: true };
        const started = Date.now();
        const reply = await toolContext.askUser(item.request, { signal });
        const ms = Date.now() - started;
        waitedMs += ms;
        questions += 1;
        const { data, checks } = answerResult(item.request, reply);
        yield { type: 'tool_end', id, name: ASK_USER, status: 'ok', args, data, checks, error: null, ms, cached: false, pageParams: null };
        return reply;
    }

    if (canAsk() && questionText) {
        const item = clarifier.beforeTurn(questionText);
        if (item) {
            const reply = yield* askFromApp(item);
            const { note } = resolveClarification(item, reply);
            questionText = `${questionText}\n\n(The Workbench asked: "${item.request.question}" — ${note})`;
            messages[questionIndex] = { ...messages[questionIndex], content: questionText };
        }
    }

    const keepWithinTurnBudget = (currentRoundStart) => {
        let total = toolMessages.reduce((sum, index) => sum + messages[index].content.length, 0);
        for (const index of toolMessages) {
            if (total <= turnCharBudget) break;
            if (index >= currentRoundStart) break;
            const before = messages[index].content.length;
            const compact = compactToolMessage(messages[index].content);
            if (compact.length < before) {
                messages[index] = { ...messages[index], content: compact };
                total -= before - compact.length;
            }
        }
    };

    for (;;) {
        const outOfTime = rounds > 0 && Date.now() - startedAt - waitedMs > maxTurnMs;
        const finalRound = useTools && (rounds >= maxRounds || outOfTime || stalled);
        // The last-round note is added once (a reminder can repeat that round).
        if (finalRound && !finalNoted) {
            finalNoted = true;
            addNote(messages, stalled ? noProgressNote(noProgressRounds) : toolBudgetNote(rounds));
        }
        yield { type: 'round', round: rounds };

        let text = '';
        let reasoning = '';
        let looped = false;
        let checkedAt = 0;
        let finishReason = null;
        const fragments = [];
        // The request's own signal: aborted with the turn, or by a deadline.
        const request = new AbortController();
        const onTurnAbort = () => request.abort();
        signal?.addEventListener('abort', onTurnAbort, { once: true });
        try {
            const stream = callModel({
                messages,
                tools: useTools ? requestTools : undefined,
                toolChoice: useTools ? (finalRound ? 'none' : 'auto') : undefined,
                signal: request.signal,
                ...(warmer ? { antiRepeat: true } : {})
            });
            for await (const chunk of withDeadlines(stream, { stallMs, maxMs: requestMaxMs, cancel: () => request.abort() })) {
                if (chunk.content) {
                    text += chunk.content;
                    yield { type: 'text', delta: chunk.content };
                }
                if (chunk.reasoning) {
                    reasoning += chunk.reasoning;
                    yield { type: 'reasoning', delta: chunk.reasoning };
                }
                if (chunk.toolCall) accumulateToolCall(fragments, chunk.toolCall);
                if (chunk.finishReason) finishReason = chunk.finishReason;
                if (text.length + reasoning.length - checkedAt >= LOOP_CHECK_CHARS) {
                    checkedAt = text.length + reasoning.length;
                    if (repeatsItself(reasoning) || repeatsItself(text)) {
                        looped = true;
                        request.abort();
                        break;
                    }
                }
            }
        } catch (error) {
            // Only the turn's first request falls back: later ones already
            // carried tools successfully, so their failure is a real one.
            if (useTools && rounds === 0 && !text && isToolsUnsupportedError(error)) {
                useTools = false;
                // The tool rules would only confuse a model that has none.
                if (noToolsSystemPrompt && messages[0]?.role === 'system') {
                    messages[0] = { ...messages[0], content: noToolsSystemPrompt };
                }
                yield { type: 'tools_unsupported', error: error.message };
                continue;
            }
            throw error;
        } finally {
            signal?.removeEventListener('abort', onTurnAbort);
        }

        if (looped) {
            // A looping reply is no answer: its text goes, and it is asked again.
            if (text) yield { type: 'retract' };
            text = '';
            fragments.length = 0;
        }
        let calls = useTools && !finalRound ? finalizeToolCalls(fragments, `${rounds}_`) : [];
        // A call written as JSON text (small local models do that) is run.
        if (!calls.length && useTools && !finalRound && text.trim()) {
            const written = textToolCalls(text, [...byName.keys()]);
            if (written.length) {
                calls = written.map((call, index) => ({ id: `call_text_${rounds}_${index}`, name: call.name, arguments: call.arguments, extraContent: null }));
                text = '';
                yield { type: 'retract' };
            }
        }
        const cut = finishReason === 'length';
        const why = looped ? 'loop' : (cut ? 'cut' : 'empty');
        if (!calls.length) {
            if (!text.trim() && emptyNudges < 1) {
                // A reasoning-only, empty or looping reply: one reminder, then
                // ask again (warmer after a loop).
                emptyNudges += 1;
                warmer = warmer || looped;
                addNote(messages, emptyNudge(why));
                continue;
            }
            if (!text.trim()) {
                // Ollama ignores tool_choice 'none': a last round can come back
                // as tool calls only. Say so instead of ending on nothing.
                if (finalRound) {
                    const final = stalled ? EMPTY_STALLED : EMPTY_FINAL;
                    yield { type: 'text', delta: final };
                    yield {
                        type: 'done',
                        text: final,
                        outcome: {
                            status: 'not achieved',
                            note: stalled ? 'the analyses stopped adding anything new and no answer was written' : 'no answer was written within the analysis budget',
                            body: final
                        },
                        rounds,
                        stoppedBy: stalled ? 'no_progress' : 'max_rounds',
                        cut,
                        audit: []
                    };
                } else {
                    yield {
                        type: 'done', text: '', outcome: null, rounds, stoppedBy: 'empty', cut, audit: [],
                        reason: EMPTY_REASONS[why]
                    };
                }
                return;
            }
            // The evidence gate: a conclusion whose analysis did not run goes
            // back to the model, once per rule, while tools can still run.
            const gate = evidence && useTools && !finalRound && gates.length < maxEvidenceGates
                ? evidence.check(text, { ran: ranCalls(), fired: gatesFired })
                : null;
            if (gate) {
                gatesFired.add(gate.id);
                gates.push({ round: rounds, id: gate.id, missing: gate.missing });
                yield { type: 'gate', id: gate.id, what: gate.what, missing: gate.missing, notice: gate.notice };
                yield { type: 'retract' };
                messages.push({ role: 'assistant', content: assistantContent(text) });
                messages.push({ role: 'user', content: gate.nudge });
                continue;
            }
            // Numbers by reference, then the checks on the filled answer, then
            // the verdict from the checks: the model's own grade never decides.
            const filled = fillReferences(text, references);
            const answer = filled.text;
            const findings = [
                ...audit(answer),
                ...(filled.unknown.length ? [{
                    rule: 'reference',
                    text: `The answer cites ${filled.unknown.map((name) => `{{${name}}}`).join(', ')}, which `
                        + `${filled.unknown.length === 1 ? 'names' : 'name'} no value this question's analyses returned: those numbers are missing.`
                }] : []),
                ...(evidence && useTools ? evidence.unmet(answer, ranCalls()).map((rule) => ({
                    rule: 'evidence',
                    id: rule.id,
                    missing: rule.missing,
                    text: `The answer concludes ${rule.what}, but ${rule.missing.join(', ')} did not run in this question.`
                })) : [])
            ];
            const stoppedBy = finalRound ? (stalled ? 'no_progress' : 'max_rounds') : 'answer';
            const outcome = computeOutcome({ text: answer, steps: stepLog, findings, stoppedBy, modelOutcome: parseOutcome(answer) });
            yield {
                type: 'done',
                text: answer,
                ...(answer !== text ? { rawText: text } : {}),
                outcome,
                rounds,
                stoppedBy,
                cut,
                audit: findings,
                references: { written: filled.written, filled: filled.filled, unknown: filled.unknown },
                gates
            };
            return;
        }

        rounds += 1;
        // The assistant message must list every call that the tool messages
        // below answer, and every listed call gets exactly one answer. Its
        // arguments are echoed as valid JSON (unparseable text as '{}': some
        // providers re-parse them and would reject the whole request); the
        // tool message says what was wrong. Empty content is null, never ''.
        // Provider data on a call (Gemini's thought_signature in
        // extra_content) is echoed unchanged, or the next request is refused.
        messages.push({
            role: 'assistant',
            content: assistantContent(text),
            tool_calls: calls.map((call) => {
                const parsed = parseToolArguments(call.arguments);
                const echoed = {
                    id: call.id,
                    type: 'function',
                    function: { name: call.name, arguments: parsed.ok ? JSON.stringify(parsed.value) : '{}' }
                };
                if (call.extraContent) echoed.extra_content = call.extraContent;
                return echoed;
            })
        });
        const roundStart = messages.length;
        let askedThisRound = false;
        // Whether this round brought anything new: a result that was not an
        // error and not repeated, or an answer from the user.
        let progressed = false;

        for (const [position, call] of calls.entries()) {
            callCounter += 1;
            const stepId = `${call.id}#${callCounter}`;
            const tool = byName.get(call.name);
            let parsed = parseToolArguments(call.arguments);
            let runArguments = call.arguments;
            let appNote = null;
            // The user answered a question for this call: news, even if the
            // call it settles was run before.
            let answered = false;
            // Before an analysis runs in a window that cuts a shell, the
            // Workbench asks which window (clarify.js beforeCall) — once per
            // window in a turn: a choice already made, or one the question
            // makes, is applied without asking; an identical call already
            // answered from the cache asks nothing.
            if (clarifier && tool && !tool.interactive && parsed.ok && position < maxCallsPerRound) {
                const validated = validateToolArguments(tool.parameters, parsed.value);
                const repeat = validated.ok && !tool.view && cache.has(`${call.name}:${canonicalJson(validated.value)}`);
                const item = validated.ok && !repeat ? clarifier.beforeCall(tool, validated.value, { questionText, session: clarifySession }) : null;
                if (item?.settled) {
                    runArguments = JSON.stringify(item.args);
                    parsed = { ok: true, value: item.args };
                    appNote = `the Workbench applied the user's choice: ${item.note}`;
                } else if (item && !askedThisRound && canAsk()) {
                    askedThisRound = true;
                    const reply = yield* askFromApp(item);
                    const resolved = resolveClarification(item, reply, { args: validated.value, tool, session: clarifySession });
                    runArguments = JSON.stringify(resolved.args);
                    parsed = { ok: true, value: resolved.args };
                    appNote = `the Workbench asked the user about the window: ${resolved.note}`;
                    answered = true;
                }
            }
            const shownArgs = parsed.ok ? parsed.value : call.arguments;
            yield {
                type: 'tool_start',
                id: stepId,
                name: call.name,
                label: tool ? toolLabel(tool, shownArgs) : call.name,
                args: shownArgs,
                page: tool?.page || null,
                cards: tool && parsed.ok ? toolCards(tool, parsed.value) : [],
                interactive: Boolean(tool?.interactive),
                reading: Boolean(tool?.keep)
            };

            let result;
            let note = null;
            let cached = false;
            if (position >= maxCallsPerRound) {
                result = { args: shownArgs, data: null, checks: [], ms: 0,
                    error: `skipped: at most ${maxCallsPerRound} tool calls per step — call it again in the next step if it is still needed` };
            } else if (!tool) {
                result = { args: shownArgs, data: null, checks: [], ms: 0,
                    error: `unknown tool "${call.name}"; available: ${[...byName.keys()].join(', ') || 'none'}` };
            } else if (tool.interactive && askedThisRound) {
                result = { args: shownArgs, data: null, checks: [], ms: 0,
                    error: 'one question at a time: ask the next one after this answer, if it is still needed' };
            } else if (tool.interactive && questions >= maxQuestions) {
                result = { args: shownArgs, data: null, checks: [], ms: 0,
                    error: `already asked ${questions} question${questions === 1 ? '' : 's'} for this request: `
                        + 'proceed with your default and say in the answer what you assumed' };
            } else if (tool.interactive) {
                // A question is a real exchange every time: never cached, and
                // no time limit while the user thinks. Only a question the
                // user saw counts toward the limit (not malformed arguments).
                askedThisRound = true;
                result = await executeToolCall(tool, call.arguments, { signal, timeoutMs: null, context: toolContext });
                waitedMs += result.ms;
                if (!result.error) questions += 1;
            } else {
                // Keyed on the validated arguments, so "2" and 2 are one call.
                const validated = parsed.ok ? validateToolArguments(tool.parameters, parsed.value) : null;
                // A view tool (tool.view) changes the screen, so a repeat is a
                // real request (point at the card again), never a cached result.
                const key = validated?.ok && !tool.view ? `${call.name}:${canonicalJson(validated.value)}` : null;
                if (key && cache.has(key)) {
                    result = { ...cache.get(key), ms: 0 };
                    cached = true;
                    note = 'identical call earlier in this question — the same result is repeated; change the arguments or answer';
                } else {
                    result = await executeToolCall(tool, runArguments, { signal, timeoutMs: toolTimeoutMs, context: toolContext });
                    if (key && !result.error) cache.set(key, result);
                }
            }
            if (appNote) note = note ? `${appNote}; ${note}` : appNote;
            if (answered || (!result.error && !cached)) progressed = true;

            if (result.facts && !result.error) facts.push(result.facts);
            if (!tool?.keep) toolMessages.push(messages.length);
            const budget = Number.isFinite(tool?.resultBudget) && tool.resultBudget > 0 ? tool.resultBudget : resultCharBudget;
            // A result the model can quote by reference: its ref, kept with
            // the data it was sent.
            const ref = result.error ? null : `r${references.size() + 1}`;
            const sent = toolMessageContent({ ...result, note }, budget, ref);
            if (ref) references.put(ref, sent.data);
            stepLog.push({
                name: call.name, args: result.args, status: stepStatus(result), error: result.error || null,
                interactive: Boolean(tool?.interactive), reading: Boolean(tool?.keep), ref
            });
            messages.push({ role: 'tool', tool_call_id: call.id, content: sent.content });
            keepWithinTurnBudget(roundStart);
            yield {
                type: 'tool_end',
                id: stepId,
                name: call.name,
                status: stepStatus(result),
                args: result.args,
                // What the model read: the result fitted to its budget.
                data: sent.data,
                checks: result.checks,
                error: result.error,
                ms: result.ms,
                cached,
                ref,
                pageParams: result.error ? null : pageParamsFor(tool, result.args)
            };
        }
        noProgressRounds = progressed ? 0 : noProgressRounds + 1;
        if (noProgressRounds >= maxNoProgressRounds) stalled = true;
    }
}
