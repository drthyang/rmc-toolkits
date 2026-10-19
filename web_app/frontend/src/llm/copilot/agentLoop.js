// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { asModelPort, ToolsUnsupportedError } from '../provider/model';
import { ASK_USER, MAX_QUESTIONS, answerResult } from './askUser';
import { resolveClarification } from './clarify';
import { beforeStall, keyNumbers, repeatsItself, textToolCalls, withDeadlines } from './robustness';
import {
    accumulateToolCall, canonicalJson, finalizeToolCalls, parseToolArguments, validateToolArguments
} from './toolCalls';
import {
    executeToolCall, fitToBudget, stepStatus, toolCards, toolLabel, toRequestTools, usableTools
} from './tools';
import { createTurnPolicy } from './turnPolicy';
import { createVerification } from './verification';

// One Copilot turn: plan → act → observe → assess. The model is called with
// the tools; when it asks for some, each call is validated, run and answered
// with a `tool` message (result + deterministic checks), and the model is
// called again — until it answers without tools, or the round budget is used
// and a last request (tool_choice 'none') asks for the answer. An async
// generator of UI events; no React. The turn module (turn.js) is its one
// production caller: it assembles the messages, the tools and the checks this
// loop runs with. `model` is a model port (provider/model.js: prepare(intent)
// → sent, stream(sent)), so tests script the model's replies.
//
// What ends the turn, and why each request is made, is the turn policy's
// (turnPolicy.js, its seam): the loop tells it what happened — the turn began,
// a reply came back, a round of calls ran, the model refused its tools — and
// does what its step says: ask (with the notes, the tool choice and the retry
// the step names), run the reply's calls, or end on the step's `done`. The loop
// keeps the stream, the calls and the messages.
//
// Events:
//   { type: 'round', round }                       a model request starts
//   { type: 'request', round, why, sent }          the request as sent (after the port's
//                                                  prepare), before its reply; why: 'first',
//                                                  'tools' (results of the round before),
//                                                  'reminder' (after an empty or looping
//                                                  reply), 'gate' (an answer sent back),
//                                                  'final' (the budget is used: tool_choice
//                                                  'none'), 'no_tools' (the retry without tools).
//                                                  A request the port refuses in prepare (the
//                                                  context guard, a spend cap) has none
//   { type: 'text', delta } / { type: 'reasoning', delta }
//   { type: 'tool_start', id, name, label, args, page, cards, interactive }
//   { type: 'tool_end', id, name, status, args, data, checks, error, ms, cached, ref, analysis, pageParams }
//                                                  ref: the result's reference (null on an error);
//                                                  analysis: the step is an analysis
//                                                  (verification isAnalysis; the footer counts these)
//   { type: 'tools_unsupported', error }           the port's ToolsUnsupportedError on the
//                                                  first request: retried without tools
//   { type: 'gate', id, what, missing, notice }    the evidence gate sent the answer back
//                                                  (a 'retract' withdraws its text)
//   { type: 'done', text, rawText, outcome, rounds, stoppedBy, cut, reason, audit, references, gates }
//                                                  the turn policy's, one shape on every path;
//                                                  stoppedBy: 'answer' | 'max_rounds' |
//                                                  'no_progress' | 'empty' (turnPolicy.js
//                                                  STOPPED_BY); an answer's fields are the
//                                                  verification's (verification/ judge):
//                                                  text: the answer with its {{references}} filled in
//                                                  (rawText: as the model wrote it)
//                                                  cut: the last reply stopped at the length limit
//                                                  reason: why an empty reply was empty, else null
//                                                  audit: findings on the answer (answer checks,
//                                                  unfilled references, missing evidence)
//                                                  outcome: computed from the checks (verification/verdict.js)
// Errors (HTTP, network) and aborts are thrown to the caller.
//
// Budgets: at most `maxRounds` model requests with tool calls and
// `maxCallsPerRound` calls in each; one call may take `toolTimeoutMs`; after
// `maxTurnMs` no new round of tools starts. The rounds, the turn's time and
// the rounds without progress end a turn: the policy holds them. The calls
// per round and the questions only refuse a call: the loop holds them, with
// the per-call time and the character budgets. Each result is fitted to
// `resultCharBudget` characters, and the turn's tool messages together to
// `turnCharBudget` — older results are compacted to their checks first — so a
// small local model's context window keeps the run context and the rules.
// A tool may carry its own `resultBudget` (a skill is longer than an analysis
// result) and `keep: true` (its message is never compacted and does not count
// toward the turn budget: the agent follows a skill for the whole question).
// An interactive tool (ask_user) waits for the user: it has no time limit, the
// wait does not count toward `maxTurnMs`, and a turn asks at most
// `maxQuestions` questions, one per step. Turn time is read from `clock.now()`
// (the turn's clock; Date.now by default) by the policy, which starts the turn
// when it is created; the loop reports each wait for the user to it.
//
// `clarifier` (clarify.js) lets the Workbench ask itself where the run data
// show that a choice changes the result: before the first request
// (`beforeTurn`: the answer joins the user's question for the model) and
// before an analysis runs in a window that cuts a shell (`beforeCall`: it runs
// in the window the user picks). These questions are ask_user steps
// (`app: true`) and count toward `maxQuestions`.
//
// Evidence, numbers and the verdict (THY-18) are the verification's
// (verification/: one per turn, the turn's `verification`; without one, a
// verification with no claims). The loop records each tool step with it,
// which names the step's `ref` (r1, r2, …: the model quotes values by it,
// {{r2.peak.significance_sigma}}) and keeps the data the message carries.
// The policy asks it to judge each answer (its `judge`): whether a send-back
// is allowed is the policy's (only while tools can still run), the send-back
// the verification's (the gate: once per claim, at most its maxSendBacks times
// a turn); otherwise the `done` — the answer with its references filled, its
// findings and the Outcome computed from the checks; the model's Outcome line
// only explains it.
//
// On a local model (robustness.js): a model request that sends nothing for
// `stallMs`, or runs past `requestMaxMs`, is cancelled and explained; a reply
// with neither an answer nor a call (a reasoning-only reply), or one stopped
// because it kept writing the same passage (`repeatsItself`), gets one
// reminder, asked a little warmer after a loop (the intent's `retry: 'warmer'`:
// the model port raises the temperature, as greedy decoding is what loops),
// before the turn ends on an explained empty reply; a call written
// as JSON text is run; and two rounds in a row that add nothing new (every
// call repeated or failed) end the analyses: the next request asks for the
// answer (`stoppedBy: 'no_progress'`), so repetition never spends the round
// budget. The loop finds the loop, the calls and whether a round added
// anything; the reminder, the warmer retry and the no-progress rule are the
// policy's.

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
    maxNoProgressRounds: 2
};

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

// The turn's clock when none is given (ms).
const DEFAULT_CLOCK = Object.freeze({ now: () => Date.now() });

export async function* runCopilotTurn({
    // The model port (provider/model.js), or a bare stream function.
    model,
    messages: initialMessages,
    tools = [],
    limits = {},
    signal,
    toolsEnabled = true,
    noToolsSystemPrompt = null,
    // Passed to every tool's run() beside its signal (the card-focus
    // controller for view tools).
    toolContext = null,
    // The Workbench's own questions (clarify.js buildClarifier).
    clarifier = null,
    // The turn's verification (verification/ createVerification): the
    // references, the evidence gate, the checks on the answer and the
    // Outcome. Null: one with no claims (no gate; references, checks and the
    // Outcome as ever).
    verification = null,
    // The turn's clock: its time limit and the user's waits.
    clock = DEFAULT_CLOCK
}) {
    const port = asModelPort(model);
    const now = () => clock.now();
    const verify = verification || createVerification({ offered: [] });
    const {
        maxRounds, maxCallsPerRound, toolTimeoutMs, maxTurnMs, resultCharBudget, turnCharBudget, maxQuestions,
        stallMs, requestMaxMs, maxNoProgressRounds
    } = { ...DEFAULT_LIMITS, ...limits };
    let questions = 0;
    const available = toolsEnabled ? usableTools(tools) : [];
    const byName = new Map(available.map((tool) => [tool.name, tool]));
    const requestTools = available.length ? toRequestTools(available) : null;
    const messages = [...initialMessages];
    const cache = new Map();
    let callCounter = 0;
    // What ends the turn (turnPolicy.js); the turn's time starts here.
    const policy = createTurnPolicy({
        limits: { maxRounds, maxTurnMs, maxNoProgressRounds }, clock, tools: Boolean(requestTools), judge: verify.judge
    });
    // Indices (into `messages`) of this turn's tool messages, oldest first.
    const toolMessages = [];

    // The user's question: the last user message as the turn starts.
    const questionIndex = messages.map((message) => message.role).lastIndexOf('user');
    // With the user's answers to the Workbench's questions appended: a window
    // the user chose there is the user's, as if written in the question.
    let questionText = questionIndex >= 0 && typeof messages[questionIndex].content === 'string' ? messages[questionIndex].content : '';
    // Asked before the first request, or before a call (calls run only while
    // tools are on).
    const canAsk = () => Boolean(clarifier && requestTools && typeof toolContext?.askUser === 'function' && questions < maxQuestions);
    let appQuestions = 0;
    // The window choices made in this turn (clarify.js), by (centre, partner, window).
    const clarifySession = { windows: new Map() };
    // One question the Workbench asks itself, shown as an ask_user step: the
    // user's reply and the time the user took (off the turn's time).
    async function* askFromApp(item) {
        appQuestions += 1;
        const id = `app_ask_${appQuestions}`;
        const args = {
            question: item.request.question,
            options: item.request.options,
            default: item.request.options[item.request.defaultIndex]?.label
        };
        yield { type: 'tool_start', id, name: ASK_USER, label: 'Question for you', args, page: null, cards: [], interactive: true, app: true };
        const started = now();
        const reply = await toolContext.askUser(item.request, { signal });
        const ms = now() - started;
        questions += 1;
        const { data, checks } = answerResult(item.request, reply);
        yield { type: 'tool_end', id, name: ASK_USER, status: 'ok', args, data, checks, error: null, ms, cached: false, pageParams: null };
        return { reply, ms };
    }

    // The time the Workbench's question before the first request waited.
    let startWaitedMs = 0;
    if (canAsk() && questionText) {
        const item = clarifier.beforeTurn(questionText);
        if (item) {
            const { reply, ms } = yield* askFromApp(item);
            startWaitedMs = ms;
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

    // The policy's step: the next request, the calls to run, or the end.
    let step = policy.next({ kind: 'start', waitedMs: startWaitedMs });
    for (;;) {
        if (step.next === 'end') {
            if (step.say) yield { type: 'text', delta: step.say };
            yield step.done;
            return;
        }
        // The notes for this request (a reminder, the last-request note).
        step.notes.forEach((note) => addNote(messages, note));
        yield { type: 'round', round: step.round };

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
            // The intent (a snapshot of the messages: the loop keeps
            // appending to its own list), prepared by the port into the
            // request as sent, which the caller sees before the reply. A
            // request the port refuses to send (the context guard, a spend
            // cap) throws here, so it has no `request` event. Preparing may
            // probe the server: it runs under the stall clock too.
            const intent = {
                messages: [...messages],
                // Tools are off exactly when the step names no tool choice.
                tools: step.toolChoice ? requestTools : undefined,
                toolChoice: step.toolChoice,
                retry: step.retry
            };
            const sent = await beforeStall(port.prepare(intent, { signal: request.signal }), { stallMs, cancel: () => request.abort() });
            yield { type: 'request', round: step.round, why: step.why, sent };
            const stream = port.stream(sent, { signal: request.signal });
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
            // The port says, by type, that the model cannot take tools (it
            // reads each provider's failure; provider/model.js). Only the
            // turn's first request falls back: later ones already carried
            // tools successfully, so their failure is a real one.
            if (step.toolChoice && step.round === 0 && !text && error instanceof ToolsUnsupportedError) {
                // The tool rules would only confuse a model that has none.
                if (noToolsSystemPrompt && messages[0]?.role === 'system') {
                    messages[0] = { ...messages[0], content: noToolsSystemPrompt };
                }
                yield { type: 'tools_unsupported', error: error.message };
                step = policy.next({ kind: 'no_tools' });
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
        // The last request's calls are never run. The port enforces tool_choice
        // 'none' where the provider has no such field (Ollama: no tools are
        // sent); this guards the providers that take the field but are not
        // verified to honour it (LM Studio, Gemini, Anthropic).
        const callsAllowed = step.toolChoice === 'auto';
        let calls = callsAllowed ? finalizeToolCalls(fragments, `${step.round}_`) : [];
        // A call written as JSON text (small local models do that) is run.
        if (!calls.length && callsAllowed && text.trim()) {
            const written = textToolCalls(text, [...byName.keys()]);
            if (written.length) {
                calls = written.map((call, index) => ({ id: `call_text_${step.round}_${index}`, name: call.name, arguments: call.arguments, extraContent: null }));
                text = '';
                yield { type: 'retract' };
            }
        }
        // The policy reads the reply: a reminder, a send-back, the end, or the
        // calls to run. A last request that came back as tool calls only (a
        // provider that does not honour tool_choice 'none'), or with nothing
        // after its reminder, ends on an answer that says so.
        const after = policy.next({ kind: 'reply', text, calls, looped, cut: finishReason === 'length' });
        if (after.next !== 'run') {
            if (after.sendBack) {
                // The evidence gate sent the answer back: it is withdrawn, and
                // the model reads it with the app's note.
                const { id, what, missing, notice, nudge } = after.sendBack;
                yield { type: 'gate', id, what, missing, notice, text };
                yield { type: 'retract' };
                messages.push({ role: 'assistant', content: assistantContent(text) });
                messages.push({ role: 'user', content: nudge });
            }
            step = after;
            continue;
        }

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
        // The time this round waited for the user's answers.
        let waitedMs = 0;

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
                    const { reply, ms } = yield* askFromApp(item);
                    waitedMs += ms;
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
            // Whether the tool really ran (or its identical earlier run is
            // repeated): not a skipped call, an unknown name or a refused question.
            let executed = false;
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
                executed = true;
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
                    executed = true;
                    note = 'identical call earlier in this question — the same result is repeated; change the arguments or answer';
                } else {
                    result = await executeToolCall(tool, runArguments, { signal, timeoutMs: toolTimeoutMs, context: toolContext });
                    executed = true;
                    if (key && !result.error) cache.set(key, result);
                }
            }
            if (appNote) note = note ? `${appNote}; ${note}` : appNote;
            if (answered || (!result.error && !cached)) progressed = true;

            if (!tool?.keep) toolMessages.push(messages.length);
            const budget = Number.isFinite(tool?.resultBudget) && tool.resultBudget > 0 ? tool.resultBudget : resultCharBudget;
            // The step, recorded with the verification: a result the model
            // can quote by reference gets its ref, kept with the data it was
            // sent; its facts are for the checks on the answer.
            const sent = verify.record({
                name: call.name, args: result.args, status: stepStatus(result), error: result.error || null,
                interactive: Boolean(tool?.interactive), reading: Boolean(tool?.keep), view: Boolean(tool?.view), executed,
                facts: result.facts
            }, (ref) => toolMessageContent({ ...result, note }, budget, ref));
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
                ref: sent.ref,
                analysis: sent.analysis,
                pageParams: result.error ? null : pageParamsFor(tool, result.args)
            };
        }
        // The next request carries this round's results.
        step = policy.next({ kind: 'round', progressed, waitedMs });
    }
}
