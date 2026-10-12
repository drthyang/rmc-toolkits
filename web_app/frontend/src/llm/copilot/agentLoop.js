// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { toolBudgetNote } from '../prompts/copilot';
import { auditAnswer, lowerOutcome } from './answerAudit';
import { ASK_USER, MAX_QUESTIONS, answerResult } from './askUser';
import { resolveClarification } from './clarify';
import { parseOutcome } from './outcome';
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
//   { type: 'done', text, outcome, rounds, stoppedBy: 'answer' | 'max_rounds' | 'empty', cut, audit }
//                                                  cut: the last reply stopped at the length limit
//                                                  audit: findings on the answer (answerAudit.js);
//                                                  any finding lowers the outcome to 'not verified'
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

export const DEFAULT_LIMITS = {
    maxRounds: 6,
    maxCallsPerRound: 4,
    toolTimeoutMs: 120000,
    maxTurnMs: 600000,
    resultCharBudget: 3000,
    turnCharBudget: 12000,
    maxQuestions: MAX_QUESTIONS
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
        const { checks, error } = JSON.parse(content);
        return JSON.stringify({
            compacted: 'this earlier result was shortened to fit the context; its checks are kept',
            ...(error ? { error } : {}),
            ...(checks ? { checks } : {})
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

// The tool message the model reads: the result (fitted to the budget) and
// the checks, or the error. Checks are never trimmed: they are the verdict.
const toolMessageContent = ({ data, checks, error, note }, budget) => {
    const message = {};
    if (error) message.error = error;
    if (note) message.note = note;
    if (checks?.length) message.checks = checks;
    const head = JSON.stringify(message).length;
    if (!error) message.data = fitToBudget(data, Math.max(400, budget - head));
    return JSON.stringify(message);
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
    clarifier = null
}) {
    const {
        maxRounds, maxCallsPerRound, toolTimeoutMs, maxTurnMs, resultCharBudget, turnCharBudget, maxQuestions
    } = { ...DEFAULT_LIMITS, ...limits };
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
        const finalRound = useTools && (rounds >= maxRounds || outOfTime);
        if (finalRound) messages.push(toolBudgetNote(rounds));
        yield { type: 'round', round: rounds };

        let text = '';
        let finishReason = null;
        const fragments = [];
        try {
            const stream = callModel({
                messages,
                tools: useTools ? requestTools : undefined,
                toolChoice: useTools ? (finalRound ? 'none' : 'auto') : undefined,
                signal
            });
            for await (const chunk of stream) {
                if (chunk.content) {
                    text += chunk.content;
                    yield { type: 'text', delta: chunk.content };
                }
                if (chunk.reasoning) yield { type: 'reasoning', delta: chunk.reasoning };
                if (chunk.toolCall) accumulateToolCall(fragments, chunk.toolCall);
                if (chunk.finishReason) finishReason = chunk.finishReason;
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
        }

        const calls = useTools && !finalRound ? finalizeToolCalls(fragments, `${rounds}_`) : [];
        const cut = finishReason === 'length';
        if (!calls.length) {
            if (!text.trim()) {
                // Ollama ignores tool_choice 'none': a last round can come back
                // as tool calls only. Say so instead of ending on nothing.
                if (finalRound) {
                    yield { type: 'text', delta: EMPTY_FINAL };
                    yield {
                        type: 'done',
                        text: EMPTY_FINAL,
                        outcome: { status: 'not achieved', note: 'no answer was written within the analysis budget', body: EMPTY_FINAL },
                        rounds,
                        stoppedBy: 'max_rounds',
                        cut,
                        audit: []
                    };
                } else {
                    yield { type: 'done', text: '', outcome: null, rounds, stoppedBy: 'empty', cut, audit: [] };
                }
                return;
            }
            // A computed contradiction outranks the model's own verdict.
            const findings = audit(text);
            const outcome = lowerOutcome(parseOutcome(text), findings, text);
            yield { type: 'done', text, outcome, rounds, stoppedBy: finalRound ? 'max_rounds' : 'answer', cut, audit: findings };
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

        for (const [position, call] of calls.entries()) {
            callCounter += 1;
            const stepId = `${call.id}#${callCounter}`;
            const tool = byName.get(call.name);
            let parsed = parseToolArguments(call.arguments);
            let runArguments = call.arguments;
            let appNote = null;
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

            if (result.facts && !result.error) facts.push(result.facts);
            if (!tool?.keep) toolMessages.push(messages.length);
            const budget = Number.isFinite(tool?.resultBudget) && tool.resultBudget > 0 ? tool.resultBudget : resultCharBudget;
            messages.push({
                role: 'tool',
                tool_call_id: call.id,
                content: toolMessageContent({ ...result, note }, budget)
            });
            keepWithinTurnBudget(roundStart);
            yield {
                type: 'tool_end',
                id: stepId,
                name: call.name,
                status: stepStatus(result),
                args: result.args,
                data: result.data,
                checks: result.checks,
                error: result.error,
                ms: result.ms,
                cached,
                pageParams: result.error ? null : pageParamsFor(tool, result.args)
            };
        }
    }
}
