// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { answerResult, replyChoice } from './askUser';
import { keyNumbers } from './robustness';
import { parseToolArguments, validateToolArguments } from './toolCalls';

// The Copilot's tool contract: what a tool is, how one call is executed, and
// how every call ends. A tool is a plain object the host (or this module, for
// tools over its own props) supplies:
//
//   {
//     name: 'compute_bond_angles',           // /^[a-z][a-z0-9_]{0,63}$/
//     description: '…',                      // what it computes, for the model
//     parameters: { type: 'object', properties: {…}, required: [...] },
//     dataShape: { sites: [{ ref: 'number', … }], … },   // the shape of `data`
//                                            // (dataShape.js): the paths the model
//                                            // quotes by reference. Internal, never
//                                            // sent; a tool that returns data declares it
//     label: (args) => 'Bond angles Se–Ta–Se',   // step-card title (optional)
//     page: 'geometry',                      // page key for "Open page" (optional)
//     pageParams: (args) => ({ … }),         // what that page should load to show
//                                            // the same analysis (optional)
//     cards: (args) => ['card-id'],          // the host cards this call is about:
//                                            // they breathe while it runs (optional)
//     view: true,                            // changes only what is on screen (a
//                                            // card's focus, a 3D camera): never
//                                            // answered from the call cache (optional)
//     interactive: true,                     // waits for the user (ask_user): never
//                                            // cached, no time limit (optional)
//     resultBudget: 9000,                    // its own result budget, characters (optional)
//     keep: true,                            // its result is never compacted later in
//                                            // the question (read_skill) (optional)
//     run: async (args, { signal, focus, askUser }) => ({ data, checks, facts, choice })
//   }
//
// `data` is the compact JSON the model reads (never a raw engine payload);
// `checks` are deterministic verdicts on the full result, each
// `{ level: 'ok' | 'note' | 'warn' | 'fail', text }`: a `fail` check is the
// tool's own refusal ("this does not answer as asked"). `facts` (optional, not
// sent to the model) are machine-readable findings the final answer must not
// contradict: one key per row of the verification's facts catalogue
// (verification/facts.js FACTS: each key, the tools that return it, its shape
// and the check that reads it). `choice` (an interactive tool's, optional:
// ask_user's) is the user's reply, `{ source: 'option' | 'other' |
// 'skipped', index }` (askUser.js replyChoice). A thrown error becomes an
// error result whose message the model reads, so it can correct its
// arguments. `focus` (when the chat has a page to point at) is the card-focus
// controller (focus/focus.js): view tools call focus.show(ids) and
// focus.onClear(undo); analysis tools leave it alone. `askUser(request, {
// signal })` (when the conversation can show a question card) is how an
// interactive tool asks the user and waits for the reply (askUser.js).
//
// How a call ended is decided here, never in the loop that forwards it. Every
// call ends as one ToolResult:
//
//   { resultStatus, args, data, checks, error, facts?, ms, cached, note,
//     choice?: { source, index, confirmed } }
//
//   resultStatus (RESULT_STATUSES, closed):
//     'ok'            it ran and answered;
//     'refused'       it ran and refused (no error, a `fail` check);
//     'error'         it could not run as asked: arguments that are not JSON or
//                     that the schema refuses, a thrown error, the time limit;
//     'unknown_tool'  no tool has that name;
//     'skipped'       the loop did not run it (SKIP_REASONS: the calls per round,
//                     one question at a time, the questions per request).
//   args     the validated arguments (the parsed ones when they failed the
//            schema, the raw text when they did not parse);
//   error    the message the model reads instead of data: non-empty exactly for
//            'error', 'unknown_tool' and 'skipped', which have no data and no
//            facts;
//   cached   repeated from an identical earlier call (only 'ok' or 'refused');
//   note     what the app adds for the model beside the result (the Workbench's
//            note, the repeat's); for a Workbench question, which sends no
//            message of its own, the text the model reads of it instead;
//   choice   a question the user answered: the reply's source and option index,
//            and `confirmed`, whether the window the user chose is their own
//            (clarify.js; null when the question is not about a window).
//
// The constructors (toolResult, skippedCall, unknownToolCall, repeatedCall,
// withNote, questionResult, and executeToolCall itself) check these
// invariants and throw a TypeError on a result that breaks one. The tool
// message the model reads depends only on the ref, the error, the note, the
// checks and the data, in that key order (toolMessage); the status and the
// choice are never sent. severityOf gives the step's status as the dock and
// the Outcome read it ('ok' | 'warn' | 'failed'), ran whether the tool ran.

export const TOOL_NAME_PATTERN = /^[a-z][a-z0-9_]{0,63}$/;
export const CHECK_LEVELS = ['ok', 'note', 'warn', 'fail'];

/** How a call ended (see the header): a closed set. */
export const RESULT_STATUSES = Object.freeze(['ok', 'refused', 'error', 'unknown_tool', 'skipped']);
/** The statuses whose result is an error the model reads, with no data. */
export const ERROR_STATUSES = Object.freeze(['error', 'unknown_tool', 'skipped']);
// Where a question's reply came from (askUser.js).
export const CHOICE_SOURCES = Object.freeze(['option', 'other', 'skipped']);

/** Why the loop skipped a call: each reason's message, given the limit or the questions asked. */
export const SKIP_REASONS = Object.freeze({
    call_limit: (n) => `skipped: at most ${n} tool calls per step — call it again in the next step if it is still needed`,
    one_question: () => 'one question at a time: ask the next one after this answer, if it is still needed',
    question_limit: (n) => `already asked ${n} question${n === 1 ? '' : 's'} for this request: `
        + 'proceed with your default and say in the answer what you assumed'
});
// The reasons whose message names a number.
const COUNTED_REASONS = Object.freeze(['call_limit', 'question_limit']);

/** The note on a call answered from the cache. */
export const REPEAT_NOTE = 'identical call earlier in this question — the same result is repeated; change the arguments or answer';

/** Tools whose shape is usable: a valid unique name, a description and a run(). */
export const usableTools = (tools) => {
    const seen = new Set();
    return (Array.isArray(tools) ? tools : []).filter((tool) => {
        if (!tool || typeof tool.run !== 'function' || typeof tool.description !== 'string') return false;
        if (!TOOL_NAME_PATTERN.test(tool.name || '') || seen.has(tool.name)) return false;
        seen.add(tool.name);
        return true;
    });
};

/** The `tools` request field (OpenAI function-calling format). */
export const toRequestTools = (tools) => tools.map((tool) => ({
    type: 'function',
    function: {
        name: tool.name,
        description: tool.description,
        parameters: tool.parameters || { type: 'object', properties: {} }
    }
}));

const normalizeChecks = (checks) => (Array.isArray(checks) ? checks : [])
    .filter((check) => check && typeof check.text === 'string' && check.text)
    .map((check) => ({ level: CHECK_LEVELS.includes(check.level) ? check.level : 'note', text: check.text }));

const isFail = (check) => check?.level === 'fail';

const invalid = (where, message) => {
    throw new TypeError(`${where}: ${message}`);
};

// A choice is { source, index, confirmed }: a known source, an integer index
// or null, and a confirmed that is true, false or null.
const checkChoice = (where, choice) => {
    if (!choice || typeof choice !== 'object' || Array.isArray(choice)) invalid(where, 'choice must be { source, index, confirmed }');
    if (!CHOICE_SOURCES.includes(choice.source)) invalid(where, `choice.source must be one of ${CHOICE_SOURCES.join(', ')}, got ${JSON.stringify(choice.source)}`);
    if (!(choice.index === null || Number.isInteger(choice.index))) invalid(where, 'choice.index must be an integer or null');
    if (![true, false, null].includes(choice.confirmed)) invalid(where, 'choice.confirmed must be true, false or null');
};

/**
 * A ToolResult (see the header), checked: a known status, an error exactly
 * when the status is 'error', 'unknown_tool' or 'skipped' (then no data and
 * no facts), a fail check exactly when it is 'refused' (or beside an error),
 * a repeat only of 'ok' or 'refused', a choice only on 'ok'. Frozen (shallow).
 * Every field but the status has a default, so `toolResult({ facts })` is a
 * step that ran and answered.
 */
export const toolResult = ({
    resultStatus = 'ok', args = {}, data = null, checks = [], error = null, facts, ms = 0, cached = false, note = null, choice
} = {}) => {
    const where = 'toolResult';
    if (!RESULT_STATUSES.includes(resultStatus)) invalid(where, `resultStatus must be one of ${RESULT_STATUSES.join(', ')}, got ${JSON.stringify(resultStatus)}`);
    const failed = ERROR_STATUSES.includes(resultStatus);
    if (failed && !(typeof error === 'string' && error)) invalid(where, `a '${resultStatus}' result names its error`);
    if (!failed && error !== null) invalid(where, `an '${resultStatus}' result has no error, got ${JSON.stringify(error)}`);
    if (!Array.isArray(checks)) invalid(where, 'checks must be an array');
    if (!failed && (resultStatus === 'refused') !== checks.some(isFail)) {
        invalid(where, resultStatus === 'refused' ? 'a \'refused\' result carries the tool\'s fail check' : 'an \'ok\' result has no fail check: it is \'refused\'');
    }
    if (failed && data !== null) invalid(where, `a '${resultStatus}' result has no data`);
    if (failed && facts !== undefined) invalid(where, `a '${resultStatus}' result has no facts`);
    if (cached && resultStatus !== 'ok' && resultStatus !== 'refused') invalid(where, 'only an \'ok\' or \'refused\' result is repeated from the cache');
    if (!(note === null || (typeof note === 'string' && note))) invalid(where, 'note must be a non-empty string or null');
    if (choice !== undefined) {
        if (resultStatus !== 'ok') invalid(where, `a '${resultStatus}' result has no choice: only a question the user answered`);
        checkChoice(where, choice);
    }
    return Object.freeze({
        resultStatus, args, data, checks, error,
        ...(facts !== undefined ? { facts } : {}),
        ms, cached: Boolean(cached), note,
        ...(choice !== undefined ? { choice: Object.freeze({ ...choice }) } : {})
    });
};

/** A call the loop did not run, for `reason` (SKIP_REASONS); `n` is the limit or the questions asked. */
export const skippedCall = (args, reason, n) => {
    if (!Object.hasOwn(SKIP_REASONS, reason)) invalid('skippedCall', `unknown reason ${JSON.stringify(reason)} (reasons: ${Object.keys(SKIP_REASONS).join(', ')})`);
    if (COUNTED_REASONS.includes(reason) && !(Number.isInteger(n) && n >= 0)) invalid('skippedCall', `'${reason}' names its number, got ${JSON.stringify(n)}`);
    return toolResult({ resultStatus: 'skipped', args, error: SKIP_REASONS[reason](n) });
};

/** A call to a tool the turn does not offer: `names` are the ones it does. */
export const unknownToolCall = (args, name, names) => toolResult({
    resultStatus: 'unknown_tool', args, error: `unknown tool "${name}"; available: ${(names || []).join(', ') || 'none'}`
});

/** An identical earlier call's result, repeated: no time, `cached`, and the repeat's note. */
export const repeatedCall = (cached) => {
    if (!cached || (cached.resultStatus !== 'ok' && cached.resultStatus !== 'refused')) {
        invalid('repeatedCall', `only an 'ok' or 'refused' result is repeated, got ${JSON.stringify(cached?.resultStatus)}`);
    }
    return toolResult({ ...cached, ms: 0, cached: true, note: REPEAT_NOTE });
};

/** `result` with the app's `note` before its own (the Workbench's, then the repeat's); unchanged without one. */
export const withNote = (result, note) => (note ? toolResult({ ...result, note: result.note ? `${note}; ${result.note}` : note }) : result);

/**
 * The result of a question the user answered (`reply` from the ask-user
 * port, `request` the card's): answerResult's data and checks, and the choice
 * with `confirmed` (clarify.js: whether a window question's answer keeps the
 * user's own window; null otherwise). `args` and `ms` are the step's.
 */
export const questionResult = (request, reply, { args = {}, ms = 0, confirmed = null } = {}) => {
    const { data, checks } = answerResult(request, reply);
    return toolResult({ resultStatus: 'ok', args, data, checks, ms, choice: { ...replyChoice(reply), confirmed } });
};

/** A step's severity as the dock and the Outcome read it, from its error and checks: failed > warn > ok. */
export const severityOf = ({ error, checks }) => {
    if (error || checks?.some(isFail)) return 'failed';
    if (checks?.some((check) => check.level === 'warn')) return 'warn';
    return 'ok';
};

/** Whether the tool ran (or its identical earlier run was repeated): not an unknown tool, not a skipped call. */
export const ran = (result) => RESULT_STATUSES.includes(result?.resultStatus)
    && result.resultStatus !== 'unknown_tool' && result.resultStatus !== 'skipped';

// Trim the largest array found anywhere in `value` to half its length,
// recording the cut as `<key>_omitted` beside it (or a trailing note for an
// array nested in an array). Returns false when no array is left to trim.
const trimLargestArray = (root) => {
    let best = null;
    const visit = (node, parent, key) => {
        if (Array.isArray(node)) {
            const size = JSON.stringify(node).length;
            if (node.length > 1 && (!best || size > best.size)) best = { node, parent, key, size };
            node.forEach((child, index) => visit(child, node, index));
        } else if (node && typeof node === 'object') {
            Object.entries(node).forEach(([childKey, child]) => visit(child, node, childKey));
        }
    };
    visit(root, null, null);
    if (!best) return false;
    const keep = Math.max(1, Math.floor(best.node.length / 2));
    const omitted = best.node.length - keep;
    const kept = best.node.slice(0, keep);
    if (best.parent && !Array.isArray(best.parent)) {
        best.parent[best.key] = kept;
        const omittedKey = `${best.key}_omitted`;
        best.parent[omittedKey] = (Number(best.parent[omittedKey]) || 0) + omitted;
    } else if (best.parent) {
        best.parent[best.key] = [...kept, `… ${omitted} more omitted`];
    } else {
        return false;
    }
    return true;
};

/**
 * Fit `value` within `budget` characters of JSON: halve the largest array
 * (recording `*_omitted` counts) until it fits, then cut the text hard as a
 * last resort, flagged `truncated: true`. Never silent. Returns a new value.
 */
export const fitToBudget = (value, budget) => {
    let json = JSON.stringify(value ?? null);
    if (json.length <= budget) return value ?? null;
    const copy = JSON.parse(json);
    if (copy && typeof copy === 'object' && !Array.isArray(copy)) {
        while (json.length > budget && trimLargestArray(copy)) {
            json = JSON.stringify(copy);
        }
        if (json.length <= budget) return copy;
    }
    return { truncated: true, text: json.slice(0, Math.max(0, budget - 40)) };
};

export const budgetJson = (value, budget) => JSON.stringify(fitToBudget(value, budget));

/**
 * The tool message the model reads for a result: its ref (with no error),
 * the error, the note, the checks and the data fitted to `budget` characters
 * less the rest, in that key order. Checks are never trimmed: they are the
 * verdict. Returns the message text and the data as fitted (null with an
 * error), which the step shows as the result sent to the model and the ref
 * resolves to. The status and the choice are never sent.
 */
export const toolMessage = ({ data, checks, error, note }, { budget, ref = null } = {}) => {
    if (!(Number.isFinite(budget) && budget > 0)) invalid('toolMessage', `budget must be a positive number of characters, got ${JSON.stringify(budget)}`);
    const message = {};
    if (ref && !error) message.ref = ref;
    if (error) message.error = error;
    if (note) message.note = note;
    if (checks?.length) message.checks = checks;
    const head = JSON.stringify(message).length;
    if (!error) message.data = fitToBudget(data, Math.max(400, budget - head));
    return { content: JSON.stringify(message), data: error ? null : message.data };
};

/**
 * An earlier tool message shortened to fit the turn's budget: its ref (which
 * still resolves: the app keeps the result it sent), its error, its checks
 * and the headline numbers of its data (keyNumbers). It reads its own JSON
 * back, so the bytes it keeps are the message's; a message already
 * compacted is returned as it is (compacting again would lose its numbers).
 */
export const compactToolMessage = (content) => {
    try {
        const { ref, checks, error, data, compacted } = JSON.parse(content);
        if (compacted) return content;
        const key = error ? null : keyNumbers(data);
        return JSON.stringify({
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

const timeoutError = (ms) => new Error(`the analysis did not finish within ${Math.round(ms / 1000)} s`);
const abortError = () => new DOMException('Aborted', 'AbortError');
// A thrown error with no text of its own still names an error.
const NO_MESSAGE = 'the tool failed without saying why';

// An interactive tool's reply, as the result's choice (confirmed: null, the
// loop's to set for the Workbench's own questions); dropped when malformed.
const readChoice = (tool, choice) => (tool?.interactive && choice && CHOICE_SOURCES.includes(choice.source)
    ? { source: choice.source, index: Number.isInteger(choice.index) ? choice.index : null, confirmed: null }
    : undefined);

/**
 * Parse, validate and run one call. Resolves (never rejects, except on abort
 * of the turn) to a ToolResult: 'ok', 'refused' (the tool answered with a
 * fail check) or 'error' (arguments that do not parse or that the schema
 * refuses, a thrown error, the time limit); `args` is the validated object
 * the tool ran with (or the raw text when it could not be parsed). The
 * tool's `run` receives its own `signal`, aborted when the turn is stopped
 * or the call times out, so an engine that can be cancelled (the Copilot's
 * worker, an HTTP request) stops working instead of running on unseen.
 * Without a finite positive `timeoutMs` (an interactive tool waiting for the
 * user) only the turn's own signal stops it.
 */
export const executeToolCall = async (tool, rawArguments, { signal, timeoutMs, context = null }) => {
    const started = performance.now();
    const elapsed = () => performance.now() - started;
    const failed = (args, error, checks = []) => toolResult({ resultStatus: 'error', args, error: error || NO_MESSAGE, checks, ms: elapsed() });
    const parsed = parseToolArguments(rawArguments);
    if (!parsed.ok) return failed(rawArguments, parsed.error);
    const validated = validateToolArguments(tool.parameters, parsed.value);
    const notes = validated.notes.map((text) => ({ level: 'note', text }));
    if (!validated.ok) return failed(parsed.value, validated.errors.join('; '), notes);
    if (signal?.aborted) throw abortError();

    const call = new AbortController();
    let timer = null;
    const onTurnAbort = () => call.abort();
    signal?.addEventListener('abort', onTurnAbort, { once: true });
    const limits = new Promise((_, reject) => {
        if (Number.isFinite(timeoutMs) && timeoutMs > 0) {
            timer = setTimeout(() => {
                reject(timeoutError(timeoutMs));
                call.abort();
            }, timeoutMs);
        }
        call.signal.addEventListener('abort', () => reject(signal?.aborted ? abortError() : timeoutError(timeoutMs)), { once: true });
    });
    let answer;
    try {
        answer = await Promise.race([
            Promise.resolve().then(() => tool.run(validated.value, { ...(context || {}), signal: call.signal })),
            limits
        ]);
    } catch (error) {
        if (signal?.aborted) throw abortError();
        return failed(validated.value, error?.message || String(error), notes);
    } finally {
        clearTimeout(timer);
        signal?.removeEventListener('abort', onTurnAbort);
        // The limit promise never settles once the tool has answered (its
        // timer is cleared); if it lost the race by rejecting, swallow that.
        limits.catch(() => {});
    }
    // The tool answered: 'refused' when one of its checks fails, else 'ok'.
    const checks = [...notes, ...normalizeChecks(answer?.checks)];
    const resultStatus = checks.some(isFail) ? 'refused' : 'ok';
    const choice = resultStatus === 'ok' ? readChoice(tool, answer?.choice) : undefined;
    return toolResult({
        resultStatus,
        args: validated.value,
        data: answer?.data ?? null,
        checks,
        ...(answer?.facts && typeof answer.facts === 'object' ? { facts: answer.facts } : {}),
        ms: elapsed(),
        ...(choice ? { choice } : {})
    });
};

/** The host cards a call is about (tool.cards(args)), [] when none or on error. */
export const toolCards = (tool, args) => {
    if (!tool || typeof tool.cards !== 'function' || !args || typeof args !== 'object') return [];
    try {
        const ids = tool.cards(args);
        return Array.isArray(ids) ? ids.filter((id) => typeof id === 'string' && id) : [];
    } catch {
        return [];
    }
};

/** A step-card title: the tool's own label, else its name. */
export const toolLabel = (tool, args) => {
    if (tool && typeof tool.label === 'function' && args && typeof args === 'object') {
        try {
            const text = tool.label(args);
            if (typeof text === 'string' && text) return text;
        } catch {
            // A label is cosmetic; fall back to the name.
        }
    }
    return tool?.name || 'tool';
};
