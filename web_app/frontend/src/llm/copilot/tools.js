// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { parseToolArguments, validateToolArguments } from './toolCalls';

// The Copilot's tool contract and how one call is executed. A tool is a plain
// object the host (or this module, for tools over its own props) supplies:
//
//   {
//     name: 'compute_bond_angles',           // /^[a-z][a-z0-9_]{0,63}$/
//     description: '…',                      // what it computes, for the model
//     parameters: { type: 'object', properties: {…}, required: [...] },
//     label: (args) => 'Bond angles Se–Ta–Se',   // step-card title (optional)
//     page: 'geometry',                      // page key for "Open page" (optional)
//     pageParams: (args) => ({ … }),         // what that page should load to show
//                                            // the same analysis (optional)
//     run: async (args, { signal }) => ({ data, checks })
//   }
//
// `data` is the compact JSON the model reads (never a raw engine payload);
// `checks` are deterministic verdicts on the full result, each
// `{ level: 'ok' | 'note' | 'warn' | 'fail', text }`. A thrown error becomes a
// failed step whose message the model reads, so it can correct its arguments.

export const TOOL_NAME_PATTERN = /^[a-z][a-z0-9_]{0,63}$/;
export const CHECK_LEVELS = ['ok', 'note', 'warn', 'fail'];

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

/** A step's status from its error and checks: failed > warn > ok. */
export const stepStatus = ({ error, checks }) => {
    if (error || checks?.some((check) => check.level === 'fail')) return 'failed';
    if (checks?.some((check) => check.level === 'warn')) return 'warn';
    return 'ok';
};

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

const timeoutError = (ms) => new Error(`the analysis did not finish within ${Math.round(ms / 1000)} s`);
const abortError = () => new DOMException('Aborted', 'AbortError');

/**
 * Parse, validate and run one call. Resolves (never rejects, except on abort
 * of the turn) to `{ args, data, checks, error, ms }`; `args` is the validated
 * object the tool ran with (or the raw text when it could not be parsed).
 * The tool's `run` receives its own `signal`, aborted when the turn is stopped
 * or the call times out, so an engine that can be cancelled (the Copilot's
 * worker, an HTTP request) stops working instead of running on unseen.
 */
export const executeToolCall = async (tool, rawArguments, { signal, timeoutMs }) => {
    const started = performance.now();
    const finish = (fields) => ({ ms: performance.now() - started, data: null, checks: [], error: null, ...fields });
    const parsed = parseToolArguments(rawArguments);
    if (!parsed.ok) return finish({ args: rawArguments, error: parsed.error });
    const validated = validateToolArguments(tool.parameters, parsed.value);
    const notes = validated.notes.map((text) => ({ level: 'note', text }));
    if (!validated.ok) {
        return finish({ args: parsed.value, error: validated.errors.join('; '), checks: notes });
    }
    if (signal?.aborted) throw abortError();

    const call = new AbortController();
    let timer = null;
    const onTurnAbort = () => call.abort();
    signal?.addEventListener('abort', onTurnAbort, { once: true });
    const limits = new Promise((_, reject) => {
        timer = setTimeout(() => {
            reject(timeoutError(timeoutMs));
            call.abort();
        }, timeoutMs);
        call.signal.addEventListener('abort', () => reject(signal?.aborted ? abortError() : timeoutError(timeoutMs)), { once: true });
    });
    try {
        const result = await Promise.race([
            Promise.resolve().then(() => tool.run(validated.value, { signal: call.signal })),
            limits
        ]);
        return finish({
            args: validated.value,
            data: result?.data ?? null,
            checks: [...notes, ...normalizeChecks(result?.checks)]
        });
    } catch (error) {
        if (signal?.aborted) throw abortError();
        return finish({ args: validated.value, error: error?.message || String(error), checks: notes });
    } finally {
        clearTimeout(timer);
        signal?.removeEventListener('abort', onTurnAbort);
        // The limit promise never settles once the tool has answered (its
        // timer is cleared); if it lost the race by rejecting, swallow that.
        limits.catch(() => {});
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
