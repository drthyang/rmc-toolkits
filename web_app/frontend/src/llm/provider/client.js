// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Minimal OpenAI-compatible chat client for local LLM servers (Ollama's /v1,
// LM Studio, or anything speaking the same dialect). Hand-rolled on fetch — no
// SDK, no SSE library — so the whole provider surface stays readable in one
// file and adds zero runtime dependencies.

import { isAnthropicUrl } from './presets';

const trimBase = (baseUrl) => (baseUrl || '').replace(/\/+$/, '');

// Cloud providers (OpenAI, Gemini, …) authenticate with a Bearer token; local
// servers need none, so the header is only added when a key is present.
const authHeaders = (apiKey) => (apiKey ? { Authorization: `Bearer ${apiKey}` } : {});

// Anthropic's API blocks browser (CORS) requests unless the caller opts in with
// `anthropic-dangerous-direct-browser-access`, and its endpoints expect an API
// version. Its OpenAI-compatible /chat/completions still authenticates with a
// Bearer token, so only these two extra headers are Anthropic-specific (the
// host test is presets.js isAnthropicUrl; which models take no temperature is
// sampling.js).

// Headers for one request: Bearer auth plus, for Anthropic, its browser-access
// and version headers. Keyed off the base URL so no per-call plumbing is needed.
const requestHeaders = (baseUrl, apiKey) => ({
    ...authHeaders(apiKey),
    ...(isAnthropicUrl(baseUrl)
        ? {
            'anthropic-version': '2023-06-01',
            'anthropic-dangerous-direct-browser-access': 'true'
        }
        : {})
});

// A failed fetch to localhost surfaces as a bare TypeError both when the server
// is not running and when the browser blocked the response for CORS, so the
// hint has to name both causes — the user cannot tell them apart from the page.
// The fixes live in the setup guide ("How to connect a model"), which carries
// the copy-ready OLLAMA_ORIGINS command with this page's origin, the LM Studio
// CORS step and the Safari note; the settings drawer opens it on this failure.
const unreachableHint = (baseUrl) => (
    `Could not reach ${trimBase(baseUrl)}: the server is not running, or CORS blocks this page. `
    + 'See “How to connect a model”.'
);

const readErrorBody = async (response) => {
    try {
        return await response.json();
    } catch {
        // Non-JSON error bodies are fine; the status code is enough.
        return null;
    }
};

const errorMessage = (status, payload) => {
    const detail = payload?.error?.message || payload?.error || '';
    return `HTTP ${status}${detail ? `: ${detail}` : ''}`;
};

export const describeHttpError = async (response) => errorMessage(response.status, await readErrorBody(response));

/**
 * An HTTP failure as an Error: the message (`HTTP <status>: <detail>`), its
 * `status`, and `body`, the parsed error payload (null when it was not JSON),
 * so the model port can tell a failure's kind by its shape
 * (provider/model.js toolsUnsupportedFailure) rather than by its text.
 */
export const httpError = async (response) => {
    const body = await readErrorBody(response);
    const error = new Error(errorMessage(response.status, body));
    error.status = response.status;
    error.body = body;
    return error;
};

const httpHint = (status) => {
    if (status === 401) return 'Got 401 — the provider rejected the API key. Check that your key is correct and active.';
    if (status === 404) return 'Got 404 — check that the base URL ends in /v1 (e.g. http://localhost:11434/v1).';
    if (status === 403) return 'Got 403 — the server rejected this origin or key. Check its CORS/allowed-origins or key permissions.';
    if (status === 429) return 'Got 429 — the provider is rate-limiting or your quota is exhausted.';
    return null;
};

export const listModels = async (baseUrl, { signal, apiKey } = {}) => {
    const response = await fetch(`${trimBase(baseUrl)}/models`, { signal, headers: requestHeaders(baseUrl, apiKey) });
    if (!response.ok) throw await httpError(response);
    const payload = await response.json();
    return (payload.data || []).map((entry) => entry.id).filter(Boolean);
};

// Probe the server and translate failures into actionable setup hints.
// Returns { ok, models, error, hint } and never throws; a failure also says
// whether the server was unreachable (network / CORS) rather than answering.
export const checkConnection = async (baseUrl, { signal, apiKey } = {}) => {
    try {
        const models = await listModels(baseUrl, { signal, apiKey });
        return { ok: true, models, error: null, hint: null };
    } catch (error) {
        if (error?.name === 'AbortError') throw error;
        const isNetworkError = error instanceof TypeError;
        return {
            ok: false,
            models: [],
            error: error.message || 'Connection failed',
            hint: isNetworkError ? unreachableHint(baseUrl) : httpHint(error.status),
            unreachable: isNetworkError
        };
    }
};

// `tools` / `toolChoice` join the body only when given, so plain chat and
// watchdog requests are byte-identical to servers that predate tool calling.
// The body is built from what it is given: a null temperature is left out
// (the model's own default sampling; sampling.js says which models take
// none, and the model port and the watchdog send them null), and
// `max_tokens` joins only when `maxTokens` is a number (never by default).
// An HTTP failure carries its status and its parsed body (httpError), which
// lets the model port tell "this model cannot call tools" from a dead server.
const postChat = async ({ baseUrl, model, messages, temperature, maxTokens, stream, signal, apiKey, tools, toolChoice }) => {
    const body = { model, messages, ...(temperature === null ? {} : { temperature }), stream };
    if (Number.isFinite(maxTokens)) body.max_tokens = maxTokens;
    if (Array.isArray(tools) && tools.length) {
        body.tools = tools;
        if (toolChoice) body.tool_choice = toolChoice;
    }
    const response = await fetch(`${trimBase(baseUrl)}/chat/completions`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json', ...requestHeaders(baseUrl, apiKey) },
        body: JSON.stringify(body),
        signal
    });
    if (!response.ok) throw await httpError(response);
    return response;
};

// Stream a chat completion, yielding `{ content, reasoning }` deltas as they
// arrive (each chunk carries one or the other). The SSE body is `data: {json}`
// lines terminated by `data: [DONE]`; chunks can split mid-line, so incomplete
// tail lines are buffered across reads. Reasoning models (e.g. qwen3 via
// Ollama) stream their chain-of-thought in a separate `reasoning` field before
// the answer arrives in `content` — surfaced so the UI can show a thinking
// indicator instead of an empty bubble.
//
// With `tools`, a model that wants one streams `delta.tool_calls` fragments:
// each is yielded as `{ toolCall: { index, id, name, arguments } }` (any part
// may be missing; copilot/toolCalls.js assembles them; `extraContent` carries
// a call's provider-specific `extra_content`, which Gemini requires back), and the choice's
// `finish_reason`, when a chunk carries one, as `{ finishReason }`.
//
// `temperature` null sends none (an absent one is 0.2); `maxTokens` caps the
// reply as `max_tokens`, only when given.
export async function* streamChat({ baseUrl, model, messages, temperature = 0.2, maxTokens, signal, apiKey, tools, toolChoice }) {
    const response = await postChat({ baseUrl, model, messages, temperature, maxTokens, stream: true, signal, apiKey, tools, toolChoice });
    if (!response.body) throw new Error('The server returned no response body to stream');
    const reader = response.body.getReader();
    const decoder = new TextDecoder();
    let buffer = '';
    try {
        for (;;) {
            const { done, value } = await reader.read();
            if (done) break;
            buffer += decoder.decode(value, { stream: true });
            const lines = buffer.split('\n');
            buffer = lines.pop() ?? '';
            for (const line of lines) {
                if (!line.startsWith('data:')) continue;
                const data = line.slice(5).trim();
                if (!data) continue;
                if (data === '[DONE]') return;
                let parsed;
                try {
                    parsed = JSON.parse(data);
                } catch {
                    continue;
                }
                // A failure after the 200 arrives as an event of its own,
                // `{ error: … }` (LM Studio on a backend crash such as a GPU
                // out-of-memory, Ollama, OpenAI): surface it, or the turn
                // would end as an empty reply with no reason given.
                if (parsed?.error) {
                    const detail = typeof parsed.error === 'string'
                        ? parsed.error
                        : (parsed.error.message || JSON.stringify(parsed.error));
                    throw new Error(`The server stopped the reply: ${detail}`);
                }
                const choice = parsed.choices?.[0];
                const delta = choice?.delta;
                if (delta) {
                    if (delta.content) yield { content: delta.content };
                    // OpenAI-compatible servers name this `reasoning` (Ollama) or
                    // `reasoning_content` (DeepSeek and others).
                    const reasoning = delta.reasoning ?? delta.reasoning_content;
                    if (reasoning) yield { reasoning };
                    if (Array.isArray(delta.tool_calls)) {
                        for (const call of delta.tool_calls) {
                            if (!call || typeof call !== 'object') continue;
                            yield {
                                toolCall: {
                                    index: Number.isInteger(call.index) ? call.index : null,
                                    id: typeof call.id === 'string' && call.id ? call.id : null,
                                    name: typeof call.function?.name === 'string' ? call.function.name : null,
                                    arguments: typeof call.function?.arguments === 'string'
                                        ? call.function.arguments
                                        // A few servers send already-parsed arguments.
                                        : (call.function?.arguments && typeof call.function.arguments === 'object'
                                            ? JSON.stringify(call.function.arguments)
                                            : null),
                                    // Provider data that must travel back with the call
                                    // (Gemini's thought_signature lives here).
                                    extraContent: call.extra_content && typeof call.extra_content === 'object'
                                        ? call.extra_content
                                        : null
                                }
                            };
                        }
                    }
                }
                if (choice?.finish_reason) yield { finishReason: choice.finish_reason };
            }
        }
    } finally {
        reader.cancel().catch(() => {});
    }
}

// Non-streaming completion, used where the whole reply is parsed at once
// (watchdog status lines, report narrative fallback).
export const completeChat = async ({ baseUrl, model, messages, temperature = 0, signal, apiKey }) => {
    const response = await postChat({ baseUrl, model, messages, temperature, stream: false, signal, apiKey });
    const payload = await response.json();
    const content = payload.choices?.[0]?.message?.content;
    if (typeof content !== 'string') throw new Error('The model returned no message content');
    return content;
};
