// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Copilot never sends a request its model would cut. A local server runs
// a model with a fixed context window; a prompt longer than the window is not
// refused but shortened (Ollama keeps the system message and the newest ones
// and drops the rest; a llama.cpp runner also shifts out the prompt's start
// when the reply fills the window), so the model answers without the run
// context or the rules and nobody is told. The Copilot's prompt is long — the
// rules, ~2,000 tokens of run context, the tool schemas, then the analyses —
// so a 4,096-token default (Ollama's on a machine without a usable GPU, older
// releases everywhere) cuts it on the first request.
//
// Before each request this module estimates the prompt's size and:
//   - Ollama: sends the request through the native API (ollama.js) with a
//     `num_ctx` that holds the prompt and room for the reply, and a
//     `num_predict` that keeps the reply inside it. It reuses the context the
//     model is loaded with when that is large enough, else asks for a fitting
//     one; it refuses the turn when even the model's own maximum is too small.
//   - LM Studio: the context is fixed when the model is loaded; when the loaded
//     model's window is too small it refuses the turn and says what to change.
//   - Anything else (cloud, other servers) is sent unchanged (client.js).

import { streamChat } from './client';
import { serverOrigin, streamOllamaChat } from './ollama';
import { isLocalUrl } from './presets';

// A deliberately low characters-per-token ratio, so the estimate errs high:
// the prompt is English, JSON and numbers (~3–4 characters per token).
export const CHARS_PER_TOKEN = 2.5;
// Template tokens around each message.
export const MESSAGE_OVERHEAD_TOKENS = 8;
// Room for the reply (and a reasoning model's thinking) beyond the prompt.
export const REPLY_RESERVE_TOKENS = 8192;
// LM Studio cannot be told to stop a reply early: it needs at least this much
// room beyond the prompt.
export const LMSTUDIO_REPLY_MIN_TOKENS = 4096;
// A context asked for is a multiple of this, and at least MIN_CONTEXT, so a
// growing conversation rarely changes it (each change reloads the model).
export const CONTEXT_STEP = 8192;
export const MIN_CONTEXT = 16384;

const OLLAMA_PORT = '11434';

export class ContextLimitError extends Error {
    constructor(message, { needed = null, available = null } = {}) {
        super(message);
        this.name = 'ContextLimitError';
        this.needed = needed;
        this.available = available;
    }
}

const textLength = (value) => (typeof value === 'string' ? value.length : (value == null ? 0 : JSON.stringify(value).length));

/** A conservative token count for the request: messages, tool calls, tool schemas. */
export const estimatePromptTokens = ({ messages = [], tools = null } = {}) => {
    let chars = 0;
    for (const message of messages || []) {
        chars += textLength(message.content);
        for (const call of message.tool_calls || []) {
            chars += textLength(call.function?.name) + textLength(call.function?.arguments);
        }
    }
    if (Array.isArray(tools) && tools.length) chars += JSON.stringify(tools).length;
    return Math.ceil(chars / CHARS_PER_TOKEN) + MESSAGE_OVERHEAD_TOKENS * (messages?.length || 0);
};

const formatTokens = (count) => Math.round(count).toLocaleString('en-US');

// Only a server on this machine, or on Ollama's port, is probed: a cloud
// provider is never sent a request it was not meant to get.
const probeable = (baseUrl) => {
    try {
        return isLocalUrl(baseUrl) || new URL(baseUrl).port === OLLAMA_PORT;
    } catch {
        return false;
    }
};

const servers = new Map();
const maxContexts = new Map();

/** Forget what each server and model is (tests; a server can be replaced on its port). */
export const forgetServers = () => {
    servers.clear();
    maxContexts.clear();
};

const getJson = async (url, init) => {
    const response = await fetch(url, init);
    if (!response.ok) return null;
    try {
        return await response.json();
    } catch {
        return null;
    }
};

/**
 * What kind of server a base URL points at: `{ kind: 'ollama', version }`
 * (GET /api/version answers a version), `{ kind: 'lmstudio' }` (its REST API
 * lists models with their load state), else `{ kind: 'openai' }`. Remembered
 * per base URL; an unreachable server is asked again next time.
 */
export const detectServer = async (baseUrl, { signal } = {}) => {
    if (!probeable(baseUrl)) return { kind: 'openai' };
    if (servers.has(baseUrl)) return servers.get(baseUrl);
    const origin = serverOrigin(baseUrl);
    let result = { kind: 'openai' };
    try {
        const version = await getJson(`${origin}/api/version`, { signal });
        if (typeof version?.version === 'string') {
            result = { kind: 'ollama', version: version.version };
        } else {
            const models = await getJson(`${origin}/api/v0/models`, { signal });
            if (Array.isArray(models?.data) && models.data.some((entry) => entry && 'state' in entry)) result = { kind: 'lmstudio' };
        }
    } catch (error) {
        if (error?.name === 'AbortError') throw error;
        return result;
    }
    servers.set(baseUrl, result);
    return result;
};

const withTag = (name) => (typeof name === 'string' && name && !name.includes(':') ? `${name}:latest` : name);

/** The context an Ollama model is loaded with now (GET /api/ps), or null. */
export const ollamaLoadedContext = async (baseUrl, model, { signal } = {}) => {
    try {
        const payload = await getJson(`${serverOrigin(baseUrl)}/api/ps`, { signal });
        const entry = (payload?.models || []).find((item) => withTag(item.name) === withTag(model) || withTag(item.model) === withTag(model));
        return Number.isFinite(entry?.context_length) ? entry.context_length : null;
    } catch (error) {
        if (error?.name === 'AbortError') throw error;
        return null;
    }
};

/** The longest context the model was trained for (POST /api/show), or null. */
export const ollamaMaxContext = async (baseUrl, model, { signal } = {}) => {
    const key = `${baseUrl}|${model}`;
    if (maxContexts.has(key)) return maxContexts.get(key);
    try {
        const payload = await getJson(`${serverOrigin(baseUrl)}/api/show`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ model }),
            signal
        });
        const entry = Object.entries(payload?.model_info || {}).find(([name, value]) => name.endsWith('.context_length') && Number.isFinite(value));
        const max = entry ? entry[1] : null;
        if (payload) maxContexts.set(key, max);
        return max;
    } catch (error) {
        if (error?.name === 'AbortError') throw error;
        return null;
    }
};

/**
 * The context for one request: `{ numCtx, numPredict, needed }`, or
 * `{ refuse: true, needed, max }` when the model cannot hold it. Prefers the
 * context the model is loaded with, then the one this conversation asked for
 * before, then a fitting multiple of CONTEXT_STEP — so it changes (and the
 * model reloads) only when the request outgrows it.
 */
export const planContext = ({ estimate, loaded = null, max = null, previous = null }) => {
    const needed = estimate + REPLY_RESERVE_TOKENS;
    if (Number.isFinite(max) && needed > max) return { refuse: true, needed, max };
    let numCtx;
    if (Number.isFinite(loaded) && loaded >= needed) numCtx = loaded;
    else if (Number.isFinite(previous) && previous >= needed) numCtx = previous;
    else numCtx = Math.max(MIN_CONTEXT, Math.ceil(needed / CONTEXT_STEP) * CONTEXT_STEP);
    if (Number.isFinite(max)) numCtx = Math.min(numCtx, max);
    return { numCtx, numPredict: numCtx - estimate, needed };
};

/** LM Studio's loaded context for the model (GET /api/v0/models), or null when not loaded. */
export const lmStudioLoadedContext = async (baseUrl, model, { signal } = {}) => {
    try {
        const payload = await getJson(`${serverOrigin(baseUrl)}/api/v0/models`, { signal });
        const entry = (payload?.data || []).find((item) => item?.id === model);
        return entry?.state === 'loaded' && Number.isFinite(entry.loaded_context_length) ? entry.loaded_context_length : null;
    } catch (error) {
        if (error?.name === 'AbortError') throw error;
        return null;
    }
};

const WHAT_IT_HOLDS = 'the Copilot\'s instructions, the run context, the analyses so far and room for the reply';

/**
 * The Copilot's model call: one streamed chat that the model receives whole.
 * Yields client.js fragments. `session` (an object the conversation keeps)
 * remembers the context asked for per model. Throws ContextLimitError instead
 * of sending a request the model would cut.
 */
export async function* streamCopilotChat({
    baseUrl, model, temperature, apiKey, messages, tools, toolChoice, signal, session = null
}) {
    const server = await detectServer(baseUrl, { signal });
    const estimate = server.kind === 'openai' ? null : estimatePromptTokens({ messages, tools });
    if (server.kind === 'ollama') {
        const [loaded, max] = await Promise.all([
            ollamaLoadedContext(baseUrl, model, { signal }),
            ollamaMaxContext(baseUrl, model, { signal })
        ]);
        const key = `${baseUrl}|${model}`;
        const plan = planContext({ estimate, loaded, max, previous: session?.[key] ?? null });
        if (plan.refuse) {
            throw new ContextLimitError(
                `This request needs about ${formatTokens(plan.needed)} tokens of context (${WHAT_IT_HOLDS}), but ${model} `
                + `holds at most ${formatTokens(plan.max)}. Start a new conversation, turn the analyses off, or pick a model `
                + 'with a longer context. Nothing was sent: the model would have answered without part of its instructions.',
                { needed: plan.needed, available: plan.max }
            );
        }
        if (session) session[key] = plan.numCtx;
        yield* streamOllamaChat({
            baseUrl, model, messages, tools, temperature, numCtx: plan.numCtx, numPredict: plan.numPredict, signal
        });
        return;
    }
    if (server.kind === 'lmstudio') {
        const loaded = await lmStudioLoadedContext(baseUrl, model, { signal });
        const needed = estimate + LMSTUDIO_REPLY_MIN_TOKENS;
        if (Number.isFinite(loaded) && loaded < needed) {
            const suggest = Math.max(MIN_CONTEXT, Math.ceil((estimate + REPLY_RESERVE_TOKENS) / CONTEXT_STEP) * CONTEXT_STEP);
            throw new ContextLimitError(
                `LM Studio runs ${model} with a ${formatTokens(loaded)}-token context, but this request needs about `
                + `${formatTokens(needed)} (${WHAT_IT_HOLDS}). Reload the model in LM Studio with a context length of at least `
                + `${formatTokens(suggest)}. Nothing was sent: the model would have answered without part of its instructions.`,
                { needed, available: loaded }
            );
        }
    }
    yield* streamChat({ baseUrl, model, temperature, apiKey, messages, tools, toolChoice, signal });
}
