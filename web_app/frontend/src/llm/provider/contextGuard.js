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
// Before each request this module estimates the prompt's size and decides
// (planCopilotChat, which the model port calls in its `prepare`, so the
// decision is in the request as sent, or refuses it before anything is
// streamed; streamCopilotChat then sends under that plan):
//   - Ollama: sends the request through the native API (ollama.js) with a
//     `num_ctx` that holds the prompt and room for the reply, and a
//     `num_predict` that keeps the reply inside it. It reuses the context the
//     model is loaded with when that is large enough, else asks for a fitting
//     one; it refuses the turn when even the model's own maximum is too small.
//   - LM Studio: the context is fixed when the model is loaded; when the loaded
//     model's window is too small it refuses the turn and says what to change.
//     A model not loaded yet is loaded by the request itself (JIT) with
//     LM Studio's own default context: its window is checked once it is
//     loaded, before any of the reply is shown.
//   - Anything else (cloud, other servers) is sent unchanged (client.js).

import { streamChat } from './client';
import { serverOrigin, streamOllamaChat } from './ollama';
import { isLocalUrl } from './presets';

// A deliberately low characters-per-token ratio, so the estimate errs high:
// the prompt is English, JSON and numbers (~3–4 characters per token).
export const CHARS_PER_TOKEN = 2.5;
// Tool results and call arguments are JSON dense with numbers, which
// tokenizes far finer: 1.64 characters per token measured on gemma4:26b
// (the 2026-10-12 review). Counted at this ratio, so they err high too.
export const JSON_CHARS_PER_TOKEN = 1.4;
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
// Ollama reports how many prompt tokens the model read (prompt_eval_count:
// the whole prompt, cached or not; after a cut, what was left). A count below
// this share of the estimate means the prompt was cut. Measured: an uncut
// request reads 0.74 of the estimate (gemma4:26b on the Demo, 8,437 of
// 11,386); one cut to a 4,096-token window read 0.18 (2,051).
export const READ_SHARE_LIMIT = 0.4;

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

/**
 * A conservative token count for the request: messages and tool schemas at
 * CHARS_PER_TOKEN, tool results and call arguments (JSON with numbers) at
 * JSON_CHARS_PER_TOKEN.
 */
export const estimatePromptTokens = ({ messages = [], tools = null } = {}) => {
    let chars = 0;
    let jsonChars = 0;
    for (const message of messages || []) {
        if (message.role === 'tool') jsonChars += textLength(message.content);
        else chars += textLength(message.content);
        for (const call of message.tool_calls || []) {
            chars += textLength(call.function?.name);
            jsonChars += textLength(call.function?.arguments);
        }
    }
    if (Array.isArray(tools) && tools.length) chars += JSON.stringify(tools).length;
    return Math.ceil(chars / CHARS_PER_TOKEN + jsonChars / JSON_CHARS_PER_TOKEN) + MESSAGE_OVERHEAD_TOKENS * (messages?.length || 0);
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

/** The key a conversation's `session` keeps a model's Ollama context under. */
export const sessionKey = (baseUrl, model) => `${baseUrl}|${model}`;

// LM Studio's refusal: the loaded window cannot hold the prompt and
// LMSTUDIO_REPLY_MIN_TOKENS of reply. `sent`: the request had already gone
// out (a model loaded by it, JIT) and its reply was stopped unseen.
const lmStudioTooSmall = ({ model, estimate, needed, loaded, sent }) => {
    const suggest = Math.max(MIN_CONTEXT, Math.ceil((estimate + REPLY_RESERVE_TOKENS) / CONTEXT_STEP) * CONTEXT_STEP);
    return new ContextLimitError(
        `LM Studio runs ${model} with a ${formatTokens(loaded)}-token context, but this request needs about `
        + `${formatTokens(needed)} (${WHAT_IT_HOLDS}). Reload the model in LM Studio with a context length of at least `
        + `${formatTokens(suggest)}. ${sent
            ? 'LM Studio loaded it for this request with that context, so the reply was stopped unseen'
            : 'Nothing was sent'}: the model would have answered without part of its instructions.`,
        { needed, available: loaded }
    );
};

/**
 * The guard's decision for one request, before anything is sent: the plan
 * the request goes out under, or a ContextLimitError when the model cannot
 * hold it. `server` is detectServer's answer (detected here when null).
 *   - Ollama: `{ server: 'ollama', version, estimate, needed, numCtx,
 *     numPredict, loadedContext, maxContext }` (planContext; the context
 *     this conversation asked before is read from `session`);
 *   - LM Studio: `{ server: 'lmstudio', estimate, needed, loadedContext }`,
 *     needed = estimate + LMSTUDIO_REPLY_MIN_TOKENS; loadedContext null when
 *     the model is not loaded (the request loads it, JIT, and its window is
 *     checked at the reply's first fragment);
 *   - anything else: `{ server: 'openai' }` (sent unchanged, never probed).
 * The plan is plain JSON: no key, no session.
 */
export const planCopilotChat = async ({ baseUrl, model, messages, tools }, { server = null, session = null, signal } = {}) => {
    const detected = server ?? await detectServer(baseUrl, { signal });
    if (detected.kind === 'ollama') {
        const estimate = estimatePromptTokens({ messages, tools });
        const [loaded, max] = await Promise.all([
            ollamaLoadedContext(baseUrl, model, { signal }),
            ollamaMaxContext(baseUrl, model, { signal })
        ]);
        const plan = planContext({ estimate, loaded, max, previous: session?.[sessionKey(baseUrl, model)] ?? null });
        if (plan.refuse) {
            throw new ContextLimitError(
                `This request needs about ${formatTokens(plan.needed)} tokens of context (${WHAT_IT_HOLDS}), but ${model} `
                + `holds at most ${formatTokens(plan.max)}. Start a new conversation, turn the analyses off, or pick a model `
                + 'with a longer context. Nothing was sent: the model would have answered without part of its instructions.',
                { needed: plan.needed, available: plan.max }
            );
        }
        return {
            server: 'ollama',
            version: detected.version ?? null,
            estimate,
            needed: plan.needed,
            numCtx: plan.numCtx,
            numPredict: plan.numPredict,
            loadedContext: loaded,
            maxContext: max
        };
    }
    if (detected.kind === 'lmstudio') {
        const estimate = estimatePromptTokens({ messages, tools });
        const needed = estimate + LMSTUDIO_REPLY_MIN_TOKENS;
        const loaded = await lmStudioLoadedContext(baseUrl, model, { signal });
        if (Number.isFinite(loaded) && loaded < needed) throw lmStudioTooSmall({ model, estimate, needed, loaded, sent: false });
        return { server: 'lmstudio', estimate, needed, loadedContext: Number.isFinite(loaded) ? loaded : null };
    }
    return { server: 'openai' };
};

/**
 * The Copilot's model call: one streamed chat that the model receives whole.
 * Yields client.js fragments. `plan` is planCopilotChat's decision for this
 * request (the model port makes it in `prepare`); without one it is made
 * here. `session` (an object the conversation keeps) remembers the context
 * asked for per model, written as the request goes out. Throws
 * ContextLimitError instead of sending a request the model would cut; after
 * sending, when Ollama read only part of the prompt, or a model LM Studio
 * loaded for the request (JIT) has a window too small (the reply stopped
 * unseen).
 */
export async function* streamCopilotChat({
    baseUrl, model, temperature, apiKey, messages, tools, toolChoice, signal, session = null
}, { plan = null } = {}) {
    const planned = plan ?? await planCopilotChat({ baseUrl, model, messages, tools }, { session, signal });
    if (planned.server === 'ollama') {
        const { estimate, numCtx, numPredict } = planned;
        if (session) session[sessionKey(baseUrl, model)] = numCtx;
        for await (const chunk of streamOllamaChat({
            baseUrl, model, messages, tools, temperature, apiKey, numCtx, numPredict, signal
        })) {
            // A backstop: the request was sized to fit, so a cut means the
            // server did not honour the context asked for.
            const read = chunk.usage?.promptTokens;
            if (Number.isFinite(read) && read < READ_SHARE_LIMIT * estimate) {
                throw new ContextLimitError(
                    `Ollama read only ${formatTokens(read)} of about ${formatTokens(estimate)} prompt tokens: it cut the prompt, `
                    + 'so this reply did not see all of its instructions and the run context. Do not rely on it; '
                    + `check the server's context settings (num_ctx ${formatTokens(numCtx)} was asked for).`,
                    { needed: estimate, available: read }
                );
            }
            yield chunk;
        }
        return;
    }
    if (planned.server === 'lmstudio' && !Number.isFinite(planned.loadedContext)) {
        // Not loaded: this request loads it (JIT) with LM Studio's default
        // context. Its window is read once it is loaded, at the first
        // fragment, and a window too small stops the reply unseen.
        const { estimate, needed } = planned;
        const request = new AbortController();
        const onAbort = () => request.abort();
        signal?.addEventListener('abort', onAbort, { once: true });
        try {
            let checked = false;
            for await (const chunk of streamChat({ baseUrl, model, temperature, apiKey, messages, tools, toolChoice, signal: request.signal })) {
                if (!checked) {
                    checked = true;
                    const now = await lmStudioLoadedContext(baseUrl, model, { signal });
                    if (Number.isFinite(now) && now < needed) {
                        request.abort();
                        throw lmStudioTooSmall({ model, estimate, needed, loaded: now, sent: true });
                    }
                }
                yield chunk;
            }
        } finally {
            signal?.removeEventListener('abort', onAbort);
        }
        return;
    }
    yield* streamChat({ baseUrl, model, temperature, apiKey, messages, tools, toolChoice, signal });
}
