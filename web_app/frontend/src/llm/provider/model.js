// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The model port: the one seam between a Copilot turn (copilot/turn.js and
// its loop, copilot/agentLoop.js) and whatever answers it. A port has two
// operations:
//
//   prepare(intent, { signal }) → sent       optional; the identity when absent
//   stream(sent, { signal })    → AsyncIterable<fragment>
//
// The intent is what the turn asks for: { messages, tools, toolChoice, retry }
// (`tools` and `toolChoice` undefined when the request carries no tools;
// `retry` null, or 'warmer' for a request asked again after a reply that kept
// repeating itself, since greedy decoding is what loops). `sent` is the
// request as it goes out: the turn shows it in its `request` event, and the
// validation keeps it in each record. `stream` sends it and yields client.js
// fragments ({ content }, { reasoning }, { toolCall }, { finishReason }). A
// bare function is a port whose prepare is the identity: fn(intent, { signal }).
// Either operation may return a promise.
//
// The port's typed failures, the only ones the loop reacts to by kind:
//   ContextLimitError      the request would not fit the model's window
//                          (contextGuard.js). Thrown by prepare, before
//                          anything is sent, so the turn emits no `request`
//                          for it; thrown by stream only after sending (Ollama
//                          read part of the prompt; a model LM Studio loaded
//                          for the request has too small a window)
//   ToolsUnsupportedError  the model or its server cannot take the request's
//                          tools. Thrown by stream before any fragment; the
//                          loop retries the turn's first request without tools
// Any other error is the provider's own (HTTP, network, a stopped reply).
//
// The adapters: providerModel below, the production one over the Copilot's
// provider (contextGuard.js); validation/spend.js meteredModel, a decorator
// that prices each request before it is sent; and the tests' scripted models.
// The provider quirks live in providerModel and nowhere above it: the server
// kind, the form each server takes tool_choice in (TOOL_CHOICE_FORMS), the
// context guard's plan, and what "tools unsupported" looks like per provider.

import { ContextLimitError, detectServer, planCopilotChat, streamCopilotChat } from './contextGuard';

export { ContextLimitError };

// The temperature a request asked again after a loop runs at, at least.
export const ANTI_REPEAT_TEMPERATURE = 0.6;

/** The model (or its server) cannot take the request's tools: retry without them. */
export class ToolsUnsupportedError extends Error {
    constructor(message, { status = null, server = null, cause = null } = {}) {
        super(message);
        this.name = 'ToolsUnsupportedError';
        this.status = status;
        this.server = server;
        if (cause) this.cause = cause;
    }
}

// How a request's tool choice goes out (`sent.plan.toolChoice`), by form:
export const TOOL_CHOICE_FORMS = Object.freeze({
    // OpenAI, Gemini, Anthropic (their OpenAI-compatible endpoints), LM Studio:
    // "tool_choice" in the /chat/completions body, beside the tools.
    body: 'tool_choice in the /chat/completions body, beside the tools',
    // Ollama, 'auto': its native /api/chat has no tool_choice; its default is auto.
    omitted: 'no tool_choice field (Ollama\'s native /api/chat has none; its default is auto), the tools sent',
    // Ollama, 'none': enforced by sending no tools, so no call can come back.
    tools_withheld: 'no tools sent (Ollama\'s native /api/chat has no tool_choice), so the model cannot call one'
});

/** The form a tool choice takes on a server kind ('ollama', 'lmstudio', 'openai'). */
export const toolChoiceForm = (kind, toolChoice) => {
    if (kind === 'ollama') return toolChoice === 'none' ? 'tools_withheld' : 'omitted';
    return toolChoice ? 'body' : 'omitted';
};

const hasTools = (tools) => Array.isArray(tools) && tools.length > 0;

// Ollama's capability check, on its native /api/chat (and /v1): HTTP 400,
// { "error": "<model> does not support tools" }.
const OLLAMA_NO_TOOLS = /\bdoes not support tools$/;
// The fallback, the one place a failure's text decides: a 4xx saying the
// model or server cannot take a `tools` field, for the providers whose shape
// is not known (LM Studio, Gemini, Anthropic). Deliberately narrow: a 400
// about a tool's schema is a real error and must be shown, not hidden by a
// retry without tools.
const TOOLS_UNSUPPORTED_TEXT = /(does not|doesn't|do not) support (tools|tool calling|tool use|function calling)|(tools|tool calling|tool use|function calling) (is |are )?not supported/i;

/**
 * Whether a provider failure says the model cannot take tools: by the
 * provider's own shape (an Error from client.js / ollama.js httpError, with
 * `status` and the parsed `body`), then by the fallback text. `server` is
 * the plan's server kind, or null when unknown (every shape is tried).
 *   - Ollama: 400, `{ error: "<model> does not support tools" }`;
 *   - OpenAI: 400, `{ error: { param: 'tools', code: 'unsupported_parameter' } }`;
 *   - any 4xx whose message matches TOOLS_UNSUPPORTED_TEXT (the fallback).
 */
export const isToolsUnsupported = (error, { server = null } = {}) => {
    if (error instanceof ToolsUnsupportedError) return true;
    const status = error?.status;
    if (!(Number.isInteger(status) && status >= 400 && status < 500)) return false;
    const detail = error.body?.error;
    if ((server === 'ollama' || server == null) && status === 400 && typeof detail === 'string' && OLLAMA_NO_TOOLS.test(detail)) {
        return true;
    }
    if (server !== 'ollama' && detail && typeof detail === 'object' && detail.param === 'tools' && detail.code === 'unsupported_parameter') {
        return true;
    }
    return TOOLS_UNSUPPORTED_TEXT.test(error.message || '');
};

// One request through the provider under its plan; a failure before any
// fragment that says the model cannot take the tools becomes a
// ToolsUnsupportedError (the provider's message kept).
async function* send(sent, { apiKey, session, signal }) {
    const { plan = null, ...request } = sent || {};
    let streamed = false;
    try {
        for await (const chunk of streamCopilotChat({ ...request, apiKey, session, signal }, { plan })) {
            streamed = true;
            yield chunk;
        }
    } catch (error) {
        const server = plan?.server ?? null;
        if (!streamed && hasTools(request.tools) && error?.name !== 'AbortError' && isToolsUnsupported(error, { server })) {
            throw new ToolsUnsupportedError(error.message, { status: error.status ?? null, server, cause: error });
        }
        throw error;
    }
}

const identity = (intent) => intent;

/** A model port from a port ({ prepare?, stream }) or a bare stream function; TypeError otherwise. */
export const asModelPort = (model) => {
    if (typeof model === 'function') return { prepare: identity, stream: model };
    if (model && typeof model.stream === 'function') {
        if (model.prepare != null && typeof model.prepare !== 'function') {
            throw new TypeError('the model port\'s prepare must be a function');
        }
        return {
            prepare: model.prepare ? (intent, options) => model.prepare(intent, options) : identity,
            stream: (sent, options) => model.stream(sent, options)
        };
    }
    throw new TypeError('model must be a model port ({ prepare?, stream }) or a stream function');
};

/**
 * The production port over the Copilot's provider (contextGuard.js), with
 * the user's settings ({ baseUrl, model, temperature, apiKey }), read when the
 * port is made. `session` is the object a conversation keeps for the context
 * asked per model.
 *
 * prepare (async) detects the server (contextGuard.js detectServer: only
 * localhost or Ollama's port is probed), puts the tool choice in that
 * server's form, and asks the context guard for its plan, so `sent` is
 *   { baseUrl, model, temperature, messages, tools, toolChoice, plan }
 * exactly as it goes out: a 'warmer' retry at ANTI_REPEAT_TEMPERATURE or
 * above; on Ollama, tool_choice 'none' as no tools (the native API has no
 * tool_choice); `plan` the guard's decision (planCopilotChat: the server kind,
 * Ollama's numCtx / numPredict, LM Studio's loaded context, the estimate) and
 * the tool choice's form (TOOL_CHOICE_FORMS). A request the guard refuses
 * throws ContextLimitError here, before anything is sent. `sent` never holds
 * the API key, which joins the request only as it is sent, so no request
 * event or record carries it.
 *
 * stream sends `sent` under its plan; a model that cannot take tools fails
 * with a ToolsUnsupportedError.
 */
export const providerModel = (settings, { session = null } = {}) => {
    const { baseUrl, model, temperature, apiKey } = settings || {};
    return {
        prepare: async (intent, { signal } = {}) => {
            const server = await detectServer(baseUrl, { signal });
            const asked = {
                baseUrl,
                model,
                temperature: intent?.retry === 'warmer'
                    ? Math.max(Number(temperature) || 0, ANTI_REPEAT_TEMPERATURE)
                    : temperature,
                messages: intent?.messages,
                tools: intent?.tools,
                toolChoice: intent?.toolChoice
            };
            const form = hasTools(asked.tools) ? toolChoiceForm(server.kind, asked.toolChoice) : null;
            const request = form === 'tools_withheld' ? { ...asked, tools: undefined } : asked;
            const plan = await planCopilotChat(request, { server, session, signal });
            return { ...request, plan: form ? { ...plan, toolChoice: form } : plan };
        },
        stream: (sent, { signal } = {}) => send(sent, { apiKey, session, signal })
    };
};
