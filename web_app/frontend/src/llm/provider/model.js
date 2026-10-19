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
//
// The adapters: providerModel below, the production one over the Copilot's
// provider (contextGuard.js streamCopilotChat); validation/spend.js
// meteredModel, a decorator that prices each request before it is sent; and
// the tests' scripted models. Provider quirks (tool_choice forms, typed
// "tools unsupported" failures, the context plan in `sent`) belong behind
// this interface, in providerModel.

import { streamCopilotChat } from './contextGuard';

// The temperature a request asked again after a loop runs at, at least.
export const ANTI_REPEAT_TEMPERATURE = 0.6;

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
 * The production port: each request through streamCopilotChat (the context
 * guard, Ollama's native API, provider/contextGuard.js) with the user's
 * settings ({ baseUrl, model, temperature, apiKey }), read when the port is
 * made. `session` is the object a conversation keeps for the context asked
 * per model. A 'warmer' retry runs at ANTI_REPEAT_TEMPERATURE or above. `sent`
 * is { baseUrl, model, temperature, messages, tools, toolChoice }: it never
 * holds the API key, which joins the request only as it is sent, so no
 * request event or record carries it.
 */
export const providerModel = (settings, { session = null } = {}) => {
    const { baseUrl, model, temperature, apiKey } = settings || {};
    return {
        prepare: (intent) => ({
            baseUrl,
            model,
            temperature: intent?.retry === 'warmer'
                ? Math.max(Number(temperature) || 0, ANTI_REPEAT_TEMPERATURE)
                : temperature,
            messages: intent?.messages,
            tools: intent?.tools,
            toolChoice: intent?.toolChoice
        }),
        stream: (sent, { signal } = {}) => streamCopilotChat({ ...sent, apiKey, session, signal })
    };
};
