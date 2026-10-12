// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Ollama's native chat API (`POST /api/chat`), used for the Copilot's own
// requests to an Ollama server because it is the only Ollama endpoint that
// takes a context length. The OpenAI-compatible /v1 endpoint runs a model at
// the server's default context (OLLAMA_CONTEXT_LENGTH, or one chosen from the
// GPU memory: 4096 on a machine without a usable GPU), and a longer prompt is
// cut to fit without an error — the model then answers without the rules or
// the run context it was given. contextGuard.js sizes `num_ctx` from the
// request; this file speaks the dialect.
//
// The stream is newline-delimited JSON: `{ message: { content, thinking,
// tool_calls }, done: false }` chunks, then `{ done: true, done_reason,
// prompt_eval_count, eval_count, … }`, or `{ error }`. It yields the same
// fragments as client.js streamChat (`{ content }`, `{ reasoning }`,
// `{ toolCall }`, `{ finishReason }`), plus `{ usage }` from the last chunk.

import { describeHttpError } from './client';

/** The server root of an OpenAI-style base URL: `http://host:11434/v1` → `http://host:11434`. */
export const serverOrigin = (baseUrl) => (baseUrl || '').replace(/\/+$/, '').replace(/\/v1$/, '');

const parseArguments = (text) => {
    if (text && typeof text === 'object') return text;
    try {
        const value = JSON.parse(text || '{}');
        return value && typeof value === 'object' && !Array.isArray(value) ? value : {};
    } catch {
        return {};
    }
};

/**
 * OpenAI-style messages → Ollama's: content is always a string, a call's
 * arguments are an object, and a tool result names its tool (tool_name) as
 * well as the call it answers.
 */
export const toOllamaMessages = (messages) => {
    const namesById = new Map();
    return (messages || []).map((message) => {
        if (message.role === 'assistant' && Array.isArray(message.tool_calls)) {
            message.tool_calls.forEach((call) => namesById.set(call.id, call.function?.name));
            return {
                role: 'assistant',
                content: message.content || '',
                tool_calls: message.tool_calls.map((call) => ({
                    ...(call.id ? { id: call.id } : {}),
                    function: { name: call.function?.name, arguments: parseArguments(call.function?.arguments) }
                }))
            };
        }
        if (message.role === 'tool') {
            const name = namesById.get(message.tool_call_id);
            return {
                role: 'tool',
                content: message.content || '',
                ...(message.tool_call_id ? { tool_call_id: message.tool_call_id } : {}),
                ...(name ? { tool_name: name } : {})
            };
        }
        return { role: message.role, content: message.content || '' };
    });
};

/**
 * Stream one chat from Ollama's native API. `numCtx` is the context the model
 * must run with; `numPredict` caps the reply so it cannot outgrow it.
 * `tool_choice` has no native equivalent and is not sent (Ollama ignores it
 * on /v1 too).
 */
export async function* streamOllamaChat({
    baseUrl, model, messages, tools, temperature = 0.2, apiKey = '', numCtx, numPredict, signal
}) {
    const body = {
        model,
        messages: toOllamaMessages(messages),
        stream: true,
        options: {
            temperature,
            ...(Number.isFinite(numCtx) ? { num_ctx: numCtx } : {}),
            ...(Number.isFinite(numPredict) ? { num_predict: numPredict } : {})
        }
    };
    if (Array.isArray(tools) && tools.length) body.tools = tools;
    // An Ollama behind a proxy that wants a key gets it, as on /v1.
    const headers = { 'Content-Type': 'application/json' };
    if (apiKey) headers.Authorization = `Bearer ${apiKey}`;
    const response = await fetch(`${serverOrigin(baseUrl)}/api/chat`, {
        method: 'POST',
        headers,
        body: JSON.stringify(body),
        signal
    });
    if (!response.ok) {
        const error = new Error(await describeHttpError(response));
        error.status = response.status;
        throw error;
    }
    if (!response.body) throw new Error('The server returned no response body to stream');
    const reader = response.body.getReader();
    const decoder = new TextDecoder();
    let buffer = '';
    let calls = 0;
    const handle = function* (line) {
        const text = line.trim();
        if (!text) return;
        let chunk;
        try {
            chunk = JSON.parse(text);
        } catch {
            return;
        }
        if (chunk?.error) {
            const detail = typeof chunk.error === 'string' ? chunk.error : (chunk.error.message || JSON.stringify(chunk.error));
            throw new Error(`The server stopped the reply: ${detail}`);
        }
        const message = chunk?.message;
        if (message?.thinking) yield { reasoning: message.thinking };
        if (message?.content) yield { content: message.content };
        if (Array.isArray(message?.tool_calls)) {
            for (const call of message.tool_calls) {
                if (!call?.function?.name) continue;
                const args = call.function.arguments;
                yield {
                    toolCall: {
                        // Ollama sends each call whole: one index per call.
                        index: calls,
                        id: typeof call.id === 'string' && call.id ? call.id : null,
                        name: call.function.name,
                        arguments: typeof args === 'string' ? args : JSON.stringify(args ?? {}),
                        extraContent: null
                    }
                };
                calls += 1;
            }
        }
        if (chunk?.done) {
            yield {
                usage: {
                    promptTokens: Number.isFinite(chunk.prompt_eval_count) ? chunk.prompt_eval_count : null,
                    outputTokens: Number.isFinite(chunk.eval_count) ? chunk.eval_count : null,
                    doneReason: chunk.done_reason || null
                }
            };
            yield { finishReason: chunk.done_reason === 'length' ? 'length' : (calls ? 'tool_calls' : 'stop') };
        }
    };
    try {
        for (;;) {
            const { done, value } = await reader.read();
            if (done) break;
            buffer += decoder.decode(value, { stream: true });
            const lines = buffer.split('\n');
            buffer = lines.pop() ?? '';
            for (const line of lines) yield* handle(line);
        }
        buffer += decoder.decode();
        yield* handle(buffer);
    } finally {
        reader.cancel().catch(() => {});
    }
}
