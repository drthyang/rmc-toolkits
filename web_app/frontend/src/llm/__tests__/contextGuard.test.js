// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { isToolsUnsupportedError } from '../copilot/agentLoop';
import {
    CHARS_PER_TOKEN, CONTEXT_STEP, ContextLimitError, MIN_CONTEXT, REPLY_RESERVE_TOKENS,
    detectServer, estimatePromptTokens, forgetServers, planContext, streamCopilotChat
} from '../provider/contextGuard';
import { serverOrigin, streamOllamaChat, toOllamaMessages } from '../provider/ollama';

const encoder = new TextEncoder();
const streamOf = (chunks) => new ReadableStream({
    start(controller) {
        chunks.forEach((chunk) => controller.enqueue(encoder.encode(chunk)));
        controller.close();
    }
});
const ndjson = (objects) => objects.map((object) => `${JSON.stringify(object)}\n`).join('');
const json = (body, status = 200) => ({ ok: status < 400, status, json: async () => body });

const collect = async (iterator) => {
    const parts = [];
    for await (const part of iterator) parts.push(part);
    return parts;
};

// A request the size of the Copilot's first one on a run: rules, ~6,000
// characters of run context, ten tool schemas, a question.
const bigRequest = () => ({
    messages: [
        { role: 'system', content: 'r'.repeat(5200) },
        { role: 'user', content: `Context:\n${'c'.repeat(6000)}` },
        { role: 'assistant', content: 'Understood.' },
        { role: 'user', content: 'What are the Se–Ta–Se angles?' }
    ],
    tools: Array.from({ length: 10 }, (_, index) => ({
        type: 'function',
        function: { name: `tool_${index}`, description: 'd'.repeat(600), parameters: { type: 'object', properties: {} } }
    }))
});

// An Ollama server: /api/version, /api/ps, /api/show, /api/chat (NDJSON).
const ollamaServer = ({ loaded = null, max = 262144, reply = null } = {}) => {
    const calls = [];
    const fetchMock = vi.fn(async (url, options = {}) => {
        const path = new URL(url).pathname;
        calls.push({ path, body: options.body ? JSON.parse(options.body) : null });
        if (path === '/api/version') return json({ version: '0.40.2' });
        if (path === '/api/ps') {
            return json({ models: loaded ? [{ name: 'gemma4:26b', model: 'gemma4:26b', context_length: loaded }] : [] });
        }
        if (path === '/api/show') return json({ model_info: { 'general.architecture': 'gemma4', 'gemma4.context_length': max } });
        if (path === '/api/chat') {
            return { ok: true, status: 200, body: streamOf([reply ?? ndjson([
                { message: { role: 'assistant', content: 'Hello' }, done: false },
                { message: { role: 'assistant', content: '' }, done: true, done_reason: 'stop', prompt_eval_count: 5000, eval_count: 3 }
            ])]) };
        }
        return json({ error: 'not found' }, 404);
    });
    vi.stubGlobal('fetch', fetchMock);
    return { calls, fetchMock };
};

beforeEach(() => forgetServers());
afterEach(() => vi.unstubAllGlobals());

describe('the prompt estimate', () => {
    it('counts every message, tool call and schema, and errs high', () => {
        const request = bigRequest();
        const chars = request.messages.reduce((sum, message) => sum + message.content.length, 0) + JSON.stringify(request.tools).length;
        const estimate = estimatePromptTokens(request);
        expect(estimate).toBe(Math.ceil(chars / CHARS_PER_TOKEN) + 8 * 4);
        // At least as many tokens as a 3.5-characters-per-token tokenizer makes.
        expect(estimate).toBeGreaterThan(chars / 3.5);
        const withCall = estimatePromptTokens({
            messages: [{ role: 'assistant', content: null, tool_calls: [{ function: { name: 'echo', arguments: '{"value":1}' } }] }]
        });
        expect(withCall).toBe(Math.ceil(15 / CHARS_PER_TOKEN) + 8);
    });
});

describe('planContext', () => {
    it('keeps the context the model is loaded with when the request fits it', () => {
        expect(planContext({ estimate: 10000, loaded: 262144, max: 262144 })).toEqual({
            numCtx: 262144, numPredict: 252144, needed: 10000 + REPLY_RESERVE_TOKENS
        });
    });

    it('asks for a fitting context when the loaded one would cut the prompt (a 4,096 default)', () => {
        const plan = planContext({ estimate: 10000, loaded: 4096, max: 262144 });
        expect(plan.numCtx).toBe(Math.ceil((10000 + REPLY_RESERVE_TOKENS) / CONTEXT_STEP) * CONTEXT_STEP);
        expect(plan.numCtx).toBeGreaterThanOrEqual(10000 + REPLY_RESERVE_TOKENS);
        expect(plan.numPredict).toBe(plan.numCtx - 10000);
        expect(planContext({ estimate: 100, loaded: null, max: null }).numCtx).toBe(MIN_CONTEXT);
    });

    it('keeps the context this conversation asked for while it still fits', () => {
        expect(planContext({ estimate: 11000, loaded: null, max: 262144, previous: 24576 }).numCtx).toBe(24576);
        expect(planContext({ estimate: 20000, loaded: null, max: 262144, previous: 24576 }).numCtx).toBe(32768);
    });

    it('never asks for more than the model holds, and refuses what it cannot hold', () => {
        expect(planContext({ estimate: 25000, loaded: null, max: 40960 }).numCtx).toBe(40960);
        expect(planContext({ estimate: 7000, loaded: null, max: 8192 })).toEqual({ refuse: true, needed: 7000 + REPLY_RESERVE_TOKENS, max: 8192 });
    });
});

describe('detectServer', () => {
    it('knows Ollama by its version endpoint', async () => {
        ollamaServer();
        await expect(detectServer('http://localhost:11434/v1')).resolves.toEqual({ kind: 'ollama', version: '0.40.2' });
    });

    it('knows LM Studio by its model list, which says what is loaded', async () => {
        vi.stubGlobal('fetch', vi.fn(async (url) => (String(url).endsWith('/api/version')
            ? json({ error: 'Unexpected endpoint or method. (GET /api/version)' })
            : json({ data: [{ id: 'qwen/qwen3.8-27b', state: 'loaded', loaded_context_length: 80000 }] }))));
        await expect(detectServer('http://localhost:1234/v1')).resolves.toEqual({ kind: 'lmstudio' });
    });

    it('never probes a cloud provider, and asks an unreachable server again', async () => {
        const fetchMock = vi.fn(async () => { throw new TypeError('Failed to fetch'); });
        vi.stubGlobal('fetch', fetchMock);
        await expect(detectServer('https://api.openai.com/v1')).resolves.toEqual({ kind: 'openai' });
        expect(fetchMock).not.toHaveBeenCalled();
        await detectServer('http://localhost:11434/v1');
        await detectServer('http://localhost:11434/v1');
        expect(fetchMock).toHaveBeenCalledTimes(2);
    });
});

describe('the Copilot on Ollama', () => {
    it('sends the native request with a context that holds the prompt, when the loaded one would cut it', async () => {
        const { calls } = ollamaServer({ loaded: 4096 });
        const request = bigRequest();
        const estimate = estimatePromptTokens(request);
        const parts = await collect(streamCopilotChat({
            baseUrl: 'http://localhost:11434/v1', model: 'gemma4:26b', temperature: 0.2, ...request, toolChoice: 'auto'
        }));
        expect(parts).toEqual([
            { content: 'Hello' },
            { usage: { promptTokens: 5000, outputTokens: 3, doneReason: 'stop' } },
            { finishReason: 'stop' }
        ]);
        const chat = calls.find((call) => call.path === '/api/chat').body;
        expect(chat.options.num_ctx).toBeGreaterThanOrEqual(estimate + REPLY_RESERVE_TOKENS);
        expect(chat.options.num_ctx % CONTEXT_STEP).toBe(0);
        expect(chat.options.num_predict).toBe(chat.options.num_ctx - estimate);
        expect(chat.options.temperature).toBe(0.2);
        expect(chat.tools).toHaveLength(10);
        expect('tool_choice' in chat).toBe(false);
        expect(calls.some((call) => call.path === '/v1/chat/completions')).toBe(false);
    });

    it('keeps a loaded context that fits, so the model is not reloaded', async () => {
        const { calls } = ollamaServer({ loaded: 262144 });
        await collect(streamCopilotChat({ baseUrl: 'http://localhost:11434/v1', model: 'gemma4:26b', ...bigRequest() }));
        expect(calls.find((call) => call.path === '/api/chat').body.options.num_ctx).toBe(262144);
    });

    it('keeps the context it asked for across the conversation', async () => {
        const { calls } = ollamaServer({ loaded: null });
        const session = {};
        await collect(streamCopilotChat({ baseUrl: 'http://localhost:11434/v1', model: 'gemma4:26b', session, ...bigRequest() }));
        const first = calls.filter((call) => call.path === '/api/chat').at(-1).body.options.num_ctx;
        const longer = bigRequest();
        longer.messages.push({ role: 'tool', tool_call_id: 'x', content: 't'.repeat(500) });
        await collect(streamCopilotChat({ baseUrl: 'http://localhost:11434/v1', model: 'gemma4:26b', session, ...longer }));
        expect(calls.filter((call) => call.path === '/api/chat').at(-1).body.options.num_ctx).toBe(first);
        expect(session['http://localhost:11434/v1|gemma4:26b']).toBe(first);
        // It moves up a step only when the conversation outgrows it.
        longer.messages.push({ role: 'tool', tool_call_id: 'y', content: 't'.repeat(6000) });
        await collect(streamCopilotChat({ baseUrl: 'http://localhost:11434/v1', model: 'gemma4:26b', session, ...longer }));
        expect(calls.filter((call) => call.path === '/api/chat').at(-1).body.options.num_ctx).toBe(first + CONTEXT_STEP);
    });

    it('refuses, and sends nothing, when the model cannot hold the request', async () => {
        const { calls } = ollamaServer({ max: 8192 });
        const run = collect(streamCopilotChat({ baseUrl: 'http://localhost:11434/v1', model: 'gemma4:26b', ...bigRequest() }));
        await expect(run).rejects.toBeInstanceOf(ContextLimitError);
        await expect(collect(streamCopilotChat({ baseUrl: 'http://localhost:11434/v1', model: 'gemma4:26b', ...bigRequest() })))
            .rejects.toThrow(/needs about [\d,]+ tokens of context .* but gemma4:26b holds at most 8,192\. .* Nothing was sent/);
        expect(calls.some((call) => call.path === '/api/chat')).toBe(false);
    });
});

describe('the Copilot on LM Studio', () => {
    const lmStudio = (loadedContext) => {
        const fetchMock = vi.fn(async (url) => {
            const path = new URL(url).pathname;
            if (path === '/api/version') return json({ error: 'Unexpected endpoint or method. (GET /api/version)' });
            if (path === '/api/v0/models') return json({ data: [{ id: 'qwen/qwen3.8-27b', state: 'loaded', loaded_context_length: loadedContext }] });
            if (path === '/v1/chat/completions') {
                return { ok: true, status: 200, body: streamOf([`data: ${JSON.stringify({ choices: [{ delta: { content: 'ok' } }] })}\n\ndata: [DONE]\n\n`]) };
            }
            return json({ error: 'not found' }, 404);
        });
        vi.stubGlobal('fetch', fetchMock);
        return fetchMock;
    };

    it('refuses when the loaded model\'s window is too small, and says how to fix it', async () => {
        const fetchMock = lmStudio(4096);
        await expect(collect(streamCopilotChat({ baseUrl: 'http://localhost:1234/v1', model: 'qwen/qwen3.8-27b', ...bigRequest() })))
            .rejects.toThrow(/LM Studio runs qwen\/qwen3\.8-27b with a 4,096-token context, but .* Reload the model in LM Studio with a context length of at least/);
        expect(fetchMock.mock.calls.some(([url]) => String(url).endsWith('/chat/completions'))).toBe(false);
    });

    it('sends the request unchanged when it fits', async () => {
        lmStudio(80000);
        const parts = await collect(streamCopilotChat({ baseUrl: 'http://localhost:1234/v1', model: 'qwen/qwen3.8-27b', ...bigRequest() }));
        expect(parts).toEqual([{ content: 'ok' }]);
    });
});

describe('Ollama native chat', () => {
    it('converts OpenAI-style messages: string content, object arguments, named tool results', () => {
        expect(serverOrigin('http://localhost:11434/v1/')).toBe('http://localhost:11434');
        expect(toOllamaMessages([
            { role: 'system', content: 'S' },
            { role: 'assistant', content: null, tool_calls: [{ id: 'c1', type: 'function', function: { name: 'find_gr_peaks', arguments: '{"pair":"Ta-Se"}' } }] },
            { role: 'tool', tool_call_id: 'c1', content: '{"data":1}' }
        ])).toEqual([
            { role: 'system', content: 'S' },
            { role: 'assistant', content: '', tool_calls: [{ id: 'c1', function: { name: 'find_gr_peaks', arguments: { pair: 'Ta-Se' } } }] },
            { role: 'tool', content: '{"data":1}', tool_call_id: 'c1', tool_name: 'find_gr_peaks' }
        ]);
    });

    it('streams thinking, content and whole tool calls from NDJSON split anywhere', async () => {
        const text = ndjson([
            { message: { role: 'assistant', content: '', thinking: 'Which shell…' }, done: false },
            { message: { role: 'assistant', content: 'Checking.' }, done: false },
            { message: { role: 'assistant', content: '', tool_calls: [
                { function: { name: 'find_gr_peaks', arguments: { pair: 'Ta-Se' } } },
                { id: 'call_7', function: { name: 'list_sites', arguments: {} } }
            ] }, done: false },
            { message: { role: 'assistant', content: '' }, done: true, done_reason: 'stop', prompt_eval_count: 812, eval_count: 40 }
        ]);
        vi.stubGlobal('fetch', vi.fn(async () => ({ ok: true, status: 200, body: streamOf([text.slice(0, 17), text.slice(17, 90), text.slice(90)]) })));
        const parts = await collect(streamOllamaChat({ baseUrl: 'http://localhost:11434/v1', model: 'm', messages: [] }));
        expect(parts).toEqual([
            { reasoning: 'Which shell…' },
            { content: 'Checking.' },
            { toolCall: { index: 0, id: null, name: 'find_gr_peaks', arguments: '{"pair":"Ta-Se"}', extraContent: null } },
            { toolCall: { index: 1, id: 'call_7', name: 'list_sites', arguments: '{}', extraContent: null } },
            { usage: { promptTokens: 812, outputTokens: 40, doneReason: 'stop' } },
            { finishReason: 'tool_calls' }
        ]);
    });

    it('says when the reply stopped at the length limit, and surfaces errors', async () => {
        vi.stubGlobal('fetch', vi.fn(async () => ({ ok: true, status: 200, body: streamOf([ndjson([
            { message: { content: 'Partial' }, done: false },
            { message: { content: '' }, done: true, done_reason: 'length', prompt_eval_count: 9, eval_count: 99 }
        ])]) })));
        const parts = await collect(streamOllamaChat({ baseUrl: 'http://x:11434', model: 'm', messages: [] }));
        expect(parts.at(-1)).toEqual({ finishReason: 'length' });

        vi.stubGlobal('fetch', vi.fn(async () => ({ ok: true, status: 200, body: streamOf([ndjson([{ error: 'model runner has unexpectedly stopped' }])]) })));
        await expect(collect(streamOllamaChat({ baseUrl: 'http://x:11434', model: 'm', messages: [] })))
            .rejects.toThrow('The server stopped the reply: model runner has unexpectedly stopped');

        // A model without tool calling: the agent loop's fallback must recognise it.
        vi.stubGlobal('fetch', vi.fn(async () => json({ error: 'registry.ollama.ai/library/gemma3:latest does not support tools' }, 400)));
        const failure = await collect(streamOllamaChat({ baseUrl: 'http://x:11434', model: 'gemma3', messages: [], tools: [{}] })).catch((error) => error);
        expect(failure.status).toBe(400);
        expect(isToolsUnsupportedError(failure)).toBe(true);
    });
});
