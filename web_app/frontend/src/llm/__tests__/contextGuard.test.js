// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import {
    CHARS_PER_TOKEN, CONTEXT_STEP, ContextLimitError, JSON_CHARS_PER_TOKEN, LMSTUDIO_REPLY_MIN_TOKENS, MIN_CONTEXT, REPLY_RESERVE_TOKENS,
    detectServer, estimatePromptTokens, forgetServers, planContext, planCopilotChat, streamCopilotChat
} from '../provider/contextGuard';
import { isToolsUnsupported } from '../provider/model';
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
        expect(withCall).toBe(Math.ceil(4 / CHARS_PER_TOKEN + 11 / JSON_CHARS_PER_TOKEN) + 8);
    });

    it('errs high on tool results too, so the reply cap never exceeds the room left (review 2026-10-12 #8)', () => {
        // Tool JSON dense with numbers: 1.64 characters per token on gemma4:26b.
        const result = JSON.stringify({ data: { rows: Array.from({ length: 800 }, (_, i) => [i / 7, (i * 3) / 11]) } });
        const real = result.length / 1.64;
        const estimate = estimatePromptTokens({ messages: [{ role: 'tool', tool_call_id: 'c1', content: result }] });
        expect(estimate).toBeGreaterThan(real);
        // So numPredict = numCtx − estimate stays inside what the prompt leaves.
        const plan = planContext({ estimate, loaded: 32768, max: 262144 });
        expect(plan.numPredict).toBeLessThanOrEqual(plan.numCtx - real);
        // And an uncut request still reads far above the backstop's 0.4 share.
        expect(real / estimate).toBeGreaterThan(0.75);
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

describe('planCopilotChat: the guard decides before anything is sent (the model port\'s prepare)', () => {
    it('Ollama: the plan the request goes out under, read from the probes only, the session untouched', async () => {
        const { calls } = ollamaServer({ loaded: 4096 });
        const request = bigRequest();
        const estimate = estimatePromptTokens(request);
        const session = {};
        const plan = await planCopilotChat({ baseUrl: 'http://localhost:11434/v1', model: 'gemma4:26b', ...request }, { session });
        const numCtx = Math.max(MIN_CONTEXT, Math.ceil((estimate + REPLY_RESERVE_TOKENS) / CONTEXT_STEP) * CONTEXT_STEP);
        expect(plan).toEqual({
            server: 'ollama', version: '0.40.2', estimate, needed: estimate + REPLY_RESERVE_TOKENS,
            numCtx, numPredict: numCtx - estimate, loadedContext: 4096, maxContext: 262144
        });
        expect(calls.map((call) => call.path)).toEqual(['/api/version', '/api/ps', '/api/show']);
        // The session keeps what was sent: a plan alone is not a request.
        expect(session).toEqual({});
        // Sent under that plan: no probe again, the plan's num_ctx and num_predict.
        calls.length = 0;
        await collect(streamCopilotChat({ baseUrl: 'http://localhost:11434/v1', model: 'gemma4:26b', session, ...request }, { plan }));
        expect(calls.map((call) => call.path)).toEqual(['/api/chat']);
        expect(calls[0].body.options).toMatchObject({ num_ctx: plan.numCtx, num_predict: plan.numPredict });
        expect(session['http://localhost:11434/v1|gemma4:26b']).toBe(plan.numCtx);
    });

    it('Ollama: refuses a request the model cannot hold, with the same message, and sends nothing', async () => {
        const { calls } = ollamaServer({ max: 8192 });
        const failure = await planCopilotChat({ baseUrl: 'http://localhost:11434/v1', model: 'gemma4:26b', ...bigRequest() }).catch((error) => error);
        expect(failure).toBeInstanceOf(ContextLimitError);
        expect(failure.message).toMatch(/needs about [\d,]+ tokens of context .* but gemma4:26b holds at most 8,192\. .* Nothing was sent/);
        expect(failure).toMatchObject({ available: 8192 });
        expect(calls.some((call) => call.path === '/api/chat')).toBe(false);
    });

    it('LM Studio: the loaded window and the 4,096-token room; too small refused; not loaded left to the JIT check', async () => {
        const request = bigRequest();
        const estimate = estimatePromptTokens(request);
        const studio = (entry) => vi.stubGlobal('fetch', vi.fn(async (url) => (new URL(url).pathname === '/api/v0/models'
            ? json({ data: [{ id: 'qwen/qwen3.8-27b', ...entry }] })
            : json({ error: 'Unexpected endpoint or method. (GET /api/version)' }))));
        studio({ state: 'loaded', loaded_context_length: 80000 });
        await expect(planCopilotChat({ baseUrl: 'http://localhost:1234/v1', model: 'qwen/qwen3.8-27b', ...request }))
            .resolves.toEqual({ server: 'lmstudio', estimate, needed: estimate + LMSTUDIO_REPLY_MIN_TOKENS, loadedContext: 80000 });
        forgetServers();
        studio({ state: 'loaded', loaded_context_length: 4096 });
        await expect(planCopilotChat({ baseUrl: 'http://localhost:1234/v1', model: 'qwen/qwen3.8-27b', ...request }))
            .rejects.toThrow(/LM Studio runs qwen\/qwen3\.8-27b with a 4,096-token context, .* Nothing was sent/);
        forgetServers();
        studio({ state: 'not-loaded' });
        await expect(planCopilotChat({ baseUrl: 'http://localhost:1234/v1', model: 'qwen/qwen3.8-27b', ...request }))
            .resolves.toMatchObject({ server: 'lmstudio', loadedContext: null });
    });

    it('any other server: sent unchanged, never probed', async () => {
        const fetchMock = vi.fn(async () => { throw new Error('probed'); });
        vi.stubGlobal('fetch', fetchMock);
        await expect(planCopilotChat({ baseUrl: 'https://api.openai.com/v1', model: 'gpt', ...bigRequest() })).resolves.toEqual({ server: 'openai' });
        expect(fetchMock).not.toHaveBeenCalled();
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

    it('stops the turn when the server reports it read only part of the prompt', async () => {
        // Seen with a 4,096-token window: Ollama cut the 8,437-token first request to 2,051.
        ollamaServer({
            loaded: 262144,
            reply: ndjson([
                { message: { role: 'assistant', content: 'Looking at the view.' }, done: false },
                { message: { role: 'assistant', content: '' }, done: true, done_reason: 'stop', prompt_eval_count: 2051, eval_count: 30 }
            ])
        });
        const parts = [];
        const run = (async () => {
            for await (const part of streamCopilotChat({ baseUrl: 'http://localhost:11434/v1', model: 'gemma4:26b', ...bigRequest() })) parts.push(part);
        })();
        await expect(run).rejects.toThrow(/Ollama read only 2,051 of about [\d,]+ prompt tokens: it cut the prompt/);
        expect(parts).toEqual([{ content: 'Looking at the view.' }]);
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

    // LM Studio loads a model on its first request (JIT), with its own
    // default context: not loaded when asked, loaded by the time it answers.
    const jitStudio = (loadedContext) => {
        let started = false;
        let chatSignal = null;
        vi.stubGlobal('fetch', vi.fn(async (url, options = {}) => {
            const path = new URL(url).pathname;
            if (path === '/api/version') return json({ error: 'Unexpected endpoint or method. (GET /api/version)' });
            if (path === '/api/v0/models') {
                return json({ data: [started
                    ? { id: 'qwen/qwen3.8-27b', state: 'loaded', loaded_context_length: loadedContext }
                    : { id: 'qwen/qwen3.8-27b', state: 'not-loaded', max_context_length: 262144 }] });
            }
            if (path === '/v1/chat/completions') {
                started = true;
                chatSignal = options.signal;
                return { ok: true, status: 200, body: streamOf([`data: ${JSON.stringify({ choices: [{ delta: { content: 'ok' } }] })}\n\ndata: [DONE]\n\n`]) };
            }
            return json({ error: 'not found' }, 404);
        }));
        return { signal: () => chatSignal };
    };

    it('checks a model it loads for the request, and stops a reply its window cut, unseen (review 2026-10-12 #8)', async () => {
        const server = jitStudio(4096);
        const parts = [];
        await expect((async () => {
            for await (const part of streamCopilotChat({ baseUrl: 'http://localhost:1234/v1', model: 'qwen/qwen3.8-27b', ...bigRequest() })) parts.push(part);
        })()).rejects.toThrow(/LM Studio runs qwen\/qwen3\.8-27b with a 4,096-token context, .* LM Studio loaded it for this request with that context, so the reply was stopped unseen/);
        expect(parts).toEqual([]);
        expect(server.signal().aborted).toBe(true);
        forgetServers();
        jitStudio(80000);
        await expect(collect(streamCopilotChat({ baseUrl: 'http://localhost:1234/v1', model: 'qwen/qwen3.8-27b', ...bigRequest() })))
            .resolves.toEqual([{ content: 'ok' }]);
    });
});

describe('Ollama native chat', () => {
    it('sends the API key, as the /v1 path does (review 2026-10-12 #8)', async () => {
        const fetchMock = vi.fn(async () => ({ ok: true, status: 200, body: streamOf([ndjson([{ message: { content: '' }, done: true, done_reason: 'stop' }])]) }));
        vi.stubGlobal('fetch', fetchMock);
        await collect(streamOllamaChat({ baseUrl: 'http://gpu-box:11434/v1', model: 'gemma4:26b', messages: [{ role: 'user', content: 'q' }], apiKey: 'k-1' }));
        expect(fetchMock.mock.calls[0][1].headers).toEqual({ 'Content-Type': 'application/json', Authorization: 'Bearer k-1' });
        await collect(streamOllamaChat({ baseUrl: 'http://localhost:11434/v1', model: 'gemma4:26b', messages: [{ role: 'user', content: 'q' }] }));
        expect(fetchMock.mock.calls[1][1].headers).toEqual({ 'Content-Type': 'application/json' });
    });

    it('leaves a null temperature out of the options; a number goes in as before', async () => {
        const fetchMock = vi.fn(async () => ({ ok: true, status: 200, body: streamOf([ndjson([{ message: { content: '' }, done: true, done_reason: 'stop' }])]) }));
        vi.stubGlobal('fetch', fetchMock);
        const base = { baseUrl: 'http://localhost:11434/v1', model: 'm', messages: [{ role: 'user', content: 'q' }], numCtx: 16384, numPredict: 8000 };
        await collect(streamOllamaChat({ ...base, temperature: null }));
        await collect(streamOllamaChat({ ...base, temperature: 0.2 }));
        const [omitted, sent] = fetchMock.mock.calls.map(([, options]) => JSON.parse(options.body).options);
        expect(omitted).toEqual({ num_ctx: 16384, num_predict: 8000 });
        expect(sent).toEqual({ temperature: 0.2, num_ctx: 16384, num_predict: 8000 });
    });

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

        // A model without tool calling: the failure keeps its status and body,
        // which the model port reads as "tools unsupported" (provider/model.js).
        vi.stubGlobal('fetch', vi.fn(async () => json({ error: 'registry.ollama.ai/library/gemma3:latest does not support tools' }, 400)));
        const failure = await collect(streamOllamaChat({ baseUrl: 'http://x:11434', model: 'gemma3', messages: [], tools: [{}] })).catch((error) => error);
        expect(failure.status).toBe(400);
        expect(failure.body).toEqual({ error: 'registry.ollama.ai/library/gemma3:latest does not support tools' });
        expect(failure.message).toBe('HTTP 400: registry.ollama.ai/library/gemma3:latest does not support tools');
        expect(isToolsUnsupported(failure, { server: 'ollama' })).toBe(true);
    });
});
