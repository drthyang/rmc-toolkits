// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { afterEach, describe, expect, it, vi } from 'vitest';
import { checkConnection, completeChat, listModels, streamChat } from '../provider/client';

// Build a Response-like object whose body streams the given chunks, so the SSE
// parser is exercised against realistic reads — including chunks that split a
// `data:` line in the middle.
const streamResponse = (chunks) => {
    const encoder = new TextEncoder();
    const body = new ReadableStream({
        start(controller) {
            chunks.forEach((chunk) => controller.enqueue(encoder.encode(chunk)));
            controller.close();
        }
    });
    return { ok: true, status: 200, body };
};

const sseLine = (content) => `data: ${JSON.stringify({ choices: [{ delta: { content } }] })}\n\n`;

const collect = async (iterator) => {
    const parts = [];
    for await (const part of iterator) parts.push(part);
    return parts;
};

afterEach(() => {
    vi.unstubAllGlobals();
});

describe('streamChat', () => {
    it('yields content deltas and stops at [DONE]', async () => {
        vi.stubGlobal('fetch', vi.fn(async () => streamResponse([
            sseLine('Hel'),
            sseLine('lo'),
            'data: [DONE]\n\n',
            sseLine('never seen')
        ])));
        const parts = await collect(streamChat({ baseUrl: 'http://localhost:11434/v1', model: 'm', messages: [] }));
        expect(parts).toEqual([{ content: 'Hel' }, { content: 'lo' }]);
    });

    it('yields reasoning deltas separately from content', async () => {
        const reasoningLine = (r) => `data: ${JSON.stringify({ choices: [{ delta: { reasoning: r } }] })}\n\n`;
        vi.stubGlobal('fetch', vi.fn(async () => streamResponse([
            reasoningLine('think'),
            reasoningLine('ing'),
            sseLine('answer'),
            'data: [DONE]\n'
        ])));
        const parts = await collect(streamChat({ baseUrl: 'http://x/v1', model: 'm', messages: [] }));
        expect(parts).toEqual([{ reasoning: 'think' }, { reasoning: 'ing' }, { content: 'answer' }]);
    });

    it('reassembles lines split across chunk boundaries', async () => {
        const line = sseLine('split-delta');
        vi.stubGlobal('fetch', vi.fn(async () => streamResponse([
            line.slice(0, 12),
            line.slice(12),
            'data: [DONE]\n'
        ])));
        const parts = await collect(streamChat({ baseUrl: 'http://x/v1', model: 'm', messages: [] }));
        expect(parts).toEqual([{ content: 'split-delta' }]);
    });

    it('skips malformed JSON lines and empty deltas', async () => {
        vi.stubGlobal('fetch', vi.fn(async () => streamResponse([
            'data: {not json}\n',
            `data: ${JSON.stringify({ choices: [{ delta: {} }] })}\n`,
            sseLine('ok'),
            'data: [DONE]\n'
        ])));
        const parts = await collect(streamChat({ baseUrl: 'http://x/v1', model: 'm', messages: [] }));
        expect(parts).toEqual([{ content: 'ok' }]);
    });

    it('throws on an error event sent mid-stream, after the 200', async () => {
        // LM Studio reports a backend crash (here a GPU out-of-memory) this way.
        const oom = 'Encountered fatal exception in the backend generation thread: Insufficient Memory';
        vi.stubGlobal('fetch', vi.fn(async () => streamResponse([
            sseLine('Par'),
            `data: ${JSON.stringify({ error: { message: oom } })}\n\n`,
            sseLine('never seen')
        ])));
        const parts = [];
        await expect((async () => {
            for await (const part of streamChat({ baseUrl: 'http://x/v1', model: 'm', messages: [] })) parts.push(part);
        })()).rejects.toThrow(`The server stopped the reply: ${oom}`);
        expect(parts).toEqual([{ content: 'Par' }]);

        vi.stubGlobal('fetch', vi.fn(async () => streamResponse([`data: ${JSON.stringify({ error: 'model unloaded' })}\n\n`])));
        await expect(collect(streamChat({ baseUrl: 'http://x/v1', model: 'm', messages: [] })))
            .rejects.toThrow('The server stopped the reply: model unloaded');
    });

    it('throws with the server error message on HTTP failure', async () => {
        vi.stubGlobal('fetch', vi.fn(async () => ({
            ok: false,
            status: 404,
            json: async () => ({ error: { message: 'model not found' } })
        })));
        await expect(collect(streamChat({ baseUrl: 'http://x/v1', model: 'm', messages: [] })))
            .rejects.toThrow('HTTP 404: model not found');
    });

    it('sends the expected chat-completions payload', async () => {
        const fetchMock = vi.fn(async () => streamResponse(['data: [DONE]\n']));
        vi.stubGlobal('fetch', fetchMock);
        const messages = [{ role: 'user', content: 'hi' }];
        await collect(streamChat({ baseUrl: 'http://x/v1/', model: 'llama3.2', messages, temperature: 0.7 }));
        const [url, options] = fetchMock.mock.calls[0];
        expect(url).toBe('http://x/v1/chat/completions');
        expect(JSON.parse(options.body)).toEqual({ model: 'llama3.2', messages, temperature: 0.7, stream: true });
        // Non-Anthropic hosts must not receive the Anthropic-only headers.
        expect(options.headers['anthropic-version']).toBeUndefined();
        expect(options.headers['anthropic-dangerous-direct-browser-access']).toBeUndefined();
    });

    it('leaves a null temperature out of the body and adds max_tokens only when a cap is given', async () => {
        const fetchMock = vi.fn(async () => streamResponse(['data: [DONE]\n']));
        vi.stubGlobal('fetch', fetchMock);
        const messages = [{ role: 'user', content: 'hi' }];
        await collect(streamChat({ baseUrl: 'https://api.anthropic.com/v1', model: 'claude-sonnet-5-5', messages, temperature: null, maxTokens: 8192 }));
        await collect(streamChat({ baseUrl: 'http://x/v1', model: 'm', messages }));
        await collect(streamChat({ baseUrl: 'http://x/v1', model: 'm', messages, temperature: 0, maxTokens: undefined }));
        const [capped, plain, zero] = fetchMock.mock.calls.map(([, options]) => options.body);
        // The bytes: no temperature key at all (not null), the cap as max_tokens.
        expect(capped).toBe(JSON.stringify({ model: 'claude-sonnet-5-5', messages, stream: true, max_tokens: 8192 }));
        expect(capped).not.toContain('temperature');
        // Without a cap nothing is added; an absent temperature is the default 0.2, a zero is sent.
        expect(JSON.parse(plain)).toEqual({ model: 'm', messages, temperature: 0.2, stream: true });
        expect(JSON.parse(zero)).toEqual({ model: 'm', messages, temperature: 0, stream: true });
    });

    it('adds Anthropic browser-access and version headers for api.anthropic.com', async () => {
        const fetchMock = vi.fn(async () => streamResponse(['data: [DONE]\n']));
        vi.stubGlobal('fetch', fetchMock);
        await collect(streamChat({
            baseUrl: 'https://api.anthropic.com/v1',
            model: 'claude-opus-4-8',
            messages: [{ role: 'user', content: 'hi' }],
            apiKey: 'sk-ant-test'
        }));
        const [url, options] = fetchMock.mock.calls[0];
        expect(url).toBe('https://api.anthropic.com/v1/chat/completions');
        expect(options.headers.Authorization).toBe('Bearer sk-ant-test');
        expect(options.headers['anthropic-version']).toBe('2023-06-01');
        expect(options.headers['anthropic-dangerous-direct-browser-access']).toBe('true');
    });

    it('sends tools and tool_choice only when tools are given', async () => {
        const fetchMock = vi.fn(async () => streamResponse(['data: [DONE]\n']));
        vi.stubGlobal('fetch', fetchMock);
        const tools = [{ type: 'function', function: { name: 'list_sites', description: 'd', parameters: { type: 'object', properties: {} } } }];
        await collect(streamChat({ baseUrl: 'http://x/v1', model: 'm', messages: [], tools, toolChoice: 'auto' }));
        await collect(streamChat({ baseUrl: 'http://x/v1', model: 'm', messages: [], tools: [], toolChoice: 'auto' }));
        const withTools = JSON.parse(fetchMock.mock.calls[0][1].body);
        expect(withTools.tools).toEqual(tools);
        expect(withTools.tool_choice).toBe('auto');
        const without = JSON.parse(fetchMock.mock.calls[1][1].body);
        expect('tools' in without).toBe(false);
        expect('tool_choice' in without).toBe(false);
    });

    it('yields streamed tool-call fragments and the finish reason', async () => {
        const toolLine = (calls, finish = null) => `data: ${JSON.stringify({
            choices: [{ delta: { tool_calls: calls }, finish_reason: finish }]
        })}\n`;
        const first = toolLine([{ index: 0, id: 'call_1', type: 'function', function: { name: 'compute_bond_angles', arguments: '{"end1":' } }]);
        vi.stubGlobal('fetch', vi.fn(async () => streamResponse([
            first.slice(0, 30),
            first.slice(30),
            toolLine([{ index: 0, function: { arguments: '"Se"}' } }]),
            // Some servers send parsed arguments; they are re-serialized.
            toolLine([{ index: 1, id: 'call_2', function: { name: 'list_sites', arguments: { a: 1 } } }]),
            `data: ${JSON.stringify({ choices: [{ delta: {}, finish_reason: 'tool_calls' }] })}\n`,
            'data: [DONE]\n'
        ])));
        const parts = await collect(streamChat({ baseUrl: 'http://x/v1', model: 'm', messages: [], tools: [{}] }));
        expect(parts).toEqual([
            { toolCall: { index: 0, id: 'call_1', name: 'compute_bond_angles', arguments: '{"end1":', extraContent: null } },
            { toolCall: { index: 0, id: null, name: null, arguments: '"Se"}', extraContent: null } },
            { toolCall: { index: 1, id: 'call_2', name: 'list_sites', arguments: '{"a":1}', extraContent: null } },
            { finishReason: 'tool_calls' }
        ]);
    });

    it('passes a Gemini-style whole call (no index) with its extra_content through', async () => {
        const call = {
            id: 'function-call-1',
            type: 'function',
            function: { name: 'list_sites', arguments: '{}' },
            extra_content: { google: { thought_signature: 'c2ln' } }
        };
        vi.stubGlobal('fetch', vi.fn(async () => streamResponse([
            `data: ${JSON.stringify({ choices: [{ delta: { role: 'assistant', tool_calls: [call] }, finish_reason: 'stop' }] })}\n`,
            'data: [DONE]\n'
        ])));
        const parts = await collect(streamChat({ baseUrl: 'http://x/v1', model: 'm', messages: [], tools: [{}] }));
        expect(parts).toEqual([
            { toolCall: { index: null, id: 'function-call-1', name: 'list_sites', arguments: '{}', extraContent: { google: { thought_signature: 'c2ln' } } } },
            { finishReason: 'stop' }
        ]);
    });

    it('attaches the HTTP status to a failed request', async () => {
        vi.stubGlobal('fetch', vi.fn(async () => ({
            ok: false,
            status: 400,
            json: async () => ({ error: 'llama-tiny does not support tools' })
        })));
        await expect(collect(streamChat({ baseUrl: 'http://x/v1', model: 'm', messages: [], tools: [{}] })))
            .rejects.toMatchObject({ status: 400, message: 'HTTP 400: llama-tiny does not support tools' });
    });
});

describe('completeChat', () => {
    it('returns the message content', async () => {
        vi.stubGlobal('fetch', vi.fn(async () => ({
            ok: true,
            status: 200,
            json: async () => ({ choices: [{ message: { content: 'STATUS: improving — dropping' } }] })
        })));
        const reply = await completeChat({ baseUrl: 'http://x/v1', model: 'm', messages: [] });
        expect(reply).toBe('STATUS: improving — dropping');
    });

    it('rejects when the reply has no content', async () => {
        vi.stubGlobal('fetch', vi.fn(async () => ({ ok: true, status: 200, json: async () => ({}) })));
        await expect(completeChat({ baseUrl: 'http://x/v1', model: 'm', messages: [] }))
            .rejects.toThrow('no message content');
    });
});

describe('listModels / checkConnection', () => {
    it('lists model ids', async () => {
        vi.stubGlobal('fetch', vi.fn(async () => ({
            ok: true,
            status: 200,
            json: async () => ({ data: [{ id: 'llama3.2' }, { id: 'qwen2.5' }] })
        })));
        expect(await listModels('http://x/v1')).toEqual(['llama3.2', 'qwen2.5']);
    });

    it('maps a network TypeError to the server-or-CORS hint', async () => {
        vi.stubGlobal('fetch', vi.fn(async () => {
            throw new TypeError('Failed to fetch');
        }));
        const result = await checkConnection('http://localhost:11434/v1');
        expect(result.ok).toBe(false);
        expect(result.hint).toContain('How to connect a model');
        expect(result.hint).toContain('CORS');
        expect(result.unreachable).toBe(true);
    });

    it('maps HTTP 404 to the /v1 base-URL hint', async () => {
        vi.stubGlobal('fetch', vi.fn(async () => ({ ok: false, status: 404, json: async () => ({}) })));
        const result = await checkConnection('http://localhost:11434');
        expect(result.ok).toBe(false);
        expect(result.hint).toContain('/v1');
    });

    it('reports success with the model list', async () => {
        vi.stubGlobal('fetch', vi.fn(async () => ({
            ok: true,
            status: 200,
            json: async () => ({ data: [{ id: 'llama3.2' }] })
        })));
        expect(await checkConnection('http://x/v1')).toEqual({
            ok: true,
            models: ['llama3.2'],
            error: null,
            hint: null
        });
    });
});
