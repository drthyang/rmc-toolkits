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
        expect(parts).toEqual(['Hel', 'lo']);
    });

    it('reassembles lines split across chunk boundaries', async () => {
        const line = sseLine('split-delta');
        vi.stubGlobal('fetch', vi.fn(async () => streamResponse([
            line.slice(0, 12),
            line.slice(12),
            'data: [DONE]\n'
        ])));
        const parts = await collect(streamChat({ baseUrl: 'http://x/v1', model: 'm', messages: [] }));
        expect(parts).toEqual(['split-delta']);
    });

    it('skips malformed JSON lines and empty deltas', async () => {
        vi.stubGlobal('fetch', vi.fn(async () => streamResponse([
            'data: {not json}\n',
            `data: ${JSON.stringify({ choices: [{ delta: {} }] })}\n`,
            sseLine('ok'),
            'data: [DONE]\n'
        ])));
        const parts = await collect(streamChat({ baseUrl: 'http://x/v1', model: 'm', messages: [] }));
        expect(parts).toEqual(['ok']);
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
        expect(result.hint).toContain('OLLAMA_ORIGINS');
        expect(result.hint).toContain('CORS');
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
