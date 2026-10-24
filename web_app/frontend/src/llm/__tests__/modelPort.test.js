// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The production model port (provider/model.js providerModel) at its
// interface, prepare and stream, over fake servers behind fetch (Ollama, LM
// Studio, a cloud provider): nothing is mocked above fetch and no request
// leaves the test. It pins what the port guarantees:
//   - prepare returns the request as it goes out, with the context guard's
//     plan and the tool choice's form per server (TOOL_CHOICE_FORMS), or
//     throws ContextLimitError before anything is sent;
//   - on Ollama, tool_choice 'none' goes out as no tools (the native API has
//     no tool_choice): the model cannot call one;
//   - stream sends exactly that plan, and a model that cannot take tools
//     fails with a ToolsUnsupportedError, read from the provider's response
//     (its shape first, the fallback text last), and only then;
//   - the loop reacts to that type: the turn's first request is asked again
//     without tools.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { runCopilotTurn } from '../copilot/agentLoop';
import { ContextLimitError, estimatePromptTokens, forgetServers } from '../provider/contextGuard';
import {
    ContextLimitError as PortContextLimitError, TOOL_CHOICE_FORMS, ToolsUnsupportedError, isToolsUnsupported, providerModel, toolChoiceForm
} from '../provider/model';

const encoder = new TextEncoder();
const streamOf = (text) => new ReadableStream({
    start(controller) {
        controller.enqueue(encoder.encode(text));
        controller.close();
    }
});
const json = (body, status = 200) => ({ ok: status < 400, status, json: async () => body });
const sse = (text) => ({ ok: true, status: 200, body: streamOf(`data: ${JSON.stringify({ choices: [{ delta: { content: text } }] })}\n\n`
    + `data: ${JSON.stringify({ choices: [{ delta: {}, finish_reason: 'stop' }] })}\n\ndata: [DONE]\n\n`) });
const ndjson = (lines) => ({ ok: true, status: 200, body: streamOf(lines.map((line) => `${JSON.stringify(line)}\n`).join('')) });
const ollamaText = (text) => ndjson([
    { message: { role: 'assistant', content: text }, done: false },
    { message: { role: 'assistant', content: '' }, done: true, done_reason: 'stop', prompt_eval_count: 4000, eval_count: 9 }
]);
const ollamaCall = (name, args) => ndjson([
    { message: { role: 'assistant', content: '', tool_calls: [{ function: { name, arguments: args } }] }, done: false },
    { message: { role: 'assistant', content: '' }, done: true, done_reason: 'stop', prompt_eval_count: 4000, eval_count: 9 }
]);

/**
 * A fake server: 'ollama' (/api/version, /api/ps, /api/show, /api/chat),
 * 'lmstudio' (/api/v0/models, /v1/chat/completions) or 'cloud' (only
 * /chat/completions; a probe fails the test). `chat` answers each chat
 * request in order: a response, or a function of the parsed body.
 */
const fakeServer = (kind, { max = 262144, loaded = 49152, chat = [] } = {}) => {
    const calls = [];
    const queue = [...chat];
    vi.stubGlobal('fetch', vi.fn(async (url, options = {}) => {
        const { pathname } = new URL(url);
        const body = options.body ? JSON.parse(options.body) : null;
        calls.push({ path: pathname, body, signal: options.signal });
        if (kind === 'cloud' && pathname.startsWith('/api/')) throw new Error(`a cloud provider was probed: ${url}`);
        if (pathname === '/api/version') return kind === 'ollama' ? json({ version: '0.40.2' }) : json({ error: 'Unexpected endpoint' });
        if (pathname === '/api/v0/models') return json({ data: [{ id: 'm', state: 'loaded', loaded_context_length: loaded }] });
        if (pathname === '/api/ps') return json({ models: [] });
        if (pathname === '/api/show') return json({ model_info: { 'm.context_length': max } });
        if (pathname === '/api/chat' || pathname.endsWith('/chat/completions')) {
            const next = queue.shift();
            if (!next) return kind === 'ollama' ? ollamaText('Done.') : sse('Done.');
            return typeof next === 'function' ? next(body) : next;
        }
        return json({ error: 'not found' }, 404);
    }));
    return { calls, chats: () => calls.filter((call) => call.path === '/api/chat' || call.path.endsWith('/chat/completions')) };
};

const URLS = { ollama: 'http://localhost:11434/v1', lmstudio: 'http://localhost:1234/v1', cloud: 'https://api.openai.com/v1' };
const portOn = (kind, session = {}) => providerModel({ baseUrl: URLS[kind], model: 'm', temperature: 0.2, apiKey: kind === 'cloud' ? 'test-key' : '' }, { session });
const TOOLS = [{ type: 'function', function: { name: 'echo', description: 'Echo.', parameters: { type: 'object', properties: { value: { type: 'number' } } } } }];
const MESSAGES = [{ role: 'system', content: 'Rules.' }, { role: 'user', content: 'Echo 1.' }];
const intent = (toolChoice, extra = {}) => ({ messages: MESSAGES, tools: TOOLS, toolChoice, retry: null, ...extra });
const collect = async (iterator) => {
    const parts = [];
    for await (const part of iterator) parts.push(part);
    return parts;
};

beforeEach(() => forgetServers());
afterEach(() => vi.unstubAllGlobals());

describe('prepare: the request as it goes out, with the guard\'s plan', () => {
    it('a cloud provider: never probed; the tool choice in the body, beside the tools, in every round', async () => {
        const { calls } = fakeServer('cloud');
        const port = portOn('cloud');
        for (const choice of ['auto', 'none']) {
            expect(await port.prepare(intent(choice))).toEqual({
                baseUrl: URLS.cloud, model: 'm', temperature: 0.2, messages: MESSAGES, tools: TOOLS, toolChoice: choice,
                plan: { server: 'openai', toolChoice: 'body' }
            });
        }
        expect(calls).toEqual([]);
        // The key joins only as the request goes out: never in `sent`.
        expect(JSON.stringify(await port.prepare(intent('auto')))).not.toContain('test-key');
    });

    it('LM Studio: the loaded window in the plan; the tool choice in the body', async () => {
        fakeServer('lmstudio', { loaded: 80000 });
        const sent = await portOn('lmstudio').prepare(intent('none'));
        const estimate = estimatePromptTokens({ messages: MESSAGES, tools: TOOLS });
        expect(sent.tools).toBe(TOOLS);
        expect(sent.toolChoice).toBe('none');
        expect(sent.plan).toEqual({ server: 'lmstudio', estimate, needed: estimate + 4096, loadedContext: 80000, toolChoice: 'body' });
    });

    it('Ollama: num_ctx and num_predict in the plan; \'auto\' omitted (its default), \'none\' as no tools', async () => {
        fakeServer('ollama');
        const port = portOn('ollama');
        const auto = await port.prepare(intent('auto'));
        expect(auto.tools).toBe(TOOLS);
        expect(auto.plan).toMatchObject({ server: 'ollama', version: '0.40.2', loadedContext: null, maxContext: 262144, toolChoice: 'omitted' });
        expect(auto.plan.numPredict).toBe(auto.plan.numCtx - auto.plan.estimate);
        const final = await port.prepare(intent('none'));
        expect(final.tools).toBeUndefined();
        expect(final.toolChoice).toBe('none');
        expect(final.plan.toolChoice).toBe('tools_withheld');
        // The estimate counts what is sent: no tool schemas.
        expect(final.plan.estimate).toBe(estimatePromptTokens({ messages: MESSAGES }));
        expect(final.plan.estimate).toBeLessThan(auto.plan.estimate);
    });

    it('a request without tools has no tool-choice form', async () => {
        fakeServer('ollama');
        const sent = await portOn('ollama').prepare({ messages: MESSAGES, tools: undefined, toolChoice: undefined, retry: null });
        expect(sent.plan.server).toBe('ollama');
        expect('toolChoice' in sent.plan).toBe(false);
    });

    it('refuses before anything is sent: ContextLimitError from prepare (Ollama\'s maximum, LM Studio\'s loaded window)', async () => {
        const ollama = fakeServer('ollama', { max: 8192 });
        const big = intent('auto', { messages: [...MESSAGES, { role: 'user', content: 'x'.repeat(30000) }] });
        await expect(portOn('ollama').prepare(big)).rejects.toBeInstanceOf(ContextLimitError);
        expect(ollama.chats()).toEqual([]);
        forgetServers();
        const studio = fakeServer('lmstudio', { loaded: 4096 });
        await expect(portOn('lmstudio').prepare(big)).rejects.toThrow(/LM Studio runs m with a 4,096-token context, .* Nothing was sent/);
        expect(studio.chats()).toEqual([]);
        // One error type, re-exported by the port.
        expect(PortContextLimitError).toBe(ContextLimitError);
    });

    it('a reply cap (maxTokens): in the request as sent; on Ollama a ceiling on the plan\'s numPredict; none by default', async () => {
        fakeServer('ollama');
        const capped = providerModel({ baseUrl: URLS.ollama, model: 'm', temperature: 0.2, maxTokens: 4096 }, { session: {} });
        const sent = await capped.prepare(intent('auto'));
        expect(sent.maxTokens).toBe(4096);
        expect(sent.plan.numPredict).toBe(Math.min(4096, sent.plan.numCtx - sent.plan.estimate));
        expect(sent.plan.numPredict).toBe(4096);
        forgetServers();
        fakeServer('ollama');
        const plain = await portOn('ollama').prepare(intent('auto'));
        expect('maxTokens' in plain).toBe(false);
        expect(plain.plan.numPredict).toBe(plain.plan.numCtx - plain.plan.estimate);
        // A cap that is not a positive integer is refused when the port is made, never dropped.
        for (const bad of [0, -1, 8192.5, '8192', Number.NaN]) {
            expect(() => providerModel({ baseUrl: URLS.cloud, model: 'm', maxTokens: bad })).toThrow(RangeError);
        }
    });

    it('the forms, by server kind', () => {
        expect(['openai', 'lmstudio', 'ollama'].map((kind) => [toolChoiceForm(kind, 'auto'), toolChoiceForm(kind, 'none')])).toEqual([
            ['body', 'body'], ['body', 'body'], ['omitted', 'tools_withheld']
        ]);
        expect(Object.keys(TOOL_CHOICE_FORMS)).toEqual(['body', 'omitted', 'tools_withheld']);
    });
});

describe('stream: sends the plan', () => {
    it('Ollama: the plan\'s num_ctx and num_predict; the final request carries no tools and no tool_choice', async () => {
        const { chats } = fakeServer('ollama');
        const session = {};
        const port = portOn('ollama', session);
        const sent = await port.prepare(intent('none'));
        expect(await collect(port.stream(sent))).toContainEqual({ content: 'Done.' });
        const [chat] = chats();
        expect(chat.path).toBe('/api/chat');
        expect('tools' in chat.body).toBe(false);
        expect('tool_choice' in chat.body).toBe(false);
        expect(chat.body.options).toMatchObject({ num_ctx: sent.plan.numCtx, num_predict: sent.plan.numPredict });
        expect(session).toEqual({ 'http://localhost:11434/v1|m': sent.plan.numCtx });
    });

    it('Ollama with a reply cap: the capped num_predict is what the body carries', async () => {
        const { chats } = fakeServer('ollama');
        const port = providerModel({ baseUrl: URLS.ollama, model: 'm', temperature: 0.2, maxTokens: 4096 }, { session: {} });
        const sent = await port.prepare(intent('auto'));
        await collect(port.stream(sent));
        expect(chats()[0].body.options).toEqual({ temperature: 0.2, num_ctx: sent.plan.numCtx, num_predict: 4096 });
        expect('max_tokens' in chats()[0].body).toBe(false);
    });

    it('a model that takes no temperature (null): none in the request as sent or the body, a warmer retry included', async () => {
        for (const kind of ['cloud', 'lmstudio']) {
            forgetServers();
            const { chats } = fakeServer(kind);
            const port = providerModel({ baseUrl: URLS[kind], model: 'm', temperature: null, apiKey: 'test-key' }, { session: {} });
            for (const retry of [null, 'warmer']) {
                const sent = await port.prepare(intent('auto', { retry }));
                expect('temperature' in sent).toBe(false);
                await collect(port.stream(sent));
            }
            expect(chats().map((chat) => chat.body)).toEqual([0, 1].map(() => ({ model: 'm', messages: MESSAGES, stream: true, tools: TOOLS, tool_choice: 'auto' })));
        }
        forgetServers();
        const { chats } = fakeServer('ollama');
        const port = providerModel({ baseUrl: URLS.ollama, model: 'm', temperature: null }, { session: {} });
        await collect(port.stream(await port.prepare(intent('auto', { retry: 'warmer' }))));
        expect('temperature' in chats()[0].body.options).toBe(false);
    });

    it('LM Studio and the cloud: the body as before, tool_choice "none" beside the tools', async () => {
        for (const kind of ['lmstudio', 'cloud']) {
            forgetServers();
            const { chats } = fakeServer(kind);
            const port = portOn(kind);
            await collect(port.stream(await port.prepare(intent('none'))));
            expect(chats()[0].body).toEqual({ model: 'm', messages: MESSAGES, temperature: 0.2, stream: true, tools: TOOLS, tool_choice: 'none' });
        }
    });
});

describe('stream: "tools unsupported" as a typed failure', () => {
    const failing = (kind, response) => {
        fakeServer(kind, { chat: [response] });
        const port = portOn(kind);
        return port.prepare(intent('auto')).then((sent) => collect(port.stream(sent))).catch((error) => error);
    };

    it('Ollama: a 400 whose error says the model does not support tools', async () => {
        const failure = await failing('ollama', json({ error: 'registry.ollama.ai/library/gemma3:latest does not support tools' }, 400));
        expect(failure).toBeInstanceOf(ToolsUnsupportedError);
        expect(failure).toMatchObject({
            name: 'ToolsUnsupportedError', status: 400, server: 'ollama',
            message: 'HTTP 400: registry.ollama.ai/library/gemma3:latest does not support tools'
        });
        expect(failure.cause.body).toEqual({ error: 'registry.ollama.ai/library/gemma3:latest does not support tools' });
    });

    it('OpenAI: a 400 naming the tools parameter as unsupported, by its shape (the text alone would not tell)', async () => {
        const body = { error: { message: 'Unsupported parameter: \'tools\' is not supported with this model.', type: 'invalid_request_error', param: 'tools', code: 'unsupported_parameter' } };
        const failure = await failing('cloud', json(body, 400));
        expect(failure).toBeInstanceOf(ToolsUnsupportedError);
        expect(failure).toMatchObject({ status: 400, server: 'openai' });
        expect(isToolsUnsupported(Object.assign(new Error(`HTTP 400: ${body.error.message}`), { status: 400 }))).toBe(false);
    });

    it('LM Studio (no shape known): the fallback text', async () => {
        const failure = await failing('lmstudio', json({ error: 'tools are not supported by this model' }, 400));
        expect(failure).toBeInstanceOf(ToolsUnsupportedError);
        expect(failure.server).toBe('lmstudio');
    });

    it('never for a schema error, a 5xx, a request without tools, a failure after the reply began, or an abort', async () => {
        const schema = await failing('cloud', json({ error: { message: 'Invalid schema for function \'echo\'', param: 'tools[0].function.parameters', code: 'invalid_function_parameters' } }, 400));
        expect(schema).not.toBeInstanceOf(ToolsUnsupportedError);
        expect(schema).toMatchObject({ status: 400, message: 'HTTP 400: Invalid schema for function \'echo\'' });
        forgetServers();
        expect(await failing('ollama', json({ error: 'gemma3 does not support tools' }, 500))).not.toBeInstanceOf(ToolsUnsupportedError);

        forgetServers();
        fakeServer('ollama', { chat: [json({ error: 'gemma3 does not support tools' }, 400)] });
        const plain = portOn('ollama');
        const noTools = await plain.prepare({ messages: MESSAGES, tools: undefined, toolChoice: undefined, retry: null })
            .then((sent) => collect(plain.stream(sent))).catch((error) => error);
        expect(noTools).not.toBeInstanceOf(ToolsUnsupportedError);
        expect(noTools.status).toBe(400);

        forgetServers();
        const midway = await failing('cloud', { ok: true, status: 200, body: streamOf(`data: ${JSON.stringify({ choices: [{ delta: { content: 'Par' } }] })}\n\n`
            + `data: ${JSON.stringify({ error: { message: 'tools are not supported' } })}\n\n`) });
        expect(midway).not.toBeInstanceOf(ToolsUnsupportedError);
        expect(midway.message).toBe('The server stopped the reply: tools are not supported');

        forgetServers();
        fakeServer('cloud', { chat: [async () => { throw new DOMException('Aborted', 'AbortError'); }] });
        const port = portOn('cloud');
        const aborted = await collect(port.stream(await port.prepare(intent('auto')))).catch((error) => error);
        expect(aborted.name).toBe('AbortError');
    });

    it('isToolsUnsupported: the shapes per server, the fallback text, never outside a 4xx', () => {
        const err = (message, status = 400, body = null) => Object.assign(new Error(message), { status, body });
        expect(isToolsUnsupported(err('HTTP 400: m does not support tools', 400, { error: 'm does not support tools' }), { server: 'ollama' })).toBe(true);
        expect(isToolsUnsupported(err('HTTP 400', 400, { error: { param: 'tools', code: 'unsupported_parameter' } }), { server: 'openai' })).toBe(true);
        // A shape is read only on its own server.
        expect(isToolsUnsupported(err('HTTP 400', 400, { error: { param: 'tools', code: 'unsupported_parameter' } }), { server: 'ollama' })).toBe(false);
        // The fallback text, on any server.
        expect(isToolsUnsupported(err('HTTP 400: tools are not supported by this model'))).toBe(true);
        expect(isToolsUnsupported(err('HTTP 400: function calling is not supported'))).toBe(true);
        expect(isToolsUnsupported(err('HTTP 500: model does not support tools', 500))).toBe(false);
        expect(isToolsUnsupported(new Error('model does not support tools'))).toBe(false);
        // A schema problem is a real error, never "this model has no tools".
        expect(isToolsUnsupported(err("HTTP 400: Invalid schema for function 'list_sites'"))).toBe(false);
        expect(isToolsUnsupported(err('HTTP 400: tools.0.function.parameters is invalid'))).toBe(false);
        expect(isToolsUnsupported(new ToolsUnsupportedError('typed'))).toBe(true);
    });
});

describe('the loop over the port (fake servers)', () => {
    const echo = {
        name: 'echo',
        description: 'Echo.',
        parameters: { type: 'object', properties: { value: { type: 'number' } }, required: ['value'] },
        run: async ({ value }) => ({ data: { value }, checks: [] })
    };

    it('Ollama: the final round is sent without tools, so its reply is the answer (behaviour change of the model port)', async () => {
        const { chats } = fakeServer('ollama', { chat: [ollamaCall('echo', { value: 1 }), ollamaText('The value is 1.\n\nOutcome: achieved.')] });
        const events = [];
        for await (const event of runCopilotTurn({ model: portOn('ollama'), messages: MESSAGES, tools: [echo], limits: { maxRounds: 1 } })) events.push(event);
        const requests = events.filter((event) => event.type === 'request');
        expect(requests.map((event) => [event.why, event.sent.plan.toolChoice, Boolean(event.sent.tools)])).toEqual([
            ['first', 'omitted', true], ['final', 'tools_withheld', false]
        ]);
        expect(chats().map((chat) => 'tools' in chat.body)).toEqual([true, false]);
        // The history's call and result still travel as messages.
        expect(chats()[1].body.messages.map((message) => message.role)).toEqual(['system', 'user', 'assistant', 'tool', 'user']);
        expect(chats()[1].body.messages.at(-1).content).toMatch(/used the analysis budget/);
        expect(events.at(-1)).toMatchObject({ type: 'done', stoppedBy: 'max_rounds', text: 'The value is 1.\n\nOutcome: achieved.' });
    });

    it('a model without tools: the port\'s typed failure, then the first request asked again without tools', async () => {
        const noTools = (body) => (body.tools ? json({ error: 'registry.ollama.ai/library/gemma3:latest does not support tools' }, 400) : ollamaText('Plain answer.'));
        const { chats } = fakeServer('ollama', { chat: [noTools, noTools] });
        const events = [];
        for await (const event of runCopilotTurn({ model: portOn('ollama'), messages: MESSAGES, tools: [echo], noToolsSystemPrompt: 'PLAIN' })) events.push(event);
        expect(events.filter((event) => ['request', 'tools_unsupported', 'done'].includes(event.type)).map((event) => event.type))
            .toEqual(['request', 'tools_unsupported', 'request', 'done']);
        expect(events.find((event) => event.type === 'tools_unsupported').error).toMatch(/does not support tools/);
        expect(chats().map((chat) => 'tools' in chat.body)).toEqual([true, false]);
        expect(events.at(-1)).toMatchObject({ type: 'done', text: 'Plain answer.' });
    });
});
