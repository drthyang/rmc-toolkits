// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The run policy (plan item 1.2, validation/runPolicy.js), hermetic: every
// server is a fake fetch (LM Studio, Ollama, the cloud providers), and no
// request leaves the test. It pins:
//   - the pre-run check: an LM Studio model loaded below the minimum (48k),
//     not loaded, or a server of another kind or none, is refused before
//     anything is sent, with what was found, what is needed and how to reload;
//     Ollama is not checked, only recorded; a cloud provider is never probed;
//   - what a record, its transcript and the report gain: the policy, the
//     window the run ran under (LM Studio's loaded context before and after,
//     the num_ctx each Ollama request was sent), and the columns side by side;
//   - the scripted check of the final round: what the model port sends each
//     provider for tool_choice 'none' (FINAL_ROUND_FORMS): the field beside
//     the tools on the OpenAI-compatible servers, no tools on Ollama (since
//     protocol 2026-10-19.7). Acceptance of the field by each provider is not
//     testable without a live call; the doc says which are known to honour it.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { DEFAULT_LIMITS, runCopilotTurn } from '../../src/llm/copilot/agentLoop';
import * as guard from '../../src/llm/provider/contextGuard';
import { ANTI_REPEAT_TEMPERATURE, providerModel } from '../../src/llm/provider/model';
import { PROVIDER_PRESETS } from '../../src/llm/provider/presets';
import { buildCase, DEV_CASES } from '../cases';
import { renderReport } from '../report';
import { liveRun } from '../run';
import { MODELS, resolveModel } from '../models';
import {
    buildRunPolicy, CONTEXT_GUARD, FINAL_ROUND_FORMS, MIN_CONTEXT_TOKENS, noteWindowAfter, preflight, RunPolicyError
} from '../runPolicy';
import { renderTranscript } from '../transcript';

const encoder = new TextEncoder();
const streamOf = (text) => new ReadableStream({
    start(controller) {
        controller.enqueue(encoder.encode(text));
        controller.close();
    }
});
const json = (body, status = 200) => ({ ok: status < 400, status, json: async () => body });
const streamed = (text) => ({ ok: true, status: 200, body: streamOf(text) });

// One reply in each dialect: { call: [name, args] } or { text }.
const sse = (reply) => {
    const events = reply.call
        ? [{ choices: [{ delta: { tool_calls: [{ index: 0, id: 'c1', function: { name: reply.call[0], arguments: JSON.stringify(reply.call[1]) } }] } }] },
            { choices: [{ delta: {}, finish_reason: 'tool_calls' }] }]
        : [{ choices: [{ delta: { content: reply.text } }] }, { choices: [{ delta: {}, finish_reason: 'stop' }] }];
    return `${events.map((event) => `data: ${JSON.stringify(event)}\n\n`).join('')}data: [DONE]\n\n`;
};
const ndjson = (reply) => [
    reply.call
        ? { message: { role: 'assistant', content: '', tool_calls: [{ function: { name: reply.call[0], arguments: reply.call[1] } }] }, done: false }
        : { message: { role: 'assistant', content: reply.text }, done: false },
    { message: { role: 'assistant', content: '' }, done: true, done_reason: 'stop' }
].map((line) => `${JSON.stringify(line)}\n`).join('');

/**
 * A fake server behind fetch. `kind`: 'lmstudio' (answers /api/v0/models with
 * the model's load state), 'ollama' (/api/version, /api/ps, /api/show,
 * /api/chat; a request loads the model at its num_ctx, as Ollama does), or
 * 'cloud' (only /chat/completions; any probe fails the test). `replies` are
 * served in order to the chat requests. `down`: every request fails as a
 * refused connection does.
 */
const fakeServer = ({ kind, id = 'qwen3.8-27b', loaded = 49152, state = 'loaded', max = 262144, replies = [], down = false }) => {
    const calls = [];
    const chats = [];
    const queue = [...replies];
    const ollama = { loaded };
    vi.stubGlobal('fetch', vi.fn(async (url, options = {}) => {
        if (down) throw new TypeError('fetch failed');
        const { pathname } = new URL(url);
        const body = options.body ? JSON.parse(options.body) : null;
        calls.push({ url, path: pathname, body });
        if (kind === 'cloud' && pathname.startsWith('/api/')) throw new Error(`a cloud provider was probed: ${url}`);
        if (pathname === '/api/version') return kind === 'ollama' ? json({ version: '0.40.2' }) : json({ error: 'not found' }, 404);
        if (pathname === '/api/v0/models') {
            return kind === 'lmstudio'
                ? json({ data: [{ id, state, ...(state === 'loaded' ? { loaded_context_length: loaded } : {}) }] })
                : json({ error: 'not found' }, 404);
        }
        if (pathname === '/api/ps') return json({ models: ollama.loaded ? [{ name: id, model: id, context_length: ollama.loaded }] : [] });
        if (pathname === '/api/show') return json({ model_info: { 'general.architecture': 'gemma4', 'gemma4.context_length': max } });
        if (pathname === '/api/chat' || pathname.endsWith('/chat/completions')) {
            chats.push({ url, path: pathname, body });
            const reply = queue.shift() ?? { text: 'Done.' };
            if (pathname === '/api/chat') {
                ollama.loaded = body.options?.num_ctx ?? ollama.loaded;
                return streamed(ndjson(reply));
            }
            return streamed(sse(reply));
        }
        return json({ error: 'not found' }, 404);
    }));
    return { calls, chats };
};

beforeEach(() => guard.forgetServers());
afterEach(() => vi.unstubAllGlobals());

const qwen = () => resolveModel('qwen3.8-27b', {});
const gemma = () => ({ ...resolveModel('gemma4:26b', {}) });
const built = (id) => buildCase(DEV_CASES.find((item) => item.id === id));
// An answer that concludes nothing, so the evidence gate sends nothing back.
const ANSWER = 'The Nb site holds its average position.\n\nOutcome: answered.';

describe('the policy\'s constants', () => {
    it('the LM Studio minimum is 48k (48 × 1,024, six context steps); Ollama and the cloud have none', () => {
        expect(MIN_CONTEXT_TOKENS).toEqual({ lmstudio: 49152, ollama: null, openai: null });
        expect(MIN_CONTEXT_TOKENS.lmstudio).toBe(48 * 1024);
        expect(MIN_CONTEXT_TOKENS.lmstudio % guard.CONTEXT_STEP).toBe(0);
    });

    it('states the guard\'s decisions as the guard has them, unchanged', () => {
        expect(CONTEXT_GUARD).toMatchObject({
            on_overflow: 'refuse',
            chars_per_token: guard.CHARS_PER_TOKEN,
            json_chars_per_token: guard.JSON_CHARS_PER_TOKEN,
            message_overhead_tokens: guard.MESSAGE_OVERHEAD_TOKENS,
            reply_reserve_tokens: guard.REPLY_RESERVE_TOKENS,
            lmstudio_reply_min_tokens: guard.LMSTUDIO_REPLY_MIN_TOKENS,
            context_step: guard.CONTEXT_STEP,
            min_context: guard.MIN_CONTEXT,
            read_share_limit: guard.READ_SHARE_LIMIT
        });
        // The decision of 2026-10-19 keeps them: 2.5 / 1.4 characters a token, 8,192 reserved.
        expect([guard.CHARS_PER_TOKEN, guard.JSON_CHARS_PER_TOKEN, guard.REPLY_RESERVE_TOKENS]).toEqual([2.5, 1.4, 8192]);
    });

    it('every model names its server kind, and every kind has a minimum, a sizing rule, a reply cap and a final-round form', () => {
        expect(Object.fromEntries(Object.entries(MODELS).map(([name, model]) => [name, model.server]))).toEqual({
            'qwen3.8-27b': 'lmstudio', 'gemma4:26b': 'ollama', 'meta/muse-glimmer': 'lmstudio', 'claude-sonnet-5.5': 'openai'
        });
        for (const kind of ['lmstudio', 'ollama', 'openai']) {
            const policy = buildRunPolicy({ model: qwen(), server: { kind, version: null, loaded_tokens: null, max_tokens: null } });
            expect(policy.context.sizing).toEqual(expect.any(String));
            expect(policy.reply_cap.rule).toEqual(expect.any(String));
            expect(policy.final_round).toEqual(FINAL_ROUND_FORMS[kind]);
        }
    });

    it('a policy keeps the limits as run, the temperature and its warmer floor; an API model its spend cap', () => {
        const policy = buildRunPolicy({ model: qwen(), limits: { maxRounds: 3 } });
        expect(policy.limits).toEqual({ ...DEFAULT_LIMITS, maxRounds: 3 });
        expect(policy.temperature).toEqual({ base: 0.2, warmer_retry_min: ANTI_REPEAT_TEMPERATURE });
        expect(policy.server).toEqual({ kind: null, expected: 'lmstudio', version: null, base_url: 'http://localhost:1234/v1' });
        expect(policy.history).toEqual({ earlier_turns: 0, reasoning_resent: false });
        expect(policy.spend).toBeNull();
        const sonnet = buildRunPolicy({ model: resolveModel('claude-sonnet-5.5', { ANTHROPIC_API_KEY: 'k' }) });
        expect(sonnet.spend).toEqual({ cap_usd: 12, priced_reply_tokens: 8192 });
        expect(JSON.stringify(sonnet)).not.toContain('"k"');
    });
});

describe('the pre-run check', () => {
    it('lets an LM Studio model loaded at 48k or more run, and reads its window as the guard does', async () => {
        const { calls } = fakeServer({ kind: 'lmstudio', loaded: 49152 });
        await expect(preflight(qwen())).resolves.toEqual({ kind: 'lmstudio', version: null, loaded_tokens: 49152, max_tokens: null });
        expect(calls.map((call) => call.path)).toEqual(['/api/version', '/api/v0/models', '/api/v0/models']);
        guard.forgetServers();
        fakeServer({ kind: 'lmstudio', loaded: 131072 });
        await expect(preflight(qwen())).resolves.toMatchObject({ loaded_tokens: 131072 });
    });

    it('refuses an LM Studio model loaded below 48k: what was found, what is needed, how to reload it', async () => {
        const { chats } = fakeServer({ kind: 'lmstudio', loaded: 32768 });
        const failure = await preflight(qwen()).catch((error) => error);
        expect(failure).toBeInstanceOf(RunPolicyError);
        expect(failure.message).toBe(
            'LM Studio runs qwen3.8-27b with a 32,768-token context, but the run policy needs at least 49,152 (48k) for an '
            + 'LM Studio model, so the run was not started. In LM Studio, eject qwen3.8-27b and load it again with a Context '
            + 'Length of 49,152 or more (its load settings), or from a terminal: lms unload qwen3.8-27b && lms load qwen3.8-27b '
            + '--context-length 49152. Then run the same command again: the runs already recorded are kept.'
        );
        expect(failure).toMatchObject({ needed: 49152, server: { kind: 'lmstudio', loaded_tokens: 32768 } });
        expect(chats).toEqual([]);
    });

    it('refuses an LM Studio model that is not loaded: a request would load it at LM Studio\'s default', async () => {
        fakeServer({ kind: 'lmstudio', state: 'not-loaded' });
        await expect(preflight(qwen())).rejects.toThrow(
            /^LM Studio has no model loaded as 'qwen3\.8-27b' \(the id the harness sends\).*lms load qwen3\.8-27b --context-length 49152/
        );
        guard.forgetServers();
        fakeServer({ kind: 'lmstudio', id: 'qwen/qwen3.8-27b' });
        await expect(preflight(qwen())).rejects.toThrow(/no model loaded as 'qwen3\.8-27b'.*lms ps lists the loaded models/);
    });

    it('refuses when no server, or a server of another kind, answers at the model\'s URL', async () => {
        fakeServer({ kind: 'lmstudio', down: true });
        await expect(preflight(qwen())).rejects.toThrow(
            'qwen3.8-27b runs on LM Studio at http://localhost:1234/v1, but no LM Studio server answered there, so the run was '
            + 'not started. Start LM Studio\'s server and load qwen3.8-27b, then run the same command again.'
        );
        guard.forgetServers();
        fakeServer({ kind: 'ollama' });
        await expect(preflight(qwen())).rejects.toThrow(/but the server there is Ollama, so the run was not started/);
        guard.forgetServers();
        fakeServer({ kind: 'lmstudio', down: true });
        await expect(preflight(gemma())).rejects.toThrow(/^gemma4:26b runs on Ollama at http:\/\/localhost:11434\/v1, but no Ollama server answered there/);
    });

    it('does not check an Ollama model\'s window (the guard sizes num_ctx), only reads it', async () => {
        fakeServer({ kind: 'ollama', id: 'gemma4:26b', loaded: 4096, max: 262144 });
        await expect(preflight(gemma())).resolves.toEqual({ kind: 'ollama', version: '0.40.2', loaded_tokens: 4096, max_tokens: 262144 });
    });

    it('never probes a cloud provider', async () => {
        const { calls } = fakeServer({ kind: 'cloud' });
        await expect(preflight(resolveModel('claude-sonnet-5.5', { ANTHROPIC_API_KEY: 'k' })))
            .resolves.toEqual({ kind: 'openai', version: null, loaded_tokens: null, max_tokens: null });
        expect(calls).toEqual([]);
    });
});

describe('a live run under the policy (run.js liveRun, fake servers)', () => {
    it('LM Studio below 48k: refused before any request, and no record', async () => {
        const { chats } = fakeServer({ kind: 'lmstudio', loaded: 32768, replies: [{ text: ANSWER }] });
        const model = qwen();
        await expect(liveRun({ built: built('dev-offcentre'), arm: 'C', model, runIndex: 0, base: providerModel(model, { session: {} }), session: {} }))
            .rejects.toBeInstanceOf(RunPolicyError);
        expect(chats).toEqual([]);
    });

    it('LM Studio at 48k: the record keeps the policy and the window, read before and after; the transcript and the report show it', async () => {
        const { chats } = fakeServer({ kind: 'lmstudio', loaded: 49152, replies: [{ call: ['list_sites', {}] }, { text: ANSWER }] });
        const model = qwen();
        const session = {};
        const record = await liveRun({ built: built('dev-offcentre'), arm: 'C', model, runIndex: 0, base: providerModel(model, { session }), session });
        expect(record.error).toBeNull();
        expect(chats.map((chat) => chat.url)).toEqual(['http://localhost:1234/v1/chat/completions', 'http://localhost:1234/v1/chat/completions']);
        expect(record.policy).toMatchObject({
            server: { kind: 'lmstudio', expected: 'lmstudio', version: null, base_url: 'http://localhost:1234/v1' },
            temperature: { base: 0.2, warmer_retry_min: 0.6 },
            limits: DEFAULT_LIMITS,
            context: { min_tokens: 49152, loaded_tokens: 49152, loaded_tokens_after: 49152, model_max_tokens: null, num_ctx_used: [] },
            reply_cap: { sent: null },
            final_round: FINAL_ROUND_FORMS.lmstudio,
            context_guard: CONTEXT_GUARD,
            warnings: []
        });
        // No num_ctx is sent to LM Studio.
        expect(record.requests.map((request) => request.num_ctx)).toEqual([null, null]);

        const transcript = renderTranscript(record);
        expect(transcript).toContain('- Server: LM Studio · context window: 49,152 loaded (minimum 49,152)');
        expect(transcript).toContain('- Reply cap: none sent; the guard needs 4,096 tokens of room beyond the prompt · final round: tool_choice "none", tools still sent');
        expect(transcript).toContain('- Temperature 0.2 (a request asked again after a loop: 0.6 or more) · limits: 6 rounds × 4 calls, 600 s a turn');
        expect(transcript).toContain('- Context guard: refuses a request that would not fit, never cuts it (2.5 characters a token, 1.4 for tool JSON; 8,192-token reply reserve)');

        const report = renderReport([record]);
        expect(report).toContain('## What a comparison across models also compares');
        expect(report).toMatch(/\| qwen3\.8-27b · C \| LM Studio \| 49,152 loaded \(minimum 49,152\) \| none sent; .* \| tool_choice "none" \(tools sent\) \| 0\.2 \(retry ≥ 0\.6\) \| 6 rounds × 4 calls, .* \| 0 \|/);
        expect(report).toContain('The limits, the temperature and the context guard are the same in every run with a policy.');
        expect(report).toContain('- the tokenizer: the guard\'s estimate of a request is the same for every model');
    });

    it('LM Studio reloaded during the run: the window read after differs, and the record says so', async () => {
        fakeServer({ kind: 'lmstudio', loaded: 49152, replies: [{ text: ANSWER }] });
        const model = qwen();
        const session = {};
        const record = await liveRun({ built: built('dev-offcentre'), arm: 'C', model, runIndex: 0, base: providerModel(model, { session }), session });
        // As if LM Studio had loaded the model again at its default context.
        noteWindowAfter(record.policy, { kind: 'lmstudio', loaded_tokens: 4096 });
        expect(record.policy.warnings).toEqual([expect.stringMatching(/^LM Studio's loaded context changed during the run: 49,152 tokens before, 4,096 after/)]);
        expect(renderTranscript(record)).toContain('- **Run policy warning:** LM Studio\'s loaded context changed during the run');
        expect(renderReport([record])).toMatch(/49,152 loaded, 4,096 after \(minimum 49,152\) \|.*\| 1 \|$/m);
    });

    it('Ollama: no window check; each request keeps the num_ctx the guard sent it, and the policy lists them', async () => {
        const { chats } = fakeServer({ kind: 'ollama', id: 'gemma4:26b', loaded: null, max: 262144, replies: [{ call: ['list_sites', {}] }, { text: ANSWER }] });
        const model = gemma();
        const session = {};
        const record = await liveRun({ built: built('dev-offcentre'), arm: 'C', model, runIndex: 0, base: providerModel(model, { session }), session });
        expect(record.error).toBeNull();
        const sent = chats.map((chat) => chat.body.options.num_ctx);
        expect(sent).toHaveLength(2);
        expect(sent.every((value) => value >= guard.MIN_CONTEXT && value % guard.CONTEXT_STEP === 0)).toBe(true);
        expect(record.requests.map((request) => request.num_ctx)).toEqual(sent);
        // Read from the guard's plan, in each request as sent (the model port's prepare).
        expect(record.requests.map((request) => request.sent.plan)).toEqual(chats.map((chat) => expect.objectContaining({
            server: 'ollama', version: '0.40.2', numCtx: chat.body.options.num_ctx, numPredict: chat.body.options.num_predict, maxContext: 262144
        })));
        expect(record.policy.context).toMatchObject({
            min_tokens: null, loaded_tokens: null, model_max_tokens: 262144, num_ctx_used: [...new Set(sent)]
        });
        expect('loaded_tokens_after' in record.policy.context).toBe(false);
        expect(record.policy).toMatchObject({ server: { kind: 'ollama', version: '0.40.2' }, reply_cap: { sent: 'num_predict' }, final_round: FINAL_ROUND_FORMS.ollama });
        // num_predict keeps the reply inside the num_ctx asked for.
        expect(chats.every((chat) => chat.body.options.num_predict < chat.body.options.num_ctx)).toBe(true);
        expect(renderTranscript(record)).toContain(`- Server: Ollama 0.40.2 · context window: num_ctx ${[...new Set(sent)].map((value) => value.toLocaleString('en-US')).join(' → ')} (sized per request; the model holds 262,144)`);
    });

    it('Ollama, a model that cannot hold the turn: the guard refuses in prepare, so no request is recorded; the run\'s error is the guard\'s', async () => {
        const { chats } = fakeServer({ kind: 'ollama', id: 'gemma4:26b', loaded: null, max: 8192, replies: [{ text: ANSWER }] });
        const model = gemma();
        const record = await liveRun({ built: built('dev-offcentre'), arm: 'C', model, runIndex: 0, base: providerModel(model, { session: {} }) });
        expect(chats).toEqual([]);
        expect(record.requests).toEqual([]);
        expect(record.error).toMatchObject({ name: 'ContextLimitError', message: expect.stringMatching(/but gemma4:26b holds at most 8,192\. .* Nothing was sent/) });
        expect(record.policy.context.num_ctx_used).toEqual([]);
    });

    it('the report sets two columns\' policies side by side and flags limits that differ', async () => {
        fakeServer({ kind: 'lmstudio', loaded: 65536, replies: [{ text: ANSWER }] });
        const lmStudio = qwen();
        const a = await liveRun({ built: built('dev-offcentre'), arm: 'C', model: lmStudio, runIndex: 0, base: providerModel(lmStudio, { session: {} }), session: {} });
        guard.forgetServers();
        fakeServer({ kind: 'ollama', id: 'gemma4:26b', loaded: null, replies: [{ text: ANSWER }] });
        const ollama = gemma();
        const session = {};
        const b = await liveRun({ built: built('dev-offcentre'), arm: 'C', model: ollama, runIndex: 0, base: providerModel(ollama, { session }), session, limits: { maxRounds: 3 } });
        const report = renderReport([a, b]);
        expect(report).toMatch(/\| gemma4:26b · C \| Ollama 0\.40\.2 \| num_ctx [\d,]+, sized per request \(the model holds 262,144\) \| num_ctx minus the prompt estimate, per request \| no tools sent \| 0\.2 \(retry ≥ 0\.6\) \| 3 rounds × 4 calls/);
        expect(renderTranscript(b)).toContain('· final round: no tools sent (the native API has no tool_choice)');
        // An Ollama record from before protocol 2026-10-19.7 keeps its own form.
        const before = { ...b, policy: { ...b.policy, final_round: { tool_choice: null, tools: 'sent', how: 'not sent' } } };
        expect(renderReport([before])).toMatch(/\| no tool_choice \(tools sent\) \|/);
        expect(renderTranscript(before)).toContain('· final round: no tool_choice, tools still sent');
        expect(report).toMatch(/\| qwen3\.8-27b · C \| LM Studio \| 65,536 loaded \(minimum 49,152\) \|/);
        expect(report).toContain('**Not the same in every run: the limits.**');
        // A record from before the policy is shown as such, not dropped.
        const { policy, ...old } = a;
        void policy;
        expect(renderReport([old])).toMatch(/\| qwen3\.8-27b · C \| not recorded \| not recorded \|/);
        expect(renderTranscript(old)).toContain('- Run policy: not recorded');
    });
});

// The scripted check of plan item 1.2: what the final round (the analysis
// budget used: the loop's tool_choice 'none') sends, through the production
// port (providerModel → streamCopilotChat → client.js / ollama.js) to each
// provider preset. The model port decides the form (provider/model.js
// toolChoiceForm): the field beside the tools on the OpenAI-compatible
// servers; on Ollama, whose native API has no tool_choice, no tools (since
// protocol 2026-10-19.7: before it, the tools were still offered and a final
// reply could come back as calls only).
describe('the final round\'s tool_choice \'none\', as the model port sends it to each provider', () => {
    const echo = {
        name: 'echo',
        description: 'Echo the value.',
        parameters: { type: 'object', properties: { value: { type: 'number' } }, required: ['value'] },
        run: async ({ value }) => ({ data: { value }, checks: [] })
    };
    const KINDS = { ollama: 'ollama', lmstudio: 'lmstudio', openai: 'cloud', gemini: 'cloud', anthropic: 'cloud' };

    // What each preset was sent in the first and the final round, and the
    // form each request as sent names (its plan).
    const sentForms = async (preset) => {
        const kind = KINDS[preset.id];
        const { chats, calls } = fakeServer({ kind, id: 'm', loaded: 49152, replies: [{ call: ['echo', { value: 1 }] }, { text: 'Done.' }] });
        const port = providerModel({ baseUrl: preset.baseUrl, model: 'm', temperature: 0.2, apiKey: preset.cloud ? 'test-key' : '' }, { session: {} });
        const whys = [];
        const forms = [];
        for await (const event of runCopilotTurn({
            model: port, messages: [{ role: 'system', content: 'Rules.' }, { role: 'user', content: 'Echo 1.' }], tools: [echo], limits: { maxRounds: 1 }
        })) {
            if (event.type === 'request') {
                whys.push(event.why);
                forms.push(event.sent.plan.toolChoice);
            }
            if (event.type === 'tool_end') expect(event.error).toBeFalsy();
            if (event.type === 'done') expect(event.text).toBe('Done.');
        }
        const form = (body) => ({ tool_choice: body.tool_choice ?? null, tools: Array.isArray(body.tools) && body.tools.length ? 'sent' : 'not sent' });
        return { whys, forms, chats, calls, first: form(chats[0].body), final: form(chats[1].body) };
    };

    it.each(PROVIDER_PRESETS.map((preset) => [preset.label, preset]))('%s', async (_, preset) => {
        const { whys, forms, chats, calls, first, final } = await sentForms(preset);
        expect(whys).toEqual(['first', 'final']);
        const kind = { ollama: 'ollama', lmstudio: 'lmstudio' }[preset.id] ?? 'openai';
        // The final round's form is the policy's, for the kind the guard detects.
        expect(final).toEqual({ tool_choice: FINAL_ROUND_FORMS[kind].tool_choice, tools: FINAL_ROUND_FORMS[kind].tools });
        if (kind === 'ollama') {
            // The native API: no tool_choice in any round; the tools offered in
            // the first request, none in the final one (the model cannot call one).
            expect(chats.map((chat) => chat.url)).toEqual(['http://localhost:11434/api/chat', 'http://localhost:11434/api/chat']);
            expect(first).toEqual({ tool_choice: null, tools: 'sent' });
            expect(final).toEqual({ tool_choice: null, tools: 'not sent' });
            expect(chats.every((chat) => !('tool_choice' in chat.body))).toBe(true);
            expect(forms).toEqual(['omitted', 'tools_withheld']);
        } else {
            expect(chats.map((chat) => chat.url)).toEqual([`${preset.baseUrl}/chat/completions`, `${preset.baseUrl}/chat/completions`]);
            expect(first).toEqual({ tool_choice: 'auto', tools: 'sent' });
            // The string form, beside the tools that are still sent: the body as before.
            expect(chats[1].body.tool_choice).toBe('none');
            expect(forms).toEqual(['body', 'body']);
        }
        // A cloud provider is never probed: only its chat endpoint is called.
        if (preset.cloud) expect(calls.map((call) => call.url)).toEqual(chats.map((chat) => chat.url));
    });
});
