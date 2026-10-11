// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it, vi } from 'vitest';
import { DEFAULT_LIMITS, assistantContent, isToolsUnsupportedError, runCopilotTurn } from '../copilot/agentLoop';
import { ASK_USER, buildAskUserTool } from '../copilot/askUser';
import { budgetJson, fitToBudget, stepStatus, usableTools } from '../copilot/tools';
import { ASK_RULES, buildCopilotMessages } from '../prompts/copilot';

// A scripted model: each request consumes the next reply. A reply is either
// { text } or { calls: [{ name, args }] } (optionally with text), or an
// Error to throw, or a function of the request returning one of those (a
// model that reads what it was sent). Every request's arguments are recorded.
const scriptedModel = (replies) => {
    const requests = [];
    const callModel = (request) => {
        requests.push({ ...request, messages: request.messages.map((message) => ({ ...message })) });
        const next = replies.shift();
        const reply = typeof next === 'function' ? next(request) : next;
        return (async function* stream() {
            if (!reply) throw new Error('script exhausted');
            if (reply instanceof Error) throw reply;
            if (reply.text) {
                for (const piece of reply.text.match(/.{1,7}/gs)) yield { content: piece };
            }
            for (const [index, call] of (reply.calls || []).entries()) {
                const args = typeof call.args === 'string' ? call.args : JSON.stringify(call.args ?? {});
                yield {
                    toolCall: {
                        index,
                        id: call.id ?? `c${requests.length}_${index}`,
                        name: call.name,
                        arguments: args.slice(0, 5),
                        extraContent: call.extraContent ?? null
                    }
                };
                yield { toolCall: { index, id: null, name: null, arguments: args.slice(5) } };
            }
            yield { finishReason: reply.calls ? 'tool_calls' : 'stop' };
        }());
    };
    return { callModel, requests };
};

const collect = async (generator) => {
    const events = [];
    for await (const event of generator) events.push(event);
    return events;
};

const echoTool = (overrides = {}) => ({
    name: 'echo',
    description: 'Echo the value.',
    parameters: { type: 'object', properties: { value: { type: 'number', maximum: 10 } }, required: ['value'] },
    label: ({ value }) => `Echo ${value}`,
    page: 'dashboard',
    run: vi.fn(async ({ value }) => ({ data: { value, list: [1, 2, 3] }, checks: [{ level: 'ok', text: 'fine' }] })),
    ...overrides
});

const baseMessages = [{ role: 'system', content: 'SYS' }, { role: 'user', content: 'question' }];

describe('runCopilotTurn', () => {
    it('runs a tool, feeds the result back, and parses the Outcome line', async () => {
        const tool = echoTool();
        const { callModel, requests } = scriptedModel([
            { text: 'Checking.', calls: [{ name: 'echo', args: { value: 4 } }] },
            { text: 'The value is 4.\nOutcome: achieved — echoed 4.' }
        ]);
        const events = await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [tool] }));

        expect(tool.run).toHaveBeenCalledWith({ value: 4 }, { signal: expect.any(AbortSignal) });
        const start = events.find((event) => event.type === 'tool_start');
        expect(start).toMatchObject({ name: 'echo', label: 'Echo 4', args: { value: 4 }, page: 'dashboard' });
        const end = events.find((event) => event.type === 'tool_end');
        expect(end).toMatchObject({ status: 'ok', data: { value: 4, list: [1, 2, 3] }, error: null, cached: false });

        // The first request offers the tools; the second carries the exchange.
        expect(requests[0].tools[0]).toMatchObject({ type: 'function', function: { name: 'echo' } });
        expect(requests[0].toolChoice).toBe('auto');
        const second = requests[1].messages;
        expect(second.at(-2)).toEqual({
            role: 'assistant',
            content: 'Checking.',
            tool_calls: [{ id: 'c1_0', type: 'function', function: { name: 'echo', arguments: '{"value":4}' } }]
        });
        expect(second.at(-1).role).toBe('tool');
        expect(second.at(-1).tool_call_id).toBe('c1_0');
        expect(JSON.parse(second.at(-1).content)).toEqual({
            checks: [{ level: 'ok', text: 'fine' }],
            data: { value: 4, list: [1, 2, 3] }
        });

        const done = events.at(-1);
        expect(done).toMatchObject({ type: 'done', rounds: 1, stoppedBy: 'answer' });
        expect(done.outcome).toEqual({ status: 'achieved', note: 'echoed 4.', body: 'The value is 4.' });
    });

    it('returns argument errors to the model, which can correct them', async () => {
        const tool = echoTool();
        const { callModel, requests } = scriptedModel([
            { calls: [{ name: 'echo', args: { value: 99 } }] },
            { calls: [{ name: 'echo', args: { value: 9 } }] },
            { text: 'Done.' }
        ]);
        const events = await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [tool] }));
        const ends = events.filter((event) => event.type === 'tool_end');
        expect(ends.map((event) => event.status)).toEqual(['failed', 'ok']);
        expect(ends[0].error).toBe('value must be ≤ 10, got 99');
        expect(JSON.parse(requests[1].messages.at(-1).content)).toEqual({ error: 'value must be ≤ 10, got 99' });
        expect(tool.run).toHaveBeenCalledTimes(1);
        // The assistant message is null-content, never '' (Anthropic rejects empty text).
        expect(requests[1].messages.at(-2).content).toBeNull();
    });

    it('answers unknown tools, invalid JSON and thrown errors as tool errors', async () => {
        const failing = echoTool({ name: 'boom', run: async () => { throw new Error('engine said no'); } });
        const { callModel, requests } = scriptedModel([
            { calls: [
                { name: 'nope', args: {} },
                { name: 'echo', args: '{value:1}' },
                { name: 'boom', args: { value: 1 } }
            ] },
            { text: 'ok' }
        ]);
        const events = await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [echoTool(), failing] }));
        const errors = events.filter((event) => event.type === 'tool_end').map((event) => event.error);
        expect(errors[0]).toBe('unknown tool "nope"; available: echo, boom');
        expect(errors[1]).toMatch(/not valid JSON/);
        expect(errors[2]).toBe('engine said no');
        // Every listed call is answered, and unparseable arguments are echoed as '{}'.
        const assistant = requests[1].messages.find((message) => message.tool_calls);
        expect(assistant.tool_calls.map((call) => call.function.arguments)).toEqual(['{}', '{}', '{"value":1}']);
        const answered = requests[1].messages.filter((message) => message.role === 'tool').map((message) => message.tool_call_id);
        expect(answered).toEqual(assistant.tool_calls.map((call) => call.id));
    });

    it('reuses the result of an identical call and tells the model', async () => {
        const tool = echoTool();
        const { callModel, requests } = scriptedModel([
            { calls: [{ name: 'echo', args: { value: 2 } }] },
            { calls: [{ name: 'echo', args: { value: '2' } }] },
            { text: 'ok' }
        ]);
        const events = await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [tool] }));
        expect(tool.run).toHaveBeenCalledTimes(1);
        expect(events.filter((event) => event.type === 'tool_end').map((event) => event.cached)).toEqual([false, true]);
        expect(JSON.parse(requests[2].messages.at(-1).content).note).toMatch(/identical call/);
    });

    it('answers calls beyond the per-step limit as skipped', async () => {
        const tool = echoTool();
        const { callModel, requests } = scriptedModel([
            { calls: [1, 2, 3].map((value) => ({ name: 'echo', args: { value } })) },
            { text: 'ok' }
        ]);
        const events = await collect(runCopilotTurn({
            callModel, messages: baseMessages, tools: [tool], limits: { maxCallsPerRound: 2 }
        }));
        expect(tool.run).toHaveBeenCalledTimes(2);
        expect(events.filter((event) => event.type === 'tool_end').at(-1).error).toMatch(/^skipped/);
        expect(requests[1].messages.filter((message) => message.role === 'tool')).toHaveLength(3);
    });

    it('asks for the answer with tool_choice none once the round budget is used', async () => {
        const tool = echoTool();
        const { callModel, requests } = scriptedModel([
            { calls: [{ name: 'echo', args: { value: 1 } }] },
            { calls: [{ name: 'echo', args: { value: 2 } }] },
            // The last request ignores tool calls even if the server sends some.
            { text: 'Best effort.\nOutcome: partial — budget used.', calls: [{ name: 'echo', args: { value: 3 } }] }
        ]);
        const events = await collect(runCopilotTurn({
            callModel, messages: baseMessages, tools: [tool], limits: { maxRounds: 2 }
        }));
        expect(requests).toHaveLength(3);
        expect(requests[2].toolChoice).toBe('none');
        expect(requests[2].messages.at(-1).content).toMatch(/used the analysis budget/);
        expect(tool.run).toHaveBeenCalledTimes(2);
        expect(events.at(-1)).toMatchObject({ type: 'done', stoppedBy: 'max_rounds', rounds: 2 });
        expect(events.at(-1).outcome.status).toBe('partial');
    });

    it('retries once without tools when the model cannot call them', async () => {
        const unsupported = Object.assign(new Error('HTTP 400: registry.ollama.ai/library/gemma3 does not support tools'), { status: 400 });
        const { callModel, requests } = scriptedModel([unsupported, { text: 'Plain answer.' }]);
        const events = await collect(runCopilotTurn({
            callModel, messages: baseMessages, tools: [echoTool()], noToolsSystemPrompt: 'PLAIN'
        }));
        expect(events.some((event) => event.type === 'tools_unsupported')).toBe(true);
        expect(requests[1].tools).toBeUndefined();
        expect(requests[1].toolChoice).toBeUndefined();
        expect(requests[1].messages[0].content).toBe('PLAIN');
        expect(events.at(-1)).toMatchObject({ type: 'done', text: 'Plain answer.' });
    });

    it('does not mistake other HTTP failures for missing tool support', async () => {
        const quota = Object.assign(new Error('HTTP 429: quota'), { status: 429 });
        const { callModel } = scriptedModel([quota]);
        await expect(collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [echoTool()] })))
            .rejects.toThrow('HTTP 429: quota');
        const err = (message, status = 400) => Object.assign(new Error(message), { status });
        expect(isToolsUnsupportedError(err('HTTP 500: model does not support tools', 500))).toBe(false);
        expect(isToolsUnsupportedError(err('HTTP 400: tools are not supported by this model'))).toBe(true);
        expect(isToolsUnsupportedError(err('HTTP 400: function calling is not supported'))).toBe(true);
        // A schema problem is a real error, never "this model has no tools".
        expect(isToolsUnsupportedError(err("HTTP 400: Invalid schema for function 'list_sites'"))).toBe(false);
        expect(isToolsUnsupportedError(err('HTTP 400: tools.0.function.parameters is invalid'))).toBe(false);
    });

    it('surfaces a tool-schema 400 instead of falling back', async () => {
        const schema = Object.assign(new Error("HTTP 400: Invalid schema for function 'echo'"), { status: 400 });
        const { callModel, requests } = scriptedModel([schema]);
        await expect(collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [echoTool()] })))
            .rejects.toThrow('Invalid schema');
        expect(requests).toHaveLength(1);
    });

    it('sends no tools when they are switched off', async () => {
        const { callModel, requests } = scriptedModel([{ text: 'hi' }]);
        await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [echoTool()], toolsEnabled: false }));
        expect(requests[0].tools).toBeUndefined();
    });

    it('stops waiting for a tool when the turn is aborted', async () => {
        const controller = new AbortController();
        const slow = echoTool({ run: () => new Promise(() => {}) });
        const { callModel } = scriptedModel([{ calls: [{ name: 'echo', args: { value: 1 } }] }]);
        const generator = runCopilotTurn({ callModel, messages: baseMessages, tools: [slow], signal: controller.signal });
        const seen = [];
        const run = (async () => {
            for await (const event of generator) {
                seen.push(event.type);
                if (event.type === 'tool_start') controller.abort();
            }
        })();
        await expect(run).rejects.toMatchObject({ name: 'AbortError' });
        expect(seen).toContain('tool_start');
        expect(seen).not.toContain('tool_end');
    });

    it('times out a tool that never resolves', async () => {
        const slow = echoTool({ run: () => new Promise(() => {}) });
        const { callModel } = scriptedModel([{ calls: [{ name: 'echo', args: { value: 1 } }] }, { text: 'ok' }]);
        const events = await collect(runCopilotTurn({
            callModel, messages: baseMessages, tools: [slow], limits: { toolTimeoutMs: 20 }
        }));
        expect(events.find((event) => event.type === 'tool_end').error).toMatch(/did not finish within/);
    });

    it('fits large results to the budget without dropping the checks', async () => {
        const big = echoTool({
            run: async () => ({
                data: { rows: Array.from({ length: 500 }, (_, index) => ({ index, value: index * 1.5 })) },
                checks: [{ level: 'warn', text: 'big' }]
            })
        });
        const { callModel, requests } = scriptedModel([{ calls: [{ name: 'echo', args: { value: 1 } }] }, { text: 'ok' }]);
        await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [big], limits: { resultCharBudget: 1200 } }));
        const content = requests[1].messages.at(-1).content;
        expect(content.length).toBeLessThanOrEqual(1300);
        const parsed = JSON.parse(content);
        expect(parsed.checks).toEqual([{ level: 'warn', text: 'big' }]);
        expect(parsed.data.rows_omitted).toBeGreaterThan(400);
        expect(parsed.data.rows.length + parsed.data.rows_omitted).toBe(500);
    });

    it('echoes a call\'s extra_content (Gemini thought_signature) and sends no index', async () => {
        const signature = { google: { thought_signature: 'c2ln' } };
        const { callModel, requests } = scriptedModel([
            { calls: [{ name: 'echo', args: { value: 1 }, extraContent: signature }] },
            { text: 'ok' }
        ]);
        await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [echoTool()] }));
        const [call] = requests[1].messages.at(-2).tool_calls;
        expect(call.extra_content).toEqual(signature);
        expect('index' in call).toBe(false);
    });

    it('sends whitespace-only text before a tool call as null content', async () => {
        const { callModel, requests } = scriptedModel([
            { text: '\n\n ', calls: [{ name: 'echo', args: { value: 1 } }] },
            { text: 'ok' }
        ]);
        await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [echoTool()] }));
        expect(requests[1].messages.at(-2).content).toBeNull();
        expect(assistantContent(' x ')).toBe(' x ');
        expect(assistantContent('')).toBeNull();
    });

    it('says so when the last round comes back with tool calls only (Ollama ignores tool_choice)', async () => {
        const { callModel } = scriptedModel([
            { calls: [{ name: 'echo', args: { value: 1 } }] },
            { calls: [{ name: 'echo', args: { value: 2 } }] }
        ]);
        const events = await collect(runCopilotTurn({
            callModel, messages: baseMessages, tools: [echoTool()], limits: { maxRounds: 1 }
        }));
        const done = events.at(-1);
        expect(done.stoppedBy).toBe('max_rounds');
        expect(done.text).toMatch(/analysis budget/);
        expect(done.outcome.status).toBe('not achieved');
        expect(events.filter((event) => event.type === 'text').map((event) => event.delta).join('')).toBe(done.text);
    });

    it('marks an empty reply', async () => {
        const { callModel } = scriptedModel([{ text: '  ' }]);
        const events = await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [echoTool()] }));
        expect(events.at(-1)).toMatchObject({ type: 'done', text: '', outcome: null, stoppedBy: 'empty' });
    });

    it('compacts earlier results to their checks when the turn outgrows its budget', async () => {
        const big = echoTool({
            run: async ({ value }) => ({
                data: { value, rows: Array.from({ length: 200 }, (_, index) => index) },
                checks: [{ level: 'ok', text: `round ${value}` }]
            })
        });
        const { callModel, requests } = scriptedModel([
            { calls: [{ name: 'echo', args: { value: 1 } }] },
            { calls: [{ name: 'echo', args: { value: 2 } }] },
            { calls: [{ name: 'echo', args: { value: 3 } }] },
            { text: 'ok' }
        ]);
        await collect(runCopilotTurn({
            callModel, messages: baseMessages, tools: [big], limits: { resultCharBudget: 900, turnCharBudget: 1500 }
        }));
        const toolContents = requests[3].messages.filter((message) => message.role === 'tool').map((message) => JSON.parse(message.content));
        expect(toolContents[0]).toEqual({
            compacted: 'this earlier result was shortened to fit the context; its checks are kept',
            checks: [{ level: 'ok', text: 'round 1' }]
        });
        // The newest result is never compacted.
        expect(toolContents[2].data.value).toBe(3);
        const total = requests[3].messages.filter((message) => message.role === 'tool').reduce((sum, message) => sum + message.content.length, 0);
        expect(total).toBeLessThanOrEqual(1500);
    });

    it('aborts the tool\'s own signal when it times out', async () => {
        let seen = null;
        const slow = echoTool({ run: (args, { signal }) => { seen = signal; return new Promise(() => {}); } });
        const { callModel } = scriptedModel([{ calls: [{ name: 'echo', args: { value: 1 } }] }, { text: 'ok' }]);
        await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [slow], limits: { toolTimeoutMs: 10 } }));
        expect(seen.aborted).toBe(true);
    });

    it('starts no new round of tools after the turn time limit', async () => {
        const { callModel, requests } = scriptedModel([
            { calls: [{ name: 'echo', args: { value: 1 } }] },
            { text: 'Late.\nOutcome: partial — out of time.' }
        ]);
        const events = await collect(runCopilotTurn({
            callModel, messages: baseMessages, tools: [echoTool()], limits: { maxTurnMs: -1 }
        }));
        expect(requests[0].toolChoice).toBe('auto');
        expect(requests[1].toolChoice).toBe('none');
        expect(events.at(-1).stoppedBy).toBe('max_rounds');
    });

    it('gives each step what its page should load, from the validated arguments', async () => {
        const tool = echoTool({ pageParams: ({ value }) => ({ shown: value }) });
        const failing = echoTool({ name: 'boom', pageParams: () => ({ never: true }), run: async () => { throw new Error('no'); } });
        const throwing = echoTool({ name: 'odd', pageParams: () => { throw new Error('bad mapping'); } });
        const { callModel } = scriptedModel([
            { calls: [{ name: 'echo', args: { value: '3' } }, { name: 'boom', args: { value: 1 } }, { name: 'odd', args: { value: 1 } }] },
            { text: 'ok' }
        ]);
        const events = await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [tool, failing, throwing] }));
        expect(events.filter((event) => event.type === 'tool_end').map((event) => event.pageParams)).toEqual([{ shown: 3 }, null, null]);
    });

    it('uses sensible default limits', () => {
        expect(DEFAULT_LIMITS).toEqual({
            maxRounds: 6, maxCallsPerRound: 4, toolTimeoutMs: 120000, maxTurnMs: 600000, resultCharBudget: 3000, turnCharBudget: 12000,
            maxQuestions: 2
        });
    });
});

describe('tool helpers', () => {
    it('keeps only well-formed, uniquely named tools', () => {
        const ok = echoTool();
        expect(usableTools([ok, { ...ok }, { ...ok, name: 'Bad Name' }, { name: 'x', description: 'd' }, null])).toEqual([ok]);
    });

    it('derives a step status from the error and checks', () => {
        expect(stepStatus({ error: 'x', checks: [] })).toBe('failed');
        expect(stepStatus({ error: null, checks: [{ level: 'fail', text: 'a' }] })).toBe('failed');
        expect(stepStatus({ error: null, checks: [{ level: 'warn', text: 'a' }, { level: 'ok', text: 'b' }] })).toBe('warn');
        expect(stepStatus({ error: null, checks: [{ level: 'note', text: 'a' }] })).toBe('ok');
    });

    it('trims nested arrays and hard-cuts as a last resort, never silently', () => {
        const value = { a: { b: Array.from({ length: 100 }, (_, i) => i) }, c: 'x' };
        const fitted = fitToBudget(value, 120);
        expect(JSON.stringify(fitted).length).toBeLessThanOrEqual(120);
        expect(fitted.a.b_omitted + fitted.a.b.length).toBe(100);
        const cut = JSON.parse(budgetJson({ s: 'y'.repeat(500) }, 100));
        expect(cut.truncated).toBe(true);
        expect(fitToBudget(null, 10)).toBeNull();
    });
});

describe('view tools', () => {
    it('are never answered from the cache, get the tool context, and name their cards', async () => {
        const focus = { show: vi.fn() };
        const view = {
            name: 'point',
            view: true,
            description: 'Point at a card.',
            parameters: { type: 'object', properties: { card: { type: 'string' } }, required: ['card'] },
            cards: ({ card }) => [card],
            run: vi.fn(async ({ card }, context) => {
                context.focus.show([card]);
                return { data: { card }, checks: [{ level: 'ok', text: 'shown' }] };
            })
        };
        const { callModel } = scriptedModel([
            { calls: [{ name: 'point', args: { card: 'a' } }] },
            { calls: [{ name: 'point', args: { card: 'a' } }] },
            { text: 'Done.\nOutcome: achieved — pointed.' }
        ]);
        const events = await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [view], toolContext: { focus } }));
        expect(view.run).toHaveBeenCalledTimes(2);
        expect(view.run).toHaveBeenLastCalledWith({ card: 'a' }, { focus, signal: expect.any(AbortSignal) });
        expect(focus.show).toHaveBeenCalledTimes(2);
        expect(events.filter((event) => event.type === 'tool_end').map((event) => event.cached)).toEqual([false, false]);
        expect(events.find((event) => event.type === 'tool_start').cards).toEqual(['a']);
    });

    it('lists no cards for a tool without them, or whose cards() throws', async () => {
        const broken = echoTool({ name: 'broken', cards: () => { throw new Error('x'); } });
        const { callModel } = scriptedModel([
            { calls: [{ name: 'echo', args: { value: 1 } }, { name: 'broken', args: { value: 2 } }] },
            { text: 'ok\nOutcome: achieved — x' }
        ]);
        const events = await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [echoTool(), broken] }));
        expect(events.filter((event) => event.type === 'tool_start').map((event) => event.cards)).toEqual([[], []]);
    });
});


// Replays of the clarifying-question flow (ask_user). The numbers are the
// Demo run's (GaTa4Se8, 250 K): the first Ta–Se peak window 2.33–2.71 Å holds
// 3.86 Se per Ta, the whole first shell 2.33–2.96 Å (peaks at 2.52 and 2.74 Å)
// 5.89 — the octahedron the spot check's first window cut in two.
describe('ask_user: the agent asks instead of guessing a parameter', () => {
    const SHELL = [2.33, 2.96];
    const FIRST_PEAK = [2.33, 2.71];
    const grPeaks = () => ({
        name: 'find_gr_peaks',
        description: 'Partial g(r) peaks and the first shell.',
        parameters: { type: 'object', properties: { pair: { type: 'string' } }, required: ['pair'] },
        run: vi.fn(async ({ pair }) => ({
            data: {
                pair,
                peaks_A: [{ r: 2.52, fwhm: 0.19 }, { r: 2.74, fwhm: 0.2 }],
                first_peak_window_A: FIRST_PEAK,
                first_shell: { window_A: SHELL, peaks_in_shell_A: [2.52, 2.74], cn: { Se_per_Ta: 5.89 } }
            },
            checks: [{ level: 'note', text: 'The first shell holds 2 peaks: a split shell.' }]
        }))
    });
    const bondAngles = () => ({
        name: 'compute_bond_angles',
        description: 'Bond angles A–B–C.',
        parameters: {
            type: 'object',
            properties: {
                end1: { type: 'string' }, apex: { type: 'string' }, end2: { type: 'string' },
                r12_min: { type: 'number' }, r12_max: { type: 'number' }
            },
            required: ['end1', 'apex', 'end2', 'r12_min', 'r12_max']
        },
        run: vi.fn(async (args) => ({ data: { window: [args.r12_min, args.r12_max] }, checks: [{ level: 'ok', text: 'angles' }] }))
    });
    // The last tool result the model was sent, parsed.
    const lastResult = (request) => JSON.parse(request.messages.filter((message) => message.role === 'tool').at(-1).content);

    // The model of an ambiguous question: it finds the shells, offers them as
    // options, and computes with whatever window the answer names.
    const ambiguousScript = () => [
        { text: 'First, the Ta–Se shells.', calls: [{ name: 'find_gr_peaks', args: { pair: 'Ta-Se' } }] },
        (request) => {
            const { data } = lastResult(request);
            return {
                text: 'The shell is split; which window should the angles use?',
                calls: [{
                    name: ASK_USER,
                    args: {
                        question: 'Which Ta–Se window should the Se–Ta–Se angles use?',
                        options: [
                            { label: 'Whole first shell', description: `${data.first_shell.window_A.join('–')} Å, both Ta–Se peaks (CN 5.89)` },
                            { label: 'First peak only', description: `${data.first_peak_window_A.join('–')} Å, the short bonds (CN ~3.9)` }
                        ],
                        default: 'Whole first shell'
                    }
                }]
            };
        },
        (request) => {
            const { data } = lastResult(request);
            const windows = { 'Whole first shell': SHELL, 'First peak only': FIRST_PEAK };
            // A typed answer: the model reads the window off it.
            const typed = /(\d+(?:\.\d+)?)\D+(\d+(?:\.\d+)?)/.exec(data.answer || '');
            const window = windows[data.answer] || (typed ? [Number(typed[1]), Number(typed[2])] : SHELL);
            return { calls: [{ name: 'compute_bond_angles', args: { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: window[0], r12_max: window[1] } }] };
        },
        { text: 'The Se–Ta–Se angles …\nOutcome: achieved — computed in the chosen window.' }
    ];

    const run = async (replies, askUser, extra = {}) => {
        const tools = [grPeaks(), bondAngles(), buildAskUserTool()];
        const { callModel, requests } = scriptedModel(replies);
        const events = await collect(runCopilotTurn({
            callModel, messages: baseMessages, tools, toolContext: askUser ? { askUser } : null, ...extra
        }));
        return { events, requests, tools };
    };

    it.each([
        ['the default option', { answer: 'Whole first shell', source: 'option', index: 0 }, SHELL],
        ['the other option', { answer: 'First peak only', source: 'option', index: 1 }, FIRST_PEAK],
        ['an answer typed in Other', { answer: 'use 2.3 to 3.1 Å', source: 'other', index: null }, [2.3, 3.1]],
        ['Skip (the default stands)', { answer: null, source: 'skipped', index: null }, SHELL]
    ])('an ambiguous question asks which shell, and %s decides the window', async (_, reply, window) => {
        const askUser = vi.fn(async () => reply);
        const { events, requests, tools } = await run(ambiguousScript(), askUser);

        // The question the card shows: two options with their meaning, the default marked.
        expect(askUser).toHaveBeenCalledTimes(1);
        const [request, { signal }] = askUser.mock.calls[0];
        expect(request).toEqual({
            question: 'Which Ta–Se window should the Se–Ta–Se angles use?',
            options: [
                { label: 'Whole first shell', description: '2.33–2.96 Å, both Ta–Se peaks (CN 5.89)' },
                { label: 'First peak only', description: '2.33–2.71 Å, the short bonds (CN ~3.9)' }
            ],
            allowOther: true,
            defaultIndex: 0
        });
        expect(signal).toBeInstanceOf(AbortSignal);
        const start = events.find((event) => event.type === 'tool_start' && event.name === ASK_USER);
        expect(start).toMatchObject({ interactive: true, label: 'Question for you' });

        // The answer reaches the model as the tool result …
        const sent = JSON.parse(requests[2].messages.at(-1).content);
        expect(sent.checks[0].text).toMatch(reply.source === 'skipped' ? /skipped .* default, "Whole first shell"/ : /Use exactly this choice|Follow it exactly/);
        // … and the analysis runs in the window the user chose.
        const angles = tools[1].run;
        expect(angles).toHaveBeenCalledTimes(1);
        expect(angles.mock.calls[0][0]).toMatchObject({ r12_min: window[0], r12_max: window[1] });
        expect(events.at(-1)).toMatchObject({ type: 'done', stoppedBy: 'answer' });
    });

    it('a clear question runs the analysis without asking', async () => {
        const askUser = vi.fn();
        const { events, tools } = await run([
            { calls: [{ name: 'compute_bond_angles', args: { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.33, r12_max: 2.96 } }] },
            { text: 'The angles …\nOutcome: achieved — computed.' }
        ], askUser);
        expect(askUser).not.toHaveBeenCalled();
        expect(events.some((event) => event.interactive)).toBe(false);
        expect(tools[1].run).toHaveBeenCalledTimes(1);
        expect(events.at(-1).stoppedBy).toBe('answer');
    });

    it('tells the model when and when not to ask, in the system prompt', () => {
        const [system] = buildCopilotMessages({ context: {}, userText: 'q', askUser: true });
        expect(system.content).toContain(ASK_RULES);
        expect(ASK_RULES).toMatch(/which coordination shell or r-window, which site or element\s+pair, which symmetry tolerance, which run or dataset/);
        expect(ASK_RULES).toMatch(/Never ask when the question is clear, when the\s+context or an analysis answers it/);
        expect(ASK_RULES).toMatch(/One question at a time, at most two per request/);
        expect(ASK_RULES).toMatch(/default/);
        expect(buildCopilotMessages({ context: {}, userText: 'q' })[0].content).not.toContain('ask_user');
    });

    it('stops cleanly when the turn is cancelled while the question waits', async () => {
        const controller = new AbortController();
        let seenSignal = null;
        const askUser = vi.fn((request, { signal }) => new Promise((_, reject) => {
            seenSignal = signal;
            signal.addEventListener('abort', () => reject(new DOMException('Aborted', 'AbortError')), { once: true });
        }));
        const tools = [grPeaks(), bondAngles(), buildAskUserTool()];
        const { callModel, requests } = scriptedModel(ambiguousScript());
        const seen = [];
        const turn = (async () => {
            for await (const event of runCopilotTurn({ callModel, messages: baseMessages, tools, signal: controller.signal, toolContext: { askUser } })) {
                seen.push(event);
                if (event.type === 'tool_start' && event.interactive) setTimeout(() => controller.abort(), 5);
            }
        })();
        await expect(turn).rejects.toMatchObject({ name: 'AbortError' });
        expect(seenSignal.aborted).toBe(true);
        // No answer was sent, no analysis ran after the question, no further request.
        expect(seen.filter((event) => event.type === 'tool_end').map((event) => event.name)).toEqual(['find_gr_peaks']);
        expect(tools[1].run).not.toHaveBeenCalled();
        expect(requests).toHaveLength(2);
    });

    it('waits without a time limit, and the wait does not count toward the turn limit', async () => {
        const askUser = vi.fn(() => new Promise((resolve) => setTimeout(() => resolve({ answer: 'First peak only', source: 'option', index: 1 }), 60)));
        const { requests, tools } = await run(ambiguousScript(), askUser, { limits: { toolTimeoutMs: 10, maxTurnMs: 40 } });
        // After a 60 ms answer under a 40 ms turn limit, the model still gets its tools.
        expect(requests[2].toolChoice).toBe('auto');
        expect(tools[1].run.mock.calls[0][0]).toMatchObject({ r12_min: 2.33, r12_max: 2.71 });
    });

    it('asks one question per step and at most two per request, and never from the cache', async () => {
        const askUser = vi.fn(async () => ({ answer: 'A', source: 'option', index: 0 }));
        const ask = { name: ASK_USER, args: { question: 'Which?', options: ['A', 'B'] } };
        const { events, requests } = await run([
            { calls: [ask, { ...ask, args: { ...ask.args, question: 'And?' } }] },
            { calls: [ask] },
            { calls: [ask] },
            { text: 'Done.\nOutcome: partial — assumed the default.' }
        ], askUser);
        // Asked in steps 1 and 2 (the same question twice: a real exchange
        // each time); the second call of step 1 and the third question refused.
        expect(askUser).toHaveBeenCalledTimes(2);
        const ends = events.filter((event) => event.type === 'tool_end');
        expect(ends.map((event) => event.cached)).toEqual([false, false, false, false]);
        expect(ends[1].error).toMatch(/one question at a time/);
        expect(ends[3].error).toMatch(/already asked 2 questions .* proceed with your default/);
        expect(JSON.parse(requests[3].messages.at(-1).content).error).toMatch(/already asked 2 questions/);
    });

    it('reads plain strings as options, and says so when nobody can be asked', async () => {
        const { events } = await run([
            { calls: [{ name: ASK_USER, args: { question: 'Which run?', options: ['250 K', '5 K'] } }] },
            { text: 'Assumed 250 K.\nOutcome: partial — assumed the run.' }
        ], null);
        const end = events.find((event) => event.type === 'tool_end');
        expect(end.status).toBe('failed');
        expect(end.error).toMatch(/cannot be asked here: proceed with your default/);
        expect(end.args.options).toEqual([{ label: '250 K' }, { label: '5 K' }]);
    });
});
