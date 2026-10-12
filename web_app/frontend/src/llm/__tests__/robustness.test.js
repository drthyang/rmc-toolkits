// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { runCopilotTurn } from '../copilot/agentLoop';
import {
    ModelTimeoutError, compactAnswer, explainModelError, keyNumbers, repeatsItself, textToolCalls, withDeadlines
} from '../copilot/robustness';
import { applyTurnEvent, newAssistantTurn } from '../copilot/turnState';
import { conversationMemory } from '../useCopilot';

const collect = async (generator) => {
    const events = [];
    for await (const event of generator) events.push(event);
    return events;
};

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

// A model whose replies are scripted: { text, reasoning, calls } or a stream
// function (request) => async iterable.
const scriptedModel = (replies) => {
    const requests = [];
    const callModel = (request) => {
        requests.push({ ...request, messages: request.messages.map((message) => ({ ...message })) });
        const reply = replies.shift();
        if (typeof reply === 'function') return reply(request);
        return (async function* stream() {
            if (!reply) throw new Error('script exhausted');
            if (reply.reasoning) yield { reasoning: reply.reasoning };
            if (reply.text) yield { content: reply.text };
            for (const [index, call] of (reply.calls || []).entries()) {
                yield { toolCall: { index, id: `c${requests.length}_${index}`, name: call.name, arguments: JSON.stringify(call.args ?? {}) } };
            }
            yield { finishReason: reply.finish || (reply.calls ? 'tool_calls' : 'stop') };
        }());
    };
    return { callModel, requests };
};

const echoTool = (overrides = {}) => ({
    name: 'echo',
    description: 'Echo the value.',
    parameters: { type: 'object', properties: { value: { type: 'number' } }, required: ['value'] },
    label: (args) => `Echo ${args.value}`,
    run: async ({ value }) => ({ data: { value, mean_A: 2.52 }, checks: [] }),
    ...overrides
});

const baseMessages = [{ role: 'system', content: 'SYS' }, { role: 'user', content: 'question' }];

// A passage a looping reasoning model writes again and again.
const LOOP = 'Wait, let me reconsider the Ta–Se window once more: the first shell spans 2.4 to 2.9 Å, '
    + 'so the coordination number should be six, but the cut at 2.65 Å splits it, so let me check again. ';

describe('withDeadlines', () => {
    it('passes a stream through', async () => {
        const stream = (async function* s() { yield 1; yield 2; }());
        const seen = [];
        for await (const value of withDeadlines(stream, { stallMs: 1000, maxMs: 1000 })) seen.push(value);
        expect(seen).toEqual([1, 2]);
    });

    it('stops a stream that sends nothing, and cancels it', async () => {
        let cancelled = false;
        const stream = (async function* s() { yield 1; await sleep(200); yield 2; }());
        const seen = [];
        await expect((async () => {
            for await (const value of withDeadlines(stream, { stallMs: 30, maxMs: 5000, cancel: () => { cancelled = true; } })) seen.push(value);
        })()).rejects.toMatchObject({ name: 'ModelTimeoutError', kind: 'stall' });
        expect(seen).toEqual([1]);
        expect(cancelled).toBe(true);
    });

    it('stops a reply that keeps streaming past the total limit', async () => {
        const stream = (async function* s() { for (;;) { await sleep(5); yield 'x'; } }());
        await expect((async () => {
            for await (const value of withDeadlines(stream, { stallMs: 1000, maxMs: 40 })) void value;
        })()).rejects.toMatchObject({ name: 'ModelTimeoutError', kind: 'max' });
    });
});

describe('textToolCalls', () => {
    const names = ['list_sites', 'compute_bond_angles'];

    it('reads a call written as JSON, fenced or in <tool_call> tags', () => {
        expect(textToolCalls('{"name": "list_sites", "arguments": {}}', names)).toEqual([{ name: 'list_sites', arguments: '{}' }]);
        expect(textToolCalls('```json\n{"name":"compute_bond_angles","parameters":{"center":"Ta"}}\n```', names))
            .toEqual([{ name: 'compute_bond_angles', arguments: '{"center":"Ta"}' }]);
        expect(textToolCalls('<tool_call>{"name":"list_sites","arguments":{}}</tool_call>\n<tool_call>{"name":"compute_bond_angles","arguments":"{}"}</tool_call>', names))
            .toEqual([{ name: 'list_sites', arguments: '{}' }, { name: 'compute_bond_angles', arguments: '{}' }]);
    });

    it('never reads prose, an unknown tool or a partial call as a call', () => {
        expect(textToolCalls('I will call list_sites next.', names)).toEqual([]);
        expect(textToolCalls('{"name": "delete_files", "arguments": {}}', names)).toEqual([]);
        expect(textToolCalls('Here: {"name": "list_sites", "arguments": {}}', names)).toEqual([]);
        expect(textToolCalls('{"name": "list_sites", "arguments": [1]}', names)).toEqual([]);
        expect(textToolCalls('', names)).toEqual([]);
    });
});

describe('repeatsItself', () => {
    it('catches a passage written again and again, not ordinary reasoning', () => {
        expect(repeatsItself(LOOP.repeat(5))).toBe(true);
        expect(repeatsItself(LOOP.repeat(2))).toBe(false);
        const varied = Array.from({ length: 40 }, (_, index) => `Step ${index}: the distance ${2 + index / 100} Å is inside the window. `).join('');
        expect(repeatsItself(varied)).toBe(false);
    });
});

describe('keyNumbers and compactAnswer', () => {
    it('keeps the headline of a result', () => {
        expect(keyNumbers({
            site: 'Ta1', mean_A: 2.52, cn: 6, centrosymmetric: false, window: [2.3, 2.9], note: 'a long explanation',
            rows: [1, 2, 3, 4, 5, 6], nested: { a: 1 }
        })).toEqual({ site: 'Ta1', mean_A: 2.52, cn: 6, centrosymmetric: false, window: [2.3, 2.9] });
        expect(keyNumbers([1, 2])).toBeNull();
        expect(keyNumbers({ rows: [1, 2, 3, 4, 5] })).toBeNull();
    });

    it('shortens a long answer to its sentences with numbers and its verdict', () => {
        const filler = 'This paragraph explains the method in words only. '.repeat(60);
        const answer = `${filler}The Ta–Se bond is 2.52 Å. The coordination number is 6.0. ${filler}\nOutcome: achieved — six Se around Ta.`;
        const short = compactAnswer(answer, 500);
        expect(short.length).toBeLessThanOrEqual(500);
        expect(short).toContain('The Ta–Se bond is 2.52 Å.');
        expect(short).toContain('The coordination number is 6.0.');
        expect(short).toMatch(/\nOutcome: achieved — six Se around Ta\.$/);
        expect(compactAnswer('Short. Outcome: achieved.', 500)).toBe('Short. Outcome: achieved.');
    });

    it('is applied to long answers in the conversation memory', () => {
        const long = `${'Words without figures. '.repeat(200)}The mean is 2.52 Å.\nOutcome: achieved.`;
        const memory = conversationMemory([{ role: 'user', content: 'q?' }, { role: 'assistant', text: long }]);
        expect(memory[1].content.length).toBeLessThanOrEqual(2000);
        expect(memory[1].content).toContain('2.52 Å');
    });
});

describe('explainModelError', () => {
    it('adds what to do for the known failures', () => {
        expect(explainModelError(new Error('Insufficient Memory (kIOGPUCommandBufferCallbackErrorOutOfMemory)')).hint).toMatch(/ran out of memory/);
        expect(explainModelError(new Error('Context length exceeded: 40000 > 32768')).hint).toMatch(/no longer fits/);
        expect(explainModelError(new Error('No models loaded. Please load a model')).hint).toMatch(/not loaded/);
        expect(explainModelError(new TypeError('Failed to fetch')).hint).toMatch(/connection/);
        expect(explainModelError(new ModelTimeoutError('stall', 180000))).toEqual({
            message: 'The model sent nothing for 3 min, so the request was stopped.',
            hint: expect.stringMatching(/busy with another job/)
        });
        expect(explainModelError(new Error('something else'))).toEqual({ message: 'something else', hint: null });
    });
});

describe('runCopilotTurn on a local model', () => {
    it('stops a request that stalls and explains it', async () => {
        let requestSignal = null;
        const callModel = ({ signal }) => {
            requestSignal = signal;
            return (async function* stream() { yield { reasoning: 'Thinking' }; await sleep(500); yield { content: 'late' }; }());
        };
        await expect(collect(runCopilotTurn({
            callModel, messages: baseMessages, tools: [echoTool()], limits: { stallMs: 30 }
        }))).rejects.toBeInstanceOf(ModelTimeoutError);
        expect(requestSignal.aborted).toBe(true);
    });

    it('answers after one reminder when the first reply is reasoning only', async () => {
        const { callModel, requests } = scriptedModel([
            { reasoning: 'Let me think about the bonds…', finish: 'stop' },
            { text: 'The bond is 2.52 Å.\n\nOutcome: achieved — measured.' }
        ]);
        const events = await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [echoTool()] }));
        expect(requests).toHaveLength(2);
        // The reminder joins the question: never two user messages in a row.
        expect(requests[1].messages).toHaveLength(2);
        expect(requests[1].messages[1].content).toMatch(/^question\n\nYour last reply had no answer/);
        expect(events.at(-1)).toMatchObject({ stoppedBy: 'answer', outcome: { status: 'achieved' } });
    });

    it('runs a call the model wrote as JSON text, and drops that text', async () => {
        const { callModel, requests } = scriptedModel([
            { text: '{"name": "echo", "arguments": {"value": 4}}' },
            { text: 'The value is 4.\n\nOutcome: achieved.' }
        ]);
        const events = await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [echoTool()] }));
        expect(events.filter((event) => event.type === 'tool_end')).toHaveLength(1);
        expect(events.some((event) => event.type === 'retract')).toBe(true);
        expect(JSON.parse(requests[1].messages.at(-1).content).data).toEqual({ value: 4, mean_A: 2.52 });
        expect(events.at(-1)).toMatchObject({ text: 'The value is 4.\n\nOutcome: achieved.', stoppedBy: 'answer' });
    });

    it('clears the shown text of a call written as text', () => {
        let turn = newAssistantTurn(0);
        for (const event of [{ type: 'round', round: 0 }, { type: 'text', delta: '{"name": "echo"}' }, { type: 'retract' }]) {
            turn = applyTurnEvent(turn, event, 0);
        }
        expect(turn.items.filter((item) => item.kind === 'text').map((item) => item.text)).toEqual(['']);
    });

    it('stops a looping reasoning stream early and asks again, warmer', async () => {
        let pieces = 0;
        const looping = () => (async function* stream() {
            for (; pieces < 1000; pieces += 1) yield { reasoning: LOOP };
        }());
        const { callModel, requests } = scriptedModel([looping, { text: 'Six Se at 2.52 Å.\n\nOutcome: achieved.' }]);
        const events = await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [echoTool()] }));
        expect(pieces).toBeLessThan(40);
        expect(requests[0].signal.aborted).toBe(true);
        expect(requests[1].antiRepeat).toBe(true);
        expect(requests[1].messages.at(-1).content).toMatch(/kept repeating the same passage/);
        expect(events.at(-1)).toMatchObject({ stoppedBy: 'answer' });
    });

    it('ends with a reason when the model loops again after the reminder', async () => {
        const looping = () => (async function* stream() {
            for (let index = 0; index < 1000; index += 1) yield { content: LOOP };
        }());
        const { callModel } = scriptedModel([looping, looping]);
        const events = await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [echoTool()] }));
        expect(events.at(-1)).toMatchObject({
            type: 'done', text: '', stoppedBy: 'empty', reason: 'the model kept repeating the same passage, even after a reminder'
        });
    });

    it('ends the analyses when repeated calls add nothing, long before the round budget', async () => {
        const same = { calls: [{ name: 'echo', args: { value: 1 } }] };
        const { callModel, requests } = scriptedModel([
            same, same, same,
            { text: 'The value is 1.\n\nOutcome: achieved.' }
        ]);
        const events = await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [echoTool()], limits: { maxRounds: 6 } }));
        expect(requests).toHaveLength(4);
        expect(requests[3].toolChoice).toBe('none');
        expect(requests[3].messages.at(-1).content).toMatch(/added nothing new/);
        expect(events.filter((event) => event.type === 'tool_end')).toHaveLength(3);
        expect(events.at(-1)).toMatchObject({ stoppedBy: 'no_progress', outcome: { status: 'achieved' } });
    });

    it('ends a model that repeats a window the user already settled (review 2026-10-12 #4 with the loop guard)', async () => {
        // The Workbench applies the remembered choice without asking: the
        // repeats are cached results, so they add nothing and end the analyses.
        const shellTool = echoTool({
            shellWindows: (args) => [{ centre: 'Ta', partner: 'Se', window: [args.value, 9], apply: ([lo]) => ({ ...args, value: lo }) }]
        });
        const clarifier = {
            beforeTurn: () => null,
            beforeCall: (tool, args) => (args.value === 1 ? { kind: 'cut', settled: true, args: { value: 2 }, note: 'the user already chose 2' } : null)
        };
        const same = { calls: [{ name: 'echo', args: { value: 1 } }] };
        const { callModel, requests } = scriptedModel([same, same, same, { text: 'Done.\n\nOutcome: achieved.' }]);
        const events = await collect(runCopilotTurn({
            callModel, messages: baseMessages, tools: [shellTool], clarifier, toolContext: { askUser: async () => ({}) }
        }));
        expect(requests).toHaveLength(4);
        expect(events.at(-1)).toMatchObject({ stoppedBy: 'no_progress' });
    });

    it('counts failing calls as no progress too', async () => {
        const bad = { calls: [{ name: 'echo', args: { value: 'x' } }] };
        const { callModel, requests } = scriptedModel([bad, bad, { text: 'I could not run it.\n\nOutcome: not achieved.' }]);
        const events = await collect(runCopilotTurn({ callModel, messages: baseMessages, tools: [echoTool()] }));
        expect(requests).toHaveLength(3);
        expect(events.at(-1)).toMatchObject({ stoppedBy: 'no_progress' });
    });
});
