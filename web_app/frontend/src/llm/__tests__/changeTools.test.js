// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Change tools at the tool contract's and the turn's interface: what a tool
// does to the Workbench is its effect (copilot/effects.js, one list for every
// tool); a change tool asks the user through the turn's approve port before
// it changes anything, runs with no time limit while the card waits, is never
// answered from the cache and is not an analysis, and returns a change
// record (applied, declined, unchanged, stale) the step keeps and the model
// never reads. A scripted model port; no model runs. The host's change tool
// on the Bond Geometry page: src/__tests__/copilotChangeTools.test.jsx.

import { describe, expect, it, vi } from 'vitest';
import {
    COPILOT_TOOL_EFFECTS, DOCK_EFFECTS, TOOL_EFFECTS, effectOf, isChangeTool, toolsWithEffects
} from '../copilot/effects';
import { CHANGE_STATUSES, executeToolCall, toolMessage, toolResult } from '../copilot/tools';
import { describeTurn, runTurn } from '../copilot/turn';
import { applyTurnEvent, changeReply, newAssistantTurn } from '../copilot/turnState';

const abortError = () => new DOMException('Aborted', 'AbortError');

// A scripted model port: each request takes the next reply ({ text, calls }).
const scriptedPort = (replies) => {
    const requests = [];
    return {
        requests,
        port: {
            prepare: (intent) => {
                const sent = { ...intent };
                requests.push(sent);
                return sent;
            },
            stream: (sent, { signal } = {}) => {
                const reply = replies.shift();
                return (async function* stream() {
                    if (signal?.aborted) throw abortError();
                    if (!reply) throw new Error('script exhausted');
                    if (reply.text) yield { content: reply.text };
                    for (const [index, call] of (reply.calls || []).entries()) {
                        yield { toolCall: { index, id: `c${requests.length}_${index}`, name: call.name, arguments: JSON.stringify(call.args ?? {}) } };
                    }
                    yield { finishReason: reply.calls ? 'tool_calls' : 'stop' };
                }());
            }
        }
    };
};

const collect = async (generator) => {
    const events = [];
    for await (const event of generator) events.push(event);
    return events;
};

const RUN = Object.freeze({
    runName: 'test-run', plotFiles: [], rValueFile: null, structure: null, symmetry: null,
    runSettings: null, pcaSites: null, liveData: false
});

const echo = () => ({
    name: 'echo',
    description: 'Echo the value.',
    parameters: { type: 'object', properties: { value: { type: 'number' } }, required: ['value'] },
    run: vi.fn(async ({ value }) => ({ data: { value }, checks: [{ level: 'ok', text: 'echoed' }] }))
});

// A change tool as a test declares one (its name is not on the list): it
// asks the approve port, then reports what it did.
const setter = ({ run } = {}) => ({
    name: 'set_value',
    effect: 'change',
    description: 'Set the value on the page.',
    parameters: { type: 'object', properties: { value: { type: 'number' } }, required: ['value'], additionalProperties: false },
    run: vi.fn(run || (async ({ value }, { approve, signal }) => {
        const line = `value ${value}`;
        const { approved, auto } = await approve({ tool: 'set_value', title: 'Set the value', line, details: [`0 → ${value}`] }, { signal });
        if (!approved) {
            return { data: { declined: true }, checks: [{ level: 'note', text: 'The user declined this change.' }], change: { status: 'declined', line, details: [], auto, id: null } };
        }
        return { data: { applied: true, value }, checks: [{ level: 'ok', text: `Applied: ${line}.` }], change: { status: 'applied', line, details: [], auto, id: 'v1' } };
    }))
});

const answer = { text: 'Done.\n\nOutcome: achieved — set it.' };
const spec = (overrides = {}) => ({ question: 'Set the value to 3.', run: RUN, tools: [echo(), setter()], ...overrides });

describe('the one list of tool effects (copilot/effects.js)', () => {
    it('names every effect a tool can have, and every tool of the Copilot with one of them', () => {
        expect(TOOL_EFFECTS).toEqual(['compute', 'view', 'change']);
        expect(DOCK_EFFECTS).toEqual(TOOL_EFFECTS);
        expect(Object.values(COPILOT_TOOL_EFFECTS).every((effect) => TOOL_EFFECTS.includes(effect))).toBe(true);
        expect(Object.entries(COPILOT_TOOL_EFFECTS).filter(([, effect]) => effect !== 'compute')).toEqual([
            ['highlight_cards', 'view'], ['orient_3d_view', 'view'], ['current_view', 'view'], ['set_bond_geometry', 'change']
        ]);
        expect(Object.isFrozen(COPILOT_TOOL_EFFECTS)).toBe(true);
    });

    it('a listed name takes the list\'s effect; another tool its own declaration: effect, then view, else compute', () => {
        expect(effectOf({ name: 'compute_bond_angles', effect: 'change' })).toBe('compute');
        expect(effectOf({ name: 'set_bond_geometry' })).toBe('change');
        expect(effectOf({ name: 'highlight_cards' })).toBe('view');
        expect(effectOf({ name: 'probe', effect: 'change' })).toBe('change');
        expect(effectOf({ name: 'probe', view: true })).toBe('view');
        expect(effectOf({ name: 'probe', effect: 'nonsense' })).toBe('compute');
        expect(effectOf({ name: 'toString' })).toBe('compute');
        expect(effectOf(null)).toBe('compute');
        expect(isChangeTool({ name: 'set_bond_geometry' })).toBe(true);
        expect(isChangeTool({ name: 'compute_bond_angles' })).toBe(false);
    });

    it('toolsWithEffects keeps the tools of the effects asked for, in order, and refuses an unknown effect', () => {
        const tools = [echo(), { name: 'highlight_cards', view: true }, null, setter(), { name: 'list_sites' }];
        expect(toolsWithEffects(tools, ['compute']).map((tool) => tool.name)).toEqual(['echo', 'list_sites']);
        expect(toolsWithEffects(tools, DOCK_EFFECTS).map((tool) => tool.name)).toEqual(['echo', 'highlight_cards', 'set_value', 'list_sites']);
        expect(toolsWithEffects(tools, ['change']).map((tool) => tool.name)).toEqual(['set_value']);
        expect(() => toolsWithEffects(tools, ['write'])).toThrow(TypeError);
        expect(() => toolsWithEffects(tools, [])).toThrow(TypeError);
        expect(toolsWithEffects(null, ['compute'])).toEqual([]);
    });
});

describe('a change record on the tool result (tools.js)', () => {
    const record = (status, extra = {}) => ({ status, line: 'value 3', details: ['0 → 3'], auto: false, id: null, ...extra });
    const fail = [{ level: 'fail', text: 'Not applied: the page changed.' }];

    it('rides on an ok result (applied, declined, unchanged) and on a refused one (stale), frozen', () => {
        expect(CHANGE_STATUSES).toEqual(['applied', 'declined', 'unchanged', 'stale']);
        for (const status of ['applied', 'declined', 'unchanged']) {
            const result = toolResult({ change: record(status) });
            expect(result.change).toEqual(record(status));
            expect(Object.isFrozen(result.change)).toBe(true);
            expect(Object.isFrozen(result.change.details)).toBe(true);
        }
        expect(toolResult({ resultStatus: 'refused', checks: fail, change: record('stale') }).change.status).toBe('stale');
    });

    it('the constructor refuses a record that breaks an invariant', () => {
        expect(() => toolResult({ resultStatus: 'error', error: 'bad', change: record('applied') })).toThrow(/no change record/);
        expect(() => toolResult({ change: record('stale') })).toThrow(/stale change is a 'refused' result/);
        expect(() => toolResult({ resultStatus: 'refused', checks: fail, change: record('applied') })).toThrow(/'applied' change is an 'ok' result/);
        expect(() => toolResult({ change: record('done') })).toThrow(/change.status must be one of/);
        expect(() => toolResult({ change: record('applied', { line: '' }) })).toThrow(/change.line/);
        expect(() => toolResult({ change: record('applied', { details: [3] }) })).toThrow(/change.details/);
        expect(() => toolResult({ change: record('applied', { auto: 'yes' }) })).toThrow(/change.auto/);
        expect(() => toolResult({ change: record('applied', { id: '' }) })).toThrow(/change.id/);
    });

    it('executeToolCall keeps a change tool\'s record, drops another tool\'s and a malformed one, and never sends it', async () => {
        const approve = vi.fn(async () => ({ approved: true, auto: false }));
        const applied = await executeToolCall(setter(), '{"value":3}', { timeoutMs: 1000, context: { approve } });
        expect(applied.resultStatus).toBe('ok');
        expect(applied.change).toEqual({ status: 'applied', line: 'value 3', details: [], auto: false, id: 'v1' });
        const message = toolMessage(applied, { budget: 3000, ref: 'r1' }).content;
        expect(message).not.toContain('"change"');
        expect(message).not.toContain('v1');

        const compute = { ...echo(), run: async () => ({ data: { value: 1 }, checks: [], change: record('applied') }) };
        expect((await executeToolCall(compute, '{"value":1}', { timeoutMs: 1000 })).change).toBeUndefined();
        const malformed = setter({ run: async () => ({ data: { applied: true }, checks: [], change: { status: 'applied' } }) });
        expect((await executeToolCall(malformed, '{"value":1}', { timeoutMs: 1000 })).change).toBeUndefined();
        // A stale record without its refusal does not fit how the call ended: dropped.
        const unrefused = setter({ run: async () => ({ data: { stale: true }, checks: [], change: record('stale') }) });
        expect((await executeToolCall(unrefused, '{"value":1}', { timeoutMs: 1000 })).change).toBeUndefined();
    });

    it('a strict schema refuses an unknown argument before the tool runs', async () => {
        const tool = setter();
        const result = await executeToolCall(tool, '{"value":3,"window":2}', { timeoutMs: 1000, context: { approve: vi.fn() } });
        expect(result.resultStatus).toBe('error');
        expect(result.error).toBe('unknown argument "window"; the arguments are value');
        expect(tool.run).not.toHaveBeenCalled();
    });
});

describe('the turn and a change tool', () => {
    it('is offered only with an approve port, which reaches its run() beside the signal', async () => {
        const without = describeTurn({ run: RUN, tools: [echo(), setter()] });
        expect(without.offered.tools.map((tool) => tool.name)).toEqual(['echo']);
        const approve = vi.fn(async () => ({ approved: true, auto: false }));
        const withPort = describeTurn({ run: RUN, tools: [echo(), setter()], approve });
        expect(withPort.offered.tools.map((tool) => tool.name)).toEqual(['echo', 'set_value']);
        expect(() => describeTurn({ run: RUN, approve: 'yes' })).toThrow(/approve must be a function or null/);

        const tool = setter();
        const { port } = scriptedPort([{ calls: [{ name: 'set_value', args: { value: 3 } }] }, answer]);
        await collect(runTurn(spec({ model: port, tools: [tool], approve })));
        expect(tool.run).toHaveBeenCalledWith({ value: 3 }, { approve, signal: expect.any(AbortSignal) });
        expect(approve).toHaveBeenCalledWith(
            { tool: 'set_value', title: 'Set the value', line: 'value 3', details: ['0 → 3'] }, { signal: expect.any(AbortSignal) }
        );
    });

    it('waits for the approval with no time limit, and the step ends with its effect and change record, not as an analysis', async () => {
        // The card takes longer than the per-call limit: the change still lands.
        const approve = () => new Promise((resolve) => { setTimeout(() => resolve({ approved: true, auto: false }), 60); });
        const { port, requests } = scriptedPort([{ calls: [{ name: 'set_value', args: { value: 3 } }] }, answer]);
        const events = await collect(runTurn(spec({ model: port, approve, config: { limits: { toolTimeoutMs: 20 } } })));
        const start = events.find((event) => event.type === 'tool_start');
        const end = events.find((event) => event.type === 'tool_end');
        expect(start).toMatchObject({ name: 'set_value', effect: 'change' });
        expect(end).toMatchObject({ name: 'set_value', resultStatus: 'ok', status: 'ok', analysis: false, error: null });
        expect(end.change).toEqual({ status: 'applied', line: 'value 3', details: [], auto: false, id: 'v1' });
        // The model read the data and the checks, never the record.
        const message = requests[1].messages.at(-1);
        expect(JSON.parse(message.content)).toEqual({ ref: 'r1', checks: [{ level: 'ok', text: 'Applied: value 3.' }], data: { applied: true, value: 3 } });
        // A turn whose only step is a change answered from the run context.
        const done = events.at(-1);
        expect(done.type).toBe('done');
        expect(done.outcome.note).toMatch(/answered from the run context/);
        // The dock's step: its effect and record kept, its reading.
        let turn = newAssistantTurn(0);
        for (const event of events) turn = applyTurnEvent(turn, event, 0);
        const step = turn.items.find((item) => item.kind === 'tool');
        expect(step).toMatchObject({ effect: 'change', change: { status: 'applied', id: 'v1' }, analysis: false });
        expect(changeReply(step)).toEqual({ state: 'applied', text: 'value 3' });
    });

    it('is never answered from the cache: the same call twice runs twice', async () => {
        const approve = vi.fn(async () => ({ approved: true, auto: true }));
        const tool = setter();
        const compute = echo();
        const call = { calls: [{ name: 'set_value', args: { value: 3 } }, { name: 'echo', args: { value: 1 } }] };
        const { port } = scriptedPort([call, { ...call }, answer]);
        const events = await collect(runTurn(spec({ model: port, tools: [compute, tool], approve })));
        expect(tool.run).toHaveBeenCalledTimes(2);
        expect(compute.run).toHaveBeenCalledTimes(1);
        const ends = events.filter((event) => event.type === 'tool_end');
        expect(ends.map((event) => [event.name, event.cached])).toEqual([['set_value', false], ['echo', false], ['set_value', false], ['echo', true]]);
        expect(ends.filter((event) => event.name === 'set_value').every((event) => event.change.auto)).toBe(true);
    });

    it('a declined change returns {"declined": true} and its record; Stop while the card waits ends the turn as stopped', async () => {
        const decline = vi.fn(async () => ({ approved: false, auto: false }));
        const { port, requests } = scriptedPort([{ calls: [{ name: 'set_value', args: { value: 3 } }] }, answer]);
        const events = await collect(runTurn(spec({ model: port, approve: decline })));
        const end = events.find((event) => event.type === 'tool_end');
        expect(end).toMatchObject({ resultStatus: 'ok', data: { declined: true }, change: { status: 'declined' } });
        expect(JSON.parse(requests[1].messages.at(-1).content).data).toEqual({ declined: true });

        const controller = new AbortController();
        const waiting = (request, { signal }) => new Promise((resolve, reject) => {
            signal.addEventListener('abort', () => reject(abortError()), { once: true });
            setTimeout(() => controller.abort(), 10);
        });
        const stopped = await collect(runTurn(spec({
            model: scriptedPort([{ calls: [{ name: 'set_value', args: { value: 3 } }] }, answer]).port, approve: waiting, signal: controller.signal
        })));
        expect(stopped.at(-1)).toMatchObject({ type: 'stopped', stoppedBy: 'stopped' });
        expect(stopped.some((event) => event.type === 'tool_end')).toBe(false);
    });

    it('changeReply reads every way a change step can stand', () => {
        const step = (fields) => ({ kind: 'tool', status: 'ok', error: null, checks: [], change: null, ...fields });
        expect(changeReply(step({ status: 'running' }))).toEqual({ state: 'running', text: '' });
        expect(changeReply(step({ status: 'stopped' })).state).toBe('stopped');
        expect(changeReply(step({ status: 'failed', error: 'no run is open' }))).toEqual({ state: 'not_applied', text: 'no run is open' });
        expect(changeReply(step({ change: { status: 'declined', line: 'value 3' } }))).toEqual({ state: 'declined', text: 'value 3' });
        expect(changeReply(step({ change: { status: 'unchanged', line: 'value 3' } })).state).toBe('unchanged');
        expect(changeReply(step({
            status: 'failed', checks: [{ level: 'fail', text: 'Not applied: the inputs changed.' }], change: { status: 'stale', line: 'value 3' }
        }))).toEqual({ state: 'stale', text: 'Not applied: the inputs changed.' });
    });
});
