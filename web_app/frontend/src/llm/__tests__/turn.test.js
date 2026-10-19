// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The turn module (copilot/turn.js) at its interface: runTurn and
// describeTurn, with a scripted model port, a fake ask-user and a fake clock.
// No model runs and no server is contacted: where providerModel is under
// test, streamCopilotChat is replaced and fetch is a fake (the port's prepare
// probes the server and asks the context guard). The first requests of a dock-shaped and a
// harness-shaped spec against the model-visible snapshot files need the
// host's run module and the validation cases, so they sit beside those files
// (src/__tests__/modelVisible/turnRequests.test.js): this folder imports
// nothing from the host.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { ASK_USER, buildAskUserTool } from '../copilot/askUser';
import { READ_SKILL, buildReadSkillTool } from '../copilot/skills';
import { MEMORY_ANSWER_CHARS, TURN_METHODS, describeTurn, runTurn } from '../copilot/turn';
import { configurationBreaks, turnHistory } from '../useCopilot';
import {
    ASK_RULES, COPILOT_PROMPT, COPILOT_PROMPT_NO_TOOLS, COPILOT_PROMPT_VIEW, HISTORY_MESSAGES, copilotSystemPrompt
} from '../prompts/copilot';
import { ContextLimitError, forgetServers } from '../provider/contextGuard';
import { ANTI_REPEAT_TEMPERATURE, ToolsUnsupportedError, providerModel } from '../provider/model';

const provider = vi.hoisted(() => ({ calls: [], replies: [] }));
vi.mock('../provider/contextGuard', async (importOriginal) => ({
    ...(await importOriginal()),
    // The provider providerModel sends through: records each call (frozen at
    // the call) and streams the next scripted reply.
    streamCopilotChat: (request) => {
        provider.calls.push({ ...request, messages: request.messages.map((message) => ({ ...message })) });
        const reply = provider.replies.shift();
        return typeof reply === 'function' ? reply(request) : (async function* stream() {
            if (reply?.text) yield { content: reply.text };
            yield { finishReason: 'stop' };
        }());
    }
}));

const abortError = () => new DOMException('Aborted', 'AbortError');
const TERMINAL = ['done', 'stopped', 'failed'];

// A scripted model port: each request takes the next reply — { text,
// reasoning, calls, finish }, an Error to throw, or a function of the sent
// request returning one — and every request is recorded with its intent, the
// request as sent and its signal. `prepare` (optional) maps an intent to the
// request as sent.
const scriptedPort = (replies, { prepare = null } = {}) => {
    const requests = [];
    const port = {
        prepare: (intent, { signal } = {}) => {
            const sent = prepare ? prepare(intent) : { ...intent };
            requests.push({ intent, sent, signal });
            return sent;
        },
        stream: (sent, { signal } = {}) => {
            const request = requests.find((entry) => entry.sent === sent);
            if (request) request.streamed = true;
            const next = replies.shift();
            const reply = typeof next === 'function' ? next(sent) : next;
            return (async function* stream() {
                if (signal?.aborted) throw abortError();
                if (!reply) throw new Error('script exhausted');
                if (reply instanceof Error) throw reply;
                if (reply.reasoning) yield { reasoning: reply.reasoning };
                if (reply.text) yield { content: reply.text };
                for (const [index, call] of (reply.calls || []).entries()) {
                    yield { toolCall: { index, id: `c${requests.length}_${index}`, name: call.name, arguments: JSON.stringify(call.args ?? {}) } };
                }
                yield { finishReason: reply.finish || (reply.calls ? 'tool_calls' : 'stop') };
            }());
        }
    };
    return { port, requests };
};

// A clock that ticks 1 ms each time it is read.
const fakeClock = (start = 1000) => {
    let time = start;
    return { now: () => { time += 1; return time; } };
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

const echoTool = (overrides = {}) => ({
    name: 'echo',
    description: 'Echo the value.',
    parameters: { type: 'object', properties: { value: { type: 'number' } }, required: ['value'] },
    run: vi.fn(async ({ value }) => ({ data: { value }, checks: [{ level: 'ok', text: 'echoed' }] })),
    ...overrides
});
const viewTool = () => ({
    name: 'point_at', description: 'Point at a card.', view: true,
    parameters: { type: 'object', properties: {} }, run: async () => ({ data: { shown: true }, checks: [] })
});
const modulations = () => echoTool({ name: 'analyze_modulations', parameters: { type: 'object', properties: {} }, run: async () => ({ data: { ordered: [] }, checks: [] }) });

const skillFile = (name) => [
    '---', `name: ${name}`, `description: How to judge ${name}.`, '---', '', `# ${name}`, '', '## Steps', '', '1. Read the run.'
].join('\n');
const SKILLS = Object.freeze({ files: { '.claude/skills/fit-check/SKILL.md': skillFile('fit-check') }, references: [], load: async () => '' });

const answer = (text = 'The value is 4.') => ({ text: `${text}\n\nOutcome: achieved — read it.` });
const spec = (overrides = {}) => ({ question: 'What is the value?', run: RUN, tools: [echoTool()], clock: fakeClock(), ...overrides });
const types = (events) => events.map((event) => event.type);
const requestsOf = (events) => events.filter((event) => event.type === 'request');

// A TaSe6-like split Ta–Se shell (peaks at 2.52 and 2.74 Å) in a 10 Å cubic
// cell, as the partials file and the structure carry it: the Workbench asks
// which window a Se–Ta–Se question means.
const grid = Array.from({ length: 500 }, (_, i) => i * 0.02);
const shell = (r) => 25 * Math.exp(-((r - 2.52) ** 2) / 0.0072) + 12 * Math.exp(-((r - 2.74) ** 2) / 0.0288) + (r > 4.5 ? 1 : 0);
const SPLIT_RUN = Object.freeze({
    ...RUN,
    plotFiles: [{ name: 'run_PDFpartials.csv', plotKind: 'pdf_partials', plotData: { series: [{ label: 'Ta-Se', x: grid, y: grid.map(shell) }] } }],
    structure: { elementCounts: { Ta: 16, Se: 32 }, latticeVectors: [[10, 0, 0], [0, 10, 0], [0, 0, 10]], supercell: [1, 1, 1] }
});

describe('runTurn: the events and their order', () => {
    it('starts with start, then each request\'s round, request and reply, and ends with exactly one done', async () => {
        const tool = echoTool();
        const { port, requests } = scriptedPort(
            [{ text: 'Checking.', calls: [{ name: 'echo', args: { value: 4 } }] }, answer()],
            { prepare: (intent) => ({ ...intent, temperature: 0.3 }) }
        );
        const events = await collect(runTurn(spec({ model: port, tools: [tool] })));

        expect(types(events).filter((type) => type !== 'text')).toEqual([
            'start', 'round', 'request', 'tool_start', 'tool_end', 'round', 'request', 'done'
        ]);
        expect(events.filter((event) => TERMINAL.includes(event.type))).toHaveLength(1);
        // Every event is stamped by the turn's clock, in order.
        const stamps = events.map((event) => event.at);
        expect(stamps.every((at, index) => Number.isFinite(at) && (index === 0 || at > stamps[index - 1]))).toBe(true);
        // A request is the port's prepared request, sent as it was shown.
        const sent = requestsOf(events);
        expect(sent.map((event) => [event.round, event.why])).toEqual([[0, 'first'], [1, 'tools']]);
        expect(sent.map((event) => event.sent)).toEqual(requests.map((request) => request.sent));
        expect(sent[0].sent).toBe(requests[0].sent);
        expect(sent[0].sent).toMatchObject({ temperature: 0.3, toolChoice: 'auto', retry: null });
        expect(requests.every((request) => request.streamed)).toBe(true);
        // The tool pair shares its id; the done event is the loop's.
        const [start, end] = events.filter((event) => event.type.startsWith('tool_'));
        expect(end.id).toBe(start.id);
        expect(events.at(-1)).toMatchObject({ type: 'done', stoppedBy: 'answer', outcome: { status: 'achieved' } });
    });

    it('start is describeTurn for the same spec, stamped', async () => {
        const base = spec({ askUser: vi.fn(), skills: SKILLS });
        const { port } = scriptedPort([answer()]);
        const [start] = await collect(runTurn({ ...base, model: port }));
        const { type, at, ...rest } = start;
        expect(type).toBe('start');
        expect(Number.isFinite(at)).toBe(true);
        const description = describeTurn(base);
        expect(rest.offered.tools.map((tool) => tool.name)).toEqual(description.offered.tools.map((tool) => tool.name));
        expect({ ...rest, offered: { ...rest.offered, tools: null } }).toEqual({ ...description, offered: { ...description.offered, tools: null } });
        expect(rest).toMatchObject({
            prompt: 'tools', method: 'full', toolsActive: true, skillsError: null,
            offered: { skills: ['fit-check'], askUser: true, view: false },
            shown: { run: 'test-run', live_mode: false }
        });
    });

    it('a Workbench question comes right after start, as an ask_user pair, before any request', async () => {
        const askUser = vi.fn(async (request) => ({ answer: request.options[1].label, source: 'option', index: 1 }));
        const { port, requests } = scriptedPort([answer('Computed.')]);
        const events = await collect(runTurn(spec({
            question: 'What are the Se–Ta–Se bond angles?', run: SPLIT_RUN, tools: [], askUser, model: port
        })));
        expect(types(events).slice(0, 5)).toEqual(['start', 'tool_start', 'tool_end', 'round', 'request']);
        expect(events[1]).toMatchObject({ name: ASK_USER, app: true, interactive: true });
        expect(events[2]).toMatchObject({ id: events[1].id, name: ASK_USER, status: 'ok' });
        expect(askUser).toHaveBeenCalledTimes(1);
        expect(askUser.mock.calls[0][0].question).toBe('Which Ta–Se window should the analysis use?');
        // The answer joins the question in the request as sent.
        expect(requests[0].sent.messages.at(-1).content).toMatch(
            /^What are the Se–Ta–Se bond angles\?\n\n\(The Workbench asked: "Which Ta–Se window should the analysis use\?" — the user chose "First peak only"/
        );
        expect(events.at(-1).type).toBe('done');
    });

    it('a gate is always followed by a retract, and the next request says why it was made', async () => {
        const { port } = scriptedPort([
            { text: 'Ti forms a superstructure.\nOutcome: sure.' },
            { calls: [{ name: 'analyze_modulations', args: {} }] },
            answer('No order between cells.')
        ]);
        const events = await collect(runTurn(spec({ tools: [modulations()], model: port })));
        const gate = events.findIndex((event) => event.type === 'gate');
        expect(gate).toBeGreaterThan(0);
        expect(events[gate + 1].type).toBe('retract');
        expect(requestsOf(events).map((event) => event.why)).toEqual(['first', 'gate', 'tools']);
        expect(events.at(-1)).toMatchObject({ type: 'done', gates: [expect.objectContaining({ id: 'order-between-cells' })] });
    });

    it('a model without tools: tools_unsupported once, before any step, then a round and request under the no-tools prompt', async () => {
        // The port's typed failure (provider/model.js reads it from the provider's response).
        const unsupported = new ToolsUnsupportedError('HTTP 400: qwen does not support tools', { status: 400 });
        const { port } = scriptedPort([unsupported, answer('Plain answer.')]);
        const events = await collect(runTurn(spec({ model: port })));
        expect(types(events).filter((type) => type !== 'text')).toEqual(['start', 'round', 'request', 'tools_unsupported', 'round', 'request', 'done']);
        const [first, retry] = requestsOf(events);
        expect(first).toMatchObject({ why: 'first', sent: { toolChoice: 'auto' } });
        expect(retry.why).toBe('no_tools');
        expect(retry.sent.tools).toBeUndefined();
        expect(retry.sent.toolChoice).toBeUndefined();
        expect(retry.sent.messages[0]).toEqual({ role: 'system', content: COPILOT_PROMPT_NO_TOOLS });
    });

    it('the last request asks for the answer: why "final", tool_choice none', async () => {
        const call = { calls: [{ name: 'echo', args: { value: 1 } }] };
        const { port } = scriptedPort([call, { calls: [{ name: 'echo', args: { value: 2 } }] }, answer()]);
        const events = await collect(runTurn(spec({ model: port, config: { limits: { maxRounds: 2 } } })));
        expect(requestsOf(events).map((event) => [event.why, event.sent.toolChoice])).toEqual([
            ['first', 'auto'], ['tools', 'auto'], ['final', 'none']
        ]);
    });

    it('a reply with no answer gets a reminder: why "reminder"', async () => {
        const { port } = scriptedPort([{ reasoning: 'Thinking about it…', finish: 'stop' }, answer()]);
        const events = await collect(runTurn(spec({ model: port })));
        expect(requestsOf(events).map((event) => event.why)).toEqual(['first', 'reminder']);
        expect(requestsOf(events)[1].sent.retry).toBeNull();
    });

    it('keeps no state between calls: the same spec gives the same events', async () => {
        const run = async () => {
            const { port } = scriptedPort([{ calls: [{ name: 'echo', args: { value: 4 } }] }, answer()]);
            const events = await collect(runTurn(spec({ model: port, clock: fakeClock() })));
            return JSON.parse(JSON.stringify(events.map(({ ms, ...event }) => event)));
        };
        expect(await run()).toEqual(await run());
    });
});

describe('runTurn: how a turn ends', () => {
    it('stopped: an abort while the model streams ends with one stopped event, never a throw', async () => {
        const controller = new AbortController();
        const port = {
            stream: (sent, { signal }) => (async function* stream() {
                yield { reasoning: 'Thinking' };
                await new Promise((resolve, reject) => {
                    signal.addEventListener('abort', () => reject(abortError()), { once: true });
                    setTimeout(() => controller.abort(), 5);
                });
                yield { content: 'never' };
            }())
        };
        const events = await collect(runTurn(spec({ model: port, signal: controller.signal })));
        expect(events.at(-1)).toMatchObject({ type: 'stopped', error: { name: 'AbortError' } });
        expect(events.filter((event) => TERMINAL.includes(event.type))).toHaveLength(1);
    });

    it('stopped: Stop while the user is asked rejects the question and ends the turn', async () => {
        const controller = new AbortController();
        let seen = null;
        const askUser = vi.fn((request, { signal }) => new Promise((_, reject) => {
            seen = signal;
            signal.addEventListener('abort', () => reject(abortError()), { once: true });
        }));
        const { port, requests } = scriptedPort([
            { calls: [{ name: ASK_USER, args: { question: 'Which window?', options: ['Whole shell', 'First peak'] } }] }
        ]);
        const events = [];
        for await (const event of runTurn(spec({ model: port, askUser, signal: controller.signal }))) {
            events.push(event);
            if (event.type === 'tool_start' && event.interactive) setTimeout(() => controller.abort(), 5);
        }
        expect(seen.aborted).toBe(true);
        expect(requests).toHaveLength(1);
        expect(events.at(-1).type).toBe('stopped');
    });

    it('failed: a model error ends with its message and what to do', async () => {
        const { port } = scriptedPort([new Error('Insufficient Memory (kIOGPUCommandBufferCallbackErrorOutOfMemory)')]);
        const events = await collect(runTurn(spec({ model: port })));
        const last = events.at(-1);
        expect(last).toMatchObject({ type: 'failed', message: 'Insufficient Memory (kIOGPUCommandBufferCallbackErrorOutOfMemory)' });
        expect(last.hint).toMatch(/ran out of memory/);
        expect(last.error).toBeInstanceOf(Error);
        expect(events.filter((event) => TERMINAL.includes(event.type))).toHaveLength(1);
    });

    it('failed: a context-limit refusal, with its own message and no hint', async () => {
        const refusal = new ContextLimitError('This request needs about 30,000 tokens of context, but the model holds at most 8,192. Nothing was sent.');
        const { port } = scriptedPort([refusal]);
        const last = (await collect(runTurn(spec({ model: port })))).at(-1);
        expect(last).toMatchObject({ type: 'failed', error: refusal, message: refusal.message, hint: null });
    });

    it('failed: a request the port refuses to send (a spend cap) has no request event', async () => {
        const capped = Object.assign(new Error('spend cap reached'), { name: 'SpendCapError' });
        const port = { prepare: () => { throw capped; }, stream: () => { throw new Error('never sent'); } };
        const events = await collect(runTurn(spec({ model: port })));
        expect(types(events)).toEqual(['start', 'round', 'failed']);
        expect(events.at(-1).error).toBe(capped);
    });

    it('invalid input throws a TypeError on the first next(), never on the call', async () => {
        const { port } = scriptedPort([answer()]);
        const bad = [
            {},
            { ...spec(), model: port, question: '   ' },
            { ...spec(), model: null },
            { ...spec(), model: port, run: null },
            { ...spec(), model: port, tools: [buildAskUserTool()] },
            { ...spec(), model: port, tools: [buildReadSkillTool([{ name: 'x', description: 'x', body: 'x', references: [] }])] },
            { ...spec(), model: port, askUser: 'yes' },
            { ...spec(), model: port, skills: { list: 'none' } },
            { ...spec(), model: port, history: ['q'] },
            { ...spec(), model: port, config: { method: 'off' } },
            { ...spec(), model: port, config: { toolsEnabled: 'yes' } }
        ];
        for (const input of bad) {
            const turn = runTurn(input);
            await expect(turn.next()).rejects.toBeInstanceOf(TypeError);
        }
        expect(TURN_METHODS).toEqual(['full']);
    });
});

describe('the tools a turn offers, and the prompt for them', () => {
    const firstRequest = async (overrides) => {
        const { port, requests } = scriptedPort([answer()]);
        const events = await collect(runTurn(spec({ model: port, ...overrides })));
        return { events, start: events[0], sent: requests[0].sent };
    };
    const names = (sent) => (sent.tools || []).map((tool) => tool.function.name);

    it('ask_user and read_skill come after the caller\'s tools, in that order, only with tools on', async () => {
        const askUser = vi.fn();
        const on = await firstRequest({ tools: [echoTool()], askUser, skills: SKILLS });
        expect(names(on.sent)).toEqual(['echo', ASK_USER, READ_SKILL]);
        expect(on.start).toMatchObject({ toolsActive: true, prompt: 'tools', offered: { askUser: true, skills: ['fit-check'], view: false } });
        expect(on.sent.messages[0].content).toBe(copilotSystemPrompt({ askUser: true, skills: '- fit-check: How to judge fit-check.' }));

        const off = await firstRequest({ tools: [echoTool()], askUser, skills: SKILLS, config: { toolsEnabled: false } });
        expect(off.sent.tools).toBeUndefined();
        expect(off.start).toMatchObject({ toolsActive: false, prompt: 'no_tools', offered: { askUser: false, skills: [] } });
        expect(off.start.offered.tools.map((tool) => tool.name)).toEqual(['echo']);
        expect(off.sent.messages[0].content).toBe(COPILOT_PROMPT_NO_TOOLS);
        expect(off.sent.messages[2].content).not.toContain('analyses I run');
    });

    it('tools are active exactly when they are on and at least one is offered', async () => {
        const none = await firstRequest({ tools: [], askUser: null, skills: null });
        expect(none.start).toMatchObject({ toolsActive: false, prompt: 'no_tools' });
        expect(none.sent.tools).toBeUndefined();
        const asking = await firstRequest({ tools: [], askUser: vi.fn() });
        expect(asking.start.toolsActive).toBe(true);
        expect(names(asking.sent)).toEqual([ASK_USER]);
    });

    it('tells the model when and when not to ask, only when it can ask (the "when to ask" rules)', async () => {
        const asks = await firstRequest({ askUser: vi.fn() });
        expect(asks.sent.messages[0].content).toContain(ASK_RULES);
        expect(ASK_RULES).toMatch(/which coordination shell or r-window, which site or element\s+pair, which symmetry tolerance, which run or dataset/);
        expect(ASK_RULES).toMatch(/Never ask when the question is clear, when the\s+context or an analysis answers it/);
        expect(ASK_RULES).toMatch(/One question at a time, at most two per request/);
        expect(ASK_RULES).toMatch(/default/);
        const silent = await firstRequest({ askUser: null });
        expect(silent.sent.messages[0].content).toBe(COPILOT_PROMPT);
        expect(silent.sent.messages[0].content).not.toContain('ask_user');
        expect(names(silent.sent)).toEqual(['echo']);
    });

    it('without an askUser port: no ask_user and no Workbench question, even where it would ask', async () => {
        const { events, sent } = await firstRequest({ question: 'What are the Se–Ta–Se bond angles?', run: SPLIT_RUN, tools: [], askUser: null });
        expect(events.some((event) => event.type === 'tool_start')).toBe(false);
        expect(sent.messages.at(-1).content).toBe('What are the Se–Ta–Se bond angles?');
    });

    it('with a view tool, asks to point at cards and cite them; with tools off, no view rules either', async () => {
        const view = await firstRequest({ tools: [echoTool(), viewTool()] });
        expect(view.start).toMatchObject({ prompt: 'view', offered: { view: true } });
        expect(view.sent.messages[0].content).toBe(COPILOT_PROMPT_VIEW);
        expect(COPILOT_PROMPT_VIEW).toContain('[short label](card:<id>)');
        const off = await firstRequest({ tools: [echoTool(), viewTool()], config: { toolsEnabled: false } });
        expect(off.sent.messages[0].content).toBe(COPILOT_PROMPT_NO_TOOLS);
    });

    it('a malformed skill is reported, not thrown: the turn runs without skills', async () => {
        const broken = { files: { '.claude/skills/x/SKILL.md': 'no frontmatter' }, references: [] };
        expect(describeTurn({ run: RUN, skills: broken }).skillsError).toMatch(/has no frontmatter/);
        const { start, sent } = await firstRequest({ skills: broken });
        expect(start.skillsError).toMatch(/has no frontmatter/);
        expect(start.offered.skills).toEqual([]);
        expect(names(sent)).toEqual(['echo']);
        // Already parsed skills are taken as given.
        const parsed = describeTurn({ run: RUN, skills: { list: [{ name: 'fit-check', description: 'd', body: 'b', references: [] }] } });
        expect(parsed.offered.skills).toEqual(['fit-check']);
    });

    it('describeTurn is pure and needs no question or model', () => {
        const input = { run: RUN, tools: [echoTool()], askUser: vi.fn(), skills: SKILLS };
        const first = describeTurn(input);
        expect(first.offered.tools.map((tool) => tool.name)).toEqual(['echo', ASK_USER, READ_SKILL]);
        expect(first).toMatchObject({ prompt: 'tools', method: 'full', toolsActive: true, shown: { run: 'test-run', live_mode: false } });
        expect(JSON.stringify(describeTurn(input).shown)).toBe(JSON.stringify(first.shown));
        expect(() => describeTurn({ run: RUN, tools: [buildAskUserTool()] })).toThrow(TypeError);
    });

    it('a tool\'s run() gets the caller\'s tool context, the ask-user port and its signal', async () => {
        const focus = { show: vi.fn() };
        const askUser = vi.fn();
        const tool = echoTool();
        const { port } = scriptedPort([{ calls: [{ name: 'echo', args: { value: 1 } }] }, answer()]);
        await collect(runTurn(spec({ model: port, tools: [tool], askUser, toolContext: { focus } })));
        expect(tool.run).toHaveBeenCalledWith({ value: 1 }, { focus, askUser, signal: expect.any(AbortSignal) });
    });
});

describe('the history window', () => {
    const firstMessages = async (history) => {
        const { port, requests } = scriptedPort([answer()]);
        await collect(runTurn(spec({ model: port, history })));
        // system, context, acknowledgement … the question.
        return requests[0].sent.messages;
    };

    it('puts the context before the conversation and the question last', async () => {
        const messages = await firstMessages([]);
        expect(messages[0]).toEqual({ role: 'system', content: COPILOT_PROMPT });
        expect(messages[1].content).toContain('"run": "test-run"');
        expect(messages[2].content).toBe('Understood. I will answer questions about this RMC modeling run from the context above and the analyses I run.');
        expect(messages.at(-1)).toEqual({ role: 'user', content: 'What is the value?' });
    });

    it('keeps only the newest messages: 8, four exchanges', async () => {
        const history = Array.from({ length: 10 }, (_, index) => ({ question: `q${index}`, answer: `a${index}` }));
        const memory = (await firstMessages(history)).slice(3, -1);
        expect(HISTORY_MESSAGES).toBe(8);
        expect(memory).toHaveLength(HISTORY_MESSAGES);
        expect(memory.map((message) => message.content)).toEqual(['q6', 'a6', 'q7', 'a7', 'q8', 'a8', 'q9', 'a9']);
    });

    it('pairs each answered question with its answer and drops unanswered ones (the dock\'s turns as history)', async () => {
        const turns = [
            { role: 'user', content: 'q1' },
            { role: 'assistant', text: 'a1\nOutcome: achieved — x' },
            { role: 'user', content: 'q2 (stopped before an answer)' },
            { role: 'assistant', text: '   ' },
            { role: 'user', content: 'q3 (failed: no assistant turn)' },
            { role: 'user', content: 'q4' },
            { role: 'assistant', text: 'a4' }
        ];
        expect(turnHistory(turns)).toEqual([
            { question: 'q1', answer: 'a1\nOutcome: achieved — x' },
            { question: 'q2 (stopped before an answer)', answer: '   ' },
            { question: 'q3 (failed: no assistant turn)', answer: '' },
            { question: 'q4', answer: 'a4' }
        ]);
        expect((await firstMessages(turnHistory(turns))).slice(3, -1)).toEqual([
            { role: 'user', content: 'q1' },
            { role: 'assistant', content: 'a1\nOutcome: achieved — x' },
            { role: 'user', content: 'q4' },
            { role: 'assistant', content: 'a4' }
        ]);
    });

    it('drops the answers about an earlier configuration of the run (configKey, a Live Data save)', async () => {
        const turns = [
            { role: 'user', content: 'q1', configKey: 'epoch|0' },
            { role: 'assistant', text: 'a1' },
            { role: 'user', content: 'q2', configKey: 'epoch|0' },
            { role: 'assistant', text: 'a2' },
            { role: 'user', content: 'q3', configKey: 'epoch|1' },
            { role: 'assistant', text: 'a3' }
        ];
        // The questions asked under the current configuration, and only them.
        expect(turnHistory(turns, 'epoch|1')).toEqual([{ question: 'q3', answer: 'a3' }]);
        expect(turnHistory(turns.slice(0, 4), 'epoch|0')).toEqual([
            { question: 'q1', answer: 'a1' },
            { question: 'q2', answer: 'a2' }
        ]);
        // A new configuration before any question about it: no memory, so the
        // first request after the save is a fresh conversation's.
        expect(turnHistory(turns, 'epoch|2')).toEqual([]);
        expect(await firstMessages(turnHistory(turns, 'epoch|2'))).toEqual(await firstMessages([]));
        // A configuration seen before is still a new one: only the questions
        // since the last change count.
        expect(turnHistory([...turns, { role: 'user', content: 'q4', configKey: 'epoch|0' }, { role: 'assistant', text: 'a4' }], 'epoch|0'))
            .toEqual([{ question: 'q4', answer: 'a4' }]);
        expect(turnHistory(turns, 'epoch|0')).toEqual([]);
    });

    it('marks where the conversation crosses a new configuration (configurationBreaks)', () => {
        expect(configurationBreaks([], 'c1')).toEqual({ before: [], after: false });
        // No host key: one configuration throughout, as before configKey.
        expect(configurationBreaks([{ role: 'user', content: 'q' }, { role: 'assistant', text: 'a' }]))
            .toEqual({ before: [], after: false });
        const turns = [
            { role: 'user', content: 'q1', configKey: 'c1' },
            { role: 'assistant', text: 'a1' },
            { role: 'user', content: 'q2 (failed: no assistant turn)', configKey: 'c1' },
            { role: 'user', content: 'q3', configKey: 'c2' },
            { role: 'assistant', text: 'a3' }
        ];
        expect(configurationBreaks(turns, 'c2')).toEqual({ before: [3], after: false });
        expect(configurationBreaks(turns, 'c3')).toEqual({ before: [3], after: true });
    });

    it('compacts a long answer to its sentences with numbers and its verdict', async () => {
        const long = `${'Words without figures. '.repeat(200)}The mean is 2.52 Å.\nOutcome: achieved.`;
        const [, answerMessage] = (await firstMessages([{ question: 'q?', answer: long }])).slice(3, -1);
        expect(MEMORY_ANSWER_CHARS).toBe(2000);
        expect(answerMessage.content.length).toBeLessThanOrEqual(MEMORY_ANSWER_CHARS);
        expect(answerMessage.content).toContain('2.52 Å');
        expect(answerMessage.content).toMatch(/\nOutcome: achieved\.$/);
    });
});

describe('the production model port (providerModel)', () => {
    // A passage a looping reasoning model writes again and again.
    const LOOP = 'Wait, let me reconsider the Ta–Se window once more: the first shell spans 2.4 to 2.9 Å, '
        + 'so the coordination number should be six, but the cut at 2.65 Å splits it, so let me check again. ';
    const settings = { baseUrl: 'http://localhost:1234/v1', model: 'local-model', temperature: 0.2, apiKey: 'secret-key' };
    const json = (body, status = 200) => ({ ok: status < 400, status, json: async () => body });
    // No server answers at the base URL: the port finds none, so the request
    // goes out as an OpenAI-compatible one, unchanged (kind 'openai').
    const noServer = () => vi.stubGlobal('fetch', vi.fn(async () => { throw new TypeError('fetch failed'); }));
    // An Ollama whose model holds `max` tokens (the guard's probes only:
    // streamCopilotChat is replaced, so /api/chat is never reached).
    const ollama = (max) => {
        const paths = [];
        vi.stubGlobal('fetch', vi.fn(async (url) => {
            const path = new URL(url).pathname;
            paths.push(path);
            if (path === '/api/version') return json({ version: '0.40.2' });
            if (path === '/api/ps') return json({ models: [] });
            if (path === '/api/show') return json({ model_info: { 'gemma4.context_length': max } });
            return json({ error: 'not found' }, 404);
        }));
        return paths;
    };
    const onOllama = { baseUrl: 'http://localhost:11434/v1', model: 'gemma4:26b', temperature: 0.2, apiKey: '' };

    beforeEach(() => {
        forgetServers();
        provider.calls.length = 0;
    });
    afterEach(() => vi.unstubAllGlobals());

    it('asks again after a loop at a temperature of 0.6 or more, and never shows the API key', async () => {
        noServer();
        let pieces = 0;
        provider.replies.splice(0, provider.replies.length,
            () => (async function* stream() { for (; pieces < 1000; pieces += 1) yield { reasoning: LOOP }; }()),
            { text: 'Six Se at 2.52 Å.\n\nOutcome: achieved.' });
        const session = {};
        const events = await collect(runTurn(spec({ model: providerModel(settings, { session }) })));

        expect(pieces).toBeLessThan(40);
        const [first, again] = requestsOf(events);
        expect(first).toMatchObject({
            why: 'first', sent: { temperature: 0.2, model: 'local-model', baseUrl: 'http://localhost:1234/v1', plan: { server: 'openai', toolChoice: 'body' } }
        });
        expect(again.why).toBe('reminder');
        expect(again.sent.temperature).toBe(ANTI_REPEAT_TEMPERATURE);
        expect(again.sent.temperature).toBeGreaterThanOrEqual(0.6);
        expect(again.sent.messages.at(-1).content).toMatch(/kept repeating the same passage/);
        // The provider got the request as shown, with the key and the session added.
        expect(provider.calls.map((call) => call.temperature)).toEqual([0.2, 0.6]);
        expect(provider.calls[0]).toMatchObject({ apiKey: 'secret-key', session });
        expect(provider.calls[0].signal.aborted).toBe(true);
        expect(JSON.stringify(requestsOf(events))).not.toContain('secret-key');
        expect(events.at(-1)).toMatchObject({ type: 'done', stoppedBy: 'answer' });
    });

    it('keeps a temperature above 0.6 as it is, and maps the intent field by field, with the guard\'s plan', async () => {
        noServer();
        const port = providerModel({ ...settings, temperature: 0.9 });
        const intent = { messages: [{ role: 'user', content: 'q' }], tools: [{ type: 'function' }], toolChoice: 'none', retry: 'warmer' };
        expect(await port.prepare(intent)).toEqual({
            baseUrl: settings.baseUrl, model: 'local-model', temperature: 0.9, messages: intent.messages, tools: intent.tools, toolChoice: 'none',
            plan: { server: 'openai', toolChoice: 'body' }
        });
        expect((await providerModel(settings).prepare({ ...intent, retry: null })).temperature).toBe(0.2);
    });

    it('the request event carries the guard\'s plan; the provider gets the request\'s own fields, the plan beside them', async () => {
        const paths = ollama(262144);
        provider.replies.splice(0, provider.replies.length, answer());
        const events = await collect(runTurn(spec({ model: providerModel(onOllama, { session: {} }) })));
        const [request] = requestsOf(events);
        const { plan } = request.sent;
        expect(plan).toMatchObject({ server: 'ollama', version: '0.40.2', loadedContext: null, maxContext: 262144, toolChoice: 'omitted' });
        expect(plan.numCtx).toBeGreaterThanOrEqual(plan.needed);
        expect(plan.numPredict).toBe(plan.numCtx - plan.estimate);
        // What the provider is called with keeps its fields (the model-visible snapshots capture them).
        expect(Object.keys(provider.calls[0]).sort()).toEqual(['apiKey', 'baseUrl', 'messages', 'model', 'session', 'signal', 'temperature', 'toolChoice', 'tools']);
        expect(paths).toEqual(['/api/version', '/api/ps', '/api/show']);
        expect(events.at(-1)).toMatchObject({ type: 'done', stoppedBy: 'answer' });
    });

    it('a request the context guard refuses is never sent and has no request event (as a spend cap)', async () => {
        const paths = ollama(8192);
        const events = await collect(runTurn(spec({ model: providerModel(onOllama, { session: {} }) })));
        expect(types(events)).toEqual(['start', 'round', 'failed']);
        expect(events.at(-1).error).toBeInstanceOf(ContextLimitError);
        expect(events.at(-1).message).toMatch(/but gemma4:26b holds at most 8,192\. .* Nothing was sent/);
        expect(provider.calls).toEqual([]);
        expect(paths).not.toContain('/api/chat');
    });
});
