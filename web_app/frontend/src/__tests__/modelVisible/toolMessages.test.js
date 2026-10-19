// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// What the model sees, the paths no other snapshot covers: the exact `tool`
// messages runCopilotTurn (llm/copilot/agentLoop.js) answers a round of calls
// with — a call that did not run, a question refused, asked or skipped,
// arguments the loop could not read, a tool that threw or ran out of time, a
// repeat, the Workbench's notes, a result fitted to its budget, an earlier
// result compacted, and the view tool's three answers — as raw strings, byte
// for byte and in order. One file per path in __snapshots__/tool-messages/;
// each holds every tool message of the request named, the model's own text
// left out. Also the question as the model reads it after the Workbench asked
// before the turn (its `user` message): the same note, joined to the question.
//
// The toy tools run here (no engine), except the Workbench's paths, which run
// on the Demo run with its real tools, clarifier and engines in-thread, as
// tool-results/ does. tool-results/ pins one result of every real tool
// (parsed); this file pins the paths around them.
//
// One substitution: a JSON.parse error message is the JavaScript engine's own
// wording, which differs between Node versions; it is written as
// <JSON.parse message> in the file and checked against this runtime's.
//
// A change to these files is a change to what the model sees: call it out in
// the PR.

import { beforeAll, describe, expect, it } from 'vitest';
import { runCopilotTurn } from '../../llm/copilot/agentLoop';
import { buildAskUserTool } from '../../llm/copilot/askUser';
import { buildClarifier } from '../../llm/copilot/clarify';
import { buildRunTools } from '../../llm/copilot/runTools';
import { buildFocusTools } from '../../llm/copilot/viewTools';
import { createVerification } from '../../llm/copilot/verification';
import { buildRunContext } from '../../llm/context/runContext';
import { buildCopilotMessages } from '../../llm/prompts/copilot';
import { buildCopilotTools } from '../../copilotTools';
import { inProcessEngine } from '../../workers/pcaKdeWorker';
import { loadDemoAssistantRun } from './helpers';

// --- a scripted model -------------------------------------------------------------

// Each request takes the next reply (the last repeats): { calls: [[name,
// args], …], text? } or { text }. Arguments given as a string are sent as
// that text (unparseable JSON); otherwise as JSON. Every request is kept with
// a copy of its messages as sent.
const scripted = (replies) => {
    const requests = [];
    const stream = (intent) => {
        requests.push({ ...intent, messages: intent.messages.map((message) => ({ ...message })) });
        const reply = replies[Math.min(requests.length - 1, replies.length - 1)];
        return (async function* fragments() {
            if (reply.text) yield { content: reply.text };
            for (const [index, [name, args]] of (reply.calls || []).entries()) {
                yield {
                    toolCall: {
                        index,
                        id: `call_${requests.length}_${index}`,
                        name,
                        arguments: typeof args === 'string' ? args : JSON.stringify(args)
                    }
                };
            }
            yield { finishReason: reply.calls ? 'tool_calls' : 'stop' };
        }());
    };
    return { stream, requests };
};

const ANSWER = { text: 'Recorded.\n\nOutcome: a recorded tool message.' };
const BASE = [{ role: 'system', content: 'SYS' }, { role: 'user', content: 'Run the analysis.' }];

// One turn: its events and every request as sent.
const turn = async ({ replies, tools, messages = BASE, ...rest }) => {
    const model = scripted(replies);
    const events = [];
    for await (const event of runCopilotTurn({ model: model.stream, messages, tools, ...rest })) events.push(event);
    return { events, requests: model.requests };
};

// --- the toy tools ----------------------------------------------------------------

const echo = (overrides = {}) => ({
    name: 'echo',
    description: 'Echo the value.',
    parameters: { type: 'object', properties: { value: { type: 'number', maximum: 10 } }, required: ['value'] },
    run: async ({ value }) => ({ data: { value, list: [1, 2, 3] }, checks: [{ level: 'ok', text: 'fine' }] }),
    ...overrides
});
const boom = echo({ name: 'boom', description: 'Always fails.', run: async () => { throw new Error('engine said no'); } });
const slow = echo({ name: 'slow', description: 'Never answers.', run: () => new Promise(() => {}) });
const rows = echo({
    name: 'rows',
    description: 'Many rows.',
    run: async ({ value }) => ({
        data: { value, rows: Array.from({ length: 300 }, (_, index) => ({ index, r: Number((1.5 + index * 0.01).toFixed(2)) })) },
        checks: [{ level: 'warn', text: `round ${value}: many rows` }]
    })
});
const prose = echo({
    name: 'prose',
    description: 'One long text.',
    run: async () => ({ data: { text: 'The first coordination shell. '.repeat(60) }, checks: [{ level: 'note', text: 'long' }] })
});

// --- the rendering ----------------------------------------------------------------

// Written as <JSON.parse message> in the file: the engine's own wording.
const jsonParseMessage = (text) => {
    try {
        JSON.parse(text);
    } catch (error) {
        return error.message;
    }
    throw new Error(`${text} parses`);
};

/**
 * Every tool message of `request` (and, with `user`, its question as the
 * model reads it), in order, raw: a header line naming each, then its
 * content exactly as sent.
 */
const render = (title, how, request, { user = false, substitute = null } = {}) => {
    const lines = [`# The tool messages the model reads: ${title}`, `# ${how}`];
    request.messages.forEach((message, index) => {
        if (message.role === 'tool') {
            lines.push(`## messages[${index}] · tool · ${message.tool_call_id}`, message.content);
        } else if (user && message.role === 'user' && index === request.messages.map((item) => item.role).lastIndexOf('user')) {
            lines.push(`## messages[${index}] · user (the question as the model reads it)`, message.content);
        }
    });
    const text = `${lines.join('\n')}\n`;
    return substitute ? text.split(substitute).join('<JSON.parse message>') : text;
};

const file = (name) => `__snapshots__/tool-messages/${name}.txt`;

// --- the paths with toy tools ------------------------------------------------------

describe('tool messages the model reads (toy tools)', () => {
    it('a call that runs: the result, its ref and its checks', async () => {
        const { requests } = await turn({ replies: [{ calls: [['echo', { value: 4 }]] }, ANSWER], tools: [echo()] });
        await expect(render('a call that runs', 'echo({ value: 4 }): its result with its ref', requests[1]))
            .toMatchFileSnapshot(file('ok'));
    });

    it('an unknown tool', async () => {
        const { requests } = await turn({ replies: [{ calls: [['nope', {}]] }, ANSWER], tools: [echo(), boom] });
        await expect(render('an unknown tool', 'a call to "nope", which the turn does not offer (it offers echo and boom)', requests[1]))
            .toMatchFileSnapshot(file('unknown-tool'));
    });

    it('the fifth call of a round, skipped', async () => {
        const calls = [1, 2, 3, 4, 5].map((value) => ['echo', { value }]);
        const { requests } = await turn({ replies: [{ calls }, ANSWER], tools: [echo()] });
        await expect(render('the fifth call of a round, skipped', 'five echo calls in one round: four per round run (DEFAULT_LIMITS), the fifth is answered skipped', requests[1]))
            .toMatchFileSnapshot(file('fifth-call-skipped'));
    });

    it('unparseable JSON', async () => {
        const { requests } = await turn({ replies: [{ calls: [['echo', '{value:1}']] }, ANSWER], tools: [echo()] });
        await expect(render('unparseable JSON', 'echo called with the text {value:1}, which is not JSON', requests[1], { substitute: jsonParseMessage('{value:1}') }))
            .toMatchFileSnapshot(file('unparseable-json'));
    });

    it('arguments the schema refuses, with the notes on the keys it dropped', async () => {
        const { requests } = await turn({ replies: [{ calls: [['echo', { value: 99, colour: 'red', size: 2 }]] }, ANSWER], tools: [echo()] });
        await expect(render('arguments the schema refuses', 'echo({ value: 99, colour, size }): value above its maximum, two unknown keys dropped with a note each', requests[1]))
            .toMatchFileSnapshot(file('schema-invalid'));
    });

    it('a tool that throws', async () => {
        const { requests } = await turn({ replies: [{ calls: [['boom', { value: 1 }]] }, ANSWER], tools: [boom] });
        await expect(render('a tool that throws', 'boom({ value: 1 }) throws "engine said no"', requests[1]))
            .toMatchFileSnapshot(file('thrown-error'));
    });

    it('a tool that runs out of time', async () => {
        const { requests } = await turn({ replies: [{ calls: [['slow', { value: 1 }]] }, ANSWER], tools: [slow], limits: { toolTimeoutMs: 20 } });
        await expect(render('a tool that runs out of time', 'slow({ value: 1 }) never answers; toolTimeoutMs 20', requests[1]))
            .toMatchFileSnapshot(file('timeout'));
    });

    it('an identical call repeated from the cache', async () => {
        const { requests } = await turn({
            replies: [{ calls: [['echo', { value: 2 }]] }, { calls: [['echo', { value: '2' }]] }, ANSWER], tools: [echo()]
        });
        await expect(render('an identical call repeated', 'echo({ value: 2 }), then echo({ value: "2" }): the same validated arguments, answered from the cache with a note', requests[2]))
            .toMatchFileSnapshot(file('cached-repeat'));
    });

    it('a result fitted to its budget: the largest list halved, with <key>_omitted', async () => {
        const { requests } = await turn({ replies: [{ calls: [['rows', { value: 1 }]] }, ANSWER], tools: [rows], limits: { resultCharBudget: 1200 } });
        await expect(render('a result fitted to its budget', 'rows({ value: 1 }): 300 rows, resultCharBudget 1200', requests[1]))
            .toMatchFileSnapshot(file('omitted'));
    });

    it('a result cut hard as a last resort', async () => {
        const { requests } = await turn({ replies: [{ calls: [['prose', { value: 1 }]] }, ANSWER], tools: [prose], limits: { resultCharBudget: 600 } });
        await expect(render('a result cut hard', 'prose({ value: 1 }): one long text and no list to halve, resultCharBudget 600', requests[1]))
            .toMatchFileSnapshot(file('hard-truncation'));
    });

    it('earlier results compacted to their checks and headline numbers', async () => {
        const { requests } = await turn({
            replies: [{ calls: [['rows', { value: 1 }]] }, { calls: [['rows', { value: 2 }]] }, { calls: [['rows', { value: 3 }]] }, ANSWER],
            tools: [rows],
            limits: { resultCharBudget: 900, turnCharBudget: 1000 }
        });
        await expect(render('earlier results compacted', 'rows three rounds running: resultCharBudget 900, turnCharBudget 1000; the newest is never compacted', requests[3]))
            .toMatchFileSnapshot(file('compacted-key-numbers'));
    });

    it('an error with a note, then compacted: its error and checks kept, its note dropped', async () => {
        // A window the Workbench settled (a stub clarifier) on a call that then
        // fails: an error with a note, the one error message that compacts.
        const clarifier = {
            beforeTurn: () => null,
            beforeCall: (tool, args) => (tool.name === 'boom'
                ? { kind: 'cut', settled: true, args, note: 'the user already chose 2.21–3.23 Å for the Ta–Se window 2.3–2.7 Å in this question: run in it.' }
                : null)
        };
        const { requests } = await turn({
            replies: [{ calls: [['boom', { value: 1 }]] }, { calls: [['rows', { value: 2 }]] }, { calls: [['rows', { value: 3 }]] }, ANSWER],
            tools: [boom, rows],
            clarifier,
            limits: { resultCharBudget: 900, turnCharBudget: 1000 }
        });
        await expect(render('an error with a note', 'boom with a window the Workbench settled (a stub clarifier) fails: the error and the note', requests[1]))
            .toMatchFileSnapshot(file('error-with-note'));
        await expect(render('an earlier error compacted', 'boom with a settled window (a stub clarifier) fails, then rows twice: turnCharBudget 1000', requests[3]))
            .toMatchFileSnapshot(file('compacted-error'));
    });
});

// --- ask_user ---------------------------------------------------------------------

describe('ask_user\'s messages', () => {
    const ask = (args) => ['ask_user', { question: 'Which Ta–Se window?', options: ['Whole first shell', 'First peak only'], ...args }];
    const replying = (reply) => ({ askUser: async () => reply });

    it('one question at a time: a second question in the same round is refused', async () => {
        const { requests } = await turn({
            replies: [{ calls: [ask({}), ask({ question: 'And the Ta–Ta window?' })] }, ANSWER],
            tools: [buildAskUserTool()],
            toolContext: replying({ answer: 'Whole first shell', source: 'option', index: 0 })
        });
        await expect(render('one question at a time', 'two ask_user calls in one round: the first asked (the user picks option 1), the second refused', requests[1]))
            .toMatchFileSnapshot(file('ask-one-at-a-time'));
    });

    it('already asked two: a third question in a request is refused', async () => {
        const { requests } = await turn({
            replies: [{ calls: [ask({})] }, { calls: [ask({})] }, { calls: [ask({})] }, ANSWER],
            tools: [buildAskUserTool()],
            toolContext: replying({ answer: 'First peak only', source: 'option', index: 1 })
        });
        await expect(render('already asked two', 'ask_user in three rounds: two asked (option 2 each time), the third refused', requests[3]))
            .toMatchFileSnapshot(file('ask-already-asked-two'));
    });

    it('nobody to ask: no ask-user port', async () => {
        const { requests } = await turn({ replies: [{ calls: [ask({})] }, ANSWER], tools: [buildAskUserTool()] });
        await expect(render('nobody to ask', 'ask_user with no askUser in the tools\' context', requests[1]))
            .toMatchFileSnapshot(file('ask-no-port'));
    });

    it('a malformed question: fewer than two distinct options', async () => {
        const { requests } = await turn({
            replies: [{ calls: [['ask_user', { question: 'Which window?', options: ['Whole shell', 'whole SHELL'] }]] }, ANSWER],
            tools: [buildAskUserTool()],
            toolContext: replying({ answer: 'Whole shell', source: 'option', index: 0 })
        });
        await expect(render('a malformed question', 'ask_user with two options that are one label: refused before it reaches the user', requests[1]))
            .toMatchFileSnapshot(file('ask-malformed'));
    });

    it('a skip, with a default', async () => {
        const { requests } = await turn({
            replies: [{ calls: [ask({ options: [{ label: 'Whole first shell', description: '2.21–3.23 Å, both Ta–Se peaks' }, 'First peak only'], default: 'Whole first shell' })] }, ANSWER],
            tools: [buildAskUserTool()],
            toolContext: replying({ answer: null, source: 'skipped', index: null })
        });
        await expect(render('a skip, with a default', 'ask_user with a default and its description; the user skips', requests[1]))
            .toMatchFileSnapshot(file('ask-skip-default'));
    });

    it('a skip, without a default', async () => {
        const { requests } = await turn({
            replies: [{ calls: [ask({})] }, ANSWER],
            tools: [buildAskUserTool()],
            toolContext: replying({ answer: null, source: 'skipped', index: null })
        });
        await expect(render('a skip, without a default', 'ask_user with no default; the user skips', requests[1]))
            .toMatchFileSnapshot(file('ask-skip-no-default'));
    });

    it('an answer typed in Other', async () => {
        const { requests } = await turn({
            replies: [{ calls: [ask({})] }, ANSWER],
            tools: [buildAskUserTool()],
            toolContext: replying({ answer: ' use 2.3 to 3.1 Å ', source: 'other', index: null })
        });
        await expect(render('an answer typed in Other', 'ask_user; the user types " use 2.3 to 3.1 Å "', requests[1]))
            .toMatchFileSnapshot(file('ask-typed-other'));
    });
});

// --- highlight_cards --------------------------------------------------------------

describe('highlight_cards\' messages', () => {
    const cards = [
        { id: 'plot:rmc_run_PDFpartials.csv', page: 'dashboard', title: 'Partial g(r)', about: 'the partial pair distribution functions', marks: true },
        { id: 'symmetry', page: 'dashboard', title: 'Detected SG', about: 'the space group at the card\'s tolerance' }
    ];
    const pages = { dashboard: 'Dashboard' };
    const point = (focus, args) => turn({
        replies: [{ calls: [['highlight_cards', args]] }, ANSWER],
        tools: buildFocusTools({ cards, pages }),
        toolContext: { focus }
    });

    it('shown', async () => {
        const { requests } = await point({ show: async () => ({}) }, { cards: ['plot:rmc_run_PDFpartials.csv'], marks: [{ x: 2.52, label: 'Ta–Se 2.52 Å' }] });
        await expect(render('highlight_cards, shown', 'the partial g(r) card with one mark; the focus controller shows it', requests[1]))
            .toMatchFileSnapshot(file('highlight-ok'));
    });

    it('not on screen', async () => {
        const { requests } = await point({ show: async (ids) => ({ missing: ids }) }, { cards: ['symmetry'] });
        await expect(render('highlight_cards, not on screen', 'the Detected SG card; the focus controller reports it missing', requests[1]))
            .toMatchFileSnapshot(file('highlight-not-on-screen'));
    });

    it('an unknown card', async () => {
        const { requests } = await point({ show: async () => ({}) }, { cards: ['bond-angles'] });
        await expect(render('highlight_cards, an unknown card', 'a card id the host does not have', requests[1]))
            .toMatchFileSnapshot(file('highlight-unknown-card'));
    });
});

// --- the Workbench's notes (the Demo run) -------------------------------------------

describe('the Workbench\'s notes on the Demo run', () => {
    const CUT = { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.3, r12_max: 2.7 };
    let setup;

    beforeAll(async () => {
        const { run, rmc6fText } = await loadDemoAssistantRun();
        const tools = [
            ...buildRunTools(run, { cards: [] }),
            ...buildCopilotTools({ request: inProcessEngine(rmc6fText), structure: run.structure, elements: Object.keys(run.structure.elementCounts || {}).sort() }),
            buildAskUserTool()
        ];
        const context = buildRunContext(run);
        setup = { tools, context, clarifier: buildClarifier(run, { noiseFloor: context.symmetry?.noise_floor_A ?? null }) };
    });

    // A turn on the Demo, as toolResults.test.js runs one: the dock's tools
    // without the view tools, the Workbench's clarifier, the verification.
    const demoTurn = async ({ question = 'Run the analysis.', replies, reply }) => {
        const { tools, context, clarifier } = setup;
        const asked = [];
        const messages = buildCopilotMessages({ context, history: [], userText: question, toolsEnabled: true, viewTools: false, askUser: true, skills: '' });
        const result = await turn({
            replies,
            tools,
            messages,
            toolContext: { askUser: async (request) => { asked.push(request.question); return reply; } },
            clarifier,
            verification: createVerification({ context, catalog: clarifier.catalog, offered: tools.map((tool) => tool.name) })
        });
        return { ...result, asked };
    };

    it('asked before the call: the user keeps the window as asked', async () => {
        const { requests, asked } = await demoTurn({ replies: [{ calls: [['compute_bond_angles', CUT]] }, ANSWER], reply: { answer: 'Keep 2.3–2.7 Å', source: 'option', index: 1 } });
        expect(asked).toHaveLength(1);
        await expect(render('the Workbench asked: Keep', 'compute_bond_angles in 2.3–2.7 Å, which cuts the Ta–Se shell: the Workbench asks, the user keeps the window', requests[1]))
            .toMatchFileSnapshot(file('workbench-asked-keep'));
    });

    it('asked before the call: the user types a window', async () => {
        const { requests, asked } = await demoTurn({ replies: [{ calls: [['compute_bond_angles', CUT]] }, ANSWER], reply: { answer: '2.2 to 3.0', source: 'other', index: null } });
        expect(asked).toHaveLength(1);
        await expect(render('the Workbench asked: a typed window', 'compute_bond_angles in 2.3–2.7 Å: the Workbench asks, the user types "2.2 to 3.0"', requests[1]))
            .toMatchFileSnapshot(file('workbench-asked-typed-window'));
    });

    it('asked before the call: analyze_clusters in a window that cuts the Ta–Ta shell, the default kept', async () => {
        const { requests, asked } = await demoTurn({
            replies: [{ calls: [['analyze_clusters', { element: 'Ta', r_min: 2.65, r_max: 2.95 }]] }, ANSWER],
            reply: { answer: null, source: 'skipped', index: null }
        });
        expect(asked).toHaveLength(1);
        await expect(render('the Workbench asked: clusters', 'analyze_clusters for Ta in 2.65–2.95 Å, which cuts the Ta–Ta shell: the Workbench asks, the user skips', requests[1]))
            .toMatchFileSnapshot(file('workbench-asked-clusters'));
    });

    it('settled without asking: the question names the first peak', async () => {
        const { requests, asked } = await demoTurn({
            question: 'What are the Se–Ta–Se angles in the first Ta–Se peak only?',
            replies: [{ calls: [['compute_bond_angles', CUT]] }, ANSWER],
            reply: { answer: null, source: 'skipped', index: null }
        });
        expect(asked).toEqual([]);
        await expect(render('the Workbench settled: the first peak named', 'compute_bond_angles in 2.3–2.7 Å for a question that names the first Ta–Se peak: no question, the window kept as the user\'s', requests[1]))
            .toMatchFileSnapshot(file('workbench-settled'));
    });

    it('settled without asking, and answered from the cache: the window chosen earlier in the question', async () => {
        const { requests, asked } = await demoTurn({
            replies: [{ calls: [['compute_bond_angles', CUT]] }, { calls: [['compute_bond_angles', CUT]] }, ANSWER],
            reply: { answer: 'Whole first shell', source: 'option', index: 0 }
        });
        expect(asked).toHaveLength(1);
        await expect(render('the Workbench settled, then the cache', 'compute_bond_angles in 2.3–2.7 Å twice: asked once (the whole shell chosen); the repeat runs in the chosen window unasked, answered from the cache', requests[2]))
            .toMatchFileSnapshot(file('workbench-settled-cached'));
    });

    it('asked before the turn: the answer joins the question the model reads', async () => {
        const { requests, asked } = await demoTurn({
            question: 'Compute the bond angles around Ta.',
            replies: [ANSWER],
            reply: { answer: null, source: 'skipped', index: null }
        });
        expect(asked).toHaveLength(1);
        await expect(render('the Workbench asked before the turn', 'a bond-angle question naming only Ta: the Workbench asks which neighbours, the user skips; no tool message, the question as read', requests[0], { user: true }))
            .toMatchFileSnapshot(file('workbench-before-turn'));
    });
});
