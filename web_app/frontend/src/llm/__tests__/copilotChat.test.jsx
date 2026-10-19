// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The Copilot dock end to end against a mocked OpenAI-compatible server:
// nothing is fetched until the rail opens it; the model asks for a tool,
// the dock runs it and shows a step card, the result reaches the next
// request, and the answer ends in an Outcome chip. "Show in <page>" hands the
// host the step's page and parameters. Collapsing keeps the conversation, an
// answer that lands while collapsed marks the rail unread, and the open state
// and width are remembered. The agent points at host cards (highlight_cards):
// the card breathes, its page opens, a click elsewhere clears it, and the
// answer's card chip brings it back. No live model is involved.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import React, { act } from 'react';
import { createRoot } from 'react-dom/client';
import CopilotChat from '../components/CopilotChat';
import { DOCK_STORAGE_KEY } from '../components/dockLayout';
import { forgetServers } from '../provider/contextGuard';
import { saveSettings } from '../settings';

const sse = (chunks) => {
    const encoder = new TextEncoder();
    const body = new ReadableStream({
        start(controller) {
            chunks.forEach((chunk) => controller.enqueue(encoder.encode(`data: ${JSON.stringify(chunk)}\n\n`)));
            controller.enqueue(encoder.encode('data: [DONE]\n\n'));
            controller.close();
        }
    });
    return { ok: true, status: 200, body };
};
const delta = (fields, finish = null) => ({ choices: [{ delta: fields, finish_reason: finish }] });

// The mocked server is OpenAI-compatible only: the Copilot's probes of an
// Ollama or LM Studio API (provider/contextGuard.js) get a 404.
const notFound = { ok: false, status: 404, json: async () => ({ error: 'not found' }) };
const isServerProbe = (url) => /\/api\//.test(String(url));
const modelList = { ok: true, status: 200, json: async () => ({ data: [{ id: 'tool-model' }] }) };

const replies = () => [
    sse([
        delta({ content: 'Let me compute the angles.' }),
        delta({ tool_calls: [{ index: 0, id: 'call_a', type: 'function', function: { name: 'compute_bond_angles', arguments: '' } }] }),
        delta({ tool_calls: [{ index: 0, function: { arguments: '{"end1":"Ta","apex":"Ta","end2":"Ta","r12_min":2.7,"r12_max":3.3}' } }] }, 'tool_calls')
    ]),
    sse([
        delta({ content: 'The Ta–Ta–Ta angles peak at **60°**.\n\nOutcome: achieved — the Ta4 cluster angle.' }, 'stop')
    ])
];

const settle = () => act(async () => {
    for (let i = 0; i < 10; i += 1) await new Promise((resolve) => setTimeout(resolve, 0));
});

// Settle until `ready()` holds (a multi-step answer takes several hops).
const settleUntil = async (ready, rounds = 20) => {
    for (let i = 0; i < rounds && !ready(); i += 1) await settle();
};

const setTextarea = (textarea, value) => {
    Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype, 'value').set.call(textarea, value);
    textarea.dispatchEvent(new Event('input', { bubbles: true }));
};

describe('CopilotChat', () => {
    let container;
    let root;
    let fetchMock;
    let chatBodies;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        forgetServers();
        localStorage.removeItem(DOCK_STORAGE_KEY);
        saveSettings({ baseUrl: 'http://localhost:11434/v1', model: '', apiKey: '', copilotTools: true });
        chatBodies = [];
        const queue = replies();
        fetchMock = vi.fn(async (url, options = {}) => {
            if (isServerProbe(url)) return notFound;
            if (String(url).endsWith('/v1/models')) return modelList;
            chatBodies.push(JSON.parse(options.body));
            return queue.shift();
        });
        vi.stubGlobal('fetch', fetchMock);
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.unstubAllGlobals();
    });

    const bondAngles = {
        name: 'compute_bond_angles',
        description: 'Bond angles.',
        parameters: {
            type: 'object',
            properties: {
                end1: { type: 'string' }, apex: { type: 'string' }, end2: { type: 'string' },
                r12_min: { type: 'number' }, r12_max: { type: 'number' }
            },
            required: ['end1', 'apex', 'end2', 'r12_min', 'r12_max']
        },
        label: ({ end1, apex, end2 }) => `Bond angles ${end1}–${apex}–${end2}`,
        page: 'geometry',
        pageParams: (args) => ({ end1: args.end1, apex: args.apex, end2: args.end2, r12Min: args.r12_min, r12Max: args.r12_max }),
        run: vi.fn(async () => ({ data: { modes: [{ angle_deg: 59.5 }] }, checks: [{ level: 'ok', text: '41336 angles.' }] }))
    };

    const chat = (props = {}) => (
        <CopilotChat
            runName="demo"
            plotFiles={[{ plotKind: 'neutron_sq', plotData: { title: 'S(Q)', metrics: { rwp: 0.04 } } }]}
            tools={[bondAngles]}
            pages={{ geometry: 'Bond Geometry' }}
            runKey="run-1"
            {...props}
        />
    );

    const openChat = async () => {
        act(() => container.querySelector('.llm-chat-launcher').click());
        await settle();
    };

    const ask = async (text) => {
        act(() => setTextarea(container.querySelector('.llm-composer-input'), text));
        act(() => container.querySelector('.llm-composer').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true })));
        await settle();
    };

    it('is a rail until opened, and contacts no server before that', async () => {
        const onOpenChange = vi.fn();
        act(() => root.render(chat({ onOpenChange })));
        await settle();
        const dock = container.querySelector('.llm-dock');
        expect(dock.tagName).toBe('ASIDE');
        expect(dock.getAttribute('aria-label')).toBe('AI Copilot');
        // jsdom has no matchMedia: docked beside the pages.
        expect(dock.dataset.mode).toBe('docked');
        expect(dock.dataset.open).toBe('false');
        const rail = container.querySelector('.llm-dock-rail');
        expect(rail.textContent).toBe('AI Copilot');
        expect(rail.getAttribute('aria-expanded')).toBe('false');
        expect(container.querySelector('.llm-chat-panel')).toBeNull();
        expect(fetchMock).not.toHaveBeenCalled();

        await openChat();
        expect(onOpenChange).toHaveBeenLastCalledWith(true);
        expect(dock.dataset.open).toBe('true');
        const panel = container.querySelector('.llm-chat-panel');
        expect(container.querySelector('.llm-dock-rail')).toBeNull();
        expect(panel.textContent).toContain('Connected');
        // A resize handle on its left edge.
        const handle = container.querySelector('.llm-dock-resize');
        expect(handle.getAttribute('role')).toBe('separator');
        expect(handle.getAttribute('aria-orientation')).toBe('vertical');
        // The composer takes focus when the dock opens; the open state is kept.
        expect(document.activeElement).toBe(container.querySelector('.llm-composer-input'));
        expect(JSON.parse(localStorage.getItem(DOCK_STORAGE_KEY)).open).toBe(true);
    });

    it('runs a tool, shows the step, ends with the Outcome chip, and shows the analysis on its page', async () => {
        const onOpenPage = vi.fn();
        act(() => root.render(chat({ onOpenPage })));
        await openChat();
        await ask('What are the Ta–Ta–Ta angles?');

        expect(chatBodies).toHaveLength(2);
        expect(chatBodies[0].tools.map((tool) => tool.function.name)).toEqual(['read_run_data', 'compute_bond_angles', 'ask_user']);
        const toolMessage = chatBodies[1].messages.at(-1);
        expect(toolMessage).toMatchObject({ role: 'tool', tool_call_id: 'call_a' });
        expect(JSON.parse(toolMessage.content).data).toEqual({ modes: [{ angle_deg: 59.5 }] });

        const step = container.querySelector('.llm-step');
        expect(step.className).toContain('is-ok');
        expect(step.querySelector('.llm-step-label').textContent).toBe('Bond angles Ta–Ta–Ta');
        expect(step.querySelector('.llm-step-headline').textContent).toBe('41336 angles.');
        const outcome = container.querySelector('.llm-outcome');
        expect(outcome.className).toContain('is-achieved');
        expect(outcome.textContent).toContain('the Ta4 cluster angle.');
        expect(container.querySelector('.llm-turn-footer').textContent).toBe('1 analysis');

        // "Show in" hands the host the page and what to load there.
        const show = step.querySelector('.llm-step-open');
        expect(show.textContent).toBe('Show in Bond Geometry');
        act(() => show.click());
        expect(onOpenPage).toHaveBeenCalledWith('geometry', { end1: 'Ta', apex: 'Ta', end2: 'Ta', r12Min: 2.7, r12Max: 3.3 });
        // The dock stays open beside the page.
        expect(container.querySelector('.llm-chat-panel').className).not.toContain('is-closed');
    });

    it('collapses on Escape, keeps the conversation, and reopens where it was', async () => {
        const onOpenChange = vi.fn();
        act(() => root.render(chat({ onOpenChange })));
        await openChat();
        await ask('angles?');
        expect(container.querySelectorAll('.llm-message')).toHaveLength(2);

        const panel = container.querySelector('.llm-chat-panel');
        act(() => panel.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true })));
        expect(onOpenChange).toHaveBeenLastCalledWith(false);
        // Hidden, not unmounted: the conversation survives.
        expect(container.querySelector('.llm-chat-panel').className).toContain('is-closed');
        const rail = container.querySelector('.llm-dock-rail');
        expect(document.activeElement).toBe(rail);
        expect(JSON.parse(localStorage.getItem(DOCK_STORAGE_KEY)).open).toBe(false);

        await openChat();
        expect(container.querySelectorAll('.llm-message')).toHaveLength(2);
        expect(fetchMock.mock.calls.filter(([url]) => String(url).endsWith('/v1/models'))).toHaveLength(1);
    });

    it('starts a new conversation when another run opens', async () => {
        act(() => root.render(chat()));
        await openChat();
        await ask('q');
        expect(container.querySelectorAll('.llm-message')).toHaveLength(2);

        act(() => root.render(chat({ runKey: 'run-2' })));
        await settle();
        expect(container.querySelectorAll('.llm-message')).toHaveLength(0);
        expect(container.textContent).toContain('What should I check?');
    });

    it('survives an empty reply and an error sent mid-stream', async () => {
        // Seen with LM Studio: a GPU out-of-memory ends the stream with an
        // error event and no content. That used to end the turn as an empty
        // reply, which was dropped, and the chat then crashed rendering.
        // The empty reply gets one reminder, which is empty too.
        const queue = [
            sse([delta({}, 'stop')]),
            sse([delta({}, 'stop')]),
            sse([{ error: { message: 'Insufficient Memory (kIOGPUCommandBufferCallbackErrorOutOfMemory)' } }])
        ];
        vi.stubGlobal('fetch', vi.fn(async (url) => {
            if (isServerProbe(url)) return notFound;
            return String(url).endsWith('/v1/models') ? modelList : queue.shift();
        }));
        act(() => root.render(chat()));
        await openChat();

        await ask('first?');
        expect(container.querySelectorAll('.llm-message')).toHaveLength(2);
        expect(container.querySelector('.llm-turn-footer').textContent)
            .toBe('The model returned an empty reply: the model wrote no answer, even after a reminder.');

        await ask('second?');
        expect(container.querySelector('[role="alert"]').textContent)
            .toBe('The server stopped the reply: Insufficient Memory (kIOGPUCommandBufferCallbackErrorOutOfMemory) '
                + 'The model server ran out of memory. Close other models (or other jobs using the GPU), load this one '
                + 'with a smaller context, or pick a smaller model.');
        expect(container.querySelector('.llm-composer-input').disabled).toBe(false);
    });

    it('sends the earlier question and its answer with the next question, and forgets them in a new conversation', async () => {
        stubReplies([
            sse([delta({ content: 'Not converged.\n\nOutcome: achieved — still falling.' }, 'stop')]),
            sse([delta({ content: 'Six Se.\n\nOutcome: achieved — counted.' }, 'stop')]),
            sse([delta({ content: 'Fresh.\n\nOutcome: achieved — x' }, 'stop')])
        ]);
        act(() => root.render(chat()));
        await openChat();
        await ask('Has it converged?');
        await ask('And around Ta?');
        expect(chatBodies).toHaveLength(2);
        // system, context, acknowledgement, then the memory, then the question.
        expect(chatBodies[1].messages.slice(3)).toEqual([
            { role: 'user', content: 'Has it converged?' },
            { role: 'assistant', content: 'Not converged.\n\nOutcome: achieved — still falling.' },
            { role: 'user', content: 'And around Ta?' }
        ]);
        act(() => container.querySelector('.llm-link-btn').click());
        await ask('Again?');
        expect(chatBodies[2].messages.slice(3)).toEqual([{ role: 'user', content: 'Again?' }]);
    });

    describe('a new configuration of the run (configKey: a Live Data save)', () => {
        const DIVIDER = 'New configuration saved — earlier answers were about the previous one';
        // The conversation log, one word per entry.
        const logEntries = () => [...container.querySelector('.llm-chat-log').children].map((element) => {
            if (element.classList.contains('llm-config-divider')) return 'divider';
            return element.classList.contains('is-user') ? 'user' : 'assistant';
        });

        it('keeps the conversation on screen, ends it with a divider, and sends none of it as memory', async () => {
            stubReplies([
                sse([delta({ content: 'Not converged.\n\nOutcome: achieved — still falling.' }, 'stop')]),
                sse([delta({ content: 'Six Se.\n\nOutcome: achieved — counted.' }, 'stop')]),
                sse([delta({ content: 'Still six.\n\nOutcome: achieved — counted.' }, 'stop')])
            ]);
            act(() => root.render(chat({ configKey: 'epoch|0' })));
            await openChat();
            await ask('Has it converged?');

            // The same configuration (a re-render, a save of the plots only): nothing changes.
            act(() => root.render(chat({ configKey: 'epoch|0' })));
            await settle();
            expect(logEntries()).toEqual(['user', 'assistant']);

            // A new configuration: the answers stay, the divider ends them, and
            // the screen reader hears it.
            act(() => root.render(chat({ configKey: 'epoch|1' })));
            await settle();
            expect(logEntries()).toEqual(['user', 'assistant', 'divider']);
            expect(container.querySelector('.llm-config-divider').textContent).toBe(DIVIDER);
            expect(container.querySelector('.llm-sr-only[aria-live]').textContent).toContain(DIVIDER);
            // A second save before any question: still one divider.
            act(() => root.render(chat({ configKey: 'epoch|2' })));
            await settle();
            expect(logEntries()).toEqual(['user', 'assistant', 'divider']);

            // The next question is sent as in a fresh conversation: no memory
            // of the answer about the previous configuration.
            await ask('And around Ta?');
            expect(chatBodies[1].messages.slice(3)).toEqual([{ role: 'user', content: 'And around Ta?' }]);
            expect(logEntries()).toEqual(['user', 'assistant', 'divider', 'user', 'assistant']);

            // Memory holds the exchanges about the current configuration, and only them.
            await ask('Still?');
            expect(chatBodies[2].messages.slice(3)).toEqual([
                { role: 'user', content: 'And around Ta?' },
                { role: 'assistant', content: 'Six Se.\n\nOutcome: achieved — counted.' },
                { role: 'user', content: 'Still?' }
            ]);

            // Another run still starts a new conversation, divider and all.
            act(() => root.render(chat({ runKey: 'run-2', configKey: 'epoch|2' })));
            await settle();
            expect(container.querySelector('.llm-chat-log')).toBeNull();
            expect(container.querySelector('.llm-config-divider')).toBeNull();
        });

        it('puts the divider after an answer still running at the save: it was asked about the previous one', async () => {
            let release;
            stubReplies([
                () => new Promise((resolve) => {
                    release = () => resolve(sse([delta({ content: 'Not converged.\n\nOutcome: achieved — x' }, 'stop')]));
                }),
                sse([delta({ content: 'Six Se.\n\nOutcome: achieved — counted.' }, 'stop')])
            ]);
            act(() => root.render(chat({ configKey: 'epoch|0' })));
            await openChat();
            await ask('Has it converged?');
            await settleUntil(() => typeof release === 'function');
            expect(container.querySelector('.llm-send.is-stop')).not.toBeNull();

            act(() => root.render(chat({ configKey: 'epoch|1' })));
            await settle();
            expect(logEntries()).toEqual(['user', 'assistant', 'divider']);

            await act(async () => { release(); });
            await settleUntil(() => !container.querySelector('.llm-send.is-stop'));
            expect(logEntries()).toEqual(['user', 'assistant', 'divider']);
            expect(container.textContent).toContain('Not converged.');

            await ask('And around Ta?');
            expect(chatBodies[1].messages.slice(3)).toEqual([{ role: 'user', content: 'And around Ta?' }]);
            expect(logEntries()).toEqual(['user', 'assistant', 'divider', 'user', 'assistant']);
        });
    });

    it('says why it cannot ask when the host hands it a tool the turn adds itself', async () => {
        const shadow = { ...bondAngles, name: 'ask_user', description: 'A host ask_user.' };
        act(() => root.render(chat({ tools: [shadow] })));
        await openChat();
        // The panel stays up and says what is wrong; asking fails with the same reason.
        expect(container.querySelector('[role="alert"]').textContent).toMatch(/tools must not contain ask_user/);
        await ask('angles?');
        expect(chatBodies).toHaveLength(0);
        expect(container.querySelector('[role="alert"]').textContent).toMatch(/tools must not contain ask_user/);
        expect(container.querySelector('.llm-composer-input').disabled).toBe(false);
    });

    it('says when analyses are off, and asks for a run when none is open', async () => {
        saveSettings({ copilotTools: false });
        act(() => root.render(chat()));
        await openChat();
        expect(container.querySelector('.llm-composer-status').textContent).toMatch(/Analyses are off/);

        act(() => root.render(chat({ plotFiles: [], structure: null })));
        await settle();
        expect(container.querySelector('.llm-chat-norun').textContent).toMatch(/Open a run folder/);
    });
    // A queued reply may be a function of the request body: a model that
    // reads what it was sent.
    const stubReplies = (queue) => vi.stubGlobal('fetch', vi.fn(async (url, options = {}) => {
        if (isServerProbe(url)) return notFound;
        if (String(url).endsWith('/v1/models')) return modelList;
        const body = JSON.parse(options.body);
        chatBodies.push(body);
        const next = queue.shift();
        return typeof next === 'function' ? next(body) : next;
    }));

    it('remembers the open state and the width, and resizes from the keyboard', async () => {
        const width = Object.getOwnPropertyDescriptor(window, 'innerWidth');
        Object.defineProperty(window, 'innerWidth', { value: 1600, configurable: true });
        try {
            localStorage.setItem(DOCK_STORAGE_KEY, JSON.stringify({ open: true, width: 500 }));
            const onOpenChange = vi.fn();
            act(() => root.render(chat({ onOpenChange })));
            await settle();
            const dock = container.querySelector('.llm-dock');
            expect(dock.dataset.open).toBe('true');
            expect(onOpenChange).toHaveBeenCalledWith(true);
            expect(dock.style.getPropertyValue('--llm-dock-w')).toBe('500px');

            const handle = container.querySelector('.llm-dock-resize');
            expect(handle.getAttribute('aria-valuenow')).toBe('500');
            act(() => handle.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowLeft', bubbles: true })));
            expect(dock.style.getPropertyValue('--llm-dock-w')).toBe('524px');
            expect(JSON.parse(localStorage.getItem(DOCK_STORAGE_KEY)).width).toBe(524);
            // End: the widest dock, 760px here (a 1600px window keeps the pages
            // 840px); on a 1200px window the pages' 760px would cap it at 440.
            act(() => handle.dispatchEvent(new KeyboardEvent('keydown', { key: 'End', bubbles: true })));
            expect(dock.style.getPropertyValue('--llm-dock-w')).toBe('760px');
            Object.defineProperty(window, 'innerWidth', { value: 1200, configurable: true });
            act(() => window.dispatchEvent(new Event('resize')));
            expect(dock.style.getPropertyValue('--llm-dock-w')).toBe('440px');

            // A drag moves the edge live and is remembered on release.
            const pointer = (type, clientX) => {
                const event = new MouseEvent(type, { bubbles: true, clientX, button: 0 });
                Object.defineProperty(event, 'pointerId', { value: 1 });
                return event;
            };
            act(() => handle.dispatchEvent(pointer('pointerdown', 800)));
            act(() => handle.dispatchEvent(pointer('pointermove', 850)));
            expect(dock.style.getPropertyValue('--llm-dock-w')).toBe('390px');
            act(() => handle.dispatchEvent(pointer('pointerup', 850)));
            expect(JSON.parse(localStorage.getItem(DOCK_STORAGE_KEY)).width).toBe(390);
            expect(handle.getAttribute('aria-valuenow')).toBe('390');
        } finally {
            if (width) Object.defineProperty(window, 'innerWidth', width);
        }
    });

    it('marks the rail unread when an answer lands while collapsed', async () => {
        act(() => root.render(chat()));
        await openChat();
        act(() => setTextarea(container.querySelector('.llm-composer-input'), 'angles?'));
        act(() => container.querySelector('.llm-composer').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true })));
        // Collapse while the answer is still running.
        act(() => container.querySelector('.llm-chat-panel').dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true })));
        await settle();
        const rail = container.querySelector('.llm-dock-rail');
        expect(rail.className).toContain('has-unread');
        expect(rail.getAttribute('aria-label')).toBe('AI Copilot — new answer');
        await openChat();
        expect(container.querySelectorAll('.llm-message')).toHaveLength(2);
        act(() => container.querySelector('.llm-chat-panel').dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true })));
        expect(container.querySelector('.llm-dock-rail').className).not.toContain('has-unread');
    });

    describe('the Workbench asks first', () => {
        // A TaSe6-like split Ta–Se shell (peaks at 2.52 and 2.74 Å) in a
        // 10 Å cubic cell, as the partials file and the structure carry it.
        const grid = Array.from({ length: 500 }, (_, i) => i * 0.02);
        const shell = (r) => 25 * Math.exp(-((r - 2.52) ** 2) / 0.0072) + 12 * Math.exp(-((r - 2.74) ** 2) / 0.0288) + (r > 4.5 ? 1 : 0);
        const partials = { name: 'run_PDFpartials.csv', plotKind: 'pdf_partials', plotData: { series: [{ label: 'Ta-Se', x: grid, y: grid.map(shell) }] } };
        const structure = { elementCounts: { Ta: 16, Se: 32 }, latticeVectors: [[10, 0, 0], [0, 10, 0], [0, 0, 10]], supercell: [1, 1, 1] };

        it('shows its question before the model is asked, and the model gets the choice', async () => {
            stubReplies([sse([delta({ content: 'Computed.\n\nOutcome: achieved — x' }, 'stop')])]);
            act(() => root.render(chat({ plotFiles: [partials], structure })));
            await openChat();
            await ask('What are the Se–Ta–Se bond angles?');

            // Waiting on the card: no request has gone to the model yet.
            expect(chatBodies).toHaveLength(0);
            const card = container.querySelector('.llm-ask.is-waiting');
            expect(card.querySelector('.llm-ask-question').textContent).toBe('Which Ta–Se window should the analysis use?');
            expect(card.querySelector('.llm-step-meta').textContent).toBe('The Workbench asks before the analysis');
            const options = [...card.querySelectorAll('.llm-ask-option')];
            expect(options.map((button) => button.querySelector('.llm-ask-option-label').firstChild.textContent))
                .toEqual(['Whole first shell', 'First peak only']);

            act(() => options[1].click());
            await settleUntil(() => container.querySelector('.llm-outcome'));
            expect(chatBodies).toHaveLength(1);
            expect(chatBodies[0].messages.at(-1).content).toMatch(
                /^What are the Se–Ta–Se bond angles\?\n\n\(The Workbench asked: "Which Ta–Se window should the analysis use\?" — the user chose "First peak only" for the Ta–Se bonds: the window /
            );
            expect(container.querySelector('.llm-ask').className).toContain('is-answered');
            // The focus is not dropped to the page's start: the answered card
            // keeps it while the turn runs, then the composer (review 2026-10-12 #12).
            expect(document.activeElement).toBe(container.querySelector('.llm-composer-input'));
            // What the user sees as their question is what they typed.
            expect(container.querySelector('.llm-message.is-user .llm-bubble').textContent).toBe('What are the Se–Ta–Se bond angles?');
        });
    });

    describe('the answer check', () => {
        // A still-falling ln(χ²) history: the watchdog says improving.
        const rValueFile = { name: 'run.log', plotData: { chiColumn: 'X_ray_(R)1', series: [{ y: Array.from({ length: 60 }, (_, i) => 5 - i * 0.05) }] } };
        const convergenceCall = () => sse([delta({ tool_calls: [{ index: 0, id: 'call_c', type: 'function', function: { name: 'assess_convergence', arguments: '{}' } }] }, 'tool_calls')]);

        it('lowers a contradicted answer to "Not verified" and says why', async () => {
            stubReplies([
                convergenceCall(),
                sse([delta({ content: 'The run has converged — a plateau.\n\nOutcome: achieved — flat history.' }, 'stop')])
            ]);
            act(() => root.render(chat({ rValueFile })));
            await openChat();
            await ask('Has it converged?');
            const audit = container.querySelector('.llm-audit');
            expect(audit.getAttribute('aria-label')).toBe('Checks on this answer');
            expect(audit.textContent).toMatch(/^CheckThe answer says the run has converged, but the χ² history is still improving/);
            // The model's own "achieved" does not stand.
            const chip = container.querySelector('.llm-outcome');
            expect(chip.className).toContain('is-not-verified');
            expect(chip.querySelector('strong').textContent).toBe('Not verified');
            expect(chip.textContent).toContain('a check below contradicts the answer. The model said "achieved": flat history.');
            expect(container.querySelector('.llm-sr-only[aria-live]').textContent).toMatch(/^Outcome: Not verified\. .* Check: The answer says the run has converged/);
        });

        it('keeps "Achieved" on an answer the checks agree with', async () => {
            stubReplies([
                convergenceCall(),
                sse([delta({ content: 'Not converged: ln χ² is still improving.\n\nOutcome: achieved — the watchdog says improving.' }, 'stop')])
            ]);
            act(() => root.render(chat({ rValueFile })));
            await openChat();
            await ask('Has it converged?');
            expect(container.querySelector('.llm-audit')).toBeNull();
            const chip = container.querySelector('.llm-outcome');
            expect(chip.className).toContain('is-achieved');
            expect(chip.querySelector('strong').textContent).toBe('Achieved');
        });
    });

    describe('skills', () => {
        const skillText = [
            '---', 'name: rmc-fit-assessment', 'description: Judge a run before trusting it.', '---', '',
            '# RMC fit assessment', '', '**Status: Draft: awaiting the user\'s interview.**', '',
            '## Steps', '', '1. Call assess_convergence (to confirm with the user).'
        ].join('\n');
        const skills = { files: { '.claude/skills/rmc-fit-assessment/SKILL.md': skillText }, references: [], load: vi.fn() };
        // A still-falling ln(χ²) history: the watchdog says improving.
        const rValueFile = { name: 'run.log', plotData: { chiColumn: 'X_ray_(R)1', series: [{ y: Array.from({ length: 60 }, (_, i) => 5 - i * 0.05) }] } };

        it('lists the skills in the prompt, reads one as a step that is not an analysis, and checks the answer', async () => {
            stubReplies([
                sse([delta({ tool_calls: [{ index: 0, id: 'call_s', type: 'function', function: { name: 'read_skill', arguments: '{"name":"rmc-fit-assessment"}' } }] }, 'tool_calls')]),
                sse([delta({ content: 'The run has converged — a plateau.\n\nOutcome: achieved — flat history.' }, 'stop')])
            ]);
            act(() => root.render(chat({ skills, rValueFile })));
            await openChat();
            expect(container.querySelector('.llm-composer-status').textContent).toMatch(/following 1 skill\.$/);
            await ask('Has it converged?');

            const system = chatBodies[0].messages[0].content;
            expect(system).toContain('- rmc-fit-assessment: Judge a run before trusting it.');
            expect(system).not.toContain('Call assess_convergence (to confirm');
            expect(chatBodies[0].tools.map((tool) => tool.function.name)).toContain('read_skill');
            expect(JSON.parse(chatBodies[1].messages.at(-1).content).data.text).toContain('Call assess_convergence');

            expect(container.querySelector('.llm-step-label').textContent).toBe('Reading the rmc-fit-assessment skill');
            expect(container.querySelector('.llm-turn-footer')).toBeNull();
            // The answer contradicts the watchdog: a check under it says so.
            const audit = container.querySelector('.llm-audit');
            expect(audit.getAttribute('aria-label')).toBe('Checks on this answer');
            expect(audit.textContent).toMatch(/^CheckThe answer says the run has converged, but the χ² history is still improving/);
        });

        it('says so when a skill file cannot be read', async () => {
            act(() => root.render(chat({ skills: { files: { '.claude/skills/x/SKILL.md': 'no frontmatter' }, references: [] } })));
            await openChat();
            expect(container.querySelector('[role="alert"]').textContent).toMatch(/The Copilot's skills could not be read: .*has no frontmatter/);
        });
    });

    describe('the context guard on Ollama', () => {
        const ndjson = (objects) => {
            const encoder = new TextEncoder();
            return { ok: true, status: 200, body: new ReadableStream({
                start(controller) {
                    controller.enqueue(encoder.encode(objects.map((object) => `${JSON.stringify(object)}\n`).join('')));
                    controller.close();
                }
            }) };
        };
        // Ollama with gemma4:26b loaded at a 4,096-token default context.
        const ollama = ({ max = 262144, replies = [] } = {}) => {
            const chats = [];
            const fetchMock = vi.fn(async (url, options = {}) => {
                const path = new URL(String(url)).pathname;
                if (path === '/v1/models') return { ok: true, status: 200, json: async () => ({ data: [{ id: 'gemma4:26b' }] }) };
                if (path === '/api/version') return { ok: true, status: 200, json: async () => ({ version: '0.40.2' }) };
                if (path === '/api/ps') return { ok: true, status: 200, json: async () => ({ models: [{ name: 'gemma4:26b', model: 'gemma4:26b', context_length: 4096 }] }) };
                if (path === '/api/show') return { ok: true, status: 200, json: async () => ({ model_info: { 'gemma4.context_length': max } }) };
                if (path === '/api/chat') {
                    chats.push(JSON.parse(options.body));
                    return replies.shift();
                }
                return notFound;
            });
            vi.stubGlobal('fetch', fetchMock);
            return { chats, fetchMock };
        };

        it('runs the turn through the native API with a context that holds the whole prompt', async () => {
            const angles = { ...bondAngles, run: vi.fn(async () => ({ data: { modes: [{ angle_deg: 60 }] }, checks: [{ level: 'ok', text: 'angles' }] })) };
            const { chats, fetchMock } = ollama({ replies: [
                ndjson([
                    { message: { role: 'assistant', content: 'Computing.', tool_calls: [{ function: { name: 'compute_bond_angles', arguments: { end1: 'Ta', apex: 'Ta', end2: 'Ta', r12_min: 2.7, r12_max: 3.3 } } }] }, done: false },
                    { message: { content: '' }, done: true, done_reason: 'stop', prompt_eval_count: 3000, eval_count: 30 }
                ]),
                ndjson([
                    { message: { content: 'They peak at 60°.\n\nOutcome: achieved — x' }, done: false },
                    { message: { content: '' }, done: true, done_reason: 'stop', prompt_eval_count: 3300, eval_count: 20 }
                ])
            ] });
            saveSettings({ model: 'gemma4:26b' });
            act(() => root.render(chat({ tools: [angles] })));
            await openChat();
            await ask('What are the Ta–Ta–Ta angles?');
            await settleUntil(() => container.querySelector('.llm-outcome'));

            expect(chats).toHaveLength(2);
            // Not the 4,096 the model is loaded with: room for the prompt and the reply.
            for (const body of chats) {
                expect(body.options.num_ctx).toBeGreaterThanOrEqual(16384);
                expect(body.options.num_predict).toBeGreaterThan(8000);
            }
            expect(chats[1].options.num_ctx).toBe(chats[0].options.num_ctx);
            expect(chats[1].messages.at(-1)).toMatchObject({ role: 'tool', tool_name: 'compute_bond_angles' });
            expect(angles.run).toHaveBeenCalledTimes(1);
            expect(fetchMock.mock.calls.some(([url]) => String(url).endsWith('/v1/chat/completions'))).toBe(false);
            expect(container.querySelector('.llm-outcome').className).toContain('is-achieved');
        });

        it('refuses the turn, with a reason, when the model cannot hold the prompt', async () => {
            const { chats } = ollama({ max: 8192 });
            saveSettings({ model: 'gemma4:26b' });
            act(() => root.render(chat()));
            await openChat();
            await ask('Has the run converged?');
            expect(chats).toHaveLength(0);
            expect(container.querySelector('[role="alert"]').textContent)
                .toMatch(/needs about [\d,]+ tokens of context .* but gemma4:26b holds at most 8,192\. .*Nothing was sent/);
            expect(container.querySelector('.llm-composer-input').disabled).toBe(false);
        });

        it('says when an answer stopped at the length limit', async () => {
            ollama({ replies: [ndjson([
                { message: { content: 'The angles are' }, done: false },
                { message: { content: '' }, done: true, done_reason: 'length', prompt_eval_count: 3000, eval_count: 13000 }
            ])] });
            saveSettings({ model: 'gemma4:26b' });
            act(() => root.render(chat()));
            await openChat();
            await ask('Angles?');
            const footer = container.querySelector('.llm-turn-footer');
            expect(footer.className).toContain('is-cut');
            expect(footer.textContent).toContain('the reply reached the model\'s length limit and is cut short');
        });
    });

    describe('asking the user (ask_user)', () => {
        const askCall = (args) => sse([
            delta({ content: 'The Ta–Se shell is split.' }),
            delta({ tool_calls: [{ index: 0, id: 'call_q', type: 'function', function: { name: 'ask_user', arguments: JSON.stringify(args) } }] }, 'tool_calls')
        ]);
        const shellQuestion = {
            question: 'Which Ta–Se window should the Se–Ta–Se angles use?',
            options: [
                { label: 'Whole first shell', description: '2.33–2.96 Å, both Ta–Se peaks' },
                { label: 'First peak only', description: '2.33–2.71 Å, the short bonds' }
            ],
            default: 'Whole first shell'
        };
        const WINDOWS = { 'Whole first shell': [2.33, 2.96], 'First peak only': [2.33, 2.71] };
        // The model computes in the window the answer names.
        const computeFromAnswer = (body) => {
            const { data } = JSON.parse(body.messages.at(-1).content);
            const typed = /(\d+(?:\.\d+)?)\D+(\d+(?:\.\d+)?)/.exec(data.answer || '');
            const [min, max] = WINDOWS[data.answer] || [Number(typed[1]), Number(typed[2])];
            return sse([delta({ tool_calls: [{ index: 0, id: 'call_b', type: 'function', function: {
                name: 'compute_bond_angles',
                arguments: JSON.stringify({ end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: min, r12_max: max })
            } }] }, 'tool_calls')]);
        };
        const answered = () => sse([delta({ content: 'The Se–Ta–Se angles …\n\nOutcome: achieved — in the chosen window.' }, 'stop')]);
        const freshAngles = () => ({ ...bondAngles, run: vi.fn(async () => ({ data: { modes: [{ angle_deg: 90 }] }, checks: [{ level: 'ok', text: 'angles' }] })) });

        const card = () => container.querySelector('.llm-ask');
        const options = () => [...container.querySelectorAll('.llm-ask-option')];

        it('shows the question as a card, waits, and sends the chosen option back to the model', async () => {
            stubReplies([askCall(shellQuestion), computeFromAnswer, answered()]);
            const angles = freshAngles();
            act(() => root.render(chat({ tools: [angles] })));
            await openChat();
            await ask('What are the Se–Ta–Se angles?');

            // The prompt carries the rules, and the tool is offered.
            expect(chatBodies[0].messages[0].content).toContain('Ask instead of guessing a parameter that changes the science');
            expect(chatBodies[0].tools.map((tool) => tool.function.name)).toContain('ask_user');
            // The turn waits on the card: one request so far.
            expect(chatBodies).toHaveLength(1);
            expect(card().className).toContain('is-waiting');
            expect(card().getAttribute('role')).toBe('group');
            expect(card().querySelector('.llm-ask-question').textContent).toBe(shellQuestion.question);
            expect(options().map((button) => button.textContent)).toEqual([
                '1Whole first shelldefault2.33–2.96 Å, both Ta–Se peaks',
                '2First peak only2.33–2.71 Å, the short bonds'
            ]);
            expect(container.querySelector('.llm-ask-other-input').placeholder).toBe('Other…');
            // The empty composer gives focus to the default option; arrows move it.
            expect(document.activeElement).toBe(options()[0]);
            act(() => options()[0].dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true })));
            expect(document.activeElement).toBe(options()[1]);
            expect(container.querySelector('.llm-sr-only[aria-live]').textContent).toBe(`The Copilot asks: ${shellQuestion.question}`);
            expect(container.querySelector('.llm-send').getAttribute('aria-label')).toBe('Stop');

            // Picking "First peak only" with its number key answers it.
            act(() => options()[1].dispatchEvent(new KeyboardEvent('keydown', { key: '2', bubbles: true })));
            await settleUntil(() => container.querySelector('.llm-outcome'));
            const result = JSON.parse(chatBodies[1].messages.at(-1).content);
            expect(result.data).toEqual({ question: shellQuestion.question, answer: 'First peak only', means: '2.33–2.71 Å, the short bonds' });
            expect(angles.run).toHaveBeenCalledTimes(1);
            expect(angles.run.mock.calls[0][0]).toMatchObject({ r12_min: 2.33, r12_max: 2.71 });

            // The card keeps the question and the reply; the answer follows.
            expect(card().className).toContain('is-answered');
            expect(card().querySelector('.llm-ask-reply').textContent).toBe('YouFirst peak only');
            expect(container.querySelector('.llm-outcome').className).toContain('is-achieved');
            // The question is not counted as an analysis.
            expect(container.querySelector('.llm-turn-footer').textContent).toBe('1 analysis');
        });

        it('takes an answer typed in Other', async () => {
            stubReplies([askCall(shellQuestion), computeFromAnswer, answered()]);
            const angles = freshAngles();
            act(() => root.render(chat({ tools: [angles] })));
            await openChat();
            await ask('Se–Ta–Se angles?');
            const input = container.querySelector('.llm-ask-other-input');
            expect(container.querySelector('.llm-ask-send').disabled).toBe(true);
            act(() => {
                Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set.call(input, '2.3 to 3.1 Å');
                input.dispatchEvent(new Event('input', { bubbles: true }));
            });
            act(() => container.querySelector('.llm-ask-other').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true })));
            await settleUntil(() => container.querySelector('.llm-outcome'));
            expect(JSON.parse(chatBodies[1].messages.at(-1).content).data).toEqual({ question: shellQuestion.question, answer: '2.3 to 3.1 Å', typed: true });
            expect(angles.run.mock.calls[0][0]).toMatchObject({ r12_min: 2.3, r12_max: 3.1 });
        });

        it('cancels cleanly: Stop while waiting ends the turn, and the next question works', async () => {
            stubReplies([askCall(shellQuestion), sse([delta({ content: 'Fine.\n\nOutcome: achieved — x' }, 'stop')])]);
            act(() => root.render(chat({ tools: [freshAngles()] })));
            await openChat();
            await ask('Se–Ta–Se angles?');
            expect(card().className).toContain('is-waiting');
            act(() => container.querySelector('.llm-send').click());
            await settle();
            expect(chatBodies).toHaveLength(1);
            expect(card().className).toContain('is-stopped');
            expect(card().textContent).toContain('Not answered');
            expect(container.querySelector('.llm-ask-option')).toBeNull();
            expect(container.querySelector('.llm-turn-footer').textContent).toBe('Stopped.');
            expect(container.querySelector('[role="alert"]')).toBeNull();
            await ask('Anything else?');
            expect(chatBodies).toHaveLength(2);
            expect(container.querySelector('.llm-outcome').className).toContain('is-achieved');
        });

        it('drops a waiting question when another run opens', async () => {
            stubReplies([askCall(shellQuestion)]);
            act(() => root.render(chat({ tools: [freshAngles()] })));
            await openChat();
            await ask('Se–Ta–Se angles?');
            expect(card()).not.toBeNull();
            act(() => root.render(chat({ tools: [freshAngles()], runKey: 'run-2' })));
            await settle();
            expect(card()).toBeNull();
            expect(container.querySelectorAll('.llm-message')).toHaveLength(0);
            expect(container.querySelector('[role="alert"]')).toBeNull();
            expect(chatBodies).toHaveLength(1);
        });

        it('marks the rail when a question arrives while the Copilot is collapsed', async () => {
            stubReplies([askCall(shellQuestion), computeFromAnswer, answered()]);
            act(() => root.render(chat({ tools: [freshAngles()] })));
            await openChat();
            act(() => setTextarea(container.querySelector('.llm-composer-input'), 'Se–Ta–Se angles?'));
            act(() => container.querySelector('.llm-composer').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true })));
            act(() => container.querySelector('.llm-chat-panel').dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true })));
            await settle();
            const rail = container.querySelector('.llm-dock-rail');
            expect(rail.className).toContain('has-unread');
            expect(rail.getAttribute('aria-label')).toBe('AI Copilot — a question for you');
            await openChat();
            act(() => options()[0].click());
            await settleUntil(() => container.querySelector('.llm-outcome'));
            expect(JSON.parse(chatBodies[1].messages.at(-1).content).data.answer).toBe('Whole first shell');
            expect(container.querySelector('.llm-outcome')).not.toBeNull();
        });
    });

    describe('pointing at the host cards', () => {
        let host;
        const cards = [
            { id: 'chi2-history', page: 'dashboard', title: 'χ² history', kind: 'r_value', marks: true, about: 'ln χ² against time' },
            { id: 'bond-angles', page: 'geometry', title: 'Bond-angle distribution', about: 'the angle histogram' }
        ];

        beforeEach(() => {
            host = document.createElement('div');
            host.innerHTML = '<div id="dash"><article data-copilot-card="chi2-history">chart</article></div>'
                + '<div id="geom" aria-hidden="true"><article data-copilot-card="bond-angles">angles</article></div>';
            document.body.appendChild(host);
        });
        afterEach(() => host.remove());

        const wait = (ms) => act(() => new Promise((resolve) => setTimeout(resolve, ms)));

        it('highlights a card with marks, cites it as a chip, and clears on a click elsewhere', async () => {
            stubReplies([
                sse([
                    delta({ content: 'Here is the history.' }),
                    delta({ tool_calls: [{ index: 0, id: 'call_h', type: 'function', function: { name: 'highlight_cards', arguments: '{"cards":["chi2-history"],"marks":[{"x":120,"label":"plateau"}]}' } }] }, 'tool_calls')
                ]),
                sse([delta({ content: 'The [χ² history](card:chi2-history) still falls; see [[bond-angles]].\n\nOutcome: achieved — read the card.' }, 'stop')])
            ]);
            const onMarks = vi.fn();
            const onOpenPage = vi.fn((page) => { if (page === 'geometry') host.querySelector('#geom').removeAttribute('aria-hidden'); });
            act(() => root.render(chat({ cards, onMarks, onOpenPage })));
            await openChat();
            await ask('Has it converged?');

            // The view tools are offered, and the prompt says how to point.
            expect(chatBodies[0].tools.map((tool) => tool.function.name)).toContain('highlight_cards');
            expect(chatBodies[0].messages[0].content).toContain('[short label](card:<id>)');
            const card = host.querySelector('[data-copilot-card="chi2-history"]');
            expect(card.getAttribute('data-copilot-focus')).toBe('primary');
            expect(onMarks).toHaveBeenLastCalledWith({ 'chi2-history': [{ x: 120, label: 'plateau' }] });
            expect(container.querySelector('.llm-step-label').textContent).toBe('Pointing at χ² history');
            const step = container.querySelector('.llm-step');
            expect(step.className).toContain('is-ok');

            const chips = [...container.querySelectorAll('.llm-card-chip')];
            expect(chips.map((chip) => chip.textContent)).toEqual(['χ² history', 'Bond-angle distribution']);

            // A click on the page (not in the Copilot) ends the highlight.
            act(() => host.querySelector('#dash').dispatchEvent(new Event('click', { bubbles: true })));
            expect(card.hasAttribute('data-copilot-focus')).toBe(false);
            expect(onMarks).toHaveBeenLastCalledWith({});

            // The chip brings it back; a chip for a card on another page opens it.
            act(() => chips[0].click());
            await wait(20);
            expect(card.getAttribute('data-copilot-focus')).toBe('primary');
            act(() => chips[1].click());
            await wait(150);
            expect(onOpenPage).toHaveBeenCalledWith('geometry', null);
            expect(host.querySelector('[data-copilot-card="bond-angles"]').getAttribute('data-copilot-focus')).toBe('primary');
            expect(card.hasAttribute('data-copilot-focus')).toBe(false);
            // Collapsing the Copilot ends the highlight too.
            act(() => container.querySelector('.llm-chat-panel').dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true })));
            expect(host.querySelector('[data-copilot-focus]')).toBeNull();
        });

        it('lets the card of a running analysis breathe until the step ends', async () => {
            let release;
            const gate = new Promise((resolve) => { release = resolve; });
            const convergence = {
                name: 'slow_convergence',
                description: 'Convergence.',
                parameters: { type: 'object', properties: {} },
                cards: () => ['chi2-history'],
                run: vi.fn(async () => { await gate; return { data: { status: 'improving' }, checks: [] }; })
            };
            stubReplies([
                sse([delta({ tool_calls: [{ index: 0, id: 'call_c', type: 'function', function: { name: 'slow_convergence', arguments: '{}' } }] }, 'tool_calls')]),
                sse([delta({ content: 'Improving.\n\nOutcome: achieved — x' }, 'stop')])
            ]);
            act(() => root.render(chat({ cards, tools: [convergence] })));
            await openChat();
            act(() => setTextarea(container.querySelector('.llm-composer-input'), 'converged?'));
            act(() => container.querySelector('.llm-composer').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true })));
            await settle();
            const card = host.querySelector('[data-copilot-card="chi2-history"]');
            expect(card.getAttribute('data-copilot-focus')).toBe('running');
            release();
            await settle();
            expect(card.hasAttribute('data-copilot-focus')).toBe(false);
        });
    });
});
