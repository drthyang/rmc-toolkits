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
        localStorage.removeItem(DOCK_STORAGE_KEY);
        saveSettings({ baseUrl: 'http://localhost:11434/v1', model: '', apiKey: '', copilotTools: true });
        chatBodies = [];
        const queue = replies();
        fetchMock = vi.fn(async (url, options = {}) => {
            if (String(url).endsWith('/models')) {
                return { ok: true, status: 200, json: async () => ({ data: [{ id: 'tool-model' }] }) };
            }
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
        expect(chatBodies[0].tools.map((tool) => tool.function.name)).toEqual(['read_run_data', 'compute_bond_angles']);
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
        expect(fetchMock.mock.calls.filter(([url]) => String(url).endsWith('/models'))).toHaveLength(1);
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
        const queue = [
            sse([delta({}, 'stop')]),
            sse([{ error: { message: 'Insufficient Memory (kIOGPUCommandBufferCallbackErrorOutOfMemory)' } }])
        ];
        vi.stubGlobal('fetch', vi.fn(async (url) => (String(url).endsWith('/models')
            ? { ok: true, status: 200, json: async () => ({ data: [{ id: 'tool-model' }] }) }
            : queue.shift())));
        act(() => root.render(chat()));
        await openChat();

        await ask('first?');
        expect(container.querySelectorAll('.llm-message')).toHaveLength(2);
        expect(container.querySelector('.llm-turn-footer').textContent).toBe('The model returned an empty reply.');

        await ask('second?');
        expect(container.querySelector('[role="alert"]').textContent)
            .toBe('The server stopped the reply: Insufficient Memory (kIOGPUCommandBufferCallbackErrorOutOfMemory)');
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
    const stubReplies = (queue) => vi.stubGlobal('fetch', vi.fn(async (url, options = {}) => {
        if (String(url).endsWith('/models')) {
            return { ok: true, status: 200, json: async () => ({ data: [{ id: 'tool-model' }] }) };
        }
        chatBodies.push(JSON.parse(options.body));
        return queue.shift();
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
