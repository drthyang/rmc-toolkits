// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The Copilot chatbox end to end against a mocked OpenAI-compatible server:
// nothing is fetched until the launcher opens it; the model asks for a tool,
// the chatbox runs it and shows a step card, the result reaches the next
// request, and the answer ends in an Outcome chip. "Show in <page>" hands the
// host the step's page and parameters. Closing keeps the conversation. No
// live model is involved.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import React, { act } from 'react';
import { createRoot } from 'react-dom/client';
import CopilotChat from '../components/CopilotChat';
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

    it('is a launcher until opened, and contacts no server before that', async () => {
        const onOpenChange = vi.fn();
        act(() => root.render(chat({ onOpenChange })));
        await settle();
        const launcher = container.querySelector('.llm-chat-launcher');
        expect(launcher.textContent).toBe('AI Copilot');
        expect(launcher.getAttribute('aria-expanded')).toBe('false');
        expect(container.querySelector('.llm-chat-panel')).toBeNull();
        expect(fetchMock).not.toHaveBeenCalled();

        await openChat();
        expect(onOpenChange).toHaveBeenLastCalledWith(true);
        const panel = container.querySelector('.llm-chat-panel');
        expect(panel.getAttribute('role')).toBe('dialog');
        expect(container.querySelector('.llm-chat-launcher')).toBeNull();
        expect(panel.textContent).toContain('Connected');
        // The composer takes focus when the chatbox opens.
        expect(document.activeElement).toBe(container.querySelector('.llm-composer-input'));
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
        // The chatbox stays open over the page.
        expect(container.querySelector('.llm-chat-panel').className).not.toContain('is-closed');
    });

    it('closes on Escape, keeps the conversation, and reopens where it was', async () => {
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
        const launcher = container.querySelector('.llm-chat-launcher');
        expect(document.activeElement).toBe(launcher);

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

    it('says when analyses are off, and asks for a run when none is open', async () => {
        saveSettings({ copilotTools: false });
        act(() => root.render(chat()));
        await openChat();
        expect(container.querySelector('.llm-composer-status').textContent).toMatch(/Analyses are off/);

        act(() => root.render(chat({ plotFiles: [], structure: null })));
        await settle();
        expect(container.querySelector('.llm-chat-norun').textContent).toMatch(/Open a run folder/);
    });
});
