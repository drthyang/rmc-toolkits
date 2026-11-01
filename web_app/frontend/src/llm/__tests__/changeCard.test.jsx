// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// A change tool through the dock (CopilotChat → useCopilot → the turn), end
// to end against a mocked OpenAI-compatible server: the approval card (one
// line, the details behind "?", Apply / Decline) and the turn waiting on it;
// Ask first by default and the drawer's Auto, which a reload resets; a
// decline the model reads as {"declined": true}; Stop while the card waits;
// a change the page refused after the card (stale); Undo on the applied
// card, and an Undo refused. The change tool here is a test's (effect
// 'change'); the real one on the Bond Geometry page is tested in
// src/__tests__/copilotChangeTools.test.jsx. No live model.

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
const notFound = { ok: false, status: 404, json: async () => ({ error: 'not found' }) };
const modelList = { ok: true, status: 200, json: async () => ({ data: [{ id: 'tool-model' }] }) };

const settle = () => act(async () => {
    for (let i = 0; i < 10; i += 1) await new Promise((resolve) => setTimeout(resolve, 0));
});
const settleUntil = async (ready, rounds = 30) => {
    for (let i = 0; i < rounds && !ready(); i += 1) await settle();
};
const setTextarea = (textarea, value) => {
    Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype, 'value').set.call(textarea, value);
    textarea.dispatchEvent(new Event('input', { bubbles: true }));
};

const LINE = 'Se–Ta–Se · 2.21–2.65 Å · 1° bins';
const DETAILS = ['A–B window: 2.00–3.00 → 2.21–2.65 Å'];

// A change tool as the contract has it: asks the approve port, then reports
// what it did; `stale` makes the page move while the card waits.
const changeTool = ({ stale = false, undo = null } = {}) => ({
    name: 'set_window',
    effect: 'change',
    description: 'Set the window on the page.',
    parameters: { type: 'object', properties: { r12_max: { type: 'number' } }, additionalProperties: false },
    label: () => 'Set the window',
    page: 'geometry',
    run: vi.fn(async (args, { approve, signal }) => {
        const { approved, auto } = await approve({ tool: 'set_window', title: 'Change the window', line: LINE, details: DETAILS }, { signal });
        if (!approved) {
            return {
                data: { declined: true },
                checks: [{ level: 'note', text: 'The user declined this change.' }],
                change: { status: 'declined', line: LINE, details: DETAILS, auto, id: null }
            };
        }
        if (stale) {
            return {
                data: { stale: true },
                checks: [{ level: 'fail', text: 'Not applied: the page changed while the card waited.' }],
                change: { status: 'stale', line: LINE, details: DETAILS, auto, id: null }
            };
        }
        return {
            data: { applied: true, computed: true },
            checks: [{ level: 'ok', text: `Applied: ${LINE}.` }],
            change: { status: 'applied', line: LINE, details: DETAILS, auto, id: 'w1' }
        };
    }),
    undo: undo || vi.fn(async () => ({ ok: true, computed: false }))
});

describe('a change in the dock: the approval card', () => {
    let container;
    let root;
    let chatBodies;

    const callChange = () => sse([
        delta({ content: 'I will set the window.' }),
        delta({ tool_calls: [{ index: 0, id: 'call_c', type: 'function', function: { name: 'set_window', arguments: '{"r12_max":2.65}' } }] }, 'tool_calls')
    ]);
    const answered = () => sse([delta({ content: 'Set.\n\nOutcome: achieved — the page shows it.' }, 'stop')]);
    const stubReplies = (queue) => vi.stubGlobal('fetch', vi.fn(async (url, options = {}) => {
        if (/\/api\//.test(String(url))) return notFound;
        if (String(url).endsWith('/v1/models')) return modelList;
        const body = JSON.parse(options.body);
        chatBodies.push(body);
        const next = queue.shift();
        return typeof next === 'function' ? next(body) : next;
    }));

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        forgetServers();
        localStorage.removeItem(DOCK_STORAGE_KEY);
        saveSettings({ baseUrl: 'http://localhost:11434/v1', model: '', apiKey: '', copilotTools: true });
        chatBodies = [];
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.unstubAllGlobals();
    });

    const chat = (props = {}) => (
        <CopilotChat
            runName="demo"
            plotFiles={[{ plotKind: 'neutron_sq', plotData: { title: 'S(Q)', metrics: { rwp: 0.04 } } }]}
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
    const card = () => container.querySelector('.llm-change');
    const button = (text) => [...card().querySelectorAll('button')].find((item) => item.textContent === text);
    const openSettings = () => act(() => container.querySelector('[aria-label="Connection settings"]').click());
    const mode = (label) => [...container.querySelectorAll('.llm-approval-option')].find((item) => item.textContent === label);

    it('Ask first: the turn waits on the card — one line, the details behind "?" — and Apply makes the change', async () => {
        stubReplies([callChange(), answered()]);
        const tool = changeTool();
        const onOpenPage = vi.fn();
        act(() => root.render(chat({ tools: [tool], onOpenPage })));
        await openChat();
        // Offered with the analyses, and said so in the status line.
        expect(container.querySelector('.llm-composer-status').textContent).toMatch(/changes its inputs once you approve\.$/);
        await ask('Set the Ta–Se window to the whole shell.');

        expect(chatBodies[0].tools.map((item) => item.function.name)).toContain('set_window');
        expect(chatBodies).toHaveLength(1);
        expect(card().className).toContain('is-waiting');
        expect(card().getAttribute('role')).toBe('group');
        expect(card().querySelector('.llm-ask-kicker').textContent).toBe('Bond Geometry changeNeeds your approval');
        expect(card().querySelector('.llm-change-text').textContent).toBe(LINE);
        // The details only behind "?".
        expect(card().querySelector('.llm-change-details')).toBeNull();
        act(() => card().querySelector('.llm-change-help').click());
        expect([...card().querySelectorAll('.llm-change-details li')].map((item) => item.textContent)).toEqual(DETAILS);
        expect(container.querySelector('.llm-sr-only[aria-live]').textContent).toBe(`The Copilot asks to apply a change: ${LINE}`);
        // The waiting card took the focus from the empty composer.
        expect(document.activeElement).toBe(card());
        expect(container.querySelector('.llm-send').getAttribute('aria-label')).toBe('Stop');

        act(() => button('Apply').click());
        await settleUntil(() => container.querySelector('.llm-outcome'));
        expect(tool.run).toHaveBeenCalledTimes(1);
        expect(JSON.parse(chatBodies[1].messages.at(-1).content)).toEqual({
            ref: 'r1', checks: [{ level: 'ok', text: `Applied: ${LINE}.` }], data: { applied: true, computed: true }
        });
        expect(card().className).toContain('is-applied');
        expect(card().querySelector('.llm-ask-kicker').textContent).toBe('Bond Geometry changeApplied by the Copilot');
        expect(button('Undo')).toBeDefined();
        act(() => button('Open Bond Geometry').click());
        expect(onOpenPage).toHaveBeenCalledWith('geometry', null);
        // A change is not an analysis.
        expect(container.querySelector('.llm-turn-footer')).toBeNull();
        expect(container.querySelector('.llm-outcome').className).toContain('is-achieved');
    });

    it('Decline: the model reads {"declined": true}, and the card says nothing changed', async () => {
        stubReplies([callChange(), answered()]);
        act(() => root.render(chat({ tools: [changeTool()] })));
        await openChat();
        await ask('Set the window.');
        act(() => button('Decline').click());
        await settleUntil(() => container.querySelector('.llm-outcome'));
        expect(JSON.parse(chatBodies[1].messages.at(-1).content).data).toEqual({ declined: true });
        expect(card().className).toContain('is-declined');
        expect(card().textContent).toContain('Declined: nothing changed.');
        expect(button('Undo')).toBeUndefined();
    });

    it('Auto, from the drawer: no card to wait on, and the card says so; a reload is back at Ask first', async () => {
        stubReplies([callChange(), answered()]);
        act(() => root.render(chat({ tools: [changeTool()] })));
        await openChat();
        openSettings();
        expect(mode('Ask first').getAttribute('aria-checked')).toBe('true');
        expect(mode('Auto').getAttribute('aria-checked')).toBe('false');
        act(() => mode('Auto').click());
        expect(mode('Auto').getAttribute('aria-checked')).toBe('true');
        openSettings();
        expect(container.querySelector('.llm-composer-status').textContent).toMatch(/changes its inputs without asking \(Auto\)\.$/);
        await ask('Set the window.');
        await settleUntil(() => container.querySelector('.llm-outcome'));
        expect(chatBodies).toHaveLength(2);
        expect(card().className).toContain('is-applied');
        expect(card().querySelector('.llm-ask-kicker').textContent).toBe('Bond Geometry changeApplied by the Copilot · Auto');

        // A reload: a new dock, back at Ask first (Auto is never saved).
        act(() => root.unmount());
        root = createRoot(container);
        act(() => root.render(chat({ tools: [changeTool()] })));
        await settle();
        if (container.querySelector('.llm-chat-launcher')) await openChat();
        openSettings();
        expect(mode('Ask first').getAttribute('aria-checked')).toBe('true');
        // Nothing saved holds the mode.
        const saved = Array.from({ length: localStorage.length }, (_, index) => `${localStorage.key(index)}=${localStorage.getItem(localStorage.key(index))}`);
        expect(saved.length).toBeGreaterThan(0);
        expect(saved.join('\n')).not.toMatch(/approve|auto/i);
    });

    it('Stop while the card waits ends the turn: nothing applied, the card says stopped', async () => {
        stubReplies([callChange()]);
        const tool = changeTool();
        act(() => root.render(chat({ tools: [tool] })));
        await openChat();
        await ask('Set the window.');
        expect(card().className).toContain('is-waiting');
        act(() => container.querySelector('.llm-send').click());
        await settle();
        expect(card().className).toContain('is-stopped');
        expect(card().textContent).toContain('Stopped before the change was made.');
        expect(button('Apply')).toBeUndefined();
        expect(chatBodies).toHaveLength(1);
    });

    it('a change the page refused after the card (stale) says it was not applied', async () => {
        stubReplies([callChange(), answered()]);
        act(() => root.render(chat({ tools: [changeTool({ stale: true })] })));
        await openChat();
        await ask('Set the window.');
        act(() => button('Apply').click());
        await settleUntil(() => container.querySelector('.llm-outcome'));
        expect(card().className).toContain('is-stale');
        expect(card().textContent).toContain('Not applied: the page changed while this waited for you.');
        expect(button('Undo')).toBeUndefined();
    });

    it('Undo on the applied card puts the page back through the tool; a refused Undo says why and stays', async () => {
        stubReplies([callChange(), answered(), callChange(), answered()]);
        const undo = vi.fn()
            .mockResolvedValueOnce({ ok: false, reason: 'the inputs changed on the page since' })
            .mockResolvedValueOnce({ ok: true, computed: false });
        const tool = changeTool({ undo });
        act(() => root.render(chat({ tools: [tool] })));
        await openChat();
        await ask('Set the window.');
        act(() => button('Apply').click());
        await settleUntil(() => container.querySelector('.llm-outcome'));

        await act(async () => { button('Undo').click(); });
        await settle();
        expect(undo).toHaveBeenCalledWith(expect.objectContaining({ status: 'applied', id: 'w1', line: LINE }));
        expect(card().textContent).toContain('Not undone: the inputs changed on the page since.');
        expect(button('Undo')).toBeDefined();

        await act(async () => { button('Undo').click(); });
        await settle();
        expect(undo).toHaveBeenCalledTimes(2);
        expect(card().textContent).toContain('Undone: the earlier inputs are back.');
        expect(card().querySelector('.llm-ask-kicker').textContent).toBe('Bond Geometry changeUndone');
        expect(button('Undo')).toBeUndefined();
    });

    it('says in one short line why an applied change did not compute, and why a change was refused', async () => {
        stubReplies([callChange(), answered(), callChange(), answered()]);
        const notComputed = changeTool();
        notComputed.run = vi.fn(async (args, { approve, signal }) => {
            const { auto } = await approve({ tool: 'set_window', title: 'Change the window', line: LINE, details: DETAILS }, { signal });
            return {
                data: { applied: true, computed: false, not_computed: 'the angle budget is exceeded' },
                checks: [{ level: 'warn', text: `Applied: ${LINE}, but the page did not compute it: the angle budget is exceeded.` }],
                change: { status: 'applied', line: LINE, details: DETAILS, auto, id: 'w2', note: 'The page did not compute it: the angle budget is exceeded.' }
            };
        });
        act(() => root.render(chat({ tools: [notComputed] })));
        await openChat();
        await ask('Set the window.');
        act(() => button('Apply').click());
        await settleUntil(() => container.querySelector('.llm-outcome'));
        const lines = () => [...card().querySelectorAll('.llm-step-headline')].map((item) => item.textContent);
        expect(lines()).toEqual(['The page did not compute it: the angle budget is exceeded.']);

        // A change the tool refused before any card: its reason, capitalised.
        const refusing = changeTool();
        refusing.run = vi.fn(async () => { throw new Error('no run is open: open a run folder, or load the Demo'); });
        act(() => root.render(chat({ tools: [refusing], runKey: 'run-2' })));
        await settle();
        await ask('Set the window.');
        await settleUntil(() => container.querySelector('.llm-outcome'));
        expect(card().className).toContain('is-not_applied');
        expect(lines()).toEqual(['No run is open: open a run folder, or load the Demo']);
    });

    it('marks the rail when a change waits while the Copilot is collapsed', async () => {
        stubReplies([callChange(), answered()]);
        act(() => root.render(chat({ tools: [changeTool()] })));
        await openChat();
        act(() => setTextarea(container.querySelector('.llm-composer-input'), 'Set the window.'));
        act(() => container.querySelector('.llm-composer').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true })));
        act(() => container.querySelector('.llm-chat-panel').dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true })));
        await settle();
        const rail = container.querySelector('.llm-dock-rail');
        expect(rail.className).toContain('has-unread');
        expect(rail.getAttribute('aria-label')).toBe('AI Copilot — a change waits for your approval');
        await openChat();
        act(() => button('Apply').click());
        await settleUntil(() => container.querySelector('.llm-outcome'));
        expect(card().className).toContain('is-applied');
    });
});
