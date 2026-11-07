// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The context meter in the dock (THY-76), against a scripted LM Studio
// server whose loaded context the test sets: the meter shows the guard's
// numbers (the requests as sent, counted by estimatePromptTokens, against
// the loaded window with LM Studio's 4,096 tokens of reply), turns amber and
// red, offers Start fresh and Summarize and continue, and shows a refused
// request. The chat says when the history drops a question or shortens an
// answer. No live model is involved.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import React, { act } from 'react';
import { createRoot } from 'react-dom/client';
import CopilotChat from '../components/CopilotChat';
import { DOCK_STORAGE_KEY } from '../components/dockLayout';
import { LMSTUDIO_REPLY_MIN_TOKENS, estimatePromptTokens, forgetServers } from '../provider/contextGuard';
import { saveSettings } from '../settings';
import { SUMMARY_ANSWER_CHARS } from '../useCopilot';

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
const answer = (text) => sse([{ choices: [{ delta: { content: text }, finish_reason: 'stop' }] }]);
const json = (body, status = 200) => ({ ok: status < 400, status, json: async () => body });

const settle = () => act(async () => {
    for (let i = 0; i < 10; i += 1) await new Promise((resolve) => setTimeout(resolve, 0));
});

const setTextarea = (textarea, value) => {
    Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype, 'value').set.call(textarea, value);
    textarea.dispatchEvent(new Event('input', { bubbles: true }));
};

// A long answer: sentences with numbers, beyond the turn's 2,000-character memory.
const LONG = `${Array.from({ length: 80 }, (_, index) => `Shell ${index} holds ${index + 2}.1 Se at 2.${index} Å.`).join(' ')}\n\nOutcome: achieved — counted.`;

describe('the context meter in the dock', () => {
    let container;
    let root;
    let chatBodies;
    let loaded;
    let replies;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        forgetServers();
        localStorage.removeItem(DOCK_STORAGE_KEY);
        saveSettings({ baseUrl: 'http://localhost:1234/v1', model: 'tool-model', apiKey: '', copilotTools: true });
        chatBodies = [];
        loaded = 1_000_000;
        replies = [];
        vi.stubGlobal('fetch', vi.fn(async (url, options = {}) => {
            const path = new URL(String(url)).pathname;
            if (path === '/api/v0/models') return json({ data: [{ id: 'tool-model', state: 'loaded', loaded_context_length: loaded }] });
            if (path.startsWith('/api/')) return json({ error: 'not found' }, 404);
            if (path.endsWith('/v1/models')) return json({ data: [{ id: 'tool-model' }] });
            chatBodies.push(JSON.parse(options.body));
            return replies.shift() || answer('Fine.\n\nOutcome: achieved — x.');
        }));
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.unstubAllGlobals();
    });

    const open = async () => {
        act(() => root.render(
            <CopilotChat runName="demo" plotFiles={[{ plotKind: 'neutron_sq', plotData: { title: 'S(Q)', metrics: { rwp: 0.04 } } }]} runKey="run-1" />
        ));
        await settle();
        act(() => container.querySelector('.llm-chat-launcher').click());
        await settle();
    };
    const ask = async (text) => {
        act(() => setTextarea(container.querySelector('.llm-composer-input'), text));
        act(() => container.querySelector('.llm-composer').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true })));
        await settle();
        await settle();
    };
    const meter = () => container.querySelector('.llm-meter');
    const details = () => meter().querySelector('.llm-meter-details').textContent;
    const button = (name) => [...container.querySelectorAll('button')].find((element) => element.textContent === name) || null;
    // What the meter says a request uses (its details' first line).
    const shownUsed = () => Number(details().match(/^([\d,]+) of/)[1].replace(/,/g, ''));
    // The guard's count of a request as it went out, with LM Studio's reply room.
    const usedBy = (body) => estimatePromptTokens({ messages: body.messages, tools: body.tools }) + LMSTUDIO_REPLY_MIN_TOKENS;

    it('shows the next request against the loaded window, in the guard\'s numbers', async () => {
        await open();
        expect(meter().dataset.level).toBe('neutral');
        expect(meter().textContent).toMatch(/^Context 1%/);
        expect(details()).toMatch(/of 1,000,000 tokens: the prompt [\d,]+ \+ the reply reserve 4,096\./);
        expect(details()).toContain('Window: the context LM Studio loaded the model with.');
        expect(details()).toContain('Estimated high');
        // Local: no cost line.
        expect(details()).not.toContain('Cost');
        const before = shownUsed();

        await ask('Has it converged?');
        // The request as sent is what the meter previewed (the question aside).
        const question = Math.ceil('Has it converged?'.length / 2.5);
        expect(Math.abs(usedBy(chatBodies[0]) - before)).toBeLessThanOrEqual(question + 1);
        expect(details()).toMatch(/Next question starts at [\d,]+; the last answer peaked at [\d,]+\./);
        expect(shownUsed()).toBeGreaterThanOrEqual(usedBy(chatBodies[0]));
    });

    it('turns amber at 70 %, red above 90 %, and offers Start fresh and Summarize and continue', async () => {
        replies = [answer(LONG)];
        await open();
        await ask('Which shells?');
        const used = shownUsed();

        // The window the next request fills to 80 %: amber (re-read after an answer).
        loaded = Math.round(used / 0.8);
        await ask('And the next?');
        expect(meter().dataset.level).toBe('amber');
        expect(meter().textContent).toMatch(/Context (7\d|8\d|90)%/);
        expect(button('Start fresh')).not.toBeNull();
        // The long first answer can go shorter: Summarize and continue is offered.
        expect(button('Summarize and continue')).not.toBeNull();
        // At amber the meter's Start fresh replaces New conversation.
        expect(button('New conversation')).toBeNull();

        act(() => button('Summarize and continue').click());
        await settle();
        expect(container.querySelector('.llm-chat-log').textContent).toContain('Summarized — the answers above go to the model as their sentences with numbers');
        await ask('Third?');
        const memory = chatBodies.at(-1).messages.slice(3, -1);
        expect(memory[1].role).toBe('assistant');
        expect(memory[1].content.length).toBeLessThanOrEqual(SUMMARY_ANSWER_CHARS);
        expect(memory[1].content).toContain('Outcome: achieved — counted.');

        // Red: the next request at ~95 % of the window.
        loaded = Math.round(shownUsed() / 0.95);
        await ask('Fourth?');
        expect(meter().dataset.level).toBe('red');

        // Start fresh leaves the prompt alone: on a window this small it stays
        // red, and the details say what makes room.
        const full = shownUsed();
        expect(details()).toContain('Reload the model in LM Studio with a longer context');
        act(() => button('Start fresh').click());
        await settle();
        expect(container.querySelectorAll('.llm-message')).toHaveLength(0);
        expect(shownUsed()).toBeLessThan(full);
    });

    it('shows a request the guard refused, over 100 %, and nothing was sent', async () => {
        await open();
        await ask('First?');
        const sent = chatBodies.length;
        loaded = Math.round(shownUsed() / 1.25);
        await ask('Second?');
        expect(chatBodies).toHaveLength(sent);
        expect(meter().dataset.level).toBe('red');
        expect(meter().textContent).toMatch(/Context 1[23]\d%/);
        expect(details()).toMatch(/^Refused: this request needs about [\d,]+ tokens; the model holds [\d,]+\. Nothing was sent\./);
        expect(container.querySelector('[role="alert"]').textContent).toMatch(/LM Studio runs tool-model with a [\d,]+-token context/);
    });

    it('says in the chat when the oldest question stops being sent and when an answer goes shortened', async () => {
        replies = [answer(LONG)];
        await open();
        await ask('q1');
        await ask('q2');
        const notes = () => [...container.querySelectorAll('.llm-memory-note')].map((element) => element.textContent);
        expect(notes()).toEqual(['Memory: the last answer is sent shortened to its sentences with numbers.']);
        await ask('q3');
        await ask('q4');
        await ask('q5');
        expect(notes()).toHaveLength(1);
        await ask('q6');
        expect(notes().at(-1)).toBe('Memory: the earliest question is no longer sent to the model (it keeps the last 4).');
        // What was sent agrees: q1 is not in the last request.
        expect(chatBodies.at(-1).messages.some((message) => message.content === 'q1')).toBe(false);
        expect(chatBodies.at(-1).messages.some((message) => message.content === 'q2')).toBe(true);
    });
});
