// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The AI Copilot dock's state lifecycle, App end to end (plan 1.4; the smoke
// test the manager asked for instead of a manual browser check): the real
// App, Dashboard, Detected SG card, run module and dock, on the Demo run's
// files, with no model and no server.
//
//   1. the Demo's files picked as a folder (static build, File System Access,
//      Live Data on), the Copilot opened, a question asked and answered;
//   2. RMCProfile saves a new configuration into the folder: the Live Data
//      poll sees the .rmc6f signature change (configurationSignature, App's
//      configKey);
//   3. the dock draws the divider, and the next request carries no earlier
//      question or answer: its system prompt and context are a fresh
//      conversation's, for the new configuration; memory resumes after it;
//   4. a Detected SG pick is kept across another save of the same run (card
//      and Copilot context), and back at the card's default
//      (DEFAULT_VIEW.toleranceA) when another run (the Demo button) opens.
//
// What is not the app's own code: the environment (static mode with the File
// System Access API, the folder handle over public/demo, an in-thread stand-in
// for the Dashboard's structure worker, which jsdom lacks), and the model: a
// recording fake at streamCopilotChat (provider/contextGuard.js), where the
// dock's model port sends each request (as modelVisible/dock.test.jsx). fetch
// serves the Demo button's files from disk and gives the connection probe one
// model; axios refuses every call. Nothing is fetched from a network.

import { afterEach, beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const model = vi.hoisted(() => ({ current: null }));
const seen = vi.hoisted(() => ({ chat: null, axios: [], fetch: [] }));

// No server: static mode makes no axios call, and any call would fail here.
vi.mock('axios', () => {
    const refuse = (method) => vi.fn(async (url) => {
        seen.axios.push(`${method} ${url}`);
        throw new Error(`no server in this test (${method} ${url})`);
    });
    return { default: { get: refuse('GET'), post: refuse('POST'), isCancel: () => false } };
});

// The static build in a Chromium browser: a picked folder, watched by Live Data.
vi.mock('../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => true,
    supportsFileSystemAccess: () => true,
}));

// The model: each request recorded where the dock's model port sends it.
vi.mock('../llm/provider/contextGuard', async (importOriginal) => ({
    ...(await importOriginal()),
    streamCopilotChat: (request) => model.current.callModel(request),
}));

// The dock as App renders it; its props are read on the way in.
vi.mock('../llm', async (importOriginal) => {
    const real = await importOriginal();
    const { createElement } = await import('react');
    return {
        ...real,
        CopilotChat: function RecordedCopilotChat(props) {
            seen.chat = props;
            return createElement(real.CopilotChat, props);
        },
    };
});

const { default: App } = await import('../App');
const { WATCH_INTERVAL_MS } = await import('../browserData');
const { COPILOT_RUN_PROPS, DEFAULT_VIEW } = await import('../copilotRun');
const { DOCK_STORAGE_KEY } = await import('../llm/components/dockLayout');
const { DEFAULT_SETTINGS, saveSettings } = await import('../llm/settings');
const { forgetServers } = await import('../llm/provider/contextGuard');
const { buildRunContext } = await import('../llm/context/runContext');
const { contextBlock } = await import('../llm/prompts/copilot');
const { skillIndex } = await import('../llm/copilot/skills');
const { handleStructureMessage } = await import('../workers/localStructureWorker');
const { bundledSkills, DEMO_TIME, demoEntries, promptVariants, recordingModel } = await import('./modelVisible/helpers');

const HERE = dirname(fileURLToPath(import.meta.url));
const DEMO = join(HERE, '../../public/demo');
const RMC6F = 'GTS_250K.rmc6f';
const FOLDER = 'GTS_250K_run';
const MODEL_ID = 'lifecycle-model';
const DIVIDER = 'New configuration saved — earlier answers were about the previous one';
const MOVES_SAVED = 1560869;

// jsdom has no Web Worker: the Dashboard's structure worker runs its handler
// in-thread, answering with the worker's own messages (localStructureWorker.js).
class InThreadWorker {
    constructor(url) {
        this.url = String(url);
        this.onmessage = null;
        this.onerror = null;
        this.terminated = false;
    }

    postMessage(data) {
        if (!this.url.includes('localStructureWorker')) throw new Error(`no in-thread worker for ${this.url}`);
        const reply = (message) => {
            if (!this.terminated) this.onmessage?.({ data: { id: data?.id, ...message } });
        };
        handleStructureMessage(data).then(
            (result) => reply({ result }),
            (error) => reply({ error: error.message || 'Browser structure parser failed' })
        );
    }

    terminate() {
        this.terminated = true;
    }
}

// The run folder on "disk": the Demo's files, by name. RMCProfile's save
// writes a new .rmc6f (a later move count, one more prior save, a later time).
const folder = new Map();
let rmc6fText = '';
let saves = 0;
const saveConfiguration = (movesGenerated) => {
    saves += 1;
    const text = rmc6fText
        .replace(/(Number of moves generated:\s*)\d+/, `$1${movesGenerated}`)
        .replace(/(Number of prior configuration saves:\s*)\d+/, `$1${saves}`);
    folder.set(RMC6F, new File([text], RMC6F, { lastModified: DEMO_TIME + saves * 600_000 }));
};
// A FileSystemDirectoryHandle over the folder: getFile() reads it as it is now.
const folderHandle = () => ({
    kind: 'directory',
    name: FOLDER,
    async *entries() {
        for (const [name, file] of folder) yield [name, { kind: 'file', getFile: async () => file }];
    },
});

const notFound = { ok: false, status: 404, json: async () => ({ error: 'not found' }) };
const fetchStub = async (url) => {
    const href = String(url);
    seen.fetch.push(href);
    // The Demo button's files (browserData.loadDemoRun), from public/demo.
    if (href.includes('/demo/')) {
        const bytes = readFileSync(join(DEMO, decodeURIComponent(href.split('/').pop())));
        return { ok: true, status: 200, blob: async () => new Blob([bytes]) };
    }
    // The connection probe: an OpenAI-compatible server with one model; its
    // Ollama / LM Studio probes (contextGuard.js) find nothing.
    if (!/\/api\//.test(href) && href.endsWith('/v1/models')) {
        return { ok: true, status: 200, json: async () => ({ data: [{ id: MODEL_ID }] }) };
    }
    return notFound;
};

// The fake model answers each question by its number, so memory can be read back.
const ANSWERS = [
    'Answer one: the χ² is still falling.\n\nOutcome: achieved — read from the χ² history.',
    'Answer two: about three accepted moves per atom.\n\nOutcome: achieved — read from the move counters.',
    'Answer three: about twelve per cent.\n\nOutcome: achieved — read from the move counters.',
    'Answer four: more saves are needed.\n\nOutcome: achieved — read from the move counters.',
    'Answer five: the χ² is still falling.\n\nOutcome: achieved — read from the χ² history.',
];
const QUESTIONS = [
    'Has this run converged?',
    'How far has the sampling come?',
    'And the acceptance ratio?',
    'Is the sampling enough now?',
    'Has the Demo run converged?',
];

const settle = () => act(async () => {
    for (let i = 0; i < 10; i += 1) await new Promise((resolve) => setTimeout(resolve, 0));
});
// Settle until `ready()` holds; a structure parse or the symmetry finder can
// take a few hundred milliseconds of wall time.
const waitFor = async (ready, what, timeoutMs = 60_000) => {
    const start = Date.now();
    while (!ready()) {
        if (Date.now() - start > timeoutMs) throw new Error(`timed out waiting for ${what}`);
        await settle();
    }
};

// The eight run props App hands the dock, as they are now.
const runProps = () => Object.fromEntries(COPILOT_RUN_PROPS.map((key) => [key, seen.chat[key]]));
// The context JSON in a request's context message.
const contextOf = (request) => JSON.parse(request.messages[1].content.match(/```json\n([\s\S]*)\n```/)[1]);
const roundSig = (value) => Number(value.toPrecision(2));

describe('the Copilot dock\'s state lifecycle in App (no model, no server)', () => {
    let container;
    let root;

    beforeAll(() => {
        rmc6fText = readFileSync(join(DEMO, RMC6F), 'utf8');
    });

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        // Only the Live Data poll's clock is faked: the dock's own timers run.
        vi.useFakeTimers({ toFake: ['setInterval', 'clearInterval'] });
        vi.stubGlobal('fetch', vi.fn(fetchStub));
        vi.stubGlobal('Worker', InThreadWorker);
        forgetServers();
        localStorage.removeItem(DOCK_STORAGE_KEY);
        saveSettings({ ...DEFAULT_SETTINGS, model: '' });
        model.current = recordingModel(ANSWERS.map((text) => ({ text })));
        seen.chat = null;
        seen.axios = [];
        seen.fetch = [];
        folder.clear();
        saves = 0;
        for (const { file } of demoEntries()) folder.set(file.name, file);
        window.showDirectoryPicker = async () => folderHandle();
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        delete window.showDirectoryPicker;
        vi.unstubAllGlobals();
        vi.useRealTimers();
        // Hermetic: no server was asked, and fetch only served the Demo's files
        // and the stubbed model server.
        expect(seen.axios).toEqual([]);
        expect(seen.fetch.every((href) => href.includes('/demo/') || href.startsWith(DEFAULT_SETTINGS.baseUrl.replace(/\/v1$/, '')))).toBe(true);
    });

    const click = async (element) => {
        expect(element).toBeTruthy();
        await act(async () => { element.click(); });
    };
    const button = (label) => [...container.querySelectorAll('button')]
        .find((element) => element.textContent.trim() === label);
    const liveData = () => container.querySelector('label.ui-switch-outline input[type="checkbox"]');
    const poll = async () => {
        await act(async () => { vi.advanceTimersByTime(WATCH_INTERVAL_MS); });
        await settle();
    };
    const dividers = () => [...container.querySelectorAll('.llm-config-divider')];
    // The Detected SG card's bricks, the selected one, and one by its group.
    const activeBrick = () => container.querySelector('.sym-brick.is-active .sym-brick-label')?.textContent ?? null;
    const brick = (group) => [...container.querySelectorAll('.sym-brick')]
        .find((element) => element.querySelector('.sym-brick-label')?.textContent === group);
    // The run as the Copilot is handed it is complete: parsed plots, the χ²
    // history, the .dat, and the symmetry of the structure the Dashboard holds.
    const runReady = (movesGenerated) => {
        const chat = seen.chat;
        return Boolean(chat?.symmetry && chat.runSettings && chat.rValueFile
            && chat.plotFiles.length && chat.plotFiles.every((file) => file.plotData || file.parseError)
            && chat.structure?.moves?.generated === movesGenerated);
    };
    // Ask through the composer; the request it sends and the props it was asked with.
    const ask = async (index) => {
        const sent = model.current.requests.length;
        const props = runProps();
        const textarea = container.querySelector('.llm-composer-input');
        expect(textarea.disabled).toBe(false);
        act(() => {
            Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype, 'value').set.call(textarea, QUESTIONS[index]);
            textarea.dispatchEvent(new Event('input', { bubbles: true }));
        });
        act(() => container.querySelector('.llm-composer').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true })));
        await waitFor(
            () => model.current.requests.length > sent && !container.querySelector('.llm-send.is-stop')
                && container.querySelector('.llm-chat-log')?.textContent.includes(ANSWERS[index].split('\n')[0]),
            `the answer to question ${index + 1}`
        );
        // One request per question: the fake answer is never sent back.
        expect(model.current.requests).toHaveLength(sent + 1);
        return { request: model.current.requests[sent], props };
    };
    // A request's messages after its system prompt hold none of these
    // questions or answers (memory carries each verbatim).
    const expectNoMemoryOf = (request, indices) => {
        const conversation = request.messages.slice(1);
        for (const index of indices) {
            expect(conversation.some(({ content }) => content === QUESTIONS[index])).toBe(false);
            expect(conversation.some(({ content }) => content.includes(ANSWERS[index].split('\n')[0]))).toBe(false);
        }
    };

    it('ends the model\'s memory at a new configuration, keeps the Detected SG pick across a save, and resets it with another run', async () => {
        const systemPrompt = promptVariants(skillIndex(bundledSkills())).get('tools-view-ask-skills');

        // --- 1. the Demo's files picked as a folder, Live Data on, the Copilot open ---
        await act(async () => { root.render(<App />); });
        await click(button('Select Folder'));
        await click(liveData());
        expect(liveData().checked).toBe(true);
        await click(container.querySelector('.llm-chat-launcher'));
        await waitFor(() => runReady(MOVES_SAVED) && container.querySelector('.llm-composer-input:not(:disabled)'),
            'the run published to the open Copilot');
        const { runKey, configKey } = seen.chat;
        expect(seen.chat.runName).toBe(FOLDER);
        expect(configKey).toMatch(/^local\|/);
        // The card's default, on the card and in the Copilot's context.
        expect(activeBrick()).toBe('F-43m');
        expect(seen.chat.symmetry.toleranceA).toBe(DEFAULT_VIEW.toleranceA);

        const first = await ask(0);
        // A fresh conversation's first request: the prompt, the context, the
        // acknowledgement, the question.
        expect(first.request.model).toBe(MODEL_ID);
        expect(first.request.messages.map((message) => message.role)).toEqual(['system', 'user', 'assistant', 'user']);
        expect(first.request.messages[0].content).toBe(systemPrompt);
        expect(first.request.messages[1].content).toBe(contextBlock(buildRunContext(first.props)));
        expect(first.request.messages[3]).toEqual({ role: 'user', content: QUESTIONS[0] });
        expect(contextOf(first.request).live_mode).toBe(true);
        expect(contextOf(first.request).symmetry.tolerance_A).toBe(DEFAULT_VIEW.toleranceA);
        expect(dividers()).toHaveLength(0);

        // --- 2. RMCProfile saves a new configuration; the Live Data poll sees it ---
        saveConfiguration(1_600_000);
        await poll();
        await waitFor(() => runReady(1_600_000), 'the new configuration published to the Copilot');
        expect(seen.chat.runKey).toBe(runKey);
        expect(seen.chat.configKey).not.toBe(configKey);

        // --- 3. the divider, then a request with no memory of the earlier answer ---
        expect(dividers().map((element) => element.textContent)).toEqual([DIVIDER]);
        expect(container.querySelector('.llm-chat-log').textContent).toContain(ANSWERS[0].split('\n')[0]);
        const second = await ask(1);
        expect(second.request.messages.map((message) => message.role)).toEqual(['system', 'user', 'assistant', 'user']);
        expectNoMemoryOf(second.request, [0]);
        // As for a fresh conversation: the same prompt, tools and acknowledgement,
        // and the context of the run as App hands it now, which is the first
        // context with the new configuration's move count, and only that.
        expect(second.request.messages[0].content).toBe(systemPrompt);
        expect(second.request.messages[1].content).toBe(contextBlock(buildRunContext(second.props)));
        expect(second.request.messages[1].content).toBe(first.request.messages[1].content
            .replace(`"moves_generated": ${MOVES_SAVED}`, '"moves_generated": 1600000'));
        expect(second.request.messages[2]).toEqual(first.request.messages[2]);
        expect(second.request.messages[3]).toEqual({ role: 'user', content: QUESTIONS[1] });
        expect(second.request.tools).toEqual(first.request.tools);
        // The divider stays between the two configurations' exchanges.
        expect(dividers()).toHaveLength(1);

        // Memory resumes within the new configuration, and only from the divider.
        const third = await ask(2);
        expect(third.request.messages.slice(3)).toEqual([
            { role: 'user', content: QUESTIONS[1] },
            { role: 'assistant', content: ANSWERS[1] },
            { role: 'user', content: QUESTIONS[2] },
        ]);
        expectNoMemoryOf(third.request, [0]);

        // --- 4a. a Detected SG pick is kept across a Live Data save of the same run ---
        await click(brick('P-42_1m'));
        await waitFor(() => seen.chat.symmetry?.spaceGroup === 'P-42_1m', 'the pick in the Copilot\'s run');
        const picked = seen.chat.symmetry.toleranceA;
        expect(picked).not.toBe(DEFAULT_VIEW.toleranceA);
        expect(activeBrick()).toBe('P-42_1m');

        saveConfiguration(1_700_000);
        await poll();
        await waitFor(() => runReady(1_700_000), 'the second save published to the Copilot');
        expect(seen.chat.runKey).toBe(runKey);
        expect(activeBrick()).toBe('P-42_1m');
        expect(seen.chat.symmetry.toleranceA).toBe(picked);
        expect(dividers()).toHaveLength(2);

        const fourth = await ask(3);
        expect(fourth.request.messages.map((message) => message.role)).toEqual(['system', 'user', 'assistant', 'user']);
        expectNoMemoryOf(fourth.request, [0, 1, 2]);
        expect(fourth.request.messages[1].content).toBe(contextBlock(buildRunContext(fourth.props)));
        expect(contextOf(fourth.request).symmetry.tolerance_A).toBe(roundSig(picked));
        expect(contextOf(fourth.request).symmetry.space_group).toMatch(/^P-42_1m /);
        expect(contextOf(fourth.request).configuration_optimization.moves_generated).toBe(1_700_000);

        // --- 4b. another run (the Demo button): the pick is back at the card's default ---
        await click(container.querySelector('button.demo-button'));
        await waitFor(() => seen.chat.runName === 'Demo' && runReady(MOVES_SAVED), 'the Demo run published to the Copilot');
        expect(seen.chat.runKey).not.toBe(runKey);
        expect(seen.chat.symmetry.toleranceA).toBe(DEFAULT_VIEW.toleranceA);
        expect(activeBrick()).toBe('F-43m');
        // A new conversation: nothing of the last run's on screen.
        expect(container.querySelector('.llm-chat-log')).toBeNull();
        expect(dividers()).toHaveLength(0);

        const fifth = await ask(4);
        expect(fifth.request.messages.map((message) => message.role)).toEqual(['system', 'user', 'assistant', 'user']);
        expect(fifth.request.messages[0].content).toBe(systemPrompt);
        expect(fifth.request.messages[1].content).toBe(contextBlock(buildRunContext(fifth.props)));
        expect(contextOf(fifth.request).run).toBe('Demo');
        expect(contextOf(fifth.request).symmetry.tolerance_A).toBe(DEFAULT_VIEW.toleranceA);
        expect(contextOf(fifth.request).symmetry.space_group).toMatch(/^F-43m /);
    }, 180_000);

    // The open run as the header, the dock and the Copilot's run show it: its
    // name, identity, configuration, move count and the folder of its plots.
    const openRun = () => ({
        header: container.querySelector('.local-file-bar .ui-fieldbar__value')?.textContent,
        runName: seen.chat.runName,
        runKey: seen.chat.runKey,
        configKey: seen.chat.configKey,
        moves: seen.chat.structure?.moves?.generated,
        folders: [...new Set(seen.chat.plotFiles.map((file) => file.path.split('/')[0]))],
    });

    it('stops watching the picked folder when the Demo opens, and watches the folder picked next', async () => {
        model.current = recordingModel([{ text: ANSWERS[4] }]);
        await act(async () => { root.render(<App />); });
        await click(button('Select Folder'));
        await click(liveData());
        await click(container.querySelector('.llm-chat-launcher'));
        await waitFor(() => seen.chat?.runName === FOLDER && runReady(MOVES_SAVED), 'the picked folder published to the Copilot');

        // The Demo opens with Live Data still on, and gets a conversation.
        await click(container.querySelector('button.demo-button'));
        await waitFor(() => seen.chat.runName === 'Demo' && runReady(MOVES_SAVED)
            && container.querySelector('.llm-composer-input:not(:disabled)'), 'the Demo run published to the Copilot');
        const demo = openRun();
        expect(demo).toMatchObject({ header: 'Demo', runName: 'Demo', moves: MOVES_SAVED, folders: ['Demo'] });
        await ask(4);

        // RMCProfile saves into the folder that was open before: the Demo stays.
        saveConfiguration(1_600_000);
        await poll();
        await poll();
        for (let i = 0; i < 5; i += 1) await settle();
        expect(openRun()).toEqual(demo);
        expect(dividers()).toHaveLength(0);
        expect(container.querySelector('button.demo-button').getAttribute('aria-pressed')).toBe('true');
        expect(liveData().checked).toBe(true);
        expect(container.querySelector('.ui-status.is-error')).toBeNull();

        // Leaving the Demo for the picked folder: another run, read as it is
        // now, and watched again.
        await click(button('Select Folder'));
        await waitFor(() => seen.chat.runName === FOLDER && runReady(1_600_000), 'the folder picked again published to the Copilot');
        const picked = openRun();
        expect(picked.runKey).not.toBe(demo.runKey);
        expect(picked).toMatchObject({ header: FOLDER, folders: [FOLDER] });
        saveConfiguration(1_700_000);
        await poll();
        await waitFor(() => runReady(1_700_000), 'the save in the folder picked again', 10_000);
        expect(openRun()).toMatchObject({ header: FOLDER, runKey: picked.runKey, folders: [FOLDER] });
        expect(seen.chat.configKey).not.toBe(picked.configKey);
    }, 120_000);

    it('watches the folder picked after the Demo when Live Data was switched on with the Demo open', async () => {
        await act(async () => { root.render(<App />); });
        await click(container.querySelector('.llm-chat-launcher'));
        await click(container.querySelector('button.demo-button'));
        await waitFor(() => seen.chat?.runName === 'Demo' && runReady(MOVES_SAVED), 'the Demo run published to the Copilot');
        await click(liveData());
        expect(liveData().checked).toBe(true);

        await click(button('Select Folder'));
        await waitFor(() => seen.chat.runName === FOLDER && runReady(MOVES_SAVED), 'the picked folder published to the Copilot');
        const picked = openRun();
        saveConfiguration(1_600_000);
        await poll();
        await waitFor(() => runReady(1_600_000), 'the save in the picked folder', 10_000);
        expect(openRun()).toMatchObject({ header: FOLDER, runKey: picked.runKey, folders: [FOLDER] });
        expect(seen.chat.configKey).not.toBe(picked.configKey);
    }, 120_000);
});
