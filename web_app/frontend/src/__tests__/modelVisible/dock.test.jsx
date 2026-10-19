// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// What the model sees, 2-3 of 4, the app's path: the dock (CopilotChat →
// CopilotView → useCopilot → the turn module) on the Demo run, props as
// App.jsx passes them, a question typed and sent. The first request is
// captured where the dock's model port (providerModel, provider/model.js)
// sends it, streamCopilotChat (provider/contextGuard.js), replaced here by a
// recording fake; everything before it is the app's own code. No model runs
// and no server is contacted (the connection probe gets a stubbed model
// list). The file's header line still says "as useCopilot.js calls it", its
// wording before the turn module, so the file stays byte-identical.
//
//   dock/first-request.txt   fields, tool names, messages
//   dock/tools.json          the `tools` field exactly as sent
//
// The run props and the host's tools come from the run module's core
// (copilotRun.js) as App.jsx gets them through useCopilotRun, except the
// engine: the browser adapter runs the engines in a Web Worker (none in
// jsdom), so the run gets the same handler in-thread (pcaKdeWorker's
// handlePcaMessage); the tool definitions do not depend on it.
// useCopilotRun.test.jsx checks that the hook hands App what the core gives.
//
// A change to these files is a change to what the model sees: call it out in
// the PR.

import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import React, { act } from 'react';
import { createRoot } from 'react-dom/client';
import CopilotChat from '../../llm/components/CopilotChat';
import { DOCK_STORAGE_KEY } from '../../llm/components/dockLayout';
import { DEFAULT_SETTINGS, saveSettings } from '../../llm/settings';
import { skillIndex } from '../../llm/copilot/skills';
import { copilotSkills } from '../../copilotSkills';
import { buildCopilotCards, COPILOT_PAGES } from '../../copilotCards';
import { buildViewTools } from '../../copilotViewTools';
import {
    bundledSkills, jsonFile, loadDemoAssistantRun, promptVariants, recordingModel, renderRequest
} from './helpers';

const model = vi.hoisted(() => ({ current: null }));
vi.mock('../../llm/provider/contextGuard', async (importOriginal) => ({
    ...(await importOriginal()),
    streamCopilotChat: (request) => model.current.callModel(request)
}));

const QUESTION = 'Has this run converged?';
const MODEL_ID = 'snapshot-model';

const settle = () => act(async () => {
    for (let i = 0; i < 10; i += 1) await new Promise((resolve) => setTimeout(resolve, 0));
});

describe('the dock\'s first request (the Demo run)', () => {
    let demo;
    let container;
    let root;

    beforeAll(async () => {
        demo = await loadDemoAssistantRun();
    });

    afterEach(() => {
        act(() => root?.unmount());
        container?.remove();
        vi.unstubAllGlobals();
    });

    it('is the app\'s own request', async () => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        localStorage.removeItem(DOCK_STORAGE_KEY);
        saveSettings({ ...DEFAULT_SETTINGS, model: '' });
        // The connection probe: an OpenAI-compatible server with one model.
        vi.stubGlobal('fetch', vi.fn(async (url) => (String(url).endsWith('/models')
            ? { ok: true, status: 200, json: async () => ({ data: [{ id: MODEL_ID }] }) }
            : { ok: false, status: 404, json: async () => ({ error: 'not found' }) })));
        model.current = recordingModel();

        const { run, copilot } = demo;
        const cards = buildCopilotCards({ plotFiles: run.plotFiles, rValueFile: run.rValueFile, symmetry: run.symmetry });
        const tools = [
            ...copilot.tools,
            ...buildViewTools({ cards, pages: COPILOT_PAGES, activePage: () => 'dashboard' })
        ];

        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
        act(() => root.render(
            <CopilotChat
                className="app-copilot"
                {...run}
                tools={tools}
                pages={COPILOT_PAGES}
                onOpenPage={() => {}}
                runKey="Demo|1"
                cards={cards}
                onMarks={() => {}}
                skills={copilotSkills}
            />
        ));
        act(() => container.querySelector('.llm-chat-launcher').click());
        await settle();
        const textarea = container.querySelector('.llm-composer-input');
        act(() => {
            Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype, 'value').set.call(textarea, QUESTION);
            textarea.dispatchEvent(new Event('input', { bubbles: true }));
        });
        act(() => container.querySelector('.llm-composer').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true })));
        for (let i = 0; i < 20 && !model.current.requests.length; i += 1) await settle();

        const [first] = model.current.requests;
        expect(first).toBeDefined();
        expect(first.model).toBe(MODEL_ID);
        await expect(jsonFile(first.tools)).toMatchFileSnapshot('__snapshots__/dock/tools.json');
        await expect(renderRequest(first, {
            title: 'First request of the dock (CopilotChat + useCopilot), the Demo run (GaTa4Se8 250 K)',
            port: 'streamCopilotChat(request) (llm/provider/contextGuard.js), as useCopilot.js calls it',
            variants: promptVariants(skillIndex(bundledSkills())),
            toolsFile: 'dock/tools.json'
        })).toMatchFileSnapshot('__snapshots__/dock/first-request.txt');
    });
});
