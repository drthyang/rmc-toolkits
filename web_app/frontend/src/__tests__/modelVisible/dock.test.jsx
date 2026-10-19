// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// What the model sees, 2-3 of 4, the app's path: the dock (CopilotChat →
// CopilotView → useCopilot) on the Demo run, props as App.jsx passes them, a
// question typed and sent. The first request is captured at the model port
// useCopilot calls, streamCopilotChat (provider/contextGuard.js), replaced
// here by a recording fake; everything before it is the app's own code. No
// model runs and no server is contacted (the connection probe gets a stubbed
// model list).
//
//   dock/first-request.txt   fields, tool names, messages
//   dock/tools.json          the `tools` field exactly as sent
//
// The host's tools are built as App.jsx builds them, except the engine
// request: useCopilotTools runs the engines in a Web Worker (none in jsdom),
// so buildCopilotTools gets the same handler in-thread (pcaKdeWorker's
// handlePcaMessage); the tool definitions do not depend on it.
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
import { buildCopilotCards } from '../../copilotCards';
import { buildCopilotTools } from '../../copilotTools';
import { buildViewTools } from '../../copilotViewTools';
import { handlePcaMessage } from '../../workers/pcaKdeWorker';
import {
    appSource, bundledSkills, CARD_TOLERANCE_A, COPILOT_PAGES, jsonFile, loadDemoAssistantRun, PLOT_ORDER, promptVariants,
    recordingModel, renderRequest
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

    it('mirrors what App.jsx and Dashboard.jsx do not export', () => {
        const app = appSource('App.jsx');
        const pages = /const COPILOT_PAGES = \{([^}]*)\}/.exec(app)?.[1] ?? '';
        expect(Object.fromEntries([...pages.matchAll(/(\w+): '([^']*)'/g)].map((match) => [match[1], match[2]]))).toEqual(COPILOT_PAGES);
        expect(app).toContain(`const symTolState = useState(${CARD_TOLERANCE_A});`);
        const dashboard = appSource('components/Dashboard.jsx');
        expect(dashboard).toContain(`const plotOrder = [${PLOT_ORDER.map((kind) => `'${kind}'`).join(', ')}];`);
        expect(dashboard).toContain('const isDashboardPlotFile = (file) => file.plotKind && file.plotKind !== \'stog\';');
        expect(dashboard).toContain('.filter(isDashboardPlotFile)\n            .sort(comparePlotFiles);');
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

        const { run, rmc6fText } = demo;
        const request = (kind, params) => handlePcaMessage({ kind, ...params }, async () => rmc6fText);
        const cards = buildCopilotCards({ plotFiles: run.plotFiles, rValueFile: run.rValueFile, symmetry: run.symmetry });
        const tools = [
            ...buildCopilotTools({ request, structure: run.structure, elements: Object.keys(run.structure.elementCounts || {}).sort() }),
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
