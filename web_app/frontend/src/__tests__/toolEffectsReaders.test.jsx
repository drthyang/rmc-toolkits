// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The dock and the validation harness read their tools from the one list of
// tool effects (llm/copilot/effects.js): each passes the tools it built
// through toolsWithEffects — the dock with every effect (DOCK_EFFECTS), the
// live arms with 'compute' only (LIVE_ARM_EFFECTS) — and offers what comes
// back. The function is spied on here, its own behaviour kept.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import React, { act } from 'react';
import { createRoot } from 'react-dom/client';

vi.mock('../llm/copilot/effects', async (importOriginal) => {
    const real = await importOriginal();
    return { ...real, toolsWithEffects: vi.fn(real.toolsWithEffects) };
});

const { DOCK_EFFECTS, toolsWithEffects } = await import('../llm/copilot/effects');
const { default: CopilotChat } = await import('../llm/components/CopilotChat');
const { DOCK_STORAGE_KEY } = await import('../llm/components/dockLayout');
const { saveSettings } = await import('../llm/settings');
const { buildCase, DEV_CASES } = await import('../../validation/cases');
const { LIVE_ARM_EFFECTS, liveArmTools } = await import('../../validation/harness');
const { loadRun } = await import('../../validation/runFolder');

const changeTool = {
    name: 'set_bond_geometry',
    description: 'Set the Bond Geometry inputs.',
    parameters: { type: 'object', properties: {} },
    run: vi.fn(),
    undo: vi.fn()
};
const viewTool = { name: 'current_view', view: true, description: 'What is on screen.', parameters: { type: 'object', properties: {} }, run: vi.fn() };
const computeTool = { name: 'list_sites', description: 'Sites.', parameters: { type: 'object', properties: {} }, run: vi.fn() };

describe('the dock and the harness read their tools from the one list', () => {
    let container;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        toolsWithEffects.mockClear();
        localStorage.removeItem(DOCK_STORAGE_KEY);
        saveSettings({ baseUrl: 'http://localhost:11434/v1', model: '', apiKey: '', copilotTools: true });
        vi.stubGlobal('fetch', vi.fn(async (url) => (String(url).endsWith('/v1/models')
            ? { ok: true, status: 200, json: async () => ({ data: [{ id: 'm' }] }) }
            : { ok: false, status: 404, json: async () => ({}) })));
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.unstubAllGlobals();
    });

    it('the dock: every tool it built, through the list, with every effect', async () => {
        act(() => root.render(
            <CopilotChat
                runName="demo"
                plotFiles={[{ plotKind: 'neutron_sq', plotData: { title: 'S(Q)', metrics: { rwp: 0.04 } } }]}
                tools={[computeTool, viewTool, changeTool]}
                runKey="run-1"
            />
        ));
        act(() => container.querySelector('.llm-chat-launcher').click());
        await act(async () => { await new Promise((resolve) => setTimeout(resolve, 0)); });
        const call = toolsWithEffects.mock.calls.find(([tools]) => tools.includes(changeTool));
        expect(call).toBeDefined();
        expect(call[1]).toBe(DOCK_EFFECTS);
        expect(call[0]).toEqual(expect.arrayContaining([computeTool, viewTool, changeTool]));
        // What it offers is what the list returned: the change tool among them.
        expect([...container.querySelectorAll('.llm-tool-list code')].map((item) => item.textContent))
            .toEqual(expect.arrayContaining(['list_sites', 'current_view', 'set_bond_geometry']));
    });

    it('the harness: the tools it built for a case, through the list, compute only', async () => {
        const run = await loadRun(buildCase(DEV_CASES[0]));
        const injected = { ...run, tools: [...run.tools, changeTool, viewTool] };
        const tools = liveArmTools(injected);
        expect(toolsWithEffects).toHaveBeenCalledTimes(1);
        const [[built, effects]] = toolsWithEffects.mock.calls;
        expect(effects).toBe(LIVE_ARM_EFFECTS);
        expect(built).toEqual(expect.arrayContaining([changeTool, viewTool]));
        expect(tools).not.toContain(changeTool);
        expect(tools).not.toContain(viewTool);
        expect(tools.length).toBe(built.length - 2);
    });
});
