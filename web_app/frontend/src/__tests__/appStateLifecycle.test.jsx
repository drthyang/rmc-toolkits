// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// What App.jsx keeps, and for how long (plan 1.4, state lifecycle):
//
//   - the AI Copilot's run (`runKey`) and configuration (`configKey`): a Live
//     Data save that brings a new .rmc6f — or a Load of the same folder that
//     finds one — is a new configuration of the same run, never a new run;
//     a save of the plots or the .log alone is neither. Both runtimes: Flask
//     (configEpoch) and static-mode Live Data on a picked folder (the .rmc6f
//     files' signature). What the dock does with configKey is pinned in
//     src/llm/__tests__/copilotChat.test.jsx.
//   - the Detected SG pick (SymTolContext, which the card and the Copilot's
//     context both read): kept across a Live Data save and a Load of the same
//     run, back at the card's default (DEFAULT_VIEW.toleranceA) when another
//     run opens, even one that was open before.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';

const state = vi.hoisted(() => ({
    static: false,
    fsAccess: false,
    files: [],
    folder: [],
    copilot: null,
    copilotView: null,
    symTol: null,
}));

vi.mock('axios', () => ({
    default: {
        get: vi.fn(async (url) => (
            url.endsWith('/api/files') ? { data: { files: state.files } } : { data: {} }
        )),
        post: vi.fn(async () => ({ data: {} })),
        isCancel: () => false,
    },
}));

// A browser-loaded run as makeRunFromEntries shapes it, from { name, modified, size }.
const localRunOf = (name, entries) => ({
    name,
    files: entries.map((entry) => ({
        name: entry.name,
        path: `${name}/${entry.name}`,
        type: 'file',
        plotKind: entry.name.endsWith('.log') ? 'r_value' : null,
        size: entry.size,
        modified: entry.modified,
    })),
    structure: null,
    structureFile: null,
    settingsFile: null,
    structureError: null,
    diagnostics: {},
});

vi.mock('../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => state.static,
    supportsFileSystemAccess: () => state.fsAccess,
    buildLocalRunFromHandle: async (handle) => localRunOf(handle.name, state.folder),
    loadDemoRun: async () => localRunOf('Demo', [
        { name: 'GTS_250K.rmc6f', modified: 1, size: 1000 },
        { name: 'GTS_250K-01.log', modified: 1, size: 10 },
    ]),
}));

// The Dashboard reads the shared tolerance as the Detected SG card does.
vi.mock('../components/Dashboard', async () => {
    const { useContext } = await import('react');
    const { SymTolContext } = await import('../symTolContext');
    return {
        default: function TolerancePage() {
            state.symTol = useContext(SymTolContext);
            return null;
        },
    };
});
vi.mock('../components/StructurePage', () => ({ default: () => null }));
vi.mock('../components/BondGeometryPage', () => ({ default: () => null }));
vi.mock('../components/PcaKdePage', () => ({ default: () => null }));
vi.mock('../components/OrientationPage', () => ({ default: () => null }));
vi.mock('../components/AutoStogPage', () => ({ default: () => null }));
vi.mock('../llm', () => ({
    CopilotChat: ({ runKey, configKey }) => {
        state.copilot = { runKey, configKey };
        return null;
    },
}));
// The Copilot's run module: what the context is built from (its view).
vi.mock('../useCopilotRun', () => ({
    default: ({ view }) => {
        state.copilotView = view;
        return { props: { plotFiles: [], rValueFile: null, symmetry: null }, tools: [] };
    },
}));

const { default: App } = await import('../App');
const { WATCH_INTERVAL_MS } = await import('../browserData');
const { DEFAULT_VIEW } = await import('../copilotRun');

const listing = (rmc6f, log) => [
    { name: 'run.rmc6f', path: '/runs/a/run.rmc6f', type: 'file', plotKind: null, ...rmc6f },
    { name: 'run-01.log', path: '/runs/a/run-01.log', type: 'file', plotKind: 'r_value', ...log },
];

describe('App state lifecycle', () => {
    let container;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        vi.useFakeTimers();
        state.static = false;
        state.fsAccess = false;
        state.files = listing({ modified: 100, size: 5000 }, { modified: 100, size: 10 });
        state.folder = [];
        state.copilot = null;
        state.copilotView = null;
        state.symTol = null;
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.useRealTimers();
        delete window.showDirectoryPicker;
    });

    const click = async (element) => {
        await act(async () => { element.click(); });
    };
    const button = (label) => [...container.querySelectorAll('button')]
        .find((element) => element.textContent.trim() === label);
    const liveData = () => container.querySelector('label.ui-switch-outline input[type="checkbox"]');
    const poll = async () => {
        await act(async () => { await vi.advanceTimersByTimeAsync(WATCH_INTERVAL_MS); });
    };
    const pick = async (tolerance) => {
        await act(async () => { state.symTol[1](tolerance); });
    };
    // The pick as the card and as the Copilot's context read it.
    const tolerance = () => {
        expect(state.copilotView.toleranceA).toBe(state.symTol[0]);
        return state.symTol[0];
    };
    const openFolder = async (path) => {
        const input = container.querySelector('#data-path');
        await act(async () => {
            Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set.call(input, path);
            input.dispatchEvent(new Event('input', { bubbles: true }));
        });
        await act(async () => {
            container.querySelector('form.path-bar').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true }));
        });
    };

    describe('Flask mode', () => {
        it('a Live Data save of a new .rmc6f is a new configuration of the same run, and keeps the Detected SG pick', async () => {
            await act(async () => { root.render(<App />); });
            await click(liveData());
            await poll();
            const { runKey, configKey } = state.copilot;
            expect(tolerance()).toBe(DEFAULT_VIEW.toleranceA);
            await pick(0.05);
            expect(tolerance()).toBe(0.05);

            // RMCProfile saves a new configuration.
            state.files = listing({ modified: 200, size: 5100 }, { modified: 100, size: 10 });
            await poll();
            expect(state.copilot.runKey).toBe(runKey);
            expect(state.copilot.configKey).not.toBe(configKey);
            expect(tolerance()).toBe(0.05);

            // Only the .log changed: the same configuration.
            const saved = state.copilot.configKey;
            state.files = listing({ modified: 200, size: 5100 }, { modified: 300, size: 99 });
            await poll();
            expect(state.copilot).toEqual({ runKey, configKey: saved });

            // Load of the same folder: the same run, the pick kept.
            await click(button('Load'));
            expect(state.copilot).toEqual({ runKey, configKey: saved });
            expect(tolerance()).toBe(0.05);
        });

        it('a Load that finds a new .rmc6f is a new configuration too (Live Data off)', async () => {
            await act(async () => { root.render(<App />); });
            await poll();
            const { runKey, configKey } = state.copilot;
            state.files = listing({ modified: 200, size: 5100 }, { modified: 100, size: 10 });
            await click(button('Load'));
            expect(state.copilot.runKey).toBe(runKey);
            expect(state.copilot.configKey).not.toBe(configKey);
        });

        it('another run resets the Detected SG pick to the card\'s default, even a run open before', async () => {
            await act(async () => { root.render(<App />); });
            await poll();
            await pick(0.05);
            const first = state.copilot.runKey;

            await openFolder('other');
            expect(state.copilot.runKey).not.toBe(first);
            expect(tolerance()).toBe(DEFAULT_VIEW.toleranceA);

            await pick(0.4);
            await openFolder('data');
            expect(state.copilot.runKey).toBe(first);
            expect(tolerance()).toBe(DEFAULT_VIEW.toleranceA);

            // The Demo is another run as well, and so is leaving it.
            await pick(0.05);
            await click(container.querySelector('button.demo-button'));
            expect(state.copilot.runKey).not.toBe(first);
            expect(state.copilot.configKey).toMatch(/^local\|/);
            expect(tolerance()).toBe(DEFAULT_VIEW.toleranceA);
            await pick(0.05);
            await click(container.querySelector('button.demo-button'));
            expect(state.copilot.runKey).toBe(first);
            expect(tolerance()).toBe(DEFAULT_VIEW.toleranceA);
        });
    });

    describe('static mode, a picked folder with Live Data', () => {
        it('a new .rmc6f in the folder is a new configuration of the same run; picking a folder is another run', async () => {
            state.static = true;
            state.fsAccess = true;
            state.folder = [
                { name: 'run.rmc6f', modified: 100, size: 5000 },
                { name: 'run-01.log', modified: 100, size: 10 },
            ];
            window.showDirectoryPicker = async () => ({ name: 'picked' });
            await act(async () => { root.render(<App />); });
            await click(button('Select Folder'));
            await click(liveData());
            const { runKey, configKey } = state.copilot;
            expect(configKey).toMatch(/^local\|/);
            await pick(0.05);

            // Only the .log changed: the run refreshes, the configuration is the same.
            state.folder = [
                { name: 'run.rmc6f', modified: 100, size: 5000 },
                { name: 'run-01.log', modified: 200, size: 20 },
            ];
            await poll();
            expect(state.copilot).toEqual({ runKey, configKey });

            // A new configuration saved into the folder.
            state.folder = [
                { name: 'run.rmc6f', modified: 300, size: 5100 },
                { name: 'run-01.log', modified: 300, size: 30 },
            ];
            await poll();
            expect(state.copilot.runKey).toBe(runKey);
            expect(state.copilot.configKey).not.toBe(configKey);
            expect(tolerance()).toBe(0.05);

            // Picking a folder (the same one too) opens another run.
            await click(button('Select Folder'));
            expect(state.copilot.runKey).not.toBe(runKey);
            expect(tolerance()).toBe(DEFAULT_VIEW.toleranceA);
        });
    });
});
