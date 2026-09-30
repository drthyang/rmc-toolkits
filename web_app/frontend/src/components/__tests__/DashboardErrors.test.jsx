// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// Dashboard failures must reach the screen as one line: a failed "Save all
// figures" (the save menu does not await it, so an uncaught rejection used to
// vanish) says so under the Loaded-files header until the next save, and a
// run-folder listing that fails without a server message (server down, network
// error) says so instead of leaving only the "Open a run folder." prompt.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';

const state = vi.hoisted(() => ({ listError: null, zipOutcomes: [] }));

vi.mock('axios', () => ({
    default: {
        get: vi.fn(async (url) => {
            if (url.endsWith('/api/files')) {
                if (state.listError) throw state.listError;
                return {
                    data: {
                        files: [{ name: 'run_FQ1.csv', path: 'runs/a/run_FQ1.csv', type: 'file', plotKind: 'xray_sq' }],
                    },
                };
            }
            if (url.endsWith('/api/plot/metadata')) return { data: { title: 'F(Q)' } };
            if (url.endsWith('/api/structure')) throw Object.assign(new Error('404'), { response: { data: {} } });
            throw new Error(`unexpected ${url}`);
        }),
    },
}));

// Flask mode (the test environment would otherwise count as static).
vi.mock('../../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => false,
}));

vi.mock('../../figureExport', () => ({
    saveSvgFiguresAsZip: vi.fn(async () => {
        const outcome = state.zipOutcomes.shift();
        if (outcome) throw outcome;
    }),
}));

// The chart itself is not under test: a bare svg is enough for Save all.
vi.mock('../InteractivePlot', () => ({
    default: () => <div className="interactive-plot"><svg /></div>,
}));
vi.mock('../ModelSummary', () => ({ default: () => null }));
vi.mock('../../llm', () => ({ WatchdogBadge: () => null }));

const { default: Dashboard } = await import('../Dashboard');

const flush = () => act(async () => {
    for (let i = 0; i < 5; i += 1) await Promise.resolve();
});

describe('Dashboard error lines', () => {
    let container;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        state.listError = null;
        state.zipOutcomes = [];
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
    });

    const mount = async () => {
        act(() => root.render(<Dashboard directory="runs/a" localRun={null} />));
        await flush();
    };

    const saveAll = async () => {
        act(() => container.querySelector('.loaded-files-card .ui-save__trigger').click());
        await act(async () => {
            container.querySelector('.loaded-files-card [role="menuitem"]').click();
        });
        await flush();
    };

    const saveAllBanner = () => container.querySelector('.loaded-files-card .ui-banner--danger');

    it('says when Save all figures fails, until the next save', async () => {
        await mount();
        expect(container.querySelector('.loaded-files-card')).not.toBeNull();

        state.zipOutcomes = [new Error('Could not rasterize the figure')];
        await saveAll();
        expect(saveAllBanner()?.textContent).toContain('Could not rasterize the figure');
        expect(saveAllBanner()?.getAttribute('role')).toBe('alert');
        // The menu is usable again (not stuck on "Saving…").
        expect(container.querySelector('.loaded-files-card .ui-save__trigger').textContent).toContain('Save all figures');

        await saveAll();
        expect(saveAllBanner()).toBeNull();
    });

    it('drops a failed Save all when another run folder opens', async () => {
        await mount();
        state.zipOutcomes = [new Error('Could not rasterize the figure')];
        await saveAll();
        expect(saveAllBanner()).not.toBeNull();

        act(() => root.render(<Dashboard directory="runs/b" localRun={null} />));
        await flush();
        expect(saveAllBanner()).toBeNull();
    });

    it('says when the run folder cannot be listed and the server gave no reason', async () => {
        state.listError = new Error('Network Error');
        await mount();
        const alert = container.querySelector('[role="alert"]');
        expect(alert?.textContent).toContain('Could not list the run folder');
    });

    it('shows the server\'s own reason when it gives one', async () => {
        state.listError = Object.assign(new Error('403'), {
            response: { data: { error: 'Path is outside the data root' } },
        });
        await mount();
        const alert = container.querySelector('[role="alert"]');
        expect(alert?.textContent).toContain('Path is outside the data root');
        expect(alert?.textContent).not.toContain('Could not list the run folder');
    });

    it('falls back to a generic line when the save error has no message', async () => {
        await mount();
        state.zipOutcomes = [new Error('')];
        await saveAll();
        expect(saveAllBanner()?.textContent).toContain('Could not save the figures');
    });
});
