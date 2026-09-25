// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// Flask-mode Live Data must keep the analysis pages on ONE configuration.
// They fetch from the backend on demand and the backend always reads the file
// on disk, so when RMCProfile saves a new .rmc6f a page that is not refreshed
// would draw its old site table / slab points next to data its next request
// (a slider move, a site click) takes from the new configuration. App.jsx
// keys the analysis pages on a configuration epoch that changes with the
// .rmc6f signature in the /api/files listing, so they remount and re-read
// everything from the new file.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';

const state = vi.hoisted(() => ({
    files: [],
    mounts: { structure: 0, geometry: 0, ellipsoids: 0, orientation: 0 },
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

vi.mock('../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => false,
    supportsFileSystemAccess: () => false,
}));

const countingPage = async (name) => {
    const { useEffect } = await import('react');
    return {
        default: function CountingPage() {
            useEffect(() => { state.mounts[name] += 1; }, []);
            return null;
        },
    };
};

vi.mock('../components/StructurePage', () => countingPage('structure'));
vi.mock('../components/BondGeometryPage', () => countingPage('geometry'));
vi.mock('../components/PcaKdePage', () => countingPage('ellipsoids'));
vi.mock('../components/OrientationPage', () => countingPage('orientation'));
vi.mock('../components/Dashboard', () => ({ default: () => null }));
vi.mock('../components/AutoStogPage', () => ({ default: () => null }));
vi.mock('../llm', () => ({ AssistantPage: () => null }));

const { default: App } = await import('../App');
const { WATCH_INTERVAL_MS } = await import('../browserData');

const listing = (rmc6f, log) => [
    { name: 'run.rmc6f', path: '/runs/a/run.rmc6f', type: 'file', plotKind: null, ...rmc6f },
    { name: 'run-01.log', path: '/runs/a/run-01.log', type: 'file', plotKind: 'r_value', ...log },
];

describe('App Flask-mode Live Data', () => {
    let container;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        vi.useFakeTimers();
        Object.keys(state.mounts).forEach((key) => { state.mounts[key] = 0; });
        state.files = listing({ modified: 100, size: 5000 }, { modified: 100, size: 10 });
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.useRealTimers();
    });

    const click = async (element) => {
        await act(async () => { element.click(); });
    };
    const tab = (label) => [...container.querySelectorAll('nav.page-tabs button')]
        .find((button) => button.textContent.trim() === label);
    const poll = async () => {
        await act(async () => { await vi.advanceTimersByTimeAsync(WATCH_INTERVAL_MS); });
    };

    const openAnalysisPagesWithLiveData = async () => {
        await act(async () => { root.render(<App />); });
        for (const label of ['Atomic Density', 'Bond Geometry', 'PCA Ellipsoid', 'Displacement Directions']) {
            await click(tab(label));
        }
        await click(container.querySelector('label.watch-toggle input[type="checkbox"]'));
        await poll();
        expect(state.mounts).toEqual({ structure: 1, geometry: 1, ellipsoids: 1, orientation: 1 });
    };

    it('re-reads the analysis pages when the .rmc6f changes on disk', async () => {
        await openAnalysisPagesWithLiveData();

        // RMCProfile saves a new configuration.
        state.files = listing({ modified: 200, size: 5100 }, { modified: 100, size: 10 });
        await poll();
        expect(state.mounts).toEqual({ structure: 2, geometry: 2, ellipsoids: 2, orientation: 2 });

        // Nothing new on disk: no further refresh.
        await poll();
        expect(state.mounts).toEqual({ structure: 2, geometry: 2, ellipsoids: 2, orientation: 2 });
    });

    it('catches up when Live Data is switched on after the file changed', async () => {
        await act(async () => { root.render(<App />); });
        await click(tab('PCA Ellipsoid'));
        await poll();
        expect(state.mounts.ellipsoids).toBe(1);

        // Saved while Live Data was off: nothing is polled...
        state.files = listing({ modified: 200, size: 5100 }, { modified: 100, size: 10 });
        await poll();
        expect(state.mounts.ellipsoids).toBe(1);

        // ...until it is switched on, which re-reads the page at once.
        await click(container.querySelector('label.watch-toggle input[type="checkbox"]'));
        expect(state.mounts.ellipsoids).toBe(2);
    });

    it('re-checks the .rmc6f when the same folder is loaded again with Live Data off', async () => {
        await act(async () => { root.render(<App />); });
        await click(tab('PCA Ellipsoid'));
        await poll();
        expect(state.mounts.ellipsoids).toBe(1);
        const load = container.querySelector('form.path-bar button[type="submit"]');

        // Load the same folder with nothing new on disk: the page keeps its state.
        await click(load);
        expect(state.mounts.ellipsoids).toBe(1);

        // A configuration saved while Live Data is off is picked up by Load.
        state.files = listing({ modified: 200, size: 5100 }, { modified: 100, size: 10 });
        await poll();
        expect(state.mounts.ellipsoids).toBe(1);
        await click(load);
        expect(state.mounts.ellipsoids).toBe(2);
    });

    it('does not reload the analysis pages when only the plot/log files change', async () => {
        await openAnalysisPagesWithLiveData();

        state.files = listing({ modified: 100, size: 5000 }, { modified: 300, size: 99 });
        await poll();
        expect(state.mounts).toEqual({ structure: 1, geometry: 1, ellipsoids: 1, orientation: 1 });
    });
});
