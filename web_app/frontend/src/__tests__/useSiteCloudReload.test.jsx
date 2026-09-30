// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// Flask-mode Live Data reloads the site-cloud pages IN PLACE: a new
// `dataEpoch` (App.jsx's configuration epoch) re-requests the site table from
// the backend under the unchanged directory, gives `requestPca` a new identity
// (so the PCA KDE volume and the orientation histogram, whose effects are
// keyed on it, re-request too), and keeps the picked site when the new
// configuration still has it.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';

const state = vi.hoisted(() => ({ requests: [], refs: [1, 2, 3] }));

vi.mock('axios', () => ({
    default: {
        get: vi.fn(async (url, config) => {
            state.requests.push({ url, params: config?.params });
            return {
                data: {
                    elements: ['Ga'],
                    sites: state.refs.map((referenceNumber) => ({
                        referenceNumber, element: 'Ga', count: 8, copiesPerCell: 8,
                    })),
                },
            };
        }),
    },
}));

// Flask mode: a backend directory (the test environment would otherwise count
// as the static build, where no run means no request).
vi.mock('../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => false,
}));

const { default: useSiteCloud } = await import('../useSiteCloud');

describe('useSiteCloud dataEpoch (Flask Live Data)', () => {
    let container;
    let root;
    let hook;

    function Probe(props) {
        hook = useSiteCloud(props);
        return null;
    }

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        state.requests = [];
        state.refs = [1, 2, 3];
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
    });

    const render = async (props) => {
        await act(async () => { root.render(<Probe {...props} />); });
    };
    const siteRequests = () => state.requests.filter((request) => request.url.endsWith('/api/pca/sites'));

    it('re-reads the site table in place and keeps the picked site', async () => {
        await render({ directory: 'runs/a', localRun: null, dataEpoch: 0 });
        expect(siteRequests()).toHaveLength(1);
        await act(async () => { hook.setSelectedRef(3); });
        const firstRequest = hook.requestPca;

        // Re-render with nothing new: no reload.
        await render({ directory: 'runs/a', localRun: null, dataEpoch: 0 });
        expect(siteRequests()).toHaveLength(1);
        expect(hook.requestPca).toBe(firstRequest);

        // A new configuration on disk: same directory, new epoch.
        await render({ directory: 'runs/a', localRun: null, dataEpoch: 1 });
        expect(siteRequests()).toHaveLength(2);
        expect(siteRequests()[1].params.dir).toBe('runs/a');
        expect(hook.requestPca).not.toBe(firstRequest);
        expect(hook.selectedRef).toBe(3);
        expect(hook.datasetKey).toBe('runs/a');
    });

    it('falls back to a site that still exists when the pick is gone', async () => {
        await render({ directory: 'runs/a', localRun: null, dataEpoch: 0 });
        await act(async () => { hook.setSelectedRef(3); });

        state.refs = [1, 2];
        await render({ directory: 'runs/a', localRun: null, dataEpoch: 1 });
        expect(hook.selectedRef).toBe(1);
    });
});
