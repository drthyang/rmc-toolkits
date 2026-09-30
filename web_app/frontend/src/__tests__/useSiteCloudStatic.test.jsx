// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The static build has no backend. With no run open, useSiteCloud must not
// ask /api/pca/sites (that showed "Request failed with status code 404" on the
// PCA Ellipsoid, Displacement Directions and Bond Geometry pages) and must
// not be ready, so no page issues a request either.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';

const requests = vi.hoisted(() => []);
vi.mock('axios', () => ({
    default: {
        get: vi.fn(async (url) => {
            requests.push(url);
            throw Object.assign(new Error('Request failed with status code 404'), { response: { status: 404, data: {} } });
        }),
    },
}));
vi.mock('../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => true,
}));

const { default: useSiteCloud } = await import('../useSiteCloud');

describe('useSiteCloud in the static build with no run', () => {
    let container;
    let root;
    let hook;
    const Probe = (props) => {
        hook = useSiteCloud(props);
        return null;
    };

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        requests.length = 0;
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });
    afterEach(() => {
        act(() => root.unmount());
        container.remove();
    });

    it('asks nothing, reports no error and is not ready', async () => {
        await act(async () => { root.render(<Probe directory="data" localRun={null} />); });
        expect(requests).toEqual([]);
        expect(hook.sitesError).toBeNull();
        expect(hook.sites).toBeNull();
        expect(hook.loadingSites).toBe(false);
        expect(hook.ready).toBe(false);
    });

    it('drops the site pick when the run is closed, so no page asks for it', async () => {
        // A run whose file text is still loading: no worker, no request.
        const pending = { runId: 1, structureFile: { path: 'Demo/x.rmc6f', sourceFile: { text: () => new Promise(() => {}) } } };
        await act(async () => { root.render(<Probe directory="data" localRun={pending} />); });
        await act(async () => { hook.setSelectedRef(3); });
        expect(hook.selectedRef).toBe(3);
        await act(async () => { root.render(<Probe directory="data" localRun={null} />); });
        expect(hook.selectedRef).toBeNull();
        expect(hook.ready).toBe(false);
        expect(requests).toEqual([]);
    });
});
