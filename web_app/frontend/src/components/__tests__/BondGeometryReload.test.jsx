// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// Flask-mode Live Data reloads the Bond Geometry page IN PLACE: a new
// `dataEpoch` (a new .rmc6f saved under the same run directory) re-reads the
// element list, the Model information card and the partials, keeps the
// user's triplet and typed windows, and drops the computed angle
// distribution — a result from the previous configuration must never sit
// next to the new model — with a cue saying why.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';

const state = vi.hoisted(() => ({ requests: [] }));

const SITES = {
    elements: ['Ga', 'Nb', 'Se'],
    sites: [
        { referenceNumber: 1, element: 'Ga', count: 4, copiesPerCell: 4 },
        { referenceNumber: 2, element: 'Nb', count: 16, copiesPerCell: 16 },
        { referenceNumber: 3, element: 'Se', count: 32, copiesPerCell: 32 },
    ],
};

const TRIPLETS = {
    triplet: ['Se', 'Nb', 'Se'],
    bond12: [2.0, 3.0],
    bond23: [2.0, 3.0],
    sharedEnds: true,
    binCenters: [0.5, 1.5],
    sinCorrected: [0, 1],
    density: [0, 1],
    coordination: [0, 0, 1],
    apexCount: 16,
    lengths12: { uniqueBonds: 96, meanLength: 2.6 },
    angleCount: 240,
    meanAngle: 90,
    stdAngle: 5,
};

vi.mock('axios', () => ({
    default: {
        get: vi.fn(async (url, config) => {
            state.requests.push({ url, params: config?.params });
            if (url.endsWith('/api/pca/sites')) return { data: SITES };
            if (url.endsWith('/api/triplets')) return { data: TRIPLETS };
            if (url.endsWith('/api/files')) return { data: { files: [] } };
            if (url.endsWith('/api/structure')) return { data: { atomCount: 52, elements: SITES.elements } };
            throw new Error(`unexpected ${url}`);
        }),
    },
}));

// Flask mode (the dev test environment would otherwise count as static).
vi.mock('../../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => false,
}));

// The Three.js / SVG / model-card panels are not what this test is about.
vi.mock('../FoldedCellPanel', () => ({ default: () => null }));
vi.mock('../InteractivePlot', () => ({ default: () => null }));
vi.mock('../ModelSummary', () => ({ default: () => null }));

const { default: BondGeometryPage } = await import('../BondGeometryPage');

describe('BondGeometryPage dataEpoch (Flask Live Data)', () => {
    let container;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        state.requests = [];
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
    });

    const render = async (dataEpoch) => {
        await act(async () => {
            root.render(<BondGeometryPage directory="runs/a" localRun={null} dataEpoch={dataEpoch} />);
        });
    };
    const count = (suffix) => state.requests.filter((request) => request.url.endsWith(suffix)).length;
    const select = (label) => container.querySelector(`select[aria-label="${label}"]`);
    const setSelect = async (label, value) => {
        await act(async () => {
            const element = select(label);
            element.value = value;
            element.dispatchEvent(new Event('change', { bubbles: true }));
        });
    };
    const text = () => container.textContent.replace(/\s+/g, ' ');
    // The KPI rail is always there; before a result its values read "—".
    const angleKpi = () => container.querySelector('[aria-label="Triplet result"] .ui-kpi__value').textContent;

    it('reloads in place, keeps the triplet and drops the stale result', async () => {
        await render(0);
        expect(count('/api/pca/sites')).toBe(1);
        expect(count('/api/structure')).toBe(1);
        expect(count('/api/files')).toBe(1);

        // The user picks a triplet other than the default and computes it.
        await setSelect('End element A', 'Ga');
        const compute = [...container.querySelectorAll('button')].find((button) => /Compute/.test(button.textContent));
        expect(angleKpi()).toBe('—');
        await act(async () => { compute.click(); });
        expect(count('/api/triplets')).toBe(1);
        expect(angleKpi()).toBe('15.0\u2009per Nb');

        // Nothing new on disk: nothing is re-read.
        await render(0);
        expect(count('/api/pca/sites')).toBe(1);

        // RMCProfile saves a new configuration of the same run.
        await render(1);
        expect(count('/api/pca/sites')).toBe(2);
        expect(count('/api/structure')).toBe(2);
        expect(count('/api/files')).toBe(2);
        // Never recomputed unasked...
        expect(count('/api/triplets')).toBe(1);
        // ...the previous result is gone, and the page says why.
        expect(angleKpi()).toBe('—');
        expect(text()).toMatch(/New configuration — Compute again/);
        expect(text()).toMatch(/new configuration/);
        // The picks survive the reload.
        expect(select('End element A').value).toBe('Ga');
        expect(select('Central element B').value).toBe('Nb');

        // Compute again: the cue goes, the new result lands.
        await act(async () => { compute.click(); });
        expect(count('/api/triplets')).toBe(2);
        expect(angleKpi()).toBe('15.0\u2009per Nb');
        expect(text()).not.toMatch(/new configuration/i);
    });
});
