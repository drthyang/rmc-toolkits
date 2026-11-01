// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The Bond Geometry page's controller (bondGeometryControl.js), the one way
// the AI Copilot's change tool reaches the inputs: registered while the page
// is mounted, it reads the boxes and the result shown, and applies a change
// through the page's own setters — the selects and boxes then show it — and
// computes as the Compute button does. Flask mode, with axios mocked.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import { getBondGeometry } from '../../bondGeometryControl';

const state = vi.hoisted(() => ({ requests: [], sites: null }));

const SITES = {
    elements: ['Ga', 'Nb', 'Se'],
    sites: [
        { referenceNumber: 1, element: 'Ga', count: 4, copiesPerCell: 4 },
        { referenceNumber: 2, element: 'Nb', count: 16, copiesPerCell: 16 },
        { referenceNumber: 3, element: 'Se', count: 32, copiesPerCell: 32 },
    ],
};

// A triplets payload for the request asked: its triplet and windows echoed.
const triplets = (params) => ({
    triplet: [params.end1, params.apex, params.end2],
    bond12: [Number(params.r12Min), Number(params.r12Max)],
    bond23: [Number(params.r23Min ?? params.r12Min), Number(params.r23Max ?? params.r12Max)],
    sharedEnds: params.end1 === params.end2 && params.r23Min == null,
    binCenters: [0.5, 1.5],
    binWidth: Number(params.binWidth),
    sinCorrected: [0, 1],
    density: [0, 1],
    coordination: [0, 0, 1],
    coordinationSummary: { mean: 2, modal: 2, modalShare: 1, unbondedShare: 0 },
    apexCount: 16,
    lengths12: { uniqueBonds: 96, meanLength: 2.6, count: 96 },
    lengths23: { uniqueBonds: 96, meanLength: 2.6, count: 96 },
    angleCount: params.end1 === 'Ga' ? 64 : 240,
    meanAngle: 90,
    stdAngle: 5,
});

vi.mock('axios', () => ({
    default: {
        get: vi.fn(async (url, config) => {
            state.requests.push({ url, params: config?.params });
            if (url.endsWith('/api/pca/sites')) return { data: await (state.sites ?? SITES) };
            if (url.endsWith('/api/triplets')) return { data: triplets(config.params) };
            if (url.endsWith('/api/files')) return { data: { files: [] } };
            if (url.endsWith('/api/structure')) return { data: { atomCount: 52, elements: SITES.elements } };
            throw new Error(`unexpected ${url}`);
        }),
    },
}));

vi.mock('../../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => false,
}));
vi.mock('../FoldedCellPanel', () => ({ default: () => null }));
vi.mock('../InteractivePlot', () => ({ default: () => null }));
vi.mock('../ModelSummary', () => ({ default: () => null }));

const { default: BondGeometryPage } = await import('../BondGeometryPage');

const flush = () => act(async () => {
    for (let i = 0; i < 6; i += 1) await new Promise((resolve) => setTimeout(resolve, 0));
});

// A change applied: started inside act, then the page renders and computes
// (an apply resolves only once the page has re-rendered), then its outcome.
const applied = async (controller, inputs, options) => {
    let promise;
    await act(async () => { promise = controller.apply(inputs, options); });
    await flush();
    return promise;
};

describe('the Bond Geometry page\'s controller', () => {
    let container;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        state.requests = [];
        state.sites = null;
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root?.unmount());
        container.remove();
    });

    const render = async (props = {}) => {
        await act(async () => {
            root.render(<BondGeometryPage directory="runs/a" localRun={null} dataEpoch={0} {...props} />);
        });
        await flush();
    };
    const select = (label) => container.querySelector(`select[aria-label="${label}"]`);
    const box = (label) => container.querySelector(`input[aria-label="${label}"]`);
    const tripletRequests = () => state.requests.filter((request) => request.url.endsWith('/api/triplets'));
    const angleKpi = () => container.querySelector('[aria-label="Triplet result"] .ui-kpi__value').textContent;

    const NEXT = {
        end1: 'Ga', apex: 'Nb', end2: 'Se', r12Min: '2.4', r12Max: '2.9', split23: true, r23Min: '2.5', r23Max: '2.8', binWidth: '2'
    };

    it('registers while the page is mounted, and reads the page: its run, its boxes, nothing shown', async () => {
        expect(getBondGeometry()).toBeNull();
        await render();
        const controller = getBondGeometry();
        expect(controller).not.toBeNull();
        expect(controller.state()).toEqual({
            noRun: false,
            ready: true,
            sitesError: null,
            elements: ['Ga', 'Nb', 'Se'],
            // Seeded as the page seeds them: the most abundant element at the ends.
            inputs: { end1: 'Se', apex: 'Nb', end2: 'Se', r12Min: '2.00', r12Max: '3.00', split23: false, r23Min: '2.00', r23Max: '3.00', binWidth: '1.0' },
            computing: false,
            shown: null,
            result: null,
            epoch: expect.any(Number)
        });
        act(() => root.unmount());
        root = null;
        expect(getBondGeometry()).toBeNull();
    });

    it('applies a change through the page\'s own setters, then computes it as Compute does', async () => {
        await render();
        const controller = getBondGeometry();
        const outcome = await applied(controller, NEXT);
        // The controls show the change, as if the user had set them.
        expect(select('End element A').value).toBe('Ga');
        expect(select('Central element B').value).toBe('Nb');
        expect(box('A-B window minimum').value).toBe('2.4');
        expect(box('A-B window maximum').value).toBe('2.9');
        expect(container.querySelector('input[aria-label="Use a distinct B-C window"]').checked).toBe(true);
        expect(box('B-C window minimum').value).toBe('2.5');
        expect(box('Angle bin width in degrees').value).toBe('2');
        // One computation, of exactly those inputs, through the page's engine.
        expect(tripletRequests()).toHaveLength(1);
        expect(tripletRequests()[0].params).toMatchObject({
            end1: 'Ga', apex: 'Nb', end2: 'Se', r12Min: 2.4, r12Max: 2.9, r23Min: 2.5, r23Max: 2.8, binWidth: 2
        });
        expect(outcome).toMatchObject({ computed: true, result: { angleCount: 64, apexCount: 16 } });
        expect(angleKpi()).toBe('4.0 per Nb');
        expect(controller.state()).toMatchObject({ shown: { end1: 'Ga', r12Min: 2.4, r23Max: 2.8 }, result: { angleCount: 64 } });
    });

    it('a change that leaves the boxes as they are still computes them', async () => {
        await render();
        const controller = getBondGeometry();
        const outcome = await applied(controller, controller.state().inputs);
        expect(outcome.computed).toBe(true);
        expect(tripletRequests()).toHaveLength(1);
    });

    it('Undo\'s form: the earlier boxes back, the result taken off the page, nothing computed', async () => {
        await render();
        const controller = getBondGeometry();
        const before = controller.state().inputs;
        await applied(controller, NEXT);
        const outcome = await applied(controller, before, { compute: false, clear: true });
        expect(outcome).toEqual({ computed: false, reason: null });
        expect(select('End element A').value).toBe('Se');
        expect(box('A-B window maximum').value).toBe('3.00');
        expect(container.querySelector('input[aria-label="Use a distinct B-C window"]').checked).toBe(false);
        expect(tripletRequests()).toHaveLength(1);
        expect(angleKpi()).toBe('—');
        expect(controller.state()).toMatchObject({ inputs: before, shown: null, result: null });
    });

    it('the stamp moves with every edit on the page, the result shown and the run\'s configuration', async () => {
        await render();
        const controller = getBondGeometry();
        const first = controller.stamp();
        expect(controller.stamp()).toBe(first);
        await act(async () => {
            const element = select('End element A');
            element.value = 'Ga';
            element.dispatchEvent(new Event('change', { bubbles: true }));
        });
        const edited = controller.stamp();
        expect(edited).not.toBe(first);
        const compute = [...container.querySelectorAll('button')].find((button) => /Compute/.test(button.textContent));
        await act(async () => { compute.click(); });
        await flush();
        const computed = controller.stamp();
        expect(computed).not.toBe(edited);
        // A Live Data save of the run: a new configuration.
        await render({ dataEpoch: 1 });
        expect(controller.stamp()).not.toBe(computed);
        expect(controller.state().shown).toBeNull();
    });

    it('waits for the run\'s sites before it computes, and a new run settles a change still waiting', async () => {
        let release;
        state.sites = new Promise((resolve) => { release = resolve; });
        await render();
        const controller = getBondGeometry();
        expect(controller.state()).toMatchObject({ ready: false, elements: [] });
        let outcome = null;
        await act(async () => { controller.apply({ ...NEXT, split23: false }).then((value) => { outcome = value; }); });
        await flush();
        expect(outcome).toBeNull();
        expect(tripletRequests()).toHaveLength(0);
        await act(async () => { release(SITES); });
        await flush();
        expect(outcome).toMatchObject({ computed: true });
        expect(tripletRequests()).toHaveLength(1);

    });

    it('a new run settles a change still waiting to compute, as not computed', async () => {
        // The run's sites never arrive: the change waits for them.
        state.sites = new Promise(() => {});
        await render();
        const controller = getBondGeometry();
        let late = null;
        await act(async () => { controller.apply(NEXT).then((value) => { late = value; }); });
        await flush();
        expect(late).toBeNull();
        await render({ directory: 'runs/b' });
        expect(late).toEqual({ computed: false, reason: 'another run opened' });
        expect(tripletRequests()).toHaveLength(0);
    });
});
