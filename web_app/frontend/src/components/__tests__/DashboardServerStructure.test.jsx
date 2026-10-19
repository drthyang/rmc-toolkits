// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */
/* global process */

// A server folder (Flask mode): the Dashboard's structure comes from
// /api/structure, and since 0.8 its average-structure basis with it
// (rmc_toolkits/structure_basis.py). So the Detected SG card and its CIF
// button exist for a server folder as for a picked one, the structure the
// page publishes for the Copilot carries the basis (its symmetry and the two
// symmetry tools), and App's configuration epoch (`dataEpoch`: a new .rmc6f
// on disk) re-reads the structure in place -- not the plots -- keeping the
// previous complete read when the new one lands mid-write, and applying only
// the latest request's answer.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import fixture from '../../__tests__/fixtures/structure_basis_fixture.json';

const state = vi.hoisted(() => ({
    structure: null,
    calls: { files: 0, metadata: 0, structure: 0, data: 0 },
    // Set: each /api/structure call waits for its resolver in `pending`
    // (each /api/files call in `pendingFiles` with `deferFiles`).
    defer: false,
    pending: [],
    deferFiles: false,
    pendingFiles: [],
    downloads: [],
}));

// /api/plot/data for the F(Q) fit: what the published run fetches for the
// Copilot (the chart itself is a stub here).
const PLOT_DATA = vi.hoisted(() => ({
    kind: 'xray_sq', title: 'F(Q)', metrics: { rwp: 0.1 }, xLabel: 'Q (Å^{-1})', yLabel: 'F(Q)',
    series: [{ label: 'F(Q)_RMC', x: [1, 2], y: [0.1, 0.2] }, { label: 'F(Q)_Expt', x: [1, 2], y: [0.1, 0.25] }],
}));

vi.mock('axios', () => ({
    default: {
        get: vi.fn(async (url) => {
            if (url.endsWith('/api/plot/data')) {
                state.calls.data += 1;
                return { data: PLOT_DATA };
            }
            if (url.endsWith('/api/files')) {
                state.calls.files += 1;
                const listing = {
                    data: {
                        files: [
                            { name: 'GTS_250K.rmc6f', path: '/runs/a/GTS_250K.rmc6f', type: 'file', plotKind: null },
                            { name: 'GTS_250K_FQ1.csv', path: '/runs/a/GTS_250K_FQ1.csv', type: 'file', plotKind: 'xray_sq' },
                        ],
                    },
                };
                if (state.deferFiles) return new Promise((done) => state.pendingFiles.push(() => done(listing)));
                return listing;
            }
            if (url.endsWith('/api/plot/metadata')) {
                state.calls.metadata += 1;
                return { data: { kind: 'xray_sq', title: 'F(Q)', metrics: {} } };
            }
            if (url.endsWith('/api/structure')) {
                state.calls.structure += 1;
                if (state.defer) return new Promise((done) => state.pending.push((data) => done({ data })));
                return { data: state.structure };
            }
            throw new Error(`unexpected ${url}`);
        }),
    },
}));

// Flask mode (the test environment would otherwise count as static).
vi.mock('../../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => false,
}));

vi.mock('../../figureExport', async (importOriginal) => ({
    ...(await importOriginal()),
    downloadBlob: (blob, filename) => state.downloads.push({ blob, filename }),
}));

vi.mock('../InteractivePlot', () => ({ default: () => <div className="interactive-plot"><svg /></div> }));
vi.mock('../../llm', async (importOriginal) => ({ ...(await importOriginal()), WatchdogBadge: () => null }));

const { default: Dashboard } = await import('../Dashboard');
const { structureFromRmc6f } = await import('../../browserData');
const { copilotRun } = await import('../../copilotRun');
const { clearPlotDataCache } = await import('../../plotDataCache');

// /api/structure for the demo run in a server folder: the model as the parsers
// read it and the basis Flask computes (the golden's demo rows).
const RMC6F = readFileSync(resolve(process.cwd(), 'public/demo/GTS_250K.rmc6f'), 'utf8');
const { basis: _browserBasis, ...model } = structureFromRmc6f({ path: '/runs/a/GTS_250K.rmc6f', text: RMC6F });
const serverStructure = { ...model, basis: fixture.cases[0].basis };
// The same configuration with no reference numbers (a coords-only file): no basis.
const withoutBasis = { ...model, basis: [] };

const flush = () => act(async () => {
    for (let i = 0; i < 8; i += 1) await Promise.resolve();
});

describe('Dashboard: a server folder\'s structure', () => {
    let container;
    let root;
    let published;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        state.structure = serverStructure;
        state.calls = { files: 0, metadata: 0, structure: 0, data: 0 };
        clearPlotDataCache();
        state.defer = false;
        state.pending = [];
        state.deferFiles = false;
        state.pendingFiles = [];
        state.downloads = [];
        published = null;
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
    });

    const render = async (props = {}) => {
        act(() => root.render(
            <Dashboard directory="/runs/a" localRun={null} wantAssistantData onRunChange={(run) => { published = run; }} {...props} />,
        ));
        await flush();
    };
    const card = () => container.querySelector('[aria-label="Detected space group"]');
    const headline = () => card()?.querySelector('dd')?.textContent ?? '';
    const cifButton = () => container.querySelector('.sym-cif-action button');
    const modelCard = () => container.querySelector('[aria-label="Model information"]');

    it('shows the Detected SG card and downloads its CIF', async () => {
        await render();
        expect(modelCard()).not.toBeNull();
        expect(headline()).toMatch(/^F-43m/);
        expect(cifButton().textContent).toBe('⤓CIF · F-43m');
        act(() => cifButton().dispatchEvent(new MouseEvent('click', { bubbles: true })));
        expect(state.downloads).toHaveLength(1);
        expect(state.downloads[0].filename).toBe('GTS_250K_F-43m.cif');
        const text = await state.downloads[0].blob.text();
        expect(text).toMatch(/^data_GTS_250K_F-43m$/m);
        expect(text).toMatch(/_space_group_IT_number\s+216/);
    });

    it('publishes the basis for the Copilot: its symmetry and both symmetry tools', async () => {
        await render();
        expect(published.structure.basis).toBe(serverStructure.basis);
        const run = copilotRun({
            name: '/runs/a', files: published.files, structure: published.structure,
            structurePath: published.structurePath, settings: null, engine: async () => ({}),
        });
        expect(run.props.symmetry).toMatchObject({ spaceGroup: 'F-43m', toleranceA: 0.2 });
        expect(run.tools.map((tool) => tool.name)).toEqual(expect.arrayContaining(['detect_symmetry', 'measure_symmetry_lowering']));
    });

    it('re-reads the structure on a new configuration, not the plots', async () => {
        await render({ dataEpoch: 0 });
        // The plot data the published run carries for the Copilot: once.
        expect(state.calls).toEqual({ files: 1, metadata: 1, structure: 1, data: 1 });

        // Another render with the same epoch reads nothing.
        await render({ dataEpoch: 0, cardMarks: {} });
        expect(state.calls.structure).toBe(1);

        // A new configuration caught mid-write: the previous complete read stays.
        state.structure = { ...withoutBasis, parseReport: { ...model.parseReport, parsedAtoms: 30000 } };
        await render({ dataEpoch: 1 });
        expect(state.calls).toEqual({ files: 1, metadata: 1, structure: 2, data: 1 });
        expect(headline()).toMatch(/^F-43m/);
        expect(modelCard().textContent).toContain('previous read');

        // The complete file: its structure, here with no basis (no card).
        state.structure = withoutBasis;
        await render({ dataEpoch: 2 });
        expect(state.calls).toEqual({ files: 1, metadata: 1, structure: 3, data: 1 });
        expect(card()).toBeNull();
        expect(modelCard()).not.toBeNull();
        expect(modelCard().textContent).not.toContain('previous read');
        expect(published.structure).toBe(withoutBasis);
    });

    it('applies only the latest structure request\'s answer', async () => {
        await render({ dataEpoch: 0 });
        expect(headline()).toMatch(/^F-43m/);
        state.defer = true;
        await render({ dataEpoch: 1 });
        await render({ dataEpoch: 2 });
        expect(state.pending).toHaveLength(2);
        const [older, newer] = state.pending;

        await act(async () => { newer(withoutBasis); });
        await flush();
        expect(card()).toBeNull();
        // The older request answers last: dropped.
        await act(async () => { older(serverStructure); });
        await flush();
        expect(card()).toBeNull();
        expect(published.structure).toBe(withoutBasis);
    });

    it('drops the previous folder\'s structure once another folder starts loading', async () => {
        await render({ dataEpoch: 0 });
        state.defer = true;
        await render({ dataEpoch: 1 });
        const [previousFolder] = state.pending;
        // Another folder: its listing is still on its way when the previous
        // folder's re-read answers -- that answer is dropped.
        state.deferFiles = true;
        await render({ dataEpoch: 1, directory: '/runs/b' });
        await act(async () => { previousFolder(withoutBasis); });
        await flush();
        expect(headline()).toMatch(/^F-43m/);

        await act(async () => { state.pendingFiles[0](); });
        await flush();
        expect(state.pending).toHaveLength(2);
        await act(async () => { state.pending[1](withoutBasis); });
        await flush();
        expect(card()).toBeNull();
        expect(published.structure).toBe(withoutBasis);
    });
});
