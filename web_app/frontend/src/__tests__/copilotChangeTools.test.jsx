// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The Copilot's change tool, set_bond_geometry (copilotChangeTools.js), on
// the Bond Geometry page: its strict schema; refused in plain words when no
// run is open, when the page cannot read the run and when the values leave
// the page unable to compute; the approval (Ask first and Auto) and a
// decline; the stale-state check (the page moved while the card waited: not
// applied); an approved change set through the page's own controls, and the
// page computing it; Undo putting the earlier inputs back. Each call goes
// through executeToolCall, as the Copilot runs it. Flask mode, axios mocked.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import { checkShape } from '../llm/copilot/dataShape';
import { executeToolCall, toRequestTools } from '../llm/copilot/tools';
import { createChangeHistory } from '../bondGeometryControl';
import { SET_BOND_GEOMETRY_SHAPE, buildChangeTools } from '../copilotChangeTools';

const state = vi.hoisted(() => ({ requests: [], sites: null, staticMode: false }));

const SITES = {
    elements: ['Ga', 'Nb', 'Se'],
    sites: [
        { referenceNumber: 1, element: 'Ga', count: 4, copiesPerCell: 4 },
        { referenceNumber: 2, element: 'Nb', count: 16, copiesPerCell: 16 },
        { referenceNumber: 3, element: 'Se', count: 32, copiesPerCell: 32 },
    ],
};

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
    angleCount: Number(params.r12Max) < 2.2 ? 0 : 240,
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
vi.mock('../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => state.staticMode,
}));
vi.mock('../components/FoldedCellPanel', () => ({ default: () => null }));
vi.mock('../components/InteractivePlot', () => ({ default: () => null }));
vi.mock('../components/ModelSummary', () => ({ default: () => null }));

const { default: BondGeometryPage } = await import('../components/BondGeometryPage');

const settle = () => act(async () => {
    for (let i = 0; i < 8; i += 1) await new Promise((resolve) => setTimeout(resolve, 0));
});
// The tool's polls run on timers: let them, and the page's renders, go on.
const settleFor = async (ms) => {
    const end = Date.now() + ms;
    while (Date.now() < end) {
        await act(async () => { await new Promise((resolve) => setTimeout(resolve, 20)); });
    }
};

const ELEMENTS = ['Ga', 'Nb', 'Se'];
const WAITS = { page: 300, ready: 600, compute: 3000 };
const approving = (auto = false) => vi.fn(async () => ({ approved: true, auto }));

describe('set_bond_geometry on the Bond Geometry page', () => {
    let container;
    let root;
    let history;
    let tool;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        state.requests = [];
        state.sites = null;
        state.staticMode = false;
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
        history = createChangeHistory();
        [tool] = buildChangeTools({ elements: ELEMENTS, history, waits: WAITS });
    });

    afterEach(() => {
        act(() => root?.unmount());
        container.remove();
    });

    const render = async (props = {}) => {
        await act(async () => {
            root.render(<BondGeometryPage directory="runs/a" localRun={null} dataEpoch={0} {...props} />);
        });
        await settle();
    };
    // One call as the Copilot makes it: parsed, validated against the schema,
    // run with the approve port, its result as the contract decides it.
    const call = async (args, approve = approving()) => {
        let promise;
        await act(async () => { promise = executeToolCall(tool, JSON.stringify(args), { timeoutMs: null, context: { approve } }); });
        let result;
        promise.then((value) => { result = value; });
        for (let i = 0; i < 200 && !result; i += 1) await settleFor(20);
        return result;
    };
    const select = (label) => container.querySelector(`select[aria-label="${label}"]`);
    const box = (label) => container.querySelector(`input[aria-label="${label}"]`);
    const tripletRequests = () => state.requests.filter((request) => request.url.endsWith('/api/triplets'));
    const setSelect = async (label, value) => {
        await act(async () => {
            const element = select(label);
            element.value = value;
            element.dispatchEvent(new Event('change', { bubbles: true }));
        });
    };

    it('its schema mirrors the page\'s controls, strict, with the boxes\' ranges', () => {
        const [sent] = toRequestTools([tool]);
        expect(sent.function.name).toBe('set_bond_geometry');
        expect(sent.function.parameters).toEqual({
            type: 'object',
            properties: {
                end1: { type: 'string', enum: ELEMENTS, description: 'Element A.' },
                apex: { type: 'string', enum: ELEMENTS, description: 'Element B, the central atom.' },
                end2: { type: 'string', enum: ELEMENTS, description: 'Element C.' },
                r12_min: { type: 'number', minimum: 0, maximum: 15, description: 'A–B window minimum (Å).' },
                r12_max: { type: 'number', minimum: 0, maximum: 15, description: 'A–B window maximum (Å).' },
                distinct_bc: {
                    type: 'boolean',
                    description: 'The Distinct B–C switch: false makes the B–C bond use the A–B window again; r23_min with r23_max turns it on.'
                },
                r23_min: { type: 'number', minimum: 0, maximum: 15, description: 'B–C window minimum (Å); give both or neither.' },
                r23_max: { type: 'number', minimum: 0, maximum: 15, description: 'B–C window maximum (Å).' },
                bin_width: { type: 'number', minimum: 0.1, maximum: 45, description: 'Angle bin width (degrees).' }
            },
            additionalProperties: false
        });
        expect(tool.description).toMatch(/approves the change on a card first \(unless they chose Auto\)/);
        expect(tool.method).toMatch(/coordination-shells skill/);
        expect(tool.method).toMatch(/first_shell\.window_A/);
        // No run, no tool: the run's elements are what it sets.
        expect(buildChangeTools({ elements: [], history })).toEqual([]);
        expect(() => buildChangeTools({ elements: ELEMENTS })).toThrow(/history/);
    });

    it('the schema refuses before anything runs: an unknown argument, a value out of range, an element the run lacks', async () => {
        const approve = approving();
        await render();
        for (const [args, error] of [
            [{ window: [2, 3] }, 'unknown argument "window"; the arguments are end1, apex, end2, r12_min, r12_max, distinct_bc, r23_min, r23_max, bin_width'],
            [{ r12_max: 16 }, 'r12_max must be ≤ 15, got 16'],
            [{ r12_min: -1 }, 'r12_min must be ≥ 0, got -1'],
            [{ bin_width: 0.05 }, 'bin_width must be ≥ 0.1, got 0.05'],
            [{ end1: 'Ta' }, 'end1 must be one of "Ga", "Nb", "Se", got "Ta"'],
            [{ distinct_bc: 'yes' }, 'distinct_bc must be true or false, got "yes"']
        ]) {
            const result = await call(args, approve);
            expect(result, JSON.stringify(args)).toMatchObject({ resultStatus: 'error', error });
        }
        expect(approve).not.toHaveBeenCalled();
        expect(tripletRequests()).toHaveLength(0);
    });

    it('is refused in plain words: no page, no run open, a run the page cannot read, values the page cannot compute', async () => {
        const approve = approving();
        // The page never mounts.
        expect(await call({ r12_max: 2.8 }, approve)).toMatchObject({ resultStatus: 'error', error: 'the Bond Geometry page did not open' });
        // The static build with no run open.
        state.staticMode = true;
        await render();
        expect(await call({ r12_max: 2.8 }, approve)).toMatchObject({ resultStatus: 'error', error: 'no run is open: open a run folder, or load the Demo' });
        // The run's sites fail to load.
        state.staticMode = false;
        state.sites = Promise.reject(new Error('Request failed with status code 404'));
        state.sites.catch(() => {});
        await render({ directory: 'runs/broken' });
        expect(await call({ r12_max: 2.8 }, approve)).toMatchObject({
            resultStatus: 'error', error: 'the Bond Geometry page could not read the run: Request failed with status code 404'
        });
        // Values that contradict the page's own boxes, or each other.
        state.sites = null;
        await render({ directory: 'runs/a' });
        expect(await call({ r12_max: 1.5 }, approve)).toMatchObject({
            resultStatus: 'error',
            error: 'the page cannot compute that: the A–B window would be 2.00–1.50 Å: its minimum must be below its maximum'
        });
        expect(await call({ r23_min: 2.4 }, approve)).toMatchObject({ resultStatus: 'error', error: 'give r23_min and r23_max together, or neither' });
        expect(await call({ distinct_bc: true }, approve)).toMatchObject({ resultStatus: 'error', error: 'give r23_min and r23_max to turn the distinct B–C window on' });
        expect(approve).not.toHaveBeenCalled();
        expect(tripletRequests()).toHaveLength(0);
    });

    it('asks first, then sets the page\'s controls and the page computes; the change is the agent\'s, in the history', async () => {
        await render();
        const approve = approving();
        const result = await call({ end1: 'Ga', r12_min: 2.4, r12_max: 2.9, r23_min: 2.5, r23_max: 2.8, bin_width: 2 }, approve);
        expect(approve).toHaveBeenCalledTimes(1);
        expect(approve.mock.calls[0][0]).toEqual({
            tool: 'set_bond_geometry',
            title: 'Change the Bond Geometry inputs',
            line: 'Ga–Nb–Se · A–B 2.40–2.90 Å, B–C 2.50–2.80 Å · 2° bins',
            details: [
                'Triplet: Se–Nb–Se → Ga–Nb–Se',
                'A–B window: 2.00–3.00 → 2.40–2.90 Å',
                'Distinct B–C window: off → on',
                'B–C window: 2.00–3.00 → 2.50–2.80 Å',
                'Bin width: 1° → 2°',
                'Then the page computes the angle distribution; Undo puts the inputs back.'
            ]
        });
        expect(result.resultStatus).toBe('ok');
        expect(result.data).toEqual({
            applied: true,
            inputs: { end1: 'Ga', apex: 'Nb', end2: 'Se', r12_min: 2.4, r12_max: 2.9, distinct_bc: true, r23_min: 2.5, r23_max: 2.8, bin_width: 2 },
            computed: true,
            angle_count: 240
        });
        expect(checkShape(result.data, SET_BOND_GEOMETRY_SHAPE)).toEqual([]);
        expect(result.checks).toEqual([{ level: 'ok', text: 'Applied: Ga–Nb–Se · A–B 2.40–2.90 Å, B–C 2.50–2.80 Å · 2° bins; the page computed 240 angles.' }]);
        expect(result.change).toMatchObject({ status: 'applied', auto: false, id: 'bond-geometry-1', line: 'Ga–Nb–Se · A–B 2.40–2.90 Å, B–C 2.50–2.80 Å · 2° bins' });
        // The page's own controls show it, and it computed exactly that.
        expect(select('End element A').value).toBe('Ga');
        expect(box('A-B window maximum').value).toBe('2.9');
        expect(box('B-C window minimum').value).toBe('2.5');
        expect(tripletRequests().map((request) => request.params)).toEqual([
            { dir: 'runs/a', end1: 'Ga', apex: 'Nb', end2: 'Se', r12Min: 2.4, r12Max: 2.9, r23Min: 2.5, r23Max: 2.8, binWidth: 2 }
        ]);
        expect(history.entries()).toEqual([expect.objectContaining({
            id: 'bond-geometry-1', actor: 'agent', shownBefore: false, undone: false,
            before: expect.objectContaining({ end1: 'Se', r12Max: '3.00', split23: false }),
            after: expect.objectContaining({ end1: 'Ga', r12Max: '2.9', split23: true })
        })]);
    });

    it('in Auto the change is applied at once and says so; with no angles in the windows it warns', async () => {
        await render();
        const result = await call({ r12_min: 1.9, r12_max: 2.1 }, approving(true));
        expect(result.change).toMatchObject({ status: 'applied', auto: true });
        expect(result.data).toMatchObject({ applied: true, computed: true, angle_count: 0 });
        expect(result.checks).toEqual([{ level: 'warn', text: 'Applied: Se–Nb–Se · 1.90–2.10 Å · 1° bins; the page found no angles in these windows: widen them.' }]);
    });

    it('a declined change returns {"declined": true} and leaves the page as it was', async () => {
        await render();
        const result = await call({ end1: 'Ga' }, vi.fn(async () => ({ approved: false, auto: false })));
        expect(result.resultStatus).toBe('ok');
        expect(result.data).toEqual({ declined: true });
        expect(result.checks).toEqual([{ level: 'note', text: 'The user declined this change. Ask what they would prefer, or propose something else.' }]);
        expect(result.change).toMatchObject({ status: 'declined', id: null });
        expect(select('End element A').value).toBe('Se');
        expect(tripletRequests()).toHaveLength(0);
        expect(history.entries()).toEqual([]);
    });

    it('not applied when the page moved while the card waited (the stale-state check)', async () => {
        await render();
        let decide;
        const approve = vi.fn(() => new Promise((resolve) => { decide = resolve; }));
        let promise;
        await act(async () => { promise = executeToolCall(tool, JSON.stringify({ r12_max: 2.8 }), { timeoutMs: null, context: { approve } }); });
        for (let i = 0; i < 50 && !decide; i += 1) await settleFor(20);
        expect(approve).toHaveBeenCalledTimes(1);
        // The user picks another A while the card waits, then approves.
        await setSelect('End element A', 'Ga');
        decide({ approved: true, auto: false });
        let result;
        promise.then((value) => { result = value; });
        for (let i = 0; i < 50 && !result; i += 1) await settleFor(20);
        expect(result.resultStatus).toBe('refused');
        expect(result.data).toEqual({ stale: true, inputs: { end1: 'Ga', apex: 'Nb', end2: 'Se', r12_min: 2, r12_max: 3, distinct_bc: false, bin_width: 1 } });
        expect(result.checks[0]).toMatchObject({ level: 'fail', text: expect.stringMatching(/^Not applied: the Bond Geometry page changed while the card waited/) });
        expect(result.change).toMatchObject({ status: 'stale', id: null });
        // The user's edit stands, nothing was computed, nothing recorded.
        expect(select('End element A').value).toBe('Ga');
        expect(box('A-B window maximum').value).toBe('3.00');
        expect(tripletRequests()).toHaveLength(0);
        expect(history.entries()).toEqual([]);
    });

    it('asks nothing when the page already shows exactly that', async () => {
        await render();
        await call({ r12_max: 2.8 }, approving(true));
        const approve = approving();
        const result = await call({ r12_max: 2.8 }, approve);
        expect(approve).not.toHaveBeenCalled();
        expect(result.data).toEqual({ unchanged: true, inputs: expect.objectContaining({ r12_max: 2.8 }) });
        expect(result.change).toMatchObject({ status: 'unchanged' });
        expect(tripletRequests()).toHaveLength(1);
    });

    it('Undo puts the earlier inputs back: the result taken off a page that showed none, computed again on one that did', async () => {
        await render();
        const first = await call({ end1: 'Ga', r12_max: 2.9 }, approving());
        let undone;
        await act(async () => { undone = tool.undo(first.change); });
        await settleFor(100);
        expect(await undone).toEqual({ ok: true, computed: false });
        expect(select('End element A').value).toBe('Se');
        expect(box('A-B window maximum').value).toBe('3.00');
        expect(container.querySelector('[aria-label="Triplet result"] .ui-kpi__value').textContent).toBe('—');
        expect(history.entries()[0].undone).toBe(true);
        expect(await tool.undo(first.change)).toEqual({ ok: false, reason: 'it was undone already' });

        // Now the page shows a result (the user's Compute), and the agent changes it.
        const compute = [...container.querySelectorAll('button')].find((button) => /Compute|Update/.test(button.textContent));
        await act(async () => { compute.click(); });
        await settle();
        const second = await call({ r12_max: 2.7 }, approving());
        const before = tripletRequests().length;
        await act(async () => { undone = tool.undo(second.change); });
        await settleFor(100);
        expect(await undone).toEqual({ ok: true, computed: true });
        expect(box('A-B window maximum').value).toBe('3.00');
        expect(tripletRequests()).toHaveLength(before + 1);
        expect(tripletRequests().at(-1).params).toMatchObject({ end1: 'Se', r12Max: 3 });
    });

    it('Undo is refused once the user changed the inputs since, or on another run', async () => {
        await render();
        const applied = await call({ r12_max: 2.9 }, approving());
        await setSelect('End element A', 'Ga');
        expect(await tool.undo(applied.change)).toEqual({ ok: false, reason: 'the inputs changed on the page since' });
        await setSelect('End element A', 'Se');
        await render({ directory: 'runs/b' });
        expect(await tool.undo(applied.change)).toEqual({ ok: false, reason: 'another run or configuration is open' });
        expect(await tool.undo({ id: 'bond-geometry-9' })).toEqual({ ok: false, reason: 'this change is not in the history' });
        expect(history.entries()[0].undone).toBe(false);
    });
});
