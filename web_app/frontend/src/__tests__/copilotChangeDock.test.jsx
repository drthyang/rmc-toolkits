// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The change end to end: the dock (CopilotChat, with the host's change tool
// as App builds it, copilotHostTools.js) beside the Bond Geometry page. A
// scripted model asks set_bond_geometry for a new triplet and window; the
// card waits, and nothing on the page moves until Apply; then the page's own
// controls show the change and the page computes it, the model reads what the
// page holds, and Undo on the card puts the page back. Flask mode for the
// page (axios mocked); the model behind a mocked OpenAI-compatible server.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import CopilotChat from '../llm/components/CopilotChat';
import { DOCK_STORAGE_KEY } from '../llm/components/dockLayout';
import { forgetServers } from '../llm/provider/contextGuard';
import { saveSettings } from '../llm/settings';
import { createChangeHistory } from '../bondGeometryControl';
import { copilotHostTools } from '../copilotHostTools';

const state = vi.hoisted(() => ({ requests: [] }));

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
    bond23: [Number(params.r12Min), Number(params.r12Max)],
    sharedEnds: params.end1 === params.end2,
    binCenters: [0.5, 1.5],
    binWidth: Number(params.binWidth),
    sinCorrected: [0, 1],
    density: [0, 1],
    coordination: [0, 0, 1],
    coordinationSummary: { mean: 2, modal: 2, modalShare: 1, unbondedShare: 0 },
    apexCount: 16,
    lengths12: { uniqueBonds: 96, meanLength: 2.6, count: 96 },
    lengths23: { uniqueBonds: 96, meanLength: 2.6, count: 96 },
    angleCount: 64,
    meanAngle: 90,
    stdAngle: 5,
});

vi.mock('axios', () => ({
    default: {
        get: vi.fn(async (url, config) => {
            state.requests.push({ url, params: config?.params });
            if (url.endsWith('/api/pca/sites')) return { data: SITES };
            if (url.endsWith('/api/triplets')) return { data: triplets(config.params) };
            if (url.endsWith('/api/files')) return { data: { files: [] } };
            if (url.endsWith('/api/structure')) return { data: { atomCount: 52, elements: SITES.elements } };
            throw new Error(`unexpected ${url}`);
        }),
    },
}));
vi.mock('../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => false,
}));
vi.mock('../components/FoldedCellPanel', () => ({ default: () => null }));
vi.mock('../components/InteractivePlot', () => ({ default: () => null }));
vi.mock('../components/ModelSummary', () => ({ default: () => null }));

const { default: BondGeometryPage } = await import('../components/BondGeometryPage');

const sse = (chunks) => {
    const encoder = new TextEncoder();
    const body = new ReadableStream({
        start(controller) {
            chunks.forEach((chunk) => controller.enqueue(encoder.encode(`data: ${JSON.stringify(chunk)}\n\n`)));
            controller.enqueue(encoder.encode('data: [DONE]\n\n'));
            controller.close();
        }
    });
    return { ok: true, status: 200, body };
};
const delta = (fields, finish = null) => ({ choices: [{ delta: fields, finish_reason: finish }] });

// Real time passes (the tool polls the page every 50 ms) while React renders.
const wait = async (ready, ms = 4000) => {
    const end = Date.now() + ms;
    while (!ready() && Date.now() < end) {
        await act(async () => { await new Promise((resolve) => setTimeout(resolve, 20)); });
    }
};

const CARDS = [{ id: 'bond-angles', page: 'geometry', title: 'Bond-angle distribution', about: 'the angle histogram' }];

describe('a Copilot change on the Bond Geometry page, through the dock', () => {
    let container;
    let root;
    let chatBodies;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        state.requests = [];
        forgetServers();
        localStorage.removeItem(DOCK_STORAGE_KEY);
        saveSettings({ baseUrl: 'http://localhost:11434/v1', model: '', apiKey: '', copilotTools: true });
        chatBodies = [];
        const replies = [
            sse([
                delta({ content: 'Setting up Ga–Nb–Se on the page.' }),
                delta({ tool_calls: [{ index: 0, id: 'call_s', type: 'function', function: { name: 'set_bond_geometry', arguments: '{"end1":"Ga","r12_min":2.4,"r12_max":2.9}' } }] }, 'tool_calls')
            ]),
            sse([delta({ content: 'The page now shows Ga–Nb–Se in 2.40–2.90 Å.\n\nOutcome: achieved — set on the page.' }, 'stop')])
        ];
        vi.stubGlobal('fetch', vi.fn(async (url, options = {}) => {
            if (/\/api\//.test(String(url))) return { ok: false, status: 404, json: async () => ({ error: 'not found' }) };
            if (String(url).endsWith('/v1/models')) return { ok: true, status: 200, json: async () => ({ data: [{ id: 'tool-model' }] }) };
            chatBodies.push(JSON.parse(options.body));
            return replies.shift();
        }));
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.unstubAllGlobals();
    });

    it('waits for Apply, then sets the page\'s controls, the page computes, the model reads it, and Undo puts it back', async () => {
        const history = createChangeHistory();
        const tools = copilotHostTools({ elements: ['Ga', 'Nb', 'Se'], cards: CARDS, pages: { geometry: 'Bond Geometry' }, history });
        await act(async () => {
            root.render(
                <>
                    <BondGeometryPage directory="runs/a" localRun={null} dataEpoch={0} />
                    <CopilotChat
                        runName="demo"
                        plotFiles={[{ plotKind: 'neutron_sq', plotData: { title: 'S(Q)', metrics: { rwp: 0.04 } } }]}
                        tools={tools}
                        cards={CARDS}
                        pages={{ geometry: 'Bond Geometry' }}
                        onOpenPage={() => {}}
                        runKey="run-1"
                    />
                </>
            );
        });
        await wait(() => container.querySelector('select[aria-label="End element A"]')?.value === 'Se');
        const select = (label) => container.querySelector(`select[aria-label="${label}"]`);
        const box = (label) => container.querySelector(`input[aria-label="${label}"]`);
        const triplets = () => state.requests.filter((request) => request.url.endsWith('/api/triplets'));
        const card = () => container.querySelector('.llm-change');

        act(() => container.querySelector('.llm-chat-launcher').click());
        await wait(() => !container.querySelector('.llm-composer-input')?.disabled);
        const textarea = container.querySelector('.llm-composer-input');
        act(() => {
            Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype, 'value').set.call(textarea, 'Set up the Ga–Nb–Se angles on the page.');
            textarea.dispatchEvent(new Event('input', { bubbles: true }));
        });
        act(() => container.querySelector('.llm-composer').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true })));

        // The card waits; the page has not moved.
        await wait(() => card()?.className.includes('is-waiting'));
        expect(chatBodies[0].tools.map((tool) => tool.function.name)).toContain('set_bond_geometry');
        expect(card().querySelector('.llm-change-text').textContent).toBe('Ga–Nb–Se · 2.40–2.90 Å · 1° bins');
        expect(select('End element A').value).toBe('Se');
        expect(box('A-B window maximum').value).toBe('3.00');
        expect(triplets()).toHaveLength(0);
        // The card the change is about breathes while it waits.
        expect(container.querySelector('[data-copilot-card="bond-angles"]').getAttribute('data-copilot-focus')).toBeTruthy();

        act(() => [...card().querySelectorAll('button')].find((button) => button.textContent === 'Apply').click());
        await wait(() => container.querySelector('.llm-outcome'));

        // The page's own controls hold the change, and it computed exactly that.
        expect(select('End element A').value).toBe('Ga');
        expect(box('A-B window minimum').value).toBe('2.4');
        expect(box('A-B window maximum').value).toBe('2.9');
        expect(triplets()).toHaveLength(1);
        expect(triplets()[0].params).toMatchObject({ end1: 'Ga', apex: 'Nb', end2: 'Se', r12Min: 2.4, r12Max: 2.9, binWidth: 1 });
        expect(container.querySelector('[aria-label="Triplet result"] .ui-kpi__value').textContent).toBe('4.0 per Nb');
        // The model read what the page holds.
        const message = JSON.parse(chatBodies[1].messages.at(-1).content);
        expect(message.data).toEqual({
            applied: true,
            inputs: { end1: 'Ga', apex: 'Nb', end2: 'Se', r12_min: 2.4, r12_max: 2.9, distinct_bc: false, bin_width: 1 },
            computed: true,
            angle_count: 64
        });
        expect(card().className).toContain('is-applied');
        expect(history.entries()).toEqual([expect.objectContaining({ actor: 'agent', undone: false })]);

        // Undo on the card: the earlier inputs back, the agent's result off the page.
        await act(async () => { [...card().querySelectorAll('button')].find((button) => button.textContent === 'Undo').click(); });
        await wait(() => card().textContent.includes('Undone'));
        expect(select('End element A').value).toBe('Se');
        expect(box('A-B window maximum').value).toBe('3.00');
        expect(container.querySelector('[aria-label="Triplet result"] .ui-kpi__value').textContent).toBe('—');
        expect(triplets()).toHaveLength(1);
        expect(history.entries()[0].undone).toBe(true);
    });
});
