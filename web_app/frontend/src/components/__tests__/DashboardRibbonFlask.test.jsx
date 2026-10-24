// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// A server folder (Flask mode): the G(r) panel's ribbon comes from the
// `uncertainty` /api/plot/data sends with the xPDF plot
// (rmc_toolkits/ft_uncertainty.py) — the payload a picked folder gets from
// ftUncertaintyRun.js, so the same files draw the same ribbon in both modes
// (the payload here is the golden's). The chart and the caption read one
// request (plotDataCache), and Live Data asks again when a file the ribbon
// reads changes: the transform, or the F(Q), .dat or PDFpartials beside it.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import fixture from '../../__tests__/fixtures/ft_uncertainty_fixture.json';

const state = vi.hoisted(() => ({ listing: [], dataCalls: [], uncertainty: null, x: [], y: [] }));

vi.mock('axios', () => ({
    default: {
        get: vi.fn(async (url, { params } = {}) => {
            if (url.endsWith('/api/files')) return { data: { files: state.listing } };
            if (url.endsWith('/api/plot/metadata')) return { data: { kind: 'xpdf', title: 'xPDF', metrics: {} } };
            if (url.endsWith('/api/structure')) throw Object.assign(new Error('404'), { response: { data: {} } });
            if (url.endsWith('/api/plot/data')) {
                state.dataCalls.push(params.path);
                return {
                    data: {
                        kind: 'xpdf', title: 'xPDF', metrics: {}, xLabel: 'r (Å)', yLabel: 'G(r)',
                        series: [{ label: 'X_ray-calc', x: state.x, y: state.y.map((v) => 0.97 * v) }, { label: 'X_ray_exp_renorm', x: state.x, y: state.y }],
                        uncertainty: state.uncertainty
                    }
                };
            }
            throw new Error(`unexpected ${url}`);
        }),
    },
}));

vi.mock('../../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => false,
}));
vi.mock('../ModelSummary', () => ({ default: () => null }));
vi.mock('../../llm', () => ({ WatchdogBadge: () => null }));

const { default: Dashboard } = await import('../Dashboard');
const { clearPlotDataCache } = await import('../../plotDataCache');
const { IssueScope, IssueStoreProvider } = await import('../../ui');

const entryOf = (name, modified = 1) => ({ name, path: `/runs/a/${name}`, type: 'file', plotKind: name.includes('_FT_XFQ') ? 'xpdf' : null, modified, size: 100 });

const flush = () => act(async () => {
    for (let i = 0; i < 20; i += 1) await new Promise((resolve) => setTimeout(resolve, 0));
});

describe('Dashboard: a server folder\'s G(r) ribbon', () => {
    let container;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        clearPlotDataCache();
        const propagated = fixture.cases.find((entry) => entry.name === 'propagated, CONVOLVE');
        const transformText = fixture.texts[propagated.files['rmc_run_FT_XFQ1.csv']];
        const rows = transformText.trim().split('\n').slice(1).map((row) => row.split(',').map(Number));
        state.x = rows.map((row) => row[0]);
        state.y = rows.map((row) => row[2]);
        state.uncertainty = propagated.expected;
        state.dataCalls = [];
        state.listing = ['rmc_run_FT_XFQ1.csv', 'rmc_run_XFQ1.csv', 'rmc_run.dat', 'rmc_run_PDFpartials.csv'].map((name) => entryOf(name));
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
    });

    const mount = async (watchFiles = false) => {
        act(() => root.render(
            <IssueStoreProvider>
                <IssueScope page="dashboard">
                    <Dashboard directory="/runs/a" localRun={null} watchFiles={watchFiles} />
                </IssueScope>
            </IssueStoreProvider>
        ));
        await flush();
    };
    const card = () => container.querySelector('[data-copilot-card="plot:rmc_run_FT_XFQ1.csv"]');

    it('draws the server\'s propagated ribbon, one request for the chart and the caption', async () => {
        await mount();
        expect(card().querySelector('path.series-ribbon')).not.toBeNull();
        expect([...card().querySelectorAll('.plot-legend button')].map((button) => button.textContent))
            .toEqual(['X_ray_exp_renorm', 'X_ray-calc', '±1σ propagated']);
        expect(card().querySelector('.gr-ribbon-note').textContent).toContain('the r points are correlated');
        expect(state.dataCalls).toEqual(['/runs/a/rmc_run_FT_XFQ1.csv']);
    });

    it('a server payload without σ: the caption says why', async () => {
        state.uncertainty = fixture.demo.expected;
        await mount();
        expect(card().querySelector('path.series-ribbon')).toBeNull();
        expect(card().querySelector('.gr-ribbon-note').textContent).toContain(
            'No error ribbon. F(Q) carries no uncertainties: scale_ft_rmc.fq (named in GTS_250K.dat) is not in the run folder.'
        );
    });

    it('Live Data asks for the ribbon again when the F(Q) beside the transform changes', async () => {
        vi.useFakeTimers({ shouldAdvanceTime: true });
        try {
            await mount(true);
            expect(state.dataCalls).toHaveLength(1);
            // RMCProfile rewrote the run's F(Q) (not the transform): the ribbon reads it.
            state.listing = state.listing.map((entry) => (entry.name === 'rmc_run_XFQ1.csv' ? { ...entry, modified: 2 } : entry));
            state.uncertainty = { ...state.uncertainty, reason: null };
            await act(async () => { vi.advanceTimersByTime(3100); });
            await flush();
            expect(state.dataCalls).toHaveLength(2);
        } finally {
            vi.useRealTimers();
        }
    });
});
