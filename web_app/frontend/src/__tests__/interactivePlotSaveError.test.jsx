// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// A failed figure save must say so. The chart is on screen when the user
// saves it, so the error cannot take the "no plot yet" error branch (which
// replaces the whole chart): it is a one-line danger banner under the plot
// toolbar, cleared by the next save.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';

const state = vi.hoisted(() => ({ outcomes: [] }));

vi.mock('../figureExport', () => ({
    saveSvgFigure: vi.fn(async () => {
        const outcome = state.outcomes.shift();
        if (outcome) throw outcome;
    }),
}));

const { default: InteractivePlot } = await import('../components/InteractivePlot');

const PLOT = {
    title: 'F(Q)',
    xLabel: 'Q',
    yLabel: 'F(Q)',
    series: [{ label: 'F(Q)_RMC', x: [1, 2, 3], y: [0.1, 0.2, 0.3] }],
};

describe('InteractivePlot save errors', () => {
    let container;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        state.outcomes = [];
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
    });

    const render = () => act(() => {
        root.render(<InteractivePlot file={{ path: 'run_FQ1.csv', name: 'run_FQ1.csv' }} plotData={PLOT} />);
    });

    const save = async () => {
        act(() => container.querySelector('.ui-save__trigger').click());
        await act(async () => {
            container.querySelector('[role="menuitem"]').click();
        });
    };

    const banner = () => container.querySelector('.ui-banner--danger');

    it('shows the failure under the toolbar and keeps the chart', async () => {
        render();
        state.outcomes = [new Error('Could not rasterize the figure')];
        await save();
        expect(banner()?.textContent).toContain('Could not rasterize the figure');
        expect(banner()?.getAttribute('role')).toBe('alert');
        // The chart stays: the banner sits between the toolbar and the stage.
        expect(container.querySelector('.interactive-plot svg')).not.toBeNull();
        expect(banner().previousElementSibling?.classList.contains('plot-toolbar')).toBe(true);
    });

    it('falls back to a generic line when the error has no message', async () => {
        render();
        state.outcomes = [new Error('')];
        await save();
        expect(banner()?.textContent).toContain('Could not save the figure');
    });

    it('clears the banner on the next save', async () => {
        render();
        state.outcomes = [new Error('Could not encode the figure')];
        await save();
        expect(banner()).not.toBeNull();
        await save();
        expect(banner()).toBeNull();
    });
});
