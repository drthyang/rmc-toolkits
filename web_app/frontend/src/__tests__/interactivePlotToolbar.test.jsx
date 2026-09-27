// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The plot toolbar's opt-in props: `legend={false}` drops the series legend
// (a card whose title already names the curves), and `actionsTarget` moves
// the plot actions (Reset zoom, Save) into an element the caller owns — a
// card header — so a short card gives the toolbar row's height to the plot.
// Without them the toolbar renders as before (interactivePlotAxes.test pins
// that markup).

import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import InteractivePlot from '../components/InteractivePlot';

const plotData = {
    title: 'Ta-Se partial g(r)',
    xLabel: 'r (Å)',
    yLabel: 'g(r)',
    series: [
        { label: 'Ta-Se', x: [1, 2, 3, 4], y: [0, 1, 4, 1] },
        { label: 'rmin 2.0', x: [2, 2], y: [0, 4], role: 'guide', legend: false },
    ],
};

describe('InteractivePlot toolbar options', () => {
    let container;
    let target;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        container = document.createElement('div');
        target = document.createElement('span');
        document.body.append(container, target);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        target.remove();
    });

    const render = (props) => act(() => {
        root.render(<InteractivePlot file={{ path: 'p', name: 'p' }} plotData={plotData} {...props} />);
    });
    const legendLabels = () => [...container.querySelectorAll('.plot-legend button')].map((button) => button.textContent);

    it('by default the toolbar row holds the legend and the actions', () => {
        render({});
        expect(legendLabels()).toEqual(['Ta-Se']);
        expect(container.querySelector('.plot-toolbar .plot-actions .ui-save__trigger')).not.toBeNull();
    });

    it('legend={false} drops the legend and keeps the actions on the right', () => {
        render({ legend: false });
        expect(legendLabels()).toEqual([]);
        const toolbar = container.querySelector('.plot-toolbar');
        expect([...toolbar.children].map((child) => child.className)).toEqual(['plot-legend', 'plot-actions']);
        expect(toolbar.querySelector('.ui-save__trigger')).not.toBeNull();
    });

    it('with actionsTarget the actions render there and the plot has no toolbar row', () => {
        render({ legend: false, actionsTarget: target });
        expect(container.querySelector('.plot-toolbar')).toBeNull();
        expect(container.querySelector('.ui-save__trigger')).toBeNull();
        expect(target.querySelector('.plot-actions .ui-save__trigger')).not.toBeNull();
        // The svg is still the plot's own.
        expect(container.querySelector('svg[role="img"]').getAttribute('aria-label')).toBe('Ta-Se partial g(r)');
    });

    it('a null actionsTarget (not mounted yet) renders the actions nowhere, not inline', () => {
        render({ legend: false, actionsTarget: null });
        expect(container.querySelector('.plot-toolbar')).toBeNull();
        expect(document.querySelector('.ui-save__trigger')).toBeNull();
    });
});
