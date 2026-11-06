// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The angle histogram's hover hook (InteractivePlot onHoverX, what the Bond
// Geometry page highlights the 3D triplets from): the mouse reports each new
// bin once and null when it leaves, the hovered bin is shaded; a finger's tap
// pins the bin until a tap on it again or anywhere off the chart.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import InteractivePlot from '../components/InteractivePlot';

const centres = Array.from({ length: 180 }, (_, i) => i + 0.5);
const PLOT = {
    title: 'Bond angles',
    xLabel: 'θ (°)',
    yLabel: 'P(θ)',
    xDomain: [0, 180],
    xTicks: [0, 30, 60, 90, 120, 150, 180],
    yMin: 0,
    series: [{
        label: 'sin-corrected',
        x: centres,
        y: centres.map((a) => Math.exp(-((a - 90) ** 2) / 200)),
        curve: 'step',
        binWidth: 1,
        fill: true
    }],
};
// The plot area is 642 user units from 60 (as interactivePlotTouch).
const LEFT = 60;
const WIDTH = 642;
const at = (degrees) => LEFT + (degrees / 180) * WIDTH;

describe('InteractivePlot onHoverX', () => {
    let container;
    let root;
    let svg;
    let onHoverX;

    const mount = (props) => {
        act(() => {
            root.render(<InteractivePlot file={{ path: 'angles', name: 'angles' }} plotData={PLOT} {...props} />);
        });
        svg = container.querySelector('svg');
        svg.getScreenCTM = () => ({ inverse: () => ({}) });
    };

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        vi.stubGlobal('DOMPoint', class {
            constructor(x, y) { this.x = x; this.y = y; }
            matrixTransform() { return { x: this.x, y: this.y }; }
        });
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
        onHoverX = vi.fn();
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.unstubAllGlobals();
    });

    const pointer = (type, pointerType, degrees, target = svg) => {
        act(() => {
            target.dispatchEvent(new PointerEvent(type, {
                bubbles: true, cancelable: true, pointerId: 1, pointerType, clientX: at(degrees), clientY: 200,
            }));
        });
    };
    const tap = (degrees) => {
        pointer('pointerdown', 'touch', degrees);
        pointer('pointerup', 'touch', degrees);
    };
    const reported = () => onHoverX.mock.calls.map(([x]) => x).filter((x, i, all) => i === 0 || x !== all[i - 1]);
    const band = () => container.querySelector('.hover-bin');

    it('reports each bin the mouse enters once, shades it, and null when it leaves', () => {
        mount({ onHoverX });
        pointer('pointermove', 'mouse', 88.2);
        pointer('pointermove', 'mouse', 88.7);
        pointer('pointermove', 'mouse', 92.1);
        expect(onHoverX.mock.calls.filter(([x]) => x !== null).map(([x]) => x)).toEqual([88.5, 92.5]);
        expect(Number(band().getAttribute('x'))).toBeCloseTo(at(92), 6);
        expect(Number(band().getAttribute('width'))).toBeCloseTo(WIDTH / 180, 6);
        pointer('pointerout', 'mouse', 92.1);
        expect(onHoverX).toHaveBeenLastCalledWith(null);
        expect(band()).toBeNull();
    });

    it('pins a tapped bin until a tap on it again', () => {
        mount({ onHoverX });
        tap(60.3);
        expect(onHoverX).toHaveBeenLastCalledWith(60.5);
        // The finger leaving as it lifts keeps the pin.
        pointer('pointerout', 'touch', 60.3);
        expect(onHoverX).toHaveBeenLastCalledWith(60.5);
        expect(band()).not.toBeNull();
        tap(60.6);
        expect(onHoverX).toHaveBeenLastCalledWith(null);
        expect(band()).toBeNull();
        // Another bin moves the pin there.
        tap(60.3);
        tap(120.2);
        expect(reported().slice(-3)).toEqual([60.5, null, 120.5]);
    });

    it('clears a pinned bin on a tap anywhere off the chart', () => {
        mount({ onHoverX });
        tap(45.5);
        expect(onHoverX).toHaveBeenLastCalledWith(45.5);
        pointer('pointerdown', 'touch', 45.5, document.body);
        expect(onHoverX).toHaveBeenLastCalledWith(null);
        expect(container.querySelector('.plot-tooltip')).toBeNull();
    });

    it('reports null when the plot goes', () => {
        mount({ onHoverX });
        pointer('pointermove', 'mouse', 10.1);
        act(() => root.render(<div />));
        expect(onHoverX).toHaveBeenLastCalledWith(null);
    });

    it('leaves a plot without onHoverX as it was: no shaded bin, a second tap re-shows', () => {
        mount({});
        pointer('pointermove', 'mouse', 88.2);
        expect(container.querySelector('.hover-line')).not.toBeNull();
        expect(band()).toBeNull();
        tap(60.3);
        tap(60.3);
        expect(container.querySelector('.plot-tooltip')).not.toBeNull();
    });
});
