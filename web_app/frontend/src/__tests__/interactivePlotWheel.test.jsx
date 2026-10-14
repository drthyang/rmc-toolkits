// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The wheel over a chart: a plain wheel scrolls the page (the chart neither
// zooms nor cancels it); ctrl/⌘ + wheel and a trackpad pinch (ctrl + wheel
// with small pixel deltas) zoom x about the pointer and cancel the browser's
// own page zoom. React's onWheel is passive, where that cancel is ignored, so
// the listener is native; jsdom honours passive listeners as browsers do.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import InteractivePlot from '../components/InteractivePlot';

const angles = Array.from({ length: 181 }, (_, i) => i);
// A fixed, unpadded x range (the Bond Geometry angle axis) and y ticks of at
// most three characters, so the plot area is exactly 642 user units from 60.
const PLOT = {
    title: 'Bond angles',
    xLabel: 'θ (°)',
    yLabel: 'P(θ)',
    xDomain: [0, 180],
    xTicks: [0, 30, 60, 90, 120, 150, 180],
    yMin: 0,
    series: [{ label: 'Ta-Se-Ta', x: angles, y: angles.map((a) => Math.exp(-((a - 90) ** 2) / 200)) }],
};
const LEFT = 60;
const WIDTH = 642;

describe('InteractivePlot wheel', () => {
    let container;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        // jsdom has no SVG geometry: client pixels are user units here.
        vi.stubGlobal('DOMPoint', class {
            constructor(x, y) { this.x = x; this.y = y; }
            matrixTransform() { return { x: this.x, y: this.y }; }
        });
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
        act(() => {
            root.render(<InteractivePlot file={{ path: 'angles', name: 'angles' }} plotData={PLOT} />);
        });
        container.querySelector('svg').getScreenCTM = () => ({ inverse: () => ({}) });
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.unstubAllGlobals();
    });

    // The x range on screen, read back from two labelled x ticks.
    const xRange = () => {
        const ticks = [...container.querySelectorAll('text.plot-tick[text-anchor="middle"]')]
            .map((node) => ({ value: Number(node.textContent), px: Number(node.getAttribute('x')) }));
        const [a, b] = ticks;
        const perUnit = (b.px - a.px) / (b.value - a.value);
        const lo = a.value - (a.px - LEFT) / perUnit;
        return [lo, lo + WIDTH / perUnit];
    };

    // A wheel event at the middle of the plot area; true when it was cancelled.
    const wheel = (init) => {
        const event = new WheelEvent('wheel', {
            bubbles: true, cancelable: true, clientX: LEFT + WIDTH / 2, clientY: 200, ...init,
        });
        act(() => { container.querySelector('svg').dispatchEvent(event); });
        return event.defaultPrevented;
    };
    const resetPill = () => [...container.querySelectorAll('button')].find((b) => b.textContent === 'Reset zoom');

    it('leaves a plain wheel to scroll the page', () => {
        expect(wheel({ deltaY: -100 })).toBe(false);
        expect(wheel({ deltaY: 100, shiftKey: true })).toBe(false);
        expect(resetPill()).toBeUndefined();
        const [lo, hi] = xRange();
        expect(lo).toBeCloseTo(0, 6);
        expect(hi).toBeCloseTo(180, 6);
    });

    it('zooms x about the pointer on ctrl + wheel, a notch at a time', () => {
        expect(wheel({ deltaY: -100, ctrlKey: true })).toBe(true);
        expect(resetPill()).toBeDefined();
        let [lo, hi] = xRange();
        expect(hi - lo).toBeCloseTo(180 / 1.22, 6);
        expect((lo + hi) / 2).toBeCloseTo(90, 6);
        // A larger delta (or a line-mode wheel) is still one 1.22× step, and
        // out undoes in exactly.
        expect(wheel({ deltaY: 3, deltaMode: 1, ctrlKey: true })).toBe(true);
        [lo, hi] = xRange();
        expect(lo).toBeCloseTo(0, 6);
        expect(hi).toBeCloseTo(180, 6);
    });

    it('keeps the value under the pointer in place', () => {
        // A quarter of the way in: 45°, at both zoom steps.
        const clientX = LEFT + WIDTH / 4;
        wheel({ deltaY: -100, ctrlKey: true, clientX });
        let [lo, hi] = xRange();
        expect(hi - lo).toBeCloseTo(180 / 1.22, 6);
        expect(lo + (hi - lo) / 4).toBeCloseTo(45, 6);
        for (let i = 0; i < 10; i += 1) wheel({ deltaY: -2, ctrlKey: true, clientX });
        [lo, hi] = xRange();
        expect(hi - lo).toBeCloseTo(180 / 1.22 / Math.exp(0.2), 6);
        expect(lo + (hi - lo) / 4).toBeCloseTo(45, 6);
    });

    it('zooms on ⌘ + wheel too', () => {
        expect(wheel({ deltaY: -400, metaKey: true })).toBe(true);
        const [lo, hi] = xRange();
        expect(hi - lo).toBeCloseTo(180 / 1.22, 6);
    });

    it('follows a pinch in proportion to its small deltas', () => {
        // Chromium: deltaY = −100·ln(scale), so a 1.05× pinch is ≈ −4.88.
        expect(wheel({ deltaY: -100 * Math.log(1.05), ctrlKey: true })).toBe(true);
        const [lo, hi] = xRange();
        expect(hi - lo).toBeCloseTo(180 / 1.05, 6);
    });

    it('composes events that arrive before the chart re-renders', () => {
        const svg = container.querySelector('svg');
        const init = { bubbles: true, cancelable: true, clientX: LEFT + WIDTH / 2, clientY: 200, ctrlKey: true };
        act(() => {
            for (let i = 0; i < 4; i += 1) svg.dispatchEvent(new WheelEvent('wheel', { ...init, deltaY: -100 * Math.log(1.05) }));
        });
        const [lo, hi] = xRange();
        expect(hi - lo).toBeCloseTo(180 / 1.05 ** 4, 6);
    });

    it('keeps double-click as the reset', () => {
        wheel({ deltaY: -100, ctrlKey: true });
        act(() => {
            container.querySelector('svg').dispatchEvent(new MouseEvent('dblclick', { bubbles: true }));
        });
        expect(resetPill()).toBeUndefined();
        const [lo, hi] = xRange();
        expect(lo).toBeCloseTo(0, 6);
        expect(hi).toBeCloseTo(180, 6);
    });
});
