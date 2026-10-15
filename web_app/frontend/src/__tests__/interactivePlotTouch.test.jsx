// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */
/* global process */

// A finger on a chart: the chart's touch-action leaves a vertical drag (the
// page scroll) and a pinch (the page zoom) to the browser, which cancels the
// pointer when it takes the gesture. A drag that starts sideways reaches the
// chart and zooms x only; a tap shows the readout touch has no hover for. The
// mouse keeps its box zoom and hover.

import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import InteractivePlot from '../components/InteractivePlot';
import { svgToSvgString } from '../figureExport';

// Read from disk: Vitest stubs CSS imports, ?raw included.
const plotCss = readFileSync(resolve(process.cwd(), 'src/components/InteractivePlot.css'), 'utf8');

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
const at = (degrees) => LEFT + (degrees / 180) * WIDTH;

describe('InteractivePlot touch', () => {
    let container;
    let root;
    let svg;

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
        svg = container.querySelector('svg');
        svg.getScreenCTM = () => ({ inverse: () => ({}) });
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.unstubAllGlobals();
    });

    // The x range on screen, read back from two labelled x ticks and the plot
    // area (the clip rectangle: a y zoom's longer tick labels move its left).
    const xRange = () => {
        const clip = container.querySelector('clipPath rect');
        const left = Number(clip.getAttribute('x'));
        const width = Number(clip.getAttribute('width'));
        const ticks = [...container.querySelectorAll('text.plot-tick[text-anchor="middle"]')]
            .map((node) => ({ value: Number(node.textContent), px: Number(node.getAttribute('x')) }));
        const [a, b] = ticks;
        const perUnit = (b.px - a.px) / (b.value - a.value);
        const lo = a.value - (a.px - left) / perUnit;
        return [lo, lo + width / perUnit];
    };
    const yTicks = () => [...container.querySelectorAll('text.plot-tick[text-anchor="end"]')].map((node) => node.textContent);

    const pointer = (type, pointerType, clientX, clientY, target = svg) => {
        act(() => {
            target.dispatchEvent(new PointerEvent(type, {
                bubbles: true, cancelable: true, pointerId: 1, pointerType, clientX, clientY,
            }));
        });
    };
    // A press, a move to each point in turn, and a release at the last.
    const gesture = (pointerType, points) => {
        const [first, ...rest] = points;
        pointer('pointerdown', pointerType, ...first);
        rest.forEach((point) => pointer('pointermove', pointerType, ...point));
        const last = points[points.length - 1];
        pointer('pointerup', pointerType, ...last);
    };
    // React derives onPointerLeave from pointerout.
    const leave = (pointerType) => pointer('pointerout', pointerType, at(90), 200);
    const readout = () => container.querySelector('.plot-tooltip');
    const selection = () => container.querySelector('.zoom-selection');
    const resetPill = () => [...container.querySelectorAll('button')].find((b) => b.textContent === 'Reset zoom');

    it('leaves vertical drags and pinches to the page', () => {
        const rule = plotCss.match(/\.interactive-plot svg\s*\{([^}]*)\}/)[1];
        expect(rule).toMatch(/touch-action:\s*pan-y pinch-zoom;/);
    });

    it('zooms x only on a sideways swipe, however the finger drifts', () => {
        const before = yTicks();
        // 60° to 120°, drifting 40 user units down: a mouse would zoom y too.
        gesture('touch', [[at(60), 120], [at(90), 140], [at(120), 160]]);
        const [lo, hi] = xRange();
        expect(lo).toBeCloseTo(60, 6);
        expect(hi).toBeCloseTo(120, 6);
        // The peak at 90° is in the window, so the re-fitted y range is the
        // same as before.
        expect(yTicks()).toEqual(before);
        expect(selection()).toBeNull();
        expect(resetPill()).toBeDefined();
    });

    it('previews a finger selection as a band the plot height', () => {
        pointer('pointerdown', 'touch', at(60), 120);
        pointer('pointermove', 'touch', at(120), 160);
        const band = selection();
        expect(Number(band.getAttribute('x'))).toBeCloseTo(at(60), 6);
        expect(Number(band.getAttribute('width'))).toBeCloseTo(at(120) - at(60), 6);
        expect(Number(band.getAttribute('y'))).toBe(16);
        expect(Number(band.getAttribute('height'))).toBe(376);
    });

    it('keeps the mouse box zoom', () => {
        const before = yTicks();
        gesture('mouse', [[at(60), 120], [at(90), 140], [at(120), 160]]);
        const [lo, hi] = xRange();
        expect(lo).toBeCloseTo(60, 6);
        expect(hi).toBeCloseTo(120, 6);
        expect(yTicks()).not.toEqual(before);
    });

    it('zooms nothing when the browser takes the drag as a page scroll', () => {
        pointer('pointerdown', 'touch', at(90), 100);
        pointer('pointermove', 'touch', at(91), 140);
        expect(selection()).not.toBeNull();
        pointer('pointercancel', 'touch', at(91), 140);
        expect(selection()).toBeNull();
        expect(resetPill()).toBeUndefined();
        expect(readout()).toBeNull();
    });

    it('shows the readout on a tap and keeps it after the finger lifts', () => {
        // A tap may wobble, up to 16 CSS px across.
        gesture('touch', [[at(90), 200], [at(90) + 10, 205]]);
        expect(resetPill()).toBeUndefined();
        expect(readout()).not.toBeNull();
        expect(readout().querySelector('strong').textContent).toMatch(/: 93\.00$/);
        leave('touch');
        expect(readout()).not.toBeNull();
        // The next press takes it, wherever it lands.
        pointer('pointerdown', 'touch', 20, 200);
        expect(readout()).toBeNull();
    });

    it('leaves a tap readout out of a saved figure', () => {
        gesture('touch', [[at(90), 200]]);
        expect(svg.querySelector('.hover-line')).not.toBeNull();
        const figure = svgToSvgString(svg);
        expect(figure).not.toMatch(/hover-(line|dot)/);
        expect(figure).toMatch(/series-path/);
    });

    it('still hovers and clears for the mouse, and never on a finger move', () => {
        pointer('pointermove', 'mouse', at(45), 200);
        expect(readout()).not.toBeNull();
        leave('mouse');
        expect(readout()).toBeNull();
        pointer('pointermove', 'touch', at(45), 200);
        expect(readout()).toBeNull();
    });

    it('resets on a double tap and drops the readout', () => {
        gesture('touch', [[at(60), 120], [at(120), 120]]);
        expect(resetPill()).toBeDefined();
        gesture('touch', [[at(90), 200]]);
        expect(readout()).not.toBeNull();
        act(() => {
            svg.dispatchEvent(new MouseEvent('dblclick', { bubbles: true }));
        });
        expect(resetPill()).toBeUndefined();
        expect(readout()).toBeNull();
        const [lo, hi] = xRange();
        expect(lo).toBeCloseTo(0, 6);
        expect(hi).toBeCloseTo(180, 6);
    });
});
