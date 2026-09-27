// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The angle-axis options of InteractivePlot (xDomain, xTicks, xMinorStep,
// xGrid, yMin, and per-series curve / fill / width / legend) are opt-in: a
// payload that sets none of them must render exactly the markup it rendered
// before they existed. The snapshots below were written by the component as
// it was on main 28ca96f, before the options were added, for payloads shaped
// like the Dashboard's and Auto StoG's; they pin every other plot in the app.

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import InteractivePlot from '../components/InteractivePlot';

const range = (n, f) => Array.from({ length: n }, (_, i) => f(i));
const q = range(60, (i) => 0.5 + i * 0.4);
const r = range(80, (i) => 0.1 * (i + 1));

const render = (plotData, variant) => renderToStaticMarkup(
    <InteractivePlot file={{ path: 'p', name: 'p' }} plotData={plotData} variant={variant} />
);

// Dashboard: a measured/calculated pair (hollow markers + curve).
const dashboardFq = {
    kind: 'xray_sq',
    title: 'F(Q)',
    xLabel: 'Q (Å^{-1})',
    yLabel: 'F(Q)',
    series: [
        { label: 'F(Q)_RMC', x: q, y: q.map((v) => Math.sin(v) / v) },
        { label: 'F(Q)_Expt', x: q, y: q.map((v) => Math.sin(v) / v + 0.01 * Math.cos(3 * v)) },
    ],
};

// Dashboard: an R-value history with counts large enough to widen the y ticks.
const dashboardChi = {
    kind: 'r_value',
    title: 'chi²',
    xLabel: 'generated moves',
    yLabel: 'ln(chi²)',
    series: [{ label: 'chi²', x: range(40, (i) => i * 25000), y: range(40, (i) => 9000 - 120 * i) }],
};

// Auto StoG S(Q): hidden series, a neutral guide and coloured guides.
const autoStogSq = {
    title: 'S(Q) — scaled and filtered',
    xLabel: 'Q (Å^{-1})',
    yLabel: 'S(Q)',
    series: [
        { label: 'Auto-scaled a·S + b', x: q, y: q.map((v) => 1 + Math.sin(2 * v) / v) },
        { label: 'Filtered S(Q)', x: q, y: q.map((v) => 1 + 0.9 * Math.sin(2 * v) / v) },
        { label: 'Measured (unscaled)', x: q, y: q.map((v) => 3 + Math.sin(2 * v)), defaultHidden: true },
        { label: 'S → 1', x: [q[0], q[q.length - 1]], y: [1, 1], role: 'guide' },
        { label: 'Level 0.98', x: [10, 20], y: [0.98, 0.98], role: 'guide', color: '#4c7df0' },
    ],
};

// Auto StoG G_K(r): an initial y window.
const autoStogGk = {
    title: 'G_K(r) — full range',
    xLabel: 'r (Å)',
    yLabel: 'G_K(r)',
    series: [
        { label: 'G_K(r) output (RMC file)', x: r, y: r.map((v) => (v < 2 ? -0.5 : Math.sin(3 * v) / v)) },
        { label: '−⟨b⟩² theory', x: [0, 8], y: [-0.5, -0.5], role: 'guide' },
    ],
    initialYDomain: [-1.05, 1.6],
};

// The Bond Geometry partial g(r) helper as it was: one curve, two guides.
const partialGr = {
    title: 'Ta-Se partial g(r)',
    xLabel: 'r (Å)',
    yLabel: 'g(r)',
    series: [
        { label: 'Ta-Se', x: r.slice(0, 60), y: r.slice(0, 60).map((v) => 8 * Math.exp(-((v - 2.6) ** 2) / 0.01)) },
        { label: 'rmin 2.0', x: [2, 2], y: [0, 8], role: 'guide' },
        { label: 'rmax 3.0', x: [3, 3], y: [0, 8], role: 'guide' },
    ],
};

describe('InteractivePlot markup without the angle-axis options', () => {
    it('Dashboard F(Q) pair (grid card and wide variants)', () => {
        expect(render(dashboardFq)).toMatchSnapshot();
        expect(render(dashboardFq, 'wide')).toMatchSnapshot();
    });

    it('Dashboard chi² history (wide y ticks)', () => {
        expect(render(dashboardChi)).toMatchSnapshot();
    });

    it('Auto StoG S(Q) with hidden series and guides', () => {
        expect(render(autoStogSq, 'wide')).toMatchSnapshot();
    });

    it('Auto StoG G_K(r) with an initial y window', () => {
        expect(render(autoStogGk)).toMatchSnapshot();
    });

    it("the 'fit' variant before its first measurement", () => {
        expect(render(partialGr, 'fit')).toMatchSnapshot();
    });
});

// A 1° angle histogram in the Bond Geometry page's shape.
const centers = range(180, (i) => i + 0.5);
const anglePayload = (extra = {}) => ({
    title: 'Se–Ta–Se bond angles',
    xLabel: 'angle at Ta, θ (°)',
    yLabel: 'sin-corrected (random = 1)',
    xDomain: [0, 180],
    xTicks: [0, 30, 60, 90, 120, 150, 180],
    xMinorStep: 10,
    xGrid: true,
    yMin: 0,
    series: [
        { label: 'sin-corrected', x: centers, y: centers.map((c) => 1 + 5 * Math.exp(-((c - 90) ** 2) / 8)), curve: 'step', binWidth: 1, fill: true, width: 1.75 },
        { label: 'random bonds', x: [0, 180], y: [1, 1], role: 'guide' },
    ],
    ...extra,
});

const tickLabels = (markup) => [...markup.matchAll(/<text class="plot-tick" x="([\d.]+)" y="([\d.]+)" text-anchor="(middle|end)">([^<]*)<\/text>/g)]
    .map(([, x, y, anchor, label]) => ({ x: Number(x), y: Number(y), anchor, label }));

describe('InteractivePlot angle-axis options', () => {
    // Grid-card box: 720 × 450, left 60, right 18, top 16, bottom 58.
    const markup = render(anglePayload());
    const xLabels = tickLabels(markup).filter((tick) => tick.anchor === 'middle');
    const yLabels = tickLabels(markup).filter((tick) => tick.anchor === 'end');

    it('fixes the x domain to 0–180 with no padding and labels every 30°', () => {
        expect(xLabels.map((tick) => tick.label)).toEqual(['0', '30', '60', '90', '120', '150', '180']);
        expect(xLabels[0].x).toBe(60);
        expect(xLabels[6].x).toBe(702);
        expect(markup).toContain('<rect class="plot-bg" x="60" y="16" width="642" height="376">');
    });

    it('adds unlabelled minor marks every 10° and vertical grid lines at the majors', () => {
        expect(markup.match(/plot-tick-mark--minor/g)).toHaveLength(12);
        const vertical = [...markup.matchAll(/<line class="plot-grid-line" x1="([\d.]+)" x2="([\d.]+)" y1="16" y2="392">/g)];
        expect(vertical.map(([, x1]) => Number(x1))).toEqual([60, 167, 274, 381, 488, 595, 702]);
    });

    it('starts the y axis at 0, on the frame', () => {
        expect(yLabels[0].label).toBe('0');
        expect(yLabels[0].y).toBe(392 + 4.5);
    });

    it('draws the step outline across each bin, with an area to y = 0', () => {
        const step = markup.match(/<path class="series-path series-path--step" d="([^"]+)" stroke="#1f6fd6" style="stroke-width:1.75">/);
        expect(step).not.toBeNull();
        // First bin: flat from 0° to 1° (x 60 → 63.57).
        expect(step[1].startsWith('M 60.00 ')).toBe(true);
        expect(step[1]).toMatch(/^M 60\.00 ([\d.]+) L 63\.57 \1 L 63\.57 /);
        const area = markup.match(/<path class="series-area" d="([^"]+)" fill="#1f6fd6">/);
        expect(area).not.toBeNull();
        expect(area[1].startsWith('M 60.00 392.00 L 60.00 ')).toBe(true);
        expect(area[1].endsWith('L 702.00 392.00 Z')).toBe(true);
    });

    it('keeps the guide dashed and out of the palette; legend:false leaves the legend', () => {
        expect(markup).toContain('class="series-path series-path--guide" d="M 60.00 ');
        expect(markup).toContain('random bonds</button>');
        const hidden = render(anglePayload({
            series: [...anglePayload().series.slice(0, 1), { label: 'rmin 2.0', x: [2, 2], y: [0, 1], role: 'guide', legend: false }],
        }));
        expect(hidden).not.toContain('rmin 2.0</button>');
        expect(hidden).toContain('series-path--guide');
    });

    it('a step series breaks at a non-finite value', () => {
        const gap = render(anglePayload({
            series: [{ label: 's', x: [0.5, 1.5, 2.5, 3.5], y: [1, null, 2, 2], curve: 'step', binWidth: 1 }],
        }));
        const d = gap.match(/<path class="series-path series-path--step" d="([^"]+)"/)[1];
        expect(d.match(/M /g)).toHaveLength(2);
    });
});
