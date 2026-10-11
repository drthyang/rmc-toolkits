// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The AI Copilot's marks on a plot (highlight_cards → Dashboard → marks):
// a dashed guide, a ▽ and a label per mark inside the x range, nothing for a
// mark outside it, and no markup at all without marks.

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import InteractivePlot from '../components/InteractivePlot';
import { markWindow } from '../plotDomain';

const r = Array.from({ length: 200 }, (_, i) => 0.05 * (i + 1));
const partials = {
    kind: 'pdf_partials',
    title: 'Partial PDFs',
    xLabel: 'r (Å)',
    yLabel: 'g(r)',
    series: [{ label: 'Ta-Se', x: r, y: r.map((v) => Math.exp(-((v - 2.6) ** 2) / 0.01)) }]
};
const render = (marks) => renderToStaticMarkup(
    <InteractivePlot file={{ path: 'p', name: 'p' }} plotData={partials} marks={marks} />
);

describe('InteractivePlot marks', () => {
    it('draws a guide, a ▽ and a label for each mark in range', () => {
        const html = render([{ x: 2.6, label: 'Ta–Se 2.60 Å' }, { x: 8.2 }, { x: 40, label: 'outside' }]);
        expect(html).toContain('class="plot-copilot-marks"');
        expect(html.match(/class="plot-copilot-mark"/g)).toHaveLength(2);
        expect(html).toContain('Ta–Se 2.60 Å</text>');
        expect(html).not.toContain('outside');
        // The mark in the right part of the plot labels to its left (none here).
        expect(html.match(/<text[^>]*text-anchor="start"[^>]*>Ta–Se/)).not.toBeNull();
    });

    it('zooms x to the marks, in context', () => {
        const near = markWindow([{ x: 2.6 }, { x: 2.97 }], [0, 50]);
        expect(near[0]).toBe(0);
        expect(near[1]).toBeCloseTo(5.97, 12);
        expect(markWindow([{ x: 20 }, { x: 30 }], [0, 50])).toEqual([12.5, 37.5]);
        // Nothing in range, or a window that would be most of the plot: no zoom.
        expect(markWindow([{ x: 80 }], [0, 50])).toBeNull();
        expect(markWindow([{ x: 2 }, { x: 45 }], [0, 50])).toBeNull();
        expect(markWindow([], [0, 50])).toBeNull();
        // The rendered axis follows: ticks up to 6 Å, not 10.
        const ticks = [...render([{ x: 2.6, label: 'Ta–Se' }]).matchAll(/class="plot-tick" x="[^"]+" y="[^"]+" text-anchor="middle">([^<]+)</g)].map((m) => m[1]);
        expect(ticks.at(-1)).not.toBe('10');
        expect(Number(ticks.at(-1))).toBeLessThanOrEqual(6);
    });

    it('renders exactly as before without marks', () => {
        expect(render(null)).toBe(render([]));
        expect(render(null)).not.toContain('plot-copilot');
    });
});
