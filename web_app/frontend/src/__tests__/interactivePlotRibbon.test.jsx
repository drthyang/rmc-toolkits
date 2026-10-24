// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// InteractivePlot's opt-in `ribbon` ({ series, sigma, label }): a ±σ band
// around one series, under it in its colour, kept inside the y range, split at
// a gap (a non-finite y or σ), listed in the legend. Absent, nothing changes
// (the other interactivePlot* tests pin the plot without it).

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import InteractivePlot from '../components/InteractivePlot';

const file = { path: 'run_FT_XFQ1.csv', name: 'run_FT_XFQ1.csv' };
const x = [1, 2, 3, 4, 5, 6];
const plotData = {
    kind: 'xpdf', title: 'xPDF', xLabel: 'r (Å)', yLabel: 'G(r)',
    series: [
        { label: 'X_ray-calc', x, y: [0, 1, 0, -1, 0, 1] },
        { label: 'X_ray_exp_renorm', x, y: [0, 1, 0, -1, 0, 1] }
    ]
};

const ticks = (html) => [...html.matchAll(/class="plot-tick"[^>]*text-anchor="end"[^>]*>([^<]+)</g)].map((match) => Number(match[1]));

describe('InteractivePlot ribbon', () => {
    it('draws one closed band per run of finite points, under the series, in the legend', () => {
        const sigma = [0.1, 0.1, NaN, 0.1, 0.1, 0.1];
        const html = renderToStaticMarkup(
            <InteractivePlot file={file} plotData={plotData} ribbon={{ series: 'X_ray_exp_renorm', sigma, label: '±1σ propagated' }} />
        );
        const band = /<path class="series-ribbon" d="([^"]+)"/.exec(html);
        expect(band).not.toBeNull();
        // Two runs (the NaN σ splits them): [1, 2] and [4, 5, 6].
        expect(band[1].match(/M /g)).toHaveLength(2);
        expect(band[1].match(/Z/g)).toHaveLength(2);
        expect(html.indexOf('series-ribbon')).toBeLessThan(html.indexOf('series-markers'));
        expect(html).toContain('swatch-band');
        expect(html).toContain('±1σ propagated');
    });

    it('keeps the band inside the y range', () => {
        const wide = renderToStaticMarkup(
            <InteractivePlot file={file} plotData={plotData} ribbon={{ series: 'X_ray_exp_renorm', sigma: x.map(() => 3), label: '±1σ' }} />
        );
        const plain = renderToStaticMarkup(<InteractivePlot file={file} plotData={plotData} />);
        expect(Math.max(...ticks(wide))).toBeGreaterThanOrEqual(4);
        expect(Math.min(...ticks(wide))).toBeLessThanOrEqual(-4);
        expect(Math.max(...ticks(plain))).toBeLessThan(4);
    });

    it('draws nothing for an unknown series or a σ of another length', () => {
        for (const ribbon of [
            { series: 'nope', sigma: x.map(() => 0.1), label: '±1σ' },
            { series: 'X_ray_exp_renorm', sigma: [0.1, 0.1], label: '±1σ' }
        ]) {
            const html = renderToStaticMarkup(<InteractivePlot file={file} plotData={plotData} ribbon={ribbon} />);
            expect(html).not.toContain('series-ribbon');
        }
    });
});
