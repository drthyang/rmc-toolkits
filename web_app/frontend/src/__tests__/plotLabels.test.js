// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Axis labels and titles must name the function the file holds, identically to
// plots.series_titles (tests/test_plots.py::FunctionLabelTests and
// tests/test_parsers_plot_payload.py pin the same files).

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { buildRunContext } from '../llm/context/runContext';
import { detectPlotKind, plotDataFromText, plotMetadataFromFile, seriesTitles } from '../browserData';

const demoFile = (name) => ({
    name,
    path: `Demo/${name}`,
    plotKind: detectPlotKind(name),
    text: readFileSync(new URL(`../../public/demo/${name}`, import.meta.url), 'utf8'),
});

describe('function labels', () => {
    it('labels the demo *_FQ1.csv F(Q) (header F(Q)_RMC, F(Q)_Expt), not S(Q)', () => {
        const file = demoFile('GTS_250K_FQ1.csv');
        expect(plotMetadataFromFile(file).title).toBe('F(Q)');   // before parsing
        const plot = plotDataFromText(file);
        expect([plot.kind, plot.title, plot.yLabel]).toEqual(['xray_sq', 'F(Q)', 'F(Q)']);
        const context = buildRunContext({ plotFiles: [{ ...file, plotData: plot }] });
        expect(context.datasets[0].title).toBe('F(Q)');
    });

    it('labels the demo PDFpartials partial g(r), not G(r)', () => {
        const plot = plotDataFromText(demoFile('GTS_250K_PDFpartials.csv'));
        expect([plot.title, plot.yLabel, plot.xLabel]).toEqual(['Partial g(r)', 'g(r)', 'r (Å)']);
    });

    it('charts reciprocal-space datasets numbered 2 and up', () => {
        expect(detectPlotKind('run_FQ2.csv')).toBe('xray_sq');
        expect(detectPlotKind('run_SQ12.csv')).toBe('neutron_sq');
        expect(detectPlotKind('run_FQ1partials.csv')).toBeNull();
        expect(detectPlotKind('run_XFQ1.csv')).toBeNull();
        const plot = plotDataFromText({
            name: 'run_SQ2.csv', plotKind: 'neutron_sq', text: 'Q, S(Q)_RMC, S(Q)_Expt\n1.0, 0.7, 1.0\n2.0, 1.4, 2.0\n'
        });
        expect([plot.title, plot.yLabel]).toEqual(['S(Q) #2', 'S(Q)']);
        expect(plot.metrics.rwp).toBeCloseTo(0.3, 12);
    });

    it('matches plots.series_titles case by case', () => {
        expect(seriesTitles('xray_sq', 'run_FQ1.csv', ['Q', 'F(Q)_RMC', 'F(Q)_Expt'])).toEqual(['F(Q)', 'F(Q)']);
        expect(seriesTitles('xray_sq', 'run_FQ2.csv', ['Q', 'a', 'b'])).toEqual(['F(Q) #2', 'F(Q)']);
        expect(seriesTitles('neutron_sq', 'run_SQ1.csv', ['Q', 'a', 'b'])).toEqual(['S(Q)', 'S(Q)']);
        expect(seriesTitles('neutron_sq', 'run_SQ1.csv', ['Q', 'F(Q)_RMC', 'F(Q)_Expt'])).toEqual(['F(Q)', 'F(Q)']);
        expect(seriesTitles('npdf', 'run_PDF1.csv', ['r', 'D(r)_RMC', 'D(r)_Expt'])).toEqual(['PDF1', 'D(r)']);
        expect(seriesTitles('npdf', 'run_PDF2.csv', ['r', 'calc', 'expt'])).toEqual(['PDF2', 'G(r)']);
    });

    it('labels a STOG .fq file F(Q) by default', () => {
        const plot = plotDataFromText({ name: 'scale_ft_rmc.fq', plotKind: 'stog', text: '2\ntitle\n0.5 -0.9\n1.0 -0.5\n' });
        expect(plot.yLabel).toBe('F(Q)');
        expect(plotMetadataFromFile({ name: 'scale_ft_rmc.fq', plotKind: 'stog' }).title).toBe('F(Q)');
    });
});
