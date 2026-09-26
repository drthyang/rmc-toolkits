// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Flask ⟷ browser parity on RMCProfile output layouts no real run on hand
// covers: neutron *_PDFn / *_SQn, *_bragg (ToF and Q), EXAFS Q/R, partials —
// with NaN-masked regions, CRLF, trailing commas, a leading blank line and
// E-notation. The golden is the FLASK payload for each file
// (tests/generate_plot_parity_fixture.py; tests/test_parsers_plot_payload.py
// fails if it is stale); the browser parser must reproduce it.

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { detectPlotKind, plotDataFromText, plotMetadataFromFile } from '../browserData';

const { cases } = JSON.parse(readFileSync(new URL('./fixtures/plot_parity_fixture.json', import.meta.url), 'utf8'));

const summarize = (plot, metadata) => ({
    kind: plot.kind,
    title: plot.title,
    metadataTitle: metadata.title,
    xLabel: plot.xLabel,
    yLabel: plot.yLabel,
    metrics: plot.metrics,
    series: plot.series.map((series) => {
        const finite = series.y.filter(Number.isFinite);
        return {
            label: series.label,
            n: series.y.length,
            finite: finite.length,
            first: finite.length ? finite[0] : null,
            last: finite.length ? finite[finite.length - 1] : null,
        };
    }),
});

const expectClose = (actual, expected) => {
    if (expected === null) expect(actual).toBeNull();
    else expect(actual).toBeCloseTo(expected, 12);
};

describe('browser plot payloads match Flask', () => {
    it.each(cases.map((entry) => [entry.name, entry]))('%s', (name, entry) => {
        const file = { name, path: `run/${name}`, plotKind: detectPlotKind(name), text: entry.text };
        const plot = plotDataFromText(file);
        const actual = summarize(plot, plotMetadataFromFile({ ...file, plotData: plot }));
        const { expected } = entry;
        expect([actual.kind, actual.title, actual.metadataTitle, actual.xLabel, actual.yLabel])
            .toEqual([expected.kind, expected.title, expected.metadataTitle, expected.xLabel, expected.yLabel]);
        expect(Object.keys(actual.metrics).sort()).toEqual(Object.keys(expected.metrics).sort());
        Object.entries(expected.metrics).forEach(([key, value]) => expectClose(actual.metrics[key], value));
        expect(actual.series.map(({ label, n, finite }) => [label, n, finite]))
            .toEqual(expected.series.map(({ label, n, finite }) => [label, n, finite]));
        actual.series.forEach((series, index) => {
            expectClose(series.first, expected.series[index].first);
            expectClose(series.last, expected.series[index].last);
        });
    });
});

// Same cases as tests/test_parsers.py::test_read_rmc_csv_cell_rules_match_the_browser.
describe('CSV cell rules shared with read_rmc_csv', () => {
    const plotOf = (text) => plotDataFromText({ name: 'run_FQ1.csv', plotKind: 'xray_sq', text });

    it('skips blank lines and keeps NaN / **** cells as masked values', () => {
        const plot = plotOf('\nQ, F(Q)_RMC, F(Q)_Expt\n\n1.0, NaN, 0.5\n2.0, ****, 1.5D-01\n');
        expect(plot.series.map((series) => series.label)).toEqual(['F(Q)_RMC', 'F(Q)_Expt']);
        expect(plot.series[0].y.every(Number.isNaN)).toBe(true);
        expect(plot.series[1].y).toEqual([0.5, 0.15]);
    });

    it('rejects a non-numeric cell, naming the true file line', () => {
        expect(() => plotOf('Q, a, b\n\n1.0, 2.0, 3.0\n2.0, abc, 3.0\n')).toThrow("line 4: 'abc' is not a number");
    });
});
