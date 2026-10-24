// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// RMCProfile .log reads must survive Live Data polls that land mid-write, and
// must agree with Python on non-finite rows. Same cases as
// tests/test_parsers.py::ReadChiLogTests.

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { chiTermsPlotData, combineChiTerms, plotDataFromText, readChi } from '../browserData';
import { buildRunContext } from '../llm/context/runContext';
import { classifyConvergence } from '../llm/watchdog/heuristics';

const DEMO_LOG = new URL('../../public/demo/GTS_250K-02.log', import.meta.url);
const HEADER = 'Time  moves_acc moves_gen  F(Q)_1  X_ray_(R)1\nh/m/s/.th WEIGHT PARAMETERS 0.1E+01 0.1E+01\n';

// Last-column values of the complete data rows, read independently of the parser.
const lastColumn = (text) => text.split('\n').slice(2)
    .filter((line) => line.trim())
    .map((line) => Number(line.trim().split(/\s+/).pop()));

describe('readChi', () => {
    const text = readFileSync(DEMO_LOG, 'utf8');
    const complete = lastColumn(text);

    it('reads the complete demo log and names its column', () => {
        const log = readChi(text);
        expect(log.values).toEqual(complete);
        expect(log.column).toBe('X_ray_(R)1');
        expect(log.skippedRows).toBe(0);
    });

    it('never turns a half-written final line into the final chi^2', () => {
        // Cut at every byte inside the final line: a truncated move counter, a
        // lone "0." or a truncated mantissa used to become the final value.
        const finalStart = text.replace(/\n$/, '').lastIndexOf('\n') + 1;
        for (let offset = finalStart; offset < text.length; offset += 1) {
            const plot = plotDataFromText({ plotKind: 'r_value', name: 'GTS_250K-02.log', text: text.slice(0, offset) });
            expect(plot.series[0].y).toHaveLength(complete.length - 1);
            expect(plot.metrics.final_chi_r).toBe(complete[complete.length - 2]);
        }
    });

    it('skips rows whose token count differs from the header', () => {
        const log = readChi(`${HEADER}1.0 10 20 0.3E-02 0.2E-03\n2.0 20 40 0.117\n3.0 30 60 0.3E-02 0.1D-03\n`);
        expect(log.values).toEqual([0.2e-3, 0.1e-3]);
        expect(log.skippedRows).toBe(1);
    });

    it('keeps non-finite chi^2 rows as NaN (agreeing with Python)', () => {
        const log = readChi(`${HEADER}1.0 10 20 0.3E-02 0.2E-03\n2.0 20 40 0.3E-02 NaN\n3.0 30 60 0.3E-02 **********\n`);
        expect(log.values).toHaveLength(3);
        expect(log.values[0]).toBe(0.2e-3);
        expect(log.values.slice(1).every(Number.isNaN)).toBe(true);
    });

    it('falls back to the first row count for a log with no column header', () => {
        const log = readChi('header\nheader\n1 0.1 10.0\n2 0.2\n3 0.3 30.0\n');
        expect(log.values).toEqual([10, 30]);
        expect(log.column).toBeNull();
        expect(log.skippedRows).toBe(1);
        // No column names: which columns are χ² terms is unknown.
        expect(log.terms).toBeNull();
    });

    it('reads every χ² term (the columns after Time, moves_acc, moves_gen), as test_parsers.py', () => {
        const header = 'Time  moves_acc moves_gen  F(Q)_1  Curvature  X_ray_(R)1\nh/m/s/.th WEIGHT PARAMETERS 0.1E+01 0.1E+01 0.1E+01\n';
        const log = readChi(`${header}1.0 10 20 0.3E-02 0.000E+00 0.2E-03\n2.0 20 40 0.117\n3.0 30 60 0.2D-02 NaN 0.1D-03\n`);
        expect(log.terms.names).toEqual(['F(Q)_1', 'Curvature', 'X_ray_(R)1']);
        expect(log.terms.values[0]).toEqual([0.3e-2, 0.2e-2]);
        expect(log.terms.values[1][0]).toBe(0);
        expect(Number.isNaN(log.terms.values[1][1])).toBe(true);
        expect(log.terms.values[2]).toEqual(log.values);
        const plot = chiTermsPlotData(log.terms);
        expect(plot.total[0]).toBe(Math.log(0.3e-2 + 0 + 0.2e-3));
        expect(Number.isNaN(plot.total[1])).toBe(true);
        expect(plot.series[1][0]).toBe(Math.log(1e-12));
        expect(plot.segments).toEqual([2]);
    });

    it('gives the demo log\'s eight terms, the last one the charted column', () => {
        const plot = plotDataFromText({ plotKind: 'r_value', name: 'GTS_250K-02.log', text });
        expect(plot.chiTerms.names).toEqual(['F(Q)_1', ...Array(6).fill('Curvature'), 'X_ray_(R)1']);
        expect(plot.chiTerms.series.at(-1)).toEqual(plot.series[0].y);
        expect(plot.chiTerms.total).toHaveLength(complete.length);
        expect(plot.chiTerms.segments).toEqual([complete.length]);
    });
});

describe('combineChiTerms: one run\'s logs end to end', () => {
    const part = (names, series) => ({ names, series, total: series[0].map((_, i) => series.reduce((sum, s) => sum + s[i], 0)), segments: [series[0].length] });

    it('concatenates the terms when every log names the same, the total and segments always', () => {
        const a = part(['F(Q)_1', 'X_ray_(R)1'], [[1, 2], [3, 4]]);
        const b = part(['F(Q)_1', 'X_ray_(R)1'], [[5], [6]]);
        expect(combineChiTerms([a, b])).toEqual({ names: ['F(Q)_1', 'X_ray_(R)1'], series: [[1, 2, 5], [3, 4, 6]], total: [4, 6, 11], segments: [2, 1] });
        const changed = part(['F(Q)_1', 'X_ray_(R)1_new'], [[7], [8]]);
        expect(combineChiTerms([a, changed])).toEqual({ names: null, series: null, total: [4, 6, 15], segments: [2, 1] });
        // A log without terms (no column names): none for the run.
        expect(combineChiTerms([a, null])).toBeNull();
    });
});

describe('a blown-up (NaN) log tail', () => {
    // 60 improving rows, then 40 rows of NaN: the browser used to drop the NaN
    // rows, so the watchdog saw only the improving part.
    const rows = Array.from({ length: 100 }, (_, index) => (
        index < 60 ? `${index} 1 2 0.3E-02 ${(0.02 * Math.exp(-index / 20)).toExponential(3)}` : `${index} 1 2 0.3E-02 NaN`
    ));
    const plot = plotDataFromText({ plotKind: 'r_value', name: 'run-00.log', text: `${HEADER}${rows.join('\n')}\n` });

    it('keeps every row and ends on a non-finite value', () => {
        expect(plot.series[0].y).toHaveLength(100);
        expect(Number.isNaN(plot.metrics.final_chi_r)).toBe(true);
    });

    it('is classified as diverging, never improving', () => {
        expect(classifyConvergence(plot.series[0].y)).toBe('diverging');
        // The finite prefix alone is improving — what the old parser reported.
        expect(classifyConvergence(plot.series[0].y.slice(0, 60))).toBe('improving');
    });

    it('reaches the AI context as a counted non-finite tail', () => {
        const context = JSON.parse(JSON.stringify(buildRunContext({ rValueFile: { plotData: plot } })));
        expect(context.convergence.non_finite_steps).toBe(40);
        expect(context.convergence.last).toBeNull();
        expect(context.convergence.final_chi_squared).toBeUndefined();
        expect(Number.isFinite(context.convergence.min)).toBe(true);
    });
});
