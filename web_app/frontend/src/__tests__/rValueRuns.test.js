// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The chi^2 history must be ONE run's logs, like Flask's related_r_value_logs():
// static mode used to flat-map every visible log, splicing different runs into
// one "convergence" curve (and taking final_chi_r from whichever sorted last).
// It must also be labelled by its log column, not as a total "R-value".

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { chooseRValueGroup, combineRValueFiles, plotDataFromText, rValueGroupKey } from '../browserData';
import { buildRunContext } from '../llm/context/runContext';

const demo = (name) => readFileSync(new URL(`../../public/demo/${name}`, import.meta.url), 'utf8');
const HEADER = 'Time       moves_acc moves_gen     F(Q)_1   X_ray_(R)1\nh/m/s/.th    WEIGHT PARAMETERS   0.100E+01  0.100E+01\n';

const parsedLog = (path, text) => {
    const name = path.split('/').pop();
    return { name, path, plotKind: 'r_value', sourceFile: {}, plotData: plotDataFromText({ plotKind: 'r_value', name, text }) };
};

// Display order (Dashboard comparePlotFiles): lower-case stem, then sequence.
const FILES = [
    parsedLog('run/GTS_250K-00.log', demo('GTS_250K-00.log')),
    parsedLog('run/GTS_250K-01.log', demo('GTS_250K-01.log')),
    parsedLog('run/GTS_250K-02.log', demo('GTS_250K-02.log')),
    parsedLog('run/other-00.log', `${HEADER}1.0 1 2 0.5E+00 0.110E+01\n2.0 2 4 0.5E+00 0.100E+01\n`),
];
const rowsOf = (file) => file.plotData.series[0].y.length;
const GTS_ROWS = FILES.slice(0, 3).reduce((sum, file) => sum + rowsOf(file), 0);

describe('chi² logs of one run', () => {
    it('groups logs by folder and exact stem', () => {
        expect(rValueGroupKey('run/GTS_250K-01.log')).toBe('run/GTS_250K');
        expect(rValueGroupKey('run/sub/GTS_250K-01.log')).toBe('run/sub/GTS_250K');
        const { group, others } = chooseRValueGroup(FILES);
        expect(group.map((file) => file.name)).toEqual(['GTS_250K-00.log', 'GTS_250K-01.log', 'GTS_250K-02.log']);
        expect(others).toHaveLength(1);
    });

    it('concatenates only the run the structure file belongs to', () => {
        const combined = combineRValueFiles(FILES, 'run/GTS_250K.rmc6f');
        expect(combined.plotData.series[0].y).toHaveLength(GTS_ROWS);
        expect(combined.plotData.metrics.final_chi_r).toBe(FILES[2].plotData.metrics.final_chi_r);
        expect(combined.sourceNames).toEqual(['GTS_250K-00.log', 'GTS_250K-01.log', 'GTS_250K-02.log']);
        expect(combined.otherRuns).toEqual(['other']);

        const other = combineRValueFiles(FILES, 'run/other.rmc6f');
        expect(other.plotData.series[0].y).toHaveLength(2);
        expect(other.plotData.metrics.final_chi_r).toBe(1.0);
    });

    it('skips a header-only restart log silently, as read_chi_log does', () => {
        // The state right after RMCProfile starts a restart: header + WEIGHT line.
        const name = 'GTS_250K-03.log';
        let parseError = null;
        try {
            plotDataFromText({ plotKind: 'r_value', name, text: HEADER });
        } catch (error) {
            parseError = error.message;
        }
        expect(parseError).toContain('does not contain chi values');
        const restarting = { name, path: `run/${name}`, plotKind: 'r_value', sourceFile: {}, parseError };
        const combined = combineRValueFiles([...FILES.slice(0, 3), restarting], 'run/GTS_250K.rmc6f');
        expect(combined.parseError).toBe('');
        expect(combined.plotData.series[0].y).toHaveLength(GTS_ROWS);
    });

    it('names a restart whose fit term changed by both columns (read_chi_log parity)', () => {
        const changed = parsedLog('run/GTS_250K-03.log',
            `${HEADER.replace('X_ray_(R)1', 'X_ray_(R)1_new')}9.0 1 2 0.5E+00 0.110E+01\n`);
        const combined = combineRValueFiles([...FILES.slice(0, 3), changed], 'run/GTS_250K.rmc6f');
        expect(combined.plotData.title).toBe('χ² history: X_ray_(R)1 / X_ray_(R)1_new');
        expect(combined.plotData.chiColumn).toBeNull();
    });

    it('falls back to the first run (as Flask does) without a structure file', () => {
        expect(combineRValueFiles(FILES).plotData.series[0].y).toHaveLength(GTS_ROWS);
    });

    it('labels the series by the log column and tells the AI context which one', () => {
        const combined = combineRValueFiles(FILES, 'run/GTS_250K.rmc6f');
        expect(combined.plotData.title).toBe('χ² history: X_ray_(R)1');
        expect(combined.plotData.series[0].label).toBe('X_ray_(R)1');
        expect(combined.plotData.yLabel).toBe('ln(χ²)');
        const context = buildRunContext({ rValueFile: combined });
        expect(context.convergence.column).toBe('X_ray_(R)1');
        expect(context.convergence.quantity).toContain("'X_ray_(R)1'");
        expect(context.convergence.quantity).toContain('not a total');
        expect(context.convergence.n_steps).toBe(GTS_ROWS);
    });
});
