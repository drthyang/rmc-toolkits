// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The convergence badge judges every chi^2 term of the .log and their total
// when the log names its columns (the rule of 2026-10-24, the Copilot's);
// classified from ONE fit term's chi^2 (the last .log column, a log without
// column names), it names that term rather than implying the whole run.

import { describe, expect, it, vi } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';

vi.mock('../settings', () => ({ useLlmSettings: () => ({ watchdogEnabled: true }) }));

import WatchdogBadge from '../components/WatchdogBadge';

const rValueFile = (chiColumn) => ({
    path: 'run-01.log',
    plotData: { chiColumn, series: [{ label: chiColumn, x: [0, 1, 2, 3], y: [5, 4, 3, 2] }] }
});

describe('WatchdogBadge', () => {
    it('prefixes the status with the chi^2 column it was classified from', () => {
        const html = renderToStaticMarkup(<WatchdogBadge rValueFile={rValueFile('X_ray_(R)1')} />);
        expect(html).toContain('X_ray_(R)1: ');
        expect(html).toContain('not the total');
    });

    it('shows the bare status when the column is unnamed', () => {
        const html = renderToStaticMarkup(<WatchdogBadge rValueFile={rValueFile(null)} />);
        expect(html).not.toContain(': ');
        expect(html).toContain('role="status"');
    });

    it('says it judged every χ² term and their total when the log names them', () => {
        const file = rValueFile('X_ray_(R)1');
        file.plotData.chiTerms = { names: ['F(Q)_1', 'X_ray_(R)1'], series: [[1, 1, 1, 1], [5, 4, 3, 2]], total: [5.1, 4.1, 3.1, 2.1], segments: [4] };
        const html = renderToStaticMarkup(<WatchdogBadge rValueFile={file} />);
        expect(html).not.toContain('X_ray_(R)1: ');
        expect(html).toContain('every χ² term of the .log and their total');
    });
});
