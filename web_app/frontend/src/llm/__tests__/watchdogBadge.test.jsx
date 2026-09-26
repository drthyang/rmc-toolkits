// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The convergence badge is classified from ONE fit term's chi^2 (the last
// .log column); it names that term rather than implying the whole run.

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
});
