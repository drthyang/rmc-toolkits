// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A plot payload with JSON-null gaps (Flask's strict-JSON encoding of a masked
// region) and a missing axis label must render — the old AxisLabel called
// label.match() on undefined and threw, unmounting the dashboard.

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import InteractivePlot from '../components/InteractivePlot';

describe('InteractivePlot with null gaps', () => {
    it('renders the finite points and skips the gaps', () => {
        const plotData = JSON.parse(
            '{"kind":"xray_sq","title":"F(Q)","xLabel":"Q (Å^{-1})","yLabel":"F(Q)",'
            + '"series":[{"label":"F(Q)_RMC","x":[1,2,3,4],"y":[0.1,null,null,0.4]},'
            + '{"label":"F(Q)_Expt","x":[1,2,3,4],"y":[0.2,0.2,0.3,0.4]}]}'
        );
        const html = renderToStaticMarkup(<InteractivePlot file={{ path: 'run_FQ1.csv', name: 'run_FQ1.csv' }} plotData={plotData} />);
        expect(html).toContain('<svg');
        expect(html).not.toContain('NaN');
    });

    it('does not throw when the payload lacks axis labels', () => {
        const plotData = { title: 't', series: [{ label: 'a', x: [0, 1], y: [1, 2] }] };
        expect(() => renderToStaticMarkup(
            <InteractivePlot file={{ path: 'a.csv', name: 'a.csv' }} plotData={plotData} />
        )).not.toThrow();
    });
});
