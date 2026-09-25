// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { plotDataFromText, rwpColumns } from '../browserData';

// The dashboard R-factor must be normalized by the EXPERIMENT. RMCProfile writes
// fit CSVs as (x, calculated, experimental); with the calculated curve a uniform
// 0.7 × the experiment the conventional R = ‖calc − expt‖ / ‖expt‖ is exactly 0.3,
// while normalizing by the calculated column gives 0.3 / 0.7 ≈ 0.4286.
// Same cases as tests/test_plots.py::RwpColumnRoleTests.
const EXPT = [1.0, -2.0, 3.0, 0.5, -1.5];

const rwpOf = (plotKind, name, header, calcFirst = true) => {
    const rows = EXPT.map((expt, index) => {
        const calc = 0.7 * expt;
        const x = 0.1 * index;
        return (calcFirst ? [x, calc, expt] : [x, expt, calc]).map((value) => value.toFixed(7)).join(', ');
    });
    return plotDataFromText({ plotKind, name, text: [header, ...rows].join('\n') }).metrics.rwp;
};

describe('Rwp column roles', () => {
    it.each([
        ['xray_sq', 'run_FQ1.csv', 'Q, F(Q)_RMC, F(Q)_Expt'],
        ['xpdf', 'run_FT_XFQ1.csv', 'r(A), X_ray-calc, X_ray_exp_renorm'],
        ['npdf', 'run_PDF1.csv', 'r, G(r)_RMC, G(r)_Expt'],
        ['neutron_sq', 'run_SQ1.csv', 'Q, S(Q)_RMC, S(Q)_Expt'],
        ['bragg', 'run_bragg.csv', 'Flight time (us), Calculated, Experiment'],
    ])('%s divides by the experimental column', (kind, name, header) => {
        expect(rwpOf(kind, name, header)).toBeCloseTo(0.3, 12);
    });

    it('follows the RMCProfile positional order when the header names no roles', () => {
        expect(rwpOf('xray_sq', 'run_FQ1.csv', 'Q, a, b')).toBeCloseTo(0.3, 12);
    });

    it('lets header roles override the positional order', () => {
        expect(rwpOf('xray_sq', 'run_FQ1.csv', 'Q, F(Q)_Expt, F(Q)_RMC', false)).toBeCloseTo(0.3, 12);
    });

    it('resolves [calculated, experimental] indices like parsers.rwp_columns', () => {
        expect(rwpColumns(['Q', 'F(Q)_RMC', 'F(Q)_Expt'])).toEqual([1, 2]);
        expect(rwpColumns(['r(A)', 'X_ray-calc', 'X_ray_exp_renorm'])).toEqual([1, 2]);
        expect(rwpColumns(['Q', 'F(Q)_Expt', 'F(Q)_RMC'])).toEqual([2, 1]);
        expect(rwpColumns(['Q', 'observed', 'fitted'])).toEqual([2, 1]);
        expect(rwpColumns(['Q', 'a', 'b'])).toEqual([1, 2]);
        expect(rwpColumns(['Q', 'F(Q)_RMC'])).toBeNull();
    });
});
