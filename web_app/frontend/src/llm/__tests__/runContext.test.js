import { describe, expect, it } from 'vitest';
import {
    CONTEXT_CHAR_BUDGET,
    buildRunContext,
    contextToJson,
    downsampleSeries,
    recentSlope,
    seriesStats
} from '../context/runContext';

// A fixture mimicking the props Dashboard.jsx passes: parsed plot files, a
// combined R-value file, a folded structure, and a detected symmetry.
const historyOf = (length, fn) => Array.from({ length }, (_, index) => fn(index));

const fixtureProps = () => ({
    runName: '5K_try1',
    liveData: true,
    plotFiles: [
        { plotKind: 'neutron_sq', plotData: { title: 'S(Q) (neutron)', metrics: { rwp: 0.041234 } } },
        { plotKind: 'npdf', plotData: { title: 'PDF1', metrics: { rwp: 0.05331, pdf_index: 1 } } },
        { plotKind: 'stog', plotData: { title: 'scale_ft.gr', metrics: {} } },
        { plotKind: 'r_value', plotData: { title: 'R-value', metrics: { final_chi_r: 2.38 } } }
    ],
    rValueFile: {
        plotData: {
            metrics: { final_chi_r: 2.38 },
            series: [{ y: historyOf(500, (index) => 2.3 - index * 0.003) }]
        }
    },
    structure: {
        latticeVectors: [[60.12, 0, 0], [0, 60.12, 0], [0, 0, 60.12]],
        supercell: [6, 6, 6],
        totalAtoms: 12288,
        elementCounts: { Ga: 864, Nb: 3456, Se: 6912 },
        basis: [
            { el: 'Ga', frac: [0, 0, 0], dispA: 0.05 },
            { el: 'Nb', frac: [0.6, 0.6, 0.6], dispA: 0.25 },
            { el: 'Nb', frac: [0.4, 0.4, 0.4], dispA: 0.15 },
            { el: 'Se', frac: [0.37, 0.37, 0.37], dispA: 0.1 }
        ]
    },
    // The shape Dashboard passes: describeSymmetry() output + toleranceA + ladder.
    symmetry: {
        spaceGroup: 'F-43m',
        spaceGroupNumber: 216,
        pointGroup: '-43m',
        nSpace: 96,
        maxResidual: 0.142,
        toleranceA: 0.2,
        ladder: [
            { spaceGroup: 'P1', from: 0, to: 0.06, nSpace: 1 },
            { spaceGroup: 'P-42m', from: 0.06, to: 0.15, nSpace: 8 },
            { spaceGroup: 'F-43m', from: 0.15, to: 1, nSpace: 96 }
        ],
        orbits: [
            { element: 'Ga', size: 4, site: '-43m', rep: [0, 0, 0], wyckoff: 'a', members: [0] },
            { element: 'Nb', size: 16, site: '3m', rep: [0.6, 0.6, 0.6], wyckoff: 'e', members: [1, 2] },
            { element: 'Se', size: 16, site: '3m', rep: [0.37, 0.37, 0.37], wyckoff: 'e', members: [3] }
        ]
    }
});

describe('downsampleSeries', () => {
    it('keeps the first and last points and respects the cap', () => {
        const values = historyOf(1000, (index) => index * 0.5);
        const sampled = downsampleSeries(values, 48);
        expect(sampled).toHaveLength(48);
        expect(sampled[0]).toBe(0);
        expect(sampled[47]).toBe(Number((999 * 0.5).toPrecision(3)));
    });

    it('returns short series unchanged in length', () => {
        expect(downsampleSeries([1, 2, 3], 48)).toHaveLength(3);
        expect(downsampleSeries([], 48)).toEqual([]);
    });

    it('rounds to 3 significant digits', () => {
        expect(downsampleSeries([1.23456, 2.34567], 48)).toEqual([1.23, 2.35]);
    });
});

describe('seriesStats and recentSlope', () => {
    it('computes stats on the full series', () => {
        const values = historyOf(100, (index) => 100 - index);
        const stats = seriesStats(values);
        expect(stats.nSteps).toBe(100);
        expect(stats.first).toBe(100);
        expect(stats.last).toBe(1);
        expect(stats.min).toBe(1);
        expect(stats.max).toBe(100);
        expect(stats.recentSlopePerStep).toBeCloseTo(-1, 6);
    });

    it('handles empty and constant input', () => {
        expect(seriesStats([])).toBeNull();
        expect(recentSlope([5, 5, 5, 5, 5, 5])).toBe(0);
        expect(recentSlope([5])).toBe(0);
    });
});

describe('buildRunContext', () => {
    it('assembles every section from full props', () => {
        const context = buildRunContext(fixtureProps());
        expect(context.run).toBe('5K_try1');
        expect(context.live_mode).toBe(true);
        expect(context.structure.cell_A).toEqual([10.02, 10.02, 10.02]);
        expect(context.structure.angles_deg).toEqual([90, 90, 90]);
        expect(context.structure.elements.Se).toBe(6912);
        // r_value files are convergence data, not fitted datasets.
        expect(context.datasets).toHaveLength(3);
        expect(context.datasets[0]).toEqual({ kind: 'neutron_sq', title: 'S(Q) (neutron)', rwp: 0.0412 });
        expect(context.convergence.n_steps).toBe(500);
        expect(context.convergence.final_chi_squared).toBe(2.38);
        expect(context.convergence.history.length).toBeLessThanOrEqual(48);
        expect(context.convergence.quantity).toContain('ln');
    });

    it('builds the symmetry block: group, ladder, and displacement-ranked sites', () => {
        const context = buildRunContext(fixtureProps());
        const symmetry = context.symmetry;
        expect(context.structure.space_group).toBeUndefined();   // moved into symmetry
        expect(symmetry.space_group).toBe('F-43m (No. 216)');
        expect(symmetry.point_group).toBe('-43m');
        expect(symmetry.n_ops).toBe(96);
        expect(symmetry.tolerance_A).toBe(0.2);
        expect(symmetry.max_residual_A).toBe(0.14);
        expect(symmetry.ladder).toEqual([
            { sg: 'P1', holds_A: [0, 0.06], n_ops: 1 },
            { sg: 'P-42m', holds_A: [0.06, 0.15], n_ops: 8 },
            { sg: 'F-43m', holds_A: [0.15, 1], n_ops: 96 }
        ]);
        // Sites sorted by mean displacement, largest first — Nb leads.
        expect(symmetry.sites.map((site) => site.element)).toEqual(['Nb', 'Se', 'Ga']);
        expect(symmetry.sites[0]).toEqual({
            element: 'Nb',
            multiplicity: 16,
            wyckoff: '16e',
            site_sym: '3m',
            frac: [0.6, 0.6, 0.6],
            mean_disp_A: 0.2,     // (0.25 + 0.15) / 2 over the orbit's basis members
            max_disp_A: 0.25
        });
    });

    it('builds a minimal symmetry block from just a space group (Flask mode)', () => {
        const props = fixtureProps();
        props.symmetry = { spaceGroup: 'F-43m', spaceGroupNumber: 216 };
        const symmetry = buildRunContext(props).symmetry;
        expect(symmetry.space_group).toBe('F-43m (No. 216)');
        expect(symmetry.ladder).toBeUndefined();
        expect(symmetry.sites).toBeUndefined();
    });

    it('omits sections whose inputs are missing', () => {
        const context = buildRunContext({ runName: 'bare' });
        expect(context.run).toBe('bare');
        expect(context.structure).toBeUndefined();
        expect(context.symmetry).toBeUndefined();
        expect(context.pair_correlations).toBeUndefined();
        expect(context.datasets).toBeUndefined();
        expect(context.convergence).toBeUndefined();
    });

    it('tolerates a structure with no lattice (known empty-.rmc6f case)', () => {
        const context = buildRunContext({
            runName: 'partial',
            structure: { totalAtoms: 0, elementCounts: {} }
        });
        expect(context.structure).toBeUndefined();
    });
});

describe('contextToJson budget', () => {
    it('stays within the character budget for a large run', () => {
        const props = fixtureProps();
        props.plotFiles = historyOf(30, (index) => ({
            plotKind: 'npdf',
            plotData: { title: `PDF${index}`, metrics: { rwp: 0.05 + index * 0.001 } }
        }));
        props.rValueFile.plotData.series[0].y = historyOf(50000, (index) => 2.3 - index * 0.00004);
        const json = contextToJson(buildRunContext(props));
        expect(json.length).toBeLessThanOrEqual(CONTEXT_CHAR_BUDGET);
        // Shrinking the history was enough here; datasets survive intact.
        expect(JSON.parse(json).datasets).toHaveLength(30);
    });

    it('truncates the dataset list as a last resort, recording the omission', () => {
        const props = fixtureProps();
        props.plotFiles = historyOf(30, (index) => ({
            plotKind: 'npdf',
            plotData: { title: `PDF${index}`, metrics: { rwp: 0.05 + index * 0.001 } }
        }));
        const parsed = JSON.parse(contextToJson(buildRunContext(props), 1500));
        expect(parsed.datasets).toHaveLength(12);
        expect(parsed.datasets_omitted).toBe(18);
    });

    it('leaves small contexts untouched', () => {
        const context = buildRunContext({ runName: 'tiny' });
        expect(JSON.parse(contextToJson(context))).toEqual(context);
    });

    it('trims evidence in order — peaks, ladder rungs, history, sites, datasets — recording omissions', () => {
        const context = {
            symmetry: {
                ladder: historyOf(5, (index) => ({ sg: `SG${index}`, holds_A: [index, index + 1], n_ops: index + 1 })),
                sites: historyOf(10, (index) => ({ element: 'X', multiplicity: 4, mean_disp_A: 1 - index * 0.05 }))
            },
            pair_correlations: [
                { pair: 'A-B', gr_peaks_A: [{ r: 2, fwhm: 0.2 }, { r: 3, fwhm: 0.3 }] }
            ],
            convergence: { history: historyOf(48, (index) => index) },
            datasets: historyOf(20, (index) => ({ kind: 'npdf', title: `PDF${index}` }))
        };
        const parsed = JSON.parse(contextToJson(context, 10));   // force every trim step
        expect(parsed.pair_correlations[0].gr_peaks_A).toHaveLength(1);
        expect(parsed.symmetry.ladder).toHaveLength(2);
        expect(parsed.symmetry.ladder_rungs_omitted).toBe(3);
        expect(parsed.symmetry.ladder[1].sg).toBe('SG4');   // endpoints kept
        expect(parsed.convergence.history).toHaveLength(12);
        expect(parsed.symmetry.sites).toHaveLength(8);
        expect(parsed.symmetry.sites_omitted).toBe(2);
        expect(parsed.datasets).toHaveLength(12);
        expect(parsed.datasets_omitted).toBe(8);
    });
});

describe('pair correlations integration', () => {
    it('includes pair correlations when partials are parsed', () => {
        const props = fixtureProps();
        const x = historyOf(300, (index) => index * 0.02);
        props.plotFiles.push({
            plotKind: 'pdf_partials',
            plotData: {
                title: 'PDFpartials',
                metrics: {},
                series: [
                    { label: 'Nb-Nb', x, y: x.map((r) => Math.exp(-((r - 3) ** 2) / 0.02)) }
                ]
            }
        });
        const context = buildRunContext(props);
        expect(context.pair_correlations).toHaveLength(1);
        expect(context.pair_correlations[0].pair).toBe('Nb-Nb');
        expect(context.pair_correlations[0].gr_peaks_A[0].r).toBeCloseTo(3, 1);
        // Fixture basis: Nb at (.6,.6,.6) and (.4,.4,.4) → 0.2·√3 · 10.02 Å ≈ 3.47 Å.
        expect(context.pair_correlations[0].avg_structure_d_A).toBeCloseTo(3.47, 1);
    });
});
