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
        elementCounts: { Ga: 864, Nb: 3456, Se: 6912 }
    },
    symmetry: { spaceGroup: 'F-43m', spaceGroupNumber: 216 }
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
        expect(context.structure.space_group).toBe('F-43m (No. 216)');
        expect(context.structure.elements.Se).toBe(6912);
        // r_value files are convergence data, not fitted datasets.
        expect(context.datasets).toHaveLength(3);
        expect(context.datasets[0]).toEqual({ kind: 'neutron_sq', title: 'S(Q) (neutron)', rwp: 0.0412 });
        expect(context.convergence.n_steps).toBe(500);
        expect(context.convergence.final_chi_squared).toBe(2.38);
        expect(context.convergence.history.length).toBeLessThanOrEqual(48);
        expect(context.convergence.quantity).toContain('ln');
    });

    it('omits sections whose inputs are missing', () => {
        const context = buildRunContext({ runName: 'bare' });
        expect(context.run).toBe('bare');
        expect(context.structure).toBeUndefined();
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
});
