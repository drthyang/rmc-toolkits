// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// Bond Geometry Phase 1 (presentation): the controls are a form, the angle
// plot reads on a fixed 0–180° axis against the exact random-bonds line, the
// headline results sit in a KPI rail that keeps its place before Compute, a
// shown result that no longer matches the inputs says so, a validation error
// names and marks its field, and one colour system (BOND_COLORS) ties the 3D
// bonds, the window guides and the chips together.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import { BOND_COLORS } from '../../plotPalette';

const state = vi.hoisted(() => ({ requests: [], plots: {}, cell: null, partials: true }));

const SITES = {
    elements: ['Ga', 'Nb', 'Se'],
    supercell: [2, 2, 2],
    sites: [
        { referenceNumber: 1, element: 'Ga', count: 4, copiesPerCell: 4 },
        { referenceNumber: 2, element: 'Nb', count: 16, copiesPerCell: 16 },
        { referenceNumber: 3, element: 'Se', count: 32, copiesPerCell: 32 },
    ],
};

// Three 60° bins make the isotropic reference easy to check by hand.
const triplets = (params) => {
    const bond12 = [params.r12Min, params.r12Max];
    const bond23 = params.r23Min !== undefined ? [params.r23Min, params.r23Max] : bond12;
    const sharedEnds = params.end1 === params.end2 && bond12[0] === bond23[0] && bond12[1] === bond23[1];
    const lengths = { uniqueBonds: 96, count: 96, meanLength: 2.6 };
    return {
        triplet: [params.end1, params.apex, params.end2],
        bond12,
        bond23,
        sharedEnds,
        binWidth: 60,
        binCenters: [30, 90, 150],
        sinCorrected: [0.5, 2, 0.5],
        density: [0.002, 0.012, 0.002],
        coordination: [0, 0, 0, 0, 0, 0, 16],
        apexCount: 16,
        lengths12: lengths,
        lengths23: sharedEnds ? null : { ...lengths, meanLength: 2.7 },
        angleCount: 240,
        meanAngle: 90,
        stdAngle: 5,
    };
};

vi.mock('axios', () => ({
    default: {
        get: vi.fn(async (url, config) => {
            state.requests.push({ url, params: config?.params });
            if (url.endsWith('/api/pca/sites')) return { data: SITES };
            if (url.endsWith('/api/triplets')) return { data: triplets(config.params) };
            if (url.endsWith('/api/files')) {
                return { data: { files: state.partials ? [{ path: 'run/PDFpartials.csv', plotKind: 'pdf_partials' }] : [] } };
            }
            if (url.endsWith('/api/plot/data')) {
                const x = [1, 2, 3, 4, 5];
                return { data: { series: [{ label: 'Nb-Se', x, y: [0, 1, 4, 1, 0] }, { label: 'Ga-Nb', x, y: [0, 0, 1, 3, 1] }] } };
            }
            if (url.endsWith('/api/structure')) return { data: { atomCount: 52, elements: SITES.elements } };
            throw new Error(`unexpected ${url}`);
        }),
    },
}));

vi.mock('../../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => false,
}));

// Capture what the page hands its plots and its 3D panel.
vi.mock('../InteractivePlot', () => ({
    default: ({ file, plotData }) => {
        state.plots[file.path.split(':')[1]] = plotData;
        return null;
    },
}));
vi.mock('../FoldedCellPanel', () => ({
    default: (props) => {
        state.cell = props;
        return <div className="mock-cell">{props.title}</div>;
    },
}));
vi.mock('../ModelSummary', () => ({ default: () => null }));

const { default: BondGeometryPage } = await import('../BondGeometryPage');

describe('BondGeometryPage presentation (Phase 1)', () => {
    let container;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        state.requests = [];
        state.plots = {};
        state.cell = null;
        state.partials = true;
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
    });

    const render = async () => {
        await act(async () => {
            root.render(<BondGeometryPage directory="runs/a" localRun={null} />);
        });
    };
    const count = (suffix) => state.requests.filter((request) => request.url.endsWith(suffix)).length;
    const input = (label) => container.querySelector(`input[aria-label="${label}"]`);
    const select = (label) => container.querySelector(`select[aria-label="${label}"]`);
    const setSelect = async (label, value) => {
        await act(async () => {
            const element = select(label);
            element.value = value;
            element.dispatchEvent(new Event('change', { bubbles: true }));
        });
    };
    // React tracks an input's value through its own setter: go around it.
    const type = async (label, value) => {
        await act(async () => {
            const element = input(label);
            const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set;
            setter.call(element, value);
            element.dispatchEvent(new Event('input', { bubbles: true }));
        });
    };
    const submit = async () => {
        await act(async () => {
            container.querySelector('form').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true }));
        });
    };
    const runButton = () => container.querySelector('form button[type="submit"]');
    const kpis = () => [...container.querySelectorAll('[aria-label="Triplet result"] .ui-kpi')].map((tile) => ({
        label: tile.querySelector('dt').textContent,
        value: tile.querySelector('.ui-kpi__value').textContent,
        sub: tile.querySelector('.ui-kpi__sub').textContent,
    }));
    const heading = (selector) => container.querySelector(`${selector} h3 .ui-card__label`).textContent;

    it('the controls are a form: submitting it (Enter in a field) runs Compute', async () => {
        await render();
        const button = runButton();
        expect(button.closest('form')).toBe(container.querySelector('form.ui-controls'));
        expect(button.textContent).toBe('Compute');
        expect(button.className).toContain('ui-btn-primary--run');
        await submit();
        expect(count('/api/triplets')).toBe(1);
        expect(state.requests.at(-1).params).toMatchObject({ end1: 'Se', apex: 'Nb', end2: 'Se', r12Min: 2, r12Max: 3, binWidth: 1 });
    });

    it('draws the angles on a fixed 0–180° axis against the exact random-bonds line', async () => {
        await render();
        // Before Compute the same axis is there, with only the reference.
        const ghost = state.plots.angles;
        expect(ghost).toMatchObject({ xDomain: [0, 180], xTicks: [0, 30, 60, 90, 120, 150, 180], xMinorStep: 10, xGrid: true, yMin: 0 });
        expect(ghost.series).toHaveLength(1);
        expect(ghost.series[0]).toMatchObject({ label: 'random bonds', role: 'guide', x: [0, 180], y: [1, 1] });

        await submit();
        const plot = state.plots.angles;
        expect(plot).toMatchObject({ xDomain: [0, 180], xGrid: true, yMin: 0, xLabel: 'angle at Nb, θ (°)', yLabel: 'sin-corrected (random = 1)' });
        expect(plot.series[0]).toMatchObject({ y: [0.5, 2, 0.5], curve: 'step', binWidth: 60, fill: true, width: 1.75 });
        expect(plot.series[1]).toMatchObject({ label: 'random bonds', role: 'guide', y: [1, 1] });

        // Density view: the isotropic fraction of each bin per degree,
        // (cos θlo − cos θhi)/(2·w), which integrates to 1.
        await act(async () => {
            [...container.querySelectorAll('button')].find((button) => button.textContent.trim() === 'density').click();
        });
        const density = state.plots.angles;
        expect(density.yLabel).toBe('density (deg^{-1})');
        const reference = density.series[1];
        expect(reference).toMatchObject({ role: 'guide', curve: 'step', binWidth: 60 });
        expect(reference.y[0]).toBeCloseTo((1 - 0.5) / 120, 12);
        expect(reference.y[1]).toBeCloseTo((0.5 + 0.5) / 120, 12);
        expect(reference.y.reduce((acc, value) => acc + value * 60, 0)).toBeCloseTo(1, 12);
    });

    it('keeps the KPI rail in place: "—" before Compute, the results after', async () => {
        await render();
        expect(kpis().map((tile) => tile.value)).toEqual(['—', '—', '—']);
        expect(kpis().map((tile) => tile.label)).toEqual(['Angles', 'Coordination', 'Nb–Se bond']);
        await submit();
        const [angles, coordination, bond] = kpis();
        expect(angles.value).toBe('15.0 per Nb');
        expect(angles.sub).toBe('240 angles · 60.0° bins · server');
        expect(coordination.value).toBe('6.00 per Nb');
        expect(coordination.sub).toBe('6-fold 100.0% · of 16');
        expect(bond.value).toBe('2.600 Å');
        expect(bond.sub).toBe('96 bonds · 2.00–3.00 Å');
        // Two bond types: a second bond tile, in the B–C colour.
        await setSelect('End element A', 'Ga');
        await submit();
        expect(kpis().map((tile) => tile.label)).toEqual(['Angles', 'Coordination', 'Nb–Ga bond', 'Nb–Se bond']);
        expect(kpis()[3].value).toBe('2.700 Å');
    });

    it('marks a shown result whose inputs changed, until they match again', async () => {
        await render();
        await submit();
        expect(runButton().textContent).toBe('Compute');
        await type('A-B window maximum', '2.9');
        expect(runButton().textContent).toBe('Update');
        expect(runButton().className).toContain('is-stale');
        expect(container.textContent).toContain('inputs changed');
        // Same number, other spelling: not a change.
        await type('A-B window maximum', '3.0');
        expect(runButton().textContent).toBe('Compute');
        expect(container.textContent).not.toContain('inputs changed');
    });

    it('a validation error names, marks and focuses its field, in the hero prompt', async () => {
        await render();
        await type('A-B window minimum', '');
        await submit();
        expect(count('/api/triplets')).toBe(0);
        const message = container.querySelector('.geom-hero [role="alert"]');
        expect(message.textContent).toBe('A–B window minimum is empty — enter a number.');
        expect(input('A-B window minimum').getAttribute('aria-invalid')).toBe('true');
        expect(input('A-B window minimum').closest('.ui-unit-field').className).toContain('is-invalid');
        expect(document.activeElement).toBe(input('A-B window minimum'));
        // Fixing the field clears the mark and the message.
        await type('A-B window minimum', '2.1');
        expect(input('A-B window minimum').hasAttribute('aria-invalid')).toBe(false);
        expect(container.querySelector('.geom-hero [role="alert"]')).toBeNull();
    });

    it('one colour system: 3D bonds, window guides and curves wear BOND_COLORS', async () => {
        await render();
        // A lone window: neutral guides (no colour of their own), out of the legend.
        let guides = state.plots.partial.series.filter((series) => series.role === 'guide');
        expect(guides.map((guide) => guide.x[0])).toEqual([2, 3]);
        expect(guides.every((guide) => guide.legend === false && guide.color === undefined)).toBe(true);
        expect(state.plots.partial.series[0]).toMatchObject({ label: 'Nb-Se', color: BOND_COLORS.ab });

        await setSelect('End element A', 'Ga');
        await act(async () => { input('Use a distinct B-C window').click(); });
        await act(async () => { await new Promise((resolve) => setTimeout(resolve, 450)); });
        const curves = state.plots.partial.series.filter((series) => series.role !== 'guide');
        expect(curves.map((curve) => [curve.label, curve.color])).toEqual([['Ga-Nb', BOND_COLORS.ab], ['Nb-Se', BOND_COLORS.bc]]);
        guides = state.plots.partial.series.filter((series) => series.role === 'guide');
        expect(guides.map((guide) => guide.color)).toEqual([BOND_COLORS.ab, BOND_COLORS.ab, BOND_COLORS.bc, BOND_COLORS.bc]);

        await submit();
        expect(state.cell.bondSets.map((set) => [set.elements.join('-'), set.color]))
            .toEqual([['Ga-Nb', BOND_COLORS.ab], ['Nb-Se', BOND_COLORS.bc]]);
        expect(state.cell.legendEmphasis).toEqual(['Ga', 'Nb', 'Se']);
    });

    it('names the cards for the triplet, as element chips that read as text', async () => {
        await render();
        expect(heading('.geom-hero')).toMatch(/^Se–Nb–Se bond angles/);
        expect(heading('.geom-helper')).toMatch(/^Nb–Se partial g\(r\)/);
        expect(container.querySelector('.mock-cell').textContent).toBe('Nb–Se bonds');
        // The central atom is ringed; the window chip follows the inputs.
        expect(container.querySelector('.geom-hero h3 .ui-element-chip--central').textContent).toBe('Nb');
        expect(container.querySelector('.geom-helper h3 .ui-chip').textContent).toBe('2.00–3.00 Å');
        // Before Compute the 3D legend says what Compute adds.
        expect(state.cell.bondSets).toBeNull();
    });

    it('swaps the ends, and offers the swap only when they differ', async () => {
        await render();
        const swap = () => container.querySelector('button[aria-label="Swap A and C"]');
        expect(swap()).toBeNull();
        await setSelect('End element A', 'Ga');
        await act(async () => { swap().click(); });
        expect(select('End element A').value).toBe('Se');
        expect(select('End element C').value).toBe('Ga');
    });

    it('without a partials file the g(r) card is a slim row and the 3D card takes the height', async () => {
        state.partials = false;
        await render();
        expect(container.querySelector('.geom-layout').className).toContain('geom-layout--no-helper');
        const helper = container.querySelector('.geom-helper');
        expect(helper.children).toHaveLength(1);
        expect(helper.textContent).toContain('No PDFpartials.csv in this run.');
    });
});
