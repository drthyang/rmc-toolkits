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

const state = vi.hoisted(() => ({ requests: [], plots: {}, plotProps: {}, cell: null, partials: true, zeroAngles: false }));

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
    // No triplet inside the windows: the engine's all-zero curves (not NaN),
    // no mean angle, every central atom 0-fold.
    if (state.zeroAngles) {
        return {
            triplet: [params.end1, params.apex, params.end2],
            bond12,
            bond23,
            sharedEnds,
            binWidth: 60,
            binCenters: [30, 90, 150],
            sinCorrected: [0, 0, 0],
            density: [0, 0, 0],
            coordination: [16],
            apexCount: 16,
            lengths12: { uniqueBonds: 0, count: 0, meanLength: null },
            lengths23: sharedEnds ? null : { uniqueBonds: 0, count: 0, meanLength: null },
            angleCount: 0,
            meanAngle: null,
            stdAngle: null,
        };
    }
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
    default: ({ file, plotData, legend, actionsTarget }) => {
        state.plots[file.path.split(':')[1]] = plotData;
        state.plotProps[file.path.split(':')[1]] = { legend, actionsTarget };
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
        state.plotProps = {};
        state.cell = null;
        state.partials = true;
        state.zeroAngles = false;
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

    it('draws the angles on a fixed 0–180° axis, with no reference line', async () => {
        await render();
        // Before Compute the same axis is there, empty.
        const ghost = state.plots.angles;
        expect(ghost).toMatchObject({ xDomain: [0, 180], xTicks: [0, 30, 60, 90, 120, 150, 180], xMinorStep: 10, xGrid: true, yMin: 0 });
        expect(ghost.series).toEqual([]);

        await submit();
        const plot = state.plots.angles;
        expect(plot).toMatchObject({ xDomain: [0, 180], xGrid: true, yMin: 0, xLabel: 'angle at Nb, θ (°)', yLabel: 'sin-corrected (random = 1)' });
        expect(plot.series).toHaveLength(1);
        expect(plot.series[0]).toMatchObject({ y: [0.5, 2, 0.5], curve: 'step', binWidth: 60, fill: true, width: 1.75 });

        // Density view: the same single curve in its own units.
        await act(async () => {
            [...container.querySelectorAll('button')].find((button) => button.textContent.trim() === 'density').click();
        });
        const density = state.plots.angles;
        expect(density.yLabel).toBe('density (deg^{-1})');
        expect(density.series).toHaveLength(1);
        expect(density.series[0]).toMatchObject({ label: 'density', curve: 'step', binWidth: 60 });
        expect(density.series.some((series) => series.role === 'guide')).toBe(false);
    });

    it('keeps the KPI rail in place: "—" before Compute, the results after', async () => {
        await render();
        // The live region wraps the description list; the <dl> keeps its own role.
        const rail = container.querySelector('[aria-label="Triplet result"]');
        expect([rail.tagName, rail.getAttribute('role'), rail.getAttribute('aria-live')]).toEqual(['DIV', 'status', 'polite']);
        expect(rail.querySelector('dl').hasAttribute('role')).toBe(false);
        expect(kpis().map((tile) => tile.value)).toEqual(['—', '—', '—']);
        expect(kpis().map((tile) => tile.label)).toEqual(['Angles', 'Coordination', 'Nb–Se bond']);
        await submit();
        const [angles, coordination, bond] = kpis();
        expect(angles.value).toBe('15.0 per Nb');
        // The mean ± std angle is visible on the sub line (not only in a hover);
        // the hover adds the bins and where it ran.
        expect(angles.sub).toBe('240 · mean 90.0 ± 5.0°');
        expect(container.querySelector('[aria-label="Triplet result"] .ui-kpi').title)
            .toBe('240 angles, mean 90.0° ± 5.0° (std); 60.0° bins; computed in the server.');
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

    it('a result with no angles keeps the empty axis and says so, instead of a zero curve', async () => {
        state.zeroAngles = true;
        await render();
        await setSelect('End element A', 'Ga');
        await submit();
        // No series at all: the bare axis.
        const plot = state.plots.angles;
        expect(plot).toMatchObject({ xDomain: [0, 180], yMin: 0, xLabel: 'angle at Nb, θ (°)' });
        expect(plot.series).toEqual([]);
        // One short line over the (dimmed, inert) axis, with the windows used.
        const prompt = container.querySelector('.geom-hero .ui-prompt');
        expect(prompt.querySelector('p').textContent).toBe('No Ga–Nb–Se triplets in these windows.');
        expect(prompt.querySelector('.ui-chip').textContent).toBe('2.00–3.00 Å');
        expect(prompt.querySelector('button')).toBeNull();
        expect(container.querySelector('.geom-plot__frame').className).toContain('ui-dim');
        // The KPIs still report what was found.
        expect(kpis()[0]).toMatchObject({ value: '0.0\u2009per Nb', sub: '0 angles' });
        // A result with angles draws its curve again.
        state.zeroAngles = false;
        await submit();
        expect(state.plots.angles.series[0]).toMatchObject({ curve: 'step', y: [0.5, 2, 0.5] });
        expect(container.querySelector('.geom-hero .ui-prompt')).toBeNull();
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
        // One header chip names both windows, each led by its bond-role dash
        // (whose hidden text names the role for a screen reader).
        const chips = container.querySelectorAll('.geom-helper h3 .ui-chip');
        expect(chips).toHaveLength(1);
        expect(chips[0].textContent).toBe('A–B 2.00–3.00\u2003B–C 2.00–3.00 Å');
        expect([...chips[0].querySelectorAll('.ui-bond-dash')].map((dash) => dash.style.getPropertyValue('--bond')))
            .toEqual([BOND_COLORS.ab, BOND_COLORS.bc]);

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
        // The g(r) title names its curves, so the plot drops its legend row and
        // its Save goes to the card header; the hero keeps its legend.
        expect(state.plotProps.partial.legend).toBe(false);
        expect(state.plotProps.partial.actionsTarget).toBe(container.querySelector('.geom-helper h3 .ui-card__cluster'));
        expect(state.plotProps.angles).toEqual({ legend: undefined, actionsTarget: undefined });
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
