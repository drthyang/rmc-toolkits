// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The G(r) panel's error ribbon on the Dashboard (GrRibbon.jsx): drawn when the
// G(r) file has no σ and F(Q)'s σ propagates (a band in the experiment's
// colour, a legend entry, ± σ in the readout), with a caption that says the
// points are correlated and a ? help holding the reproduction residual and its
// tolerance; with the file's own σ, that one, and the propagated one as a
// toggle; no ribbon and the reason when F(Q) carries no σ; a failure also in
// the Problems section. Runs from ft_uncertainty_fixture.json, read as a
// picked folder is.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import fixture from '../../__tests__/fixtures/ft_uncertainty_fixture.json';

vi.mock('../ModelSummary', () => ({ default: () => null }));
vi.mock('../../llm', () => ({ WatchdogBadge: () => null }));

const { default: Dashboard } = await import('../Dashboard');
const { makeRunFromEntries } = await import('../../browserData');
const { clearFtUncertaintyCache } = await import('../../ftUncertaintyRun');
const { IssueScope, IssueStoreProvider } = await import('../../ui');

const filesOf = (name) => {
    const testCase = fixture.cases.find((entry) => entry.name === name);
    return Object.fromEntries(Object.entries(testCase.files).map(([file, key]) => [file, fixture.texts[key]]));
};

// The FT file with its own σ column (0.01 everywhere), the rest unchanged.
const withSigmaColumn = (files) => {
    const name = `${fixture.stem}_FT_XFQ1.csv`;
    const [header, ...rows] = files[name].trimEnd().split('\n');
    return { ...files, [name]: `${[`${header} , sigma`, ...rows.map((row) => `${row} , 0.01`)].join('\n')}\n` };
};

const runOf = async (files) => ({
    ...(await makeRunFromEntries(Object.entries(files).map(([name, text]) => ({
        path: `run/${name}`, file: new File([text], name, { lastModified: 0 })
    })))),
    runId: 1
});

const flush = () => act(async () => {
    for (let i = 0; i < 20; i += 1) await new Promise((resolve) => setTimeout(resolve, 0));
});

describe('the G(r) panel\'s error ribbon', () => {
    let container;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        clearFtUncertaintyCache();
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
    });

    const mount = async (files) => {
        const localRun = await runOf(files);
        act(() => root.render(
            <IssueStoreProvider>
                <IssueScope page="dashboard">
                    <Dashboard directory="run" localRun={localRun} />
                </IssueScope>
            </IssueStoreProvider>
        ));
        await flush();
    };

    const xpdfCard = () => container.querySelector('[data-copilot-card="plot:rmc_run_FT_XFQ1.csv"]');
    const legend = () => [...xpdfCard().querySelectorAll('.plot-legend button')].map((button) => button.textContent);
    const note = () => xpdfCard().querySelector('.gr-ribbon-note');
    const problems = () => [...container.querySelectorAll('.ui-issues .ui-issue')].map((row) => row.textContent);

    it('draws the propagated ribbon and says the points are correlated, with the tolerance in its help', async () => {
        await mount(filesOf('propagated, CONVOLVE'));
        const band = xpdfCard().querySelector('path.series-ribbon');
        expect(band).not.toBeNull();
        expect(band.getAttribute('d')).toMatch(/^M [\d.]+ [\d.]+ L/);
        expect(band.getAttribute('d')).toMatch(/Z$/);
        expect(legend()).toEqual(['X_ray_exp_renorm', 'X_ray-calc', '±1σ propagated']);
        expect(note().dataset.ribbon).toBe('propagated');
        expect(note().textContent).toContain('Error ribbon: ±1σ propagated from scale_ft_rmc.fq’s σ through RMCProfile’s transform');
        expect(note().textContent).toContain('the r points are correlated (ρ = 0.87 between neighbours)');
        const help = note().querySelector('[role="tooltip"]').textContent;
        expect(help).toContain('tolerance 1×10⁻⁶');
        expect(help).toMatch(/rmc_run_XFQ1\.csv → rmc_run_FT_XFQ1\.csv: largest deviation [\d.]+×10⁻[⁰¹²³⁴⁵⁶⁷⁸⁹]+ of max \|G\|/);
        expect(help).toContain('±1σ is the diagonal of the full covariance; the points are correlated');
        expect(help).toContain('The box convolution, left out of the propagation, changes σ by');
        expect(note().getAttribute('title')).toContain('tolerance 1×10⁻⁶');
        expect(problems()).toEqual([]);
        // The legend entry hides and shows the band.
        const entry = [...xpdfCard().querySelectorAll('.plot-legend button')].find((button) => button.textContent === '±1σ propagated');
        act(() => entry.click());
        expect(xpdfCard().querySelector('path.series-ribbon')).toBeNull();
        act(() => entry.click());
        expect(xpdfCard().querySelector('path.series-ribbon')).not.toBeNull();
    });

    it('no σ in F(Q): no ribbon, and the caption names what is missing', async () => {
        await mount(filesOf('input file missing'));
        expect(xpdfCard().querySelector('path.series-ribbon')).toBeNull();
        expect(legend()).toEqual(['X_ray_exp_renorm', 'X_ray-calc']);
        expect(note().textContent).toContain(
            'No error ribbon. F(Q) carries no uncertainties: scale_ft_rmc.fq (named in rmc_run.dat) is not in the run folder.'
        );
        // The transform is still shown to reproduce the run's G(r).
        expect(note().querySelector('[role="tooltip"]').textContent).toMatch(/largest deviation .* \(tolerance 1×10⁻⁶\) ✓/);
        expect(problems()).toEqual([]);

        act(() => root.unmount());
        root = createRoot(container);
        clearFtUncertaintyCache();
        await mount(filesOf('no error column'));
        expect(note().textContent).toContain('No error ribbon. F(Q) carries no uncertainties: scale_ft_rmc.fq has no error column.');
    });

    it('a failed reproduction: no ribbon, the reason on the card and in the Problems section', async () => {
        await mount(filesOf('G(r) by the trapezoid rule'));
        expect(xpdfCard().querySelector('path.series-ribbon')).toBeNull();
        expect(note().dataset.ribbon).toBe('failed');
        expect(note().textContent).toContain(
            'No error ribbon. RMCProfile\'s transform of rmc_run_XFQ1.csv does not reproduce rmc_run_FT_XFQ1.csv within the tolerance.'
        );
        expect(note().querySelector('[role="tooltip"]').textContent).toMatch(/\(tolerance 1×10⁻⁶\) ✗/);
        expect(problems()).toEqual([expect.stringContaining('No error ribbon. RMCProfile\'s transform of rmc_run_XFQ1.csv does not reproduce')]);
        expect(problems()[0]).toContain('rmc_run_FT_XFQ1.csv');
    });

    it('the file\'s own σ is drawn, and the propagated one is a toggle', async () => {
        await mount(withSigmaColumn(filesOf('propagated, CONVOLVE')));
        expect(legend()).toEqual(['X_ray_exp_renorm', 'X_ray-calc', '±1σ (file)']);
        expect(note().dataset.ribbon).toBe('file');
        expect(note().textContent).toContain('Error ribbon: ±1σ from the file’s own σ column.');
        const toggle = xpdfCard().querySelector('[aria-label="Error ribbon"]');
        expect([...toggle.querySelectorAll('button')].map((button) => [button.textContent, button.getAttribute('aria-pressed')]))
            .toEqual([['σ file', 'true'], ['σ propagated', 'false']]);
        act(() => toggle.querySelectorAll('button')[1].click());
        expect(legend()).toEqual(['X_ray_exp_renorm', 'X_ray-calc', '±1σ propagated']);
        expect(note().dataset.ribbon).toBe('propagated');
        expect(xpdfCard().querySelector('path.series-ribbon')).not.toBeNull();
    });

    it('Live Data: a changed file beside the transform redraws the ribbon in place', async () => {
        const files = filesOf('propagated, CONVOLVE');
        await mount(files);
        expect(note().dataset.ribbon).toBe('propagated');
        const before = xpdfCard().querySelector('path.series-ribbon');
        // The same run (runId 1), its input file rewritten without σ (a newer file).
        const changed = await makeRunFromEntries(Object.entries({ ...files, 'scale_ft_rmc.fq': filesOf('no error column')['scale_ft_rmc.fq'] })
            .map(([name, text]) => ({ path: `run/${name}`, file: new File([text], name, { lastModified: name === 'scale_ft_rmc.fq' ? 5 : 0 }) })));
        act(() => root.render(
            <IssueStoreProvider>
                <IssueScope page="dashboard">
                    <Dashboard directory="run" localRun={{ ...changed, runId: 1 }} />
                </IssueScope>
            </IssueStoreProvider>
        ));
        await flush();
        expect(before.isConnected).toBe(false);
        expect(xpdfCard().querySelector('path.series-ribbon')).toBeNull();
        expect(note().textContent).toContain('No error ribbon. F(Q) carries no uncertainties: scale_ft_rmc.fq has no error column.');
    });

    it('without the propagated one, a file σ needs no toggle and the help says why', async () => {
        await mount(withSigmaColumn(filesOf('input file missing')));
        expect(legend()).toEqual(['X_ray_exp_renorm', 'X_ray-calc', '±1σ (file)']);
        expect(xpdfCard().querySelector('[aria-label="Error ribbon"]')).toBeNull();
        expect(note().querySelector('[role="tooltip"]').textContent).toContain(
            'No propagated ribbon. F(Q) carries no uncertainties: scale_ft_rmc.fq (named in rmc_run.dat) is not in the run folder.'
        );
    });
});
