// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The run module's browser adapter (useCopilotRun.js) at its interface: what
// App hands CopilotChat. On the Demo, once the Copilot is open and the
// Dashboard has published, it equals the pure core (copilotRun.js) over the
// same run; before that it has the run's name and no data; each field keeps
// its identity until its own inputs change; a Flask-shaped run has no
// symmetry, no settings and no detect_symmetry tool; with no run open in the
// static build there is no engine. No worker starts: no tool runs here.

import { afterEach, beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import React, { act } from 'react';
import { createRoot } from 'react-dom/client';
import { copilotRun, DEFAULT_VIEW } from '../copilotRun';
import { loadDemoAssistantRun } from './modelVisible/helpers';

const mode = vi.hoisted(() => ({ static: false }));
vi.mock('../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => mode.static
}));

const { default: useCopilotRun } = await import('../useCopilotRun');

const definitions = (tools) => tools.map(({ name, description, parameters }) => ({ name, description, parameters }));

const settle = () => act(async () => {
    for (let i = 0; i < 10; i += 1) await new Promise((resolve) => setTimeout(resolve, 0));
});

describe('useCopilotRun (the browser adapter)', () => {
    let demo;
    let published;
    let container;
    let root;
    let hook;
    const Probe = (props) => {
        hook = useCopilotRun(props);
        return null;
    };
    const render = (props) => act(async () => { root.render(<Probe {...props} />); });

    beforeAll(async () => {
        demo = await loadDemoAssistantRun();
        // What the Dashboard publishes once everything is parsed.
        published = {
            files: demo.loaded.files,
            structure: demo.loaded.structure,
            structurePath: demo.loaded.structurePath,
            hiddenPlots: []
        };
    });

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        mode.static = false;
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
    });

    const demoProps = (overrides = {}) => ({
        directory: 'Demo',
        localRun: demo.localRun,
        dashboard: published,
        view: { toleranceA: DEFAULT_VIEW.toleranceA, pcaSites: null, liveData: false },
        enabled: true,
        ...overrides
    });

    it('before the Copilot is opened: the run\'s name, no data, and the engine tools', async () => {
        await render(demoProps({ enabled: false }));
        await settle();
        expect(hook.props).toEqual({
            runName: 'Demo', plotFiles: [], rValueFile: null, structure: null, symmetry: null, runSettings: null, pcaSites: null, liveData: false
        });
        expect(hook.engine).toEqual(expect.any(Function));
        expect(hook.tools.map((tool) => tool.name)).toEqual([
            'list_sites', 'analyze_displacement_directions', 'compute_bond_angles', 'analyze_clusters', 'analyze_modulations'
        ]);
        // The view's own props pass through from the start, as App always passed them.
        await render(demoProps({ enabled: false, view: { pcaSites: { sites: [] }, liveData: true } }));
        expect(hook.props).toMatchObject({ pcaSites: { sites: [] }, liveData: true, structure: null });
    });

    it('once open and published: deep-equals the core over the same run', async () => {
        await render(demoProps());
        await settle();
        const core = copilotRun(demo.loaded, DEFAULT_VIEW);
        expect(hook.props).toEqual(core.props);
        expect(hook.props.runSettings).not.toBeNull();
        expect(hook.elements).toEqual(core.elements);
        expect(hook.view).toEqual(core.view);
        expect(definitions(hook.tools)).toEqual(definitions(core.tools));
        expect(Object.isFrozen(hook.props)).toBe(true);
    });

    it('a new PCA table replaces props, never the tools, the engine or the other fields', async () => {
        await render(demoProps());
        await settle();
        const before = hook;
        const table = { sites: [{ referenceNumber: 1 }] };
        await render(demoProps({ view: { toleranceA: DEFAULT_VIEW.toleranceA, pcaSites: table, liveData: false } }));
        expect(hook.props).not.toBe(before.props);
        expect(hook.props.pcaSites).toBe(table);
        expect(hook.tools).toBe(before.tools);
        expect(hook.engine).toBe(before.engine);
        expect(hook.elements).toBe(before.elements);
        for (const key of ['plotFiles', 'rValueFile', 'structure', 'symmetry', 'runSettings']) {
            expect(hook.props[key]).toBe(before.props[key]);
        }
        expect(hook.view.hiddenPlots).toBe(before.view.hiddenPlots);
    });

    it('a tolerance changes the symmetry and keeps the ladder; a hidden log drops out of the history only', async () => {
        await render(demoProps());
        await settle();
        const before = hook;
        const rung = before.props.symmetry.ladder.find((brick) => brick.spaceGroup && brick.spaceGroup !== before.props.symmetry.spaceGroup);
        const tolerance = (rung.from + rung.to) / 2;
        await render(demoProps({ view: { toleranceA: tolerance, pcaSites: null, liveData: false } }));
        expect(hook.props.symmetry).not.toBe(before.props.symmetry);
        expect(hook.props.symmetry).toMatchObject({ spaceGroup: rung.spaceGroup, toleranceA: tolerance });
        expect(hook.props.symmetry.ladder).toBe(before.props.symmetry.ladder);
        expect(hook.props.rValueFile).toBe(before.props.rValueFile);
        expect(hook.tools).toBe(before.tools);

        const afterTolerance = hook;
        await render(demoProps({
            view: { toleranceA: tolerance, pcaSites: null, liveData: false },
            dashboard: { ...published, hiddenPlots: ['Demo/GTS_250K-01.log'] }
        }));
        expect(hook.props.rValueFile.sourceNames).toEqual(['GTS_250K-00.log', 'GTS_250K-02.log']);
        expect(hook.view.hiddenPlots).toEqual(['Demo/GTS_250K-01.log']);
        expect(hook.props.symmetry).toBe(afterTolerance.props.symmetry);
        expect(hook.props.plotFiles).toBe(afterTolerance.props.plotFiles);
        expect(hook.tools).toBe(afterTolerance.tools);
    });

    it('a Flask-shaped run: the Flask engine, no symmetry, no settings, no detect_symmetry tool', async () => {
        const { basis: _basis, ...flaskStructure } = demo.loaded.structure;
        await render({
            directory: '/runs/a',
            localRun: null,
            dashboard: {
                files: [{ name: 'run_PDFpartials.csv', path: '/runs/a/run_PDFpartials.csv', type: 'file', plotKind: 'pdf_partials' }],
                structure: flaskStructure,
                structurePath: '/runs/a/run.rmc6f',
                hiddenPlots: []
            },
            view: { toleranceA: DEFAULT_VIEW.toleranceA, pcaSites: null, liveData: true },
            enabled: true
        });
        await settle();
        expect(hook.props).toMatchObject({ runName: '/runs/a', symmetry: null, runSettings: null, liveData: true });
        expect(hook.engine).toEqual(expect.any(Function));
        const names = hook.tools.map((tool) => tool.name);
        expect(names).toContain('list_sites');
        expect(names).not.toContain('detect_symmetry');
    });

    it('the static build with no run open: named after the directory, no engine, no tools', async () => {
        mode.static = true;
        await render({ directory: 'data', localRun: null, dashboard: null, view: {}, enabled: true });
        expect(hook.props.runName).toBe('data');
        expect(hook.engine).toBeNull();
        expect(hook.tools).toEqual([]);
        expect(hook.view).toEqual(DEFAULT_VIEW);
    });
});
