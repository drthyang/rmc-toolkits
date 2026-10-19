// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A case as the harness opens it (runFolder.js loadRun → the run module's
// in-process adapter, src/openRunFolder.js): the run the Copilot is handed,
// the view each record keeps, the context sections it records (plan item
// 1.3), and the PCA-table option for a dev check (off by default). No model
// runs here.

import { describe, expect, it } from 'vitest';
import { buildRunContext } from '../../src/llm/context/runContext';
import { DEFAULT_VIEW } from '../../src/copilotRun';
import { openRunFolder } from '../../src/openRunFolder';
import { demoEntries } from '../../src/__tests__/modelVisible/helpers';
import { buildCase, DEV_CASES } from '../cases';
import { CONTEXT_SECTIONS, contextRecord, runLiveCase } from '../harness';
import { HARNESS_VIEW, loadRun, RUN_FILES, RUN_NAME } from '../runFolder';
import { renderReport } from '../report';
import { renderTranscript } from '../transcript';

const built = buildCase(DEV_CASES.find((item) => item.id === 'dev-offcentre'));
const MODEL = { name: 'scripted', label: 'scripted test model', model: 'scripted', baseUrl: 'none', temperature: 0 };
const ANSWER = 'No order.\n\nVERDICT {"features": [], "space_group": "Pm-3m", "chiral": false, "centrosymmetric": true, "polar": false}';
const answering = () => async function* call() {
    yield { content: ANSWER };
    yield { finishReason: 'stop' };
};

describe('a case opened through the run module', () => {
    it('is the run rmc_run with the default view, its partials detected as pdf_partials', async () => {
        expect(HARNESS_VIEW).toEqual(DEFAULT_VIEW);
        const run = await loadRun(built);
        expect(run.view).toEqual(DEFAULT_VIEW);
        expect(run.props.runName).toBe(RUN_NAME);
        expect(run.props.structure.source).toBe(`${RUN_NAME}/${RUN_FILES.rmc6f}`);
        // detectPlotKind now names the kind the harness used to hard-code.
        expect(run.props.plotFiles.map((file) => [file.name, file.path, file.plotKind])).toEqual([
            [RUN_FILES.partials, `${RUN_NAME}/${RUN_FILES.partials}`, 'pdf_partials']
        ]);
        expect(run.props.plotFiles[0].plotData.series.length).toBeGreaterThan(0);
        expect(run.props).toMatchObject({ rValueFile: null, runSettings: null, pcaSites: null, liveData: false });
        expect(run.props.symmetry).toMatchObject({ spaceGroup: 'Pm-3m', toleranceA: 0.2 });
        expect(run.props.symmetry.ladder.length).toBeGreaterThan(0);
        expect(run.elements).toEqual(Object.keys(run.props.structure.elementCounts).sort());
    });

    it('takes a partial view over the harness default (arm A: no ladder)', async () => {
        const run = await loadRun(built, { ladder: false });
        expect(run.view).toEqual({ ...DEFAULT_VIEW, ladder: false });
        expect(run.props.symmetry.ladder).toEqual([]);
    });

    it('the PCA-table option gives the rows list_sites reads, as the PCA Ellipsoid page would', async () => {
        const run = await loadRun(built, HARNESS_VIEW, { pcaTable: true });
        const sites = await run.engine('sites', {});
        expect(run.view.pcaSites).toEqual(sites);
        expect(run.props.pcaSites).toBe(run.view.pcaSites);
        expect(JSON.parse(JSON.stringify(run.view))).toEqual(run.view);
        expect(buildRunContext(run.props)).toHaveProperty('pca_displacements');
        expect(buildRunContext((await loadRun(built)).props)).not.toHaveProperty('pca_displacements');
    });
});

describe('the record states the view and the context it ran with (plan 1.3)', () => {
    it('by default: the default view, no PCA table, no run settings, one fresh question', async () => {
        const record = await runLiveCase({ built, arm: 'C', model: MODEL, callModel: answering() });
        expect(record.error).toBeNull();
        expect(record.view).toEqual(DEFAULT_VIEW);
        expect(record.context.sections).toEqual(['structure', 'symmetry', 'pair_correlations', 'datasets']);
        expect(record.context.absent).toEqual(['configuration_optimization', 'run_settings', 'pca_displacements', 'partial_weights', 'convergence']);
        expect(record.context.other_sections).toEqual([]);
        expect(record.context.earlier_turns).toBe(0);
        expect([...record.context.sections, ...record.context.absent].sort()).toEqual([...CONTEXT_SECTIONS].sort());
        const transcript = renderTranscript(record);
        expect(transcript).toContain('- View: Detected SG at 0.2 Å with its ladder, nothing hidden, no PCA table, Live Data off');
        expect(transcript).toMatch(/- Context sections: structure, .*\(absent: .*run_settings.*\) · one fresh question, no earlier turns/);
        expect(renderReport([record])).not.toContain('PCA table');
    });

    it('with the PCA table: recorded, in the context, and a column of its own in the report', async () => {
        const record = await runLiveCase({ built, arm: 'C', model: MODEL, callModel: answering(), pcaTable: true });
        expect(record.error).toBeNull();
        expect(record.view.pcaSites?.sites?.length).toBeGreaterThan(0);
        expect(record.context.sections).toContain('pca_displacements');
        expect(renderTranscript(record)).toContain('the PCA table, Live Data off');
        expect(renderReport([record])).toContain('scripted · C + PCA table');
    });

    it('CONTEXT_SECTIONS lists every section a run can have: the Demo with the PCA table has them all', async () => {
        const plain = await openRunFolder(demoEntries());
        const run = await openRunFolder(demoEntries(), { pcaSites: await plain.engine('sites', {}) });
        const record = contextRecord(buildRunContext(run.props));
        expect(record.sections).toEqual([...CONTEXT_SECTIONS]);
        expect(record.absent).toEqual([]);
        expect(record.other_sections).toEqual([]);
    });
});
