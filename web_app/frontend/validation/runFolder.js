// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A built case as the Workbench sees an open run: its files, written as one
// run folder, opened through the run module's in-process adapter
// (src/openRunFolder.js) exactly as the app opens a picked folder — the
// Dashboard's parses, its plot files and their kinds (detectPlotKind), the
// symmetry the Copilot is handed at the view's tolerance, and an engine that
// runs the Copilot worker's handler in-thread. Nothing here is a model.

import { DEFAULT_VIEW, copilotRun, normalizeView } from '../src/copilotRun';
import { openRunFolder, readRunFolder } from '../src/openRunFolder';

// The run's name and its files' names, the same for every case. The model
// reads the run name in its context (`run`), and a case id names the case's
// kind ("dev-zone-boundary"), so no case id goes into a name or a path here;
// the record keeps the id, and the model never sees the record.
export const RUN_NAME = 'rmc_run';
export const RUN_FILES = Object.freeze({
    rmc6f: `${RUN_NAME}.rmc6f`,
    partials: `${RUN_NAME}_PDFpartials.csv`,
    // A case with fit data (built.fitFiles): RMCProfile's X-ray file family,
    // named as RMCProfile names them, so the run module detects their kinds
    // (detectPlotKind: xray_sq, xpdf, r_value) as for a user's run.
    fq: `${RUN_NAME}_FQ1.csv`,
    gr: `${RUN_NAME}_FT_XFQ1.csv`,
    log: `${RUN_NAME}-00.log`
});

/**
 * The view every live run is given unless said otherwise: the app's default
 * (the Detected SG card at 0.2 Å with its ladder, nothing hidden, no PCA
 * table, Live Data off). Each record keeps the view it ran with.
 */
export const HARNESS_VIEW = Object.freeze({ ...DEFAULT_VIEW });

// A case file as the browser hands a picked folder's file to the app; the
// time is fixed so nothing about the moment of the run reaches it.
const entry = (name, text) => ({ path: `${RUN_NAME}/${name}`, file: new File([text], name, { lastModified: 0 }) });

/**
 * The case's run-folder files as `[{ name, text }]`: the .rmc6f, the partial
 * g(r) when the case has one, and its fit files (F(Q), its Fourier transform
 * G(r), the χ² log) when it has fit data. `run.js write` should write these.
 */
export const caseFiles = (built) => [
    { name: RUN_FILES.rmc6f, text: built.rmc6f },
    ...(built.partials ? [{ name: RUN_FILES.partials, text: built.partials }] : []),
    ...(built.fitFiles ? ['fq', 'gr', 'log'].map((role) => ({ name: RUN_FILES[role], text: built.fitFiles[role] })) : [])
];

/** The case's run folder as the browser hands a picked folder to the app (caseFiles). */
export const caseEntries = (built) => caseFiles(built).map(({ name, text }) => entry(name, text));

/**
 * The case as the Copilot sees it (a CopilotRun: { props, tools, elements,
 * engine, view }). `view` is merged over HARNESS_VIEW. `pcaTable: true` gives
 * the run the PCA table (view.pcaSites) from its own engine, the rows
 * list_sites reads, as after a visit to the PCA Ellipsoid page: a dev check of
 * that context section's effect, off by default.
 */
export const loadRun = async (built, view = HARNESS_VIEW, { pcaTable = false } = {}) => {
    const normalized = normalizeView({ ...HARNESS_VIEW, ...view });
    if (!pcaTable) return openRunFolder(caseEntries(built), normalized);
    // The table comes from the run's own engine, so the folder is read first
    // and the view built after.
    const loaded = await readRunFolder(caseEntries(built));
    if (!loaded.engine) throw new Error('the PCA table needs the run\'s engine (a .rmc6f)');
    return copilotRun(loaded, { ...normalized, pcaSites: await loaded.engine('sites', {}) });
};
