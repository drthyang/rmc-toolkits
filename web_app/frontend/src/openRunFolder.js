// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The run module's in-process adapter (copilotRun.js is its pure core), for
// Node and tests: a run folder's files read as the app reads a picked folder
// or the Demo (browserData's makeRunFromEntries), every plot file and the
// structure parsed as the Dashboard parses them, and an engine that runs the
// Copilot's worker handler in-thread. The validation harness opens every case
// through it (validation/runFolder.js), so a case's run reaches the model
// exactly as the same files would reach it in the app.
//
// No app module may import this file (eslint.config.js, no-restricted-imports
// for src/**; tests and validation/ are exempt): the workers' modules install
// `self.onmessage` whenever `self.postMessage` is a function, and on a window
// it is. They are imported here dynamically, on the first open.

import {
    DASHBOARD_MAX_POINTS, isDashboardPlotFile, makeRunFromEntries, parseDashboardPlotFile, parseRunSettings
} from './browserData';
import { copilotRun, DEFAULT_VIEW, normalizeView } from './copilotRun';

/**
 * A run folder as a LoadedRun ({ name, files, structure, structurePath,
 * settings, engine }), resolved once every plot file and the structure are
 * parsed (maxPoints 100, as on the Dashboard and in Flask). `entries` are
 * `[{ path: '<folder>/<file>', file: File }]`. Rejects with "No supported
 * RMCProfile files…" or with the .rmc6f's parse error; a plot file that fails
 * keeps its parseError, as on the Dashboard, and an unreadable .dat gives null.
 */
export const readRunFolder = async (entries) => {
    const [{ inProcessEngine }, { handleStructureMessage }] = await Promise.all([
        import('./workers/pcaKdeWorker'),
        import('./workers/localStructureWorker')
    ]);
    const run = await makeRunFromEntries(entries);
    const files = await Promise.all(run.files.map((file) => (isDashboardPlotFile(file) ? parseDashboardPlotFile(file) : file)));
    const structure = run.structureFile
        ? await handleStructureMessage({ file: run.structureFile, maxPoints: DASHBOARD_MAX_POINTS })
        : null;
    // The PCA worker's own handler in-thread, behind requestAnalysis's
    // interface (pcaKdeWorker's inProcessEngine), the text read once.
    const sourceFile = run.structureFile?.sourceFile;
    let settings = null;
    if (run.settingsFile?.sourceFile) {
        try {
            settings = parseRunSettings(await run.settingsFile.sourceFile.text());
        } catch {
            settings = null;
        }
    }
    return {
        name: run.name,
        files,
        structure,
        structurePath: run.structureFile?.path ?? null,
        settings,
        engine: run.structureFile ? inProcessEngine(() => sourceFile.text()) : null
    };
};

/** The run as the Copilot sees it, from a folder's files and a view (copilotRun over readRunFolder). */
export const openRunFolder = async (entries, view = DEFAULT_VIEW) => {
    // A bad view fails before any file is read.
    const normalized = normalizeView(view);
    return copilotRun(await readRunFolder(entries), normalized);
};
