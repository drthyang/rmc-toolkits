// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The run as the Copilot sees it: one module turns an open run and the view
// state into everything the AI Copilot is handed about it — CopilotChat's
// eight run props, the engine tools, the run's elements, the engine and the
// normalised view. This file is the module's pure core: no React, no I/O, no
// worker. Its two adapters feed it:
//   - useCopilotRun.js, the browser adapter: the run the Dashboard publishes,
//     the parsed <stem>.dat, and the engine in the Copilot's own worker or in
//     Flask;
//   - openRunFolder.js, the in-process adapter (Node and tests): a folder's
//     files through the Demo / picked-folder path, the workers' own handlers
//     in-thread. The validation harness (validation/runFolder.js) calls it.
// The seam is the LoadedRun below; whatever builds one, the Copilot is
// handed the same thing for the same files and view.
//
//   LoadedRun  { name, files, structure, structurePath, settings, engine }
//     name:          the run folder's name (the picked or Demo root, the
//                    Flask directory, the harness's root folder)
//     files:         as the Dashboard holds them (plotKind from
//                    detectPlotKind; plotData or parseError once parsed)
//     structure:     the parsed .rmc6f (structureFromRmc6f, maxPoints 100),
//                    or Flask's /api/structure (the same fields, the
//                    average-structure basis included), or null
//     structurePath: the .rmc6f the Model card describes (picks the χ² logs)
//     settings:      the parsed <stem>.dat, or null (always null in Flask mode)
//     engine:        (kind, params, { signal }) => Promise<payload>, or null
//   View       { toleranceA, hiddenPlots, pcaSites, liveData, ladder }
//   CopilotRun { props, tools, elements, engine, view }

import { chartedRValueFile, dashboardPlotFiles } from './browserData';
import { buildCopilotTools } from './copilotTools';
import { describeSymmetry, toleranceLadder } from './symmetryModel';

/**
 * The view the Copilot is handed when nothing on screen says otherwise: the
 * Detected SG card's starting tolerance (Å), no plot hidden, no PCA table (the
 * PCA Ellipsoid page not visited), Live Data off, the tolerance ladder on.
 * Frozen: the one place the card's default lives for the app and the harness.
 */
export const DEFAULT_VIEW = Object.freeze({
    toleranceA: 0.2,
    hiddenPlots: Object.freeze([]),
    pcaSites: null,
    liveData: false,
    ladder: true
});

/** CopilotChat's run props, the list in AGENTS.md ("src/llm/ import boundary"). */
export const COPILOT_RUN_PROPS = Object.freeze([
    'runName', 'plotFiles', 'rValueFile', 'structure', 'symmetry', 'runSettings', 'pcaSites', 'liveData'
]);

// The ladder spans tight → this loose (Å), as the Detected SG card's.
const LADDER_MAX_A = 1.0;

// A hidden-plot list already in normal form (frozen, strings, strictly sorted)
// is kept as it is, so normalising a normalised view keeps its identity.
const isNormalHiddenList = (list) => Array.isArray(list) && Object.isFrozen(list)
    && list.every((path, index) => typeof path === 'string' && (index === 0 || list[index - 1] < path));

/**
 * A view with every field, normalised and JSON-safe (frozen): missing fields
 * take DEFAULT_VIEW's, hiddenPlots is a sorted list of distinct paths, the
 * switches are booleans. Throws RangeError on a tolerance that is not a
 * positive number.
 */
export const normalizeView = (view = DEFAULT_VIEW) => {
    const merged = { ...DEFAULT_VIEW, ...(view || {}) };
    const toleranceA = merged.toleranceA;
    if (typeof toleranceA !== 'number' || !Number.isFinite(toleranceA) || toleranceA <= 0) {
        throw new RangeError(`the Detected SG tolerance must be a positive number of Å, got ${toleranceA}`);
    }
    const hidden = merged.hiddenPlots || [];
    return Object.freeze({
        toleranceA,
        hiddenPlots: isNormalHiddenList(hidden) ? hidden : Object.freeze([...new Set([...hidden].map(String))].sort()),
        pcaSites: merged.pcaSites ?? null,
        liveData: Boolean(merged.liveData),
        ladder: merged.ladder !== false
    });
};

// --- the fields, each a pure function of exactly its own inputs ---------------------
// copilotRun composes them; useCopilotRun memoizes each on the same inputs, so
// a field keeps its identity until its own inputs change.

/** Every Dashboard plot file in Dashboard order, hidden ones included. */
export const runPlotFiles = (files) => dashboardPlotFiles(files || []);

/** The χ² history on screen: the visible logs, for the run the structure file names. */
export const runRValueFile = (plotFiles, hiddenPlots, structurePath) => (
    chartedRValueFile(plotFiles, hiddenPlots, structurePath ?? null)
);

/** The Detected SG tolerance ladder, or [] with the ladder off or no basis to read. */
export const runLadder = (structure, ladder) => (
    ladder && structure?.basis?.length ? toleranceLadder(structure, LADDER_MAX_A) : []
);

/**
 * The Detected SG card at the view's tolerance, with that tolerance and the
 * ladder; null with no structure, no basis (a coords-only .rmc6f, a server
 * before /api/structure sent one, before load) or no result.
 */
export const runSymmetry = (structure, toleranceA, ladder) => {
    if (!structure?.basis?.length) return null;
    const detected = describeSymmetry(structure, toleranceA);
    return detected ? { ...detected, toleranceA, ladder } : null;
};

/** The run's elements: the sorted keys of the structure's element counts. */
export const runElements = (structure) => Object.keys(structure?.elementCounts || {}).sort();

/** The engine tools for the run (copilotTools.js): none without an engine, no symmetry tool without a basis. */
export const runEngineTools = (engine, structure, elements) => buildCopilotTools({ request: engine || null, structure: structure || null, elements });

/** The CopilotRun from its finished fields: `props` frozen, passed whole to CopilotChat. */
export const assembleCopilotRun = ({
    name, plotFiles, rValueFile, structure, symmetry, settings, tools, elements, engine, view
}) => ({
    props: Object.freeze({
        runName: name,
        plotFiles,
        rValueFile: rValueFile ?? null,
        structure: structure ?? null,
        symmetry: symmetry ?? null,
        runSettings: settings ?? null,
        pcaSites: view.pcaSites,
        liveData: view.liveData
    }),
    tools,
    elements,
    engine: engine || null,
    view
});

/**
 * The run as the Copilot sees it, from a LoadedRun and a view: pure and
 * synchronous. Never throws on a missing piece (no structure, no files, no
 * engine); throws RangeError on a non-positive tolerance.
 */
export const copilotRun = (loaded, view = DEFAULT_VIEW) => {
    const normalized = normalizeView(view);
    const structure = loaded?.structure ?? null;
    const engine = loaded?.engine || null;
    const plotFiles = runPlotFiles(loaded?.files);
    const elements = runElements(structure);
    return assembleCopilotRun({
        name: loaded?.name ?? '',
        plotFiles,
        rValueFile: runRValueFile(plotFiles, normalized.hiddenPlots, loaded?.structurePath),
        structure,
        symmetry: runSymmetry(structure, normalized.toleranceA, runLadder(structure, normalized.ladder)),
        settings: loaded?.settings ?? null,
        tools: runEngineTools(engine, structure, elements),
        elements,
        engine,
        view: normalized
    });
};
