// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A built case as the Workbench sees an open run: the Dashboard's parsed
// structure, its plot files (the partial g(r)), the symmetry the Dashboard
// hands the Copilot (Detected SG at 0.2 Å + the tolerance ladder), and an
// engine request that runs the browser engines in-process (the same message
// handler as the Copilot's worker). Nothing here is a model.

import { plotDataFromText, structureFromRmc6f } from '../src/browserData';
import { handlePcaMessage } from '../src/workers/pcaKdeWorker';
import { describeSymmetry, toleranceLadder } from '../src/symmetryModel';

// The Detected SG card's default tolerance (Å), as the Dashboard starts.
export const CARD_TOLERANCE_A = 0.2;

// The run's name and its files' names, the same for every case. The model
// reads the run name in its context (`run`), and a case id names the case's
// kind ("dev-zone-boundary"), so no case id goes into a name or a path here;
// the record keeps the id, and the model never sees the record.
export const RUN_NAME = 'rmc_run';
export const RUN_FILES = Object.freeze({
    rmc6f: `${RUN_NAME}.rmc6f`,
    partials: `${RUN_NAME}_PDFpartials.csv`
});

export const loadRun = (built, { ladder = true } = {}) => {
    const structure = structureFromRmc6f({ path: `${RUN_NAME}/${RUN_FILES.rmc6f}`, text: built.rmc6f });
    const plotFiles = [];
    if (built.partials) {
        const file = { name: RUN_FILES.partials, path: `${RUN_NAME}/${RUN_FILES.partials}`, plotKind: 'pdf_partials', text: built.partials };
        plotFiles.push({ ...file, text: undefined, plotData: plotDataFromText(file) });
    }
    const detected = describeSymmetry(structure, CARD_TOLERANCE_A);
    const symmetry = detected
        ? { ...detected, toleranceA: CARD_TOLERANCE_A, ladder: ladder ? toleranceLadder(structure, 1.0) : [] }
        : null;
    const request = (kind, params = {}) => handlePcaMessage({ kind, ...params }, async () => built.rmc6f);
    return {
        runName: RUN_NAME,
        structure,
        plotFiles,
        rValueFile: null,
        symmetry,
        runSettings: null,
        elements: Object.keys(structure.elementCounts || {}).sort(),
        request
    };
};
