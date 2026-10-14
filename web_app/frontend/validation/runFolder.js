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

export const loadRun = (built, { ladder = true } = {}) => {
    const structure = structureFromRmc6f({ path: `${built.id}/${built.rmc6fName}`, text: built.rmc6f });
    const plotFiles = [];
    if (built.partials) {
        const file = { name: built.partialsName, path: `${built.id}/${built.partialsName}`, plotKind: 'pdf_partials', text: built.partials };
        plotFiles.push({ ...file, text: undefined, plotData: plotDataFromText(file) });
    }
    const detected = describeSymmetry(structure, CARD_TOLERANCE_A);
    const symmetry = detected
        ? { ...detected, toleranceA: CARD_TOLERANCE_A, ladder: ladder ? toleranceLadder(structure, 1.0) : [] }
        : null;
    const request = (kind, params = {}) => handlePcaMessage({ kind, ...params }, async () => built.rmc6f);
    return {
        runName: built.id,
        structure,
        plotFiles,
        rValueFile: null,
        symmetry,
        runSettings: null,
        elements: Object.keys(structure.elementCounts || {}).sort(),
        request
    };
};
