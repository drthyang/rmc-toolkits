// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { structureFromRmc6f } from '../browserData';
import { requestNumber, requestObject } from './requestGuards.js';

// As app.MAX_STRUCTURE_POINTS and /api/structure's maxPoints rule: an integer,
// clamped to [100, 1e6]; missing means the maximum.
const MAX_STRUCTURE_POINTS = 1_000_000;

export const handleStructureMessage = async (data) => {
    const { file, maxPoints: rawMaxPoints } = requestObject(data);
    const maxPoints = requestNumber(rawMaxPoints, 'maxPoints', {
        fallback: MAX_STRUCTURE_POINTS, integer: true, clamp: [100, MAX_STRUCTURE_POINTS]
    });
    if (!file?.sourceFile) {
        throw new Error('No browser structure file available');
    }
    const text = await file.sourceFile.text();
    return structureFromRmc6f({ ...file, text }, maxPoints);
};

// Guarded so the module can be imported by tests outside a worker context.
if (typeof self !== 'undefined' && typeof self.postMessage === 'function') {
    self.onmessage = async (event) => {
        // The id is read defensively and the request validated inside the try,
        // so a null or malformed message still gets an answer.
        const data = event?.data;
        const id = data !== null && typeof data === 'object' ? data.id : undefined;
        try {
            const result = await handleStructureMessage(data);
            self.postMessage({ id, result });
        } catch (error) {
            self.postMessage({
                id,
                error: error.message || 'Browser structure parser failed'
            });
        }
    };
}
