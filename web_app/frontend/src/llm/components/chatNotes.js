// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The chat's notes on what the model was not sent whole (THY-76: never
// truncate silently).

/**
 * What an answer's memory left out that the last one sent (useCopilot.js):
 * `dropped` oldest questions no longer sent (the model keeps the last
 * `kept`), and whether the last answer goes shortened (compactAnswer).
 */
export const memoryNote = ({ dropped = 0, kept = 0, shortened = false } = {}) => {
    const parts = [
        dropped ? `the earliest ${dropped === 1 ? 'question is' : `${dropped} questions are`} no longer sent to the model (it keeps the last ${kept})` : '',
        shortened ? 'the last answer is sent shortened to its sentences with numbers' : ''
    ].filter(Boolean);
    return parts.length ? `Memory: ${parts.join('; ')}.` : '';
};

/** Earlier results of this answer sent shortened (agentLoop.js turnCharBudget). */
export const compactedNote = (count) => (count > 0
    ? `${count} earlier result${count === 1 ? ' was' : 's were'} sent shortened to ${count === 1 ? 'its' : 'their'} checks and key numbers to fit`
    : '');
