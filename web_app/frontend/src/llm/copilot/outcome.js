// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Copilot ends a reply with one verdict line,
//   Outcome: achieved|partial|not achieved — <one sentence>
// which the page shows as a chip. Models decorate it (bold, italics, a list
// bullet or heading marker, an en dash or a colon instead of the em dash),
// and some copy the backticks the prompt shows the line in, so the parse is
// lenient about decoration and strict about the three statuses. The LAST such
// line wins; a reply without one simply has no chip.

export const OUTCOME_STATUSES = ['achieved', 'partial', 'not achieved'];

const OUTCOME_LINE = /^[ \t>#*_`-]*outcome\s*[*_`]*\s*:\s*[*_`]*\s*(achieved|partial|not[ _-]achieved)\b\s*[*_`]*\s*(?:[—–:-]+\s*)?(.*?)\s*[*_`]*\s*$/gim;

/**
 * Split the verdict line out of a reply. Returns
 * `{ status, note, body }` (status one of OUTCOME_STATUSES, note may be '')
 * or null when the reply carries no verdict line.
 */
export const parseOutcome = (text) => {
    if (typeof text !== 'string' || !text) return null;
    let match = null;
    let found;
    OUTCOME_LINE.lastIndex = 0;
    while ((found = OUTCOME_LINE.exec(text)) !== null) {
        match = { index: found.index, length: found[0].length, status: found[1], note: found[2] };
        if (found[0].length === 0) OUTCOME_LINE.lastIndex += 1;
    }
    if (!match) return null;
    const status = match.status.toLowerCase().replace(/[ _-]+/, ' ');
    const body = `${text.slice(0, match.index)}${text.slice(match.index + match.length)}`
        .replace(/\n{3,}/g, '\n\n')
        .trim();
    return { status, note: match.note.replace(/[*_`]+$/, '').trim(), body };
};
