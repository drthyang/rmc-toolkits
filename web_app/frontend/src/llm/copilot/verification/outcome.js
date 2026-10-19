// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Copilot ends a reply with one Outcome line,
//   Outcome: <one sentence>
// which explains the verdict the app computes from its checks (verdict.js).
// Older prompts asked for a status first (Outcome: achieved|partial|not
// achieved — <one sentence>), and models still write one: it is read, as
// `status` (the model's claim, never the verdict), and the sentence as the
// note. Models decorate the line (bold, italics, a list bullet or heading
// marker, an en dash or a colon instead of the em dash), and some copy the
// backticks the prompt shows it in, so the parse is lenient about decoration.
// The LAST such line wins; a reply without one has no line to split out.

export const OUTCOME_STATUSES = ['achieved', 'partial', 'not achieved'];

const OUTCOME_LINE = /^[ \t>#*_`-]*outcome\s*[*_`]*\s*:\s*[*_`]*\s*(?:(achieved|partial|not[ _-]achieved)\b\s*[*_`]*\s*(?:[—–:-]+\s*)?)?(.*?)\s*[*_`]*\s*$/gim;

/**
 * Split the Outcome line out of a reply. Returns `{ status, note, body }`
 * (status one of OUTCOME_STATUSES, or null for a line with no status word;
 * note may be '') or null when the reply carries no Outcome line.
 */
export const parseOutcome = (text) => {
    if (typeof text !== 'string' || !text) return null;
    let match = null;
    let found;
    OUTCOME_LINE.lastIndex = 0;
    while ((found = OUTCOME_LINE.exec(text)) !== null) {
        match = { index: found.index, length: found[0].length, status: found[1] ?? null, note: found[2] };
        if (found[0].length === 0) OUTCOME_LINE.lastIndex += 1;
    }
    if (!match) return null;
    const status = match.status ? match.status.toLowerCase().replace(/[ _-]+/, ' ') : null;
    const body = `${text.slice(0, match.index)}${text.slice(match.index + match.length)}`
        .replace(/\n{3,}/g, '\n\n')
        .trim();
    return { status, note: match.note.replace(/[*_`]+$/, '').replace(/^[—–-]+\s*/, '').trim(), body };
};
