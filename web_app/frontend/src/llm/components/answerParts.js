// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A long answer as the dock shows it: its lead — the conclusion, which the
// prompt asks for first (prompts/copilot.js ANSWER_RULES) — then the rest
// behind "Show details", so the dock reads like the instrument's result
// first and the numbers on request (the user's report of 2026-11-01:
// "overwhelmingly long … dumping numbers"). Every analysis's own numbers
// stay on its step card. An answer as the prompt asks for it — a conclusion,
// a few numbers, a caveat — shows whole; one that runs long, or whose rest
// is a dump (a table, headings, code, more than a handful of numbers),
// folds. Nothing is ever dropped, only folded.
//
//   answerParts(markdown) → { lead, rest }   rest '' when the answer shows whole
//
// Blocks are the Markdown's own: paragraphs, lists, tables, quotes and fenced
// code, split at blank lines outside a fence. The lead is the first block,
// with the next one when the first is a heading or a short line (a title, a
// "Short answer:"). The rest folds when it holds a table, a heading or code;
// otherwise when it is at least FOLD_MIN_REST_CHARS long and the answer is
// FOLD_MIN_CHARS or longer or the rest holds more than FOLD_MAX_NUMBERS
// numbers. Pure, no React.

/** An answer this long folds its rest. */
export const FOLD_MIN_CHARS = 700;
/** A rest shorter than this shows with the lead, unless it is a table, a heading or code. */
export const FOLD_MIN_REST_CHARS = 160;
/** A rest with more numbers than this folds, however long the answer. */
export const FOLD_MAX_NUMBERS = 8;
// What a dump looks like: a table row, a heading, a fence.
const STRUCTURED = /^\s*(?:\|.*\||#{1,6}\s|```|~~~)/m;
// A number as written (not a digit inside a word: "Se6", "r2").
const NUMBER = /(?<![\w.])\d+(?:\.\d+)?/g;
// A first block this short leads with the next one too.
const SHORT_LEAD_CHARS = 90;

const FENCE = /^\s*(```|~~~)/;

/** The Markdown's blocks, in order: split at blank lines outside fenced code. */
export const markdownBlocks = (markdown) => {
    const blocks = [];
    let current = [];
    let fence = null;
    String(markdown ?? '').split('\n').forEach((line) => {
        const opens = line.match(FENCE);
        if (fence) {
            current.push(line);
            if (opens && opens[1] === fence) fence = null;
            return;
        }
        if (opens) {
            fence = opens[1];
            current.push(line);
            return;
        }
        if (!line.trim()) {
            if (current.length) blocks.push(current.join('\n'));
            current = [];
            return;
        }
        current.push(line);
    });
    if (current.length) blocks.push(current.join('\n'));
    return blocks;
};

const isHeading = (block) => block.split('\n').every((line) => /^\s{0,3}#{1,6}\s/.test(line) || /^\s*\*\*[^*]+\*\*:?\s*$/.test(line));

/** The answer's lead and the rest the dock folds (see the header). */
export const answerParts = (markdown) => {
    const text = String(markdown ?? '').trim();
    const whole = { lead: text, rest: '' };
    const blocks = markdownBlocks(text);
    if (blocks.length < 2) return whole;
    const first = blocks[0];
    const leadCount = (isHeading(first) || first.length < SHORT_LEAD_CHARS) && blocks.length > 2 ? 2 : 1;
    const rest = blocks.slice(leadCount).join('\n\n');
    const folded = { lead: blocks.slice(0, leadCount).join('\n\n'), rest };
    if (STRUCTURED.test(rest)) return folded;
    if (rest.length < FOLD_MIN_REST_CHARS) return whole;
    return text.length >= FOLD_MIN_CHARS || (rest.match(NUMBER) || []).length > FOLD_MAX_NUMBERS ? folded : whole;
};
