// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The reference repair: an answer whose numbers by reference name no value
// (references.js) goes back to the model once, with what each such path
// reaches, before anyone reads it — the user would otherwise see
// {{r5.…}} where a number belongs, and a validation verdict holding one
// does not parse. The app never fills a near miss (a key of a similar
// name, the mean of an object named whole): a wrong number is worse than a
// question asked again, so the model names the value itself.
//
//   referenceRepair(unknown, store) → { id, what, missing, nudge, notice } | null
//
// `unknown` are the references as written that named no value
// (fillReferences), `store` the turn's reference store (its explain). Null
// when there is nothing to repair. The send-back has the evidence gate's
// shape (claims.js sendBackOf), its id REFERENCE_REPAIR, and the turn
// policy asks for it as a 'repair' (turnPolicy.js). The nudge is the same
// under either method: it names no analysis. Pure, no React.

/** The id of a reference repair among a turn's send-backs (`gates`). */
export const REFERENCE_REPAIR = 'reference';

/** Repairs per turn: the answer goes back once for its references. */
export const DEFAULT_MAX_REPAIRS = 1;

// At most this many keys listed for one path; at most this many references
// explained in one nudge (the rest are named).
const KEYS_SHOWN = 12;
const EXPLAINED = 8;

const keyList = (keys) => (keys.length > KEYS_SHOWN
    ? `${keys.slice(0, KEYS_SHOWN).join(', ')}, … (${keys.length} in all)`
    : keys.join(', '));

// What one reference reaches, as a clause after its name.
const reachedText = (explained) => {
    if (!explained) return 'names no value';
    const { reason, at, keys, missing } = explained;
    if (reason === 'no_result') return `names no result of this question (its results are ${keys.length ? keyList(keys) : 'none'})`;
    if (reason === 'not_a_reference') return 'is not a reference: write the ref, then .key or .index steps, e.g. {{r2.peak.significance_sigma}}';
    if (reason === 'not_a_value') {
        return keys.length
            ? `names ${at}, which holds ${keyList(keys)}: name one value in it`
            : `names ${at}, which holds no value`;
    }
    return keys.length
        ? `names nothing at "${missing}": ${at} holds ${keyList(keys)}`
        : `names nothing at "${missing}": ${at} is a single value`;
};

/**
 * The repair send-back for an answer's unfilled references, or null when
 * there are none.
 */
export const referenceRepair = (unknown, store) => {
    if (!Array.isArray(unknown) || !unknown.length) return null;
    const lines = unknown.slice(0, EXPLAINED).map((name) => `- {{${name}}} ${reachedText(store?.explain?.(name) ?? null)}.`);
    if (unknown.length > EXPLAINED) lines.push(`- and ${unknown.slice(EXPLAINED).map((name) => `{{${name}}}`).join(', ')}.`);
    const count = unknown.length;
    const what = `${count} number${count === 1 ? '' : 's'} by reference that name${count === 1 ? 's' : ''} no value`;
    return {
        id: REFERENCE_REPAIR,
        what,
        missing: [...unknown],
        nudge: [
            `(From the app, not the user.) Your answer quotes ${what}, so the user would see ${count === 1 ? 'it' : 'them'} unfilled:`,
            ...lines,
            'Write the whole answer again, with each of these references the exact path of one value in that result\'s data '
                + '(as {{r2.peak.significance_sigma}} or {{r2.data.peak.significance_sigma}}), or without the number where no result '
                + 'holds it. Do not guess a path, and do not run an analysis for this.'
        ].join('\n'),
        notice: `The answer quoted ${what}; the app asked for ${count === 1 ? 'it' : 'them'} again.`
    };
};
