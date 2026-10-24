// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// How the verification reads an answer (private to verification/): its
// sentences and clauses, which mentions are claims, and the names in it
// (space-group symbols, element symbols). The claim catalogue (claims.js) and
// the facts catalogue (facts.js) both read through these functions, so a
// sentence is read one way for the evidence gate and for the answer checks.
//
// Deliberately narrow: a mention is a claim only when it is a clear positive
// one — not negated ("has not converged", "achiral", "no preferred
// direction"), not a hypothetical, a condition or a question ("to check
// whether it has converged", "until χ² plateaus", "I tested whether site 23
// has a preferred direction", "a converged run"), and not a JSON key
// ("chiral": false in a structured block). A negation or such a word counts
// only in the stretch of the clause just before the claim, after its last
// comma or "and"/"so": "With no remaining drift in χ², the run has
// converged" is a claim. A false alarm under a correct answer would teach the
// user to ignore the checks. Pure, no React.

import { WINDOW_CUTOFF, WINDOW_RANGES } from './phrasings';

const NEGATION = /\b(not|no|never|neither|nor|none|without|lacks?|lacking|absence|absent|cannot|non|yet|stopped|ceased)\b|n't\b/i;
// What makes a mention hypothetical, conditional or a question about it.
const NOT_ASSERTED = new RegExp('\\b(whether|if|unless|until|once|when(?:ever)?|before|would|could|might|'
    + 'to (?:check|test|see|confirm|verify|determine|assess|judge|tell|decide)|'
    + 'test(?:s|ed|ing)? (?:for|whether|if)|check(?:s|ed|ing)? (?:for|whether|if|that)|'
    + 'look(?:s|ed|ing)? for|make sure|ensure|verify that|confirm that)\\b', 'i');
// An attributive use of a state names a kind, it claims nothing: "a
// converged run" (but "has a preferred direction" is a claim).
const ATTRIBUTIVE = /\b(a|an|any|every|each)\s+$/i;
// The stretch of a clause a claim's negation or condition must sit in: after
// its last comma, "and", "so", "therefore", "thus" or "hence".
const leadOf = (text) => text.split(/,\s*|\s+(?:and|so|therefore|thus|hence)\s+/i).pop();

/** A text's sentences, trimmed, blank ones dropped. */
export const sentencesOf = (text) => String(text || '')
    .replace(/\r/g, '')
    .split(/(?<=[.!?])\s+|\n+/)
    .map((sentence) => sentence.trim())
    .filter(Boolean);

/** A sentence's clauses: a contrast ("…, but …") starts a new claim. */
export const clausesOf = (sentence) => sentence
    .split(/;\s+|:\s+|,\s+(?=(?:but|while|whereas|unlike|although|though|however)\b)/i)
    .map((clause) => clause.trim())
    .filter(Boolean);

/**
 * The asserted matches of `pattern` (global) in a clause: not in a question,
 * and not after a negation or a hypothetical word in the stretch just before
 * it (nor after an article, for a `state` such as "converged"). Group 1 of
 * the pattern is the character before the claim.
 */
export const positiveMatches = (clause, pattern, { state = false } = {}) => {
    if (/\?\s*$/.test(clause)) return [];
    const out = [];
    for (const match of clause.matchAll(pattern)) {
        // A JSON key ("chiral": false) names a field, it claims nothing.
        if (/^"\s*(?::|$)/.test(clause.slice(match.index + match[0].length))) continue;
        const lead = leadOf(clause.slice(0, match.index + (match[1] || '').length));
        if (!NEGATION.test(lead) && !NOT_ASSERTED.test(lead) && !(state && ATTRIBUTIVE.test(lead))) out.push(match);
    }
    return out;
};

/**
 * Whether a text asserts `pattern` somewhere: a positive match in one of its
 * clauses. `within(sentence)`, when given, keeps only the sentences it
 * accepts (the claims about this run, claims.js).
 */
export const asserts = (text, pattern, { within = null } = {}) => sentencesOf(text)
    .some((sentence) => (!within || within(sentence))
        && clausesOf(sentence).some((clause) => positiveMatches(clause, pattern).length > 0));

// --- space-group symbols ----------------------------------------------------------

const SUBSCRIPTS = { '₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4', '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9' };

/** A symbol as written anywhere (P-42₁m, P−4 2_1 m, $P\bar{4}2_1m$, P4̅2₁m) → p-421m. */
export const normalizeSymbol = (text) => String(text || '')
    .replace(/\\bar\{?(\d)\}?/g, '-$1')
    .replace(/(\d)̅/g, '-$1')
    .replace(/[₀-₉]/g, (digit) => SUBSCRIPTS[digit])
    .replace(/[−–—]/g, '-')
    .replace(/\(no\.?\s*\d+\)/gi, '')
    .replace(/[\s_${}\\]/g, '')
    .toLowerCase();

// Text with symbols in one spelling but words kept apart: overbars and
// subscripts written as digits, every dash a hyphen.
const symbolText = (text) => String(text || '')
    .replace(/\\bar\{?(\d)\}?/g, '-$1')
    .replace(/(\d)̅/g, '-$1')
    .replace(/[₀-₉]/g, (digit) => SUBSCRIPTS[digit])
    .replace(/[−–—]/g, '-')
    .toLowerCase();

// A normalized symbol as a pattern that allows spacing, underscores, braces
// and $ inside it (P−4 2_1 m) but needs a word boundary around it.
const symbolPattern = (key) => {
    const filler = '[\\s_{}$\\\\]*';
    const body = [...key].map((char) => char.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')).join(filler);
    return new RegExp(`(?<![a-z0-9])${body}(?![a-z0-9])`);
};

/** The known groups (normalized symbol → group) a piece of text names, each matched as a whole symbol. */
export const groupsIn = (text, groups) => {
    const prepared = symbolText(text);
    return [...groups].filter(([key]) => symbolPattern(key).test(prepared)).map(([, group]) => group);
};

// Something written like a space-group symbol (a lattice letter, then axes
// with at least one digit: P2₁2₁2, I4/mmm, R-3c), known or not.
const SYMBOL_LIKE = /(?<![a-z0-9])[pifcabr](?=[a-z/-]*\d)(?:-?\d(?:_?\d)?|\/?[mnabcde]){1,5}(?![a-z0-9])/;
/** Whether a text names something written like a space-group symbol. */
export const namesSymbol = (text) => SYMBOL_LIKE.test(symbolText(text).replace(/[\s{}$\\]/g, (char) => (/\s/.test(char) ? ' ' : '')));

// --- wavevectors ----------------------------------------------------------------

const UNICODE_FRACTIONS = { '½': 0.5, '⅓': 1 / 3, '⅔': 2 / 3, '¼': 0.25, '¾': 0.75, '⅕': 0.2, '⅙': 1 / 6, '⅛': 0.125 };
const COMPONENT = '[-−+]?(?:\\d+(?:\\.\\d+)?(?:\\/\\d+)?|\\.\\d+|[½⅓⅔¼¾⅕⅙⅛])';
// A wavevector as written after q, k or "wavevector": q = (0.1, 0, 0),
// q=[½ 0 0], wave vector (1/10, 0, 0) — never a bare [u v w] direction.
const Q_VECTOR = new RegExp(`\\b(?:q|k|wave[- ]?vectors?)\\s*(?:=|≈|~|:|of|at)?\\s*[([]\\s*(${COMPONENT})\\s*[,;\\s]\\s*(${COMPONENT})\\s*[,;\\s]\\s*(${COMPONENT})\\s*[)\\]]`, 'gi');

const componentValue = (token) => {
    const sign = /^[-−]/.test(token) ? -1 : 1;
    const body = token.replace(/^[-−+]/, '');
    if (UNICODE_FRACTIONS[body] !== undefined) return sign * UNICODE_FRACTIONS[body];
    const [top, bottom] = body.split('/');
    return sign * (bottom ? Number(top) / Number(bottom) : Number(top));
};

/** The wavevectors a text writes after q, k or "wavevector" ([q1, q2, q3] in cell fractions), in order. */
export const qVectorsIn = (text) => [...String(text || '').matchAll(Q_VECTOR)]
    .map((match) => match.slice(1, 4).map(componentValue))
    .filter((q) => q.every(Number.isFinite));

// --- windows in Å -----------------------------------------------------------------

// A window's upper edge is a bond length: at least this (Å). A displacement
// range ("0.1–0.2 Å") or a join's reach ("within 0.5 Å of the first") is
// no window.
const WINDOW_MIN_TOP_A = 1;

/**
 * The windows in Å a text states ([lo, hi], lo < hi; a cutoff is [0, hi]), in
 * the order they appear: the ranges and cutoffs of phrasings.js WINDOW_RANGES
 * and WINDOW_CUTOFF, each with its top at least WINDOW_MIN_TOP_A.
 */
export const windowsIn = (text) => {
    const source = String(text || '');
    const found = [];
    WINDOW_RANGES.forEach((pattern) => {
        for (const match of source.matchAll(pattern)) found.push({ at: match.index, window: [Number(match[1]), Number(match[2])] });
    });
    for (const match of source.matchAll(WINDOW_CUTOFF)) found.push({ at: match.index, window: [0, Number(match[1])] });
    return found
        .filter(({ window: [lo, hi] }) => Number.isFinite(lo) && Number.isFinite(hi) && lo >= 0 && hi > lo && hi >= WINDOW_MIN_TOP_A)
        .sort((a, b) => a.at - b.at)
        .map(({ window }) => window);
};

// --- element symbols --------------------------------------------------------------

/** An element symbol as written in a formula or alone (Mo, Mo3, Ta₄), not inside a word. */
export const elementPattern = (element) => new RegExp(`(?<![A-Za-z])${element}(?:[0-9₀-₉]+)?(?![a-z])`);
/** Whether a text names an element (elementPattern). */
export const namesElement = (text, element) => elementPattern(element).test(text);
