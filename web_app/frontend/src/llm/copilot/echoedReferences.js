// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A number the model typed AND quoted by reference, as the dock shows it.
// Some local models (gemma4:26b) write "2.64 Å {{r1.first_shell.mean_d_A}}":
// the value, then its reference, so the filled answer read "2.64 Å 2.64".
// Where the reference's value is the number just before it (to that number's
// decimals), the dock shows the engine's value once, with the unit or short
// label that sat between them ("3.25 Se per Ta {{r2.…}}" → "3.25 Se per Ta").
//
// Display only: the dock's transcript (turnState.js) and its streaming text
// (components/referenceDisplay.js) use it. The verification, the answer the
// model reads back in its history and the validation never do — what the
// model wrote and what is judged are unchanged. Pure, no React.

import { isReferenceName } from './references';

// A number, an optional short unit or label, then a reference in double
// braces. The label is at most three words of letters (any script: Å, σ,
// element symbols), superscripts ² ³, °, %, "/" and "per" — "Å", "deg²",
// "Se per Ta" — on one line: no digits and no punctuation that ends a clause
// (, ; : . ! ? brackets, quotes), so "6 neighbours, of which {{…}}" is left
// alone. MAX_LABEL bounds its length.
const WORD = '[\\p{L}\\p{M}²³°%/]+';
const ECHO = new RegExp(
    `(-?\\d+(?:\\.\\d+)?)((?:[ \\t]?${WORD}(?:[ \\t]+${WORD}){0,2})?)([ \\t]*)\\{\\{\\s*([^{}]{1,120}?)\\s*\\}\\}`,
    'gu'
);
const MAX_LABEL = 24;

const decimals = (literal) => (literal.includes('.') ? literal.length - literal.indexOf('.') - 1 : 0);

// Whether `value` (the reference's text) is the typed number, to its decimals.
const sameNumber = (literal, value) => {
    const typed = Number(literal);
    const engine = Number(value);
    if (!Number.isFinite(typed) || !Number.isFinite(engine)) return false;
    return engine.toFixed(decimals(literal)) === typed.toFixed(decimals(literal));
};

const echoes = (text, store) => {
    if (typeof text !== 'string' || !text.includes('{{') || !store) return [];
    const found = [];
    for (const match of text.matchAll(ECHO)) {
        const [whole, literal, unit, gap, inner] = match;
        if (unit.length > MAX_LABEL) continue;
        const name = inner.trim();
        if (!isReferenceName(name)) continue;
        const value = store.resolve(name);
        if (value == null || !sameNumber(literal, String(value))) continue;
        found.push({ whole, name, value: String(value), literal, unit, gap });
    }
    return found;
};

/**
 * The text with each echoed number dropped before its reference:
 * "2.64 Å {{r1.x}}" → "{{r1.x}} Å" (filled as usual afterwards). `store` the
 * turn's reference store (its `resolve`). Any other text is returned as given.
 */
export const collapseEchoes = (text, store) => echoes(text, store).reduce(
    (out, echo) => out.replace(echo.whole, `{{${echo.name}}}${echo.unit}`),
    text
);

/**
 * `filled` (the answer as the verification filled it) with each number that
 * `raw` (the same answer as the model wrote it) typed before its own
 * reference shown once: "2.64 Å 2.64" → "2.64 Å".
 */
export const collapseFilledEchoes = (filled, raw, store) => (typeof filled === 'string'
    ? echoes(raw, store).reduce(
        (out, echo) => out.replace(`${echo.literal}${echo.unit}${echo.gap}${echo.value}`, `${echo.value}${echo.unit}`),
        filled
    )
    : filled);
