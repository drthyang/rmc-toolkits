// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Numbers by reference (THY-18). Every tool result the model reads carries a
// `ref` — r1, r2, … in the order the turn ran them — and the model quotes a
// value as {{r2.peak.significance_sigma}}: the result's ref, then the path of
// the value in the data it was shown (an array item by its index, {{r3.sites.0.shift_A}}
// or {{r3.sites[0].shift_A}}). The run context is {{ctx.<path>}}. The app
// writes the engine's value in place of each reference where the words are
// read — the answer, the record and the score — so no number is mistyped,
// rounded otherwise than the engine rounds it, or invented. A reference that
// names nothing the turn has stays as written and is a finding on the answer,
// which the computed verdict counts (verification/: the turn's store, filled
// before the answer is judged; turnState.js fills the earlier text segments
// and a stopped turn's from the steps it shows).
//
// Values print as the tools gave them: a number as written in the data, a
// list of numbers as [a, b, c] (valid in prose and in JSON), a string as is —
// so inside JSON the writer adds its quotes ("element": "{{r1.sites.0.element}}"),
// and a list holding strings reads [a, b] in prose only. Pure, no React.

// A reference's name: a root (r<n> or ctx), then .key or [n] steps; keys as
// the tools name them (letters, digits, _ - < > +, as in "<100>").
const NAME = '(?:r\\d+|ctx)(?:\\.[\\w<>+-]+|\\[\\d+\\])*';
const NAME_ONLY = new RegExp(`^${NAME}$`);
/** {{name}} in a text. */
export const REFERENCE = new RegExp(`\\{\\{\\s*(${NAME})\\s*\\}\\}`, 'g');
// Anything in double braces, for the references that do not parse.
const BRACED = /\{\{\s*([^{}]{1,120}?)\s*\}\}/g;

const steps = (path) => {
    const out = [];
    for (const match of path.matchAll(/\.([\w<>+-]+)|\[(\d+)\]|^(r\d+|ctx)/g)) out.push(match[1] ?? match[2] ?? match[3]);
    return out;
};

/** A value as it is written in place of its reference, or null when it cannot be. */
export const formatReferenceValue = (value) => {
    if (typeof value === 'number') return Number.isFinite(value) ? String(value) : null;
    if (typeof value === 'string') return value;
    if (typeof value === 'boolean') return String(value);
    if (Array.isArray(value) && value.length && value.every((item) => ['number', 'string', 'boolean'].includes(typeof item))) {
        return `[${value.map((item) => (typeof item === 'string' ? item : String(item))).join(', ')}]`;
    }
    return null;
};

/**
 * The turn's referable values: `put(ref, data)` keeps a tool result the model
 * was shown under its ref; `resolve(path)` gives the value text or null.
 * `context` is what {{ctx.…}} reads: the run context as the model read it.
 */
export const createReferenceStore = ({ context = null } = {}) => {
    const results = new Map();
    return {
        // Under a ref the caller named (the verification names it before the
        // loop fits the data into the message that carries the ref).
        put: (ref, data) => { results.set(ref, data); },
        resolve: (path) => {
            const [root, ...rest] = steps(path);
            let value = root === 'ctx' ? context : results.get(root);
            for (const key of rest) {
                if (value == null || typeof value !== 'object') return null;
                value = value[key];
            }
            return value === undefined ? null : formatReferenceValue(value);
        },
        size: () => results.size
    };
};

/**
 * The text with every reference the store can fill written in:
 * `{ text, written, filled, unknown }` — `written` counts the references,
 * `unknown` lists those (as written) that name no value, which stay in the text.
 */
export const fillReferences = (text, store) => {
    const source = String(text ?? '');
    if (!source.includes('{{')) return { text: source, written: 0, filled: 0, unknown: [] };
    const unknown = [];
    let written = 0;
    let filled = 0;
    const out = source.replace(BRACED, (whole, inner) => {
        const name = inner.trim();
        // Not a reference at all (a LaTeX group, a template): left alone.
        if (!/^(?:r\d+|ctx)\b/.test(name)) return whole;
        written += 1;
        const value = NAME_ONLY.test(name) && store ? store.resolve(name) : null;
        if (value == null) {
            if (!unknown.includes(name)) unknown.push(name);
            return whole;
        }
        filled += 1;
        return value;
    });
    return { text: out, written, filled, unknown };
};
