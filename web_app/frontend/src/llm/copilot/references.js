// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Numbers by reference (THY-18). Every tool result the model reads carries a
// `ref` — r1, r2, … in the order the turn ran them — and the model quotes a
// value as {{r2.peak.significance_sigma}}: the result's ref, then the path of
// the value in the data it was shown (an array item by its index, {{r3.sites.0.shift_A}}
// or {{r3.sites[0].shift_A}}). The same path written from the root of the tool
// message the model read — {"ref": "r2", "checks": […], "data": {…}} —,
// {{r2.data.peak.significance_sigma}}, names the same value: most references
// gemma4:26b wrote in the validation's first live runs had that form, and so
// did a model's in the dock (2026-11-01), all of them left unfilled until
// then. That reading is exact, not a guess: every tool message holds its
// result under `data` (tools.js toolMessage), and no tool's data has a
// top-level key of that name (toolDataShapes.test.js); were one to, the
// data's own path would be read. The run context is {{ctx.<path>}}. The app
// writes the engine's value in place of each reference where the words are
// read — the answer, the record and the score — so no number is mistyped,
// rounded otherwise than the engine rounds it, or invented. A reference that
// names nothing the turn has is never filled with a near miss (a key of a
// similar name, the mean of an object named whole); `explain` says what its
// path reaches, and the verification sends the answer back once with that
// (verification/repair.js). Still unfilled, it stays as written in the text
// and is a finding on the answer, which the computed verdict counts (verification/: the turn's
// store, filled before the answer is judged; turnState.js fills the earlier
// text segments and a stopped turn's from the steps it shows).
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
// What a reference starts with: anything else in double braces (a LaTeX
// group, a template) is not one.
const REFERENCE_START = /^(?:r\d+|ctx)\b/;

/** The key a tool message holds its result's data under (tools.js toolMessage). */
export const MESSAGE_DATA_KEY = 'data';

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

const isRecord = (value) => value !== null && typeof value === 'object';

// What a value holds, as a reference steps into it: an object's keys, a
// list's indices as one "0…n-1", nothing for a scalar.
const holds = (value) => {
    if (Array.isArray(value)) return value.length ? [`0…${value.length - 1}`] : [];
    return isRecord(value) ? Object.keys(value) : [];
};

/**
 * The turn's referable values: `put(ref, data)` keeps a tool result the model
 * was shown under its ref; `resolve(path)` gives the value text or null;
 * `explain(path)` says, for a path that gives no value, what it reaches —
 * `{ reason, at, keys }`: 'no_result' (no result has that ref; `keys` the refs
 * there are), 'no_key' (the step after `at`, `missing`, names nothing; `keys`
 * what the value at `at` holds), 'not_a_value' (the path ends at an object or a list
 * that is no list of scalars; `keys` what it holds) or 'not_a_reference' (no
 * reference's grammar) — and null for a path that gives one. `at` is the
 * path as far as it led, as the model wrote it. `context` is what {{ctx.…}}
 * reads: the run context as the model read it.
 */
export const createReferenceStore = ({ context = null } = {}) => {
    const results = new Map();
    // Where a path leads: { found, value, at, keys, missing } (or
    // { noResult: true } for a ref with no result).
    const walk = (path) => {
        const [root, ...rest] = steps(path);
        if (root !== 'ctx' && !results.has(root)) return { noResult: true };
        const base = root === 'ctx' ? context : results.get(root);
        // The message-root form ({{r2.data.x}}), unless the data has a `data` key of its own.
        const fromMessage = root !== 'ctx' && rest[0] === MESSAGE_DATA_KEY
            && !(isRecord(base) && Object.hasOwn(base, MESSAGE_DATA_KEY));
        const keys = fromMessage ? rest.slice(1) : rest;
        const at = fromMessage ? [root, MESSAGE_DATA_KEY] : [root];
        let value = base;
        for (const key of keys) {
            if (!isRecord(value) || !Object.hasOwn(value, key)) return { found: false, at, keys: holds(value), missing: key };
            value = value[key];
            at.push(key);
        }
        return { found: true, value, at, keys: holds(value) };
    };
    const pathText = (at) => at.map((key, index) => (index && /^\d+$/.test(key) ? `.${key}` : `${index ? '.' : ''}${key}`)).join('');
    return {
        // Under a ref the caller named (the verification names it before the
        // loop fits the data into the message that carries the ref).
        put: (ref, data) => { results.set(ref, data); },
        resolve: (path) => {
            const reached = walk(path);
            return reached.found ? formatReferenceValue(reached.value) : null;
        },
        explain: (path) => {
            if (!NAME_ONLY.test(path)) return { reason: 'not_a_reference', at: '', keys: [] };
            const reached = walk(path);
            if (reached.noResult) return { reason: 'no_result', at: '', keys: [...results.keys()] };
            if (!reached.found) return { reason: 'no_key', at: pathText(reached.at), keys: reached.keys, missing: reached.missing };
            if (formatReferenceValue(reached.value) == null) return { reason: 'not_a_value', at: pathText(reached.at), keys: reached.keys };
            return null;
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
        if (!REFERENCE_START.test(name)) return whole;
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

/** Whether a text in double braces is a reference ({{r2.…}}, {{ctx.…}}), filled or not. */
export const isReferenceName = (inner) => REFERENCE_START.test(String(inner ?? '').trim());
