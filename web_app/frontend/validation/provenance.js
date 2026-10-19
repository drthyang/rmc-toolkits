// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Where each number in an answer came from. A checkable number is one with a
// decimal point or a unit (Å, %, σ, °, ×): counts, indices and Miller
// components are left out. Each is matched to the values the model was
// told (each tool message it read: the data as fitted, the checks, the error
// and the note; the run context as shown; the question) at the precision it
// was written, a percentage also to a fraction (28 % ↔ 0.28); an unmatched
// number is "unsourced": derived by the model, or invented. A match at low
// precision (one decimal) can be a coincidence: the record names the value
// it matched, so a reviewer can check it.

const NUMBER = /(?<![\w.])([-−]?\d+(?:\.\d+)?)(\s*(?:Å|A\b|%|σ|sigma|°|deg\b|×|x\b))?/g;

const unitOf = (raw) => {
    const unit = (raw || '').trim();
    if (!unit) return null;
    if (unit === 'A') return 'Å';
    if (unit === 'sigma') return 'σ';
    if (unit === 'deg') return '°';
    if (unit === 'x') return '×';
    return unit;
};

/** The checkable numbers of a text, with a snippet around each. */
export const answerNumbers = (text) => {
    const source = String(text || '');
    const out = [];
    for (const match of source.matchAll(NUMBER)) {
        const literal = match[1].replace('−', '-');
        const unit = unitOf(match[2]);
        const decimals = literal.includes('.') ? literal.split('.')[1].length : 0;
        if (!unit && !decimals) continue;
        const start = Math.max(0, match.index - 40);
        out.push({
            value: Number(literal),
            literal,
            decimals,
            unit,
            snippet: source.slice(start, match.index + match[0].length + 20).replace(/\s+/g, ' ').trim()
        });
    }
    return out;
};

// Every number in a JSON value or text, with where it sits.
const collectValues = (value, path, out) => {
    if (typeof value === 'number' && Number.isFinite(value)) {
        out.push({ value, path });
    } else if (typeof value === 'string') {
        for (const match of value.matchAll(/(?<![\w.])-?\d+(?:\.\d+)?(?:e-?\d+)?/gi)) out.push({ value: Number(match[0]), path });
    } else if (Array.isArray(value)) {
        value.forEach((item, index) => collectValues(item, `${path}[${index}]`, out));
    } else if (value && typeof value === 'object') {
        Object.entries(value).forEach(([key, item]) => collectValues(item, path ? `${path}.${key}` : key, out));
    }
    return out;
};

/**
 * The values the model was told: `steps` are what each tool step told it
 * ({ id, name, data, checks, error, note }: a tool message's fields; the
 * caller leaves out what the model never read), `context` the run context
 * object as shown, and the question text. In that order, each step's data,
 * checks, error, then note: a number is credited to the first that holds it.
 */
export const shownValues = ({ steps = [], context = null, question = '' }) => {
    const out = [];
    steps.forEach((step) => {
        collectValues(step.data, '', []).forEach((entry) => out.push({ ...entry, source: `tool ${step.name} (${step.id})${entry.path ? ` data.${entry.path}` : ''}` }));
        (step.checks || []).forEach((check, index) => collectValues(check.text, '', []).forEach((entry) => out.push({
            ...entry, source: `tool ${step.name} (${step.id}) check ${index + 1}`
        })));
        ['error', 'note'].forEach((key) => collectValues(step[key] ?? null, '', []).forEach((entry) => out.push({
            ...entry, source: `tool ${step.name} (${step.id}) ${key}`
        })));
    });
    if (context) collectValues(context, '', []).forEach((entry) => out.push({ ...entry, source: `run context ${entry.path}` }));
    collectValues(question, '', []).forEach((entry) => out.push({ ...entry, source: 'the question' }));
    return out;
};

// Written at `decimals` places, does `shown` round to `value`?
const sameAtPrecision = (value, shown, decimals) => {
    const half = 0.5 * 10 ** -decimals;
    return Math.abs(Math.abs(value) - Math.abs(shown)) <= half * 1.0001 + 1e-12;
};

/** Each answer number with its source (or null), and a summary. */
export const traceNumbers = (text, shown) => {
    const numbers = answerNumbers(text).map((number) => {
        const match = shown.find((entry) => sameAtPrecision(number.value, entry.value, number.decimals))
            || (number.unit === '%' ? shown.find((entry) => sameAtPrecision(number.value, entry.value * 100, number.decimals)) : null);
        return { ...number, source: match ? match.source : null };
    });
    const unsourced = numbers.filter((number) => !number.source);
    return {
        numbers,
        checkable: numbers.length,
        sourced: numbers.length - unsourced.length,
        unsourced: unsourced.map((number) => ({ literal: number.literal, unit: number.unit, snippet: number.snippet }))
    };
};
