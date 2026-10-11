// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Streamed tool calls → executable calls. OpenAI-compatible servers stream a
// call as fragments under `delta.tool_calls`: the first carries `index`, `id`
// and the function name, the rest append pieces of the JSON `arguments`
// string. Servers differ: some omit `index` and send each call whole (Gemini's
// compat layer), some reuse one index for consecutive calls. The assembler
// below keys on `index`, falls back to `id`, then to arrival order, and starts
// a new call whenever a fragment names a function with an id that differs from
// the call it would otherwise extend. Pure functions only — no React, no I/O.

const newCall = (fragment, order) => ({
    key: order,
    index: fragment.index,
    id: fragment.id || null,
    name: fragment.name || '',
    arguments: fragment.arguments || '',
    extraContent: fragment.extraContent || null
});

// Whether `fragment` opens a call distinct from `call` although it would
// otherwise extend it: it names a function, and either both carry ids that
// differ, or it carries no id and names a different function. (A server that
// repeats the same name on every fragment, without ids, is read as one call.)
const startsAnotherCall = (call, fragment) => {
    if (!fragment.name || !call.name) return false;
    if (fragment.id && call.id) return fragment.id !== call.id;
    return !fragment.id && fragment.name !== call.name;
};

/**
 * Fold one `{ index, id, name, arguments }` fragment (any part may be null)
 * into `calls` (an array this function mutates and returns).
 */
export const accumulateToolCall = (calls, fragment) => {
    if (!fragment || typeof fragment !== 'object') return calls;
    const last = calls[calls.length - 1];
    let target = null;
    if (Number.isInteger(fragment.index)) {
        // The newest call holding this index: a reused index after a finished
        // call opens a new one only when startsAnotherCall says so.
        for (let i = calls.length - 1; i >= 0; i -= 1) {
            if (calls[i].index === fragment.index) { target = calls[i]; break; }
        }
    } else if (fragment.id) {
        target = calls.find((call) => call.id === fragment.id) || null;
    } else if (!fragment.name && last) {
        // A bare argument fragment continues the newest call.
        target = last;
    }
    if (!target || startsAnotherCall(target, fragment)) {
        calls.push(newCall(fragment, calls.length));
        return calls;
    }
    if (fragment.id && !target.id) target.id = fragment.id;
    if (fragment.name && !target.name) target.name = fragment.name;
    if (fragment.arguments) target.arguments += fragment.arguments;
    if (fragment.extraContent) target.extraContent = { ...(target.extraContent || {}), ...fragment.extraContent };
    return calls;
};

/**
 * The assembled calls as `{ id, name, arguments, extraContent }`, in arrival
 * order. A call without a function name cannot be run and is dropped; a
 * missing id gets `call_<turnPrefix><n>` so every call can be answered by a
 * tool message. `extraContent` is the provider's `extra_content` (Gemini's
 * thought_signature), null when there is none; it must be echoed back.
 */
export const finalizeToolCalls = (calls, idPrefix = '') => calls
    .filter((call) => call.name)
    .map((call, position) => ({
        id: call.id || `call_${idPrefix}${position}`,
        name: call.name,
        arguments: call.arguments,
        extraContent: call.extraContent || null
    }));

/**
 * Parse a call's `arguments` string into an object. Empty → `{}`. A JSON
 * object wrapped in a Markdown code fence is accepted (small models do that).
 * Returns `{ ok: true, value }` or `{ ok: false, error }`.
 */
export const parseToolArguments = (text) => {
    if (text == null) return { ok: true, value: {} };
    if (typeof text === 'object' && !Array.isArray(text)) return { ok: true, value: text };
    let source = String(text).trim();
    if (!source) return { ok: true, value: {} };
    const fenced = /^```(?:json)?\s*([\s\S]*?)\s*```$/i.exec(source);
    if (fenced) source = fenced[1];
    let value;
    try {
        value = JSON.parse(source);
    } catch (error) {
        return { ok: false, error: `arguments are not valid JSON (${error.message})` };
    }
    if (value === null || typeof value !== 'object' || Array.isArray(value)) {
        return { ok: false, error: 'arguments must be a JSON object' };
    }
    return { ok: true, value };
};

const describeType = (value) => (value === null ? 'null' : Array.isArray(value) ? 'array' : typeof value);

// One property against its schema. Numeric strings are accepted for number /
// integer properties ("2.5" is what a small model often sends); anything else
// of the wrong type is an error the model can read and fix.
const checkProperty = (key, raw, spec) => {
    let value = raw;
    const type = spec?.type;
    if ((type === 'number' || type === 'integer') && typeof value === 'string' && value.trim() !== '') {
        const parsed = Number(value);
        if (Number.isFinite(parsed)) value = parsed;
    }
    if (type === 'number' || type === 'integer') {
        if (typeof value !== 'number' || !Number.isFinite(value)) {
            return { error: `${key} must be a finite number, got ${describeType(raw)} ${JSON.stringify(raw)}` };
        }
        if (type === 'integer' && !Number.isInteger(value)) {
            return { error: `${key} must be an integer, got ${value}` };
        }
        if (Number.isFinite(spec.minimum) && value < spec.minimum) {
            return { error: `${key} must be ≥ ${spec.minimum}, got ${value}` };
        }
        if (Number.isFinite(spec.maximum) && value > spec.maximum) {
            return { error: `${key} must be ≤ ${spec.maximum}, got ${value}` };
        }
    } else if (type === 'string') {
        if (typeof value === 'number' && Number.isFinite(value)) value = String(value);
        if (typeof value !== 'string') {
            return { error: `${key} must be a string, got ${describeType(raw)}` };
        }
    } else if (type === 'boolean') {
        if (value === 'true' || value === 'false') value = value === 'true';
        if (typeof value !== 'boolean') {
            return { error: `${key} must be true or false, got ${JSON.stringify(raw)}` };
        }
    } else if (type === 'array') {
        return checkArray(key, value, spec);
    } else if (type === 'object') {
        return checkObject(key, value, spec);
    }
    if (Array.isArray(spec?.enum) && !spec.enum.includes(value)) {
        return { error: `${key} must be one of ${spec.enum.map((item) => JSON.stringify(item)).join(', ')}, got ${JSON.stringify(value)}` };
    }
    return { value };
};

// A JSON text where an array or object belongs: small models often send the
// list as a string ("[\"a\", \"b\"]"), or a bare comma list ("3, 7").
const fromText = (value, kind) => {
    if (typeof value !== 'string') return value;
    const text = value.trim();
    try {
        const parsed = JSON.parse(text);
        if (kind === 'array' ? Array.isArray(parsed) : (parsed && typeof parsed === 'object' && !Array.isArray(parsed))) {
            return parsed;
        }
    } catch {
        // Not JSON: a bare list below, else the type error.
    }
    if (kind === 'array' && text && !/^[[{]/.test(text)) return text.split(/\s*,\s*/).filter(Boolean);
    return value;
};

// An array property: each item against `items` (one schema for every item),
// within minItems / maxItems. One bare item where a list belongs is read as
// a list of one.
function checkArray(key, raw, spec) {
    let value = fromText(raw, 'array');
    if (!Array.isArray(value) && value !== undefined && typeof value !== 'object') value = [value];
    if (!Array.isArray(value)) return { error: `${key} must be a list, got ${describeType(raw)}` };
    if (Number.isFinite(spec.minItems) && value.length < spec.minItems) {
        return { error: `${key} must hold at least ${spec.minItems} item(s), got ${value.length}` };
    }
    if (Number.isFinite(spec.maxItems) && value.length > spec.maxItems) {
        return { error: `${key} may hold at most ${spec.maxItems} item(s), got ${value.length}` };
    }
    if (!spec.items) return { value };
    const items = [];
    for (const [index, item] of value.entries()) {
        const checked = checkProperty(`${key}[${index}]`, item, spec.items);
        if (checked.error) return checked;
        items.push(checked.value);
    }
    return { value: items };
}

// An object property: its `properties` and `required`, like the arguments
// themselves (unknown keys are dropped).
function checkObject(key, raw, spec) {
    const value = fromText(raw, 'object');
    if (!value || typeof value !== 'object' || Array.isArray(value)) {
        return { error: `${key} must be an object, got ${describeType(raw)}` };
    }
    if (!spec.properties) return { value };
    const out = {};
    for (const [name, item] of Object.entries(value)) {
        if (!Object.prototype.hasOwnProperty.call(spec.properties, name) || item == null) continue;
        const checked = checkProperty(`${key}.${name}`, item, spec.properties[name]);
        if (checked.error) return checked;
        out[name] = checked.value;
    }
    const missing = (spec.required || []).find((name) => !(name in out));
    if (missing) return { error: `${key}.${missing} is required` };
    return { value: out };
}

/**
 * Check parsed arguments against the JSON-Schema subset the tools use
 * (`properties` with type / enum / minimum / maximum, arrays with `items` /
 * minItems / maxItems, objects with their own `properties`, and `required`).
 * Unknown keys are dropped with a note; `null` for an optional key means
 * "not given". Returns `{ ok, value, errors, notes }`.
 */
export const validateToolArguments = (schema, args) => {
    const properties = schema?.properties || {};
    const required = Array.isArray(schema?.required) ? schema.required : [];
    const value = {};
    const errors = [];
    const notes = [];
    Object.entries(args || {}).forEach(([key, raw]) => {
        if (!Object.prototype.hasOwnProperty.call(properties, key)) {
            notes.push(`ignored unknown argument "${key}"`);
            return;
        }
        if (raw === null || raw === undefined) return;
        const checked = checkProperty(key, raw, properties[key]);
        if (checked.error) errors.push(checked.error);
        else value[key] = checked.value;
    });
    required.forEach((key) => {
        if (!(key in value) && !errors.some((error) => error.startsWith(`${key} `))) {
            errors.push(`${key} is required`);
        }
    });
    return { ok: errors.length === 0, value, errors, notes };
};

/** Key-sorted JSON, so the same arguments in another order compare equal. */
export const canonicalJson = (value) => {
    if (Array.isArray(value)) return `[${value.map(canonicalJson).join(',')}]`;
    if (value && typeof value === 'object') {
        return `{${Object.keys(value).sort().map((key) => `${JSON.stringify(key)}:${canonicalJson(value[key])}`).join(',')}}`;
    }
    return JSON.stringify(value);
};
