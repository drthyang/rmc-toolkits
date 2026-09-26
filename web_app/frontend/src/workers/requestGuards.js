// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Request and result guards shared by the static-mode workers, mirroring the
// Flask boundary (web_app/backend/app.py) so both runtimes accept and refuse
// the same requests:
//   requestNumber      <- _number() / _query_number(): a finite number, integral
//                         where required, inside its range (blank = default)
//   hasNonFinite /
//   assertFiniteResult <- _strict_result_response(): a computed result holding
//                         NaN or Infinity is an error, never posted as success
//   requestObject      <- every worker: a null / non-object message is an error
//                         the caller receives, so its promise always settles.

const DECIMAL = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/;
const NON_FINITE_TOKEN = /^[+-]?(nan|inf|infinity)$/i;

// The value as the error names it: strings quoted, like Python's repr().
const shown = (raw) => (typeof raw === 'string' ? `'${raw}'` : String(raw));

export const isBlankRequestValue = (raw) => raw == null || (typeof raw === 'string' && !raw.trim());

/**
 * One request value as a finite number, checked like app._number(): text that
 * is not a decimal number, booleans, objects, NaN and ±Infinity are refused,
 * `integer` requires an integral value, `gt/ge/lt/le` bound it and `clamp`
 * clamps it. Missing, null or blank text returns `fallback`.
 */
export const requestNumber = (raw, name, {
    fallback, integer = false, gt, ge, lt, le, clamp
} = {}) => {
    if (isBlankRequestValue(raw)) return fallback;
    let value;
    if (typeof raw === 'number') {
        value = raw;
    } else if (typeof raw === 'string' && (DECIMAL.test(raw.trim()) || NON_FINITE_TOKEN.test(raw.trim()))) {
        value = Number(raw.trim().replace(/^([+-]?)inf$/i, '$1Infinity'));
    } else {
        throw new Error(`${name} must be a number, got ${shown(raw)}`);
    }
    if (!Number.isFinite(value)) throw new Error(`${name} must be a finite number, got ${shown(raw)}`);
    if (integer && !Number.isInteger(value)) throw new Error(`${name} must be an integer, got ${shown(raw)}`);
    if (gt !== undefined && !(value > gt)) throw new Error(`${name} must be > ${gt}, got ${value}`);
    if (ge !== undefined && !(value >= ge)) throw new Error(`${name} must be >= ${ge}, got ${value}`);
    if (lt !== undefined && !(value < lt)) throw new Error(`${name} must be < ${lt}, got ${value}`);
    if (le !== undefined && !(value <= le)) throw new Error(`${name} must be <= ${le}, got ${value}`);
    if (clamp) value = Math.min(Math.max(value, clamp[0]), clamp[1]);
    return value;
};

/** The worker message's payload object; anything else is an error. */
export const requestObject = (data) => {
    if (data === null || typeof data !== 'object' || Array.isArray(data)) {
        throw new Error(`worker request must be an object, got ${data === null ? 'null' : typeof data}`);
    }
    return data;
};

/** True when a number anywhere in `value` (objects, arrays, typed arrays) is NaN or ±Infinity. */
export const hasNonFinite = (value) => {
    if (typeof value === 'number') return !Number.isFinite(value);
    if (value === null || typeof value !== 'object') return false;
    if (ArrayBuffer.isView(value)) {
        if (value instanceof DataView) return false;
        for (let i = 0; i < value.length; i += 1) {
            if (!Number.isFinite(value[i])) return true;
        }
        return false;
    }
    if (Array.isArray(value)) return value.some(hasNonFinite);
    return Object.values(value).some(hasNonFinite);
};

// app._strict_result_response's message, word for word.
export const NON_FINITE_RESULT_MESSAGE = 'the result contains NaN or Infinity for these parameters (an extreme '
    + 'bandwidth, scale or extent?); use less extreme values';

/** Throw the Flask 400 message when a computed result holds NaN or Infinity. */
export const assertFiniteResult = (result) => {
    if (hasNonFinite(result)) throw new Error(NON_FINITE_RESULT_MESSAGE);
    return result;
};
