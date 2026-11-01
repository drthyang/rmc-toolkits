// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Bond Geometry page's inputs: what its controls let the user set, as the
// page holds them (strings, as typed in its boxes) and as the AI Copilot's
// change tool reads and sets them (set_bond_geometry, copilotChangeTools.js).
// Pure: no React, no page.
//
//   inputs = { end1, apex, end2,        the triplet A–B–C, B central (selects)
//              r12Min, r12Max,          the A–B window, Å (boxes)
//              split23,                 the Distinct B–C switch
//              r23Min, r23Max,          the B–C window, Å (boxes; read only when split)
//              binWidth }               the histogram bin, degrees (box)
//
// The tool's arguments name the same controls (toolInputs): a change sets the
// controls it names and keeps the others as the page has them, and is refused
// in plain words when the page would then hold something it cannot compute.

import { tripletRequestFromInputs } from './workers/triplets';

/**
 * The ranges of the page's own boxes (their min/max): a window bound in
 * [0, 15] Å (15 Å is the app's cap on a window maximum, the worker's and
 * /api/triplets'), a bin width in [0.1, 45]°.
 */
export const BOND_GEOMETRY_LIMITS = Object.freeze({
    windowMinA: 0,
    windowMaxA: 15,
    windowStepA: 0.05,
    binMinDeg: 0.1,
    binMaxDeg: 45,
    binStepDeg: 0.5
});

/** The page's controls before a run seeds the triplet. */
export const DEFAULT_INPUTS = Object.freeze({
    end1: '',
    apex: '',
    end2: '',
    r12Min: '2.00',
    r12Max: '3.00',
    split23: false,
    r23Min: '2.00',
    r23Max: '3.00',
    binWidth: '1.0'
});

export const INPUT_KEYS = Object.freeze(Object.keys(DEFAULT_INPUTS));

/** The inputs with every key, as the page's state holds them (strings, a boolean switch). */
export const normalizeInputs = (inputs) => Object.fromEntries(INPUT_KEYS.map((key) => {
    const value = inputs?.[key] ?? DEFAULT_INPUTS[key];
    return [key, key === 'split23' ? Boolean(value) : String(value)];
}));

/** Whether two sets of inputs are the same, box for box (what a user's edit changes). */
export const sameInputs = (one, two) => INPUT_KEYS.every((key) => normalizeInputs(one)[key] === normalizeInputs(two)[key]);

/** The engine request the inputs make (the page's), or null while a box is blank or not a number. */
export const requestOf = (inputs) => {
    try {
        return tripletRequestFromInputs(normalizeInputs(inputs));
    } catch {
        return null;
    }
};

const REQUEST_KEYS = ['end1', 'apex', 'end2', 'r12Min', 'r12Max', 'r23Min', 'r23Max', 'binWidth'];

/** Whether two engine requests are the same analysis (the page's stale test). */
export const sameRequest = (one, two) => Boolean(one && two) && REQUEST_KEYS.every((key) => one[key] === two[key]);

// A bound at 2–4 decimals, as the page labels it: 3 → "3.00", 2.205 → "2.205".
const formatBound = (value) => {
    const number = Number(value);
    if (!Number.isFinite(number)) return String(value).trim() || '…';
    const decimals = (String(value).split('.')[1] ?? '').length;
    return number.toFixed(Math.min(4, Math.max(2, decimals)));
};
// A bin width as typed, without trailing zeros: "1.0" → "1", "0.50" → "0.5".
const formatBin = (value) => {
    const number = Number(value);
    return Number.isFinite(number) ? String(Math.round(number * 1e6) / 1e6) : (String(value).trim() || '…');
};
const formatWindow = (min, max) => `${formatBound(min)}–${formatBound(max)} Å`;

/** The triplet as the page names it: "Se–Ta–Se". */
export const tripletText = (inputs) => {
    const { end1, apex, end2 } = normalizeInputs(inputs);
    return `${end1 || 'A'}–${apex || 'B'}–${end2 || 'C'}`;
};

/** The bond windows as the page names them: "2.21–2.65 Å", or "A–B 2.21–2.65 Å, B–C 2.40–2.80 Å". */
export const windowsText = (inputs) => {
    const value = normalizeInputs(inputs);
    return value.split23
        ? `A–B ${formatWindow(value.r12Min, value.r12Max)}, B–C ${formatWindow(value.r23Min, value.r23Max)}`
        : formatWindow(value.r12Min, value.r12Max);
};

/** The whole set in one plain line: "Se–Ta–Se · 2.21–2.65 Å · 1° bins". */
export const inputsLine = (inputs) => `${tripletText(inputs)} · ${windowsText(inputs)} · ${formatBin(normalizeInputs(inputs).binWidth)}° bins`;

/**
 * What a change does, for its approval card: `line`, the page after it in one
 * plain line (the result first), and `details`, one line per control it
 * changes, before → after (or that it changes none and only computes).
 */
export const describeChange = (before, after) => {
    const from = normalizeInputs(before);
    const to = normalizeInputs(after);
    const details = [];
    if (tripletText(from) !== tripletText(to)) details.push(`Triplet: ${tripletText(from)} → ${tripletText(to)}`);
    if (formatWindow(from.r12Min, from.r12Max) !== formatWindow(to.r12Min, to.r12Max)) {
        details.push(`A–B window: ${formatWindow(from.r12Min, from.r12Max).replace(' Å', '')} → ${formatWindow(to.r12Min, to.r12Max)}`);
    }
    if (from.split23 !== to.split23) details.push(`Distinct B–C window: ${from.split23 ? 'on' : 'off'} → ${to.split23 ? 'on' : 'off'}`);
    if (to.split23 && (!from.split23 || formatWindow(from.r23Min, from.r23Max) !== formatWindow(to.r23Min, to.r23Max))) {
        details.push(`B–C window: ${formatWindow(from.r23Min, from.r23Max).replace(' Å', '')} → ${formatWindow(to.r23Min, to.r23Max)}`);
    }
    if (formatBin(from.binWidth) !== formatBin(to.binWidth)) details.push(`Bin width: ${formatBin(from.binWidth)}° → ${formatBin(to.binWidth)}°`);
    if (!details.length) details.push('No input changes: the page computes these inputs.');
    return { line: inputsLine(to), details };
};

// The page's box names, for a refusal that names the box.
const BOX = {
    r12Min: 'A–B window minimum',
    r12Max: 'A–B window maximum',
    r23Min: 'B–C window minimum',
    r23Max: 'B–C window maximum',
    binWidth: 'bin width'
};
const ARGUMENT = { r12Min: 'r12_min', r12Max: 'r12_max', r23Min: 'r23_min', r23Max: 'r23_max', binWidth: 'bin_width' };

/**
 * Why the page could not compute these inputs, in plain words naming the box
 * (and the argument that sets it), or null when it can: every element one of
 * the run's, every box a number in its range, each window's minimum below its
 * maximum.
 */
export const inputsProblem = (inputs, { elements = [] } = {}) => {
    const value = normalizeInputs(inputs);
    for (const [key, label] of [['end1', 'A'], ['apex', 'B'], ['end2', 'C']]) {
        if (!value[key]) return `${label} names no element yet: give ${key}`;
        if (elements.length && !elements.includes(value[key])) {
            return `${label} is ${value[key]}, which is not an element of this run (${elements.join(', ')}): give ${key}`;
        }
    }
    const { windowMinA, windowMaxA, binMinDeg, binMaxDeg } = BOND_GEOMETRY_LIMITS;
    const numbers = {};
    const boxes = ['r12Min', 'r12Max', ...(value.split23 ? ['r23Min', 'r23Max'] : []), 'binWidth'];
    for (const key of boxes) {
        const text = value[key].trim();
        const number = text === '' ? NaN : Number(text);
        if (text === '') return `the ${BOX[key]} box is empty: give ${ARGUMENT[key]}`;
        if (!Number.isFinite(number)) return `the ${BOX[key]} box holds "${text}", not a number: give ${ARGUMENT[key]}`;
        numbers[key] = number;
    }
    for (const [lo, hi, name] of [['r12Min', 'r12Max', 'A–B'], ...(value.split23 ? [['r23Min', 'r23Max', 'B–C']] : [])]) {
        if (numbers[lo] < windowMinA || numbers[hi] > windowMaxA) {
            return `the ${name} window ${formatWindow(numbers[lo], numbers[hi])} leaves the page's range, ${windowMinA}–${windowMaxA} Å`;
        }
        if (numbers[lo] >= numbers[hi]) {
            return `the ${name} window would be ${formatWindow(numbers[lo], numbers[hi])}: its minimum must be below its maximum`;
        }
    }
    if (numbers.binWidth < binMinDeg || numbers.binWidth > binMaxDeg) {
        return `a ${formatBin(numbers.binWidth)}° bin width leaves the page's range, ${binMinDeg}–${binMaxDeg}°`;
    }
    return null;
};

/**
 * The page's inputs after a change (the tool's arguments: end1, apex, end2,
 * r12_min, r12_max, distinct_bc, r23_min, r23_max, bin_width), the controls
 * it does not name kept as `current` holds them: r23_min with r23_max turns
 * Distinct B–C on, distinct_bc false turns it off (the B–C boxes keep their
 * values, as the page's own switch does). Throws an Error naming the problem
 * for arguments that contradict each other, or for a page that could not then
 * compute (inputsProblem).
 */
export const changedInputs = (current, args, { elements = [] } = {}) => {
    const next = normalizeInputs(current);
    for (const key of ['end1', 'apex', 'end2']) if (args[key] != null) next[key] = String(args[key]);
    if (args.r12_min != null) next.r12Min = String(args.r12_min);
    if (args.r12_max != null) next.r12Max = String(args.r12_max);
    const r23Given = args.r23_min != null || args.r23_max != null;
    if ((args.r23_min == null) !== (args.r23_max == null)) throw new Error('give r23_min and r23_max together, or neither');
    if (r23Given && args.distinct_bc === false) {
        throw new Error('r23_min and r23_max set a distinct B–C window: leave distinct_bc out, or set it true');
    }
    if (r23Given) {
        next.split23 = true;
        next.r23Min = String(args.r23_min);
        next.r23Max = String(args.r23_max);
    } else if (args.distinct_bc === true && !next.split23) {
        throw new Error('give r23_min and r23_max to turn the distinct B–C window on');
    } else if (args.distinct_bc === false) {
        next.split23 = false;
    }
    if (args.bin_width != null) next.binWidth = String(args.bin_width);
    const problem = inputsProblem(next, { elements });
    if (problem) throw new Error(`the page cannot compute that: ${problem}`);
    return next;
};

/** The inputs as the tool's arguments name them, numbers as numbers (what the model reads back). */
export const toolInputs = (inputs) => {
    const value = normalizeInputs(inputs);
    const number = (text) => {
        const parsed = Number(text);
        return text.trim() !== '' && Number.isFinite(parsed) ? parsed : null;
    };
    return {
        end1: value.end1,
        apex: value.apex,
        end2: value.end2,
        r12_min: number(value.r12Min),
        r12_max: number(value.r12Max),
        distinct_bc: value.split23,
        ...(value.split23 ? { r23_min: number(value.r23Min), r23_max: number(value.r23Max) } : {}),
        bin_width: number(value.binWidth)
    };
};
