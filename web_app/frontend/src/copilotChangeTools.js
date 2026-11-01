// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The AI Copilot's change tool: set_bond_geometry, the one analysis input the
// Workbench lets a user change, set as the user would (a change tool in the
// llm tool contract: effect 'change', llm/copilot/effects.js). It reads the
// Bond Geometry page through the controller the page registers
// (bondGeometryControl.js), shows the user the change on an approval card
// (the turn's approve port: the dock's card, or at once in Auto), checks that
// the page is still as the card described it (the stale-state check), then
// applies it through the page's own setters and lets the page compute, as
// Compute does. Each applied change is the agent's entry in the change
// history, and `undo` puts the earlier inputs back while the page still holds
// the agent's.
//
// It is refused, in plain words, when no run is open, when the page cannot
// read the run, when arguments contradict each other or leave the page unable
// to compute (bondGeometryInputs.js changedInputs), and — after the approval —
// when the page moved while the card waited. The validation's live arms never
// get it: they offer compute-only tools (validation/harness.js liveArmTools).

import { getBondGeometry, waitForBondGeometry } from './bondGeometryControl';
import {
    BOND_GEOMETRY_LIMITS, changedInputs, describeChange, inputsLine, requestOf, sameInputs, sameRequest, toolInputs
} from './bondGeometryInputs';
import { elementProperty } from './copilotTools';

export const SET_BOND_GEOMETRY = 'set_bond_geometry';
// The card the change is about: the page opens on it, and it breathes while the change runs.
export const BOND_GEOMETRY_CARD = 'bond-angles';

/**
 * How long the tool waits, in ms: for the page to mount (a first visit), for
 * the run's elements to load or a running computation to end, and for the
 * page to compute the change (the loop's own per-call limit; the approval
 * itself has none).
 */
export const CHANGE_WAITS = Object.freeze({ page: 6000, ready: 30000, compute: 120000 });

const { windowMinA, windowMaxA, binMinDeg, binMaxDeg } = BOND_GEOMETRY_LIMITS;
const windowBound = (description) => ({ type: 'number', minimum: windowMinA, maximum: windowMaxA, description });

/** The shape of set_bond_geometry's data (llm copilot/dataShape.js). */
export const SET_BOND_GEOMETRY_SHAPE = Object.freeze({
    'applied?': 'boolean',
    'declined?': 'boolean',
    'unchanged?': 'boolean',
    'stale?': 'boolean',
    'inputs?': {
        end1: 'string',
        apex: 'string',
        end2: 'string',
        r12_min: 'number|null',
        r12_max: 'number|null',
        distinct_bc: 'boolean',
        'r23_min?': 'number|null',
        'r23_max?': 'number|null',
        bin_width: 'number|null'
    },
    'computed?': 'boolean',
    'angle_count?': 'number',
    'not_computed?': 'string'
});

const abortError = () => new DOMException('Aborted', 'AbortError');
const sleep = (ms) => new Promise((resolve) => { setTimeout(resolve, ms); });

// Poll `test` until it holds (true) or `timeoutMs` passes (false); an AbortError on abort.
const waitUntil = async (test, { timeoutMs, signal, pollMs = 50 }) => {
    const started = Date.now();
    for (;;) {
        if (signal?.aborted) throw abortError();
        if (test()) return true;
        if (Date.now() - started >= timeoutMs) return false;
        await sleep(pollMs);
    }
};

// The page's computation, or a "not within" outcome after `ms` (the change stands either way).
const withinTime = (promise, ms) => Promise.race([
    promise,
    sleep(ms).then(() => ({ computed: false, reason: `the page had not finished computing within ${Math.round(ms / 1000)} s` }))
]);

const seconds = (ms) => `${Math.round(ms / 1000)} s`;

// The step-card title from what the call names.
const changeLabel = ({ end1, apex, end2, r12_min: lo, r12_max: hi }) => {
    const parts = ['Set Bond Geometry'];
    if (end1 && apex && end2) parts.push(`${end1}–${apex}–${end2}`);
    if (Number.isFinite(lo) && Number.isFinite(hi)) parts.push(`${lo}–${hi} Å`);
    return parts.join(' · ');
};

// What the model reads once a change was applied, and how the page took it.
const appliedResult = (after, outcome) => {
    const line = inputsLine(after);
    if (outcome.computed) {
        const count = outcome.result?.angleCount ?? 0;
        return {
            data: { applied: true, inputs: toolInputs(after), computed: true, angle_count: count },
            checks: [count > 0
                ? { level: 'ok', text: `Applied: ${line}; the page computed ${count} angles.` }
                : { level: 'warn', text: `Applied: ${line}; the page found no angles in these windows: widen them.` }]
        };
    }
    const why = outcome.error || outcome.reason || 'it did not compute';
    return {
        data: { applied: true, inputs: toolInputs(after), computed: false, not_computed: why },
        checks: [{ level: 'warn', text: `Applied: ${line}, but the page did not compute it: ${why}.` }]
    };
};

/**
 * The change tools for one run: set_bond_geometry, offered when the run has
 * elements (a configuration to compute on). `history` is the app's change
 * history (bondGeometryControl.js createChangeHistory: one for the app, so an
 * Undo outlives a rebuilt tool list); `openPage(key)` opens a page when the
 * turn has no card focus to do it; `waits` overrides CHANGE_WAITS (tests).
 */
export const buildChangeTools = ({ elements = [], history, openPage = null, waits = {} } = {}) => {
    if (!elements.length) return [];
    if (!history || typeof history.record !== 'function') throw new TypeError('buildChangeTools: history must be a change history (createChangeHistory)');
    const wait = { ...CHANGE_WAITS, ...waits };

    return [{
        name: SET_BOND_GEOMETRY,
        description: 'Set the Bond Geometry page\'s inputs as the user would: the triplet A–B–C (B the central atom), '
            + 'the A–B bond window, a distinct B–C window, the angle bin width. Controls you leave out keep their values. '
            + 'The user approves the change on a card first (unless they chose Auto); then the page computes its angle '
            + 'distribution. Returns the inputs now on the page and whether it computed, or {"declined": true}. '
            + 'The user can undo it.',
        // The method advice (llm copilot/method.js): where the windows come
        // from (the coordination-shells skill), and which tool answers.
        method: 'Choose the windows as the coordination-shells skill says: the whole first coordination shell '
            + '(find_gr_peaks\' first_shell.window_A) unless the user asked for another. Use it when the user wants an '
            + 'analysis set up on the page; to answer a question, run compute_bond_angles. Do not ask leave with ask_user '
            + 'first: the card asks.',
        parameters: {
            type: 'object',
            properties: {
                end1: elementProperty(elements, 'Element A.'),
                apex: elementProperty(elements, 'Element B, the central atom.'),
                end2: elementProperty(elements, 'Element C.'),
                r12_min: windowBound('A–B window minimum (Å).'),
                r12_max: windowBound('A–B window maximum (Å).'),
                distinct_bc: {
                    type: 'boolean',
                    description: 'The Distinct B–C switch: false makes the B–C bond use the A–B window again; '
                        + 'r23_min with r23_max turns it on.'
                },
                r23_min: windowBound('B–C window minimum (Å); give both or neither.'),
                r23_max: windowBound('B–C window maximum (Å).'),
                bin_width: { type: 'number', minimum: binMinDeg, maximum: binMaxDeg, description: 'Angle bin width (degrees).' }
            },
            additionalProperties: false
        },
        dataShape: SET_BOND_GEOMETRY_SHAPE,
        label: changeLabel,
        page: 'geometry',
        cards: () => [BOND_GEOMETRY_CARD],

        run: async (args, { signal, focus, approve } = {}) => {
            if (typeof approve !== 'function') throw new Error('a change needs the user\'s approval, and nothing here can ask for it');
            // The page on screen first (it mounts on a first visit), so the
            // user sees what the card will change.
            if (focus) await focus.show([BOND_GEOMETRY_CARD]);
            else openPage?.('geometry');
            const controller = await waitForBondGeometry({ timeoutMs: wait.page, signal });
            if (!controller) throw new Error('the Bond Geometry page did not open');
            await waitUntil(() => {
                const page = controller.state();
                return page.noRun || Boolean(page.sitesError) || (page.ready && !page.computing);
            }, { timeoutMs: wait.ready, signal });
            const before = controller.state();
            if (before.noRun) throw new Error('no run is open: open a run folder, or load the Demo');
            if (before.sitesError) throw new Error(`the Bond Geometry page could not read the run: ${before.sitesError}`);
            if (!before.ready) throw new Error(`the Bond Geometry page had not loaded the run's elements within ${seconds(wait.ready)}`);
            if (before.computing) throw new Error(`the Bond Geometry page was still computing after ${seconds(wait.ready)}`);
            // Refused in plain words when the arguments contradict each other
            // or leave the page unable to compute.
            const after = changedInputs(before.inputs, args, { elements: before.elements });
            if (sameInputs(after, before.inputs) && sameRequest(before.shown, requestOf(after))) {
                const line = inputsLine(after);
                return {
                    data: { unchanged: true, inputs: toolInputs(after) },
                    checks: [{ level: 'ok', text: `Already on the page: ${line}.` }],
                    change: { status: 'unchanged', line, details: [], auto: false, id: null }
                };
            }

            const { line, details } = describeChange(before.inputs, after);
            const stamp = controller.stamp();
            const decision = await approve({
                tool: SET_BOND_GEOMETRY,
                title: 'Change the Bond Geometry inputs',
                line,
                details: [...details, 'Then the page computes the angle distribution; Undo puts the inputs back.']
            }, { signal });
            const auto = decision?.auto === true;
            if (!decision?.approved) {
                return {
                    data: { declined: true },
                    checks: [{ level: 'note', text: 'The user declined this change. Ask what they would prefer, or propose something else.' }],
                    change: { status: 'declined', line, details, auto, id: null }
                };
            }
            if (signal?.aborted) throw abortError();
            // The card may have waited: the page as it is now (another page
            // registered, an edit, a computation, a new run or configuration).
            const now = getBondGeometry();
            if (now !== controller || controller.stamp() !== stamp) {
                const inputs = (now || controller).state().inputs;
                return {
                    data: { stale: true, inputs: toolInputs(inputs) },
                    checks: [{
                        level: 'fail',
                        text: 'Not applied: the Bond Geometry page changed while the card waited, so what was approved no longer '
                            + 'describes it. The inputs now on the page are in the data; propose the change again if it still applies.'
                    }],
                    change: { status: 'stale', line, details, auto, id: null }
                };
            }
            const entry = history.record({ before: before.inputs, after, epoch: before.epoch, shownBefore: Boolean(before.shown) });
            const outcome = await withinTime(controller.apply(after), wait.compute);
            return { ...appliedResult(after, outcome), change: { status: 'applied', line, details, auto, id: entry.id } };
        },

        // Undo one applied change: the earlier inputs back, through the same
        // controller, while the page still holds the agent's inputs on the same
        // run and configuration; computed again when the page showed a result
        // before the change, else the agent's result taken off the page.
        // Resolves { ok: true, computed } or { ok: false, reason }.
        undo: async (change) => {
            const entry = change?.id ? history.get(change.id) : null;
            if (!entry) return { ok: false, reason: 'this change is not in the history' };
            if (entry.undone) return { ok: false, reason: 'it was undone already' };
            const controller = getBondGeometry();
            if (!controller) return { ok: false, reason: 'the Bond Geometry page is not open' };
            const page = controller.state();
            if (page.epoch !== entry.epoch) return { ok: false, reason: 'another run or configuration is open' };
            if (!sameInputs(page.inputs, entry.after)) return { ok: false, reason: 'the inputs changed on the page since' };
            if (page.computing) return { ok: false, reason: 'the page is computing: try again when it is done' };
            history.markUndone(entry.id);
            const outcome = await withinTime(
                controller.apply(entry.before, entry.shownBefore ? { compute: true } : { compute: false, clear: true }),
                wait.compute
            );
            return { ok: true, computed: Boolean(outcome.computed) };
        }
    }];
};
