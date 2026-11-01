// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Bond Geometry page's inputs as the AI Copilot's change tool reaches
// them (set_bond_geometry, copilotChangeTools.js), and the record of the
// changes it made. The page registers a controller while it is mounted
// (registerBondGeometry, from BondGeometryPage.jsx); the tool waits for it
// (waitForBondGeometry), reads the page through it, and applies a change the
// user approved through the page's own setters — as the user's typing in the
// boxes and click on Compute would. One page, one controller:
//
//   controller.state() → {
//       noRun,        the static build with no run open (nothing to compute)
//       ready,        a run's elements are loaded and the triplet is seeded
//       sitesError,   why the run's sites could not be read, or null
//       elements,     the run's elements, as the page's selects list them
//       inputs,       the boxes (bondGeometryInputs.js), as typed
//       computing,    a computation is running
//       shown,        the engine request of the result on screen, or null
//       result,       { angleCount, apexCount } of that result, or null
//       epoch         the run and configuration the page shows (a new run or
//                     a Live Data save bumps it)
//   }
//   controller.stamp() → what an approval is checked against: the epoch, the
//       boxes, the result shown and whether one computes. A card shown at one
//       stamp is not applied at another (the stale-state check).
//   controller.apply(inputs, { compute = true, clear = false })
//       → Promise<{ computed, result?, error?, reason? }>
//       sets the boxes; with `compute`, computes once they hold exactly those
//       inputs and the page can (Compute, as the user's click), and resolves
//       with the result; with `clear`, takes the result off the page (an Undo
//       back to a page that showed none). Resolves `computed: false` with a
//       reason when the boxes moved on first, another run opened or the page
//       closed.
//
// The changes are kept in a history (createChangeHistory): each one the
// agent's (`actor: 'agent'`), with the inputs before and after and the epoch,
// so an Undo can put the earlier inputs back while the page still holds the
// agent's.

let current = null;

/** Register the page's controller; returns the unregister call. */
export const registerBondGeometry = (controller) => {
    current = controller;
    return () => {
        if (current === controller) current = null;
    };
};

/** The page's controller while it is mounted, else null. */
export const getBondGeometry = () => current;

const abortError = () => new DOMException('Aborted', 'AbortError');

/** The controller once the page has mounted; null after `timeoutMs`. Rejects with an AbortError on abort. */
export const waitForBondGeometry = ({ timeoutMs = 6000, signal = null, pollMs = 50 } = {}) => new Promise((resolve, reject) => {
    const started = Date.now();
    const poll = () => {
        if (signal?.aborted) return reject(abortError());
        if (current) return resolve(current);
        if (Date.now() - started >= timeoutMs) return resolve(null);
        setTimeout(poll, pollMs);
    };
    poll();
});

/**
 * The record of the agent's changes to the inputs: `record({ before, after,
 * epoch, shownBefore })` keeps one (`{ id, actor: 'agent', at, before, after,
 * epoch, shownBefore, undone: false }`), `get(id)` finds it, `markUndone(id)`
 * marks it undone, `entries()` lists copies. Nothing else: the minimal history
 * an Undo needs.
 */
export const createChangeHistory = ({ now = () => Date.now() } = {}) => {
    const entries = [];
    return {
        record: ({ before, after, epoch, shownBefore = false }) => {
            const entry = {
                id: `bond-geometry-${entries.length + 1}`,
                actor: 'agent',
                at: now(),
                before,
                after,
                epoch,
                shownBefore: Boolean(shownBefore),
                undone: false
            };
            entries.push(entry);
            return { ...entry };
        },
        get: (id) => {
            const entry = entries.find((item) => item.id === id);
            return entry ? { ...entry } : null;
        },
        markUndone: (id) => {
            const entry = entries.find((item) => item.id === id);
            if (entry) entry.undone = true;
            return Boolean(entry);
        },
        entries: () => entries.map((entry) => ({ ...entry }))
    };
};
