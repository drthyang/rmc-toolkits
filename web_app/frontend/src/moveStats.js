// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/moveStats.js
//
// Sampling ratios from the `.rmc6f` move counters.
//
// The raw counts say little on their own — a million moves is thorough for a small
// box and negligible for a large one — so they are reported per atom:
//
//   generatedPerAtom  = N_gen / N_atom   how many moves each atom was offered
//   acceptedPerAtom   = N_acc / N_atom   how many it actually took (sampling depth)
//   acceptedPerGenerated = N_acc / N_gen the acceptance ratio
//
// Returns null when the counters are absent: pre-6f configurations and some writers
// omit those header lines, and the Flask structure endpoint only began supplying them
// alongside this change.

const isCount = (value) => Number.isFinite(value) && value >= 0;

/**
 * @param {{generated?:number, accepted?:number}|null|undefined} moves  header counters
 * @param {number} totalAtoms
 * @returns {{generatedPerAtom?:number, acceptedPerAtom?:number,
 *            acceptedPerGenerated?:number, generated?:number, accepted?:number}|null}
 */
export function moveRatios(moves, totalAtoms) {
    if (!moves || !Number.isFinite(totalAtoms) || totalAtoms <= 0) return null;
    const { generated, accepted } = moves;
    const ratios = {};
    if (isCount(generated)) ratios.generatedPerAtom = generated / totalAtoms;
    if (isCount(accepted)) ratios.acceptedPerAtom = accepted / totalAtoms;
    // Guard the denominator separately: a run can legitimately have accepted 0 moves,
    // but 0 generated leaves the acceptance ratio undefined rather than zero.
    if (isCount(accepted) && isCount(generated) && generated > 0) {
        ratios.acceptedPerGenerated = accepted / generated;
    }
    if (!Object.keys(ratios).length) return null;
    return { ...ratios, generated, accepted };
}
