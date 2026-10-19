// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The judgment subset's DEV cases (plan item 1.6, part ii;
// docs/algorithms/copilot-validation.md, "The judgment subset"): cases where
// arm A+'s fixed joining rules are ambiguous or wrong, so a model has room to
// add — two features at once, a signal at the edge of significance, two
// tools that read the same atoms differently, a choice only the user can
// make. cases.js appends them to DEV_CASES and builds them with its kinds.
//
// Paired since protocol 2026-10-19.11: each has a held-out match in the
// sealed test set (same situation, kind and option; holdoutLoader.js),
// written by a fresh author who
// never saw these parameters or the dev results: they worked from
// validation/judgment-spec.md alone (the situations, the kinds and the
// difficulty knobs, no values). Keep these parameters, and
// every result on these cases, out of that file and out of any text the
// held-out author reads; the numbers belong in the PR and the baselines.
//
// Seeds. The first three cases start from 112, 113 and 114 and take the first
// seed upward on which arm A+ errs only where the case's situation says, so
// that each case measures its own situation: a side failure would blur it.
// The situation's own failure held on every seed tried; the side failures
// seen (an engine's false positive, a misread group) are reported with the
// PR, not hidden. The octahedral case's seed was taken under the cluster
// window before protocol 2026-10-19.8; since then A+ reads the Workbench's
// own joined shell there, and the case is the subset's control. The two weak-order cases share one amplitude, set from the
// noncentral χ² of the zone-boundary power so that the modulation engine's
// expected significance is 3σ; their seeds come from a rule that takes one on
// each side of the fixed 3σ bar, so it does not depend on which side an arm
// lands: seeds 131, 132, … in order, the first whose realised significance
// at q (every element tested, as A+ and analyze_modulations test it) lies in
// [3, 4], and the first in [2.5, 3).

export const JUDGMENT_DEV_CASES = [
    // Two features at once: Jahn–Teller octahedra whose centre is also
    // off-centred.
    {
        id: 'dev-jt-off-centre', kind: 'jahn_teller', judgment: true, situation: 'two_features',
        params: { elements: ['Cs', 'Ag', 'F'], a: 8.5, short: 2.05, long: 2.4, sigma: 0.05, supercell: [6, 6, 6], offCentre: { axis: [0, 0, 1], delta: 0.15 }, seed: 114 }
    },
    // Two tools read differently: eight-site <111> off-centring.
    {
        id: 'dev-eight-site', kind: 'offcentre', judgment: true, situation: 'two_readings',
        params: { elements: ['Na', 'Nb', 'O'], a: 3.94, supercell: [8, 8, 8], axis: [1, 1, 1], delta: 0.25, wells: 8, seed: 113 }
    },
    // A choice only the user can make: elongated octahedral clusters whose
    // split M–M shell the Workbench joins and asks about.
    {
        id: 'dev-octahedral-clusters', kind: 'elongated_clusters', judgment: true, situation: 'user_choice',
        params: { element: 'W', anion: 'Br', a: 8.5, edge: 2.65, elongation: 0.5, sigma: 0.05, supercell: [8, 8, 8], seed: 115 }
    },
    // A signal at the edge of significance: one weak zone-boundary order, two
    // seeds (the rule above).
    {
        id: 'dev-weak-order-a', kind: 'zone_boundary', judgment: true, situation: 'marginal',
        params: { elements: ['Ag', 'Nb', 'O'], a: 3.95, supercell: [8, 8, 8], q: [0.5, 0.5, 0], axis: [0, 0, 1], delta: 0.018, decoy: null, weak: true, seed: 132 }
    },
    {
        id: 'dev-weak-order-b', kind: 'zone_boundary', judgment: true, situation: 'marginal',
        params: { elements: ['Ag', 'Nb', 'O'], a: 3.95, supercell: [8, 8, 8], q: [0.5, 0.5, 0], axis: [0, 0, 1], delta: 0.018, decoy: null, weak: true, seed: 133 }
    }
];
