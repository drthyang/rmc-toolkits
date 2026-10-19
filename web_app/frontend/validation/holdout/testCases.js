// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// HELD OUT: the validation's test cases (docs/algorithms/copilot-validation.md,
// "The sealed test set").
//
// Builders and agents must not open, print, grep or import this file, or
// anything else in this folder. Only the opt-in loader may import it:
// validation/holdoutLoader.js loadTestCases(), which runs on the --final path
// of run.js (arm-a, live and write with --split test) or when CI sets
// RMC_VALIDATION_HOLDOUT=1. A test that needs these cases goes through the
// loader and is skipped without the flag; split.test.js checks that no other
// module names this folder.
//
// The definitions are pinned by TEST_SET_SHA256 (validation/__tests__/split.test.js):
// a change here changes the test set, and results on the old test set are not
// comparable with the new one. Their ids are public (cases.js TEST_CASE_IDS);
// nothing else of them is written outside this folder, and no tuning ever
// reads them. arm-a-test.json beside this file holds arms A, A+ and D on these
// cases (run.js arm-a --split test --final --out validation/holdout/arm-a-test.json).

import { ACHIRAL_TRAP, b20, CHIRAL_TRAP } from '../cases';

// Sulvanite (P-43m): V at the origin, Cu on the edge centres, S at (x,x,x) and its -43m images.
const sulvanite = (x) => [
    { el: 'V', frac: [0, 0, 0] },
    { el: 'Cu', frac: [0.5, 0, 0] }, { el: 'Cu', frac: [0, 0.5, 0] }, { el: 'Cu', frac: [0, 0, 0.5] },
    ...[[x, x, x], [-x, -x, x], [-x, x, -x], [x, -x, -x]].map((v) => ({ el: 'S', frac: v.map((value) => ((value % 1) + 1) % 1) }))
];

export const TEST_CASES = [
    { id: 'test-offcentre', kind: 'offcentre', params: { elements: ['Sr', 'Ge', 'O'], a: 3.95, supercell: [8, 8, 8], axis: [1, 1, 0], delta: 0.16, seed: 201 } },
    {
        id: 'test-zone-boundary', kind: 'zone_boundary',
        params: { elements: ['Ca', 'Zr', 'O'], a: 4.1, supercell: [8, 8, 8], q: [0.5, 0.5, 0.5], axis: [1, 0, 0], delta: 0.09, decoy: { axis: [0, 0, 1], q: [0, 0.125, 0], amp: 0.12 }, seed: 202 }
    },
    { id: 'test-breathing', kind: 'breathing', params: { shape: 'tetrahedron', element: 'Fe', anion: 'Cl', a: 6.90, edge: 2.75, breathe: 0.035, sigma: 0.05, supercell: [8, 8, 8], seed: 732 } },
    { id: 'test-jahn-teller', kind: 'jahn_teller', params: { elements: ['Rb', 'Mn', 'Cl'], a: 9.5, short: 2.48, long: 2.85, sigma: 0.06, supercell: [6, 6, 6], seed: 204 } },
    { id: 'test-cell-wrap', kind: 'rigid_dimer', params: { element: 'Ti', a: 5.2, axis: 2, bond: 2.29, sigmaCommon: 0.07, sigmaOwn: 0.012, supercell: [8, 8, 8], seed: 205 } },
    { id: 'test-thermal-null', kind: 'thermal_null', params: { elements: ['Pb', 'Ti', 'O'], a: 3.97, supercell: [8, 8, 8], seed: 206 } },
    { id: 'test-symmetry-lowering', kind: 'symmetry_lowering', params: { elements: ['Ca', 'Ti', 'O'], a: 3.85, supercell: [8, 8, 8], axis: [1, 1, 1], delta: 0.06, group: 'R3m', seed: 207 } },
    {
        id: 'test-chirality-trap', kind: 'sites_only',
        params: {
            structure: 'sulvanite Cu3VS4', a: 5.39, supercell: [6, 6, 6], sigma: 0.07, sites: sulvanite(0.235), group: 'P-43m',
            facts: { chiral: false, centrosymmetric: false, polar: false }, traps: [ACHIRAL_TRAP], seed: 208
        }
    },
    {
        id: 'test-chiral', kind: 'sites_only',
        params: {
            structure: 'B20 MnSi', a: 4.558, supercell: [6, 6, 6], sigma: 0.06,
            sites: [...b20(0.138).map((frac) => ({ el: 'Mn', frac })), ...b20(0.845).map((frac) => ({ el: 'Si', frac }))],
            group: 'P2_13', facts: { chiral: true, centrosymmetric: false, polar: false }, traps: [CHIRAL_TRAP], seed: 209
        }
    },
    // The judgment subset's held-out matches (plan item 1.6 ii), one per dev
    // judgment case, written blind from judgment-spec.md alone (protocol
    // 2026-10-19.11). Each seed is the first of its own sequence that its
    // fixed rule accepts, from the case's engine readings only
    // (docs/algorithms/copilot-validation.md).
    {
        id: 'test-two-features', kind: 'jahn_teller', judgment: true, situation: 'two_features',
        params: {
            elements: ['Rb', 'Cr', 'Br'], a: 10.6, short: 2.53, long: 2.97, sigma: 0.055, supercell: [7, 7, 7],
            offCentre: { axis: [0, 1, 0], delta: 0.11 }, seed: 211
        }
    },
    {
        id: 'test-marginal-above', kind: 'zone_boundary', judgment: true, situation: 'marginal',
        params: { elements: ['Ba', 'Hf', 'O'], a: 4.17, supercell: [10, 10, 10], q: [0.5, 0, 0.5], axis: [1, 1, 0], delta: 0.0128, weak: true, decoy: null, seed: 227 }
    },
    {
        id: 'test-marginal-below', kind: 'zone_boundary', judgment: true, situation: 'marginal',
        params: { elements: ['Ba', 'Hf', 'O'], a: 4.17, supercell: [10, 10, 10], q: [0.5, 0, 0.5], axis: [1, 1, 0], delta: 0.0128, weak: true, decoy: null, seed: 223 }
    },
    {
        id: 'test-two-readings', kind: 'offcentre', judgment: true, situation: 'two_readings',
        params: { elements: ['Cs', 'Ge', 'Cl'], a: 5.43, supercell: [8, 8, 8], axis: [1, 1, 1], delta: 0.23, wells: 8, seed: 231 }
    },
    {
        id: 'test-user-choice', kind: 'elongated_clusters', judgment: true, situation: 'user_choice',
        params: { element: 'Re', anion: 'Te', a: 8.4, edge: 2.62, elongation: 0.36, sigma: 0.035, supercell: [7, 7, 7], seed: 241 }
    }
];
