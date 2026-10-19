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
    // The marginal pair: re-set at protocol 2026-10-19.20 from the updated spec
    // (its amplitude and seeds only; provenance at the end of this file).
    {
        id: 'test-marginal-above', kind: 'zone_boundary', judgment: true, situation: 'marginal',
        params: { elements: ['Ba', 'Hf', 'O'], a: 4.17, supercell: [10, 10, 10], q: [0.5, 0, 0.5], axis: [1, 1, 0], delta: 0.00853, weak: true, decoy: null, seed: 252 }
    },
    {
        id: 'test-marginal-below', kind: 'zone_boundary', judgment: true, situation: 'marginal',
        params: { elements: ['Ba', 'Hf', 'O'], a: 4.17, supercell: [10, 10, 10], q: [0.5, 0, 0.5], axis: [1, 1, 0], delta: 0.00853, weak: true, decoy: null, seed: 256 }
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

// --- Provenance, moved here at the holdout gate (2026-10-19) ---------------------------
//
// Until the gate, docs/algorithms/copilot-validation.md, docs/CHANGELOG.md, the decision
// record 2026-10-19-test-breathing-recomposed.json and the Arm A tables stated the values
// below. They were scrubbed there and kept here, behind the seal, so nothing is lost. Git
// history keeps the old text too.
//
// Protocol 2026-10-19.2 (copilot-validation.md, as written before the gate):
//
// **Protocol 2026-10-19.2 (test-breathing recomposed).** The held-out breathing case planted Nb₄
// tetrahedra with a Se anion: the skills' reference compound, whose published shells
// `.claude/skills/coordination-shells/lacunar-spinel-shells.md` gives and `read_skill` serves to
// arm
// B. That leaked the case's kind to arm B only. The fix goes in the test case; the skill is the
// user's text and is unchanged.
// - **The new composition:** Fe₄ tetrahedra with a Cl anion, a = 6.90 Å, edge 2.75 Å (Fe–Fe
//   ≈ 2.73–2.78 Å in Fe₄S₄ cubane cores). Shape, breathing (sd 0.035), thermal sd (0.05 Å) and
//   supercell (8³) are unchanged, and a/edge = 2.509 keeps the old shell separation. The shortest
//   Fe–Fe between clusters is 5.32 Å (1.94 × edge). A fresh author chose it, without seeing dev
//   results, dev parameters or arm records. No Fe, Cl, Fe–Cl, Fe₄, iron or chlorine word appears
//   in
//   any skill (main or PR #71's drafts) or in the model-visible snapshots.
// - **The seed, by a rule fixed before looking:** seeds 731, 732, … in order; the first whose
// drawn
//   ε sample sd over the 512 clusters (the generator's own `cellState` draws) lies in
//   [0.033, 0.037]; at most 20 tries. 731 gave 0.03273 (out) and 732 gave 0.03696 (in), so the
//   seed
//   is 732, used by no other case. The rule depends only on the case's own draws, not on any model
//   or
//   arm result.
// - **Difficulty, honestly:** the shorter edge makes the planted radius spread 7 % smaller
//   (0.035 R = 0.059 vs 0.063 Å). The realised spread goes the other way: seed 732's ε lie above
//   the
//   planted sd (0.03696) and the old seed's below (0.03362), so the realised radius spread is
//   0.0622 Å against the old case's 0.0607 Å (2.5 % larger). A+'s cluster engine reads the
//   breathing
//   at z = 30, against 22 before.
// - **Approved by the manager on the user's delegation, 2026-10-19**
//   (`validation/findings/decisions/2026-10-19-test-breathing-recomposed.json`). The user may
//   reverse
//   it: keep their compound in the test case and disclose it as a known arm-B advantage.
// - **Nothing is lost:** no arm B or C run has touched the test set. `TEST_SET_SHA256` changes, so
//   any result on the old test set would not be comparable.
// - **Arms A and A+ on the test set**, re-recorded with `--final` at 4ab9eab: only the three
//   test-breathing records change (readings, the verdict's element, the truth), and no score
//   moves.
//   A scores 0.50 (no page sees cluster breathing). A+ scores 1: Fe breathing found (breathing sd
//   0.107 vs 0.0635 Å for independent atoms, z = 30), no shape distortion (edge spread 0.152 vs
//   0.154 Å). D scores 1. The means stay A 0.81 and A+ 0.99. The dev records are unchanged: a
//   re-run
//   equals `arm-a-dev.json` record for record, apart from protocol and commit.
//
// The decision record's proposed fix, as written before the gate:
//
// Recompose test-breathing, not the skill (the user's text): Fe4 tetrahedra with a Cl anion (a
// 6.90 Å, edge 2.75 Å), a composition chosen by a fresh author without dev results; shape,
// breathe, sigma and supercell unchanged. The seed by a rule fixed before looking: 731, 732, ...
// the first whose drawn epsilon sample sd over the 512 clusters lies in [0.033, 0.037] (732).
// Approved by the manager on the user's delegation, 2026-10-19; the user may reverse it (keep
// their compound and disclose it as a known arm-B advantage). Protocol 2026-10-19.2.
//
// The Arm A tables' test readings, as written before the gate: offcentre, direction within 10°,
// delta 0.162 Å (planted 0.16); jahn_teller, split shell seen, Δ 0.38 Å for 0.37; batch 1's
// symmetry_lowering, the shift named, Ti 0.055 Å along [1 1 1] for 0.06. docs/CHANGELOG.md's
// suite entry: "the 0.2 Å card misses a 0.06–0.07 Å polar shift" (the test case's and the dev
// case's).
//
// Protocol 2026-10-19.11 (copilot-validation.md, as written before the gate):
//
// **Protocol 2026-10-19.11 (the judgment subset's held-out matches, plan item 1.6 ii).** Five
// judgment cases join `TEST_CASES`, one per dev judgment case, each with its situation, kind and
// option: test-two-features (two_features: Rb/Cr/Br Jahn–Teller octahedra, 4 × 2.53 + 2 × 2.97 Å,
// the Cr also ±0.11 Å along [0 1 0], thermal 0.055 Å, a = 10.6 Å, 7³), test-marginal-above and
// test-marginal-below (marginal: one Ba/Hf/O order at q = (½, 0, ½), 0.0128 Å along [1 1 0], 10³;
// the amplitude puts the expected B power at 3σ after the family-wise correction over 21 points,
// E[P] = 15.97), test-two-readings (two_readings: Ge in eight <111> wells 0.23 Å off, Cs/Ge/Cl,
// a = 5.43 Å, 8³) and test-user-choice (user_choice: Re₆ octahedra, Te anion, edge 2.62 Å,
// elongation 0.36 Å, thermal 0.035 Å, a = 8.4 Å, 7³). A fresh author wrote them from
// `validation/judgment-spec.md` alone: no dev judgment parameter, dev record, arm pin, finding,
// transcript or PR description was opened, and no arm or model ran on them. The knobs were
// explored
// on scratch seeds 9001 and up only; each case seed is the first of its own sequence that a rule
// fixed beforehand accepts, from the case's engine readings:
// - **two_features**, seeds 211, 212, …: Cr split by `isSplitSite` and the Cr–Br shell read as two
//   peaks of CN 4 and 2. Seed 211 (first try): PC1 κ = −1.32 (z = −4.99); the shell joined at
//   6.00 over a 0.13 trough, first peak 3.96 at 2.52 Å, the second at 2.96 Å.
// - **marginal**, seeds 221, 222, …: B's σ at q in the default `analyze_modulations` call (every
//   element). The first in [3, 4] is 227 (3.33σ, ratio 20.5, rank 1); the first in [2.5, 3) is 223
//   (2.58σ, ratio 14.8, rank 1). Seeds 221–227 gave 1.52, 1.73, 2.58, 1.80, 2.99, 4.35, 3.33σ.
// - **two_readings**, seeds 231, 232, …: Ge not split and PC1 unresolved, Mardia ≤ −0.3, at the
//   default tiling the map test ≥ 3σ and the peak test < 3σ, the peak ≥ 3σ near a <111> on a finer
//   tiling, the Ge–Cl shell 3 + 3. Seed 231 (first try): Mardia −0.82, orientation tensor −0.24σ;
//   ν = 2 (default) map 6.82σ, peak 0.38σ; ν = 3 peak 4.25σ exactly on <111>; peaks 2.60 and
//   2.85 Å, first peak CN 3.02, shell 6.00.
// - **user_choice**, seeds 241, 242, …: the spec's four constraints. Seed 241 (first try): first
//   minimum 2.75 Å at 0.15 of the first peak; joined at CN 4.00 (peaks 2.62 and 2.89 Å); the
//   average edge 2.79 Å beyond the trough (the first peak alone, to 2.75 Å, holds CN 1.34, a cut
//   by the Workbench's CN margin); the long edges empty from 3.11 Å, the nearest Re–Re between
//   clusters 3.98 Å.
// - **The coincidence check.** The lead compared these values with the dev cases' without showing
//   them, and two coincided: test-two-features' A-site element (now Rb, for Cs: the same spectator
//   role, a 12-coordinate ion of similar size) and test-user-choice's elongation (now 0.36 Å, 25 %
//   below the first choice, 0.48 Å, with the thermal sd from 0.045 to 0.035 Å, so the
//   elongation-to-thermal ratio and the deep trough stay). Both rules were re-run unchanged and
//   took their first seeds again (211, 241); the numbers above are the new configurations'.
//
// `TEST_SET_SHA256` changes: results on the earlier test set are not comparable with this one.

// --- Protocol 2026-10-19.20: the marginal pair re-set ----------------------------------
//
// F1's corrected null (protocol 2026-10-19.16) raised the modulation test's power, and
// judgment-spec.md §2's knob was rewritten for it (protocol 2026-10-19.19): by its formula
// the pair's 0.0128 Å now expects 5.6σ (Z = 45.5). A fresh blind author re-set
// test-marginal-above and test-marginal-below from the updated spec alone, every parameter
// kept but the amplitude and the seeds. No dev judgment parameter, dev record, arm pin,
// finding, doc, transcript or git history was opened, and no arm or model ran on them.
// - **The amplitude:** the spec's expected significance, σ at the median power, is 3σ for
//   n = 1000 cells, σ_B = 0.06 Å and F = 7 × 3 = 21 at δ = 0.0085295 Å (Z = 20.21), written
//   0.00853 Å: expected 3.0003σ (scipy's ncx2.median with the Python engine's
//   zone_boundary_sigma agrees with the JS root); at the mean power it would read 3.14σ.
// - **The engine check:** 1000 configurations on seeds 500001–501000 (outside the case
//   sequence): Hf's σ at q has median 3.06 (sampling error ≈ 0.05), 51.8 % ≥ 3σ (binomial
//   ±3.1 %), 27.6 % in [3, 4], 14.5 % in [2.5, 3), sd 1.36, nullDof 2.990; every σ ≥ 3
//   also cleared the ratio bar. A further 10 000 on seeds 600001–610000: median 2.985 ±
//   0.018, 49.6 % ≥ 3σ.
// - **The seeds, by a rule fixed before any σ was read:** 251, 252, … (at most 100); -above
//   takes the first with σ in [3, 4], -below the first in [2.5, 3) (the default call, every
//   element, the unrounded σ). Seeds 251–256 gave 1.75, 3.17, 1.37, 2.46, 2.46, 2.89σ:
//   -above is 252 (3.17σ, ratio 8.85, rank 1), -below is 256 (2.89σ, ratio 8.42, rank 1).
//   The Python engine reads both σ to 1e-14.
