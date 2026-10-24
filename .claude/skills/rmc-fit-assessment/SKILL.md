---
name: rmc-fit-assessment
description: >-
  Judge an RMCProfile run before trusting anything it says: what was fitted, whether the χ²
  history has converged, where each dataset's fit misses, which peaks the constraints shape,
  how well the configuration is sampled, signs of overfitting, and what to report. Use for
  "has it converged?", "is the fit good?", "can I trust this run?", and before any structural
  conclusion. NOT for choosing bond windows (coordination-shells) or space groups
  (local-symmetry).
---

# RMC fit assessment

**Status: Draft: awaiting the user's interview.** Written from the Workbench's own record (the
llm README's spot check of 2026-10-10 and the referee-style review that followed it), the run
dashboard's algorithm notes and the tools' checks; step 7 was added in validation batch 1 (the
2026-10-14 tool review, approved by the user that day). A step marked "(to confirm with the
user)" is not yet the user's own practice: say in the answer which choice you made.

RMCProfile moves atoms to fit the data under constraints. It is a configuration optimization, not
a Rietveld refinement: never call it a refinement, and there are no esds of refined parameters.

## Steps

1. **What was fitted.** From the run context: the datasets and their Rwp. A Fourier-transform
   file (`*_FT_*`) is the same measurement as its F(Q), shown in real space: one measurement
   fitted in both spaces counts once. Read `partial_weights` (Faber–Ziman): a partial that
   carries a few % of the signal is weakly constrained by the data (Ga–X in GaTa₄Se₈, about 5 %).
2. **Convergence.** Call `assess_convergence`, or read the context's `watchdog_status` (the
   source of truth). The history is ln(χ²) of ONE fit term, the last `.log` column named in the
   context, not the total. "improving" means still falling over the recent window, even when a
   48-point plot looks flat: quote the change over the last tenth. Say "converged" only when the
   watchdog does.
3. **Where the fit misses.** Call `compare_fit` for each dataset you judge, G(r) **and** F(Q):
   Rw, the share of the misfit per range of r or Q, the largest residual features. A misfit
   concentrated at low r points at local structure or a constraint; at the highest r or the
   lowest Q, at the box size or the scale.
4. **Constraints.** For each pair whose closest-approach constraint (`min_distances_A`) sits
   near its first peak, call `find_gr_peaks`: when g(r) is still high right above the cut, the
   constraint, not the data, shapes the peak's low-r side (the Demo's Ga–Se at 2.2 Å: g at
   8.7 % of the peak). Say which peaks are constraint-limited.
5. **Sampling.** Moves per atom from the model card and the context: about 3.5 accepted moves
   per atom (the Demo) is light sampling, and the distributions of weakly weighted atoms are the
   least reliable.
6. **Overfitting and artefacts** (to confirm with the user): an Rw below the data's noise; a
   detailed claim about a partial the data barely weight; structure at the box scale
   (`analyze_modulations` reports the strongest box-scale waves; in a 10×10×10 box they are not
   physics).
7. **When box-scale waves void a modulation claim.** A long wave (every |q_i| ≤ 0.25 of the
   cell's reciprocal axes) is a correlation RMC leaves in a finite box. It voids any claim of
   order at its own q: never report it as a superstructure, a modulation or a cell doubling. It
   does not void a zone-boundary point, or a mode at a general wavevector, that passes its own
   test (ratio ≥ 3 and ≥ 3σ against independent displacements, corrected over the points
   tested): the test already leaves long waves out of the background it compares with. A long
   wave that dominates every element says the box is small or lightly sampled for statements
   about long-range order (to confirm with the user).

## Judgment calls

- Whether "improving, but slowly" is usable for the question asked, or the run must continue
  (to confirm with the user).
- Which Rw, and which share of the misfit in one r range, counts as a good fit for x-ray and for
  neutron data (to confirm with the user).
- How many independent configurations or runs before quoting a distribution or a small effect
  (to confirm with the user).

## Acceptance checklist

- [ ] The datasets are named correctly (an F(Q) and its FT count once).
- [ ] Convergence comes from the watchdog or `assess_convergence`, with its status and the
      recent change. Never "converged" when it says improving or diverging.
- [ ] `compare_fit` ran on every dataset the answer judges.
- [ ] Constraint-limited peaks are named, with the constraint distance.
- [ ] Sampling (moves per atom) and weakly weighted partials are stated.
- [ ] No order is claimed at a box-scale wavevector (every |q_i| ≤ 0.25).
- [ ] Every number comes from a tool or the context; no invented uncertainty.

## Known failure modes

- **Called a still-falling run converged** ("converged — a plateau", the review's first round):
  the history is ln(χ²) of one term, and a slow decline looks flat. *Firm rule in the app:* an
  answer that says the run has converged while the watchdog or `assess_convergence` says
  improving or diverging gets a warning under it.
- **Read F(Q) and its FT as two datasets** (review): the context now marks the FT file as the
  same measurement.
- **Missed the Ga–Se pile-up at the 2.2 Å constraint** (review): `find_gr_peaks` now checks
  each pair's constraint.
- **Judged one space only** (the 5 K run: G(r) misfit checked, F(Q) not, still "achieved").
- **The verdict was a self-assessment:** models wrote "achieved" even over a cut coordination
  shell. *Firm rule in the app (THY-18):* the Outcome is computed from the checks; the model's
  line only explains it.

## Report

The datasets and their weights; the convergence status and the recent change; per dataset, Rw
and where the misfit sits; the constraint-limited pairs; the sampling; and what the run can and
cannot support.

## References

- [run-numbers](run-numbers.md) — the Demo's and the 5 K run's measured numbers from the 2026-10-10 review, to compare against.

Further reading for Claude Code (not loaded in the app): `docs/algorithms/run-dashboard.md`.
