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

**Status: Interviewed 2026-10-24.** The steps are the user's own answers of that day
(docs/SKILL_INTERVIEWS.md §1).

RMCProfile moves atoms to fit the data under constraints. It is a configuration optimization, not
a Rietveld refinement: never call it a refinement, and there are no esds of refined parameters.

## Steps

1. **What was fitted.** From the run context: the datasets and their Rwp. A Fourier-transform
   file (`*_FT_*`) is the same measurement as its F(Q), shown in real space: one measurement
   fitted in both spaces counts once. Read `partial_weights` (Faber–Ziman): a partial that
   carries a few % of the signal is weakly constrained by the data (Ga–X in GaTa₄Se₈, about 5 %).
2. **Convergence.** Judge every fit term and the total (`assess_convergence`, the context's
   `watchdog_status`). Say "converged" only when none is still falling beyond its noise over the
   last tenth of the run, even when a 48-point plot looks flat; quote the change over the last
   tenth. A run still improving may answer first-shell quantities (distances, CN) with a "still
   improving" caveat; for subtle effects (distortions, directions, order) the run must continue.
3. **Where the fit misses.** Call `compare_fit` for each dataset you judge, G(r) **and** F(Q):
   Rw, the share of the misfit per range of r or Q, the largest residual features. There is no
   fixed Rw cutoff, for X-ray and neutron data alike: distrust the run when most of the misfit
   sits in one r or Q range that the box size or the scale doesn't explain, or when residual
   features rise above the data's noise. A misfit concentrated at low r points at local structure
   or a constraint; at the highest r or the lowest Q, at the box size or the scale.
4. **Constraints.** For each pair whose closest-approach constraint (`min_distances_A`) sits
   near its first peak, call `find_gr_peaks`. When g(r) just above the cut is above 5 % of the
   peak, the peak is constraint-limited (the Demo's Ga–Se at 2.2 Å: g at 8.7 % of the peak):
   interpret nothing on its low-r side, name the pair, and suggest rerunning with the constraint
   loosened.
5. **Sampling.** Before quoting a distribution or a small effect, have at least 3 independent
   runs or configurations and at least 10 accepted moves per atom past convergence; below that,
   the answer says so. Moves per atom are on the model card and in the context: the Demo's 3.5
   accepted moves per atom is too light. The distributions of weakly weighted atoms are the least
   reliable.
6. **Never a finding.** Refuse to report as a finding: an Rw below the data's noise presented as
   a better fit; a detailed claim about a partial under 5 % of the signal; order at a box-scale q
   (step 7); any feature that differs between independent runs.
7. **Long waves:** a long wave is a correlation the finite box can't represent, so it says
   nothing about order at its own q; a point or mode that passes its own test at another q is
   judged on that test alone. A long wave has every |q_i| ≤ 0.25 of the cell's reciprocal axes
   (`analyze_modulations` reports the strongest). The test, as `analyze_modulations` applies it:
   ratio ≥ 3 and ≥ 3σ against independent cells, corrected over the points tested.

## Acceptance checklist

- [ ] The datasets are named correctly (an F(Q) and its FT count once).
- [ ] Convergence judged on every fit term and the total; never "converged" while one is still
      falling beyond its noise over the last tenth.
- [ ] `compare_fit` ran on every dataset the answer judges, read without a fixed Rw cutoff.
- [ ] Constraint-limited peaks are named, with the constraint distance; nothing is read on their
      low-r side.
- [ ] Sampling stated: moves per atom past convergence and the number of independent runs.
- [ ] No order is claimed at a box-scale wavevector (every |q_i| ≤ 0.25).
- [ ] Nothing that step 6 refuses is reported as a finding.
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
- **Order claimed at a box-scale wavevector** (step 7). *Firm rule in the app:* an answer that
  claims order at a box-scale wavevector (every |q_i| ≤ 0.25) gets a warning under it.
- **A distribution quoted from one configuration** (step 5). *Firm rule in the app:* an answer
  that quotes a distribution from a single configuration without saying so gets a warning under
  it.

## Report

Always: the datasets and their weights, the convergence status, each dataset's Rw and where its
misfit sits, the constraint-limited pairs, and the sampling with the number of runs. Never:
"refinement", esds as if parameters were refined, or order at a box-scale q.

## References

- [run-numbers](run-numbers.md) — the Demo's and the 5 K run's measured numbers from the 2026-10-10 review, to compare against.

Further reading for Claude Code (not loaded in the app): `docs/algorithms/run-dashboard.md`.
