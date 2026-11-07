### Validation: the 0.8.0 live gate (Claude Sonnet 5.5; Muse Glimmer in part)

- The live arms ran at commit 472d9bd (protocol 2026-11-01.1) with the turn's default limits.
  - Test, counts only:
    - Claude Sonnet 5.5, arm B, 42 runs (3 a case): 0.920, against A 0.800, A+ 0.948 and D 1.000.
    - No errors, cuts or verdict errors.
  - Spend: $19.76 by the tracker's ledger, under the $20 cap.
  - Muse Glimmer, arm B: 3 of the 14 test cases (1 run each) in the 30-minute local cap, about 8.5
    minutes a run. Its mean, 0.778, is over a subset and is not compared.
- Compared with gemma4:26b's 0.748 (protocol 2026-10-24.6), with a caveat: protocol 2026-11-01.1 now
  fills the `{{rN.data.…}}` references that zeroed gemma's verdicts, so the two protocols do not
  compare strictly.
- Every Sonnet run's Outcome is *not verified*. None of the answer checks' 104 findings matches one of
  the score's 3 contradictions, so the facts catalogue needs a review.
- Records: the test split sealed in `validation/holdout/live/`, with counts-only summaries in
  `validation/baseline/live/`. The section "Live results (2026-11-07)" in
  `docs/algorithms/copilot-validation.md`.
