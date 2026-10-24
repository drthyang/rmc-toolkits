### Copilot: the fit-assessment interview's rules in code — convergence on every χ² term, the 5 % constraint, two answer checks (protocol 2026-10-24.1)

- **Every χ² term of the log and their total.** Both runtimes read every χ² column after `Time
  moves_acc moves_gen` (`readChi` / `read_chi_log`, `CHI_TERMS_START`), each row's total and the rows
  per log, and hand them to the Copilot and the watchdog as `plotData.chiTerms` (`/api/plot/data` sends
  them in Flask mode); the per-term histories concatenate only when every log names the same terms. The
  Dashboard's chart still draws the last column.
- **Converged only when no term and not the total still falls beyond its noise** over the last tenth
  of the run (the user's rule, 2026-10-24). A series falls when the least-squares line through its last
  tenth (at least 10 rows, never back past the last restart) drops by more than 3 × the scatter of its
  points about that line — the app's 3σ bar: a flat series of independent noise passes for falling
  0.8 % of the time over 10 rows and never over the Demo's 227 (`src/llm/watchdog/heuristics.js`:
  `seriesTrend`, `judgeConvergence`; the old fixed 0.02 in ln χ² is gone). `assess_convergence` reports
  every term and the total (change, noise, trend), the context's `convergence.last_10pct` replaces
  `change_over_last_10pct` and `recent_slope_per_step`, the system prompt says what `watchdog_status`
  judges, and the watchdog badge moved to the same rule (AGENTS.md next step 6). The Demo stays
  `improving` (F(Q)_1, X_ray_(R)1 and the total fall; its curvature constraints are flat), as does
  dev-fit-unconverged.
- **Constraint-limited above 5 %.** `find_gr_peaks` calls a peak constraint-limited when g(r) just
  above the closest-approach constraint is above 5 % of the peak (`CONSTRAINT_LIMIT_FRACTION`; the
  shell separation keeps 2 %), and advises interpreting nothing on its low-r side and rerunning with
  the constraint loosened (method: arm C reads the reading alone). No reading moves on the Demo
  (Ga–Se 8.7 % at 2.2 Å) or the dev cases (no constraints).
- **Two checks on the answer:** `long-wave` flags order claimed at a box-scale wavevector (every
  |q_i| ≤ the modulation engine's 0.25, or a long wave in words): a long wave says nothing about order
  at its own q; `single-configuration` flags a distribution's numbers (bond angles, coordination
  counts, displacement directions, cluster shapes) quoted without saying they come from one
  configuration. On the answer corpus they raise one finding, on the README GIF's GaSe₄ answer
  (judged right), and no false alarm. Scripted replays on the Demo pin the user's answers 1, 2, 5 and
  6 wherever a check can catch the wrong answer.
- What the model reads: the system prompt +226 characters; `assess_convergence`'s description 274 →
  575 characters and its Demo result 435 → 1,111; the Demo context +486 characters, with the context
  budget raised from 6,000 to 6,500 characters so its second g(r) peaks stay. Arms A, A+ and D on dev
  are unchanged. Approved by the user, 2026-10-24
  (`validation/findings/decisions/2026-10-24-fit-assessment-interview.json`).
