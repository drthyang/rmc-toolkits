### Tests: a lint for the Copilot's skills (AGENT_PLAN 5.1)

- `src/__tests__/skillLint.test.js` reads every `.claude/skills` file the app bundles (it never
  edits them) and checks three things. **Tools:** every snake_case name in a skill is a tool the
  builders offer on the Demo (all 16), or else a parameter or a result or context field the model
  sees. **Gate numbers:** each threshold a skill states equals the constant of the module that
  decides it (`skillLint/tables.js` `BINDINGS`: 16 bindings to `SIGNIFICANCE_SIGMA`,
  `MODAL_COORDINATION_SHARE`, `CLUSTER_Z`, `MAP_TEST_MIN_PAIRS`, `AXIS_RESOLUTION_SIGMAS`, …, some
  with the statistic they gate and a probe of the code), and every number after a comparator is
  bound or listed as no gate. **Firm rules:** each "*Firm rule in the app:*" is mapped to the
  answer check, refusal or context that enforces it, and that code is run.
- Today's three violations, all in displacement-anisotropy, sit in an allow-list with why: the
  low-index cone has no number (the app uses 15°), "preferred direction" does not name the peak
  test, and the split site is a negative kurtosis where the tools use PC1 ≤ −0.5 (plan 5.3). A new
  violation fails, and so does an entry a fix has cleared.
- `workers/pcaKde.js` exports `AXIS_RESOLUTION_SIGMAS` (no change in behaviour).
