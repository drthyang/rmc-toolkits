### Copilot validation: arm C without the method (plan item 1.8; protocol 2026-10-19.10)

- The turn takes `config.method: 'full' | 'off'` (`src/llm/copilot/method.js`, `TURN_METHODS`), and
  `planTurn` alone decides from it what the model reads. `'full'` (the dock, arm B) sends every
  request byte for byte as before. `'off'` (arm C, decided 2026-10-19) sends none of the method: no
  `read_skill` and no skills' index, no `EVIDENCE_RULES`, no method sentence of `TOOL_RULES` (the
  calibrated-σ reading, which analysis tells a distortion, effect size with z) or `ASK_RULES` ("Settle
  what the data can first: e.g. run find_gr_peaks …"), no method advice in a tool description or a
  check. The tools, their refusals, the Workbench's questions, the evidence gate and the answer checks
  run under both; `start.method` and the harness record say which one a turn followed.
- Tool specs gain an optional `method` field: the advice split out of the description, which keeps
  what the tool does and returns; put back with `'full'` (a string after a space, `{ after, text }`
  inside a sentence), left out with `'off'`. Split out: `list_sites`, `compute_bond_angles`,
  `analyze_clusters` (two pieces), `analyze_modulations`, `detect_symmetry`,
  `measure_symmetry_lowering`, `find_gr_peaks`, `compare_fit`.
- Checks built with `methodCheck` list the advice inside their text; arm C reads the fact and the
  verdict without it, and a refusal still refuses: the cut-window refusal, no angles, the pooled
  directions note, `analyze_clusters`' window failures and split-site note, `find_gr_peaks`' two
  split-shell notes. The text everyone else reads (the step cards, the record, arm A+) is unchanged.
- The harness's `ARMS.C` sets `method: 'off'`; `ARMS.B` is unchanged. What the model sees moves for
  arm C only: `harness/first-request.arm-C.txt` (its system prompt is the `tools-ask.method-off`
  variant, 7,651 characters against 8,575), new `harness/tools.arm-C.json`, four
  `system-prompt/*.method-off.txt` variants and two `tool-results/*.method-off.json` samples; the
  dock's and arm B's snapshots and the budget pins pass unchanged. A test checks arm C's first request
  for every skill name, the skills' step headings, the method sentences, `EVIDENCE_RULES` and every
  piece of tool advice. Arms A, A+ and D are unaffected.
- The evidence gate checks the same rules under method `'off'`, but its send-back names no analysis
  and no tool: "… but this conclusion needs its evidence from an analysis run in this question. …"
  (`OFF_NUDGE`, decided by the manager on the user's delegation, 2026-10-19); arm B's wording is
  unchanged. The gate's notice and finding (the step card, the record) still name the analyses.
  Snapshots: `harness/gate-send-back.arm-B.txt` and `.arm-C.txt`.
