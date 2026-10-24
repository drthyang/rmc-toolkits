### Copilot: a refused analysis backs no claim (the user's decision A5; protocol 2026-10-24.5)

- The evidence gate counts only a call that answered. A call that ends `refused` (the tool's own
  `fail` check, such as `compute_bond_angles` refusing a window that cuts the coordination shell)
  backs no claim, under `anyOf`, `firstOffered` and `allOf` alike: a refusal of one option does not
  count, an answer of another does. Until now it ran and backed its claim, so a refused bond-angle
  window still counted its coordination. The user's words: "A5: option 2. A refused analysis
  doesn't count: the gate sends the claim back, and the model re-runs on the whole shell or asks me
  to confirm the window. Land it before the B/C runs." Approved by the user, 2026-10-24
  (`validation/findings/decisions/2026-10-24-refused-not-evidence.json`).
- The send-back says the call refused. For a window cut (a refusal whose data carry `refused`:
  `compute_bond_angles`, `polyhedral_distortion`) it offers both fixes: run it again on the whole
  shell, or, if the user wants that window, ask them (`ask_user`) and call again with
  `window_confirmed: true`. For any other refusal it says only to run the call so that it answers.
  Under method `'off'` (arm C) it names no tool, window or remedy (`OFF_REFUSED_NUDGE`). The gate's
  notice, the end-of-turn finding and the Outcome's reason ("did not run or refused") name the
  refusal too.
- Unchanged: a refused call is still an analysis (the dock's footer, the records' `analysis` and
  `refused_calls`), the Outcome still reads it as a failure, and its facts still reach the checks.
  The system prompt and the tools are unchanged, so every model-visible snapshot and the 48k
  headroom (4,241 tokens for the dock's worst turn, 7,044 for arm B's) stay. Arms A, A+ and D on dev
  are unchanged (no model, no gate).
