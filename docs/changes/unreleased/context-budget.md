### Copilot: the budget test for a 32k window (plan 5.5)

- `src/__tests__/contextBudget.test.js` runs the worst realistic turn through the real agent loop
  with a scripted model and counts it with the context guard's own estimator and constants: the
  dock's largest system prompt and full tool set, the Demo's context block, the two largest
  skills read whole (`read_skill`'s budget and keep rule), a round of analyses filling the loop's
  12,000-character turn budget, and the 8,192-token reply reserve. The same for the harness's arm
  B on the dev case with the largest first request.
- It does not fit today: the guard counts tool messages at 1.4 characters a token, so the dock's
  worst turn needs 38,525 tokens (5,757 over 32,768) and arm B's 35,081 (2,313 over). The test pins
  every term and the overflow, and keeps "fits in 32,768" as an expected failure (`it.fails`) that
  turns red once it fits. Plan 5.5's other half, moving long skill material to on-demand
  references, waits for the user's skill review; no prompt, tool, skill or snapshot changed.
