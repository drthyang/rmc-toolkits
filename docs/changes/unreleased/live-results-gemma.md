### Validation: the first live results (gemma4:26b; Muse Glimmer and Claude Sonnet 5.5 not yet run)

- The live arms ran with one local model, `gemma4:26b` (Ollama), at commit 7456af9 (protocol 2026-10-24.6).
  - Test, counts only: arm B 0.748 and arm C 0.397, against A 0.800, A+ 0.948 and D 1.000.
  - Dev smoke: arm B 0.652, against A 0.814 and A+ 0.952.
  - No run timed out or was cut.
- Three failure patterns:
  - number references written as `{{rN.data.…}}`, which zero otherwise right answers;
  - `element: "all"`;
  - about one tool call per round.
- The dock's streaming, its Live Data divider and the history reset were checked live.
- Records: the dev smoke in `validation/baseline/live/`, the test split sealed in `validation/holdout/live/`,
  and a counts-only summary. The section "Live results (2026-10-25)" in `docs/algorithms/copilot-validation.md`.
