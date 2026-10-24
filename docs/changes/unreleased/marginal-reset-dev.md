### Validation: the marginal situation re-set for the independent-cell null (protocol 2026-10-19.19)

- **Why.** The modulation engines' calibrated null (F1, protocol 2026-10-19.16) gives an isotropic
  site three components where the earlier null took one. The dev weak-order pair's amplitude had
  been set from that one-component bound to an expected 3σ, so both of its seeds read far above the
  bar, and "a signal at the edge of significance" was no longer at the edge.
- **The spec** (`validation/judgment-spec.md`, the held-out author's only description of the knob)
  gives the expected significance under the new null:
  - the B power at q is σ²/(3σ² + δ²) · χ'²₃(nδ²/σ²);
  - the null is Σ λ_k χ²₁ with λ = n/(n − 1) · (σ², σ², σ² + δ²)/(3σ² + δ²);
  - the Šidák family is every zone-boundary point of every element;
  - the expected significance is the median realised σ, the σ of the median power.

  The spec also says how to evaluate it and how to check it with the engine on fresh
  configurations. It holds no dev value, and `judgment.test.js` still scans it.
- **The dev pair.** Its amplitude is re-set by that formula to an expected 3σ, and the engine over
  fresh configurations agrees (about half clear the bar). Its seeds come from the pair's own rule,
  fixed before any σ was read at the new amplitude; every other parameter is kept. The values are
  in `validation/findings/decisions/2026-10-19-marginal-reset-dev.json`, which the held-out author
  does not read. `judgment.test.js` checks the formula with the engine's own functions, the pair's
  expected 3σ and its seed rule.
- **Arm C.** `analyze_modulations`' long-wave and in-phase notes lose their causal parentheticals
  ("(common in RMC)", "(as RMC's long waves and acoustic-like thermal motion give)"): what a pattern
  is evidence of is inference, which arm C drops while keeping definitions and limits. Both notes
  are built with `methodCheck`. Arm B and the dock read them byte for byte as before, and every
  model-visible snapshot passes unchanged. New arm C samples:
  `tool-results/analyze_modulations.method-off.json` and `analyze_modulations.element.method-off.json`.
- **Arms A, A+ and D on dev.** Only the two re-set cases' six records move. One A+ record still
  claims the weak order; the other, below the bar, no longer does. They are re-recorded in
  `validation/baseline/arm-a-dev.json`; every other record and every arm A and D score is unchanged.
  The held-out marginal matches are re-set from the spec by their author, separately.
- Approved by the manager on the user's delegation, 2026-10-19; the user may reverse it before the
  freeze.
