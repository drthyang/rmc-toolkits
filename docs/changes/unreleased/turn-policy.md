### Copilot: one turn policy decides what ends a turn (architecture candidate #10)

- `src/llm/copilot/turnPolicy.js` owns every counter that ends a turn — the rounds of tools, the
  turn's time, the rounds that added nothing, the one reminder and the warmer retry, the
  last-request note, the answer the evidence gate withdrew — and the order of the steps, behind one
  event-driven interface: `createTurnPolicy({ limits, clock, tools, judge }).next(event)`, the
  events `start`, `reply`, `round` and `no_tools`, the steps `ask` (why, tool choice, retry, notes,
  the verification's send-back), `run` and `end` (the `done`, and the text to say first). The agent
  loop reports each event and does what the step says; the policy judges each answer through the
  verification's `judge` and decides when (a send-back only with tools on, before the last request
  and with time left; the `stoppedBy` an answer is judged under; the withdrawn answer stands), the
  verification what. Turn time is read at most once per event, never in round 0, and the user's
  waits stay off it.
- One `done` shape on every path, built only by the policy: `{ type, text, rawText, outcome,
  rounds, stoppedBy, cut, reason, audit, references, gates }`. An answer gains `reason: null` and
  always has `rawText`; the "budget used" answer gains `rawText`, `references: null`, `gates: []` and
  `reason: null`, its Outcome kept to status, note and body (adding `computed` / `reasons` would
  change the harness's `record.final`); an empty reply gains `rawText: ''`, `references: null` and
  `gates: []`.
- `STOPPED_BY` is the one set of endings — answer, max_rounds, no_progress, empty, stopped, error —
  each with `cutShort` and the dock's note, and `stoppedByFailure` maps a thrown turn to 'stopped'
  (an AbortError) or 'error' (a stall included). The turn's `stopped` / `failed` events carry it,
  the dock closes a turn on it, `CopilotView` shows the row's note (its STOP_NOTES are gone, the
  sentences unchanged), the computed Outcome reads `cutShort`, and the validation harness records
  `record.stopped_by` from every terminal event.
- What left the agent loop: its rounds, start time, waits total, no-progress count, stall flag,
  reminder count, warmer flag, last-note flag and `why`, both turn-time tests, the stoppedBy
  derivation and the two inline `done` literals; the notes and fixed answers moved verbatim. The
  calls per round and the questions stay in the loop: they refuse a call and never end a turn.
- Unchanged: what the model sees (every model-visible snapshot), every scripted replay, the
  agentLoop and robustness assertions, the answer corpus, the harness's expectations and its
  `record.final` on every path (checked against the previous code on 23 scripted endings in both
  arms), the turn tests. No protocol bump: what a run sends is unchanged; `record.stopped_by` is a
  new, derived field. Tests: `turnPolicy.test.js` (rows of events → steps, the policy in the loop)
  and `src/__tests__/turnEndings.test.jsx` (every ending through the dock and the harness alike).
  Details: [the llm README](../web_app/frontend/src/llm/README.md), *The turn policy*.
