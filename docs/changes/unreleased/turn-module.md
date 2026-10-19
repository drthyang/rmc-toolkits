### Copilot: one module runs one question; the dock and the harness are its adapters

- `src/llm/copilot/turn.js` is the turn module: `runTurn(spec)` takes the question, the run props,
  a model port, the caller's tools, an ask-user port, the skills, the history, a signal, a clock
  and `{ toolsEnabled, method, limits }`, and builds everything else — the context and the context
  as the model reads it, `ask_user` and `read_skill` after the caller's tools, the prompt variant,
  the history window, the Workbench's questions, the evidence gate, the answer checks' context
  facts and the references. Its events are the loop's, each stamped with `at`, between a `start`
  (what is offered, the prompt, the method, whether tools are active, the shown context) and
  exactly one terminal event: `done`, `stopped` or `failed` (with the message and the hint).
  `describeTurn(spec)` is the same assembly, not run. Invalid input (a tool named `ask_user` or
  `read_skill`, no model) throws at the first step; after `start` it never throws.
- `src/llm/provider/model.js` is the model port (`prepare(intent) → sent`, `stream(sent)`), and
  `providerModel` its production adapter over `streamCopilotChat`: a request asked again after a
  looping reply runs at a temperature of 0.6 or more, and the request as sent never holds the API
  key. The agent loop calls the port and emits each request as sent, with why it was made.
- The dock (`useCopilot.js`) keeps its React state, Stop, the question card, the card focus, the
  models that refused tools and the run reset; CopilotChat no longer builds the clarifier, parses
  the skills or adds `ask_user` and `read_skill`. What the dock sends is unchanged (the
  model-visible snapshots pass as they were). `CHAT_HISTORY_TURNS` is renamed `HISTORY_MESSAGES`.
- The validation harness runs each question through the same module (protocol 2026-10-19.3), with
  four changes, each pinned: the warmer retry at 0.6 or more (it ran at the model's 0.2), the
  no-tools prompt for a model that refuses tools, number provenance against the context as the
  model read it, and every request kept in the record (`record.requests`). `run.js` builds the
  port with `providerModel`, and `spend.js` `meteredModel` meters it (replacing `meteredCall`; a
  request the cap refuses is never sent). Arms A and A+ are unchanged.
- Tests: `turn.test.js` at the interface (a scripted port, a fake ask-user, a fake clock);
  `modelVisible/turnRequests.test.js` (the turn's first request for a dock-shaped and a
  harness-shaped spec equals the snapshot files); the prompt-assembly, memory, warmer-retry and
  "when to ask" tests moved onto the turn; the scripted replays run through `runTurn` with their
  assertions unchanged. Details: [the llm README](../web_app/frontend/src/llm/README.md), *The
  turn*, and [copilot-validation.md](algorithms/copilot-validation.md), Protocol 2026-10-19.3.
