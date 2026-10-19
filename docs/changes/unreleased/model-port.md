### Copilot: the provider quirks behind the model port (architecture candidate 3)

- `providerModel` (`src/llm/provider/model.js`) owns what differs between providers; the turn's
  and the loop's interfaces are unchanged. Its `prepare` is now async: it detects the server
  (`detectServer`, which probes only localhost or port 11434), puts the tool choice in that
  server's form, and asks the context guard for its decision (`planCopilotChat`, split out of
  `streamCopilotChat` in `contextGuard.js`). It returns the request as it goes out, with
  `plan`: the server kind; Ollama's `numCtx`, `numPredict`, estimate, loaded and maximum context;
  LM Studio's loaded context; and the tool choice's form (`TOOL_CHOICE_FORMS`: `body`, `omitted`,
  `tools_withheld`). `stream` sends under that plan. The guard's numbers and rules are unchanged
  (2.5 and 1.4 characters a token, the 8,192-token reserve, LM Studio's 4,096 tokens of room).
- **Behaviour change for Ollama runs:** the final request (the analysis budget used, `tool_choice:
  'none'`) now carries no tools, because Ollama's native `/api/chat` has no `tool_choice`. The
  model can no longer answer it with calls only, so it no longer ends on the "budget used" fallback
  for that reason. The earlier calls and results still travel as messages. The estimate counts
  what is sent, so the final request's estimate leaves out the tool schemas. OpenAI, Gemini,
  Anthropic and LM Studio get the same body as before: `"tool_choice"` beside the tools.
- **Typed failures.** A request the guard refuses throws `ContextLimitError` from `prepare`, before
  anything is sent. The turn then emits no `request` event for it, as for a spend cap, and a record
  has no entry for it. A model that cannot take tools fails `stream` with a
  `ToolsUnsupportedError`, read from the provider's response: Ollama's 400 `{ error: "<model> does
  not support tools" }`, OpenAI's 400 `{ error: { param: 'tools', code: 'unsupported_parameter' } }`,
  then the old text pattern as the port's one fallback (LM Studio, Gemini, Anthropic). The agent
  loop retries on the type, not on the text: `isToolsUnsupportedError` and its regex left
  `agentLoop.js`. `client.js` and `ollama.js` keep the parsed error body on the Error (`httpError`).
- The loop keeps its final-round guard: calls in the final round are never run, and a reply with no
  answer ends on the fixed "budget used" answer. LM Studio, Gemini and Anthropic take the field but
  are not verified to honour it. The port's server probes now run in `prepare`, which the loop
  puts under the stall clock (`beforeStall`, `robustness.js`).
- Validation (protocol 2026-10-19.7): the run policy's Ollama final-round form is now "no tools
  sent", and the scripted check pins it. Each request's `num_ctx` is read from `sent.plan`
  (`observeContext` is removed). Ollama runs under earlier protocols are not comparable on their
  final round.
- Tests: `modelPort.test.js` (the port over fake servers: the plan and form per server, the refusal
  in `prepare`, the typed failure and every case it must not catch, the loop over the port);
  `planCopilotChat` in `contextGuard.test.js`; the refused request with no `request` event in
  `turn.test.js` and `runPolicy.test.js`. The model-visible snapshots, `arms.test.js` and the
  budget test are unchanged. Details: [the llm README](../web_app/frontend/src/llm/README.md), *The
  model port*, and [copilot-validation.md](algorithms/copilot-validation.md), Protocol 2026-10-19.7.
