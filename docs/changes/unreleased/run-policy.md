### Copilot validation: the run policy, and LM Studio models at 48k (plan 1.2)

- `validation/runPolicy.js` states what a live run runs under besides the cases, the prompt, the
  tools and the scoring: the server kind and version, the temperature and the warmer retry's
  floor, the loop's limits, the minimum window per server kind, the window read before and after
  the run, the `num_ctx` each Ollama request was sent, the reply cap, the context guard's
  decisions and the final round's `tool_choice` form. Each record keeps it (`record.policy`, and
  `num_ctx` on each request), each transcript prints it under its header, and the report adds
  *What a comparison across models also compares*: the columns' policies side by side, whether
  the limits, temperature and guard are the same in every run, and what no policy sets (sampling
  beyond the temperature, the reasoning mode, the tokenizer, the quantization).
- The decision of 2026-10-19 (the manager, on the user's delegation): the guard, its estimate, the
  8,192-token reply reserve, the budgets and the prompt are unchanged, since refusing rather than
  cutting is wanted. The live runs load each LM Studio model with at least 49,152 tokens (48k),
  which holds the worst two-skill turn of the dock and of arm B with 10,590 and 14,033 tokens to
  spare (`contextBudget.test.js`). Before every run, `run.js live` reads the server
  (`preflight`): another kind, none, a model not loaded, or one loaded below 48k stops the batch
  with what was found, what is needed and how to reload (`lms load <model> --context-length
  49152`). Nothing is sent and no record is written, so the batch resumes. Ollama still sizes
  `num_ctx` itself and is recorded only. `models.js` names each model's server kind.
- A scripted check sends a final round to each provider preset through the production port and
  pins what goes out today: `"tool_choice": "none"` with the tools for OpenAI, Gemini, Anthropic
  and LM Studio, and no `tool_choice` on Ollama's native API, where the tools are still offered.
  No provider behaviour changed; enforcing the final round on Ollama is the model port's job.
- Protocol 2026-10-19.5: what the model is sent is unchanged, but a run's conditions are now part
  of the protocol and earlier records do not state theirs. Details:
  [copilot-validation.md](algorithms/copilot-validation.md), *The run policy*.
