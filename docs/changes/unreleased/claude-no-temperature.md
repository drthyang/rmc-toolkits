### Copilot: no temperature for Claude in the dock and the watchdog

- Claude Opus 5.5 in the Copilot dock (Anthropic's OpenAI-compatible endpoint,
  `https://api.anthropic.com/v1`) failed every question with "HTTP 400: `temperature` is deprecated
  for this model": the dock always sent the drawer's temperature (0.2 by default), and the watchdog
  always sent 0. Current Claude models take only their own sampling.
- `provider/sampling.js` says which models take no temperature: `takesNoTemperature({ baseUrl,
  model })` holds for a model id that starts with `claude-` or for any model on `api.anthropic.com`.
  `temperatureFor(settings, t)` is null for those and `t` otherwise. The dock's port
  (`useCopilot.js` `copilotModel`) hands `providerModel` `temperature: null` for them, so no request
  carries a temperature key, the warmer retry and the final round included. The watchdog's request
  goes without one too. Every other model is unchanged: Ollama and LM Studio still send 0.2, 0.6 or
  more on the warmer retry, and 0 for the watchdog. The connection test only lists the models
  (`GET /models`), so it sent no temperature before either.
- The settings drawer greys out the temperature field for such a model, with one line: "Claude uses
  its own sampling". The saved value is kept for the next model that takes one.
- Tool choice was checked against the same endpoint: the loop only ever asks for `"auto"` or, on
  the final round, `"none"` (both accepted). It never sends `"required"` or a named function, which
  current Claude models refuse. Nothing changed there; a test pins it.
- No model-visible text changed (prompts, tool descriptions, skills). The validation is untouched:
  its Claude Sonnet 5.5 entry already set `temperature: null`.
- Tests (scripted, fake servers only): the dock's port on Anthropic's API sends no temperature key in
  the first request, the warmer retry or the final round, and its tool choice is only "auto" or
  "none"; LM Studio, a cloud model that is not Claude, and Ollama's native API send 0.2 and then 0.6.
  The dock end to end, through a looping first reply asked again warmer, sends no temperature for
  Claude, and 0.2 then 0.6 for a local model. The watchdog sends no temperature for Claude and 0 for
  Ollama and LM Studio. The drawer's field is off for Claude, with its note, and editable for a
  local model.
