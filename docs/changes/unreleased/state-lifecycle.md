### State lifecycle: a new configuration ends the Copilot's memory; the Detected SG pick resets with the run (plan 1.4)

- **A Live Data save that brings a new `.rmc6f` is a new configuration of the same run.** App
  passes the AI Copilot a `configKey` beside its `runKey`: `configEpoch` in Flask mode, and for a
  browser-loaded run (static-mode Live Data included) the signature of its `.rmc6f` files
  (`configurationSignature` in `browserData.js`, which the Flask check now uses too). A save of the
  plots or the `.log` alone changes neither.
- **The conversation stays, but the model no longer reads the old answers as current.** Before, a
  save kept up to four earlier exchanges about the previous configuration in the model's memory
  while the context described the new one. Now each question remembers the configuration it was
  asked under; the dock draws a one-line divider where it changed ("New configuration saved —
  earlier answers were about the previous one", announced to screen readers), after an answer that
  was still running at the save, and a turn's history holds only the questions asked under the
  current configuration since the last change (`turnHistory(turns, configKey)`,
  `configurationBreaks`). The choice keeps the user's conversation readable without misleading the
  model; starting a new conversation on every save would have thrown away what the user can still
  read.
- **What the model is sent changes in one case only:** the first question after a configuration
  change is sent with no memory, exactly as a fresh conversation's first request, and later
  questions remember only the exchanges since the change. The first request of a fresh
  conversation is unchanged (the model-visible snapshots pass as they were). The validation asks
  one question per fresh run, with no Live Data save, so its requests are unchanged: no protocol
  bump.
- **The Detected SG pick resets with the run.** The shared tolerance (`useRunTolerance` in
  `symTolContext.js`) returns to the card's default (`DEFAULT_VIEW.toleranceA`, 0.2 Å) whenever
  another run opens (`runKey`), even one that was open before, so one run's pick never sets the
  next run's card or Copilot context. It is kept across page switches and across a Live Data save or
  a Load of the same run. The analysis pages' Live Data behaviour is unchanged.
- Pinned by `src/__tests__/appStateLifecycle.test.jsx` (Flask and static-mode Live Data: run vs
  configuration, the pick kept or reset), the dock's divider and memory in
  `src/llm/__tests__/copilotChat.test.jsx`, `turnHistory` / `configurationBreaks` in
  `src/llm/__tests__/turn.test.js`, and `configurationSignature` in `browserData.test.js`.
