### Tests: an App-level smoke test of the Copilot dock's state lifecycle; the Demo no longer follows a folder's Live Data (plan 1.4)

- `src/__tests__/appDockLifecycle.test.jsx` drives the real App, Dashboard, Detected SG card, run
  module and dock on the Demo run's files, with no model (a recording fake at
  `streamCopilotChat`, where the dock's model port sends each request) and no server (axios refuses
  every call; fetch serves only the Demo's files and a stubbed model list). The Demo's files are
  picked as a folder in the static build with Live Data on; RMCProfile's save is a new `.rmc6f` in
  that folder, seen by the Live Data poll.
- It checks, end to end: after a save the dock draws the configuration divider; the next request
  carries no earlier question or answer, and its system prompt, tools and context are a fresh
  conversation's (the context is the run as App hands it now, the first context with only the new
  move count); memory resumes from the divider; a Detected SG pick (a ladder brick) is kept across
  another save of the same run, on the card and in the Copilot's context; another run (the Demo
  button) brings the card and the context back to `DEFAULT_VIEW.toleranceA` and a new conversation.
- **Fix: opening the Demo stops static-mode Live Data watching the folder picked before.** Before,
  a later save into that folder replaced the Demo's data (name, plots, structure) under the Demo's
  `runKey`, so the Copilot kept the conversation and the Detected SG pick. `handleToggleDemo` now
  clears the folder handle and its signature, and each poll tick reads the folder open at that
  moment and drops a read whose folder was replaced meanwhile. Live Data switched on with no
  folder open now watches the next folder picked (before, it did not start until switched off and
  on). Pinned by two more cases in the same test file.
- Not covered in jsdom (left for the browser check): the real File System Access handle and Web
  Worker (stand-ins here), the provider's streaming, the dock's sheet layout, and Flask-mode Live
  Data with the real dock.
