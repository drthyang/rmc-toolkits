### Copilot: result statuses on the tool contract; the event stream is the transcript (architecture candidate #11)

- `src/llm/copilot/tools.js` decides how every call ends, and the agent loop only forwards it. A
  call ends as one ToolResult — `{ resultStatus, args, data, checks, error, facts?, ms, cached, note,
  choice? }` — with `resultStatus` one of `RESULT_STATUSES`: `ok`, `refused` (the tool's own `fail`
  check, no error), `error`, `unknown_tool`, `skipped` (the loop's refusals, a closed set of
  reasons: the calls per round, one question at a time, the questions per request). The
  constructors (`executeToolCall`, `skippedCall`, `unknownToolCall`, `repeatedCall`, `withNote`,
  `questionResult`, `toolResult`) check the invariants and throw a TypeError on a result that
  breaks one; `severityOf` is the old `stepStatus`, `ran` the old `executed`; `toolMessage` and
  `compactToolMessage` moved there verbatim.
- `tool_end` carries the result: `resultStatus`, `note` (what the app added for the model; for a
  Workbench question, which sends no message, the text the model read in its place) and `choice`
  (a question the user answered: `{ source, index, confirmed }`, `confirmed` the clarifier choice's
  flag, which replaces reading "Keep …" off the label), beside the severity. No store: the dock and
  the validation harness read the same events, and the verification records the ToolResult
  itself (a `refused` step still counts as an analysis that ran).
- The dock: each step keeps the status, the note and the choice; the question card and the screen
  reader's line read a question's state from them in one wording (a skipped question is announced
  "Skipped: the Copilot uses its default, …", as the card says it, no longer "Answered: …"; a
  question that never reached the user "Not asked: …"); `ToolStep` shows the note sent to the
  model.
- The validation record: each step gains `result_status`, `note`, `choice`, `ref` and
  `analysis`; a step running when a turn stops is listed as stopped; `failed_calls` and
  `unknown_tool_calls` are read from the status (unchanged values), and `refused_calls` is added.
  The transcript prints a non-ok status and the note, and labels a Workbench question's card as not
  sent to the model.
- Protocol 2026-10-19.13: two stricter readings, approved by the manager on the user's
  delegation. "Asked before" in `unsafe_actions` counts only a question that reached the user (not
  a malformed or refused one, nor a window the Workbench settled with no card); number provenance
  credits only what the model was told (each tool message it read, its error and note included;
  never a Workbench question's card). Decision records `2026-10-19-asked-before.json` and
  `2026-10-19-provenance-told.json`.
- Unchanged: every byte, note and order the model reads (every model-visible snapshot, and a new
  corpus of raw tool messages written first, on the code before the change); `isAnalysis`, the
  evidence gate, the Outcome and the computed verdict; arms A and A+; the replays; the answer
  corpus; the scripted records (pinned as files before the change, equal after it with the added
  keys left out). Tests: `tools.test.js` (every status, the invariant table), `toolOutcomes.test.jsx`
  (one turn through the dock and the harness: the drift list), `toolMessages.test.js` (the raw
  messages), the harness's known-answer runs for both readings. Details: [the llm
  README](../web_app/frontend/src/llm/README.md), *How a call ends*.
