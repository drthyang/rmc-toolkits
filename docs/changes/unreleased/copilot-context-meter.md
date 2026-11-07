### Copilot: a context meter by the chat box; nothing cut silently (THY-76, THY-26)

- "Context 62%" beside the composer: how much of the active model's window the next request takes,
  the reply reserve counted. Neutral below 70 %, amber from 70 % to 90 %, red above. The numbers
  are the context guard's own (`provider/contextMeter.js`): the request counted by the guard's
  estimator (`estimatePromptTokens`), plus the reply room the guard keeps (8,192 tokens; 4,096 on
  LM Studio), against the line the guard refuses at: the model's maximum on Ollama, the context
  LM Studio loaded it with, a Claude model's published window. 100 % is exactly where the guard
  refuses a request. A model whose window is not known (an Ollama model that does not report its
  maximum, LM Studio before it loads the model, a cloud model the meter does not know) shows its
  tokens only, with no warning colour.
- What it shows: while the Copilot works, its latest request; between questions, the larger of the
  next question's first request and the last answer's largest request; after a refusal, the request
  that was refused (over 100 %). The details behind "?": the tokens used of the window, the prompt
  and the reply reserve, where the window comes from, the next and the last request, how accurate
  the count is (deliberately high: 2.5 characters a token, 1.4 for analysis results, about 30 %
  above gemma4's real count), and for a Claude model the cost so far (the same estimate at its
  published price; hidden thinking not included).
- At amber or red: Start fresh (a new conversation) and, when it would shrink the next request,
  Summarize and continue. Summarize re-sends each earlier answer as its sentences with numbers and
  its Outcome line, at most 600 characters; no model writes the summary, so it cannot add or change
  a number. The chat marks where it was done. The details also say what makes room on the server
  (reload the model in LM Studio with a longer context; on Ollama, a model with a longer context).
- Phones: the meter stays folded away until amber.
- Never truncated silently: the chat now says when the oldest questions stop being sent to the
  model (it keeps the last 4), when the last answer goes shortened to its sentences with numbers
  (over 2,000 characters), and when earlier results of an answer were sent shortened to their checks
  and key numbers to fit the answer's budget (a new `compacted` event from the agent loop). What is
  sent did not change; only the saying.
- InfoBadge on phones: a "?" whose popover opens upward (the context meter's, the error ribbon's)
  now shows it at the screen's foot; it used to open above the top of the screen.
- No model-visible text changed: prompts, tool descriptions, tool results and the validation's
  records are byte for byte as before; PROTOCOL_VERSION unchanged. Summarize and continue changes
  what the dock sends only when the user presses it; the validation harness never uses it.
- THY-26 audit: `ask_user` (asks before guessing a parameter that changes the science, 2–4 options
  and Other, waits and resumes, the "ask only when it changes the analysis; never what the data
  answers" rules, replay tests for ask / no ask / cancel) and the context guard (Ollama sized with
  `num_ctx`, refused rather than cut, LM Studio's loaded window checked, a read-share backstop) were
  already on main; the Displacement Directions page no longer scrolls sideways at 768–853 px since
  the device layouts.
- Tests (scripted servers only): 100 % is the guard's refusal line on Ollama and LM Studio; the
  levels' bounds; the Claude windows and prices; the previewed request is the turn's first request;
  the history window's counts; the chat's notes; the compacted event after its tool_end; the dock
  against a scripted LM Studio — the guard's numbers, amber, red, Summarize and continue, a refused
  request, the memory notes.
