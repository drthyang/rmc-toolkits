# AI Copilot (experimental LLM module)

The **AI Copilot** connects the RMCProfile Workbench to an **OpenAI-compatible** LLM: a **local**
server the user runs (Ollama, LM Studio) or a **cloud** provider (OpenAI, Gemini, Anthropic) with an
API key. It replaced the chat-only AI Assistant, which answered from one fixed
summary of the run; the Copilot also **runs the Workbench analyses itself**, reads what came back,
and judges whether that answers the question.

It is **docked beside the pages**, not over them, and it **points at the page** while it explains:

- **The dock.** On a window at least 1080 px wide (and taller than 540 px: laptops, monitors,
  tablets in landscape) the Copilot is a side
  panel in its own grid column: the pages shrink and reflow beside it (their width tiers are
  container queries), so it never covers a card. Drag its left edge (or focus the edge and use
  ←/→, Home/End; double-click or Enter resets) to resize it; the pages always keep at least
  760 px. On tablets in portrait, phones and short windows it is a **bottom sheet** in its own grid row,
  resized from its top edge; when the agent brings a card on screen, the sheet comes down until the
  card fits whole above it (to at most 50 % of the window while a card is in focus, never below
  28 %). Collapsed (the chevron, or Esc inside it) it is a slim **rail** at the
  right edge, or a bar at the bottom, with a red dot when an answer arrived while it was
  collapsed. The open state, the width and the sheet height are remembered
  (`localStorage` `rmc-llm-dock-v1`).
- **Pointing at cards.** The agent's `highlight_cards` tool gives one or more cards a breathing
  blue edge (2.4 s, the same language as MATERIA and NEBULA3D; a steady edge under
  `prefers-reduced-motion`), opens the card's page when needed and scrolls it into view; on a
  plot it can also mark x positions (the g(r) peaks it names), and the plot zooms to them. The
  card of an analysis that is running breathes softly. A highlight ends when the agent points at
  something else, a new question starts, the Copilot collapses, or the user clicks anywhere on
  the page.
- **Card chips.** The answer cites cards as `[label](card:<id>)` (or `[[id]]`); each renders as a
  chip that brings the card back on screen.
- **3D views.** The host's `orient_3d_view` turns a 3D view to look down a crystal direction
  `[u v w]` (or a, b, c, PC1–PC3), centres or rings sites in blue, selects a site as a click would,
  or resets it; `current_view` says what is on screen. So the agent can say "looking down [111],
  the displacements…" and show it.
- A step's **Show in \<page\>** button switches the workbench page and loads the same analysis
  there. Bond Geometry gets the triplet and windows and computes; Displacement Directions gets the
  site and resolution. The page shows its own result beside the conversation.
- **It asks instead of guessing** (`ask_user`). When a question leaves open a choice that changes
  the science — which coordination shell or r-window, which site or pair, which tolerance, which
  run — the agent asks one short question with 2–4 options and its default, then uses exactly the
  option the user picks. Where the run data show that such a choice changes the result, the
  Workbench asks itself, before the model guesses. See *Asking the user* below.
- **It follows skills** (`read_skill`): written methods for the questions it gets, shared with
  Claude Code (`.claude/skills/`), drafts awaiting the user's interview. See *Skills* below.
- **Its answers are checked.** Where an answer contradicts what the analyses or the run context
  established, a **Check** line says so under it and the Outcome chip drops to **Not verified**,
  whatever the model said of its own answer. See *Checks on the answer* below.

How one question runs (one turn of `copilot/turn.js`, which builds the prompt, the tools offered
and the checks from the run props; *The turn* below):

```
        question
           │
   ┌───────▼────────┐   tools offered    ┌──────────────────────────┐
   │  agent loop    │ ─────────────────▶ │ model port (model.js) →  │
   │ (agentLoop.js) │ ◀───────────────── │ contextGuard / client.js │
   └───────┬────────┘  tool calls / text └──────────────────────────┘
           │ validate → run → fit to budget
   ┌───────▼─────────────────────────────┐
   │ tools: run data (runTools.js) and   │   result + deterministic checks
   │ engines from the host (tools prop)  │ ──▶ `tool` message ──▶ next round
   └─────────────────────────────────────┘
           │ … until the model answers, ending with
           ▼   Outcome: achieved | partial | not achieved — <why>
   transcript: text, step cards (ToolStep), Outcome chip
```

The loop has four stages: **plan → act → observe → assess**. The assessment has two halves:

- every tool result carries **deterministic checks** computed from the full engine output. For
  example, "0 angles in these windows — widen them", "no significant peak (1.2σ < 3σ)", or "the
  first g(r) peak sits on the closest-approach constraint".
- the model ends its answer with a verdict line, which the chat shows as a chip.

A model that cannot call tools still works: the Copilot falls back to plain answers from the run
summary. The old Assistant behaved the same way.

The **convergence watchdog** badge on the Dashboard's χ² chart is part of the same module and is
unchanged.

## Privacy

Local providers keep the promise that run data never leaves the device. The browser talks
**directly** to `http://localhost:*`, and the only party that sees the run-derived summaries is the
server the user runs.

**Cloud** providers receive the run context and the **results of the analyses the Copilot runs**.
Both are summaries of the same kind: site statistics, angle modes, significances, never the
configuration file. Cloud use is opt-in, and the settings drawer warns that data leaves the device.
The *Let the Copilot run analyses* switch turns the tools off. Nothing is sent until the user asks a
question. API keys are stored in `localStorage` only.

**The analyses only compute; the view tools only change what is on screen.** No tool moves a
slider, changes an analysis input or the Detected-SG tolerance, or writes or downloads a file.
The view tools (`highlight_cards`, `orient_3d_view`, `current_view`) open a page, scroll, give a
card its edge, mark a plot, turn a 3D camera, ring sites, and select a site as the user's click
would — and nothing else. A step card's **Show in \<page\>** button, which loads an analysis's
inputs on its page, stays the *user's* action.

## How each seam works (reading order)

1. **`context/runContext.js` + `context/pairCorrelations.js`** build the run context.
   `buildRunContext()` distills the dashboard's parsed state into a compact JSON object:
   - cell and angles;
   - the symmetry block: space group, ladder, and Wyckoff orbits ranked by rms displacement;
   - the symmetry noise floor: below it, the ladder's group changes are sampling noise of the site
     means (rms displacement / √copies);
   - pair correlations, with each first shell's mean distance (from r²g, not the g peak) and
     coordination number;
   - move counters;
   - run settings from the `.dat` file, saying when one measurement is fitted in both spaces;
   - per-dataset Rwp, a Fourier-transform file marked as the same measurement;
   - the Faber–Ziman weight of every partial in each fitted radiation (`scatteringWeights.js`);
   - the ln(χ²) history, downsampled to ≤ 48 points with statistics computed on the full series, the
     watchdog's verdict and the change over the last tenth of the run.

   `contextToJson()` serializes it compactly (arrays and small objects on one line), enforces a
   ~6,000-character budget and records each trim with an `*_omitted` count. The context travels with
   every question, so simple questions need no tool.
2. **`prompts/system.js` + `prompts/copilot.js`**
   - `SYSTEM_PROMPT` (`systemPrompt(method)`) holds the domain rules, each part of each ground-rule
     sentence marked a definition (what a context field is), an inference (which field to use, what
     a pattern is evidence of: the method) or conduct (`GROUND_RULES`); with method `'off'` no
     inference part is sent. The rules: what Rwp means, that the history is ln(χ²) of one fit
     term, and "quote the numbers, don't guess".
   - `COPILOT_PROMPT` adds how to work with tools: say why before each call, read the `checks`,
     change the arguments on a `fail`, never repeat an identical call. It also asks for the verdict
     line. `copilotSystemPrompt()` adds `ASK_RULES` when `ask_user` is offered. `TOOL_RULES` and
     `ASK_RULES` are lists of sentences with their method sentences marked
     (`TOOL_RULE_SENTENCES`, `ASK_RULE_SENTENCES`, `METHOD_SENTENCES`): with the turn's method
     `'off'` (*The turn*) they are sent without them, and `EVIDENCE_RULES` and the skills' index
     are not sent.
   - `buildCopilotMessages()` puts the context, inside a fenced JSON block, ahead of the
     conversation. Earlier turns travel as question + final answer only (the newest
     `HISTORY_MESSAGES` = 8 messages: 4 exchanges); the turn module turns its history into them.
   - `prompts/templates.js` holds the watchdog's rigid `STATUS: …` prompt and its parser.
3. **`provider/client.js`** is a hand-rolled OpenAI-compatible client (`/models`,
   `/chat/completions`) with an SSE stream parser over `fetch` + `ReadableStream`.
   - `streamChat()` yields `{ content }`, `{ reasoning }`, `{ toolCall }` fragments and
     `{ finishReason }`.
   - `tools` / `tool_choice` join the request body only when given, so watchdog requests are
     unchanged.
   - Each fragment's provider data (`extra_content`, where Gemini puts its `thought_signature`)
     is kept, because it must be echoed back.
   - An HTTP failure carries its `status` and its parsed error body (`httpError`).
   - **`provider/contextGuard.js` + `provider/ollama.js`**: the context guard decides each request
     before it is sent (`planCopilotChat()`: the server kind, Ollama's `num_ctx` and `num_predict`,
     LM Studio's loaded context, or a refusal), and `streamCopilotChat()` sends it under that plan,
     never a request the model would cut (*The context window* below). Ollama's requests use its
     native `/api/chat` (NDJSON) with a fitting `num_ctx`.
   - **`provider/model.js`**: the model port, the one seam between a turn and whatever answers it
     (*The turn* below). `providerModel(settings, { session })` is the production port. Every
     provider quirk sits behind it: the server kind, the form each server takes `tool_choice` in,
     the guard's plan in the request as sent, and the typed failures (`ContextLimitError`,
     `ToolsUnsupportedError`).
   - `checkConnection()` turns failures into actionable hints: a bare `TypeError` means "server
     down **or** CORS blocked", and 401 / 429 name the key and quota problems.
4. **`copilot/toolCalls.js`** turns streamed fragments into calls.
   - `accumulateToolCall()` keys fragments on `index`, then `id`, then arrival order, and
     tolerates servers that send whole calls or reuse an index.
   - `parseToolArguments()` accepts `''` and fenced JSON.
   - `validateToolArguments()` checks the JSON-Schema subset the tools use: types, `enum`,
     `minimum` / `maximum`, `required`. Unknown keys are dropped with a note. Numeric strings
     are accepted for number fields.
5. **`copilot/tools.js`** holds the tool contract (below) and `executeToolCall()`: parse, validate,
   run with a timeout, and give the tool its own `AbortSignal`, which is aborted on Stop or timeout.
   It decides how every call ends (*How a call ends*, below): one ToolResult with its
   `resultStatus`, from `executeToolCall` (ok, refused, error) or the constructors for a call that
   did not run. A thrown error becomes an error result whose message the model reads.
   `toolMessage()` writes the message the model reads; `fitToBudget()` halves the largest array,
   recording `*_omitted`, until a result fits. It never trims silently.
6. **`copilot/agentLoop.js`**: `runCopilotTurn()` is a pure async generator of UI events, the turn
   module's internal module (`copilot/turn.js` is its one production caller). It calls a model
   port, so the tests script the model's replies, and before each model call it emits the request
   as sent (`request`). It records each tool step with the turn's verification and forwards its
   ToolResult as `tool_end` (the status, the note, the choice beside the severity): it picks which
   constructor answers a call, never the status, and keeps no store — the events are the
   transcript the dock and the harness both read. What ends the
   turn, and why each request is made, is the turn policy's (`copilot/turnPolicy.js`, *The turn
   policy* below): the loop reports each event (the turn began, a reply, a round of calls, the
   model refused its tools) and does what the step says; the policy judges each answer through the
   verification (*Evidence, numbers by reference, the computed verdict*). Its rules:
   - **Budgets.** ≤ 6 rounds of tool calls and no new round after 10 min (the policy's); ≤ 4 calls
     per round (extra calls are answered "skipped") and 120 s per call (the loop's: they refuse a
     call, never end a turn). Each result is fitted to 3,000
     characters, and a turn's tool messages to 12,000 in total. Older results are compacted to
     their checks first, so a small local context window keeps the run context and the rules.
   - **Duplicate calls.** An identical call (same validated arguments) returns the cached result
     with a note.
   - **Final round.** When the budget is used, a last request goes out with `tool_choice: 'none'`
     (the port enforces it on Ollama by sending no tools). Its calls are never run. If it still
     answers with tool calls only (a server that takes the field but does not honour it: LM Studio,
     Gemini and Anthropic are not verified), a fixed "budget used" answer with `Outcome: not
     achieved` closes the turn (the policy's empty-final ending).
   - **No tool support.** The port's `ToolsUnsupportedError` on the turn's first request retries
     the turn once without tools, under the tool-free prompt. The loop reacts to the type only:
     other errors, a schema 400 among them, are shown, whatever their text says.
   - **Preparing.** A request the port refuses in `prepare` (the context guard, a spend cap) is
     thrown before its `request` event. `prepare` may probe the server, so it runs under the stall
     clock (`beforeStall`).
   - **Message rules.** Every listed call gets exactly one `tool` message. Arguments are echoed as
     valid JSON (`'{}'` when unparseable). Whitespace-only assistant content is sent as `null`
     (Anthropic rejects empty text blocks).
7. **`copilot/runTools.js`** holds the tools over data the module already receives:
   - `read_run_data(section)`: one context section as built, before the budget trim;
   - `assess_convergence`: the watchdog heuristics plus tail slopes;
   - `find_gr_peaks(pair)`: partial-g(r) peaks below `r_max` (a display) and the first
     coordination shell as the run reads it (`context/shells.js`, the one shell reading the
     context, the Workbench's questions, the evidence gate, the answer checks and
     `compute_bond_angles` also read): `first_shell` — its window, peaks, mean distance,
     coordination numbers both ways and their sensitivity to the cut, whether it is separated from
     the next shell, and for two or more resolved peaks `first_peak_only` (to the minimum between
     the first two: window, CN, trough; the old `first_peak_window_A`, r ± FWHM, left on
     2026-10-19) — and whether the pair's minimum-distance constraint cuts the distribution. The
     shell does not depend on `r_max`;
   - `compare_fit(dataset)`: where the model misses the data, as Rw and the share of the misfit
     per range of r or Q, and the largest residual features.
8. **`copilot/askUser.js`**: the `ask_user` tool (*Asking the user* below).
   **`copilot/verification/`**: "is this answer verified?" — the evidence gate, the checks on the
   answer, numbers by reference and the computed Outcome behind one interface (*Evidence, numbers
   by reference, the computed verdict* below). Inside it, `outcome.js` `parseOutcome()` splits the
   verdict line from the answer: lenient about Markdown and dashes, strict about the three
   statuses, and it takes the last line.
   **`copilot/turnState.js`** folds loop events into a turn: ordered text segments (with
   reasoning) and tool steps; the footer counts the steps that ended as analyses
   (`tool_end.analysis`). Each step keeps its `resultStatus` (null while running and once
   stopped), `note` and `choice`; `questionReply(step)` reads a question's state from them
   (waiting, stopped, not asked, skipped, answered), the one wording of the card and the screen
   reader's line.
9. **`watchdog/`**: slope heuristics over the ln(χ²) history are the **source of truth** for the
   badge. The LLM only writes the tooltip note, and only when the status flips or enough new data
   arrived.
10. **`useAssistant.js`, `useCopilot.js` and `components/`**
    - `useAssistant()` holds the settings, the connection probe (with auto-connect) and the run
      context.
    - `useCopilot()` is the dock's adapter of the turn module: the conversation's React state,
      Stop, the question card (its ask-user port), the card-focus effects, which models cannot call
      tools (passed to the turn as `toolsEnabled: false`), and a new conversation when `runKey`
      changes. It hands each turn the conversation so far (`turnHistory`); the turn builds the
      memory from it (a question without an answer is dropped). Each question remembers the
      `configKey` it was asked under: after a new configuration the history holds only the
      questions asked under the current one, and `configurationBreaks` says where the dock draws
      its divider. It describes the turn
      (`describeTurn`) for the composer's status line and the tool list.
    - `CopilotChat` is the dock (`components/dockLayout.js` holds its pure size rules): the rail,
      then the panel, and the resize handle. The panel mounts on first open, so nothing is fetched
      before that, and stays mounted while collapsed, so a running answer continues. It composes
      the connection bar, the settings drawer (`ConnectionSettings`, with the analyses switch) and
      `CopilotView`. Focus goes to the composer on open (once connected) and back to the rail on
      collapse. It owns the card-focus controller.
    - `CopilotView` renders the transcript (`ToolStep` cards, `MarkdownMessage` with its card chips,
      the Outcome chip), the composer, and a polite live region that announces step and verdict
      changes.
11. **`focus/focus.js` + `focus/cardRefs.js` + `copilot/viewTools.js`** — pointing at the page.
    - `createFocusController()` sets `data-copilot-focus` (primary | secondary | running |
      preview) on the host's `[data-copilot-card]` elements; `Copilot.css` draws the breathing edge.
      `show(ids, { marks })` opens the main card's page through `onOpenPage(key, null)` when it is
      not on screen (a card under `[aria-hidden="true"]` is off screen), waits up to 3 s for it to
      mount, lowers a tall sheet, and scrolls it into view (not at all when it is already in view;
      smooth unless reduced motion). A MutationObserver re-applies the attribute to a card the host
      re-renders. `clear()` runs the tools' `onClear` undo callbacks (a 3D site ring).
    - `cardRefs.js` reads a model's spelling of a card (`card:` scheme, brackets, case, the title, a
      bare plot file name) and rewrites `[[id]]` as a `card:` link; `MarkdownMessage` lets the
      `card:` scheme through react-markdown's URL filter and the sanitizer and renders known cards as
      chip buttons (unknown ones as their label).
    - `viewTools.js` builds `highlight_cards` from the host's `cards`: the description lists every
      card for the model, and run() resolves the ids, keeps marks for plot cards (`marks: true`),
      and reports what was shown and what was not on screen.

## The turn: one question (`copilot/turn.js`)

One module runs one question. The dock (`useCopilot.js`) and the validation harness
(`validation/harness.js`) are its two adapters: each passes what is its own and gets the events
back. Neither builds a prompt, a tool the turn adds, or a check, so the two cannot drift apart
(before this module the harness lost the warmer retry and the no-tools prompt).

```js
runTurn({
  question, run,           // required; run = CopilotChat's eight run props (the run module's)
  model,                   // required: a model port (below)
  tools = [],              // the caller's tools, in order; never ask_user or read_skill (a TypeError)
  askUser = null,          // the ask-user port; null: no ask_user tool and no Workbench questions
  skills = null,           // the host's bundle { files, references, load } (parsed here) or { list, load }
  history = [],            // [{ question, answer }], oldest first
  signal = null, clock = { now: Date.now },
  toolContext = null,      // more context for a tool's run() (the dock's card-focus controller)
  config = {}              // { toolsEnabled = true, method = 'full', limits = {}, phrasings = null }
}) → AsyncIterable<event>

describeTurn(spec) → { offered: { tools, skills, askUser, view }, prompt: 'tools' | 'view' | 'no_tools',
                       method, toolsActive, shown, skillsError }   // pure: the same assembly, not run
```

- **What it builds.** The run context and the context as the model reads it (`shown`: trimmed to
  its budget); the tools offered (the caller's, then `ask_user` with an ask-user port and tools on,
  then `read_skill` with skills, tools on and method `'full'`), each as the model reads it under the
  method; the system prompt for them; the
  history window (an unanswered question is dropped, each answer compacted to
  `MEMORY_ANSWER_CHARS` = 2,000 characters, the newest `HISTORY_MESSAGES` = 8 messages kept); the
  Workbench's questions (`clarify.js`); and one verification (`createVerification`, *Evidence,
  numbers by reference, the computed verdict*): the run context, its shell catalogue, the tools
  offered, the method, the caller's claim phrasings (`config.phrasings`: the validation harness
  passes the verdict's words) and, when the limits state it, the send-backs per turn
  (`config.limits.maxEvidenceGates`; the loop's own limits do not hold it).
- **The method** (`config.method`, `copilot/method.js`, `TURN_METHODS`): `'full'` (the default:
  the dock, the validation's arm B) sends the skills and the method in the prompt, the tool
  descriptions and the checks, byte for byte as before the switch; `'off'` (the validation's arm C,
  plan item 1.8) sends none of the method. One pure function, `planTurn`, decides from it: no
  `read_skill` and no skills' index; no `EVIDENCE_RULES`, no method sentence of `TOOL_RULES` or
  `ASK_RULES` and no inference part of the domain rules (`prompts/system.js`; the tool-free prompt
  too); each tool offered through `toolAsSent` — its description without its `method`
  advice, and each check of its results without the advice it lists (*The tool contract*). The
  tools, their refusals, the Workbench's questions, the evidence gate and the answer checks run
  under both. `start.method` says which one a turn followed.
- **Tools are active** exactly when `config.toolsEnabled` holds and at least one tool is offered;
  inactive, the request carries no tools and the prompt is the no-tools one. The memory of models
  that refused tools stays with the caller (the dock passes `toolsEnabled: false` for them).
- **Events.** Every event the loop emits (`agentLoop.js`), shape unchanged, each with `at`
  (`clock.now()`), and:
  - `start { offered, prompt, method, toolsActive, shown, skillsError }`: first, exactly once
    (`describeTurn` for this question);
  - `request { round, why, sent }`: before each model call, after the port's `prepare`; `sent` is
    the request as sent (after the notes, the gate's send-backs, compaction, the Workbench's
    appended answer and the retry temperature; from `providerModel`, with the context guard's
    plan); `why` is `first`, `tools`, `reminder`, `gate`, `final` or `no_tools`. A request the
    port refuses in `prepare` (the guard, a spend cap) was never sent and has no `request` event:
    `round` is followed by `failed`;
  - exactly one terminal event, always last, each with its `stoppedBy` (a key of the turn policy's
    `STOPPED_BY`): `done` (the loop's, the policy's one shape), `stopped { error, stoppedBy }` (an
    abort: `'stopped'`) or `failed { error, message, hint, stoppedBy }` (`'error'`, a stall
    included; `explainModelError`: a context-limit refusal or a spend cap lands here too).
- **Order.** `start`, then a Workbench question (a `tool_start` / `tool_end` pair, `app: true`), then
  per request: `round`, `request`, the text and reasoning deltas, the tool pairs (sharing an id). A
  `gate` is always followed by a `retract`. `tools_unsupported` comes at most once, before any tool
  step, and a new `round` and `request` under the no-tools prompt follow it.
- **Errors.** Invalid input (a reserved tool name, no model, a blank question) throws a TypeError at
  the first `next()`; after `start` the iterator never throws. A malformed skill is `skillsError`,
  never a throw: the turn runs without skills.
- **State.** None between calls: the caller passes the history, the clock and whether tools are on.

**The model port** (`provider/model.js`): `{ prepare?(intent, { signal }) → sent, stream(sent,
{ signal }) → fragments }`, a bare function being `stream` with an identity `prepare`; either may
return a promise. The intent is `{ messages, tools, toolChoice, retry }`, `retry` null or
`'warmer'` (the request asked again after a looping reply). Every provider quirk sits behind this
port, none in the turn or the loop.

`providerModel(settings, { session })` is the production port:
- **`prepare`** (async) detects the server (`detectServer`: only `localhost` or port 11434 is
  probed), puts the tool choice in that server's form, and asks the context guard for its plan
  (`planCopilotChat`). It returns the request exactly as it goes out:
  `{ baseUrl, model, temperature, messages, tools, toolChoice, plan }`. A `'warmer'` retry runs at
  `ANTI_REPEAT_TEMPERATURE` = 0.6 or above. The plan is plain JSON:
  - Ollama: `{ server: 'ollama', version, estimate, needed, numCtx, numPredict, loadedContext,
    maxContext }`;
  - LM Studio: `{ server: 'lmstudio', estimate, needed, loadedContext }` (`loadedContext` null:
    not loaded; the request loads it and its window is checked at the reply's first fragment);
  - any other server: `{ server: 'openai' }`.

  With tools, the plan adds `toolChoice`, the form. A request the guard refuses throws
  `ContextLimitError` here, before anything is sent. `sent` never holds the API key, which joins
  the request only as it goes out, so no event or record carries it.
- **`stream`** sends `sent` under its plan, through `streamCopilotChat(request, { plan })`. Its
  first argument keeps the request's own fields, which the model-visible snapshots capture. It
  throws a `ToolsUnsupportedError` when the model cannot take the request's tools. It throws
  `ContextLimitError` only after sending: Ollama read under 0.4 of the prompt, or a model LM
  Studio loaded for the request (JIT) has too small a window, so the reply is stopped unseen.

**Typed failures.** The loop reacts to these by type, never by text:

| Error | Raised by | When | The loop |
| --- | --- | --- | --- |
| `ContextLimitError` (`needed`, `available`) | `prepare`; `stream` after sending | The request would not fit the model's window | Thrown: the turn's `failed` |
| `ToolsUnsupportedError` (`status`, `server`, `cause`; the provider's message) | `stream`, before any fragment, only for a request with tools | The provider's response says the model cannot take tools | First request: retried without tools (`tools_unsupported`); later: thrown |

`isToolsUnsupported(error, { server })` reads the provider's response by its shape first:
- Ollama: HTTP 400 with `{ "error": "<model> does not support tools" }`;
- OpenAI: HTTP 400 with `{ "error": { "param": "tools", "code": "unsupported_parameter" } }`.

The text pattern is the fallback for providers with no known shape (LM Studio, Gemini,
Anthropic), and this is its one place: a 4xx saying tools, tool calling or function calling are
not supported. A schema 400, a 5xx, a failure after the reply began, an abort and a request
without tools are never read as one.

**The tool choice, per provider** (`toolChoiceForm`, `TOOL_CHOICE_FORMS`):

| Server | `'auto'` (every round but the last) | `'none'` (the final round) |
| --- | --- | --- |
| OpenAI, Gemini, Anthropic (OpenAI-compatible), LM Studio | `body`: `"tool_choice": "auto"` beside the tools | `body`: `"tool_choice": "none"` beside the tools |
| Ollama (native `/api/chat`) | `omitted`: no field (the native API has none; its default is auto), the tools sent | `tools_withheld`: no tools sent, so the model cannot call one |

The Ollama final round changed with the model port (protocol 2026-10-19.7). Before it the tools
were still offered there, and a reply of calls only ended on the loop's "budget used" answer.
The earlier calls and results still travel as messages, and the plan's estimate counts what is
sent (no tool schemas).

The validation meters the port (`validation/spend.js` `meteredModel`): a request that could cross
the spend cap is refused in `prepare`, after the guard's plan, so it is never sent.

**What stays in the adapters.** The dock: React state, Stop, the question card (its ask-user port),
the card-focus effects (it maps events to focus), the unsupported-model memory, the `runKey`
reset and the `configKey` history window; it closes a turn on its terminal event's `stoppedBy`
and shows that ending's note. The harness: the case's run (the run module), its compute-only tools, the simulated user,
the arm's skills and method (arm C: `'off'`), the verdict's words for the gate
(`validation/verdict.js` `VERDICT_PHRASINGS`, as `config.phrasings`), its port, and the record (every request as sent, the
steps as their `tool_end`s, the score; number provenance against what the model was told — each
tool message it read and `start.shown`; the safety counts from each step's result status; the
method; `stopped_by` from every terminal event).

## The turn policy: what ends a turn (`copilot/turnPolicy.js`)

One module decides when a turn ends and why each request is made (architecture candidate #10). It
owns every counter that ends a turn — the rounds of tools, the turn's time, the rounds that added
nothing, the one reminder and the warmer retry, the last-request note, the answer the evidence gate
withdrew — and the order of the steps, behind one event-driven interface. The agent loop is its one
caller: the loop streams the model, runs the calls and keeps the messages; it tells the policy what
happened and does what the step says, so it cannot ask in the wrong order, add the last-request
note twice or end a turn twice.

```js
createTurnPolicy({ limits, clock, tools, judge }) → { next(event) → step }
STOPPED_BY               // the one set of endings: { [stoppedBy]: { cutShort, note } }
stoppedByFailure(error)  // a turn that threw: 'stopped' for an AbortError, else 'error'
```

- **Configuration.** `limits` `{ maxRounds, maxTurnMs, maxNoProgressRounds }` (numbers, `Infinity`
  allowed; the loop fills them from `DEFAULT_LIMITS`); `clock.now()` — the turn starts when the
  policy is created, and a Workbench question before the first request is a wait; `tools`, whether
  the first request carries tools; `judge`, the verification's (*Evidence, …*). The policy decides
  when: whether a send-back is allowed, the `stoppedBy` an answer is judged under, which answer
  stands. The verification decides what, and caps its own send-backs (`maxSendBacks`).
- **Events**, in order: `start`, then each request's `reply`, followed by a `round` when it had
  calls, or `no_tools` when the model refused its tools.

  | Event | Fields |
  | --- | --- |
  | `start` | `waitedMs`: the Workbench's question before the first request |
  | `reply` | `text` (none when looped), `calls` (an array), `looped` (stopped for repeating itself), `cut` (at the length limit) |
  | `round` | `progressed` (anything new: a result that was not an error or a repeat, a user's answer), `waitedMs` (the user's answers in it) |
  | `no_tools` | the first request's model cannot take tools (the port's `ToolsUnsupportedError`) |
- **Steps.**

  | Step | The loop |
  | --- | --- |
  | `{ next: 'ask', round, why, toolChoice, retry, notes, sendBack }` | adds the `notes` (`addNote`), then asks: `why` is `first`, `tools`, `reminder`, `gate`, `final` or `no_tools`; `toolChoice` `'auto'`, `'none'` (the request that asks for the answer) or undefined (no tools); `retry` `'warmer'` or null. With a `sendBack` (the verification's), it first withdraws the reply (`gate`, `retract`) and sends the nudge |
  | `{ next: 'run' }` | runs the reply's calls, then reports their `round` |
  | `{ next: 'end', done, say }` | yields `say` as text when there is one, then `done` |
- **The one `done` shape**, built only here: `{ type: 'done', text, rawText, outcome, rounds,
  stoppedBy, cut, reason, audit, references, gates }`. An answer: the verification's fields,
  `reason` null, `rawText` as the model wrote it (the text when no reference changed it). The last
  request came back empty with nothing standing: the fixed "budget used" answer (`say`), an Outcome
  of `not achieved` with its note and body only (no `computed` or `reasons`: the harness copies the
  Outcome whole into `record.final`), `rawText` = text, `references` null, `gates` []. An empty
  reply after its reminder, before the last request: text and `rawText` '', `outcome` null,
  `stoppedBy` `'empty'`, `reason` (why it was empty), `references` null, `gates` [].
- **`STOPPED_BY`**, the endings the dock, the computed Outcome and the harness read:

  | `stoppedBy` | `cutShort` | The dock's note | When |
  | --- | --- | --- | --- |
  | `answer` | no | (none) | the model answered |
  | `max_rounds` | yes | Answered after using the analysis budget for one question. | the round budget or the turn's time was spent |
  | `no_progress` | yes | Answered when repeated analyses stopped adding anything new. | two rounds in a row added nothing |
  | `empty` | no | The model returned an empty reply. (with the reason) | nothing after the reminder, before the last request |
  | `stopped` | no | Stopped. | Stop (an `AbortError`) |
  | `error` | no | Stopped by an error. | any other error, a stall (`ModelTimeoutError`) included |

  `cutShort` makes the computed Outcome partial (`verification/verdict.js`); `CopilotView` shows
  the note; `useCopilot` closes a turn on its terminal event's `stoppedBy` (`stoppedByFailure` for
  its own failure); the harness records `stopped_by` from every terminal event (`record.final`
  is unchanged).
- **Invariants.** `toolChoice` is undefined exactly when tools are off. A step's notes hold 0–2
  messages, in order (the reminder, then the last-request note); the last-request note is added
  once; the request that becomes the last gets why `final`, unless it is the first. One reminder a
  turn; after a looping reply every later request is asked warmer. Turn time = now − start − the
  waits reported, read at most once per event and never in round 0. In a reply without calls the
  order is: the reminder, then the answer the gate withdrew stands, then the empty endings, then
  the judge. The judge gets only non-empty text, and the withdrawn answer with `canSendBack` false.
  `canSendBack` = tools on, this request not the last, and turn time left (and the verification's
  own cap). Exactly one `end` a turn.
- **Errors.** A TypeError for bad configuration, an event out of order or malformed (calls on a
  request that allowed none, `no_tools` after tools ran), a send-back the policy did not allow, or a
  judge that answers otherwise than `{ sendBack }`, `{ done }` or null. The judge's own errors
  propagate (the turn's `failed`); a judge returning null gives the empty ending.
- **What stays in the loop.** The calls per round and the questions (`maxCallsPerRound`,
  `maxQuestions`: they refuse a call and never end a turn), the per-call time, the stall clock and
  the request deadline (they throw: the turn's `failed`, `'error'`), loop detection
  (`repeatsItself`), calls written as text, compaction, the messages, and `DEFAULT_LIMITS`.

## Asking the user (`ask_user`)

The spot check's worst answer was a guess: asked for "Se–Ta–Se angles against an ideal
octahedron", a model picked the first-peak window (2.33–2.71 Å) and cut the TaSe₆ shell in two,
then read the cut as "incomplete octahedra". Which window to use is the user's call when the data
offer two, so the agent asks:

- **The tool** (`copilot/askUser.js`, offered with every run's tools):
  `ask_user({ question, options: [{ label, description? }] (2–4), allow_other = true, default? })`.
  `normalizeAskRequest()` trims the question, reads plain strings as options (and splits
  `"label — description"`), drops repeated labels, and matches `default` to a label or its number;
  it refuses an empty question or fewer than two distinct options with a message the model can act
  on.
- **The card** (`components/QuestionCard.jsx`) shows in the conversation, in the docked panel and
  in the bottom sheet: the question, the options as buttons with their one-line meaning (the
  default marked), an *Other…* field and *Skip*. Keyboard: Tab, ↑/↓, 1–4 to pick, Enter in the
  field. It takes focus when it appears unless the user is typing in a field, and the log never
  scrolls its question out of view. A collapsed dock marks the rail ("a question for you").
  Answered, it keeps the question and the reply.
- **The loop** waits. `useCopilot` hands the turn its ask-user port, `askUser(request, { signal })`,
  which the turn offers `ask_user` with and puts in the tools' context; the
  card settles it, and the answer returns as the tool's result (`answerResult()`: the chosen label
  and its meaning with "use exactly this choice", a typed answer verbatim, or a skip as "use your
  default and say you assumed it"), with the reply beside it as the step's `choice`
  (`replyChoice()`: `{ source: 'option' | 'other' | 'skipped', index }`), which the card and the
  record read; the model never sees it. A question that never reached the user (malformed, nobody
  to ask, or refused by the loop) ends `error` or `skipped`, with no choice: the card reads *Not
  asked*. Stop, a new run or an unmount rejects it with an AbortError: the turn closes as
  *Stopped*, the card reads *Not answered*, and the next question works.
- **The contract.** `interactive: true` marks a tool that waits for the user: never answered from
  the duplicate-call cache, no time limit, its wait not counted toward the turn's 10 min, at most
  `maxQuestions` (2) questions per request and one per step (`agentLoop.js`). It is not a view
  tool and changes no data; the turn footer does not count it as an analysis.
- **The rules** (`ASK_RULES` in `prompts/copilot.js`, added when the tool is offered): ask only
  when the answer changes which analysis runs or its parameters and neither the context nor an
  analysis can settle it; settle what the data can first (e.g. `find_gr_peaks`, then offer its
  shells as the options); never on a clear question or for permission to run an analysis; one
  question at a time, at most two per request; offer a default; then use the chosen option and
  say which.

### With a real model (2026-10-11)

gemma4:26b (Ollama, through the context guard at `num_ctx` 24,576) on the Demo, docked at
1440 × 900, each question in a new conversation, a question card answered with its default:

| Question | Asked? | What it did |
|---|---|---|
| Has the run converged? · first Ta–Se peak, constraint-limited? · most anharmonic site and largest U_iso · Ga preferred direction? · space group at 0.2 and 0.02 Å | no | answered from the context or one analysis, correctly (improving; 2.52 Å; ref 23; no direction at −0.89σ; the P-42₁m rung below the noise floor) · 6–37 s |
| Se–Ta–Se angles against an ideal octahedron · "What are the Se–Ta–Se bond angles?" | no | used the whole first shell, 2.2–3.23 Å (CN 6), from the context or `find_gr_peaks` · 18–28 s |
| "Compute the bond angles around Ta" · coordination number of Ta · Ta–Ta–Ta angles · Ta–Se bond lengths | no | chose the main polyhedron (Se–Ta–Se) and said so, or covered every option (Se, Ta and Ga neighbours; both Ta–Se peaks and the shell mean) · 9–28 s |

It never called `ask_user` in 16 questions, even after the rules named "the bond angles around
Ta" as a case to ask about: where the context settles the choice (the whole shell) that is right,
and elsewhere it stated its assumption. The card and the loop are exercised with a scripted model
(*Tests*). One answer was wrong in the same way twice: it read the raw angle counts (modes 74.5,
88.5, 104.5°; mean 103 ± 29.6°) as "significantly distorted" octahedra, missing the trans angles
the sin θ factor depletes near 180°.

### The Workbench asks itself (2026-10-11)

gemma4 never asked; so where the computed signals show a real ambiguity that changes the result,
the Workbench shows the question card itself, with the options it computed and a default
(`copilot/clarify.js`, `buildClarifier`, over the run's shell catalogue `context/shells.js`). The card is an `ask_user` step marked `app: true`
("The Workbench asks before the analysis"); it counts toward the two questions per request, and
the model's own `ask_user` still works beside it.

- **Before the first request** (`beforeTurn`), on the user's question (`readQuestion`: the
  intent, and the run's elements, pairs and triplets it names):
  - *a bond-angle or coordination request naming only an element:* which neighbours — one
    option per bonded partner (a first shell starting within 1.4× of the nearest one), nearest
    first and the default, each with its window and CN; "Every neighbour type" for coordination.
    The Demo's "bond angles around Ta": Se–Ta–Se (2.21–3.23 Å, CN 6), Ta–Ta–Ta (2.65–3.33 Å, CN
    2.8); not Ga at 3.8 Å;
  - *a named pair or triplet whose first shell is split into peaks, or runs into the next:*
    which window — the whole first shell (default), the first peak only, or a cut 0.2 Å either
    side of the minimum — only when the coordination number changes by ≥ max(0.5, 15 %) between
    them. The Demo's "Se–Ta–Se bond angles": whole shell 2.21–3.23 Å (CN 6) or first peak
    2.21–2.65 Å (CN 3.2; to the minimum between the two Ta–Se peaks since 2026-10-19, before
    r ± FWHM: 2.33–2.71 Å, CN 3.9). The Ta–Ta shell runs into the next (6 % of the peak at the minimum) but
    a cut 0.2 Å either side changes its CN by < 0.5: not asked;
  - *displacement directions of an element with several sites whose site symmetry is not
    cubic:* one site of each orbit (default; its sites are equivalent by symmetry) or all pooled
    — pooling mixes the sites' local frames and can cancel a per-site preference, the review's
    pooled-Ta error. Ga (site symmetry -43m) is not asked about;
  - *the space group with no tolerance named:* which tolerance, when more than one group holds
    clearly above the sampling noise floor (to ≥ 1.25× it). On the Demo only F-43m does: not asked.

  The answer joins the user's question for the model: `(The Workbench asked: "…" — the user
  chose Se–Ta–Se: Se around Ta in the whole first Ta–Se shell, 2.21–3.23 Å (CN 6). Use exactly
  this.)`. The conversation still shows the question as typed. The question's step sends no
  tool message and takes no ref: its `note` is the text the model read in its place (this one, or
  before a call the note that call carries) and its `choice` the user's reply, with `confirmed`
  (whether the answer keeps the user's own window: the clarifier choice's flag, true for "Keep …",
  never read from the label).
- **Before an analysis runs** (`beforeCall`, for tools that declare `shellWindows(args)`:
  `compute_bond_angles`, `analyze_clusters`): a window that overlaps a first shell but holds less
  of it by the same margin is put to the user — the whole first shell (default) or the window as
  asked — and the analysis runs in the one chosen, its tool message saying so. The spot check's
  2.3–2.7 Å Se–Ta–Se window: "it holds 3.8 Se per Ta of the shell's 6". Once per question: the
  choice is remembered by (centre, partner, window), and a later call with that window runs in the
  chosen one unasked (a model that repeats its cut window every round is asked once); a call
  identical to one already answered from the cache asks nothing. A typed "Other" window goes to the
  window the question was about (r12 or r23).
- **Not asked:** a question that settles the choice itself — an ideal polyhedron it names that
  exactly one option's CN matches ("…compare them with an ideal octahedron": only the whole Ta–Se
  shell holds 6), the shell it names ("only the first Ta–Se peak" keeps a partial window as the
  user's, confirmed; "the whole first shell", "the entire shell", "both peaks" run a cut window in
  the whole shell), an element pair or triplet with a clear shell, a site number, a window or a
  tolerance in Å, a window the user chose earlier in the question. The six spot-check questions
  ask nothing.

## Checks on the answer

The reviews found answers that contradicted what the tools had just computed. The tools now
return `facts` beside their data (not sent to the model): `assess_convergence` its status,
`detect_symmetry` its group with the point group's facts, `analyze_displacement_directions`
whether a direction is significant (its peak ≥ 3σ), `analyze_modulations` the ordered points,
`analyze_clusters` what its clusters showed, `find_gr_peaks` and `compute_bond_angles` the
coordination they counted, `measure_symmetry_lowering` its shifts (no check reads them yet); the
run context adds the watchdog's status and the tolerance ladder's groups. When the model answers,
the verification's checks (one row per facts key, the facts catalogue FACTS,
`copilot/verification/facts.js`) read the answer, sentence by sentence and clause by clause, for a
clear positive claim they contradict:

| Check | Fires on |
|---|---|
| convergence | "converged", "reached convergence", "a plateau" while the history is still improving or diverging |
| chirality | a group called chiral whose point group is not a Sohncke group (it has mirrors, glides or rotoinversions) |
| centrosymmetry | a group called centrosymmetric whose point group has no inversion centre |
| inversion lost | an inversion centre "lost" or "broken" when no group it could be lost from has one: neither the groups the answer names nor the parent (the context's average group, or the ladder's top) |
| polar | a group called polar whose point group is not polar |
| direction | a "preferred direction" when no directions analysis of the question found a peak at ≥ 3σ (the whole-map test alone is anisotropy, which thermal ellipsoids show too: review L3) |
| superstructure | a cell doubling or zone-boundary order when no `analyze_modulations` of the question found an ordered zone-boundary point; a superstructure or incommensurate order when none found order at a zone-boundary point or a general wavevector |
| cluster | clusters of an element called distorted, breathing or rigid when every `analyze_clusters` of that element (or of every element, for "the clusters") found otherwise: distorted needs a shape distortion, a pattern (axis family, trigonality), split sites or a distorted average cluster; a polyhedron word naming neither the element nor "cluster" is left alone (review L5) |
| coordination | a coordination number for an element ("four-coordinated", "CN 4", "square-planar coordination") that every whole-shell count of it contradicts by more than half an atom (`find_gr_peaks`' first shell, `compute_bond_angles` in a whole shell), unless a window the user chose gave it; a clause about part of the shell (the first peak, the short bonds, a window in Å) is not a claim about the whole (review L5) |

A negated claim ("has not converged", "achiral", "no preferred direction"), a question, a
hypothetical or conditional mention ("to check whether it has converged", "until χ² plateaus", "I
tested whether site 23 has a preferred direction", "a converged run would…"), or a sentence that
already says the run still moves raises nothing. A negation counts only in the stretch of the
clause just before the claim, after its last comma or "and"/"so": "With no remaining drift in χ²,
the run has converged" is a claim. A chirality claim in a clause that names a group the checks do
not know ("its subgroup P2₁2₁2 is chiral") is not pinned on the sentence's other groups. Symbols
are matched however written
(`P-42₁m`, `P−4 2_1 m`, `$P\bar{4}2_1m$`). Each finding shows as a **Check** line under the
answer, and the computed Outcome is **Not verified** (below). gemma4:26b's correct answers of
2026-10-11 raise no check.

## Evidence, numbers by reference, the computed verdict (THY-18)

Three rules keep an answer tied to what the engines computed in that question. One module holds
them, the verification (`copilot/verification/`, architecture candidate #4): it answers "is this
answer verified?", one verification per turn.

```js
createVerification({ context, catalog, offered, phrasings, maxSendBacks = 2, method = 'full' }) → verification
verification.record(step, show) → { ref, analysis, ...show(ref) }
verification.judge(text, { canSendBack, rounds, stoppedBy, cut }) → { sendBack } | { done } | null
verification.claims        // the claim rows the turn checks: [{ id, what, prompt, analyses, need }]
isAnalysis(step)           // the one definition of an analysis
```

- **Who calls it.** The turn builds it (`turnChecks`: the run context, its shell catalogue, the
  tools offered, tools on or off, the method, `config.phrasings` and `config.limits.maxEvidenceGates`
  when set); the agent loop records every tool step with it, in the order the model reads them, and
  the turn policy judges every answer with its `judge`. `context` gives `{{ctx.…}}` (its shown, trimmed form only), the checks'
  context facts (the watchdog, the ladder's groups, the parents) and the gate's elements; `offered:
  []` means no claim row (no gate), while references are still filled and the checks and the
  Outcome computed (the loop's default without a verification).
- **record(step, show).** `step` is the call's ToolResult (*The tool contract*: `resultStatus`,
  `args`, `error`, `checks`, `facts`, …) with the tool's `name` and kind (`interactive`,
  `reading`, `view`); a step with no known `resultStatus` is a TypeError. The verification reads
  the step's severity (`severityOf`) and whether it ran (`ran`: not an unknown tool or a skipped
  call) from it; a `refused` step ran. `show(ref)` builds the tool message, once and
  synchronously, and the ref then resolves to exactly the data it returns (a callback: the ref is
  written into the message before the data is fitted to the budget). Refs are r1, r2, … in record
  order; a step with an error takes none, and its facts are dropped. `analysis` (`isAnalysis`:
  executed, and not a question, a skill read or a view) rides on `tool_end`; the dock footer
  counts it, as the Outcome does.
- **judge(text, …), in a fixed order.** 1. fill the references; 2. decide a send-back on the
  FILLED text: it needs `canSendBack` (tools can still run: not the final round, the turn's time
  not spent), an unmet claim not sent back yet this turn, in CLAIMS order, and fewer than
  `maxSendBacks` send-backs so far — the raw text becomes the standing answer; 3. run the checks on
  the answer, in FACTS order; 4. add the unfilled-reference finding; 5. add the missing-evidence
  findings, whether tools were on, off or unsupported; 6. compute the Outcome, which the model's
  Outcome line (read from the filled text) can only lower. Blank text judges the standing answer
  (never sent back again), or gives null when there is none; the turn policy never sends blank text
  (it keeps the withdrawn answer and judges it with `canSendBack` false), and the empty endings are
  its own. `sendBack` is `{ id, what, missing, nudge, notice }`; `done` is the turn's `done` event. One
  verification serves one turn: record or judge after `done` throws. A phrasing for an unknown
  claim, or an `asserted` pattern that is not a global RegExp, is a TypeError; model text and tool
  facts never throw.
- **The catalogues are data.** CLAIMS (`claims.js`, one row per claim, in send-back order) and
  FACTS (`facts.js`, one row per facts key, in check order: the producing tools, the shape in the
  dataShape notation, the phrasings it reads, `contradicts` — the `lowering` row has none yet, so
  the missing consumer shows), over one reading of the text (`reading.js`, `phrasings.js`).
  `EVIDENCE_RULES` stays written by hand in `prompts/copilot.js`; `verificationCatalogue.test.js`
  binds it to CLAIMS: `evidenceSentence` renders it from the rows byte for byte, each row's prompt
  names its first analysis, and `promptGaps` pins what the sentence does not say (coordination
  also takes `compute_bond_angles`; a lowering falls back to `detect_symmetry`; "no order" needs a
  per-site map and each finite cluster, and neither `find_gr_peaks` nor
  `measure_symmetry_lowering`). Generating the sentence from the catalogue is the user's call.
- **The evidence gate** (CLAIMS). An answer that concludes something needs, in
  this turn, the analysis that gives the evidence for it: order between cells, or a split site
  called random disorder, needs `analyze_modulations` (ordered antiphase displacements give a site
  the same cloud); a preferred direction `analyze_displacement_directions`; distorted, breathing or
  rigid clusters `analyze_clusters`; a coordination number, a split shell or a Jahn–Teller
  distortion `find_gr_peaks` or `compute_bond_angles` (`anyOf`); a lowering of the average structure
  `measure_symmetry_lowering`, `detect_symmetry` where it is not offered (`firstOffered`); and "no
  order beyond thermal motion" every one of them (`allOf`) — `analyze_modulations`, a per-site
  directions map, and `analyze_clusters` for every element whose own first shell holds at most 4.5
  neighbours (`context/shells.js` `clusterElementsOf`, `CLUSTER_CN_MAX`: a dimer 1, a trimer 2, a
  tetrahedron 3, an octahedron 4 per member). A model cannot be forced to call a tool (not every
  model takes `tool_choice`), so the check runs after the answer: the agent loop withdraws it and
  sends it back once per claim, at most twice a turn (`DEFAULT_MAX_SEND_BACKS`), with what to run (a
  **gate** event; the footer says so), never on the final round or once the turn's time is spent.
  An answer that still lacks its evidence stands, with an `evidence` finding under it (with tools
  off too), and a withdrawn answer is kept if the model then writes nothing. Claims are read as the answer checks read them (`asserts`: not
  negated, hypothetical or asked), and only in a sentence about this run: one that names an
  element of the run (Ti, Mo₃, TaSe₆) or the run itself (this structure, the clusters, site 3) — a
  definition ("A superstructure is order that enlarges the cell") needs no analysis; "no order" is
  always about the run. The gate never asks for a tool the turn does not offer. Under the turn's
  method `'off'` (the validation's arm C, *The turn*) it checks the same claims, but its send-back
  names no analysis and no tool, since which analysis a conclusion needs is the method: "Your answer
  concludes … , but this conclusion needs its evidence from an analysis run in this question. …"
  (`OFF_NUDGE`). The gate event's notice and the end-of-turn finding, which the model never reads,
  still name the missing analyses. The app's own phrasings read prose; a caller adds its own
  through `phrasings` (`{ [claim id]: { asserted, literal } }`): the validation harness adds the
  verdict schema's words (`validation/verdict.js` `VERDICT_PHRASINGS`: `off_centre`, the cluster
  kinds, `jahn_teller`, `symmetry_lowering`, and `"features": []` as "no order"), which the app
  itself no longer reads.
- **Numbers by reference** (filled by the verification; the store is `copilot/references.js`).
  Every tool message carries a `ref` (r1,
  r2, … in the order of the turn); the model writes `{{r2.peak.significance_sigma}}` — the ref, then
  the path in the data it was shown, array items by index (`{{r3.sites.0.shift_A}}`) — or
  `{{ctx.<path>}}` for the run context as the model read it (trimmed), and the app writes the
  engine's value in place where the answer is shown, recorded and checked. A number keeps the tool's rounding, a list of scalars reads
  `[a, b, c]`. A reference that names nothing stays as written and is a `reference` finding. The
  dock fills the earlier text segments of a turn, and a stopped one, the same way (`turnState.js`,
  `{{ctx.…}}` from the turn's `start.shown`).
- **The computed verdict** (`verification/verdict.js`). The Outcome chip is the app's: **not verified**
  when a check contradicts the answer, a reference names nothing or a conclusion's evidence did not
  run; **not achieved** when analyses ran and none succeeded; **partial** when the turn or the
  reply was cut short, an analysis failed and was not rerun successfully, or one warned;
  **achieved** otherwise. The model's `Outcome:` line explains; a grade in it is kept as `claimed`
  and named in the note when the checks disagree ("The model said 'achieved': …"); a lower grade it
  gives stands (it knows when it could not answer). Self-graded, the Outcome had been
  "achieved" nearly always.

**Facts the model reads.** `detect_symmetry` and the context's ladder state each point group's
facts (`context/pointGroups.js`: crystal system, inversion centre, Sohncke class, polarity), so the
model need not recall them (the spot check called P-42₁m chiral). `compute_bond_angles` returns the
10° histogram sin θ-corrected as well (1 for random directions: raw counts deplete 180°; the
engine's bins are shared out by overlap, so 15° or 30° bins leave no comb), and for a whole shell
of about six or four neighbours, in an A–B–A triplet with one window (a mixed or split triplet's
angles are a subset of the polyhedron's), the split against an ideal octahedron (80 % cis at
90°, 20 % trans at 180°) or tetrahedron: the Demo's TaSe₆ shell is 80 % cis at a mean 89° and 20 %
trans at 158.8°, which gemma4 had read from the raw modes as "highly distorted, 75–105°". A cut
window gets no polyhedron reading (its CN 3.9 is not a tetrahedron).

## Skills

**Skills** are Claude Code skills at the repository root, `.claude/skills/<name>/SKILL.md` (front
matter `name` and `description`), so the in-app agent and Claude Code read the same text. Four are
drafted, each opening with **"Draft: awaiting the user's interview"**, with its order of steps,
judgment calls marked "(to confirm with the user)", an acceptance checklist and known failure modes
from the spot check and the review:

| Skill | For | Its firm rules in code |
|---|---|---|
| `rmc-fit-assessment` | convergence, where the fit misses (`compare_fit`), constraint-limited peaks, sampling, overfitting, the report | the convergence check on the answer (*Checks on the answer*) |
| `coordination-shells` | bond windows from the partial g(r), coordination against the ideal polyhedron, reading angle distributions | `compute_bond_angles` **refuses** a window that cuts the shell (fewer than 60 % of the central atoms share the most common count, or most have none; the 60 % is an interview item) and returns the counts that show the cut, not the angles, unless `window_confirmed` — which the Workbench's window question sets when the user keeps their window; the polyhedron reading |
| `local-symmetry` | the tolerance ladder, what breaks between rungs, the noise floor, point-group facts | the point-group facts and the chirality / centrosymmetry / polar / inversion checks |
| `displacement-anisotropy` | ellipsoids, kurtosis, preferred directions and their σ, effect sizes | the preferred-direction check; the Workbench's one-site-or-pooled question |

- **Progressive disclosure** (`copilot/skills.js`). The system prompt lists only each skill's name
  and description (`skillIndex`, 2,150 characters for the four); the agent reads a skill
  (4,700–5,500 characters each) with `read_skill(name)` and a reference file beside it with
  `read_skill(name, reference)`. `read_skill` has its own result budget (`resultBudget` 9,000) and
  is never compacted later in the question (`keep: true` in the tool contract), so the agent can
  follow the skill throughout; its step is not counted as an analysis.
- **The host bundles them** (`src/copilotSkills.js`): the SKILL.md files with the app, the reference
  files as lazy chunks loaded when read; it passes `{ files, references, load }` as the `skills`
  prop, and the module parses them (`parseSkillFiles`: a misnamed skill or an unbundled reference
  is an error, shown in the panel). Vite reads them from the repository root
  (`server.fs.allow`); the Docker build copies them in.
- **Replays** (`src/__tests__/copilotSkillReplays.test.js`): per skill, a scripted model on the
  Demo run with the real engines; the skill's right behaviour passes, and the old wrong one
  (a cut window, "converged — a plateau", "P-42₁m is chiral and centrosymmetry is lost", a
  preferred Ga direction) is refused or flagged.
- **Interviews** turn the drafts into the user's method: `docs/SKILL_INTERVIEWS.md`.

## The context window: never cut silently

A local server runs a model with a fixed context window, and a prompt longer than the window is not
refused: it is shortened to fit, without an error, and the model answers without part of what it
was given. The Copilot's prompt is long. On the Demo the first request is ~28,400 characters (the
rules 8,500, the run context 5,000, fourteen tool schemas 14,800; measured before `read_skill` and
the skills index were added, which bring 15 tools and ~2,150 characters more), which gemma4's
tokenizer reads as ~8,600 tokens (the Ollama log recorded 8,025–8,421 for the same request before
`ask_user`);
each analysis adds up to 3,000 characters more.

- **Ollama's OpenAI-compatible endpoint cannot carry a context length.** It runs the model at the
  server's default: `OLLAMA_CONTEXT_LENGTH`, else one Ollama picks from the GPU memory it finds.
  On the Apple M5 Pro this was written on, 0.40.2 picks 262,144 (51.8 GiB of GPU memory) and the
  models' own maxima then apply (gemma4:26b 262,144, qwen3:32b 40,960), so nothing was cut there:
  Ollama logs `truncating input prompt` when it cuts, and the server logs since 2026-08-26 hold no
  such line from the Copilot. But the same Mac's log of 2026-09-27 shows Ollama starting without
  the GPU (`total_vram="0 B"`) and picking 4,096, half the Copilot's first request; older releases
  defaulted to 2,048–4,096 everywhere.
- **What a cut does (verified 2026-10-11, gemma4:26b, the Demo's first request: 28,382
  characters, 14 tools).** The OpenAI endpoint reports 8,437 prompt tokens at the 262,144 default.
  At `num_ctx` 4,096 the same request still returns 200 with no error, but Ollama first drops whole
  messages (the run context among them) and then cuts what is left to half the window, keeping
  its first 5 tokens (`truncating input prompt limit=2051 prompt=5694 keep=5`): the model reads
  2,051 tokens and, without its rules and the run, calls `current_view` instead of `find_gr_peaks`.
  At 10,240, 16,384 and 24,576 (the guard's choice) it reads all 8,437 and calls `find_gr_peaks`.
  A request whose `num_ctx` equals the context the model is loaded with reuses the runner and its
  prompt cache (no reload; 1 token evaluated).
- **So the Copilot sizes the context itself.** For an Ollama server (`GET /api/version` answers a
  version) it sends its requests to the native `/api/chat`, with `options.num_ctx` = the context
  the model is loaded with when that holds the request and 8,192 tokens of reply
  (`REPLY_RESERVE_TOKENS`), else the next multiple of 8,192 that does (at least 16,384, kept for
  the rest of the conversation, at most the model's maximum from `/api/show`), and
  `options.num_predict` = that context minus the prompt, so a long reply stops at the limit
  instead of pushing the prompt's start out of the window. The size comes from a deliberately
  high estimate (`estimatePromptTokens`: 2.5 characters per token plus 8 per message; ~32 % above
  gemma4's real count on the Demo). Tool results and call arguments, JSON dense with numbers,
  count at 1.4 characters per token (`JSON_CHARS_PER_TOKEN`; gemma4 reads them at 1.64), so the
  estimate errs high there too and `num_predict` never reaches into the prompt.
- **It refuses rather than cuts.** When the request cannot fit even the model's maximum, the turn
  stops before anything is sent, with the sizes and what to do (a new conversation, analyses off,
  a model with a longer context). The decision is made in the model port's `prepare`
  (`planCopilotChat`), so the refused request has no `request` event, and every request that is
  sent carries the plan it was sized by (`sent.plan`). LM Studio fixes the context when a model is loaded: when the
  loaded model's window (`GET /api/v0/models`) cannot hold the prompt and 4,096 tokens of reply,
  the turn stops with the context length to load it with. A model not loaded yet is loaded by the
  request itself (JIT), with LM Studio's default context: its window is read at the reply's first
  fragment, and one too small stops the reply before any of it is shown. Cloud providers and other
  servers are sent the request unchanged (only `localhost` or port 11434 is ever probed); an API key
  set for an Ollama behind a proxy is sent on the native path too.
- **A backstop.** Ollama's `prompt_eval_count` is what the model read: the whole prompt, cached
  or not, and after a cut what was left. Below 0.4 of the estimate (`READ_SHARE_LIMIT`; an uncut
  request reads 0.74 of it, the cut one above 0.18) the turn stops with an error that says the
  reply did not see all of its instructions.
- **A cut reply says so.** A reply that stops at the length limit (`done_reason` /
  `finish_reason` `length`) ends with a warning in the turn's footer instead of passing as a
  whole answer.
- **How large a window.** The worst realistic turn (two skills read whole, a full round of
  analyses, the reply reserve) does not fit 32k: at 32k the guard can refuse it
  (`src/__tests__/contextBudget.test.js`; the fix, on-demand skill references, waits for the
  user's skill review). 48k holds it. The validation's live runs keep the guard as it is and load
  every LM Studio model with at least 48k, refusing a run below it before anything is sent
  (`validation/runPolicy.js`; [copilot-validation.md](../../../../docs/algorithms/copilot-validation.md),
  *The run policy*).

## Robust on local models

A local model can stall, loop or reply with nothing: qwen3.8-27b, a dense reasoning model, took 1–5
minutes a question in the 2026-10-10 spot check and once over 15. `copilot/robustness.js` holds the
pieces the agent loop uses so a turn always ends with an answer or a reason:

| Failure | What the loop does | Limit |
| --- | --- | --- |
| No token for minutes (server busy, stopped) | Cancels the request; the error says the server may be busy | `stallMs` 180 s |
| A reply that never ends | Cancels it; the error suggests a narrower question or a faster model | `requestMaxMs` 600 s |
| Reasoning (or text) repeating a passage | Stops it when caught (`repeatsItself`: its last 200 characters 4 times in the last 8,000), asks again once warmer (`retry: 'warmer'`: the model port sends it at temperature ≥ 0.6, `ANTI_REPEAT_TEMPERATURE`; in the dock and the harness alike) | one reminder |
| Reasoning only, no answer, no call | One reminder, then `done` with `stoppedBy: 'empty'` and a `reason` the footer shows | one reminder |
| A call written as JSON text | Run as a call (`textToolCalls`: the whole reply must be the call, naming a known tool); the text is retracted | |
| Calls that repeat earlier ones or fail | Two such rounds in a row end the analyses: the next request asks for the answer (`stoppedBy: 'no_progress'`) | `maxNoProgressRounds` 2 |
| A long turn | Earlier results are compacted to their checks and top-level numbers (`keyNumbers`), once | `turnCharBudget` |
| A long conversation | Earlier answers above 2,000 characters are re-sent as their sentences with numbers and their Outcome line (`compactAnswer`, in the turn's history window) | `MEMORY_ANSWER_CHARS` (`turn.js`) |
| Out of memory, context overflow, no model loaded, connection lost | The error gets a line on what to do (`explainModelError`) | |

The reminder, the warmer retry, the no-progress rule and the endings are the turn policy's (*The
turn policy*); the loop finds the loop, the calls written as text and whether a round added
anything. Notes to the model (a reminder, the last-round note) join the last user message when that
is the last message: Gemma's chat template refuses two user messages in a row.

## The tool contract (the `tools` prop)

The host passes engine tools as plain objects; the module validates and runs them. `ask_user` and
`read_skill` are the turn's own (it adds them from its ask-user port and the `skills` prop): a host
tool with either name is refused, and the dock says so in its error line.

```js
{
  name: 'compute_bond_angles',           // /^[a-z][a-z0-9_]{0,63}$/, unique
  description: '…',                      // what it computes and returns, for the model
  method: 'Use the whole first …',       // its method advice, apart (optional; below)
  parameters: { type: 'object', properties: {…}, required: [...] },
  dataShape: { triplet: 'string', 'mean_deg?': 'number|null', … },   // the shape of `data` (below)
  label: (args) => 'Bond angles Se–Ta–Se',   // step-card title (optional)
  page: 'geometry',                      // host page key for the step's button (optional)
  pageParams: (args) => ({ … }),         // what that page loads to show the same analysis (optional)
  cards: (args) => ['bond-angles'],      // host cards this call reads: they breathe while it runs (optional)
  view: true,                            // a view tool: changes only the screen; never cached (optional)
  interactive: true,                     // waits for the user (ask_user): never cached, no time limit (optional)
  shellWindows: (args) => [{ centre, partner, window: [lo, hi], apply: ([lo, hi], { confirmed }) => args }],
                                         // its bond windows, checked against the shells before it runs (optional)
  resultBudget: 9000,                    // its own result budget in characters (read_skill) (optional)
  keep: true,                            // its result is never compacted later in the question (optional)
  run: async (args, { signal, focus, askUser }) => ({ data, checks, facts, choice }),
}
```

- `data` is a compact JSON summary built for the model, never a raw engine payload. Put the most
  important rows first: budget trimming halves arrays from the end.
- `checks` is a list of `{ level: 'ok' | 'note' | 'warn' | 'fail', text }`, computed from the
  full result. `fail` means "this does not answer as asked".
- **Method advice** (`copilot/method.js`): text that tells the model which analysis answers which
  question, where to take a parameter from, which window is right or what to run next — what the
  skills teach. The turn's method `'off'` (the validation's arm C) sends none of it, so it is kept
  apart from what a tool does:
  - a spec's `method` holds the advice split out of its `description`, which keeps what the tool
    does and returns. A string is sent after the description and a space with method `'full'`; a
    list of parts puts each back in order, a string after the description and a space, `{ after,
    text }` right after the one place `after` appears in the description (advice inside a
    sentence). With `'full'` the description goes out byte for byte as it was before the advice was
    split out; with `'off'` alone. A part whose `after` is not found exactly once is a TypeError at
    the turn's first `next()`;
  - a check built with `methodCheck(level, ...parts)` (parts: strings and `{ method: '…' }`) has the
    whole text in `text`, as the dock reads it and as the step cards, the record and arm A+ read it,
    and lists its advice in `method`; with `'off'` the turn sends `text` without it (`checkText`). A
    refusal still refuses. In this repository: `compute_bond_angles`' cut-window refusal,
    no-angles failure and octahedron note, the pooled-directions note, `analyze_clusters`' window
    failures and split-site note, `find_gr_peaks`' split-shell notes and unseparated-shell warning,
    `detect_symmetry`'s point-group note.
- `signal` aborts on Stop or timeout. Stop the work if you can.
- `focus` is the dock's card-focus controller, for view tools: `await focus.show(ids, { marks })`,
  `focus.onClear(undo)` (runs when the highlight ends), `focus.state()`.
- `askUser(request, { signal })` shows a question card and resolves with the user's reply
  (`{ answer, source: 'option' | 'other' | 'skipped', index }`), for interactive tools. An
  interactive tool returns the reply beside its result as `choice` (`{ source, index }`): kept on
  the step, never sent (another tool's `choice`, or a malformed one, is dropped).
- `facts` (optional, not sent to the model) are findings the final answer must not contradict
  (*Checks on the answer*), one key per row of the verification's facts catalogue (FACTS, with
  each key's producing tools and shape): `{ convergence }`, `{ groups: [{ symbol, point_group, centrosymmetric,
  sohncke, polar }] }`, `{ direction: { subject, significant, anisotropic, peak_sigma, map_sigma } }`,
  `{ modulation: { elements, ordered, general } }`, `{ cluster: { element, members, distorted, breathing,
  rigid, split, … } }`, `{ coordination: { centre, partner, window, cn, whole_shell, confirmed } }`,
  `{ lowering: { parent, lower, shifted } }`. A key no row names is kept and read by nothing;
  `src/__tests__/verificationFacts.test.js` checks the real tools' facts against the rows.
- The loop gives every result the model reads a `ref` (r1, r2, …) and keeps the data it sent, for
  numbers by reference (*Evidence, numbers by reference, the computed verdict*).
- `shellWindows(args)` names the bond windows a call will use, so the Workbench can ask before a
  window that cuts a coordination shell runs; `apply(window, { confirmed })` returns the
  arguments for the window the user chose (`confirmed`: they kept their own).
- Parameters may be arrays (`items`, `minItems`, `maxItems`) and objects (their own `properties` /
  `required`); a list sent as JSON text or a bare comma list is accepted.

**How a call ends (architecture candidate #11, the tool transcript).** The contract decides it;
the loop only forwards it. Every call ends as one ToolResult:

```js
{ resultStatus, args, data, checks, error, facts?, ms, cached, note,
  choice?: { source: 'option' | 'other' | 'skipped', index, confirmed: true | false | null } }

RESULT_STATUSES                                  // 'ok' | 'refused' | 'error' | 'unknown_tool' | 'skipped', closed
executeToolCall(tool, raw, { signal, timeoutMs, context }) → ToolResult   // ok | refused | error
skippedCall(args, reason, n)                     // SKIP_REASONS: 'call_limit' | 'one_question' | 'question_limit'
unknownToolCall(args, name, names) · repeatedCall(cached) · withNote(result, note)
questionResult(request, reply, { args, ms, confirmed })   // a Workbench question the user answered
toolResult(fields)                               // the checked constructor the others call
severityOf(result) → 'ok' | 'warn' | 'failed'    // the step's status (the dock, the Outcome)
ran(result)                                      // not unknown_tool, not skipped
toolMessage(result, { budget, ref }) → { content, data } · compactToolMessage(content)
```

| `resultStatus` | When | The model reads |
| --- | --- | --- |
| `ok` | the tool answered | the ref, the note, the checks, the data |
| `refused` | the tool answered with a `fail` check (its own refusal: a cut window) | the same; it ran, is an analysis and backs its claim |
| `error` | arguments that are not JSON or that the schema refuses, a throw, the time limit | the error (and the dropped-key notes) |
| `unknown_tool` | no tool has that name | the error, with the names offered |
| `skipped` | the loop did not run it: past the calls per round, a second question in a round, past the questions per request | the error, saying why |

- **Invariants.** A non-empty `error` exactly for `error`, `unknown_tool` and `skipped`, which
  then have no data, no ref and no facts; a `fail` check exactly for `refused` (beside an error a
  check may be anything); `cached` (an identical earlier call repeated, with `REPEAT_NOTE`) only
  for `ok` or `refused`; a `choice` only on `ok`. The constructors throw a TypeError on a result
  that breaks one (a thrown error with no text of its own reads "the tool failed without saying
  why"). `note` is what the app adds for the model beside the result — the Workbench's (its
  choice applied or asked, before a repeat's) — and for a Workbench question, which sends no
  message of its own, the text the model read in its place. `choice.confirmed` says whether a
  window question's answer keeps the user's own window (the clarifier choice's flag); null for
  any other question.
- **The message.** `toolMessage(result)` is the only writer of the bytes the model reads: the
  ref (with no error), the error, the note, the checks, then the data fitted to the budget less
  the rest, in that key order. The status and the choice are never sent.
  `compactToolMessage` shortens an earlier one to its ref, checks and headline numbers by reading
  its own JSON back. `src/__tests__/modelVisible/toolMessages.test.js` pins every message no other
  snapshot covers (a call that did not run, the questions, unreadable arguments, a throw, the time
  limit, a repeat, the Workbench's notes, the budget, the compaction, `highlight_cards`) as raw
  strings.
- **The event stream is the transcript.** `tool_end` carries the result as it ended: `status`
  (`severityOf`), `resultStatus`, `args`, `data` (as the model read it), `checks`, `error`, `note`,
  `choice`, `ms`, `cached`, `ref`, `analysis`, `pageParams`. Nothing re-derives a status from an
  error's text or a data field: the dock (`turnState.js`, `QuestionCard`, the screen reader's line,
  `ToolStep`, which shows the note sent to the model) and the validation harness (each record
  step's `result_status`, `note`, `choice`, `ref`, `analysis`; its safety counts) read the same
  events, and the verification records the ToolResult itself.

**The data's shape (`dataShape`).** With numbers by reference the model writes
`{{r3.sites.0.shift_A}}`: the paths of a tool's `data` are an interface it writes against, and a key
renamed or dropped breaks every reference to it. So a tool that returns data declares the shape of
that data on its own spec (`copilot/dataShape.js`, plain data that reads like the JSON):

```js
dataShape: {
  note: 'string',                        // a JSON type: number (finite), string, boolean, null
  'element?': 'string',                  // ? : the key may be absent
  sites: [{ ref: 'number', shift_A: 'number|null', direction_uvw: oneOf(['number', 'number', 'number'], 'null') }],
  cn: { '*_per_*': 'number' },           // a key with * : any keys whose names fit it
  window_A: ['number', 'number']         // [shape] a list; [shape, shape, …] a tuple
}
```

It is internal: `toRequestTools` sends a tool's name, description and parameters only. It describes
`data` as the tool returns it, before the result budget (which halves the largest list and records
`<key>_omitted` beside it). `checkShape(value, shape)` lists every mismatch by its reference path;
`src/__tests__/toolDataShapes.test.js` runs every tool on the Demo and on each dev case of the
validation and checks every result against its shape, both ways: no key the shape does not declare,
and every declared key returned by some call (the forms those runs cannot reach are listed there).
**A PR that changes what a tool returns changes its `dataShape` in the same diff**, so the change
shows; the engine tools' shapes sit beside the summarizers that build their data
(`src/copilotTools.js`), the run tools' in their specs (`copilot/runTools.js`, with
`read_run_data`'s sections). The view tools (`highlight_cards`, `orient_3d_view`, `current_view`)
describe the screen and declare none yet.

`pages` (`{ key: label }`) and `onOpenPage(key, params)` let a step offer **Show in \<page\>**.

- `params` is the tool's `pageParams(args)`, computed from the validated arguments. Without
  `pageParams`, the button reads **Open \<page\>** and `params` is null.
- `runKey` names the open run; when it changes, the conversation restarts.
- `configKey` (optional, null by default) names the configuration the run shows. When it changes
  and `runKey` does not (a Live Data save brought a new `.rmc6f`), the conversation stays on
  screen, a divider ends the earlier answers ("New configuration saved — earlier answers were about
  the previous one"; after an answer still running, which was asked about the previous one), and
  the next question is sent with only the questions asked under the current configuration as
  memory: the first one after the change is sent as in a fresh conversation. A host that passes
  none keeps one configuration per run.
- `onOpenChange(open)` reports the dock opening and collapsing.
- `cards` lists what the agent may point at: `{ id, page, title, about, kind?, file?, marks?, view3d?, sites? }`.
  Each card element carries `data-copilot-card="<id>"`. `kind` / `file` tie a card to the data the
  run tools read (a plot kind, `r_value`, `symmetry`, `model`), so `assess_convergence`,
  `find_gr_peaks`, `compare_fit` and `read_run_data` breathe the right card while they run.
- `onMarks({ [cardId]: [{ x, label }] })` hands the host the marks to draw inside a plot card.
- `onLayoutChange({ mode, open, width, sheet })` reports where the dock sits; `className` and
  `storageKey` are the host's.

In this repository the host tools live in `src/copilotTools.js`:

| Tool | What it runs |
|---|---|
| `list_sites` | the per-site PCA table, through the worker or `/api/pca/sites`: per-axis kurtosis on resolved axes, PC1 as `[u v w]`, and a split site's static displacement with its sampling error (`delta_A`, `delta_se_A`) |
| `analyze_displacement_directions` | the orientation engine with calibrated significances and `[u v w]` (Show in Displacement Directions selects the site and resolution) |
| `compute_bond_angles` | the triplets engine (Show in Bond Geometry fills the triplet and windows, then computes); the sin θ-corrected histogram, and a whole shell of ~6 or ~4 neighbours read against an ideal octahedron or tetrahedron |
| `analyze_clusters` | the cluster engine (`workers/clusters.js`, `/api/clusters`): the shape of an element's finite clusters (M₄ tetrahedra, dimers, trimers) against a shuffle null of independent atoms, and the average cluster's own distortion. The engine returns its own verdict (distorted, breathing, rigid, split member sites, …); the tool shapes the data and says each finding (`CLUSTER_SENTENCES`), and its `facts` are the verdict's answers |
| `analyze_modulations` | the modulation engine (`workers/modulations.js`, `/api/modulations`): each element's displacement power at the zone-boundary points and its strongest modes against the configuration's diffuse background — order that enlarges the cell, which the folded average cannot show, and correlation between cells at a general wavevector (incommensurate or short-range). The engine returns each point's calibrated σ (corrected over every element's points, whichever the call names), each point's polarisation (the direction the atoms move in that mode, with its share of the power) and each element's verdict (cell doubling, general order, long wave); the tool renders them (`MODULATION_SENTENCES`), the polarisation of ordered points only, and never prints a ratio or σ below its bar rounded onto it |
| `detect_symmetry` | the symmetry finder at any tolerance, without changing the card, with the point group's facts (inversion, Sohncke, polar); offered when the run's structure has an average-structure basis (a browser-loaded run, or a server folder through `/api/structure`, which above the finder's 2000-site cap sends only its size; not a coords-only file). Above the cap it answers "not analysed", in both runtimes |
| `measure_symmetry_lowering` | the size and direction of the coherent shift of a symmetry lowering: each average site's distance from the position a parent group puts it at, against the noise floor; offered as `detect_symmetry` is |
| `orient_3d_view` (view) | turns a 3D view (`src/copilotViewTools.js`): look down `[u v w]` / a, b, c / PC1–3 (c up, or b looking down c), centre or ring sites, select a site, reset |
| `current_view` (view) | the page on screen, its cards, the highlight, and each 3D view's direction and site |

The 3D views join a registry while mounted (`src/viewer3d.js`, `src/useCopilotViewer.js`,
`src/cameraPose.js`), under their card ids: `density-cell` (Atomic Density), `bond-cell` (Bond
Geometry), `pca-density` and `pca-sites` (PCA Ellipsoid), `direction-sphere` and `direction-sites`
(Displacement Directions). A turn is a 650 ms slerp of the camera (at once under reduced motion)
that lands exactly on the requested pose; in the sphere's PCA frame a crystal direction is
projected onto the site's principal axes. The host's card catalog is `src/copilotCards.js`.

A browser-loaded run is analysed in a worker of the Copilot's own, which Stop or a timeout
terminates. A server folder goes to Flask with the request cancelled on abort.

## Provider setup

- **Ollama** (local): <https://ollama.com>. Pull a model that can call tools
  (`ollama pull qwen3`, `llama3.1`, …). Start the server with this site allowed as an origin:

  ```bash
  OLLAMA_ORIGINS="https://drthyang.github.io" ollama serve
  ```

  The Copilot sets the context length of its own requests through Ollama's native API (*The
  context window* above), so `OLLAMA_CONTEXT_LENGTH` is not needed for it.
- **LM Studio** (local): <https://lmstudio.ai>. Load a tool-capable model with a large context, 48k
  or more for a turn that reads two skills (the Copilot refuses a turn the loaded context cannot
  hold, and says how large to load it), enable
  **CORS** in the Developer/server settings, and start the local server
  (default `http://localhost:1234/v1`).
- **OpenAI / Gemini / Anthropic** (cloud): pick the provider and paste an API key. There is no server
  to run, but your run context and analysis results are sent to that provider. Some providers or
  networks block direct browser calls (CORS).

In the app: click the **AI Copilot** rail at the right edge (a bar at the bottom on portrait tablets and phones), open the settings (gear),
pick a provider, press **Test**, then choose a model. Settings
persist in `localStorage` (`rmc-llm-settings-v1`; the analyses switch is `copilotTools`).

HTTPS pages may call `http://localhost` in Chrome, Edge and Firefox (the trustworthy-origin
exemption); Safari may block it. Running the app locally (`npm run dev`, or the Flask/Docker mode)
always works with local providers.

## Tested with real models

The Copilot first ran against real local models on 2026-10-10: the static dev build
(`VITE_STATIC_MODE=true npm run dev`) with the **Demo** run (GaTa₄Se₈, 250 K), default settings
(temperature 0.2, analyses on), on an Apple M5 Pro with 64 GB. A Playwright script asked six
questions, each in a new conversation, and read the step cards, the answer, the Outcome chip and the
raw stream. Other jobs shared the machine and the models, so the times are upper bounds.

✓ the answer and its verdict hold; ~ the numbers hold but the answer adds a wrong claim;
✗ wrong. Times are per question, from send to the end of the answer.

| Question | qwen/qwen3.8-27b (LM Studio) | gemma4:26b (Ollama) |
|---|---|---|
| Has the run converged? | `assess_convergence` · achieved ✓ ("improving", slope −0.00052 per row) · 91 s | `assess_convergence` · achieved ✓ · 22 s |
| First Ta–Se g(r) peak; constraint-limited? | `find_gr_peaks` · achieved ✓ (2.52 Å, constraint 2.2 Å) · 140 s | no tool: the context holds it · achieved ✓ · 8 s |
| Most anharmonic site, largest U_iso | `list_sites` · achieved ✓ (ref 23, Ta: kurtosis 5.62, 0.0175 Å²) · 55 s | `list_sites` · achieved ✓ · 8 s |
| Do the Ga atoms have a preferred direction? | `analyze_displacement_directions` · achieved ✓ (no: −0.89σ, map 0.56σ) · 63 s | same call · achieved ✓ · 11 s |
| Se–Ta–Se angles against an ideal octahedron | `find_gr_peaks`, then `compute_bond_angles` at 2.33–2.71 Å; saw coordination 3.86 and widened to 2.33–2.98 Å itself (coordination 5.91). No answer within 15 min: one reasoning step took 10.5 min | before the fixes below: one call at 2.3–2.7 Å, coordination 3.75 read as "incomplete octahedra", achieved ✗ · 48 s. After: `find_gr_peaks`, then its two-peak window 2.33–2.96 Å (coordination 5.89), achieved ~ (says there are no near-180° angles; 20 % lie at 150–170°) · 95 s |
| Space group at 0.2 Å and 0.02 Å; what breaks? | `detect_symmetry` at 0.02 Å (0.2 Å from the context) · achieved ~ (F-43m and P-42₁m are right; calls P-42₁m chiral and says centrosymmetry is lost, but F-43m has none) · 321 s | no tool: the context's tolerance ladder · achieved ✓ · 45 s |

- **Tool calling works.** All 16 calls had valid JSON arguments that passed the schemas: no
  malformed call, no unknown tool, no repeated call, no budget overrun. Both models used the
  context when it already held the answer, as the prompt asks.
- **The verdict is a self-assessment.** Neither model ever wrote anything but *achieved*, even
  over a cut coordination shell. gemma4 often ends with a bare `Outcome: achieved`, without the
  reason. The step cards' checks are the part to trust.
- **Speed.** gemma4:26b (a mixture of experts) answered in 8–95 s. qwen3.8-27b, a dense reasoning
  model at LM Studio's default reasoning effort, took 1–5 min, and over 15 min on the bond-angle
  question.
- **qwen3:32b (Ollama)** was not run: other jobs held the memory it needs.

The runs found four bugs, now fixed:

1. The verdict stayed in the answer text without a chip when the model copied the prompt's
   backticks around the line (`parseOutcome` now accepts inline code and a heading marker).
2. `find_gr_peaks`'s first-peak window (2.33–2.71 Å) cut the TaSe₆ shell, whose second Ta–Se peak
   is at 2.74 Å, and no check said so. It now adds `first_two_peaks_window_A` and a note when the
   second peak overlaps the window.
3. `compute_bond_angles` said "ok" for such a window. It now warns when fewer than 60 % of the
   central atoms share the most common bond count.
4. A server error sent mid-stream (LM Studio ran out of GPU memory) was skipped, the turn ended
   empty, and the chat crashed rendering it (`turn.items is not iterable`). The error now shows,
   and an empty reply shows its note.

### A scientific review (2026-10-10, later the same day)

The same Demo and the 5 K GaNb₄Se₈ run (`data/5K_try1`, loaded as a folder) were then put through
referee-style questions with gemma4:26b, and the answers were judged against the published
structures of both compounds (GaTa₄Se₈: Phys. Rev. Research 4, 033123 (2022); GaNb₄Se₈: Phys. Rev.
B 109, 144101 (2024)). The first round was fluent and wrong in specific, fixable ways; the fixes
(new context fields, a first-shell analysis, three new analyses — `compare_fit`, `analyze_clusters`,
`analyze_modulations` — and prompt rules) are in the changelog. The third round, after the fixes:

| Question | Before | After (tools called · answer) |
|---|---|---|
| Referee review: datasets, Rwp, convergence, sampling, constraints | "two datasets", "converged — a plateau", "no peak constraint-limited" | none (the context holds it) · one measurement fitted in both spaces; still improving (−0.11 over the last tenth); 3.5 moves/atom; Ga–Se tail shaped by the 2.2 Å constraint (g at 8.7 % of the peak); Ga weakly constrained (~5 % of the signal) ✓ |
| Is the Ta₄ cluster distorted? | "peak shifted to 2.93 Å: distorted" (no tool) | `find_gr_peaks` → `analyze_clusters` with its shell window · 23 % *less* distorted than independent atoms (z = −64), breathing +27 %, weak <110> elongation and <111> compression ✓ (it read the mean edge ladder as the range of bonds — the tool now says what the ladder is) |
| Preferred Ta displacement direction; trigonal or tetragonal? | pooled map isotropic → "no distortion" | `list_sites`, pooled directions (with the new note), `analyze_clusters` · <111> compression preferred (19 % vs 14 %, z = 11), no <001> signal ✓ |
| Is GaSe₄ regular? | window 2.21–2.63 Å (CN 3.87), "± 0.02 Å" invented | `find_gr_peaks` → `compute_bond_angles` at the shell window 2.18–2.79 Å · CN 3.99, 2.434 Å (sd 0.095 measured), 109.4° (sd 4.1°), constraint caveat ✓ |
| Where does the model miss G(r)? | "the context lacks it" (partial) | `compare_fit` · 62 % of the misfit at 2.34–5.28 Å (Rw 0.024 there), largest residuals at 2.53, 3.66, 3.82, 4.9 Å ✓ |
| What can the data say about Ga? | — | none · weak weights, light sampling, constraint-shaped Ga–Se tail; trust only the average site ✓ |
| Is the average structure still cubic? | lower rungs read as symmetry breaking | none · F-43m; the P-42₁m … P1 rungs lie below the 0.031 Å noise floor ✓ |
| 5 K: does the configuration show the reported P2₁2₁2₁ Nb₄ distortion? | — | `find_gr_peaks` → `analyze_clusters` · no: clusters 7 % less distorted than independent atoms, average cluster regular (3.025–3.035 Å), inconsistent with the report ✓ |
| 5 K: converged, and where in r and Q is the misfit? | — | `assess_convergence`, `compare_fit` (G(r) only) · improving; G(r) misfit at 10.8–15 Å ~ (did not check F(Q), yet called it achieved) |
| 5 K: cubic? A hidden superstructure? | no tool could say | `read_run_data`, `analyze_modulations` · F-43m; no zone-boundary order; the strongest modes are box-scale waves ✓ |

What a referee should take from the runs themselves: neither RMC configuration carries a
cluster distortion beyond independent atoms, nor the cell doubling diffraction finds at 5 K; both
carry box-scale displacement waves; the Ta–Ta and Ga–Se distributions show the closest-approach
and shell-gap effects above. The checks on the step cards remain the part to trust.

### The view tools with a real model (2026-10-11)

The docked Copilot's view tools were checked with **gemma4:26b** (Ollama, OpenAI-compatible
endpoint) on the static build with the Demo run, in Chrome at 1440 × 900 (docked), each question
in a new conversation. Other jobs shared the model server, so the times are upper bounds.

| Question | Tools called | Result |
|---|---|---|
| First Ta–Se g(r) peak; highlight the card and mark it | `highlight_cards` (the partials card, a mark at 2.52 Å) | ✓ the card breathes, fully in view beside the panel; the plot zooms to the marked peak; the answer cites the card as a chip · 223 s |
| Largest U_iso; cite the PCA Ellipsoid card | `list_sites` | ✓ ref 23, 0.0175 Å²; from the Dashboard, its `pca-stats` chip opened PCA Ellipsoid and brought the statistics card into view · 342 s |
| The direction sphere of Ta site 1 down [111] | `list_sites`, `analyze_displacement_directions`, `orient_3d_view` (select site 1, along `[1 1 1]`) | ✓ the page opened, the sphere turned to [1 1 1] with [-1 -1 2] up (a, b, c at 120°) · 847 s |
| Centre site 1 in the site view and ring it | `orient_3d_view` (`direction-sites`, centre and ring site 1) | ✓ site 1 at the canvas centre with the blue halo; the model wrote only the Outcome · 434 s |
| What am I looking at? | `current_view` | ✓ the page, the sphere down [1 1 1] (site 1) and the site view down [3 2 3], as on screen · 553 s |

No page error and no failed step. Two display bugs were found and fixed: the Outcome note showed
LaTeX raw (`$U_{iso}$`; it now renders as Markdown inline), and chips written as the bare
reference (`<card:direction-sphere>`) read "card:direction-sphere" (they now show the card's
title; a re-run of the last question confirmed it). The model's numbers for the sphere differ
from the drawn map because the tool lets the engine choose the resolution while the page keeps
its own (ν = 10, smoothed), as noted under *Show in* in
[ai-assistant.md Step 17](../../../../docs/algorithms/ai-assistant.md#step-17--the-copilot-tool-loop-analyses-checks-budgets-verdict).

## Import boundary and extraction recipe

This module is designed to be lifted into its own repository:

- Run data and engine tools enter **only as props** to `CopilotChat` (`runName`, `plotFiles`,
  `rValueFile`, `structure`, `symmetry`, `runSettings`, `pcaSites`, `liveData`, `tools`, `pages`,
  `onOpenPage`, `runKey`, `configKey`, `onOpenChange`, `cards`, `onMarks`, `onLayoutChange`, `className`,
  `storageKey`, `skills`) and `WatchdogBadge` (`rValueFile`). The host renders `CopilotChat` once, as the
  last child of its app grid, and places it by its `data-mode`: a column right of the pages
  (`docked`) or a row under them (`sheet`). The page cards it may point at carry
  `data-copilot-card`. Nothing in `src/llm/` imports
  `browserData.js`, `symmetryModel.js`, `useSiteCloud.js`, or any host component. Cell math is
  deliberately duplicated from `ModelSummary.jsx`.
- The host imports only `src/llm/index.js`.
- Styling: the module's own CSS lives in `components/Copilot.css` (the conversation, settings,
  watchdog badge, the breathing `[data-copilot-focus]` edge and the card chips) and
  `components/CopilotChat.css` (the dock, rail, sheet and resize handle). A kit token it reads
  keeps its old literal as a fallback (`var(--warning, #b45309)`). It uses exactly one host UI-kit
  class, **by class name only** and never kit JS (see `src/ui/README.md`): `ui-pill`, on a step's
  Show in / Open button. An extracted copy needs a host stylesheet that defines it. The sheet's
  height is a share of `--viewport-h` (100vh without it).

To extract: copy `src/llm/`, give your app grid a column and a row for the dock, render
`CopilotChat` once with props from your own data source, tag your cards with
`data-copilot-card`, and pass your own tools (and view tools, if you have 3D views).

## Tests

`__tests__/` covers everything that doesn't need a live model:

- **Context and prompts:** context assembly and the character budget; the watchdog heuristics and
  reply parsing; the prompt texts; the domain rules' split (`systemPrompt.test.js`: every part of
  every ground-rule sentence classified, the classification pinned sentence by sentence, method
  `'off'` with every definition and conduct part and no inference part).
- **The turn** (`turn.test.js`, at its interface: a scripted model port, a fake ask-user and a fake
  clock): the order of the events and the single terminal event, why each request was made, the
  tools offered and the prompt for them (the tools-on rule, the reserved names, "when to ask"), the
  method (`'full'`: each description with its advice; `'off'`: no skills, no `EVIDENCE_RULES`, no
  method sentence, descriptions and checks without their advice, a check whose advice is not in its
  text failing its call), the
  history window (the dock's turns as history, compaction), the warmer retry at 0.6 or more through
  `providerModel` (and the API key never in a request event), `stopped` and `failed`, invalid input.
  `method.test.js`: a description with its advice put back or left out, a check's text with or
  without its advice, a tool as the turn offers it.
  `src/__tests__/modelVisible/turnRequests.test.js`: the turn's first request for a dock-shaped and
  a harness-shaped spec equals the snapshot files (it needs the host's run module, so it sits
  beside them: this folder imports nothing from the host).
- **The turn policy** (`turnPolicy.test.js`, at its interface: rows of events → steps with a fake
  clock and a stub judge): every ending, the last-request note once, `first` against `final`, a
  reminder and the last-request note in one request, the warmer retry, the waits off the turn's
  time (read at most once an event), when a send-back is allowed (the last request, out of time,
  tools off, the verification's own cap), the withdrawn answer, `no_tools`, the 11 `done` keys on
  every end, the TypeErrors, `stoppedByFailure` and `STOPPED_BY`; in the loop, the three ending
  paths share one key set and the dock folds the new shape as it folded the old ones.
  `src/__tests__/turnEndings.test.jsx` reaches every ending through the dock (`CopilotView`) and the
  harness from the same scripted model: both read the same `stoppedBy`, the dock shows its note,
  a stall is `'error'` and Stop `'stopped'` in both, and the empty paths' `record.final` is pinned.
- **Streaming:** SSE parsing, including tool-call fragments split mid-arguments and Gemini-style
  whole calls with `extra_content`; tool-call assembly and argument validation.
- **The context window** (`contextGuard.test.js`): the estimate errs high; the plan keeps a loaded
  context that fits, asks for a fitting one over a 4,096 default, keeps it across the
  conversation, caps it at the model's maximum and refuses beyond it; `planCopilotChat` decides
  from the probes alone (the session written only when the request is sent, which then goes out
  under that plan); server detection never probes a cloud URL; the native NDJSON stream
  (thinking, content, whole tool calls, usage, `length`, errors, a tools-unsupported 400 with its
  body) and the message conversion; LM Studio's refusal; a reply whose server read only part of
  the prompt stops the turn.
- **The model port** (`modelPort.test.js`, `providerModel` over fake servers, nothing mocked above
  fetch): `prepare` per server (the plan, the tool choice's form, a cloud never probed, the key
  never in `sent`), Ollama's final request without tools and its estimate without their schemas,
  the refusal in `prepare` with nothing sent; `stream` under the plan (Ollama's `num_ctx` and
  `num_predict`; LM Studio and the cloud with the body as before); `ToolsUnsupportedError` from
  Ollama's and OpenAI's shapes and the fallback text, and never from a schema 400, a 5xx, a request
  without tools, a failure mid-reply or an abort; the loop over the port (an Ollama turn whose
  final reply is the answer; a model without tools asked again without them). In `turn.test.js`,
  a refused request has no `request` event.
  `copilotChat.test.jsx` runs a turn on a mocked Ollama loaded at 4,096 (native requests with
  `num_ctx` ≥ 16,384 and no `/v1` call), a refused turn, and a cut reply.
- **How a call ends** (`tools.test.js`, the tool contract at its interface): every result status
  and who decides it, with the exact texts the model reads (the skip reasons, the unknown tool, the
  repeat's note), the constructors' TypeErrors, an interactive tool's `choice` passed through and
  no other's, the tool message and its compaction, `severityOf`, `ran` and the budget helpers; and
  the invariant table over one turn that reaches every status, a repeat, the model's question and
  the Workbench's. `src/__tests__/toolOutcomes.test.jsx` runs one scripted turn on a dev case
  through the dock (`CopilotView`) and the harness alike — every status, each step's end the same
  in both — and pins the drift list: the dock's footer beside `failed_calls` and `refused_calls`,
  `unknown_tool_calls`, "asked before", a skipped question on the card and in the screen reader's
  line, the note on the step and in the record, provenance against what the model was told, a
  step running at a stop listed as stopped. The validation's scripted records are pinned as files
  (`validation/__tests__/__snapshots__/records/`, the transcript's added keys left out).
- **The agent loop, against a scripted model:**
  - the call → result → answer cycle; argument errors fed back and corrected;
  - `ask_user` replays on the Demo's Ta–Se shell: an ambiguous Se–Ta–Se question finds the shells,
    asks which window, and the analysis runs in the one chosen (an option, the default, a typed
    window, or Skip); a clear question runs without asking; Stop while the question waits; no time
    limit and the wait off the turn clock; one question per step, two per request, never cached;
  - unknown tools, thrown errors, duplicate calls, the per-step limit;
  - the final `tool_choice: 'none'` round and an empty last round;
  - the tools-unsupported fallback on the port's typed error only (an untyped error with the same
    text is thrown), and a schema 400 that must surface;
  - a request refused in `prepare` (no `request` event) and a `prepare` that hangs (the stall
    clock);
  - abort and timeout, which also abort the tool's signal;
  - per-result and per-turn budgets.
- **Skills:** skill files (front matter, references, a misnamed skill), `read_skill` (its budget,
  references on demand, never compacted) (`skills.test.js`); the repository's own skills parse,
  are drafts with every section, and reach the agent whole (`src/__tests__/copilotSkills.test.js`);
  one replay per skill on the Demo (`copilotSkillReplays.test.js`); the cut-shell refusal and its
  confirmed escape (`copilotTools.test.js`). The skill lint (`src/__tests__/skillLint.test.js`,
  tables in `skillLint/tables.js`) reads every skill and reference file and never edits them: every
  tool a skill names is one the builders offer; every gate number equals the constant that decides
  it (each number after a comparator is bound or listed as no gate); every "*Firm rule in the app:*"
  is mapped to the check or refusal that enforces it, and that code is exercised; every value a skill
  quotes as measured on the Demo that `VALUES` binds (the shell readings: `find_gr_peaks`, the shell
  catalogue, the bond-angle counts summary) equals the tool's there, rounded as the skill writes it.
  Today's violations sit in its allow-list with why (a value entry with its line, the skill's value
  and the tool's: the coordination-shells skill's old first peak, 2.33–2.71 Å and 3.9, against
  2.21–2.65 Å and 3.25); a new one fails, and so does an entry a fix cleared.
- **Checks and the Workbench's questions:** the verification at its interface, steps recorded and
  the answer judged — the answer checks with their negated and qualified phrasings
  (`answerChecks.test.js`); the evidence gate under both methods, the fixed order, numbers by
  reference, the computed verdict, the errors (`verification.test.js`); the catalogues and the
  binding of `EVIDENCE_RULES` to CLAIMS (`verificationCatalogue.test.js`); the reference store
  (`references.test.js`) and the dock's reducer and footer (`turnState.test.js`); the real tools'
  facts against FACTS (`src/__tests__/verificationFacts.test.js`); the verdict's words through
  the harness's phrasings (`validation/__tests__/verdictPhrasings.test.js`); and in the loop
  (`agentLoop.test.js`); the clarifier on synthetic shells — split,
  merged, clean, word-like element symbols, the tolerance ladder, typed windows
  (`clarify.test.js`); and on the Demo with the real engines and context
  (`src/__tests__/copilotChecksReplays.test.js`, each question through the turn module, as in the
  app; two loop-level cases said in place): a
  contradicted answer is not verified, gemma4's correct answers stay *achieved* and are never sent
  back, a cluster conclusion sent back for its analysis and its numbers filled from it, each Workbench question (an element alone, a split shell, a
  cut window kept or replaced, one site or pooled), the clear questions that never ask, the
  model's own `ask_user` beside them within two per request, and cancel. In the dock
  (`copilotChat.test.jsx`): the Workbench's card before any model request, the lowered chip.
  Real answers (`src/__tests__/answerCorpus.test.js`): the gate's claims and send-back and the
  checks that fire, pinned per answer with a judgement of each, and the false alarms counted
  ([copilot-validation.md](../../../../docs/algorithms/copilot-validation.md#the-answer-corpus-plan-item-35)).
- **Rendering:** Outcome parsing, the Markdown pipeline with its card chips, and
  `copilotChat.test.jsx`, which drives the dock end to end against a mocked server: nothing is
  fetched before it opens, Show in hands over the page parameters, Esc keeps the conversation, a
  run change restarts it, an answer that lands while collapsed marks the rail, the open state and
  width are remembered and resize from the keyboard, and the agent's highlight_cards makes a host
  card breathe (with marks), its chip brings the card back from another page, and a click on the
  page clears it. The `ask_user` card: it waits on one request, takes focus, moves with the arrow
  keys, sends the option picked by its number key back to the model, which computes in that window;
  *Other…*; Stop while waiting; a run change drops the question; a collapsed dock marks the rail.
- **What the model sees** (`src/__tests__/modelVisible/`): snapshot files of every system-prompt
  variant (with method `'off'` too), the dock's and the validation harness's first request and tool
  list (both through the turn module; arms B and C, arm C's request checked for every skill name,
  skill step, method sentence and piece of tool advice), one tool message per tool on the Demo
  (and the two whose checks carry advice, as arm C reads them), and every tool message no other
  file covers, raw (`toolMessages.test.js`: a call that did not run, the questions, unreadable
  arguments, a throw, the time limit, a repeat, the Workbench's notes, the budget, the compaction,
  `highlight_cards`). A diff there is a change to what the model sees: call it out in the PR
  ([copilot-validation.md](../../../../docs/algorithms/copilot-validation.md#what-the-model-sees-snapshot-tests)).
- **The shape of each tool's data:** the notation and its checker (`dataShape.test.js`), and every
  live tool's results on the Demo and the dev cases against its `dataShape`
  (`src/__tests__/toolDataShapes.test.js`, *The tool contract*).
- **Pointing:** the focus controller (`focus.test.js`: levels, page switch and wait, clearing,
  undo callbacks, re-rendered cards), card references (`cardRefs.test.js`), `highlight_cards`
  (`viewTools.test.js`), the dock's size rules (`dockLayout.test.js`), and list/object arguments
  (`toolCalls.test.js`).

The host tools are tested on the committed demo run with the real engines in
`src/__tests__/copilotTools.test.js`; the 3D view tools and their math in `viewer3d.test.js`,
`cameraPose.test.js`, `useCopilotViewer.test.jsx`, `copilotViewTools.test.js` and
`copilotCards.test.js` (which also checks that every card in the catalog is tagged in its page),
and the plot marks in `interactivePlotMarks.test.jsx`. Run everything with `npm test` (vitest).
