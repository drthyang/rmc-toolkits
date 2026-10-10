# AI Copilot (experimental LLM module)

The **AI Copilot** connects the RMCProfile Workbench to an **OpenAI-compatible** LLM: a **local**
server the user runs (Ollama, LM Studio) or a **cloud** provider (OpenAI, Gemini, Anthropic) with an
API key. It replaced the chat-only AI Assistant on 2026-10-09. The Assistant answered from one fixed
summary of the run. The Copilot also **runs the Workbench analyses itself**, reads what came back,
and judges whether that answers the question.

It is a **chatbox over every page**, not a page of its own:

- A labelled launcher sits in the bottom-right corner. It opens a floating panel; *Larger* widens
  it, and ✕ or Esc closes it.
- A step's **Show in \<page\>** button switches the workbench page *behind* the chat and loads the
  same analysis there. Bond Geometry gets the triplet and windows and computes; Displacement
  Directions gets the site and resolution. The page shows its own result next to the conversation.
- On a phone the panel is a full-height sheet. It closes after Show in, so the page is visible.

How one question runs:

```
        question
           │
   ┌───────▼────────┐   tools offered    ┌──────────────┐
   │  agent loop    │ ─────────────────▶ │    model     │
   │ (agentLoop.js) │ ◀───────────────── │ (client.js)  │
   └───────┬────────┘  tool calls / text └──────────────┘
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

**The tools only compute.** None switches a page, moves a slider, changes the Detected-SG
tolerance, or writes or downloads a file. A step card's **Show in \<page\>** (or **Open
\<page\>**) button is the *user's* action.

## How each seam works (reading order)

1. **`context/runContext.js` + `context/pairCorrelations.js`** build the run context.
   `buildRunContext()` distills the dashboard's parsed state into a compact JSON object:
   - cell and angles;
   - the symmetry block: space group, ladder, and Wyckoff orbits ranked by rms displacement;
   - pair correlations;
   - move counters;
   - run settings from the `.dat` file;
   - per-dataset Rwp;
   - the ln(χ²) history, downsampled to ≤ 48 points with statistics computed on the full series.

   `contextToJson()` enforces a ~4,500-character budget and records each trim with an `*_omitted`
   count. The context travels with every question, so simple questions need no tool.
2. **`prompts/system.js` + `prompts/copilot.js`**
   - `SYSTEM_PROMPT` holds the domain rules: what Rwp means, that the history is ln(χ²) of one fit
     term, and "quote the numbers, don't guess".
   - `COPILOT_PROMPT` adds how to work with tools: say why before each call, read the `checks`,
     change the arguments on a `fail`, never repeat an identical call. It also asks for the verdict
     line.
   - `buildCopilotMessages()` puts the context, inside a fenced JSON block, ahead of the
     conversation. Earlier turns travel as question + final answer only (the newest 8).
   - `prompts/templates.js` holds the watchdog's rigid `STATUS: …` prompt and its parser.
3. **`provider/client.js`** is a hand-rolled OpenAI-compatible client (`/models`,
   `/chat/completions`) with an SSE stream parser over `fetch` + `ReadableStream`.
   - `streamChat()` yields `{ content }`, `{ reasoning }`, `{ toolCall }` fragments and
     `{ finishReason }`.
   - `tools` / `tool_choice` join the request body only when given, so watchdog requests are
     unchanged.
   - Each fragment's provider data (`extra_content`, where Gemini puts its `thought_signature`)
     is kept, because it must be echoed back.
   - An HTTP failure carries its `status`.
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
   A thrown error becomes a failed step whose message the model reads. `fitToBudget()` halves the
   largest array, recording `*_omitted`, until a result fits. It never trims silently.
6. **`copilot/agentLoop.js`**: `runCopilotTurn()` is a pure async generator of UI events. The model
   call is injected, so the tests script the model's replies. Its rules:
   - **Budgets.** ≤ 6 rounds of tool calls; ≤ 4 calls per round (extra calls are answered
     "skipped"); 120 s per call; no new round after 10 min. Each result is fitted to 3,000
     characters, and a turn's tool messages to 12,000 in total. Older results are compacted to
     their checks first, so a small local context window keeps the run context and the rules.
   - **Duplicate calls.** An identical call (same validated arguments) returns the cached result
     with a note.
   - **Final round.** When the budget is used, a last request goes out with `tool_choice: 'none'`.
     If it still answers with tool calls only (Ollama ignores `tool_choice`), a fixed
     "budget used" answer with `Outcome: not achieved` closes the turn.
   - **No tool support.** A 4xx saying the model "does not support tools" retries the turn once
     without tools, under the tool-free prompt. Other 400s, such as a schema error, are shown.
   - **Message rules.** Every listed call gets exactly one `tool` message. Arguments are echoed as
     valid JSON (`'{}'` when unparseable). Whitespace-only assistant content is sent as `null`
     (Anthropic rejects empty text blocks).
7. **`copilot/runTools.js`** holds the tools over data the module already receives:
   - `read_run_data(section)`: one context section as built, before the budget trim;
   - `assess_convergence`: the watchdog heuristics plus tail slopes;
   - `find_gr_peaks(pair)`: partial-g(r) peaks, a first-peak bond window, and the pair's
     minimum-distance constraint with a "constraint-limited" check.
8. **`copilot/outcome.js`**: `parseOutcome()` splits the verdict line from the answer. It is
   lenient about Markdown and dashes, strict about the three statuses, and takes the last line.
   **`copilot/turnState.js`** folds loop events into a turn: ordered text segments (with
   reasoning) and tool steps.
9. **`watchdog/`**: slope heuristics over the ln(χ²) history are the **source of truth** for the
   badge. The LLM only writes the tooltip note, and only when the status flips or enough new data
   arrived.
10. **`useAssistant.js`, `useCopilot.js` and `components/`**
    - `useAssistant()` holds the settings, the connection probe (with auto-connect) and the run
      context.
    - `useCopilot()` holds the conversation: Stop, which models cannot call tools, and a new
      conversation when `runKey` changes. It also builds the memory (`conversationMemory`): a
      question without an answer is dropped.
    - `CopilotChat` is the chatbox: the launcher, then the panel. The panel mounts on first open,
      so nothing is fetched before that, and stays mounted while closed, so a running answer
      continues. It composes the connection bar, the settings drawer (`ConnectionSettings`, with
      the analyses switch) and `CopilotView`. Focus goes to the composer on open (once connected)
      and back to the launcher on close.
    - `CopilotView` renders the transcript (`ToolStep` cards, `MarkdownMessage`, the Outcome chip),
      the composer, and a polite live region that announces step and verdict changes.

## The tool contract (the `tools` prop)

The host passes engine tools as plain objects; the module validates and runs them:

```js
{
  name: 'compute_bond_angles',           // /^[a-z][a-z0-9_]{0,63}$/, unique
  description: '…',                      // what it computes and when to use it, for the model
  parameters: { type: 'object', properties: {…}, required: [...] },
  label: (args) => 'Bond angles Se–Ta–Se',   // step-card title (optional)
  page: 'geometry',                      // host page key for the step's button (optional)
  pageParams: (args) => ({ … }),         // what that page loads to show the same analysis (optional)
  run: async (args, { signal }) => ({ data, checks }),
}
```

- `data` is a compact JSON summary built for the model, never a raw engine payload. Put the most
  important rows first: budget trimming halves arrays from the end.
- `checks` is a list of `{ level: 'ok' | 'note' | 'warn' | 'fail', text }`, computed from the
  full result. `fail` means "this does not answer as asked".
- `signal` aborts on Stop or timeout. Stop the work if you can.

`pages` (`{ key: label }`) and `onOpenPage(key, params)` let a step offer **Show in \<page\>**.

- `params` is the tool's `pageParams(args)`, computed from the validated arguments. Without
  `pageParams`, the button reads **Open \<page\>** and `params` is null.
- `runKey` names the open run; when it changes, the conversation restarts.
- `onOpenChange(open)` reports the chatbox opening and closing.

In this repository the host tools live in `src/copilotTools.js`:

| Tool | What it runs |
|---|---|
| `list_sites` | the per-site PCA table, through the worker or `/api/pca/sites` |
| `analyze_displacement_directions` | the orientation engine with calibrated significances and `[u v w]` (Show in Displacement Directions selects the site and resolution) |
| `compute_bond_angles` | the triplets engine (Show in Bond Geometry fills the triplet and windows, then computes) |
| `detect_symmetry` | the symmetry finder at any tolerance, without changing the card (browser-loaded runs only) |

A browser-loaded run is analysed in a worker of the Copilot's own, which Stop or a timeout
terminates. A server folder goes to Flask with the request cancelled on abort.

## Provider setup

- **Ollama** (local): <https://ollama.com>. Pull a model that can call tools
  (`ollama pull qwen3`, `llama3.1`, …). Start the server with this site allowed as an origin, and
  with room for the analyses:

  ```bash
  OLLAMA_ORIGINS="https://drthyang.github.io" OLLAMA_CONTEXT_LENGTH=16384 ollama serve
  ```

  Ollama's OpenAI endpoint cannot set the context length per request. Its default window is small
  enough that a multi-step turn may lose its oldest messages.
- **LM Studio** (local): <https://lmstudio.ai>. Load a tool-capable model with a large context, enable
  **CORS** in the Developer/server settings, and start the local server
  (default `http://localhost:1234/v1`).
- **OpenAI / Gemini / Anthropic** (cloud): pick the provider and paste an API key. There is no server
  to run, but your run context and analysis results are sent to that provider. Some providers or
  networks block direct browser calls (CORS).

In the app: click the **AI Copilot** button in the bottom-right corner, open the settings (gear),
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

## Import boundary and extraction recipe

This module is designed to be lifted into its own repository:

- Run data and engine tools enter **only as props** to `CopilotChat` (`runName`, `plotFiles`,
  `rValueFile`, `structure`, `symmetry`, `runSettings`, `pcaSites`, `liveData`, `tools`, `pages`,
  `onOpenPage`, `runKey`, `onOpenChange`) and `WatchdogBadge` (`rValueFile`). The host renders
  `CopilotChat` once, over its pages. Nothing in `src/llm/` imports
  `browserData.js`, `symmetryModel.js`, `useSiteCloud.js`, or any host component. Cell math is
  deliberately duplicated from `ModelSummary.jsx`.
- The host imports only `src/llm/index.js`.
- Styling: the module's own CSS lives in `components/Copilot.css` (the conversation, settings,
  watchdog badge) and `components/CopilotChat.css` (the launcher and the floating panel, z-index
  40). A kit token it reads keeps its old literal as a fallback (`var(--warning, #b45309)`). It
  uses exactly one host UI-kit class, **by class name only** and never kit JS (see
  `src/ui/README.md`): `ui-pill`, on a step's Show in / Open button. An extracted copy needs a
  host stylesheet that defines it. On phones the panel's sheet height subtracts
  `--header-pinned-h`, the host's pinned-header strip (0 without it).

To extract: copy `src/llm/`, render `CopilotChat` once with props from your own data source, and
pass your own tools.

## Tests

`__tests__/` covers everything that doesn't need a live model:

- **Context and prompts:** context assembly and the character budget; the watchdog heuristics and
  reply parsing; the prompt builders and conversation memory.
- **Streaming:** SSE parsing, including tool-call fragments split mid-arguments and Gemini-style
  whole calls with `extra_content`; tool-call assembly and argument validation.
- **The agent loop, against a scripted model:**
  - the call → result → answer cycle; argument errors fed back and corrected;
  - unknown tools, thrown errors, duplicate calls, the per-step limit;
  - the final `tool_choice: 'none'` round and an empty last round;
  - the tools-unsupported fallback, and a schema 400 that must surface;
  - abort and timeout, which also abort the tool's signal;
  - per-result and per-turn budgets.
- **Rendering:** Outcome parsing, the Markdown pipeline, and `copilotChat.test.jsx`, which drives
  the chatbox end to end against a mocked server: nothing is fetched before it opens, Show in
  hands over the page parameters, Esc keeps the conversation, and a run change restarts it.

The host tools are tested on the committed demo run with the real engines in
`src/__tests__/copilotTools.test.js`. Run everything with `npm test` (vitest).
