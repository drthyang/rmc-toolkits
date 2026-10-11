# AI Copilot (experimental LLM module)

The **AI Copilot** connects the RMCProfile Workbench to an **OpenAI-compatible** LLM: a **local**
server the user runs (Ollama, LM Studio) or a **cloud** provider (OpenAI, Gemini, Anthropic) with an
API key. It replaced the chat-only AI Assistant on 2026-10-09. The Assistant answered from one fixed
summary of the run. The Copilot also **runs the Workbench analyses itself**, reads what came back,
and judges whether that answers the question.

It is **docked beside the pages**, not over them, and it **points at the page** while it explains:

- **The dock.** On a window at least 1080 px wide (and taller than 540 px: laptops, monitors,
  tablets in landscape) the Copilot is a side
  panel in its own grid column: the pages shrink and reflow beside it (their width tiers are
  container queries), so it never covers a card. Drag its left edge (or focus the edge and use
  ←/→, Home/End; double-click or Enter resets) to resize it; the pages always keep at least
  760 px. On tablets in portrait, phones and short windows it is a **bottom sheet** in its own grid row,
  resized from its top edge; when the agent brings a card on screen, the sheet comes down (to 28 %
  of the window at most) until the card fits whole above it. Collapsed (the chevron, or Esc inside it) it is a slim **rail** at the
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
   - `find_gr_peaks(pair)`: partial-g(r) peaks and the first coordination shell — its window,
     mean distance, coordination numbers and their sensitivity to the cut, whether it is separated
     from the next shell — and whether the pair's minimum-distance constraint cuts the distribution;
   - `compare_fit(dataset)`: where the model misses the data, as Rw and the share of the misfit
     per range of r or Q, and the largest residual features.
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
  cards: (args) => ['bond-angles'],      // host cards this call reads: they breathe while it runs (optional)
  view: true,                            // a view tool: changes only the screen; never cached (optional)
  run: async (args, { signal, focus }) => ({ data, checks }),
}
```

- `data` is a compact JSON summary built for the model, never a raw engine payload. Put the most
  important rows first: budget trimming halves arrays from the end.
- `checks` is a list of `{ level: 'ok' | 'note' | 'warn' | 'fail', text }`, computed from the
  full result. `fail` means "this does not answer as asked".
- `signal` aborts on Stop or timeout. Stop the work if you can.
- `focus` is the dock's card-focus controller, for view tools: `await focus.show(ids, { marks })`,
  `focus.onClear(undo)` (runs when the highlight ends), `focus.state()`.
- Parameters may be arrays (`items`, `minItems`, `maxItems`) and objects (their own `properties` /
  `required`); a list sent as JSON text or a bare comma list is accepted.

`pages` (`{ key: label }`) and `onOpenPage(key, params)` let a step offer **Show in \<page\>**.

- `params` is the tool's `pageParams(args)`, computed from the validated arguments. Without
  `pageParams`, the button reads **Open \<page\>** and `params` is null.
- `runKey` names the open run; when it changes, the conversation restarts.
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
| `list_sites` | the per-site PCA table, through the worker or `/api/pca/sites` |
| `analyze_displacement_directions` | the orientation engine with calibrated significances and `[u v w]` (Show in Displacement Directions selects the site and resolution) |
| `compute_bond_angles` | the triplets engine (Show in Bond Geometry fills the triplet and windows, then computes) |
| `analyze_clusters` | the cluster engine (`workers/clusters.js`, `/api/clusters`): the shape of an element's finite clusters (M₄ tetrahedra, dimers, trimers) against a shuffle null of independent atoms, and the average cluster's own distortion |
| `analyze_modulations` | the modulation engine (`workers/modulations.js`, `/api/modulations`): each element's displacement power at the zone-boundary points and its strongest modes against the configuration's diffuse background — order that enlarges the cell, which the folded average cannot show |
| `detect_symmetry` | the symmetry finder at any tolerance, without changing the card (browser-loaded runs only) |
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
its own (ν = 10, smoothed), as noted under *Show in*.

## Import boundary and extraction recipe

This module is designed to be lifted into its own repository:

- Run data and engine tools enter **only as props** to `CopilotChat` (`runName`, `plotFiles`,
  `rValueFile`, `structure`, `symmetry`, `runSettings`, `pcaSites`, `liveData`, `tools`, `pages`,
  `onOpenPage`, `runKey`, `onOpenChange`, `cards`, `onMarks`, `onLayoutChange`, `className`,
  `storageKey`) and `WatchdogBadge` (`rValueFile`). The host renders `CopilotChat` once, as the
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
- **Rendering:** Outcome parsing, the Markdown pipeline with its card chips, and
  `copilotChat.test.jsx`, which drives the dock end to end against a mocked server: nothing is
  fetched before it opens, Show in hands over the page parameters, Esc keeps the conversation, a
  run change restarts it, an answer that lands while collapsed marks the rail, the open state and
  width are remembered and resize from the keyboard, and the agent's highlight_cards makes a host
  card breathe (with marks), its chip brings the card back from another page, and a click on the
  page clears it.
- **Pointing:** the focus controller (`focus.test.js`: levels, page switch and wait, clearing,
  undo callbacks, re-rendered cards), card references (`cardRefs.test.js`), `highlight_cards`
  (`viewTools.test.js`), the dock's size rules (`dockLayout.test.js`), and list/object arguments
  (`toolCalls.test.js`).

The host tools are tested on the committed demo run with the real engines in
`src/__tests__/copilotTools.test.js`; the 3D view tools and their math in `viewer3d.test.js`,
`cameraPose.test.js`, `useCopilotViewer.test.jsx`, `copilotViewTools.test.js` and
`copilotCards.test.js` (which also checks that every card in the catalog is tagged in its page),
and the plot marks in `interactivePlotMarks.test.jsx`. Run everything with `npm test` (vitest).
