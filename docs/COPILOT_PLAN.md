# AI Copilot — implementation plan

Status: **built 2026-10-09** (uncommitted at the time of writing), then turned from a page into
a chatbox the same day (§11). The plan below was reviewed
independently before and during the build; §10 records each finding and what changed. Where §10
and the sections above disagree, §10 wins. Replaces the chat-only
**AI Assistant** page with an **AI Copilot** that runs the Workbench's analyses itself, reads what
came back, and says whether the result answers the question.

## 1. Goal and non-goals

The Assistant answers from one fixed, budget-trimmed summary of the run (`buildRunContext`, ≤ 4 500
characters). It cannot look closer: it cannot compute a bond-angle distribution, test a site's
displacement directions, or re-run the symmetry finder at another tolerance. The user has to do
that on the analysis pages and paste the numbers back.

The Copilot closes that loop: **plan → act → observe → assess**.

1. The model decides which evidence the question needs.
2. It calls **tools**: the same engines the pages use (worker in the browser, Flask in server mode).
3. Each tool returns a compact result **plus deterministic checks** (e.g. "0 angles in these
   windows", "peak not significant, 1.2σ", "parse warning: 3 lines skipped").
4. The model reads the result and the checks. It adjusts and retries (a wider window, another site)
   or moves on.
5. It ends with an answer and a machine-readable verdict: `Outcome: achieved | partial |
   not achieved — <reason>`, shown as a chip.

Non-goals for this change:

- **No model-initiated changes to the app.** Tools compute; they never switch pages, move a slider,
  change the shared Detected-SG tolerance, write or download files. A tool step offers an
  **Open page** button the *user* presses. The Copilot page stays mounted while hidden, so a
  model-driven page switch would hide the transcript mid-run. Pre-filling a page's form from a
  step (e.g. Bond Geometry's triplet + windows) is phase 2 (§9).
- No new backend routes and no server-side LLM calls. The Copilot stays browser-only, like the
  Assistant.
- No persistence of conversations across reloads.
- The convergence watchdog (`watchdog/`, `WatchdogBadge`) is unchanged.

## 2. Architecture

```
host (src/)                                       llm module (src/llm/) — imports nothing from host
───────────────────────────────────────────       ──────────────────────────────────────────────────
App.jsx ── props ──────────────────────────────▶  CopilotPage
  run data (as today: runName, plotFiles,           ├─ useAssistant   (settings, probe, run context)
  rValueFile, structure, symmetry, runSettings,     ├─ useCopilot     (turn state, Stop)
  pcaSites, liveData)                               │    └─ agentLoop (pure async generator)
  + tools  ◀── useCopilotTools (copilotTools.js)    │         ├─ provider/client.js (SSE + tool_calls)
  + onOpenPage(page)                                │         └─ copilot/tools.js (validate, run, budget)
                                                    └─ CopilotView → ToolStep, MarkdownMessage
useCopilotTools ── requestAnalysis (useSiteCloud.js) ─▶ pcaKdeWorker | /api/pca/*, /api/triplets
               └─ describeSymmetry (symmetryModel.js)
```

**Import boundary kept.** Tools enter the llm module **as a prop** (an array of tool specs with
`run` functions), exactly like run data. Nothing in `src/llm/` imports host code. Host-free tools
over data the module already receives (convergence, datasets, pair correlations, run settings) live
inside the module (`copilot/runTools.js`).

### 2.1 Tool contract (the prop interface)

```js
{
  name: 'compute_bond_angles',          // /^[a-z][a-z0-9_]{0,63}$/, unique
  description: '…',                     // for the model: what it computes, units, when to use it
  parameters: { type: 'object', properties: {…}, required: [...] },  // JSON-Schema subset (§2.3)
  label: (args) => 'Bond angles Se–Ta–Se',  // UI title of the step card
  page: 'geometry',                      // optional: the page key the Open page button opens
  run: async (args, { signal }) => ({ data, checks }),
}
```

- `data`: a compact, JSON-safe object for the model, built by the tool's own summarizer, never the
  raw engine payload (a 48³ KDE volume or 642 sphere cells are useless to a model and blow the
  budget).
- `checks`: `[{ level: 'ok' | 'note' | 'warn' | 'fail', text }]`, deterministic, computed from the
  full result. They travel to the model inside the tool message and render on the step card. This
  is the grounded half of "assess the outcome". The model's `Outcome:` line is the other half.
- A thrown error becomes `{ error: message }` for the model and a failed step card. The loop
  continues, so the model can correct its arguments.

### 2.2 Agent loop (`copilot/agentLoop.js`)

A pure async generator, `runCopilotTurn({ callModel, messages, tools, limits, signal })`. It yields
UI events and never touches React. `callModel` is injected (the real one wraps `streamChat`),
which makes it testable with scripted model replies.

```
round = 0
loop:
  stream one completion with `tools` (tool_choice 'auto')
    yield text / reasoning deltas
  if no tool calls → final; parse Outcome; done
  round += 1
  append the assistant message (content or null, tool_calls)
  for each call (≤ MAX_CALLS_PER_ROUND; extra calls → error result "skipped"):
    yield tool_start
    parse + validate args → run with timeout + abort → budget the JSON
    duplicate (same name + canonical args this turn) → cached result + note
    append { role: 'tool', tool_call_id, content: JSON }
    yield tool_end (status from checks: fail > warn > ok)
  if round == MAX_ROUNDS → one last completion with tool_choice 'none' and a
    "tool budget used, answer now" note; tool calls in it are ignored
```

Limits (constants, overridable for tests): `MAX_ROUNDS = 6`, `MAX_CALLS_PER_ROUND = 4`,
`TOOL_TIMEOUT_MS = 120 000`, `TOOL_RESULT_CHAR_BUDGET = 3 000` per result (arrays trimmed with
`*_omitted` counts first, a hard cut with `truncated: true` last), and the user's Stop button
(abort signal → stream aborted, a pending tool is abandoned, partial transcript kept).

**Models without tool support.** If the first request fails with an HTTP 400 that names tools /
function calling (Ollama: "… does not support tools"), the turn is retried once **without** tools.
The session is marked `toolsUnsupported` for that model (cleared on a model switch), and the page
says so in one line: "This model cannot run analyses — answering from the run summary only."

**Conversation memory.** Earlier turns travel as plain user / final-assistant text only (no tool
messages; the newest 8 turns, as today). Within a turn the full tool exchange is sent. The run
context block is sent at the start of every turn, as today.

### 2.3 Streaming tool calls (`provider/client.js`, `copilot/toolCalls.js`)

- `postChat` adds `tools` / `tool_choice` to the body **only when given**, so the watchdog and
  plain requests are byte-identical to today.
- `streamChat` additionally yields `{ toolCall: { index, id, name, arguments } }` fragments and
  `{ finishReason }`. Existing consumers ignore unknown keys.
- `accumulateToolCalls(state, fragment)` (pure) assembles fragments. The key is `index`; servers
  that omit it (Gemini's compat layer sends whole calls) are keyed by `id`, else by arrival order.
  `arguments` strings are concatenated. A missing id gets a generated `call_<n>`.
- Tool-call presence, not `finish_reason`, decides "the model wants tools" (Ollama reports `stop`).
- Arguments: `''` → `{}`; invalid JSON → an error result to the model naming the parse error.
  Validation covers the subset the tools use: `type` (object / number / integer / string /
  boolean), `enum`, `minimum` / `maximum`, `required`. Unknown keys are dropped with a note.
- Assistant messages carrying tool calls send `content: null` when empty. Anthropic's compat layer
  rejects an empty text block, and OpenAI / Gemini / Ollama accept null.

### 2.4 Tools (v1)

All tools are read-only computations. Availability is decided per run: a tool whose inputs are
absent is simply **not offered** to the model (`available` filtering in the host hook / runTools).

| Tool | Where | Engine | Returns (compact) | Checks (examples) |
|---|---|---|---|---|
| `read_run_data(section)` | llm | `buildRunContext` sections, **untrimmed** | the full section: `convergence`, `datasets`, `pair_correlations`, `run_settings`, `symmetry`, `pca_displacements`, `structure`, `configuration_optimization` | section absent → fail naming what the run lacks |
| `assess_convergence()` | llm | `watchdog/heuristics.js` (`classifyConvergence`, `watchdogStats`) + `seriesStats` | heuristic status, column, n steps, first/last/min, slope over the last 10 / 25 / 50 % | non-finite steps → warn; < 5 points → warn; status ≠ improving/converged → warn |
| `list_sites()` | host | `requestAnalysis('sites')` (worker or `/api/pca/sites`) | per site: ref, element (composition if mixed), count, U_iso, rms axes, anisotropy, Mardia κ, zeroSpread | `parseWarning` → warn; reconstructed sites → note; zero-spread sites → warn |
| `analyze_displacement_directions(site_ref \| element, frequency?)` | host | `requestAnalysis('orientation')` | used points, peak enhancement + **calibrated** `peakSignificance`, `mapSignificance`, antipodal asymmetry z (+ flag), Bingham anisotropy σ, orientation eigenvalues, peak direction as Cartesian unit vector **and** fractional `[u v w]` (direct-lattice direction; never called a plane normal) | used < 100 → warn; peak ≥ 3σ → ok "preferred direction"; < 3σ → note "no significant peak"; map test null → note (sparse) |
| `compute_bond_angles(end1, apex, end2, r12_min, r12_max, r23_min?, r23_max?, bin_width?)` | host | `requestAnalysis('triplets')` | angle count, mean ± sd, the 3 strongest sin-corrected modes (angle, enhancement), bond-length means + unique bonds per window, coordination histogram of the apex | 0 angles → fail "widen the windows"; > 50 % apex atoms with 0 bonds → warn; parseWarning → warn |
| `detect_symmetry(tolerance_A)` | host | `describeSymmetry(structure, tol)` | space group (+ No. only when the card would show it), point group, n ops, orbits with Wyckoff labels (`orbitLabel`) | skipped (caps) → warn with reason; undetermined / lower bound → note |

`detect_symmetry` is offered only when `structure.basis` exists (browser-loaded runs; Flask
server-folder mode has no Detected SG card either). The host tools are offered only when
`requestAnalysis` can run: a local file is loaded, or Flask mode.

Tool parameter names are snake_case for the model and mapped to the engines' request keys (`r12Min`,
`referenceNumber`, …) in the host hook. The engines' own guards (`requestNumber`, the 15 Å `rmax`
cap, `APP_MAX_ANGLES`) stay the single source of validation; the schema only catches obvious
mistakes early with a message the model can act on.

### 2.5 Prompts (`prompts/copilot.js`)

`COPILOT_PROMPT = SYSTEM_PROMPT + agent rules`:

- You can run analyses with the tools; prefer evidence over recall. Do not call a tool when the run
  context already answers.
- Before each tool call, say in one short sentence why. Read the `checks`. On `fail`, change the
  arguments or the approach; do not repeat an identical call.
- Quote numbers from tool results. Significance values are calibrated one-sided σ; `zScore`-type
  local values are not significances.
- Finish with the answer, then one final line exactly
  `Outcome: achieved|partial|not achieved — <one sentence>`.

`parseOutcome(text)` → `{ status, note, body }` (body = text without the line). It tolerates
bold / italics around the label and en/em dashes. A reply without the line renders without a chip.

### 2.6 UI

The **AI Copilot** tab replaces **AI Assistant** at the same place in the nav (page key `copilot`).
The look follows the existing module CSS (`llm-*` classes, kit classes by name only).

- **Header**: "AI Copilot" + Beta tag + the existing connection bar (status, model, gear).
- **Settings drawer**: the existing drawer plus one switch, "Run analyses (tools)", default on.
  The cloud warning gains: "…including the results of analyses the Copilot runs." Off means plain
  chat, which is today's Assistant.
- **Transcript**: user bubbles; an assistant turn is an ordered list of
  - text segments (Markdown, the existing sanitize → KaTeX pipeline, moved to `MarkdownMessage.jsx`),
  - **step cards** (`ToolStep.jsx`): status dot (running / ok / warn / failed), the tool's label,
    duration; expandable to arguments, the data sent to the model, and the checks list; an
    **Open page** pill when the tool names a page,
  - the **Outcome chip** (achieved / partial / not achieved + note) under the final answer,
  - a one-line footer: "3 analyses · 1 warning".
- **Empty state**: "What should I check?" with action-shaped suggestions ("Which site is most
  anharmonic, and does it have a preferred displacement direction?", "Compute the Se–Ta–Se bond
  angles and compare them with an ideal octahedron", "Has every fit term converged?"), plus two
  disclosures: *Tools the Copilot can run* (name + description) and *Context sent to the model*.
- Stop aborts the turn. Errors (HTTP, CORS) render as today, inside the transcript.

### 2.7 Host wiring (`App.jsx`, `copilotTools.js`, `useSiteCloud.js`)

- `useSiteCloud.js`: extract the body of `requestPca` into an exported plain function
  `requestAnalysis(kind, params, { rmc6fText, directory })`. The hook calls it. This is a
  behaviour-preserving refactor; `useSiteCloudReload` / `useSiteCloudStatic` tests must stay green.
- `copilotTools.js` (host): `useCopilotTools({ directory, localRun, structure, enabled })` returns
  the tool array. It reads the `.rmc6f` text lazily, on the first tool run (cached per file object),
  so opening the page costs nothing. The summarizers (`summarizeSites`, `summarizeOrientation`,
  `summarizeBondAngles`, `summarizeSymmetry`) and their checks are exported pure functions.
  `structure` comes from the Dashboard's published run (`assistantRun.structure`).
- `App.jsx`: tab label "AI Copilot", page key `copilot`, `IssueScope page="copilot"`, `tools` and
  `onOpenPage={handlePageChange}` props. The Dashboard props `wantAssistantData` /
  `onRunContextChange` keep their names: they feed the run context, which the Copilot still uses.

### 2.8 Files

New: `llm/copilot/{agentLoop,toolCalls,tools,runTools,outcome}.js`, `llm/prompts/copilot.js`,
`llm/useCopilot.js`, `llm/components/{CopilotPage,CopilotView,ToolStep,MarkdownMessage}.jsx`,
`llm/components/Copilot.css` (renamed from `AssistantPanel.css` + the step/outcome styles),
`src/copilotTools.js`, tests (§4).

Removed (replaced): `ChatView.jsx`, `AssistantPage.jsx`, `AssistantPanel.jsx` (an unrendered
dashboard card), `AssistantWorkspace.jsx`, `useStreamedReply.js`, `prompts/templates.js`'s
`buildChatMessages` (the watchdog builders stay). `AssistantPage.css` → `CopilotPage.css`.

`llm/index.js` exports `CopilotPage` and `WatchdogBadge`.

## 3. Privacy and safety

- Local providers: unchanged. Everything stays on the machine.
- Cloud providers: tool results are run-derived summaries of the same kind as the context (site
  statistics, angle modes, significances), but more of them. The drawer warning says so. The tools
  switch turns them off.
- The model's arguments are untrusted input. They pass the schema check, then the engines' own
  guards, with work budgets (`APP_MAX_ANGLES`, orientation `frequency` ≤ 64, smoothing ≤ 64). No tool
  writes, downloads, or navigates. The loop's round / call / time limits bound cost.
- Model text renders through the existing sanitizer (no raw `<img>`, no scripts). Tool data is shown
  as `JSON.stringify` text in a `<pre>`, never as HTML.

## 4. Tests (vitest)

- `client.test.js` (+): `tools` / `tool_choice` only in the body when given; SSE tool-call
  fragments split mid-arguments across reads; `finishReason`.
- `toolCalls.test.js`: assembly by index / id / order, Gemini whole-call chunks, argument parsing,
  schema validation (types, enum, min / max, required, unknown keys).
- `agentLoop.test.js` (scripted `callModel`): call → result → final with Outcome; unknown tool;
  invalid args fed back, then corrected; duplicate call cached; a thrown tool error; > 4 calls in
  a round; MAX_ROUNDS forces a final `tool_choice: 'none'` call; abort mid-tool; tools-unsupported
  400 → retry without tools; result budget trimming.
- `outcome.test.js`: parsing variants, missing line.
- `runTools.test.js`: `read_run_data` returns untrimmed sections and fails clearly when absent;
  `assess_convergence` on synthetic series.
- `src/__tests__/copilotTools.test.js`: the **real engines** on the committed demo run
  (`handlePcaMessage` with `public/demo` text, `structureFromRmc6f` + `describeSymmetry`) →
  summarizers → expected fields and checks (e.g. the Detected SG card's F-43m at 0.2 Å). A
  zero-angle window must produce the `fail` check.
- `copilotView.test.jsx`: the page end to end with a mocked `fetch` SSE server. A tool call
  renders a step card, its result reaches the second request's `messages`, and the Outcome chip
  renders.
- Existing suites stay green (`npm test`, `npx eslint src`, `npm run build`, Python suite untouched).

## 5. Verification in the browser

No local model server runs on this machine, so a scripted mock OpenAI-compatible server
(`/v1/models`, streaming `/v1/chat/completions` that issues a tool call and then a final answer
with an Outcome line) runs from the session scratchpad. Static-mode dev server, the Demo run, the
Copilot page: check the step card, the expand / collapse, Open page, the Outcome chip, Stop, the
tools-off switch, and the phone layout.

## 6. Docs

`src/llm/README.md` (rewritten for the Copilot: loop, tool contract, budgets, privacy),
`AGENTS.md` (architecture map, import-boundary paragraph: tools as props), `docs/algorithms/
ai-assistant.md` (a new *Copilot tool loop* section: what each tool sends, budgets, checks; data
flow table rows), `docs/ALGORITHMS.md` / `notation.md` references, `README.md` / `QuickStart.md`
(tab name, what the Copilot can run), `docs/CHANGELOG.md` (0.6.0 entry).

## 7. Build order

1. Provider: tool-call streaming + request body (+ tests).
2. Core: `toolCalls.js`, `tools.js`, `outcome.js`, `agentLoop.js`, `prompts/copilot.js` (+ tests).
3. In-module tools: `runTools.js` (+ tests).
4. Host: `requestAnalysis` refactor, `copilotTools.js` (+ demo-run tests).
5. UI: `useCopilot`, components, CSS, `App.jsx` wiring; remove the replaced files (+ view test).
6. Docs.
7. Verify: eslint, vitest, build, the browser pass with the mock server.

## 8. Risks

- **Small local models call tools badly** (wrong argument names, loops). Mitigations: few tools with
  plain snake_case arguments, schema errors fed back as text, duplicate-call caching, round limits,
  and the no-tools fallback.
- **Streaming dialect drift** across Ollama / LM Studio / OpenAI / Gemini / Anthropic compat layers.
  The assembler is tolerant (index, id or order) and unit-tested on each observed shape. Live
  provider checks need the user's keys and servers: a manual step, listed under Not done.
- **Long engine runs** (bond angles on a large box) are bounded by `APP_MAX_ANGLES` and the tool
  timeout. The worker is shared with the pages, so a Copilot request can queue behind a page
  request. That is acceptable, since the page shows its own spinner.

## 9. Phase 2 (not in this change)

- Open page **with parameters**: Bond Geometry pre-filled with the step's triplet and windows; PCA /
  Displacement Directions with the step's site selected (an App-level `pageRequest` channel).
- More tools: KDE slice statistics, partial g(r) re-binning, CIF export behind an explicit user
  confirmation, Auto StoG.
- Persisted conversations; exporting a turn (steps + data) as Markdown.

## 10. Review (2026-10-09) and what changed

An independent review of this plan, run against the code and the providers' documentation, raised
the findings below. Each line gives the change made, or why none was made.

| # | Finding | Disposition |
|---|---|---|
| 1 | Gemini 3 puts a `thought_signature` in each call's `extra_content`. Dropping it makes round 2 fail; replayed calls must not carry `index`. | **Fixed.** `streamChat` yields `extraContent`, the assembler keeps it, and the loop echoes it as `extra_content`. Echoed calls carry no `index`. Covered by tests on the Gemini whole-call shape. |
| 2 | Ollama ignores `tool_choice`, so the `'none'` round can come back as tool calls only, leaving an empty answer. | **Fixed.** An empty final round closes with a fixed "analysis budget used" answer and `Outcome: not achieved`. An empty non-final reply is marked `stoppedBy: 'empty'`. |
| 3 | Token budget vs small local context windows: Ollama's `/v1` cannot set `num_ctx`. | **Fixed.** A 12,000-character cap on all of a question's tool messages compacts older results to their checks (the newest round is never compacted). The setup guide, QuickStart and README suggest `OLLAMA_CONTEXT_LENGTH=16384` and a large LM Studio context. |
| 4a | Anthropic: Claude 4.7+ reportedly return 400 for `temperature` < 1. | **Not changed.** This predates the Copilot (the Assistant sent 0.2 too) and could not be checked without a key. Listed under Not done. |
| 4b | The tools-unsupported regex would also match a schema 400. | **Fixed.** It now matches only "does not support …" / "… is not supported" for tools, tool calling, tool use or function calling. A test pins that a schema 400 surfaces. |
| 5a | `requestAnalysis` had no `signal`; the shared serial worker keeps running a stopped or timed-out job. | **Fixed.** Each tool call gets its own `AbortSignal` (Stop or timeout). Browser-loaded runs use a **Copilot-only worker** that is terminated on abort; Flask requests pass the signal to axios. |
| 5b | `{ id, kind, text, ...params }` lets params override the transport keys. | **Fixed.** Params are spread first, in both transports. The host maps model arguments through an explicit per-tool mapping, never by spreading. |
| 6 | The page-key rename was incomplete (`visitedPages.assistant`, the test mock). | **Fixed.** `visitedPages.copilot`, `wantAssistantData={visitedPages.copilot}`, and the `appLiveData.test.jsx` mock is `CopilotPage`. |
| 7 | Whitespace-only content before a tool call is sent as a text block, which Anthropic rejects. | **Fixed** with `assistantContent()`: `text.trim() ? text : null`, tested. Whether Ollama and llama.cpp accept `null` is unverified (live check, Not done). |
| 8 | The run changing under the Copilot. | **Fixed** for run switches: `runKey` = directory + run id; a change stops the turn and restarts the conversation (tested). Tools and context are captured when the question is asked; `tools` is in no effect's dependencies. **Not changed:** a Flask Live Data save mid-question is not a run switch, so later analyses read the new file. Documented in ai-assistant.md Step 17. |
| 9 | A stopped or failed turn left two user messages in a row, or an empty assistant message, in memory. | **Fixed.** `conversationMemory()` pairs each question with its answer and drops unanswered ones. A stopped turn's partial text is kept, marked "(Interrupted before the answer was finished.)". |
| 10 | Tool details: site_ref/element exclusivity, `geometry:false`, the cell for `[u v w]`, peak ties, tolerance bounds, noisy `sinCorrected`, trimming order. | **Fixed:** one of `site_ref`/`element` enforced in `run()` (site_ref wins when both are given); `geometry:false`; the cell from the payload or the Dashboard's structure; a tied-peak note; `tolerance_A` ∈ [0.001, 1]; modes from a 5-bin running mean of the counts, not `sinCorrected`; `list_sites` ranked by \|kurtosis\| with an element filter. |
| 11 | Scope: `read_run_data` is not "untrimmed"; "every fit term" is unanswerable; no per-turn time limit. | **Fixed:** the description now says the section is as built, before the budget trim; the suggestion was reworded; `maxTurnMs` = 10 min stops new rounds. `assess_convergence` is kept: it adds the classification and the tail slopes the context lacks. |
| 12 | Accessibility and the kit: `aria-expanded`, status as text, a live region, the kit class list. | **Fixed:** the step head is a button with `aria-expanded`; status is text ("Done · 0.4 s"); a polite live region announces step and verdict changes; the class list in AGENTS.md and the llm README drops the three unused header classes. The Outcome chip and the switch use `llm-*` classes and a plain checkbox. |
| 13 | Tests: `markdown.test.jsx` copies the pipeline; `templates.test.js` tests a removed builder; missing cases. | **Fixed.** `markdown.test.jsx` renders `MarkdownMessage`; the chat-prompt cases moved to `copilotPrompt.test.js`; new cases cover `extra_content`, whitespace → null, an empty final round, compaction, abort of the tool signal on timeout, and a run change. |
| 14 | Privacy: tools default-on for cloud presets; clickable links in model output. | **Not changed.** The analysis results are summaries of the same kind as the context, which already goes to the cloud provider. The warning names them and the switch turns them off. Link rendering is the Assistant's existing behaviour; the model's input is the user's own run files. |

The review also confirmed:

- the `requestPca` extraction is behaviour-identical;
- `basis` exists only from the browser parser;
- `describeSymmetry` takes ~0.1 s on the demo;
- every payload field the summarizers read exists;
- `WatchdogBadge` *is* rendered (Dashboard), so it stays.

Claims it could not confirm, about null content and `tool_choice: 'none'` on Anthropic, Gemini and
LM Studio, are why the loop tolerates either behaviour (findings 2 and 7).

## 11. Change: a chatbox over every page (2026-10-09)

After the first build, the maintainer asked for the Copilot to be **a chatbox, not a single page**,
because it operates the analyses: the conversation should stay in view beside the page that shows
the analysis. What changed:

- **No Copilot tab.** `CopilotChat` (llm module) is rendered once by `App.jsx`, over every page:
  - a labelled launcher in the bottom-right corner opens a floating panel (28 rem × 44 rem;
    *Larger* gives 48 rem × the window height);
  - ✕ or Esc closes it, the panel stays mounted (the conversation and a running answer survive),
    and focus returns to the launcher;
  - nothing is fetched before the first open, and the Dashboard's Copilot-only work waits for it
    (`copilotOpened`);
  - on a phone the panel is a full-height sheet, and the launcher an icon-only round button.
- **Show in \<page\>** replaces §1's "Open page" and absorbs §9's first item:
  - a tool may declare `pageParams(args)`, the loop attaches them to the step, and the button calls
    `onOpenPage(page, params)`;
  - the host keeps a `pageRequest` channel, and Bond Geometry (fill + compute) and Displacement
    Directions (site + resolution) apply a `request` once per click;
  - on a phone the sheet closes after Show in.
- **Still the user's action.** §1's non-goal stands: the Copilot does not navigate or change pages
  by itself. With the chatbox, the reason given in §1 (the transcript would be hidden) no longer
  applies. Whether the model may drive pages on its own is a product decision left to the
  maintainer.

## 12. Change: docked beside the pages, pointing at the page (2026-10-10)

The maintainer asked for the highlight feature of their other two agents (MATERIA, NEBULA3D) and
for a Copilot that no longer blocks the panels. That decided §11's open question: the model may
now change **what is on screen** — never data, an analysis input or a file.

- **Docked, not floating.** The panel is a column of the app grid beside the pages (≥ 1080 × 541),
  else a bottom sheet in its own row; the pages' width tiers became container queries, so they
  reflow in the room left. It resizes from its edge, collapses to a rail with an unread dot, and
  remembers its width and open state.
- **View tools.** `highlight_cards` (the cards breathe a blue edge, a plot gets marks),
  `orient_3d_view` (look down [uvw], centre, ring or select sites, reset) and `current_view`.
  Answers cite cards as `card:` chips. Show in stays the user's action, because it loads a page's
  analysis inputs.
- The details are in the changelog and the module README.
