# AI Assistant (experimental LLM module)

An experimental layer that connects the RMCProfile Run Monitor to an
**OpenAI-compatible** LLM — a **local** server the user runs (Ollama or LM
Studio) or a **cloud** provider (OpenAI, Gemini) via an API key. It demonstrates
the standard LLM-in-a-data-pipeline pattern end to end, in the browser, with no
backend:

```
Dashboard state ──props──▶ context builder ──▶ prompt templates ──▶ provider client ──▶ chat UI
(plotFiles, rValueFile,     context/            prompts/             provider/           components/
 structure, symmetry)       runContext.js       templates.js         client.js (SSE)
```

Features: a **chat** about the loaded run (ask it to summarize, compare
datasets, judge convergence, spot local distortion, …) — with a collapsible
**"Thinking"** panel for reasoning models — and a live **convergence watchdog**
badge on the dashboard.

## Privacy

Local providers keep the promise that run data never leaves the device: the
browser talks **directly** to `http://localhost:*`, and the only party that sees
the (already summarized) run context is the server the user runs themselves.
**Cloud** providers necessarily send that context to their API — this is
opt-in, flagged in the UI with a data-leaves-your-device warning, and local
stays the default. Nothing is sent until the user acts. API keys are stored in
`localStorage` only.

## How each seam works (reading order)

1. **`context/runContext.js` + `context/pairCorrelations.js`** — the heart of
   the module. `buildRunContext()` distills the dashboard's parsed state into a
   compact JSON object: cell/angles, a **symmetry block** (space group, the
   symmetry-vs-tolerance ladder, and Wyckoff-orbit sites ranked by rms
   displacement — the local-distortion evidence), **pair correlations** (each
   partial-PDF pair's first g(r) peaks next to the average structure's
   nearest-neighbour distance), a **configuration_optimization block** (move
   counters from the `.rmc6f` header with acceptance ratio and accepted moves per
   atom — the sampling-sufficiency gauge), **run settings** (the stem-matched RMCProfile
   `.dat` control file: pair-labeled minimum-distance constraints, move sizes,
   fitted datasets, flags), per-dataset Rwp values, and the ln(χ²)
   convergence history **downsampled to ≤48 points** with summary statistics
   computed on the *full* series first. `contextToJson()` enforces a
   ~4,500-character budget with an explicit trim order (extra peaks → ladder
   rungs → history → low-displacement sites → datasets, each recorded with an
   `*_omitted` count) so small local models never overflow. Every field is
   optional — a run with no `.rmc6f`, no partials, or no logs simply omits
   those sections.
2. **`prompts/system.js` + `prompts/templates.js`** — one shared system prompt
   establishes the domain (what Rwp means, that the history is ln(χ²), "quote
   the numbers, don't guess"). `buildChatMessages()` injects the context ahead
   of the conversation; `buildWatchdogMessages()` demands a rigid
   `STATUS: … — …` one-liner that `parseWatchdogReply()` parses with a regex.
   The context always travels inside a fenced JSON block so the model can tell
   data from instructions.
3. **`provider/client.js`** — a hand-rolled OpenAI-compatible client (`/models`,
   `/chat/completions`) with an SSE stream parser over `fetch` +
   `ReadableStream`. Works identically against Ollama's `/v1`, LM Studio, and
   cloud endpoints; a Bearer API key is sent only when one is set.
   `checkConnection()` translates failures into actionable hints — a bare
   `TypeError` means "server down **or** CORS blocked", and 401/429 name the
   key/quota problems. `streamChat()` yields structured `{ content } |
   { reasoning }` chunks, so reasoning models' chain-of-thought is surfaced
   separately from the answer.
4. **`watchdog/heuristics.js` + `watchdog/useWatchdog.js`** — pure slope
   heuristics over the ln(χ²) history are the **source of truth** for the badge
   status; the LLM only writes the tooltip note. The hook observes the
   `rValueFile` prop (refreshed by the existing 3-second Live Data poll — no new
   timers) keyed off *content* stats, not object identity, and calls the model
   only when the status flips or the throttle interval elapses with a
   significant change.
5. **`useAssistant.js` + `components/`** — `useAssistant()` holds the shared
   state (settings, the connection probe + auto-connect, and the run context).
   The page (`AssistantPage`) composes `AssistantConnectionBar` (status
   indicator + model switcher + settings gear), the `ConnectionSettings` drawer
   (local + cloud providers, API key, warnings), and `ChatView` (a modern
   composer with message bubbles and the app's wave avatar; assistant replies
   render **Markdown** — GitHub-flavored tables, lists, code — via
   `react-markdown`, which builds a React tree and disallows raw HTML, so there
   is still no injection surface; reasoning models get a collapsible "Thinking"
   panel that times the thinking phase). `WatchdogBadge` renders on the
   dashboard's R-value card. `AssistantPanel` wraps the same `ChatView` workspace as a collapsible
   dashboard card, kept as a self-contained extractable surface.

## Provider setup

- **Ollama** (local) — <https://ollama.com>. Pull a model
  (`ollama pull llama3.2`), then start the server with this site allowed as an
  origin:

  ```bash
  OLLAMA_ORIGINS="https://drthyang.github.io" ollama serve   # or OLLAMA_ORIGINS="*"
  ```

- **LM Studio** (local) — <https://lmstudio.ai>. Load a model, enable **CORS**
  in the Developer/server settings, and start the local server
  (default `http://localhost:1234/v1`).
- **OpenAI / Gemini** (cloud) — pick the provider and paste an API key. No
  server to run, but your run context is sent to that provider, and some
  providers or networks may block direct browser calls (CORS).

In the app: open **AI Assistant** → pick a provider → **Test** → choose a model.
Settings persist in `localStorage` (`rmc-llm-settings-v1`).

Notes: HTTPS pages are allowed to call `http://localhost` in Chrome, Edge, and
Firefox (trustworthy-origin exemption); Safari may block it. If a connection
cannot be made, running the app locally (`npm run dev` or the Flask/Docker mode)
always works for local providers.

## Import boundary & extraction recipe

This module is designed to be lifted into its own repository:

- Run data enters **only as props** to `AssistantPage` / `AssistantPanel` /
  `WatchdogBadge` (`runName`, `plotFiles`, `rValueFile`, `structure`, `symmetry`, `runSettings`,
  `liveData`). Nothing in `src/llm/` imports `browserData.js`,
  `symmetryModel.js`, or any host component — the module is fully self-contained;
  cell math is deliberately duplicated from `ModelSummary.jsx`.
- The host imports only `src/llm/index.js`.

To extract: copy `src/llm/` and feed the seven props from your own data source.

## Tests

`__tests__/` covers everything that doesn't need a live model: downsampling
invariants, context assembly (including missing-structure/missing-log runs), the
character-budget guard, convergence heuristics on synthetic series, SSE parsing
against mocked fetch streams (chunk splits mid-line, `[DONE]`, malformed lines),
and watchdog reply parsing. Run with `npm test` (vitest).
