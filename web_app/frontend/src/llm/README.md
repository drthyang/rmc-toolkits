# AI Assistant (experimental LLM module)

An experimental layer that connects the RMCProfile Run Monitor to a **local LLM**
(Ollama or LM Studio) running on the user's own machine. It demonstrates the
standard LLM-in-a-data-pipeline pattern end to end, in the browser, with no
server and no API keys:

```
Dashboard state ──props──▶ context builder ──▶ prompt templates ──▶ provider client ──▶ streaming UI
(plotFiles, rValueFile,     context/            prompts/             provider/           components/
 structure, symmetry)       runContext.js       templates.js         client.js (SSE)
```

Features: a one-click **run summary**, **chat** about the loaded run, a
**Markdown report** generator, and a live **convergence watchdog** badge.

## Privacy contract

The hosted app's promise is that run data never leaves the device. This module
keeps that promise: the browser talks **directly** to `http://localhost:*` —
the only party that ever sees the (already summarized) run context is the LLM
server the user runs themselves. Nothing is sent until the user expands the
panel and clicks an action.

## How each seam works (reading order)

1. **`context/runContext.js`** — the heart of the module. `buildRunContext()`
   distills the dashboard's parsed state into a compact JSON object: cell and
   angles recomputed from the lattice vectors, per-dataset Rwp values, and the
   ln(χ²) convergence history **downsampled to ≤48 points** with summary
   statistics (first/last/min/max/recent slope) computed on the *full* series
   first. `contextToJson()` enforces a ~3,000-character budget so small local
   models never overflow. Every field is optional — a run with no `.rmc6f` or
   no logs simply omits those sections.
2. **`prompts/system.js` + `prompts/templates.js`** — one shared system prompt
   establishes the domain (what Rwp means, that the history is ln(χ²), "quote
   the numbers, don't guess"); per-feature builders produce the `messages`
   array. The context always travels inside a fenced JSON block so the model
   can tell data from instructions. The watchdog prompt demands a rigid
   `STATUS: … — …` one-liner that `parseWatchdogReply()` parses with a regex.
3. **`provider/client.js`** — a hand-rolled OpenAI-compatible client
   (`/models`, `/chat/completions`) with an SSE stream parser over `fetch` +
   `ReadableStream`. Works identically against Ollama's `/v1` and LM Studio.
   `checkConnection()` translates failures into actionable hints — a bare
   `TypeError` from `fetch` means "server down **or** CORS blocked", and the
   user cannot tell those apart from the page, so the hint names both fixes.
4. **`watchdog/heuristics.js` + `watchdog/useWatchdog.js`** — the layering that
   keeps the LLM honest: pure slope heuristics over the ln(χ²) history are the
   **source of truth** for the badge status; the LLM only writes the tooltip
   note. The hook observes the `rValueFile` prop (which the existing 3-second
   Live Data poll already refreshes — no new timers) keyed off *content* stats,
   not object identity, and calls the model only when the status flips or the
   throttle interval elapses with a significant change.
5. **`report/buildReport.js`** — deterministic Markdown tables built from the
   same context object, with the LLM narrative grafted in under a clearly
   labeled "AI assessment" heading with a verification disclaimer. The tables
   work with no model connected at all.
6. **`components/`** — the collapsible dashboard card (`AssistantPanel`), the
   settings drawer (`ConnectionSettings`), the three tab views, and the
   `WatchdogBadge` rendered on the R-value card. Replies render as plain
   pre-wrapped text (no HTML injection surface). The Summary tab includes a
   "Context sent to the model" inspector — open it to see exactly what the
   model sees.

## Local LLM setup

- **Ollama** — <https://ollama.com>. Pull a model (`ollama pull llama3.2`),
  then start the server with this site allowed as an origin:

  ```bash
  OLLAMA_ORIGINS="https://drthyang.github.io" ollama serve   # or OLLAMA_ORIGINS="*"
  ```

- **LM Studio** — <https://lmstudio.ai>. Load a model, open the Developer tab,
  enable **CORS** in the server settings, and start the local server
  (default `http://localhost:1234/v1`).

In the app: open the **AI Assistant** card → pick the preset → **Test
connection** → choose a model. Settings persist in `localStorage`
(`rmc-llm-settings-v1`).

Notes: HTTPS pages are allowed to call `http://localhost` in Chrome, Edge, and
Firefox (trustworthy-origin exemption); Safari may block it. Future Chrome
private-network-access rules may add restrictions — if the hosted app cannot
connect, running the app locally (`npm run dev` or the Flask/Docker mode)
always works.

## Import boundary & extraction recipe

This module is designed to be lifted into its own repository:

- Run data enters **only as props** to `AssistantPanel` / `WatchdogBadge`
  (`runName`, `plotFiles`, `rValueFile`, `structure`, `symmetry`, `liveData`).
  Nothing in `src/llm/` imports `browserData.js`, `symmetryModel.js`, or any
  host component; cell math is deliberately duplicated from `ModelSummary.jsx`.
- The **single** allowed host import is `downloadBlob` / `sanitizeFilename`
  from `../figureExport.js` (≈20 lines, trivially inlined).
- The host imports only `src/llm/index.js`.

To extract: copy `src/llm/`, inline those two helpers into
`report/buildReport.js`, and feed the six props from your own data source.

## Tests

`__tests__/` covers everything that doesn't need a live model: downsampling
invariants, context assembly (including missing-structure/missing-log runs),
the character-budget guard, convergence heuristics on synthetic series, SSE
parsing against mocked fetch streams (chunk splits mid-line, `[DONE]`,
malformed lines), watchdog reply parsing, and report snapshots. Run with
`npm test` (vitest).
