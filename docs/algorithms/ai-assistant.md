# AI Copilot — algorithm reference

> Part of the [Algorithms and Math Reference](../ALGORITHMS.md). Every step is anchored to the source; if this document and the code disagree, **the code wins**.

The page covers three things:

- what the Copilot computes from the loaded run before any model is called;
- which analyses it may run while it answers, and the deterministic checks on each result;
- exactly what data leaves the device.

The **AI Copilot** replaced the chat-only **AI Assistant** on 2026-10-09. The run context (Steps
1–12), the watchdog (Steps 14–15) and the transport (Step 16) are the Assistant's, unchanged. Step 13
describes the Copilot's messages and Step 17 its tool loop. Section anchors keep their old names.

## Contents

- [AI Assistant — run context construction, pair correlations, and convergence heuristics](#ai-assistant--run-context-construction-pair-correlations-and-convergence-heuristics)
  - [What this page shows](#what-this-page-shows)
  - [Data flow — precisely what leaves the device](#data-flow--precisely-what-leaves-the-device)
  - [Step 1 — Inputs, and when the context is built](#step-1--inputs-and-when-the-context-is-built)
  - [Step 2 — Cell and composition block (`context.structure`)](#step-2--cell-and-composition-block-contextstructure)
  - [Step 3 — Sampling counters (`context.configuration_optimization`)](#step-3--sampling-counters-contextconfiguration_optimization)
  - [Step 4 — Run-control settings (`context.run_settings`)](#step-4--run-control-settings-contextrun_settings)
  - [Step 5 — Symmetry and per-orbit displacement (`context.symmetry`)](#step-5--symmetry-and-per-orbit-displacement-contextsymmetry)
  - [Step 6 — Per-site PCA displacement summary (`context.pca_displacements`)](#step-6--per-site-pca-displacement-summary-contextpca_displacements)
  - [Step 7 — Average-structure nearest-neighbour distances](#step-7--average-structure-nearest-neighbour-distances)
  - [Step 8 — $g(r)$ peak extraction](#step-8--gr-peak-extraction)
  - [Step 9 — Assembling `context.pair_correlations`](#step-9--assembling-contextpair_correlations)
  - [Step 10 — Per-dataset residuals (`context.datasets`)](#step-10--per-dataset-residuals-contextdatasets)
  - [Step 11 — Convergence series (`context.convergence`)](#step-11--convergence-series-contextconvergence)
  - [Step 12 — The character budget (`contextToJson`)](#step-12--the-character-budget-contexttojson)
  - [Step 13 — Message assembly (what the request body actually contains)](#step-13--message-assembly-what-the-request-body-actually-contains)
  - [Step 14 — Convergence heuristics (the statistics)](#step-14--convergence-heuristics-the-statistics)
  - [Step 15 — Watchdog scheduling and badge semantics](#step-15--watchdog-scheduling-and-badge-semantics)
  - [Step 16 — Transport](#step-16--transport)
  - [Step 17 — The Copilot tool loop (analyses, checks, budgets, verdict)](#step-17--the-copilot-tool-loop-analyses-checks-budgets-verdict)
  - [Step 18 — Cluster distortion against a shuffle null (`analyze_clusters`)](#step-18--cluster-distortion-against-a-shuffle-null-analyze_clusters)
  - [Step 19 — Displacement modulations: order that enlarges the cell (`analyze_modulations`)](#step-19--displacement-modulations-order-that-enlarges-the-cell-analyze_modulations)
  - [Runtime-mode differences (static vs Flask) — read this before quoting the assistant](#runtime-mode-differences-static-vs-flask--read-this-before-quoting-the-assistant)
  - [Parameters and defaults](#parameters-and-defaults)
  - [Caveats / what this is not](#caveats--what-this-is-not)

---

## AI Assistant — run context construction, pair correlations, and convergence heuristics

### What this page shows

The **AI Copilot** is a panel docked beside every page (a bottom sheet on portrait tablets and phones,
collapsed to a slim rail) that holds a conversation about the run currently loaded in the
Workbench, and points at the page while it explains (Step 17, *View tools*).
There is also an optional **convergence watchdog** chip on the Dashboard's R-value card. Neither
feature sends your run files anywhere.

The app first *computes a compact summary* of the run: a JSON object targeted at ~4,500 characters
(a best-effort budget, not a hard cap; see Step 12). It hands that summary and your question to the
OpenAI-compatible chat endpoint you configured. When the model asks for an analysis, the Copilot runs
it in the browser (or through the Flask engines in server mode) and sends back a compact result with
deterministic checks (Step 17). By default the endpoint is a server on your own machine
(`http://localhost:11434/v1`, i.e. Ollama).

Everything documented below is browser-side JavaScript in
[`web_app/frontend/src/llm/`](../../web_app/frontend/src/llm/). The module is deliberately self-contained:
it receives run data **only as React props** and imports nothing from the rest of the app (see the
import-boundary rule in [AGENTS.md](../../AGENTS.md) and [`src/llm/README.md`](../../web_app/frontend/src/llm/README.md)).
There is no Python counterpart: none of this exists in `rmc_toolkits/`. The Flask backend is never
involved in a model request. It is involved only when a Copilot analysis of a typed server folder
runs the same `/api/pca/*` and `/api/triplets` routes the pages use. The host builds those tools in
[`copilotTools.js`](../../web_app/frontend/src/copilotTools.js) and passes them in as a prop.

Two things in this section are computed statistics that a scientist should be able to audit:

1. the **pair-correlation cue extraction** — first-peak positions and approximate widths of each partial
   $g(r)$, set against the nearest-neighbour distance of the same element pair in the *average*
   structure; and
2. the **convergence watchdog heuristics** — a tail least-squares slope on the $\ln\chi$ history with
   fixed thresholds, which is the *source of truth* for the badge (the model only writes the tooltip
   sentence).

---

### Data flow — precisely what leaves the device

This subsection is normative: it is derived from the only two `fetch()` calls in the whole module, both
in [`llm/provider/client.js`](../../web_app/frontend/src/llm/provider/client.js).

| # | When | Request | Run-derived payload |
|---|------|---------|---------------------|
| 1 | AI Copilot dock first opened (automatic, once per base URL; at load when it was left open), or **Test** pressed | `GET {baseUrl}/models` → `listModels()` | **None.** Only the API key, if one is set. |
| 2 | You send a question | `POST {baseUrl}/chat/completions` (`stream: true`) → `postChat()` via `streamChat()` | The summarized context JSON (target ≤ ~4,500 chars, best-effort; Step 12), the system prompt, the last 8 memory messages (questions and final answers only) and your question. With analyses on, the tool definitions (names, descriptions, argument schemas). |
| 2a | The model asked for analyses (repeats ≤ 6 times per question, plus one last request) | the same `POST`, same body grown by the exchange | Row 2's body, plus the model's own tool calls echoed back and one `tool` message per call: the analysis summary and its checks. Each is ≤ 3,000 chars and ≤ 12,000 per question; older results are compacted to their checks. Step 17 lists what each analysis sends. |
| 3 | Watchdog enabled *and* a model configured | `POST {baseUrl}/chat/completions` (`stream: false`) → `completeChat()` | Five numbers only: `{n_steps, first, last, min, recent_window_delta}` plus two status words. |

**Row 1 fires without you asking.** The dock's panel (`CopilotChat.jsx`) mounts on first open (or at
load, when the dock was left open: its open state persists) and
calls `useAssistant({… enabled: true})`, and
`useAssistant.js` runs `probe(false)` from an effect guarded only by `autoTestedRef` (one probe per
distinct `baseUrl`). Merely opening the Copilot therefore issues a `GET /models` to whatever
endpoint is configured — **including a cloud endpoint, with your API key attached**, before you have
typed anything. Switching the selected model does not re-probe; changing the base URL does. On success,
if the currently configured model is not in the returned list, the app silently writes `models[0]` into
`localStorage` as your selected model (`saveSettings({ model: result.models[0] })`).

`baseUrl` comes from `settings.baseUrl`, persisted in `localStorage` under the key
`rmc-llm-settings-v1` ([`llm/settings.js`](../../web_app/frontend/src/llm/settings.js)). Default:
`http://localhost:11434/v1`.

**What is never transmitted.** No file is ever uploaded. The context builder reads only already-parsed
in-memory summaries; it never touches `File` objects, file text, or the filesystem. Specifically absent
from the context object built by `buildRunContext()`:

- the per-atom rows of the `.rmc6f` configuration — no atom coordinates from the `Atoms:` section, and
  not the sampled point cloud the 3D view uses (`structure.points` is not read);
- the absolute path of the structure file (`structure.source` is not read) or of any plot file;
- the raw `.log`, `.csv`, `.dat`, or STOG file contents;
- the full $\chi$ history (only ≤48 rounded points survive), the full $g(r)$ curves (only ≤2 peaks per
  pair), and the full displacement clouds (only per-site scalars);
- in an analysis result (Step 17): never the engine payload itself. The KDE volume, the sphere-cell
  arrays and the raw angle list stay in the browser; only the summary the tool builds is sent.

**What *is* transmitted that identifies your work.** Be aware of these — they are file-derived strings
and numbers:

- `run` — the run folder name (static mode: the picked folder's name; Flask mode: the directory string
  you typed, which may be a relative path such as `data/5K_try1`). When it is empty, the literal string
  `"unnamed run"` is sent instead;
- `run_settings.title / material / phase / temperature` — verbatim from the run-control `.dat`;
- `run_settings.fitted_data[].file` — the `> FILENAME ::` value copied **verbatim** from the `.dat`
  (`parseRunSettings()` only `.trim()`s it; no basename is taken). It is usually a bare basename such as
  `scale_ft_rmc.fq`, but if the run-control file declares a path, that path is transmitted;
- `datasets[].title` — plot titles, several of which are derived from file names
  (`plotMetadataFromFile()` in [`browserData.js`](../../web_app/frontend/src/browserData.js) builds the
  `npdf`/`pdf_partials` title from the file-name stem);
- **representative fractional coordinates**: `symmetry.sites[].frac` carries, for each reported Wyckoff
  orbit (≤12 of them), the orbit representative's three fractional coordinates in the conventional cell
  at 3 s.f. These are average/circular-mean basis positions, not instantaneous atom positions, but they
  are atomic coordinates and they do leave the device (see Step 5);
- element symbols, counts, cell metrics, space group, Wyckoff labels and site symmetries, and the
  numbers listed in Steps 2–11;
- **with analyses on, more of the same kind** (Step 17), and only for the analyses the model runs:
  - per-site displacement statistics for *every* site, not ≤ 12;
  - displacement-direction peaks as Cartesian unit vectors and `[u v w]`;
  - bond-angle modes and coarse histograms;
  - `detect_symmetry` orbits with representative fractional coordinates at 4 decimals, for *every*
    orbit at any tolerance the model picks.

**Local vs cloud.** `PROVIDER_PRESETS` in
[`llm/provider/presets.js`](../../web_app/frontend/src/llm/provider/presets.js) marks Ollama
(`http://localhost:11434/v1`) and LM Studio (`http://localhost:1234/v1`) as `cloud: false`, and OpenAI,
Gemini, and Anthropic as `cloud: true`. `isLocalUrl()` classifies any base URL whose hostname is in
`{localhost, 127.0.0.1, [::1], ::1, 0.0.0.0}` as local. When the active provider is cloud, the settings
drawer renders a "Run data leaves your device — sent to *provider*, with the results of the analyses
the Copilot runs" banner (its tooltip points to a local provider), an API-key field and the *Let the
Copilot run analyses* switch (`copilotTools`, default on; off sends no tool definitions and runs no
analysis) ([`ConnectionSettings.jsx`](../../web_app/frontend/src/llm/components/ConnectionSettings.jsx)).
Cloud use is opt-in: the user must select the preset and paste a key — but note that once a cloud preset
is active, the row-1 probe above authenticates with that key on every fresh page mount.

**No rmc-toolkits server is ever contacted.** There is no telemetry, no proxy, no relay. On GitHub Pages
(static mode) the page is served from `drthyang.github.io` but the model traffic goes browser →
`localhost` (or browser → the cloud provider you chose) directly. The API key is stored only in
`localStorage` and sent only as an `Authorization: Bearer …` header, and only when non-empty
(`authHeaders()`). For `api.anthropic.com` two extra headers are added
(`anthropic-version: 2023-06-01`, `anthropic-dangerous-direct-browser-access: true`, in
`requestHeaders()`); no other provider-specific handling exists.

**You can inspect the exact payload.** There are two places to look.

- **The context, before your first message.**
  [`CopilotView.jsx`](../../web_app/frontend/src/llm/components/CopilotView.jsx) renders a
  `<details>` panel labeled *"Context sent to the model"*. It holds `contextToJson(context)`, the
  identical string that `contextBlock()` embeds in the request.
  - The panel lives in the **empty-state** branch (`hasConversation = turns.length > 0 || busy`).
    It disappears with the first question and returns only after *New conversation* or a reload.
  - A second panel there lists every analysis the Copilot may run, with the description the model
    reads.
- **Every analysis result, at any time.** Expanding a step card
  ([`ToolStep.jsx`](../../web_app/frontend/src/llm/components/ToolStep.jsx)) shows the arguments and
  the exact result object sent to the model (before the per-question compaction of Step 17). It also
re-serializes the context *as of now*: with Live Data on, the JSON shown in the panel is not necessarily
the JSON that accompanied an earlier answer.

---

### Step 1 — Inputs, and when the context is built

**Inputs.** `buildRunContext()` takes eight optional fields
([`llm/context/runContext.js`](../../web_app/frontend/src/llm/context/runContext.js)):

| Prop | Origin | Used for |
|------|--------|----------|
| `runName` | folder name / directory string | `context.run` (falls back to `'unnamed run'`) |
| `plotFiles` | `Dashboard.jsx` → `allPlotFiles`, the **unfiltered** dashboard plot list | `datasets`, `pair_correlations` |
| `rValueFile` | combined R-value log, built from the **visibility-filtered** list | `convergence` |
| `structure` | parsed `.rmc6f` summary | `structure`, `configuration_optimization`, per-site displacements, average-structure distances |
| `symmetry` | `describeSymmetry()` + `toleranceLadder()` | `symmetry` |
| `runSettings` | `parseRunSettings()` of `<stem>.dat` | `run_settings` |
| `pcaSites` | PCA Ellipsoid page's site table | `pca_displacements` |
| `liveData` | Live Data toggle | `live_mode` (boolean) |

**The visibility toggle is asymmetric.** `Dashboard.jsx` builds `allPlotFiles` (every file passing
`isDashboardPlotFile`, sorted by `comparePlotFiles`), then `plotFiles = allPlotFiles` minus
`hiddenPlotPaths`, then `rValueFiles` from `plotFiles`. The assistant is handed
`assistantRun.plotFiles = allPlotFiles` — the *unfiltered* list — but `rValueFile` from the filtered
one. Consequence: **hiding a plot on the Dashboard does not remove it from `datasets` or
`pair_correlations`, but hiding an R-value log silently removes it from `convergence`, from
`recent_slope_per_step`, and from the watchdog.**

**Gating.** `useAssistant({…, enabled})` builds the context in a `useMemo` only when `enabled` is true
([`llm/useAssistant.js`](../../web_app/frontend/src/llm/useAssistant.js)). The chatbox panel is mounted
only when the user first opens it (`copilotOpened` in
[`App.jsx`](../../web_app/frontend/src/App.jsx)), and `Dashboard.jsx` computes the symmetry finder, the
tolerance ladder, and the `.dat` parse *only* when `wantAssistantData` is set — so a Dashboard/KDE-only
session performs none of this work.

**Rounding convention.** A single helper is used throughout:

```js
roundSig = (v, digits = 3) => (Number.isFinite(v) && v !== 0 ? Number(v.toPrecision(digits)) : v)
```

i.e. round to `digits` significant figures; exact zero and non-finite values pass through unchanged.
`pairCorrelations.js` duplicates this helper verbatim to avoid a circular import. Note that
`watchdogStats()` in Step 14 does **not** use it — it calls `Number(x.toPrecision(n))` bare.

**Every section is optional.** Missing inputs simply omit their key; `buildRunContext({runName})` returns
`{run, live_mode}` and nothing else (asserted in
[`__tests__/runContext.test.js`](../../web_app/frontend/src/llm/__tests__/runContext.test.js)). The two
always-present keys are the only defaults in the builder: `run: runName || 'unnamed run'` and
`live_mode: Boolean(liveData)`.

---

### Step 2 — Cell and composition block (`context.structure`)

**Input:** `structure.latticeVectors` (3×3, rows = supercell lattice vectors $\mathbf{L}_i$ in Å),
`structure.supercell` $= (n_1,n_2,n_3)$, `structure.totalAtoms`, `structure.elementCounts`.
Returns `null` unless both `latticeVectors` and `supercell` are present.

**Math.** Conventional-cell edge lengths (Å, 4 significant figures):

$$a_i = \frac{\lVert \mathbf{L}_i \rVert}{\max(N_i,1)}, \qquad i = 1,2,3$$

The $\max(N_i,1)$ is a guard: a zero or missing supercell entry would otherwise divide by zero.

Cell angles (degrees, 4 s.f.), from the *supercell* vectors (scaling by $N_i$ does not change an angle):

$$\alpha = \angle(\mathbf{L}_2,\mathbf{L}_3),\quad \beta = \angle(\mathbf{L}_1,\mathbf{L}_3),\quad
\gamma = \angle(\mathbf{L}_1,\mathbf{L}_2), \qquad
\angle(\mathbf a,\mathbf b) = \frac{180}{\pi}\arccos\!\left(\frac{\mathbf a\cdot\mathbf b}{\max(\lVert\mathbf a\rVert\lVert\mathbf b\rVert,\,10^{-12})}\right)$$

The $\arccos$ argument is clamped to $[-1,1]$.

**Output:** `{cell_A: [a,b,c], angles_deg: [α,β,γ], supercell, total_atoms, elements}` where `elements`
is the per-element atom count of the whole supercell. Two details:

- `supercell` and `total_atoms` are **verbatim passthroughs** — no `roundSig`, no validation — as is
  `elements`, whose per-element counts are copied straight from `structure.elementCounts`. Only
  `cell_A` and `angles_deg` are rounded (4 s.f.).
- `elements` is attached only when `structure.elementCounts` exists **and** has at least one key;
  otherwise the key is absent.

**Code:** `runContext.js` → `structureContext()`, `vectorLength()`, `angleBetween()`. This math is
deliberately duplicated from `ModelSummary.jsx` to preserve the module's import boundary; both compute
the same quantities.

---

### Step 3 — Sampling counters (`context.configuration_optimization`)

**Input:** `structure.moves`, parsed from the `.rmc6f` header by `readMovesMetadata()` in
[`browserData.js`](../../web_app/frontend/src/browserData.js) in static mode, and by
`parsers.read_moves_metadata()` in Flask mode (the `/api/structure` response carries `moves` as a
keys-only dict: a counter absent from the header is simply absent). **Both runtimes** therefore produce
this block, and the same counters are also rendered on the Model information card.

**Parse gate.** Both readers scan **only the header**: the text up to the Atoms marker
(case-insensitive `/^[ \t]*atoms\b/im`, the same rule as the atom parser), or the **first 4,000
characters** when that marker is absent or at index 0.
Four regexes are applied (`Number of moves generated/tried/accepted:` and
`Accumulated time (s)…:`), and the function returns `null` unless **at least one** of them matched a
finite number. A run whose header uses different wording therefore yields no block at all rather than a
partial one.

**Outputs.** Three raw counters are emitted verbatim and **unrounded** whenever finite:
`moves_generated`, `moves_tried`, `moves_accepted`. Three derived quantities follow. With $N_{\mathrm{acc}}$
accepted, $N_{\mathrm{tried}}$ tried, $N_{\mathrm{atoms}}$ atoms, $t$ seconds:

$$\text{acceptance\_ratio} = \frac{N_{\mathrm{acc}}}{N_{\mathrm{tried}}}\ \ (\text{2 s.f.}), \qquad
\text{accepted\_moves\_per\_atom} = \frac{N_{\mathrm{acc}}}{N_{\mathrm{atoms}}}\ \ (\text{3 s.f.}), \qquad
\text{accumulated\_time\_h} = \frac{t}{3600}\ \ (\text{3 s.f.})$$

Ratios are emitted only when the denominator is $>0$; the block is dropped entirely if empty. Note the
denominator: the context's `acceptance_ratio` is accepted / **tried**, while the Model information
card's "Accepted / generated" row (`moveStats.js` → `moveRatios()`) divides by moves **generated**.

**Code:** `runContext.js` → `configurationOptimizationContext()`; `browserData.js` → `readMovesMetadata()`.

---

### Step 4 — Run-control settings (`context.run_settings`)

**Input:** the object returned by `parseRunSettings()` (in `browserData.js`) for the `.dat` file whose
stem matches the `.rmc6f` — `KEY :: value` lines with `>`-prefixed sub-lines. **Static mode only**
(Dashboard reads `localRun.settingsFile`; there is no Flask endpoint for it).

**Operations.**

- `title`, `material`, `phase`, `temperature` copied verbatim when non-empty.
- `MINIMUM_DISTANCES` (Å) and `MAXIMUM_MOVES` (Å) are number lists. They are turned into labeled maps
  when — and only when — their length matches the expected label count:
  - minimum distances use **RMCProfile's upper-triangle, row-wise element-pair order** generated by
    `pairLabels()`: for atoms $(e_1,\dots,e_m)$ the sequence is
    $e_1e_1, e_1e_2, \dots, e_1e_m, e_2e_2, \dots, e_me_m$ — i.e. $m(m+1)/2$ labels, the same order as
    the partial-PDF columns. With `ATOMS :: Ga Nb Se` and six values this yields
    `{Ga-Ga, Ga-Nb, Ga-Se, Nb-Nb, Nb-Se, Se-Se}` (pinned in `runContext.test.js`).
  - maximum moves are keyed by element (one per atom type).
  - On a length mismatch the raw array is passed through unlabeled — no guessing.
- `time_limit`, `save_period` copied as strings; `weight_optimization` becomes `true` when the key is
  present; `flags` truncated to the **first 8**.
- `fitted_data`: one entry per `*_DATA` block, `{type, file, fit}` where `type` is the block name with
  the trailing `_DATA` stripped and lower-cased (`XRAY_RECIPROCAL_SPACE_DATA` → `xray_reciprocal_space`)
  and `file` is the `> FILENAME ::` string **verbatim** (see the data-flow note above).
  `parseRunSettings()` itself caps the dataset list at 8 blocks. A block with both
  `> RECIPROCAL_SPACE_FIT` and `> REAL_SPACE_FIT` sub-lines also gets
  `fitted_in: "reciprocal and real space (one measurement: the G(r) is its Fourier transform)"`, and
  `independent_datasets` counts the blocks: an X-ray F(Q) fitted in both spaces is **one** measurement,
  which the model had been reading as two datasets (0.6.0 review).

**Code:** `runContext.js` → `runSettingsContext()`, `pairLabels()`.

---

### Step 5 — Symmetry and per-orbit displacement (`context.symmetry`)

**Input:** the `symmetry` prop = `describeSymmetry(structure, tol)` merged with
`toleranceLadder(structure, 1.0)` (`Dashboard.jsx`). Requires `structure.basis`, so this block exists in
**static mode only**.

**`tol` is a live, user-set control — not a fixed default.** It comes from `SymTolContext`
([`symTolContext.js`](../../web_app/frontend/src/symTolContext.js)), a `useState(0.2)` pair created in
`App.jsx` and shared with the "Detected SG" panel of
[`ModelSummary.jsx`](../../web_app/frontend/src/components/ModelSummary.jsx), where clicking a rung of the
*Space group vs. tolerance* ladder sets `symTol` to that rung's midpoint. **0.2 Å is only the initial
value.** Every field of this block — `space_group`, `point_group`, `n_ops`, `max_residual_A`, and the
entire Wyckoff-orbit partition that decides which sites exist and how `mean_disp_A` is averaged — is
recomputed at whatever tolerance is currently set. Only `tolerance_A` records which one was used. Two
chat sessions on the same run can therefore legitimately receive different space groups and different
site lists. The ladder itself is always computed up to a maximum tolerance of **1.0 Å**.

**Fields.**

- `space_group` — `"F-43m (No. 216)"`. The `(No. N)` suffix is appended **only** when
  `symmetry.spaceGroupNumber` is truthy, so a bare symbol such as `"F-43m"` is possible.
- `point_group`; `n_ops` (order of the space group).
- `tolerance_A` and `max_residual_A` (2 s.f., Å).
- `note` — present only for a structure the finder did **not analyse** (over the 2000-site cap or
  the 384-operation budget): the finder's `reason`, so the model reads why instead of a bare
  "not analysed". Since 0.6.0 `space_group` carries no `(No. N)` when the card shows a crystal class
  (`4/mmm class`), a lower bound (`≥ P4/mmm`), `undetermined` or `not analysed`.
- `point_group_facts` — one line per point group of the group and its ladder
  (`describePointGroup`, `context/pointGroups.js`): crystal system; centrosymmetric or not (the 11
  Laue groups); a Sohncke group (proper rotations only: 1, 2, 222, 4, 422, 3, 32, 6, 622, 23,
  432, the only ones a chiral structure can have) or achiral; polar or not (1, 2, m, mm2, 4, 4mm,
  3, 3m, 6, 6mm). The demo: `-43m: cubic; no inversion centre; achiral (not a Sohncke group);
  non-polar`.
- `ladder` — a list of rungs `{sg, pg, holds_A: [from, to], n_ops}`: the space group (and its point group) the average sites
  satisfy over each Cartesian tolerance interval in Å. **Emitted only when the ladder has more than one
  rung** (`symmetry.ladder.length > 1`); a structure whose ladder collapses to a single rung gets no
  `ladder` key at all, and nothing records that it was suppressed rather than not computed.
- `noise_floor_A` (2 s.f.) and `noise_floor_note` — present when the basis carries counts and
  displacements: the tolerance below which the ladder's verdicts are sampling noise
  (`symmetryNoiseFloor()`). Each average site is the mean of $N_s$ copies with 3-D rms displacement
  $u_s$, so its position carries a 3-D standard error $u_s/\sqrt{N_s}$; two symmetry-related means
  differ by $\sqrt2$ times that from noise alone, and the finder takes the worst case over every site
  and operation, so
  $$\text{noise\_floor\_A} = 3\sqrt2\,\max_s \frac{u_s}{\sqrt{N_s}}.$$
  On the demo (1000 copies, $u \approx 0.23$ Å) it is 0.031 Å — exactly where F-43m first holds, so
  the P1 → Cm → Cmm2 → P-42₁m rungs below it are not evidence of a distortion. The folded average also
  cannot show order that enlarges the cell (the system prompt says so).
- `sites` — one entry per Wyckoff orbit:
  `{element, multiplicity, wyckoff?, site_sym?, frac?, mean_disp_A?, max_disp_A?}`, where
  - `multiplicity` = `orbit.size`, the number of basis sites in the orbit;
  - `wyckoff` = the naming-cell multiplicity and letter
    (`` `${orbit.wyckoffMultiplicity ?? orbit.size}${orbit.wyckoff}` ``, e.g. `"16e"`), present only
    when a letter was identified. It can differ from `multiplicity` when the group is named in a
    cell other than the given one (R3m's `3a` for the Ga of a lacunar spinel kept in its F-cubic
    cell, not `4a`);
  - `site_sym` = the site-symmetry symbol (e.g. `"3m"`), present only when determined;
  - `frac` = the orbit representative's fractional coordinate triple in the conventional cell, rounded
    to **3 s.f.** — this is the atomic-coordinate content of the context (see the data-flow section);
  - `mean_disp_A` / `max_disp_A` as below.

**Per-site displacement folding — the part that is arithmetic here.** For each Wyckoff orbit,
`symmetryContext()` gathers `structure.basis[m].dispA` over the orbit's member indices `orbit.members`
and reports

$$\text{mean\_disp\_A} = \frac{1}{|M|}\sum_{m \in M} u_m, \qquad
\text{max\_disp\_A} = \max_{m \in M} u_m$$

(both 2 s.f., Å), where $M$ is the set of basis sites in that orbit and $u_m$ is that site's rms
displacement. Non-finite entries are filtered out; if none survive, the two keys are omitted.

**Where $u_m$ comes from** (`structureFromRmc6f()` in `browserData.js` — upstream of this module, but
the number is meaningless without it). `structure.basis` holds **one entry per reference-site number**
(the circular-mean position of every atom sharing that reference number), not one entry per atom. For
reference site $s$ with $N_c$ copies (atoms), each copy's within-unit-cell fraction
$w_i = \operatorname{frac}(x_i N_i)$ is offset from the site's circular mean $\bar w_i$, wrapped to the
nearest image, and mapped to Cartesian Å through the unit-cell vectors $\mathbf a_i = \mathbf L_i / N_i$:

$$\Delta\mathbf r = \sum_{i=1}^{3} d_i\,\mathbf a_i,\qquad
u_s = \sqrt{\langle|\Delta\mathbf r|^2\rangle - |\langle\Delta\mathbf r\rangle|^2},$$

the square root of the trace of the site's Cartesian displacement covariance
($\approx\sqrt{3U_\mathrm{iso}}$ of the PCA page: `dispA` divides by $n$, the PCA covariance by $n-1$, a
factor $\sqrt{(n-1)/n}$ — 0.05 % at 1000 copies, ~1.9 % for a 3×3×3 box),
which uses the **full metric** and is the same for any setting of the same cell (run-dashboard.md, Model
summary Part A Step 5). This is a **single-snapshot** spread: it mixes static disorder with thermal
motion.

**Ordering and truncation.** Sites are sorted by `mean_disp_A` **descending** (sites with no
displacement sort last, treated as $-1$), then truncated to `MAX_SITES = 12` with `sites_omitted`
recording the remainder. Ranking first means the budget trimmer in Step 12 keeps the most distorted
sites.

**Code:** `runContext.js` → `symmetryContext()`; upstream
[`symmetryModel.js`](../../web_app/frontend/src/symmetryModel.js) → `describeSymmetry()`,
`toleranceLadder()`, and [`symmetry.js`](../../web_app/frontend/src/symmetry.js) → `siteOrbits()`.

---

### Step 6 — Per-site PCA displacement summary (`context.pca_displacements`)

**Input:** `pcaSites.sites` — the per-reference-site ellipsoid table published by the **PCA Ellipsoid**
page (`onSitesChange` in [`PcaKdePage.jsx`](../../web_app/frontend/src/components/PcaKdePage.jsx)), computed
either by the static-mode worker [`workers/pcaKde.js`](../../web_app/frontend/src/workers/pcaKde.js) or by
`/api/pca/sites` → [`rmc_toolkits/pca_kde.py`](../../rmc_toolkits/pca_kde.py). This block therefore exists
**only if the user has opened the PCA Ellipsoid page at least once** in the session.

**Fields per site** (rounded to 3 s.f. unless noted):

| Key | Meaning | Units | Always present? |
|-----|---------|-------|-----------------|
| `ref` | reference-site number from the `.rmc6f` | — | yes (verbatim `site.referenceNumber`) |
| `element` | element symbol | — | yes (verbatim) |
| `U_iso_A2` | $U_\mathrm{eq}$, mean of the three PCA eigenvalues | Å² | **no** — assigned with no finite check |
| `rms_axes_A` | principal RMS amplitudes $\sqrt{\lambda_1}\ge\sqrt{\lambda_2}\ge\sqrt{\lambda_3}$ | Å | **no** — assigned `undefined` unless `site.rms` is an array |
| `anisotropy` | $\sqrt{\lambda_1/\lambda_3}$ = rms₁/rms₃ | — | **no** — assigned with no finite check; `null` for a zero-spread site |
| `non_gaussianity` | Mardia's multivariate excess kurtosis $(b_2-15)/5$, $b_2=\langle(\mathbf u^\top S^{-1}\mathbf u)^2\rangle$ — rotation-invariant; 0 for a Gaussian, **negative** for a flat-topped or bimodal cloud (a symmetric split site) | — | only when `Number.isFinite(site.nonGaussianity)` |
| `degenerate` | present only when `true` (see below) | — | conditional |
| `zero_spread` | present only when `true`: $\lambda_1 <$ `ZERO_SPREAD_VARIANCE` $=10^{-8}$ Å² (e.g. an `*AVERAGE.rmc6f`); `anisotropy` is then `null` and `non_gaussianity` absent | — | conditional |
| `mixed`, `element_counts` | present only for a mixed-occupancy site: `element` is then the majority species only, and `element_counts` the full tally | — | conditional |

**Silent field dropping.** `U_iso_A2`, `rms_axes_A` and `anisotropy` are assigned unconditionally, and
`roundSig` passes a non-finite value through unchanged: an `undefined` source deletes the key in
`JSON.stringify`, a `NaN` is serialised as `null`, and a `null` (a zero-spread site's `anisotropy`, as
both engines send it) stays `null`. There is no `*_omitted` marker, so an incomplete PCA row can arrive
at the model as `{ref, element}` alone with nothing explaining why; only a zero-spread site says so,
through `zero_spread`. `non_gaussianity` is the only field with an explicit finite guard.

**`degenerate` means rank-deficient, not isotropic.** It is set when the *smallest* eigenvalue is
negligible against the largest:

$$\texttt{degenerate} \iff \frac{\lambda_3}{\lambda_1} < \texttt{DEGENERATE\_RATIO} = 10^{-6}$$

i.e. a flat (planar) or linear displacement cloud, for which the ellipsoid and the analytic
eigen-formulas are unreliable. It does **not** mean "the eigenvalues are nearly equal". The JS worker
(`workers/pcaKde.js`, `DEGENERATE_RATIO` and `degenerate: eigenvalues[2] / largest < DEGENERATE_RATIO`)
and the Python engine (`rmc_toolkits/pca_kde.py`, `bool(ratio[index] < DEGENERATE_RATIO)`) use the same
constant and the same test — they agree exactly.

**The site partition can be synthetic.** When the loaded configuration carries no reference-site or cell
columns, `PcaKdePage.jsx` rebuilds sites by folding every atom into one unit cell and clustering atoms of
the same element within a **user-set distance** (`clusterThreshold`, a 0.4–2.5 Å slider in 0.1 Å steps,
default **1.5 Å**), and the resulting table is flagged `sites.reconstructed`. Moving that slider
re-clusters and **re-numbers** the sites. `pcaContext()` copies `site.referenceNumber` verbatim and does
**not** copy the `reconstructed` flag or the cluster distance, so the model cannot distinguish a real
reference-site partition from a reconstructed one, and the same run can produce different `ref` numbering
in two different chat sessions.

**Ordering:** descending $|$`non_gaussianity`$|$ (a missing value counts as 0), ties broken by
descending `U_iso_A2` (missing values sort as $-\infty$), then truncated to `MAX_SITES = 12` with
`sites_omitted`. The ranking is by magnitude because a symmetric split site is *negative*: before 0.6.0
the signed ranking listed split sites last, where the cap trims first. A fixed `note` string travels
with the block defining every symbol and stating the sign convention (>0 peaked / heavy-tailed or a
minority off-centre component, <0 flat-topped or bimodal), because small models otherwise misread it;
the system prompt repeats the definition and asks the model to call out large $|$`non_gaussianity`$|$
of either sign.

**Code:** `runContext.js` → `pcaContext()`. The ellipsoid math itself is documented in the PCA
Ellipsoid section; nothing is recomputed here — this step is a projection, a sort, and a truncation.
Note that the page's `probability` (ellipsoid percentile) setting does not affect any field copied here.

---

### Step 7 — Average-structure nearest-neighbour distances

**Input:** `structure.basis` (one representative site per reference number, `{el, frac}` with `frac`
fractional in the *conventional* cell), `latticeVectors`, `supercell`. Returns `null` without a basis.

**Math.** Conventional cell matrix (rows, Å): $\mathsf A_{ij} = L_{ij}/\max(N_i,1)$. For basis sites
$p,q$ with fractional coordinates $\mathbf x_p,\mathbf x_q$ and an image offset
$\mathbf n \in \{-1,0,1\}^3$,

$$d_{pq}(\mathbf n) = \left\lVert (\mathbf x_q + \mathbf n - \mathbf x_p)\,\mathsf A \right\rVert
\quad\text{(row-vector convention: } c_k = \textstyle\sum_m \Delta f_m\, \mathsf A_{mk}).$$

For each unordered element pair $(A,B)$ the reported distance is

$$d^{\,\mathrm{avg}}_{AB} = \min_{\substack{p\le q,\ \{e_p,e_q\}=\{A,B\} \\ \mathbf n \in \{-1,0,1\}^3,\ d>10^{-3}\,\text{Å}}} d_{pq}(\mathbf n),$$

with the self-term ($p=q$, $\mathbf n=\mathbf 0$) skipped and any distance $\le 10^{-3}$ Å discarded so
duplicated basis entries cannot produce a spurious zero. Keys are `"El1-El2"` with the two symbols sorted
by `Array.prototype.sort()` — UTF-16 code-unit order, which coincides with alphabetical order for
chemical symbols. Test fixture: an 8 Å box with a $2\times2\times2$ supercell, Na at $(0,0,0)$ and Cl at
$(\tfrac12,0,0)$, gives Na–Cl = 2 Å and Na–Na = Cl–Cl = 4 Å.

**Self-pairs can degenerate into a lattice repeat.** Only the $\mathbf n=\mathbf 0$ image is skipped for
$p=q$. If an element occupies **exactly one basis site**, the only surviving candidates are the 26
non-zero $\pm1$ lattice translations, so `avg_structure_d_A` for that A–A pair is simply the **shortest
conventional-cell lattice translation** — the cell edge in the fixture above — and not a coordination
distance at all. For a low-symmetry or supercell-derived conventional cell this number is meaningless as
a bond length, yet it is handed to the model beside a real $g(r)$ peak and the system prompt invites a
comparison. **Nothing in the context marks such entries.**

**Approximation:** only $\pm1$ images are searched. Since fractional coordinates lie in $[0,1)$, offsets
of $\pm1$ cover every minimum-image candidate for a reasonably compact cell, but a strongly skewed cell
could in principle place the true nearest image outside that shell.

**Code:** [`llm/context/pairCorrelations.js`](../../web_app/frontend/src/llm/context/pairCorrelations.js) →
`nearestNeighborDistances()`, `conventionalCell()`, `cartDistance()`, `pairKey()`.

---

### Step 8 — $g(r)$ peak extraction

This is the app's own peak finder — it is not RMCProfile output, and it is intentionally crude.

**Input:** one series `(x, y)` from the partials file: $x$ = $r$ in Å (assumed ascending), $y$ = the
partial pair-correlation value as parsed from the CSV (no normalization, no baseline model). Options:
`rMax = 6` Å, `maxPeaks = 2` — both hard-coded defaults, never overridden by any caller.

**Procedure (`extractPeaks()`).**

1. **Window.** $\mathrm{end} = \min\{i : x_i > r_\mathrm{max}\}$, or the full length if no such $i$ exists.
   If $\mathrm{end} < 5$ the function returns `[]`. Only $r \le 6$ Å is ever examined.
2. **Smoothing.** A 5-point moving average with truncated (not padded) edges — the window is clipped at
   **both** ends of the analysis window (`j >= 0 && j < end`):

   $$s_i = \frac{1}{|W_i|}\sum_{j \in W_i} y_j, \qquad W_i = \{\,j : |j-i| \le 2,\ 0 \le j < \mathrm{end}\,\}$$

   so the divisor runs 3, 4, 5, …, 5, 4, 3 across the window: 3 at $i=0$ and $i=\mathrm{end}-1$, 4 at
   $i=1$ and $i=\mathrm{end}-2$, 5 in the interior. This is a *uniform-index* filter: its support in Å is
   5 samples spanning $4\,\Delta r$ end to end, whatever the file's $r$ step is (0.08 Å on the test's
   0.02 Å grid). The right-edge samples are not inert — they are the
   $s_{i+1}, s_{i+2}$ comparands for the last candidate peak ($i = \mathrm{end}-3$), and the half-height
   walk runs out to $\mathrm{end}-1$. **A peak sitting just under the 6 Å cut is therefore accepted,
   rejected, and widened against edge-biased averages, and is the least reliable peak reported.**
3. **Floor.** $f = 0.15\,\max_i s_i$ over the window. If the window maximum is not $>0$ the function
   returns `[]` (so an all-zero or all-negative curve yields no peaks).
4. **Local maxima.** For $2 \le i \le \mathrm{end}-3$, index $i$ is a peak when $s_i \ge f$ and

   $$s_i \ge s_{i-1},\quad s_i \ge s_{i-2},\quad s_i > s_{i+1},\quad s_i > s_{i+2}$$

   The asymmetry — $\ge$ on the left, **strict** $>$ on the right — means a flat top is reported once at
   its **right-most** index: on a plateau $s_{i+1}$ equals $s_i$, so every plateau index except the last
   fails the strict test. The reported $r$ of a plateau is therefore biased to its **trailing edge**.
5. **Position.** The reported `r` is `x[i]`, the grid abscissa of the maximum of the **smoothed** array —
   not of the raw data, and with **no sub-grid interpolation**. Position resolution therefore equals the
   file's $r$ step (0.02 Å in the test fixtures), and for an asymmetric peak the 5-point average shifts
   the reported maximum toward the heavier flank.
6. **Width.** From the peak index, walk outward while the smoothed curve stays above **half the peak
   value**: `left` decreases while $s_\mathrm{left} > s_i/2$ (stopping at 0), `right` increases while
   $s_\mathrm{right} > s_i/2$ (stopping at $\mathrm{end}-1$); then

   $$\mathrm{FWHM} \approx x_{\mathrm{right}} - x_{\mathrm{left}}.$$

   There is **no baseline subtraction and no interpolation** — the half-height is measured from zero, and
   the crossing is snapped to the grid. Consequences, stated honestly: for overlapping peaks or a curve
   with a raised baseline the walk runs past the true half-height crossing and the width is an
   **overestimate**; in the limit where the curve never falls below half the peak, the reported FWHM is
   the whole window.
7. **Selection.** Peaks are already in ascending $r$; the function returns `peaks.slice(0, 2)` — the
   **first two by distance**, not the two tallest.

**Validation — what the tests actually assert.**
`__tests__/pairCorrelations.test.js` feeds three Gaussians at 2, 3.5, and 7 Å on a 0.02 Å grid and
asserts that exactly two peaks are found below 6 Å, at $r \approx 2$ and $r \approx 3.5$ Å (to one
decimal place). The width assertion is only a loose bracket, $0.15 < \mathrm{fwhm} < 0.4$ Å, around the
analytic $2.355\sigma = 0.2355$ Å — a $-36\,\%/+70\,\%$ window, so the test does **not** establish
accuracy better than that. A second test adds a $0.3\sin(40r)$ ripple to a peak of height 10 and asserts
only `peaks.length >= 1` and `peaks[0].r ≈ 3` — i.e. that the first reported peak is still the real one.
It does **not** assert that the 15 % floor eliminated every ripple-induced extra peak.

---

### Step 9 — Assembling `context.pair_correlations`

**Input:** the **first** entry in `plotFiles` with `plotKind === 'pdf_partials'` and a non-empty
`plotData.series`. That kind is assigned by file name: `*PDFpartials*.csv`
(`detectPlotKind()` in `browserData.js`; `detect_plot_kind()` in
[`rmc_toolkits/plots.py`](../../rmc_toolkits/plots.py) implements the identical rule for this branch — see
the last bullet of [Caveats / what this is not](#caveats--what-this-is-not) for the one branch, `stog`,
where the two do diverge). If a run has several partials files, only the first in dashboard
sort order (Step 11) is used.

**Per series:**

1. The series label must match `/^([A-Z][a-z]?)\s*-\s*([A-Z][a-z]?)$/` — two 1- or 2-character element
   symbols separated by a hyphen. This rejects totals columns, unparsed headers, **3-letter symbols**,
   and any labelled or indexed species (`"Fe1-Fe2"`, `"Nb-Nb total"`). Rejected series are skipped
   silently, with no counter.
2. `avg_structure_d_A` = the Step-7 distance for that pair, looked up under the **sorted** key
   `pairKey(el1, el2)`, 3 s.f., when a basis exists.
3. `gr_peaks_A` = the Step-8 peaks as `[{r (3 s.f.), fwhm (2 s.f.)}]`, both in Å.
4. The `pair` field is the **verbatim CSV column label**, in whatever order the file writes it. Because
   the distance is keyed on the sorted pair, the model can legitimately receive
   `{pair: "Se-Ga", avg_structure_d_A: …}` where the distance was stored under `"Ga-Se"`. The two orders
   describe the same unordered pair; only the string differs.
5. `first_shell` = `{mean_d_A, to_A, cn_<B>_per_<A>?}` from `firstShell()` (below), when a shell
   boundary is found: the shell's mean distance from $r^2 g(r)$, its outer edge, and — with the cell —
   the coordination number of B around A.
6. The entry is kept only if at least one of `avg_structure_d_A` / `gr_peaks_A` is present.

**The first coordination shell (`firstShell()`).** On the 5-point moving average $\tilde g$ of the
partial (normalised to 1 at large $r$, as RMCProfile writes it), with the first peak of Step 8 at
height $h$:

- the **onset** is the last point before the peak where $\tilde g$ rises above $0.01\,h$;
- the **boundary** is the first point after the peak that is either a local minimum below
  `SHELL_TROUGH_FRACTION` $= 0.2\,h$ — a shallower dip lies *inside* a split shell, as the demo's
  two Ta–Se distances (2.52 and 2.74 Å) of one TaSe₆ octahedron — or where $\tilde g < 0.001\,h$
  (a tail decaying into an empty gap: 3.7σ out on a Gaussian shell);
- the shell's distance statistics come from the radial distribution $r^2 g(r)$ (raw $g$, negatives
  clamped) between onset and boundary, by the trapezoid rule:
  $$\bar d = \frac{\int r^3 g\,dr}{\int r^2 g\,dr},\qquad
  \sigma_d^2 = \frac{\int (r-\bar d)^2 r^2 g\,dr}{\int r^2 g\,dr};$$
  the $g(r)$ peak itself sits *below* $\bar d$ for a broad or skewed shell, because $g$ divides by
  $r^2$ (demo Ta–Ta: peak 2.93 Å, mean 3.00 Å, average-structure distance 3.02 Å — the model had read
  the peak as a 0.09 Å contraction);
- with the partner's number density $\rho_B = N_B/V$ (`elementDensity()`: the element count over
  the supercell volume), the coordination number
  $\mathrm{CN}_{B|A} = 4\pi\rho_B\int_0^{r_b} r^2 g_{AB}(r)\,dr$, and the same integral to
  $r_b \mp 0.2$ Å (`cnRange`): how much the count depends on the cut. The reverse count is
  $\mathrm{CN}_{A|B} = \mathrm{CN}_{B|A}\,\rho_A/\rho_B$. Demo checks: Ga 3.99 Se, Ta 6.00 Se,
  Se 12.0 Se, and Ta 2.76 Ta (the intracluster 3 with the Ta–Ta gap at 3.3–4.0 Å not empty).

**Cap:** at most `maxPairs = 10` accepted entries (the loop breaks once 10 have been pushed); the function
returns `null` if none. Nothing records how many pairs were dropped at this stage.

**Why both numbers travel together.** $d^{\,\mathrm{avg}}_{AB}$ is what the *average* structure says the
first-neighbour distance is; the $g(r)$ peak is what the *instantaneous* configuration says. A peak
displaced from — or split around — $d^{\,\mathrm{avg}}$ is exactly the short-range correlation that an
average structure cannot express, and pairing them is the whole point of this block. (Except for the
single-basis-site self-pairs of Step 7, where $d^{\,\mathrm{avg}}$ is a lattice repeat and the comparison
is meaningless.)

**Code:** `pairCorrelations.js` → `pairCorrelationsContext()`.

---

### Step 10 — Per-dataset residuals (`context.datasets`)

**Input:** all `plotFiles` with a truthy `plotKind` other than `r_value` (R-value logs are convergence
data, not fitted datasets — asserted in the tests). Two upstream filters decide what can ever reach here:

- **STOG files never appear.** `Dashboard.jsx` builds its file list with
  `isDashboardPlotFile = (file) => file.plotKind && file.plotKind !== 'stog'`, so `.gr`/`.sq`/`.fq`
  datasets are dropped *before* `allPlotFiles` exists. `datasetContext()` would happily accept them —
  and the unit-test fixture contains a `stog` entry — but the running app can never produce one.
- **Hidden plots still count** (Step 1): `datasets` is built from the unfiltered list, so a plot you
  hid on the Dashboard is still described to the model.

**Output per file:** `{kind, title?, rwp?}`, where `title` is included only when it differs from `kind`,
and `rwp` is `plotData.metrics.rwp` at 3 s.f.

**Which kinds actually carry an `rwp`.** `plotDataFromText()` in `browserData.js` computes
`metrics.rwp` **only** for `kind ∈ {xpdf, npdf, xray_sq, neutron_sq, bragg}` **and** only when the parsed
CSV has at least 3 columns. RMCProfile writes those CSVs as $(x, y^{\mathrm{calc}}, y^{\mathrm{expt}})$ —
`Q, F(Q)_RMC, F(Q)_Expt`, `r(A), X_ray-calc, X_ray_exp_renorm` — so by default column 2 is the RMC
calculation and column 3 the experiment; a header that names both roles (`/exp|obs/i` vs
`/calc|rmc|fit/i`) overrides that order (`rwpColumns()`, mirrored by `rwp_columns()` in Python). `pdf_partials`,
`exafs_q`, `exafs_r` (and `r_value`, which is excluded anyway) never have one, so in a real run most
`datasets` entries are `{kind}` or `{kind, title}` with no residual — their absence is not an error.

The metric is computed by `fitRwp()` → `rwp(x, observed = y^expt, fitted = y^calc)` in `browserData.js`;
the Python equivalent is `fit_rwp()` → `rwp()` in [`rmc_toolkits/parsers.py`](../../rmc_toolkits/parsers.py).
**Both compute the same unweighted quantity over the rows finite in both columns, normalized by the
experiment:**

$$R = \sqrt{\frac{\sum_i (y^{\mathrm{calc}}_i - y^{\mathrm{expt}}_i)^2}{\sum_i (y^{\mathrm{expt}}_i)^2}}$$

with no weights $w_i$ — despite the name `rwp` and despite the system prompt describing it as "a weighted
profile residual". Treat it as a relative goodness indicator between datasets of the same run, not as a
Rietveld $R_{wp}$. **Undefined cases:** when no row is finite in both columns, or
$\sum_i (y^{\mathrm{expt}}_i)^2 = 0$, both implementations return the unavailable sentinel
`None`/`null` (never `0`, which would read as a perfect fit). `datasetContext()` keeps `rwp` only when it
is finite, so such a dataset reaches the model as `{kind, title?}` with **no `rwp` key** — an `rwp: 0`
entry cannot be produced by an undefined metric.

A file named `<stem>_FT_XFQ<n>.csv` (`xpdf`) is RMCProfile's Fourier transform of the X-ray F(Q)
it fits, and carries `same_measurement_as: "XFQ<n> (its Fourier transform)"`.

**Partial weights (`context.partial_weights`).** For each radiation the run fits — X-ray when a
`.dat` block name contains `XRAY` or an `xpdf` / `xray_sq` file exists, neutron for `NEUTRON` or
`npdf` / `neutron_sq` / `bragg` — the Faber–Ziman weight of every partial in the total,

$$w_{ab} = (2-\delta_{ab})\,\frac{c_a c_b f_a f_b}{\left(\sum_i c_i f_i\right)^2},\qquad
\sum_{a\le b} w_{ab} = 1,$$

with $c$ the atom fractions from the element counts, $f = Z$ for X-rays (the $Q\to0$ limit) and
$f = b_\mathrm{coh}$ for neutrons (`COHERENT_B_FM`, the same NIST table as
`rmc_toolkits/scattering.py`; real part for B, Cd, Dy, Eu, Gd, In, Sm), sorted by $|w|$, 2 s.f.
(`scatteringWeights.js` → `partialWeights()`). GaTa₄Se₈ under X-rays: Ta–Se 0.45, Ta–Ta 0.24,
Se–Se 0.21, Ga–Ta 0.051, Ga–Se 0.048, Ga–Ga 0.0027 — anything the run says about Ga rests on ~5 % of
the signal. An element without a tabulated amplitude drops that radiation's block.

**Code:** `runContext.js` → `datasetContext()`, `partialWeightsContext()`; `browserData.js` →
`plotDataFromText()`, `rwpColumns()`, `fitRwp()`, `rwp()`.

---

### Step 11 — Convergence series (`context.convergence`)

**Input:** `rValueFile.plotData.series[0].y`.

**What the values are.** RMCProfile's `<stem>-NN.log` files are read by `readChi()` in
`browserData.js` (mirroring `read_chi_log()` in Python): skip the first 2 lines, keep only data rows with
exactly as many tokens as line 1 names, drop an unterminated final line (a row RMCProfile is still
writing), and take the **last column** of each row — a row with a non-finite χ² (`NaN`, `Inf`, Fortran
`****`) is **kept as `NaN`**, never dropped (run-dashboard.md, Parsing Step 4d). The dashboard then stores
$y_k = \ln\!\left(\max(\chi_k, 10^{-12})\right)$ — the natural log, with a hard floor of $10^{-12}$ on the
raw $\chi$ before the log, so a $\chi$ of exactly 0 maps to $\ln(10^{-12}) \approx -27.63$ rather than
$-\infty$ — and $y_k = $ `NaN` for a non-finite $\chi_k$.

**How the logs are combined.** `combineRValueFiles()` in
[`browserData.js`](../../web_app/frontend/src/browserData.js) (called by `Dashboard.jsx`) **concatenates**
the log-transformed series of the visible R-value logs **of one run** — the logs sharing the folder and
stem of the structure file (else the first stem), as Flask's `related_r_value_logs()` does — and
re-indexes $x$ as $0,1,2,\dots$; `final_chi_r` is taken from the *last* parsed file of that run. Logs of
other runs in the folder are never spliced in (they used to be, in static mode, which fed the model
another run's `final_chi_squared` and tail slope). The order is `comparePlotFiles`: first by `plotOrder` index of `plotKind`, then — for
names matching `/^(.+)-(\d{2,})\.log$/` — by lower-cased stem and then by the **numeric log sequence
number**, and otherwise by a numeric-aware `localeCompare` of the file name.

Three early-return branches make the concatenation conditional. `combineRValueFiles()` returns
`rValueFiles[0]` **unchanged** when

1. there is only one R-value file, **or**
2. no file has any of `sourceFile` / `plotData` / `parseError`, **or**
3. **any** file has a `sourceFile` but neither `plotData` nor `parseError` — i.e. it is still being
   parsed, which is routine during a Live Data refresh.

So while a sibling log is mid-parse, the assistant's `convergence` block reflects **only the first log**,
not the concatenation, and it silently switches back to the full concatenation on the next poll. If no
file parses at all, the result carries a `parseError`, `seriesStats()` returns `null`, and the
`convergence` block is absent entirely.

So "step" = one log line (one save/print period), across restarts — it is **not** an RMCProfile move
count, and the concatenation assumes the logs are in chronological order.

**Statistics on the full series** (before any downsampling — this is deliberate):

$$\mathrm{nSteps} = N,\quad \mathrm{first} = y_0,\quad \mathrm{last} = y_{N-1},\quad
\min_k y_k,\quad \max_k y_k,\quad \mathrm{recentSlopePerStep} = m$$

with $m$ the tail least-squares slope defined in Step 14. `first`/`last` are the raw end values (a
non-finite latest value serializes as `null` — visible, not hidden); `min`, `max` and $m$ use only the
finite values, in one pass (no argument spread). When any value is non-finite the block also carries
`non_finite_steps` (the count) and a `non_finite_note` telling the model those rows are a blown-up run,
not missing data. On output, `n_steps` is the **exact integer point count and is never rounded**;
`first`, `last`, `min`, `max` are 3 s.f.; `recent_slope_per_step` is 2 s.f.

The block also carries `watchdog_status` — `classifyConvergence()` of the same series (Step 14),
the dashboard's verdict — and, with ≥ 10 finite rows, `change_over_last_10pct` = the last finite
value minus the one a tenth of the rows earlier (2 s.f.; −0.1 ≈ a 10 % χ² drop). Without them the
model read the demo's still-falling 48-point trace (−0.11 over the last tenth) as a plateau and called
the run converged.

**Downsampling (`downsampleSeries`).** Uniform-stride resampling to at most `HISTORY_POINTS = 48` points
that always retains the endpoints:

$$\mathrm{step} = \frac{N-1}{P-1}, \qquad \tilde y_i = y_{\,\mathrm{round}(i\cdot\mathrm{step})},\quad i = 0,\dots,P-1$$

for $P = 48$ (or fewer, see Step 12). Series with $N \le P$ pass through at full length. Every retained
value is rounded to 3 s.f. This is *decimation, not averaging* — no local minimum between kept indices
survives, which is why the min/max/slope statistics are computed on the full series first.
(`HISTORY_POINTS` is also the default argument of the exported `downsampleSeries`.)

**Also emitted:** `quantity` — `"ln of the chi^2 in the last .log column 'X_ray_(R)1' (natural log;
lower is better) — one fit term of the run, not a total"` (the column name comes from the log header,
`plotData.chiColumn`, and is omitted when the log has none) — `column` (that header name), and
`final_chi_squared` = `plotData.metrics.final_chi_r`, the **un-logged** last value of that column, so
$\mathrm{last} = \ln(\text{final\_chi\_squared})$ up to rounding. The series is **one term** of the
fit (in the demo run the X-ray real-space χ²; the reciprocal-space `F(Q)_1` term is not included), which
is why the context names it rather than calling it the run's goodness of fit.

> **Naming (trust the code).** Two labels attached to this data still need care:
>
> 1. **One term, named.** Every producer now labels the series $\ln(\chi^2)$ of the named log column
>    (`CHI_HISTORY_Y_LABEL = 'ln(χ²)'`, title `χ² history: X_ray_(R)1`), and the context says the same;
>    only the historical metric key `final_chi_r` keeps the old name. The column is one fit term (in the
>    demo run the X-ray real-space χ²), not the run's total χ².
> 2. **`rwp` is not weighted** (Step 10), despite its name and the system prompt's wording.
> 3. **The x-axis is not a move count.** `SYSTEM_PROMPT` in
>    [`llm/prompts/system.js`](../../web_app/frontend/src/llm/prompts/system.js) tells the model, on
>    *every single request*, that "the convergence history is the natural log of the chi-squared goodness
>    metric **versus accepted-move steps**". It is not: the index is one point per log line, concatenated
>    across restart logs (`combineRValueFiles()` sets `x = 0,1,2,…` and the plot label is
>    `'Time steps'`). The model is therefore primed to misread the axis, and **any statement it makes
>    about "moves" derived from the convergence block should be discounted.** The move counters in
>    `configuration_optimization` (Step 3) are the only real move numbers in the context.

**Code:** `runContext.js` → `convergenceContext()`, `downsampleSeries()`, `seriesStats()`;
`Dashboard.jsx` → `combineRValueFiles()`, `comparePlotFiles()`, `rValueLogParts()`.

---

### Step 12 — The character budget (`contextToJson`)

**Input:** the assembled context object and a budget, default `CONTEXT_CHAR_BUDGET = 6000` **characters**
(not tokens; ~1.8k tokens). Serialization is `compactJson()`: `JSON.stringify(context, null, 1)` with
every array of plain values, and then every object of plain values (and such arrays), collapsed onto
one line — about half the length of one value per line, and still valid JSON (a bracket or brace
followed by a line break never sits inside a JSON string). Both demo-like runs (GaTa₄Se₈ 250 K,
GaNb₄Se₈ 5 K) fit with no trim at ~5,000 characters.

**Trim ladder.** Each step is applied only if the current serialization still exceeds the budget:

| Order | Condition | Action | Bookkeeping |
|-------|-----------|--------|-------------|
| 1 | any pair has >1 `gr_peaks_A` | keep only the first peak per pair | **none** |
| 2 | `symmetry.ladder.length > 2` | keep first + last rung | `ladder_rungs_omitted` |
| 3 | `convergence.history.length > 24` | re-downsample to 24 points | **none** |
| 4 | `convergence.history.length > 12` | re-downsample to 12 points | **none** |
| 5 | `pca_displacements.sites.length > 6` | keep top 6 (most non-Gaussian) | `sites_omitted` (accumulated) |
| 6 | `pca_displacements.note` present | drop the explanatory note | **none** |
| 7 | `symmetry.sites.length > 8` | keep top 8 (most displaced) | `sites_omitted` (accumulated) |
| 8 | `datasets.length > 12` | keep first 12 | `datasets_omitted` |

**Only half the ladder is recorded.** The trims that discard whole *ranked entries* — ladder rungs, PCA
sites, symmetry sites, datasets — write an `*_omitted` count. Dropping a pair's second $g(r)$ peak,
re-downsampling the history (twice), and dropping the PCA `note` record **nothing**: the context simply
arrives with less evidence and no marker. (The function's own header comment claims every truncation is
recorded; the code disagrees, and the code is what runs.)

**The history re-downsample is nested.** Rows 3 and 4 call
`downsampleSeries(c.convergence.history, …)` on the array `convergenceContext()` already produced — the
48-point, already-3-s.f.-rounded trace, not the source series. The resulting 24- and 12-point traces are
therefore a decimation of a decimation: a strict subset of the 48 kept indices, generally *not* the same
indices a single uniform resample of the full $N$-point series to 24 or 12 would pick, and the values
carry the earlier rounding.

The rationale is visible in the order: the *least* essential redundancy goes first (a second $g(r)$ peak,
intermediate ladder rungs), the ranked evidence blocks are trimmed from the bottom of their own ranking
(rows 5 and 7 exploit the sorts of Steps 5–6), and the flat dataset list is the last resort. Rows 3–4
contain the function's only early return: if the budget is met before either history trim, the current
JSON is returned immediately.

**This is best-effort, not a guarantee.** Several sections are *never* trimmed — `structure`,
`configuration_optimization`, `run_settings` (including a long `min_distances_A` map for a many-element
system), the `convergence` summary statistics, and the pair list itself (only extra peaks are dropped,
never whole pairs). A pathological run can therefore exceed 6,000 characters after all eight steps. The
tests assert the budget holds for a 30-dataset, 50,000-step fixture and that a forced budget of 10
characters triggers every trim step in order, but no test asserts the budget is met unconditionally.

**Code:** `runContext.js` → `contextToJson()`.

---

### Step 13 — Message assembly (what the request body actually contains)

**Copilot** (`buildCopilotMessages()` in
[`llm/prompts/copilot.js`](../../web_app/frontend/src/llm/prompts/copilot.js)), in order:

1. `system`: `COPILOT_PROMPT`, the fixed `SYSTEM_PROMPT` from
   [`llm/prompts/system.js`](../../web_app/frontend/src/llm/prompts/system.js) plus the tool rules
   and the Outcome rule (Step 17). Without tools it is `COPILOT_PROMPT_NO_TOOLS` (the Outcome rule
   only). There is no run data in it, but it carries hard numeric interpretation rules (below).
2. `user` — `"Context for this RMC modeling run:"` followed by ```` ```json … ``` ```` containing
   `contextToJson(context)`. The fence exists so the model can distinguish data from instructions.
3. `assistant` — a fixed one-line acknowledgement.
4. the last `CHAT_HISTORY_TURNS = 8` memory messages (`slice(-8)`: 8 *messages*, i.e. 4 exchanges).
5. `user`: the question you typed.

Then, within the question, the tool exchange (Step 17): one `assistant` message with `tool_calls` per
round, and one `tool` message per call.

**What "memory" contains.** The turns live in `useCopilot`'s state: never persisted, and cleared on
reload, on *New conversation*, and when another run opens (`runKey`).
`conversationMemory()` ([`llm/useCopilot.js`](../../web_app/frontend/src/llm/useCopilot.js)) sends
each earlier question with its **final answer text** only.

- Tool calls, tool results and chain-of-thought (`reasoning` / `reasoningMs`) are **never re-sent**.
  The model re-runs an analysis it needs again.
- A question whose turn produced no answer (an empty or failed reply) is dropped together with it, so
  memory never holds two user messages in a row.
- A stopped turn keeps its partial text, marked "(Interrupted before the answer was finished.)".

The full context is **re-sent with every request**, which is why old turns can be dropped cheaply.
The request body is `{model, messages, temperature, stream: true}`, plus `tools` and `tool_choice`
when analyses are on. `temperature` comes from settings (default **0.2**).

#### Numeric rules the system prompt asserts (no code backs them)

These shape every answer but appear nowhere in the pipeline above, and nothing in the codebase enforces
or cites them. Quoted from `SYSTEM_PROMPT`:

| Assertion | Prompt text (abridged) | Code backing |
|---|---|---|
| Rwp acceptability band | "typical acceptable RMC fits land around 0.01-0.10, but the threshold is dataset-dependent" | **none** — no threshold exists in code; and the quantity is unweighted (Step 10) |
| `accepted_moves_per_atom` rubric | "single digits mean the configuration is lightly sampled … tens or more indicate a well-sampled run" | **none** — the value is only computed, never compared |
| Ladder → displacement heuristic | "the tolerance where a higher-symmetry group first holds approximates the mean static displacement" | **none** — an interpretive rule of thumb, not a derived relation |

Treat all three as prompt-level assertions from the app author, not as validated results. They are listed
in the parameters table below for completeness.

**Watchdog** (`buildWatchdogMessages()`): the same system prompt plus one user message containing
`JSON.stringify(stats)` (Step 14) — introduced as the ln χ² of the named fit term (`chiColumn`, e.g.
`X_ray_(R)1`, "one .log column, not the total") when the column is known — the heuristic label, the
previous label if any, and a demand for
exactly one line of the form `` `STATUS: improving|converged|stalled|diverging — <one short sentence
citing a number>` ``. **This payload is not fenced** — unlike the chat context, the stats JSON is
interpolated into the middle of an English sentence, so the data/instruction separation the chat path
argues for is not applied on the watchdog request path. (The payload is five numbers, so the risk is
small, but the asymmetry is real.) The message also discloses both the current heuristic label and the
previous one to the model. Sent non-streaming at `temperature: 0`. The reply is parsed by

```js
/STATUS:\s*(improving|converged|stalled|diverging)\s*[—–-]?\s*(.*)/i
```

returning `{status, note}` or `null`. **The parsed `status` is not used to set the badge** — only the
`note` is (see Step 15).

---

### Step 14 — Convergence heuristics (the statistics)

Source of truth: [`llm/watchdog/heuristics.js`](../../web_app/frontend/src/llm/watchdog/heuristics.js).
These are pure functions over the $\ln\chi$ array; they work with **no model connected at all**.

**Recent window.** For a series of $N$ points,

$$n = \max\!\big(\min(N,5),\ \lceil 0.2N \rceil\big)$$

i.e. the last 20 % of the series, but never more than $N$ and never fewer than 5 **once $N \ge 5$**. For
$5 \le N \le 25$ this is exactly 5 points; for $N < 5$ it is the whole series ($n = N$), which is
reachable because the watchdog runs from $N \ge 2$; for $N = 500$ it is 100 points.

**Tail slope (`recentSlope`).** Ordinary least squares of the last $n$ values against their index
$i = 0,\dots,n-1$:

$$m = \frac{\sum_i (i - \bar\imath)(y_i - \bar y)}{\sum_i (i - \bar\imath)^2},
\qquad \bar\imath = \frac{n-1}{2},\qquad \bar y = \frac1n\sum_i y_i$$

with $m = 0$ when the denominator vanishes ($n<2$) or the input has fewer than 2 points. Units:
$\ln\chi$ per log step.

**Window delta — the actual decision variable.**

$$\Delta_{\mathrm{win}} = m \,(n-1)$$

the model-predicted change in $\ln\chi$ across the whole recent window. Because the series is a
logarithm, an additive shift is a *relative* change in $\chi$:
$\Delta_{\mathrm{win}} = \ln(\chi_{\mathrm{end}}/\chi_{\mathrm{start}})$, so the threshold is
magnitude-independent — the same rule works for $\chi \sim 10^{-2}$ and $\chi \sim 10^{2}$.

> **$m$ and $\Delta_{\mathrm{win}}$ are different numbers, and both are shipped.** The chat context carries
> the *per-step* slope, `recent_slope_per_step` $= m$ (Step 11). The watchdog payload carries the
> *window* delta, `recent_window_delta` $= m\,(n-1)$. They differ by the factor $n-1$, which for a long
> run is $0.2N - 1$. Do not compare the two fields as if they were the same quantity.

**Thresholds (module constants).**

- `WINDOW_DELTA_EPSILON = 0.02` — $|\Delta_{\mathrm{win}}| \le 0.02$ counts as flat.
  $e^{0.02}-1 = 2.02\,\%$, so this is "the fit moved less than ~2 % across the recent window".
- `MIN_TOTAL_DROP = 0.1` — a run must have improved by at least 0.1 in $\ln\chi$ *overall*
  ($e^{0.1}-1 = 10.5\,\%$) to earn "converged" rather than "stalled".

**Classification (`classifyConvergence`).**

Non-finite values (the NaN rows above, or `null`) are handled first: **a non-finite latest value is
`diverging`** — the run produced NaN/Inf χ² — and otherwise the rule below runs on the **finite values
only** (an isolated earlier NaN does not break the trend; `unknown` if fewer than 2 remain). The
browser used to drop such rows, so a log whose last 40 of 100 rows were NaN was classified `improving`
from its finite prefix. `detectDivergence()`/`detectStall()` apply the same handling.

$$\mathrm{status} = \begin{cases}
\texttt{unknown} & N < 2\ \text{or input not an array}\\[2pt]
\texttt{diverging} & \Delta_{\mathrm{win}} > 0.02\\[2pt]
\texttt{improving} & \Delta_{\mathrm{win}} < -0.02\\[2pt]
\texttt{stalled} & |\Delta_{\mathrm{win}}| \le 0.02 \ \wedge\ (y_0 - y_{N-1}) < 0.1\\[2pt]
\texttt{converged} & |\Delta_{\mathrm{win}}| \le 0.02 \ \wedge\ (y_0 - y_{N-1}) \ge 0.1
\end{cases}$$

`detectDivergence()` and `detectStall()` are the same predicates exposed separately and are tested to
agree with the classifier. Note the total-drop test uses the **first** point of the (possibly
concatenated) history, not the maximum — a run that rose before falling is judged on
$y_{\mathrm{first}} \to y_{\mathrm{last}}$.

**Test coverage** (`__tests__/heuristics.test.js`), on synthetic 200–300 point series: a steady
$-0.01$/step ramp → `improving`; 150 steps of ramp followed by 150 flat steps → `converged`; a
$-0.0001$/step ramp (total drop 0.02 < 0.1) → `stalled`; a $-0.01$ ramp followed by a $+0.01$ ramp →
`diverging`.

**Watchdog payload (`watchdogStats`).** Five fields:
`{n_steps, first, last, min, recent_window_delta}`, plus `non_finite_steps` when the history holds
non-finite values. `n_steps` is the exact integer point count (never rounded); `first`, `last`, `min`
are 3 s.f. and `recent_window_delta` 2 s.f., each `null` when not finite (a non-finite latest value is
sent as `last: null`). The full history is *never* sent to the watchdog. Two implementation details
worth knowing:

- It **recomputes** `windowDelta()` over the finite values — a second `recentSlope` pass over the tail —
  instead of reusing `stats.recentSlopePerStep` from the `seriesStats()` call it already made.
- `max` and `recentSlopePerStep` from `seriesStats()` are deliberately dropped.

**Re-ask gate (`significantChange`).** Given the stats of the last LLM call and the current stats:

$$\mathrm{significant} \iff
(\,n_{\mathrm{steps}}^{\mathrm{new}} - n_{\mathrm{steps}}^{\mathrm{old}} \ge 200\,)
\ \vee\
\big|\,y^{\mathrm{new}}_{\mathrm{last}} - y^{\mathrm{old}}_{\mathrm{last}}\big| \ge \ln(1 + 0.02) \approx 0.0198$$

i.e. 200 new log lines, or a ≥2 % move in $\chi$. It returns `true` when there are no previous stats and
`false` when there are no current stats. **Non-finite fallback:** if either `last` is not a finite number
the $\ln(1.02)$ test is skipped entirely and the function returns the strict inequality
`prevStats.last !== nextStats.last`. Since `watchdogStats` sends a non-finite `last` as `null`, two
consecutive non-finite polls compare equal and do not force a model call. Unlike `WINDOW_DELTA_EPSILON` and `MIN_TOTAL_DROP` (module constants), the
two re-ask thresholds are per-call options (`{relativeDelta = 0.02, stepDelta = 200}`); no caller
overrides them — `useWatchdog` calls `significantChange(prev, next)` with no third argument.

---

### Step 15 — Watchdog scheduling and badge semantics

[`llm/watchdog/useWatchdog.js`](../../web_app/frontend/src/llm/watchdog/useWatchdog.js) runs only when
`settings.watchdogEnabled` is true (default **false** — the badge renders nothing otherwise) and the
history has ≥2 points.

- **`'off'` is a distinct status from `'unknown'`.** When the watchdog is disabled or the history has
  fewer than 2 points the hook returns the literal `{status: 'off', source: 'heuristic', note: null,
  lastCheckedAt: null}` — that is what the badge renders against. `'unknown'` is what
  `classifyConvergence()` returns for a too-short array and is only ever seen as the hook's *initial*
  state.
- **No timers of its own.** It observes the `rValueFile` prop, which the existing 3-second Live Data poll
  already refreshes. The effect is keyed off *content* — its dependency array is
  `[enabled, nSteps, lastValue, baseUrl, model, apiKey, watchdogIntervalMin, column]` — not object identity,
  because Live Data produces a fresh object on every poll. **Blind spot:** a history whose *length and
  final value* are unchanged but whose interior changed does not re-classify. During Live Data this is a
  real (if brief) staleness window.
- **The badge names its χ² term.** The history is the last `.log` column only (one fit term, e.g.
  `X_ray_(R)1`), so the badge reads `X_ray_(R)1: Improving` and its tooltip adds "χ² of X_ray_(R)1
  only, not the total" (`WatchdogBadge.jsx`, from `rValueFile.plotData.chiColumn`). Classifying on
  every χ² term or on the total is deferred beyond 0.6.0.
- **The badge status is always the heuristic.** `classifyConvergence()` runs on every change and sets the
  badge for free. The LLM's parsed `status` is discarded; only `parsed.note` (or, if the format was
  ignored, the first 200 characters of the raw reply) becomes the tooltip, with `source` flipping to
  `'llm'`. When the heuristic status *flips*, the existing LLM note is cleared to `null` and `source`
  reverts to `'heuristic'` until the next model reply lands.
- **The LLM branch needs a configured endpoint.** Beyond `watchdogEnabled`, it requires both `baseUrl`
  **and** `model` to be non-empty, plus no request already in flight.
- **Throttle.** A model call is made only when the heuristic status *changed* since the last call, **or**
  when `watchdogIntervalMin` (default **5** minutes) has elapsed *and* `significantChange()` is true. A
  single in-flight guard prevents overlap, and the attempt timestamp is recorded even on failure so a
  dead server is retried on the throttle schedule rather than every poll.
- **Failures are silent by design** — an unreachable server leaves the heuristic badge in place. The
  in-flight request is aborted **only on unmount** (`useEffect(() => () => abortRef.current?.abort(), [])`);
  nothing cancels it when the settings change mid-request.

---

### Step 16 — Transport

[`llm/provider/client.js`](../../web_app/frontend/src/llm/provider/client.js) is a hand-rolled
OpenAI-compatible client (no SDK, no SSE library).

- `listModels()` → `GET {base}/models`, returns `payload.data[].id`. `checkConnection()` wraps it and
  translates failures into `{ok: false, error, hint}` — a bare `TypeError` means "server down **or** CORS
  blocked" (indistinguishable from the page), flagged `unreachable: true`; its hint names both causes
  and points to the drawer's setup guide ("How to connect a model"), which the drawer opens after a
  failed Test and which carries the copy-ready `OLLAMA_ORIGINS` command with this page's
  `window.location.origin`, the LM Studio CORS step and the Safari note; 401/403/404/429 get specific
  messages. **One exception:** an
  `AbortError` (the caller's own abort signal) is **re-thrown**, not translated. The function's own
  comment says it "never throws"; the code disagrees, and the code is what runs.
- `streamChat()` parses the SSE body manually: split on `\n`, keep only lines starting with `data:`,
  buffer the incomplete tail across reads (chunks can split mid-line), stop at `data: [DONE]`, ignore
  unparseable lines. It yields deltas of several kinds; the reader is always cancelled in a `finally`.
  - `{content}`.
  - `{reasoning}`: reasoning models report chain-of-thought in `delta.reasoning` (Ollama) or
    `delta.reasoning_content` (DeepSeek and others). It is surfaced separately in the collapsible
    "Thinking" panel.
  - `{toolCall: {index, id, name, arguments, extraContent}}`: fragments of `delta.tool_calls`,
    assembled by `copilot/toolCalls.js`. `extraContent` is a call's `extra_content`; Gemini puts its
    `thought_signature` there and requires it back.
  - `{finishReason}`.

  `tools` / `tool_choice` join the body only when given. An HTTP failure carries its `status` (the
  Copilot's tools-unsupported test reads it).
- `completeChat()` is the non-streaming variant used by the watchdog.
- Rendering of replies is Markdown via `react-markdown` with `rehypeRaw → rehypeSanitize → rehypeKatex`;
  `<img>` is stripped from the sanitize schema so model output cannot trigger an external load
  ([`MarkdownMessage.jsx`](../../web_app/frontend/src/llm/components/MarkdownMessage.jsx)). Analysis
  results on the step cards are shown as `JSON.stringify` text in a `<pre>`, never as HTML.

---

### Step 17 — The Copilot tool loop (analyses, checks, budgets, verdict)

**The loop** ([`llm/copilot/agentLoop.js`](../../web_app/frontend/src/llm/copilot/agentLoop.js),
`runCopilotTurn()`):

1. The model is called with the tools (`tool_choice: 'auto'`).
2. If it returns tool calls, each call is checked and run in order:
   - its arguments are parsed (`''` → `{}`, a Markdown-fenced object accepted);
   - they are validated against the tool's schema subset (`type`, `enum`, `minimum` / `maximum`,
     `required`; numeric strings coerced; unknown keys dropped with a note);
   - the tool runs and its result is fitted to the budget.
3. The answers go back as `tool` messages, and the model is called again.
4. The turn ends when the model answers without tool calls.

What the loop guarantees:

- **Message rules.**
  - The assistant message that lists the calls carries `content: null` when its text is only
    whitespace (Anthropic rejects empty text blocks).
  - Its arguments are re-serialized as valid JSON (`'{}'` when unparseable); its `extra_content` is
    echoed; there is no `index`.
  - Every listed call gets exactly one `tool` message, including calls that were skipped, unknown or
    invalid.
- **Errors become information.** An unknown tool, invalid arguments or a thrown engine error
  becomes a `tool` message `{error}` that the model reads, so it can correct itself. The step shows
  as failed.
- **Duplicate calls.** A call identical to an earlier one in the same question (same validated
  arguments, key-sorted) returns the cached result with the note "identical call … change the
  arguments or answer".
- **Budgets** (`DEFAULT_LIMITS`):

  | Limit | Value |
  |---|---|
  | rounds of tool calls | 6 |
  | calls per round | 4 (extra calls are answered "skipped") |
  | one call | 120 s |
  | no new round after | 10 min |
  | one result | 3,000 chars |
  | all results of one question | 12,000 chars |
  | questions to the user (`ask_user`) | 2 per question, 1 per round; no time limit, and the wait does not count toward the 10 min |

  - A result over 3,000 chars has its largest array halved repeatedly, each cut recorded as
    `<key>_omitted`; a hard cut flagged `truncated: true` is the last resort.
  - Over 12,000 chars per question, the oldest results are replaced by
    `{compacted: …, checks}`. The newest round is never compacted, and checks are never trimmed.
- **When the budget is used**, a user-role note says so and the last request is sent with
  `tool_choice: 'none'`. Tool calls in that reply are ignored. If it brings no text (Ollama ignores
  `tool_choice`), a fixed "analysis budget used" answer closes the turn with
  `Outcome: not achieved`.
- **Models without tool support.** A 4xx whose message says the model or server "does not support"
  (or "is not supported" for) tools, tool calling, tool use or function calling retries the
  question once without tools, under `COPILOT_PROMPT_NO_TOOLS`. The model is then remembered for the
  session. Any other HTTP error is shown, including a schema 400 that merely mentions a function.
- **Stop / timeout.** Each tool run gets its own `AbortSignal`, aborted on Stop or timeout. The
  Copilot's worker (one of its own, not the pages' shared worker) is then terminated, and a Flask
  request is cancelled. `describeSymmetry` runs synchronously and cannot be interrupted (≈ 0.1 s on
  the demo).
- **One run per conversation.** The tools and the context are captured when the question is asked.
  Opening another run (`runKey` = directory + run id) stops the turn and starts a new
  conversation. In Flask mode a Live Data save mid-question is *not* a new run: a later analysis in
  that question reads the newer `.rmc6f`.

**The analyses.** Each tool returns `{data, checks}`: a summary built for the model, and verdicts
computed from the full engine result. A check is `ok`, `note`, `warn` or `fail`; `fail` means "does
not answer as asked". A tool is offered only when its inputs exist.

| Tool | Offered when | Engine | `data` (what is sent) | Checks (thresholds) |
|---|---|---|---|---|
| `read_run_data(section)` | the context has the section | `buildRunContext` | one context section as built (the builder's caps apply: ≤ 12 sites, 48 history points), before `contextToJson`'s trims | — |
| `assess_convergence()` | an R-value series exists | `classifyConvergence` (Step 14), `seriesStats` | column, status, n, first/last/min, total drop, least-squares slope per row over the last 10 / 25 / 50 % (≥ 3 points each) | note: judged on one `.log` column, not the total; ok: converged; note: improving; warn: stalled/diverging; warn: non-finite rows; warn: < 5 rows; fail: < 2 finite rows |
| `find_gr_peaks(pair, r_max = 6, max_peaks = 3)` | a parsed partials file | `extractPeaks` (Step 8), `firstShell` (Step 9) | peaks (r, FWHM); first-peak window `[r − FWHM, r + FWHM]`; `first_shell`: its window (onset − one grid step to the boundary), the peaks inside it, mean distance and sd from $r^2g$, the CN in both directions and the CN with the cut 0.2 Å earlier / later, the trough $\tilde g/h$ at the boundary; the pair's `MINIMUM_DISTANCES` value, matched order-free, and `g_at_constraint_fraction` | warn: the first shell is not separated from the next ($\tilde g/h > 0.02$ at its boundary: distances fill the gap, the CN depends on the cut); note: the shell holds several peaks (a split shell: use its window); note: no boundary below max(r_max, r₁ + 3 Å); warn: the closest-approach constraint cuts the distribution — raw $g$ in $(d_{\min}, d_{\min} + 2.5\,\Delta r]$ above 2 % of the peak (`CONSTRAINT_PILEUP_FRACTION`; demo Ga–Se 9 % at 2.2 Å), or the old $r - \mathrm{FWHM}/2 \le d_{\min} + 0.05$ Å; ok: peak, shell edge, mean, CNs; fail: no peak |
| `list_sites(element?)` | a browser-loaded run or Flask mode | worker `sites` / `/api/pca/sites` | per site: ref, label (composition if mixed), count, U_iso, rms axes, anisotropy, Mardia kurtosis; sorted by \|kurtosis\| then U_iso; top-3 refs by \|kurtosis\|, U_iso and anisotropy | warn: parse warning; note: reconstructed sites; warn: zero-spread sites; fail: no site (holding that element) |
| `analyze_displacement_directions(site_ref \| element, frequency?)` | same | worker `orientation` / `/api/pca/orientation`, `geometry: false` | used / rejected points; peak enhancement, count, expected and calibrated `peakSignificance`; peak direction as a Cartesian unit vector, as `[u v w]` (direct lattice, $M^{-T}d$ scaled to max \|component\| = 1, sign kept, since the map is not antipodally folded) and as angles to a, b, c; `mapSignificance`; antipodal asymmetry, its z and flag; orientation-tensor eigenvalues, anisotropy and its Bingham significance; the principal axis `[u v w]` (`crystalOrientationRows`) | warn: < 100 vectors; ok: peak ≥ 3σ; note: peak < 3σ; note: map test undefined; ok: map ≥ 3σ; note: tied peak cells; note: significant ± asymmetry; warn: parse warning; note: an element pooled over > 1 site (the site count from the basis) — where their local axes differ, as for the four vertices of a tetrahedral cluster, a per-site preference cancels, so an isotropic pooled map is not evidence against a per-site or cluster distortion |
| `compute_bond_angles(end1, apex, end2, r12_min, r12_max, r23_min?, r23_max?, bin_width = 1)` | same | worker `triplets` / `/api/triplets` (`APP_MAX_ANGLES` applies) | angle count, central atoms, mean ± sd; ≤ 3 modes (local maxima of a 5-bin running mean of the counts, each with the share of angles within ± 7.5°); the distribution in 18 bins of 10°, raw and sin θ-corrected ($S_j = (N_j/N)\,/\,[(\cos\theta_j - \cos\theta_{j+1})/2]$, 1 for random directions); for a window that does not cut the shell and a mean coordination within 0.5 of 6, the share and mean angle below and above 135° against an ideal octahedron (80 % cis at 90°, 20 % trans at 180°), within 0.5 of 4, the share within 15° of 109.47°; mean bond length, its sd from the engine's 40-bin length histogram, and unique bonds per window; coordination histogram and mean | fail: 0 angles ("widen the window"); warn: > 50 % of central atoms have no A–B bond in the window; otherwise warn: < 60 % of central atoms have the most common A–B bond count (the window may cut through a coordination shell); note: raw counts deplete 180° (with the octahedron split); warn: parse warning; ok otherwise. Refused before running: min ≥ max, an r23 bound without its pair |
| `detect_symmetry(tolerance_A ∈ [0.001, 1])` | the structure has a `basis` (browser-loaded runs) | `describeSymmetry` (the Detected SG card's finder), *not* changing the card | group, number (only when certain), point group and its facts (`point_group_facts`, as in the context), n ops, max residual, every orbit with `orbitLabel` Wyckoff label, site symmetry and representative fractional coordinates (4 decimals) | ok: group holds; note: the point group's facts ("take inversion, chirality and polarity from this"); note: undetermined; note: "≥" lower bound; warn: not analysed (caps); fail: no basis |
| `compare_fit(dataset, n_ranges = 6)` | a fit file with a model and a measured column (`xpdf`, `npdf`, `xray_sq`, `neutron_sq`, `bragg`; repeated kinds keyed `npdf_2`, …) | `fitTool` in `runTools.js` (roles by label `/exp\|obs/` vs `/calc\|rmc\|fit/`, else column order — `rwpColumns`, duplicated) | the dataset's $R$ (Step 10) over its range; per equal-width range of r or Q: $R$ there and its share of $\sum (y^\mathrm{calc}-y^\mathrm{expt})^2$; the 5 largest local extrema of the 5-point mean residual (position, model − data, relative to the data's rms) | ok (the card's headline): the overall $R$ and where the misfit lives — one range holding ≥ 40 % of the squared residual (`MISFIT_CONCENTRATION`), else spread; note: another range's own $R$ ≥ 3× the dataset's (`RELATIVE_MISFIT_FACTOR`) — weak signal poorly fitted, e.g. high-Q F(Q); fail: < 10 points |
| `analyze_clusters(element, r_min, r_max)` | a browser-loaded run with site/cell columns, or Flask mode | worker `clusters` / `/api/clusters` (Step 18) | the cluster kind (members, edges, count, templates; other sizes); the average cluster's sorted edges, spread, median edge sampling error and spread/SE; per statistic obs / null mean / z / change = obs/null − 1: the edge ladder, ladder width, shape rms, breathing sd, mean edge; tetrahedra: E and T2 rms and T2 trigonality; the fractions of long / short gyration axes within 15° of <100>, <110>, <111>; bonds outside the window, contacts between clusters per atom | fail: extended network; fail: no clusters; warn: other sizes hold > 10 % of clusters; note: the average cluster is distorted (spread ≥ 5 SE, `STATIC_SPREAD_SE`) or regular; ok: ladder width departs from the null (\|z\| ≥ 3 **and** \|change\| ≥ 5 %: more distorted — a correlated distortion — or less — rigid units) or does not; note: breathing departs (same rule); note: trigonality departs (\|z\| ≥ 3, \|Δ\| ≥ 0.02); note: an axis family enhanced (z ≥ 3, obs/null ≥ 1.1); warn: > 5 % of bonds outside the window; warn: > 0.05 inter-cluster contacts per atom |
| `analyze_modulations(element?)` | same | worker `modulations` / `/api/modulations` (Step 19) | per element: the zone-boundary points as `[q, ratio to the median, percentile]` rows, the three strongest modes the same way; the supercell, wavevector count and untestable points | ok: ordered displacements at a zone-boundary point (ratio ≥ 3 and top 1 % or strongest); ok: no cell doubling (the highest point named); note: a long-wavelength wave dominates (≥ 5×, every \|q_i\| ≤ 0.25); note: untestable points (odd axes); fail: no displacements |

**Show in \<page\>.** A tool may declare `page` and `pageParams(args)`. The loop computes the
parameters from the validated arguments (null after an error), and the step's button hands them to
the host: `onOpenPage(page, params)` → App `pageRequest` → the page's `request` prop, applied once
per click.

- `compute_bond_angles` → Bond Geometry `{end1, apex, end2, r12Min, r12Max, r23Min, r23Max,
  binWidth}`. The page fills its boxes (`tripletInputsFromRequest`) and computes when the boxes
  hold exactly that request, through the same engine, so its numbers equal the step's.
- `analyze_displacement_directions` with a `site_ref` → Displacement Directions
  `{referenceNumber, frequency}`, with frequency `'auto'` when the Copilot let the engine choose.
  The page's smoothing, weight and frame are left as the user set them: the drawn enhancement can
  differ (the Copilot's is unsmoothed), while the significance tests use raw counts. Element
  pooling has no page view, so that step only opens the page.
- Other steps open their page without parameters.

The Copilot never calls this itself; only the user's click does.

**View tools** (2026-10-10). Three tools change what is on screen and nothing else (`view: true`;
never answered from the duplicate-call cache):

- `highlight_cards(cards, marks?)` (module, `copilot/viewTools.js`) gives the host cards it names
  (`data-copilot-card`, catalog `src/copilotCards.js`) a breathing edge through the dock's focus
  controller (`focus/focus.js`): it opens the main card's page, waits ≤ 3 s for it to mount, and
  scrolls it into view; `marks` are x positions drawn on a plot card, which zooms to them
  (`markWindow`: the marks' span plus ¾ of it, ≥ 6 % of the range per side). It ends at the next
  highlight, a new question, a click outside the Copilot, or a collapse.
- `orient_3d_view(card, along?, center_site?, highlight_sites?, select_site?, reset?)` (host,
  `src/copilotViewTools.js`): `along` is read by `parseDirection` (`viewer3d.js`) and becomes the
  Cartesian line of sight $\hat d = (u\mathbf a + v\mathbf b + w\mathbf c)/\lVert\cdot\rVert$
  (a, b, c: the cell rows; PC1–3: the site's principal axes). The camera goes to
  target + $\hat d$·distance (the direction points at the viewer), screen-up the first of c, b, a
  more than 18° from $\hat d$ with its component along $\hat d$ removed. On the sphere's PCA
  frame $\hat d$ is projected onto the PC axes. A selected site is awaited (≤ 15 s) before the
  view turns, because a new site re-frames its view.
- `current_view()` reports the page, its cards, the highlight, and each 3D view's line of sight as
  the simplest integer $[uvw]$ (|index| ≤ 4) within 2° (`nearestLowIndex`), else its components
  scaled to a largest of 1 (marked ≈).

**Checks on the answer** (2026-10-11; `copilot/answerAudit.js`). A tool may return `facts` beside
`data` and `checks` (kept by the loop, not sent to the model): `assess_convergence` its status,
`detect_symmetry` its group with the point-group facts, `analyze_displacement_directions`
`{subject, significant (peak or map ≥ 3σ), peak_sigma, map_sigma}`. The run context adds the
watchdog's status and every ladder group (`contextAuditFacts`). When the model answers,
`auditAnswer(text)` splits it into sentences and clauses (a new clause at `;`, `:` and before
"but", "while", "whereas", "unlike", "although", "though", "however") and looks for positive
claims: one in a question, or with a negation (not, no, never, neither, nor, none, without,
lack, absence, cannot, non, yet, n't) in the 60 characters before it, does not count.

- *convergence:* "converged", "reached convergence", "at convergence" or "plateau" while the
  latest status is improving or diverging, unless the sentence also says improving, diverging or
  "still falling";
- *chirality / centrosymmetry:* "chiral" (not "achiral", "non-chiral") or "centrosymmetric" (not
  "non-centrosymmetric") in a clause whose named groups (matched as whole symbols, however written)
  all lack the property;
- *inversion lost:* "lose / break / remove … inversion (centrosymmetry, centre of symmetry)" or
  "… inversion is lost" when none of the groups the answer names (else the context's group) has one;
- *direction:* "preferred (preferential) direction (orientation, axis)" when the question's
  directions analyses found none at ≥ 3σ.

The findings (`{rule, text}`, one per rule) reach the turn as `audit` and show under the answer.
`lowerOutcome` then replaces a verdict of *achieved* or *partial* (or none) by **not verified**
(`claimed` keeps the model's, the note says what it claimed); *not achieved* stays.

**The Workbench's own questions** (2026-10-11; `copilot/clarify.js`). `buildShellCatalog` reads the
partials file: per ordered pair (centre, partner) the first shell as `firstShell` finds it (Step 9),
its window rounded to 0.01 Å, the peaks in it, the CN $4\pi\rho_\text{partner}\int r^2 g\,dr$ over
the window (linear interpolation between grid points), the first-peak window $[r_1 - w, r_1 + w]$
(w the FWHM) and its CN, and the CN with the window's top 0.2 Å earlier and later. A window
*changes the answer* when its CN differs from the whole shell's by ≥ max(0.5, 0.15 × CN)
(`CN_CHANGE_ABS`, `CN_CHANGE_REL`). `readQuestion` reads the user's question for its intent
(bond angles; coordination, polyhedra, neighbours; directions; space group), the run's elements
(element symbols that are English words — In, As, S… — only after "around", "of", "for" or
before "atoms", "sites"), pairs and triplets, a polyhedron (octahedron 6, tetrahedron 4, trigonal
prism 6, square planar 4, cuboctahedron 12), a length or tolerance in Å, and a site number.

- `beforeTurn(question)`, in order: an angle or coordination request with one element and no pair
  or triplet → which bonded partner (shells starting within 1.4× the nearest one's start,
  `BONDED_REACH`; at most 4 options; "Every neighbour type" for coordination); not asked when only
  one partner or when the named polyhedron's CN matches exactly one partner. A named pair or
  triplet (A–B–C: the B–A and B–C shells) → which window, when `ambiguousWindows` finds a split
  shell whose first peak changes the CN, or a shell whose trough is above 2 % of the peak and
  whose 0.2 Å cuts change it; not asked when the named polyhedron's CN matches exactly one option.
  A directions question with one element and no site number → one site per orbit or all pooled,
  when the element has ≥ 2 sites and an orbit whose site symmetry is not cubic (23, m-3, 432,
  -43m, m-3m). A space-group question with no length in Å → which group, when ≥ 2 ladder rungs
  reach ≥ 1.25× the noise floor (`NOISE_MARGIN`).
- `beforeCall(tool, args)`, for a tool with `shellWindows(args)`: a window overlapping its first
  shell whose CN is below the whole shell's by the margin is asked about (whole shell, or keep),
  unless both its bounds appear in the question (or in an answer to an earlier Workbench question).
- `resolveClarification`: the chosen option's note (appended to the user's question as
  `(The Workbench asked: "…" — …)`) or the chosen window's arguments; a skip takes the default; a
  typed answer is passed on verbatim, or for a window read as two numbers.

The loop asks only when tools are on, an `askUser` exists, and fewer than `maxQuestions` questions
were asked; each is an `ask_user` step with `app: true`.

**Asking the user** (`ask_user`, 2026-10-11; module, `copilot/askUser.js`). An *interactive*
tool (`interactive: true`): it is never answered from the duplicate-call cache, has no time limit,
and its wait is subtracted from the turn's clock. `normalizeAskRequest()` turns the arguments into
`{question, options: 2–4 × {label, description?}, allowOther, defaultIndex}` (plain strings read as
labels, `"label — description"` split at a spaced dash, labels deduplicated case-insensitively, the
default matched by label or by number 1–4) and throws on an empty question or fewer than two
distinct options. The tool calls `askUser(request, { signal })` from the tools' context
(`useCopilot`), which shows the question card and resolves with
`{answer, source: 'option' | 'other' | 'skipped', index}`; `answerResult()` turns that into the
tool's data and one check ("use exactly this choice", the typed answer verbatim, or "use your
default and say you assumed it"). A second call in the same round, or a third in the question, is
answered with an error and never shown. The turn's abort signal rejects the wait (Stop, a run
change, an unmount). `ASK_RULES` join the prompt when the tool is offered (`copilotSystemPrompt`).

**The context window** (2026-10-11; `provider/contextGuard.js`, `provider/ollama.js`). Every
Copilot request goes through `streamCopilotChat()`. A server is probed only when it is local or on
port 11434: `GET /api/version` answering a version is Ollama; `GET /api/v0/models` listing load
states is LM Studio; anything else is sent unchanged (`client.js`).

- **The estimate** (`estimatePromptTokens`): the characters of every message's content and tool
  calls plus the tool schemas' JSON, divided by `CHARS_PER_TOKEN` = 2.5, plus
  `MESSAGE_OVERHEAD_TOKENS` = 8 per message. Deliberately high: the Demo's first request is 28,382
  characters → 11,386 estimated, 8,437 measured (gemma4:26b).
- **Ollama** (`planContext`): needed $= 	ext{estimate} + 8192$ (`REPLY_RESERVE_TOKENS`). Refused
  when needed exceeds the model's maximum (`/api/show` `*.context_length`). Otherwise `num_ctx`
  = the loaded context (`/api/ps`) when it holds the request, else this conversation's earlier
  value when it still does, else $\max(16384, 8192\lceil 	ext{needed}/8192
ceil)$, capped at the
  maximum; `num_predict` = `num_ctx` − estimate. The request goes to the native `/api/chat`
  (NDJSON; messages converted: string content, object arguments, `tool_name` on tool results; no
  `tool_choice`). After the reply, `prompt_eval_count` below 0.4 × the estimate
  (`READ_SHARE_LIMIT`) stops the turn: the server cut the prompt anyway.
- **LM Studio**: when the model is loaded with `loaded_context_length` < estimate + 4096
  (`LMSTUDIO_REPLY_MIN_TOKENS`), the turn is refused with the context length to load it with.
- **Why** (measured on Ollama 0.40.2): a prompt longer than `num_ctx` returns 200, after Ollama has
  dropped whole messages and cut the rest to half the window keeping its first 5 tokens
  (`truncating input prompt limit=2051 prompt=5694 keep=5` at 4,096); the model then answers
  without the rules or the run.
- A reply that stops at the length limit (`done_reason` / `finish_reason` `length`) sets `cut` on
  the turn, and its footer says so.

**The verdict.** The prompt asks every answer to end with
`Outcome: achieved|partial|not achieved — <one sentence>`.
[`copilot/outcome.js`](../../web_app/frontend/src/llm/copilot/outcome.js) `parseOutcome()` reads
the line:

- the **last** such line wins;
- decoration is tolerated (bold or italics, inline-code backticks, a list bullet or heading
  marker, an en dash or a colon). The prompt shows the line in backticks, and some models copy
  them;
- `not_achieved` is accepted for "not achieved";
- the three statuses are matched strictly; any other word is not a verdict.

The line is removed from the answer text and shown as the chip. A reply without it has no chip. The
verdict is the model's own self-assessment. The checks are the deterministic part: read both.

---

### Step 18 — Cluster distortion against a shuffle null (`analyze_clusters`)

**Why.** A g(r) peak cannot tell a cluster distortion from thermal spread, and a displacement-direction
map pooled over a cluster's vertices cancels it (each vertex has its own local frame). The question a
lacunar-spinel referee asks — *are the M₄ tetrahedra distorted, how, and is that more than the atoms'
own motion?* — needs the clusters themselves. Reference: `cluster_distortion()` in
[`rmc_toolkits/clusters.py`](../../rmc_toolkits/clusters.py); browser port
[`workers/clusters.js`](../../web_app/frontend/src/workers/clusters.js), parity-tested against
`clusters_fixture.json` (`tests/generate_clusters_fixture.py`; integers exact, floats to $10^{-9}$; the
demo's 4000 clusters agree to $6\times10^{-11}$ relative).

**Input.** The per-site displacement clouds of the PCA engine (`load_site_displacements` /
`siteDisplacementsFromRmc6f`), which now keep each atom's cell indices; an element; a window
$[r_\min, r_\max]$ (Å, ≤ 8 at the boundaries). A coordinates-only file (no site/cell columns) is refused.

**Templates (average structure).** Sites whose label is the element are bonded when the distance
between site means, over unit-cell image shifts $|n_i| \le \lceil r_\max/w_i\rceil + 1$ ($w_i$ the
perpendicular widths), lies in the window (± `WINDOW_TOL`). A breadth-first walk from each unvisited site
assigns integer cell offsets; a site bonded to its own image, or reached twice at different offsets,
makes the network **extended** (chains or a net — no clusters). Each finite component is a template:
member sites in site order, their offsets, and its bonds as sorted member pairs. The dominant kind
(members, bonds) by cluster count is analysed; other kinds are listed.

**Clusters (configuration).** One cluster per template per cell $c$: member $m$ is the atom of site
$s_m$ in cell $c + o_m$ (mod the supercell); a cluster with a missing atom is skipped and counted.
Member positions are $\mathbf b_m + \mathbf u$, with $\mathbf b_m = (\mathbf f_{s_m} + \mathbf o_m)\,A$ the
template geometry and $\mathbf u$ the atom's displacement from its site mean — so a bond's length is
measured however far it stretches.

**Statistics** (pooled over the dominant kind's templates; $L_{ce}$ = length of bond $e$ in cluster $c$):

| Key | Definition |
|---|---|
| `meanEdge` | $\langle L\rangle$ over all bonds |
| `ladder` | mean over clusters of each cluster's $k$-th shortest bond |
| `ladderWidth` | $\langle \max_e L_{ce} - \min_e L_{ce}\rangle_c$ |
| `shapeRms` | $\sqrt{\langle (L_{ce} - \bar L_c)^2\rangle}$ — shape distortion |
| `breathingSd` | population sd of $\bar L_c$ — size fluctuation |
| `eRms`, `t2Rms` (4 members, 6 bonds) | with opposite pairs $(e,e')$ = (01,23), (02,13), (03,12): sums $p_k$ and differences $d_k$; $\sqrt{\langle\sum_k (p_k-\bar p)^2\rangle/2}$ (E, tetragonal-type) and $\sqrt{\langle\sum_k d_k^2\rangle/3}$ (T2) |
| `t2Trigonality` | $\langle \sum d_k^4 / (\sum d_k^2)^2\rangle$ over clusters with $\sum d_k^2 > 0$: $1/3$ when the three opposite pairs differ alike (three bonds at one vertex long or short — trigonal, C₃ᵥ), 1 when one pair alone does (one bond against its opposite); a uniformly random direction averages 3/5 |
| `longAxis`, `shortAxis` | fraction of clusters whose gyration-tensor axis of largest / smallest eigenvalue (`eigh`; the port's Jacobi `eigenDecomposition`) is within 15° (`AXIS_CONE_DEG`) of a lattice direction of the family <100> (a, b, c), <110> (a ± b, …) or <111> (a ± b ± c, 4); for a dimer the bond itself is the long axis |

**The shuffle null.** For shuffle $k = 0 … 23$ (`DEFAULT_SHUFFLES`), each member $m$ of template $t$
takes its displacement from the cluster $\pi_{k,t,m}(c)$ instead of $c$, with $\pi$ the stable argsort
of the keys $\mathrm{mix32}(h_{k,t,m} \oplus c)$, $h$ chained lowbias32 hashes of the seed (`0x5EED`),
$k+1$, $t+1$, $m+1$ — the same integers in both engines. This keeps every site's displacement
distribution (spread, anisotropy, non-Gaussianity, bimodality) and destroys only the correlation between
the atoms of one cluster. RMC moves atoms one at a time, so without data a configuration *is* this null;
a departure is a correlation the data put in. Each statistic is reported as observed, null mean, null sd
($n-1$) and $z = (\text{obs} - \text{mean})/\text{sd}$ (`null` when sd = 0), elementwise for the
ladder and the axis families. With thousands of clusters a 1 % effect reaches $z > 3$, so the Copilot's
checks also require an effect size (Step 17).

**The average cluster.** The template's own bonds $|\mathbf b_j - \mathbf b_i|$ (a *static*
distortion shows here, not in the shuffle test, which keeps the site means) with each bond's sampling
error $\sqrt{\sum_{s\in\{i,j\}} \langle (\mathbf u\cdot\hat e)^2\rangle_s / N_s}$ from the two sites'
displacement spreads along it; `ideal` gives the sorted bonds averaged over templates, the largest
spread, the median sampling error and the largest spread / error. Six independent normal bonds span
~2.5σ on average and exceed 5σ in under 1 % of cases.

**Outliers.** `brokenEdgeFraction` — bonds of the analysed clusters whose length in the configuration
lies outside the window; `extraContacts` — pairs of the element within the window (the triplets
engine's linked-cell count over every atom of the element, positions rebuilt as
$(\mathbf f_s + \text{cell})/N + \mathbf u A^{-1}$) minus the template bonds inside it: atoms crossing
the gap between clusters; per atom, $2\times$ that over the element's atom count.

**Synthetic checks** (`tests/test_clusters.py`, 6×6×6 cells of an Nb₄ tetrahedron, a Ga₂ dimer and a
Se): independent noise — every \|z\| < 3; one random vertex per cluster pushed 0.3 Å out along its
3-fold axis — average cluster regular, trigonality z ≈ −8, long axes along <111> z ≈ 9, breathing
lower than the null's binomial count of pushed vertices; the same vertex pushed everywhere — average
cluster spread ≈ 45 SE, local statistics at the null.

**The demo (GaTa₄Se₈, 250 K, Ta 2.65–3.33 Å, 4000 clusters).** Ladder width 0.295 Å against 0.386 Å
in the null (−23 %, z = −64): the RMC clusters are *more rigid* than independent atoms; breathing +27 %;
short axes along <111> in 19 % of clusters against 14 % (z = 11), long axes along <110> 23 % against
20 % (z = 4) — the axis of the monoclinic distortion small-box PDF fits report for this compound
([Yang et al., Phys. Rev. Research 4, 033123 (2022)](https://doi.org/10.1103/PhysRevResearch.4.033123));
average cluster 3.031–3.054 Å, 5.5 SE.

**Cost.** Templates: (sites of the element)² × image shifts. Clusters and each shuffle: linear in the
atoms of the element. The demo: 1.5 s in Python (load included), 0.7 s in the port after a 0.2 s parse.

---

### Step 19 — Displacement modulations: order that enlarges the cell (`analyze_modulations`)

**Why.** The average structure (and so `detect_symmetry` and the noise floor of Step 5) folds every
cell into one: it cannot show order whose period is longer than one cell — the cell doubling of
GaNb₄Se₈ below 29 K (P2₁2₁2₁, $c' = 2c$; Yang et al., Phys. Rev. B 109, 144101 (2024)) or of
GaTa₄Se₈ below 50 K, tilts, charge-density waves. The configuration can. Reference:
`displacement_modulations()` in [`rmc_toolkits/modulations.py`](../../rmc_toolkits/modulations.py);
port [`workers/modulations.js`](../../web_app/frontend/src/workers/modulations.js) (a separable DFT
for numpy's `fftn`), parity-tested against `modulations_fixture.json`
(`tests/generate_modulations_fixture.py`; the 5 K run agrees to $10^{-15}$).

**Method.** For each reference site $s$ of an element, its displacements $\mathbf u_s(\mathbf c)$ from the
site mean are placed on the $N_1\times N_2\times N_3$ cell grid (a missing copy is a zero, and
counted) and transformed,

$$\mathbf U_s(\mathbf q) = \sum_{\mathbf c} \mathbf u_s(\mathbf c)\,e^{-2\pi i\,\mathbf q\cdot\mathbf c},\qquad
\mathbf q = (h/N_1, k/N_2, l/N_3),\qquad
P_s(\mathbf q) = \frac{|\mathbf U_s(\mathbf q)|^2}{n_s\,\mathrm{tr}\,C_s},$$

so independent displacements give $\langle P_s\rangle = 1$ at every $\mathbf q \ne 0$ ($C_s$ the
site's population covariance). The element's power is the mean over its sites (zero-spread sites
skipped). Every $\mathbf q \ne 0$ of the grid is the comparison set — the configuration's own
diffuse background, which keeps RMC's long-wavelength correlations a white-noise null ignores (on
both real runs a χ² test against independent cells gave χ²/dof ≈ 0.4–0.9, i.e. *less* zone-boundary
power than white noise, because neighbouring cells move together). Reported per element:

- the **zone-boundary points** with components in $\{0, \tfrac12\}$ (seven; only those whose halved
  axes have an even $N_i$, the rest listed as `untestable`);
- the five **strongest modes**, $\pm\mathbf q$ counted once as the member with the smaller grid
  index (their powers agree to rounding; folded components in $(-\tfrac12, \tfrac12]$);

each with its ratio to the median power over $\mathbf q\ne0$ and its percentile (the share of
$\mathbf q \ne 0$ with power below $P(1 - 10^{-9})$, so Friedel partners count alike in both engines).

**Checks** (`summarizeModulations`): ok — ordered displacements at a zone-boundary point when its
ratio ≥ 3 (`MODULATION_RATIO`) **and** it is in the top 1 % or is the element's strongest mode (a box
with < 100 wavevectors cannot reach the 99th percentile); otherwise ok — no cell doubling, naming the
highest point; note — an element's strongest mode within 0.25 of the zone centre on every axis at
≥ 5× the median is a long-wavelength (box-scale) wave, not a superstructure; note — untestable points.
A planted doubling (Nb ± 0.1 Å along z in alternate cells, 6³ cells) gives 204× at (0,0,½), the
strongest of 215 wavevectors. Both real runs: no zone-boundary point above the 85th percentile; the
strongest modes are box-scale waves — GaNb₄Se₈ 5 K: q = (0.1, −0.2, −0.2) for all three elements
(5–14×); GaTa₄Se₈ 250 K: (0.1, 0.2, −0.1) and (0.2, 0, 0) (13–33×). The 5 K configuration therefore
carries neither the orthorhombic cell doubling nor a cluster distortion beyond independent atoms
(Step 18), where diffraction finds both: the RMC started from, and stayed near, the cubic model.

**Cost.** One 3-D transform per site and component: 0.02 s (numpy) and 0.04 s (port) for 52 sites
on a 10³ grid.

---

### Runtime-mode differences (static vs Flask) — read this before quoting the assistant

The context is far richer in **static mode** (the browser parses the run itself) than in **Flask mode**,
because the Flask endpoints return figures and metadata rather than series and basis data. What survives
in each mode:

| Context block | Static mode (browser parsing) | Flask mode (backend) |
|---|---|---|
| `run`, `live_mode` | yes | yes |
| `structure` (cell, angles, composition) | yes | yes (`/api/structure`) |
| `configuration_optimization` | yes | yes — `/api/structure` returns the header's move counters (`moves`) |
| `run_settings` | yes | **no** — the `.dat` is read only from a locally picked folder |
| `symmetry` (+ per-site displacements, `frac`) | yes | **no** — `/api/structure` returns no `basis`, so `describeSymmetry()` returns `null` |
| `pair_correlations` | yes (needs parsed partials) | **no** — plot files carry no `plotData.series` |
| `datasets` | `kind` + `title` + `rwp` (never `stog`) | `kind` only (never `stog`) |
| `convergence` + watchdog badge | yes | **no** — no parsed R-value series, so the badge stays hidden |
| `pca_displacements` | after opening the PCA Ellipsoid page | after opening the PCA Ellipsoid page (`/api/pca/sites`) |
| `partial_weights` | yes (element counts + the `.dat` blocks or plot kinds) | from the plot kinds only (no `.dat`) |
| `analyze_clusters`, `analyze_modulations` | yes (worker; needs site/cell columns) | yes (`/api/clusters`, `/api/modulations`) |
| `compare_fit`, `find_gr_peaks` | yes | **no** — no parsed series |

Two cross-cutting notes: STOG datasets (`.gr`, `.sq`, `.fq`) are filtered out by `isDashboardPlotFile`
before the assistant sees anything, in **both** modes; and the Dashboard visibility toggle affects
`convergence`/watchdog but not `datasets`/`pair_correlations` (Step 1).

A code comment in `symmetryContext()` says the block also works from a minimal `{spaceGroup}` object "in
Flask mode"; that code path exists and is tested, but no caller currently supplies it — in Flask mode
`Dashboard.jsx` passes `symmetry = null` and the block is absent entirely.

---

### Parameters and defaults

| Name | Value | Units | Where | Meaning |
|------|-------|-------|-------|---------|
| `CONTEXT_CHAR_BUDGET` | 6000 | characters | `runContext.js` | Target size of the serialized context (`compactJson`) |
| `HISTORY_POINTS` | 48 | points | `runContext.js` | Convergence history cap; also `downsampleSeries`'s default argument |
| history trim ladder | 24, then 12 | points | `contextToJson()` | Successive re-downsamples of the *already downsampled* 48-point array |
| `MAX_SITES` | 12 | sites | `runContext.js` | Cap on symmetry sites and PCA sites at build time |
| budget site caps | 6 (PCA), 8 (symmetry), 12 (datasets) | entries | `contextToJson()` | Caps applied only under budget pressure |
| `roundSig` default | 3 | significant figures | both context files | Global rounding; 2 s.f. for displacements, residuals, tolerances, slope; 4 s.f. for cell lengths/angles; `n_steps`, `supercell`, `total_atoms` and the three move counters are **not** rounded |
| chi floor | $10^{-12}$ | raw $\chi$ (before the natural log) | `browserData.js` | $\chi \le 10^{-12}$ maps to $\ln(10^{-12}) \approx -27.6$, not $-\infty$ |
| recent-window length | $\max(\min(N,5), \lceil 0.2N\rceil)$ | points | `runContext.js`, `heuristics.js` | Tail used for the slope ($=N$ when $N<5$) |
| `WINDOW_DELTA_EPSILON` | 0.02 | $\ln\chi$ (≈2 % in $\chi$) | `heuristics.js` | Flat-vs-moving threshold on $\Delta_\mathrm{win}$ |
| `MIN_TOTAL_DROP` | 0.1 | $\ln\chi$ (≈10.5 % in $\chi$) | `heuristics.js` | Converged-vs-stalled threshold on the total drop |
| `significantChange.relativeDelta` | 0.02 | fraction of $\chi$ | `heuristics.js` | Re-ask threshold, compared as $\ln(1+0.02)$; a call-site **option**, never overridden |
| `significantChange.stepDelta` | 200 | log lines | `heuristics.js` | Re-ask threshold on new points; also a call-site option |
| `watchdogIntervalMin` | 5 (editable 1–120, step 1) | minutes | `settings.js`, `ConnectionSettings.jsx` | Minimum spacing between watchdog LLM calls; clamped to ≥1 on save |
| `watchdogEnabled` | `false` | — | `settings.js` | Badge is off until enabled |
| `rMax` | 6 | Å | `pairCorrelations.js` | Upper $r$ limit for peak finding |
| `maxPeaks` | 2 | peaks | `pairCorrelations.js` | Peaks kept per pair (first two by $r$) |
| smoothing width | 5 | grid points | `pairCorrelations.js` | Moving-average window (truncated at **both** edges: divisors 3/4/5/4/3) |
| peak floor | 0.15 | × window max | `pairCorrelations.js` | Minimum smoothed height for a peak |
| coincidence cutoff | $10^{-3}$ | Å | `pairCorrelations.js` | Distances below this are ignored |
| image shell | $\pm1$ | cells | `pairCorrelations.js` | Periodic images searched for nearest neighbours |
| `maxPairs` | 10 | pairs | `pairCorrelations.js` | Cap on `pair_correlations` entries |
| `CHAT_HISTORY_TURNS` | 8 | messages | `prompts/copilot.js` | Memory messages re-sent (questions + final answers, in-memory only) |
| `maxRounds` | 6 | rounds | `copilot/agentLoop.js` | Model requests with tool calls per question |
| `maxCallsPerRound` | 4 | calls | `copilot/agentLoop.js` | Extra calls are answered "skipped" |
| `toolTimeoutMs` | 120 000 | ms | `copilot/agentLoop.js` | One analysis; its signal is aborted |
| `maxTurnMs` | 600 000 | ms | `copilot/agentLoop.js` | No new round of tools after this |
| `resultCharBudget` | 3 000 | characters | `copilot/agentLoop.js` | One `tool` message |
| `turnCharBudget` | 12 000 | characters | `copilot/agentLoop.js` | All `tool` messages of one question; older ones compacted to their checks |
| `maxQuestions` (`MAX_QUESTIONS`) | 2 | questions | `copilot/askUser.js`, `copilot/agentLoop.js` | `ask_user` calls per question (one per round), the Workbench's included |
| `CN_CHANGE_ABS`, `CN_CHANGE_REL` | 0.5, 0.15 | atoms, of the shell's CN | `copilot/clarify.js` | A window changes the answer when its CN differs from the whole shell's by at least the larger |
| `BONDED_REACH` | 1.4 | × the nearest shell's start | `copilot/clarify.js` | Partners offered for "around <element>" |
| `NOISE_MARGIN` | 1.25 | × the noise floor | `copilot/clarify.js` | A ladder group must hold to this to make the tolerance a question |
| `POLYHEDRON_SPLIT_DEG` | 135 | ° | `copilotTools.js` | Splits cis from trans for a whole shell of 6 ± 0.5 neighbours (demo TaSe₆: 80 % at a mean 89°, 20 % at 158.8°) |
| `MIN_OPTIONS`, `MAX_OPTIONS` | 2, 4 | options | `copilot/askUser.js` | Choices a question offers |
| `CHARS_PER_TOKEN`, `MESSAGE_OVERHEAD_TOKENS` | 2.5, 8 | chars / token, tokens | `provider/contextGuard.js` | The prompt estimate (errs high) |
| `REPLY_RESERVE_TOKENS` | 8 192 | tokens | `provider/contextGuard.js` | Room for the reply in an Ollama `num_ctx` |
| `CONTEXT_STEP`, `MIN_CONTEXT` | 8 192, 16 384 | tokens | `provider/contextGuard.js` | An asked-for `num_ctx` is a multiple of the step, at least the minimum |
| `LMSTUDIO_REPLY_MIN_TOKENS` | 4 096 | tokens | `provider/contextGuard.js` | Room an LM Studio context must leave beyond the prompt |
| `READ_SHARE_LIMIT` | 0.4 | of the estimate | `provider/contextGuard.js` | Ollama `prompt_eval_count` below this: the prompt was cut; the turn stops |
| `copilotTools` | `true` | — | `settings.js`, `ConnectionSettings.jsx` | *Let the Copilot run analyses* |
| `MIN_ORIENTATION_POINTS` | 100 | vectors | `copilotTools.js` | Below: "the map is noisy" |
| `SIGNIFICANCE_SIGMA` | 3 | σ (calibrated, one-sided) | `copilotTools.js`, `prompts/copilot.js` | Peak / map significance threshold for the checks |
| `UNBONDED_APEX_FRACTION` | 0.5 | fraction | `copilotTools.js` | Central atoms with no A–B bond in the window above this: warn |
| `MODAL_COORDINATION_SHARE` | 0.6 | fraction | `copilotTools.js` | Central atoms with the most common A–B bond count below this share: warn (a spread coordination, as when the window cuts through a shell; demo Se–Ta–Se at 2.3–2.7 Å: 0.47, whole shells 0.83–0.96) |
| `SHELL_TROUGH_FRACTION` | 0.2 | of the first peak | `pairCorrelations.js` | A deeper local minimum closes the first shell; a shallower dip is inside a split shell |
| `SHELL_ONSET_FRACTION`, `SHELL_EMPTY_FRACTION` | 0.01, 0.001 | of the first peak | `pairCorrelations.js` | Shell onset; a tail below this closes the shell |
| `CONSTRAINT_PILEUP_FRACTION` | 0.02 | of the first peak | `runTools.js` | $g$ just above a closest-approach constraint above this: the constraint cuts the distribution |
| `SHELL_OVERLAP_FRACTION` | 0.02 | of the first peak | `runTools.js` | $\tilde g$ at the shell boundary above this: shells not separated |
| `MISFIT_CONCENTRATION`, `RELATIVE_MISFIT_FACTOR` | 0.4, 3 | share; × the dataset's R | `runTools.js` | `compare_fit`'s concentration and weak-signal notes |
| `DEFAULT_SHUFFLES`, `SHUFFLE_SEED`, `AXIS_CONE_DEG` | 24, `0x5EED`, 15° | — | `clusters.py`, `clusters.js` | The cluster null and the axis families |
| `MAX_CLUSTER_WINDOW` | 8 | Å | `app.py`, `pcaKdeWorker.js` | Largest `rMax` of a cluster request |
| `CLUSTER_Z`, `CLUSTER_EFFECT`, `TRIGONALITY_EFFECT`, `AXIS_ENHANCEMENT` | 3, 5 %, 0.02, 1.1 | — | `copilotTools.js` | A cluster statistic is a finding only when both its z and its size clear these |
| `STATIC_SPREAD_SE` | 5 | sampling errors | `copilotTools.js` | The average cluster counts as distorted above this spread |
| `CLUSTER_OUTLIER_FRACTION` | 0.05 | fraction / per atom | `copilotTools.js` | Bonds outside the window, contacts between clusters |
| `MODULATION_RATIO`, `MODULATION_PERCENTILE` | 3, 99 | × median, percentile | `copilotTools.js` | A zone-boundary point is order above both (or as the strongest mode) |
| `LONG_WAVE_RATIO`, `LONG_WAVE_Q` | 5, 0.25 | × median, r.l.u. | `copilotTools.js` | A strongest mode this strong and this close to the zone centre is a box-scale wave |
| constraint margin | 0.05 | Å | `copilot/runTools.js` | First-peak low side within this of the closest-approach constraint: warn |
| `temperature` (chat) | 0.2 (editable 0–2, step 0.1) | — | `settings.js`, `ConnectionSettings.jsx` | Sampling temperature |
| `temperature` (watchdog) | 0 | — | `useWatchdog.js` | Fixed |
| `baseUrl` | `http://localhost:11434/v1` | — | `settings.js` | Default endpoint (Ollama) |
| `model`, `apiKey` | `''` (empty) | — | `settings.js` | `model` may be auto-set to `models[0]` by the mount probe |
| symmetry tolerance | 0.2 **initial**, then user-set | Å | `App.jsx` (`SymTolContext`), `ModelSummary.jsx` | Shared "Detected SG" tolerance; changes the whole `symmetry` block |
| ladder max tolerance | 1.0 | Å | `Dashboard.jsx` | `toleranceLadder(structure, 1.0)` |
| PCA cluster distance | 1.5 (slider 0.4–2.5, step 0.1) | Å | `PcaKdePage.jsx` | Only for reconstructed sites; re-numbers `ref` when moved |
| `DEGENERATE_RATIO` | $10^{-6}$ | — | `pcaKde.js`, `pca_kde.py` | $\lambda_3/\lambda_1$ below this flags a planar/linear cloud |
| Rwp acceptability band | 0.01–0.10 | — | `prompts/system.js` | **Prompt assertion only** — no code enforces or cites it |
| `accepted_moves_per_atom` rubric | single digits = lightly sampled; tens+ = well sampled | moves/atom | `prompts/system.js` | **Prompt assertion only** |
| ladder-onset heuristic | onset tolerance ≈ mean static displacement | Å | `prompts/system.js` | **Prompt assertion only** |

---

### Caveats / what this is not

- **The assistant sees summary statistics, never your data.** It cannot re-fit anything, cannot inspect a
  curve, and cannot verify its own claims. Any statement it makes beyond the numbers in the context block
  is unsupported — the system prompt says as much, but a prompt is not a guarantee. Verify against the
  plots.
- **"No atom coordinates" is not quite true.** Per-atom rows never leave the device, but the
  representative fractional coordinate of every reported Wyckoff orbit (`symmetry.sites[].frac`, up to 12
  triples at 3 s.f.) does. Those are average-structure positions, not instantaneous ones — but they are
  coordinates.
- **Cloud providers receive the context, and the analysis results.** Selecting OpenAI, Gemini, or
  Anthropic sends the summarized run context (including run/dataset names, composition, and the
  `frac` triples above) to that provider under your API key. With analyses on, the results of every
  analysis the model runs go too: Step 17, which includes every orbit's coordinates from
  `detect_symmetry`. Turn the analyses off to send the context only. Merely opening the page issues an authenticated `GET /models` before you type
  anything. This is opt-in and warned about in the UI, but it is a real disclosure. Local providers keep
  everything on the machine.
- **The context inspector is only available before your first question.** It lives in the empty
  state, so once a conversation starts you cannot re-read the context without *New conversation*
  (which clears the conversation). Analysis results stay inspectable on their step cards.
- **The Outcome chip is the model's self-assessment.** A small model may write "achieved" over a
  failed check. The step cards' checks are computed, not generated: read them. It also re-serializes the *current* context, which under Live Data may differ from
  what accompanied an earlier answer.
- **The peak finder is a cue extractor, not a PDF refinement.** No baseline model, no profile function,
  no interpolation, no uncertainty. The reported `r` is the grid abscissa of the *smoothed* local
  maximum — quantised to the file's $r$ step and pulled toward the heavier flank of an asymmetric peak —
  and on a flat top it reports the plateau's trailing edge. The FWHM is a grid-snapped, zero-baseline
  half-height walk and is an overestimate whenever peaks overlap or the curve sits above a nonzero
  background. Peaks near the 6 Å cut are additionally tested and widened against edge-truncated
  smoothing averages. It looks only below 6 Å and keeps only the first two peaks by distance. Never quote
  these positions as refined peak centres or these widths as measured correlation widths.
- **Nearest-neighbour distances are average-structure distances**, computed from the circular-mean basis
  sites of the folded configuration over a $\pm1$ image shell, not from a coordination-shell analysis of
  the instantaneous configuration. For an element with a single basis site, the A–A entry degenerates to
  the shortest conventional-cell lattice translation — a lattice repeat, not a bond — and nothing in the
  context flags it.
- **`mean_disp_A` / `max_disp_A` mix static disorder with thermal motion**: they are the single-snapshot
  rms displacement of each reference site about its mean (Cartesian, through the full cell metric).
- **The whole `symmetry` block is tolerance-dependent, and the tolerance is a live UI control.** Space
  group, point group, operation count, residual, and the entire orbit partition (hence which sites exist
  and how their displacements are averaged) change with the shared "Detected SG" tolerance. Only
  `tolerance_A` records the value used. Two sessions on the same run can get different answers.
- **PCA sites may be reconstructed.** If the configuration has no reference-site columns, the site
  partition is a distance clustering at a user-set threshold, and `ref` numbering changes when that
  slider moves. The context carries neither the `reconstructed` flag nor the threshold, so the model
  cannot tell.
- **A PCA row can silently degrade to `{ref, element}`.** Missing or non-finite `U_iso_A2`,
  `rms_axes_A`, `anisotropy` are dropped by JSON serialization with no marker.
- **`degenerate` means planar/linear, not isotropic** — $\lambda_3/\lambda_1 < 10^{-6}$, a rank-deficient
  cloud whose ellipsoid statistics should not be trusted.
- **The convergence "steps" axis is log lines, concatenated across restart logs**, not RMCProfile move
  counts, and the concatenation assumes the logs sort chronologically. While any sibling log is mid-parse
  (routine during Live Data), the concatenation is skipped and only the first log is described.
- **The system prompt actively mislabels that axis** as "accepted-move steps" on every request. Discount
  anything the model says about "moves" that is derived from the convergence block.
- **`rwp` is unweighted** (Step 10) despite its name and the system prompt's description, and exists for
  only five dataset kinds. It is normalized by the experimental column; when it is undefined (no finite
  pair, or an all-zero experiment) it is `null` and the key is simply absent from the dataset entry.
- **The convergence series is one χ² term, not the total** (Step 11): the last `.log` column, named in
  `convergence.column`. A run whose other terms stall while that one improves is still reported as
  improving.
- **The character budget is best-effort.** Several blocks are never trimmed; an unusual run can exceed
  4,500 characters. And only four of the eight trim steps leave an `*_omitted` counter — a dropped second
  $g(r)$ peak, a re-downsampled history, and a dropped PCA `note` are invisible in the delivered JSON.
- **Downsampling is decimation, not averaging.** A spike between kept indices disappears from the history
  shown to the model, though it still affects `min`/`max` (computed on the full series). Under budget
  pressure the 24- and 12-point traces are decimations *of the 48-point trace*, not fresh resamples of
  the source.
- **`seriesStats()` uses `Math.min(...values)` / `Math.max(...values)`**, a spread over the entire
  history. Tests exercise 50,000 points; a much longer history could exceed the JS engine's argument
  limit and throw. No guard exists. The same unbounded-spread pattern appears twice more, though on
  much smaller arrays: `Math.max(...disps)` over one Wyckoff orbit's basis members in `symmetryContext()`
  (bounded by the number of reference sites in the conventional cell, **not** by the supercell atom
  count, because `structure.basis` holds one entry per reference number), and `Math.max(...smooth)` over
  the peak-finding window in `extractPeaks()` (bounded by the number of grid points below 6 Å).
- **The watchdog thresholds are heuristics, not statistics.** There is no noise model, no significance
  test, and no dependence on the number of datasets or the acceptance rate. `stalled` vs `converged` is
  decided by a single 0.1-in-$\ln\chi$ total-drop rule that a run resumed from an already-good
  configuration will fail by construction (it never dropped 10 % *in this log*). The watchdog also
  re-classifies only when the history length or its last value changes.
- **The LLM never overrides the badge.** If the model is wrong, the status chip is still the heuristic;
  only the tooltip sentence comes from the model, and a malformed reply falls back to the first 200
  characters of whatever it said.
- **Chat history is ephemeral and can be lopsided.** Turns live in component state, are cleared by a
  reload, and an aborted or reasoning-only reply leaves your question in the re-sent history with no
  answer beside it. Reasoning text is never re-sent.
- **Three of the system prompt's numeric rules have no code behind them** (Step 13): the 0.01–0.10 Rwp
  band, the accepted-moves-per-atom rubric, and the ladder-onset ≈ mean-static-displacement heuristic.
  They are author assertions in a prompt string, not validated thresholds.
- **No Python equivalent exists for this pipeline.** Nothing in this section is computed in
  `rmc_toolkits/`. Two upstream helpers do exist in both languages, and they do **not** agree equally
  well:
  - `rwp` — `browserData.js` and `rmc_toolkits/parsers.py` implement the identical unweighted formula,
    the identical column-role resolution (`rwpColumns` / `rwp_columns`) *and* the identical `null`/`None`
    sentinel for the undefined cases. **Exact match** (to floating-point round-off).
  - `detect_plot_kind` — the `pdf_partials` (`*PDFpartials*.csv`), `npdf`, `xpdf`, `xray_sq`,
    `neutron_sq`, `bragg`, `exafs_q`, `exafs_r` and `r_value` (`*-NN.log`) branches match exactly, which
    covers everything this section depends on. The **`stog` branch differs**: the JS returns `'stog'` for
    any name matching `/\.(gr|sq|fq)$/i`, while the Python returns `'stog'` only for the three literal
    default names `{scale_ft.gr, scale_ft.sq, scale_ft_rmc.fq}`. (This has no effect on the assistant,
    which never sees `stog` files at all.)
  - `DEGENERATE_RATIO` in the PCA engines (Step 6) also matches exactly between `pcaKde.js` and
    `pca_kde.py`.
