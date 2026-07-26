# AI Assistant — algorithm reference

> Part of the [Algorithms and Math Reference](../ALGORITHMS.md). Every step is anchored to the source; if this document and the code disagree, **the code wins**.

What the assistant computes from the loaded run before any model is called — and exactly what data leaves the device.

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
  - [Runtime-mode differences (static vs Flask) — read this before quoting the assistant](#runtime-mode-differences-static-vs-flask--read-this-before-quoting-the-assistant)
  - [Parameters and defaults](#parameters-and-defaults)
  - [Caveats / what this is not](#caveats--what-this-is-not)

---

## AI Assistant — run context construction, pair correlations, and convergence heuristics

### What this page shows

The **AI Assistant** page is a chat window over the run that is currently loaded in the Workbench, plus
an optional **convergence watchdog** chip that appears on the Dashboard's R-value card. Neither feature
sends your run files anywhere. What the app does is *compute a compact summary* of the run — a JSON
object targeted at ~4,500 characters (a best-effort budget, not a hard cap — see Step 12) — and hand
that summary, plus your question, to whatever
OpenAI-compatible chat endpoint you configured. By default that endpoint is a server on your own
machine (`http://localhost:11434/v1`, i.e. Ollama).

Everything documented below is browser-side JavaScript in
[`web_app/frontend/src/llm/`](../../web_app/frontend/src/llm/). The module is deliberately self-contained:
it receives run data **only as React props** and imports nothing from the rest of the app (see the
import-boundary rule in [AGENTS.md](../../AGENTS.md) and [`src/llm/README.md`](../../web_app/frontend/src/llm/README.md)).
There is no Python counterpart — none of this exists in `rmc_toolkits/`, and the Flask backend is never
involved in an assistant request.

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
| 1 | AI Assistant page mounted (automatic, once per base URL), or **Test** pressed | `GET {baseUrl}/models` → `listModels()` | **None.** Only the API key, if one is set. |
| 2 | You press send in the chat | `POST {baseUrl}/chat/completions` (`stream: true`) → `postChat()` via `streamChat()` | The summarized context JSON (target ≤ ~4,500 chars, best-effort — Step 12) + the system prompt + the last 8 chat messages + your question. |
| 3 | Watchdog enabled *and* a model configured | `POST {baseUrl}/chat/completions` (`stream: false`) → `completeChat()` | Five numbers only: `{n_steps, first, last, min, recent_window_delta}` plus two status words. |

**Row 1 fires without you asking.** `AssistantPage.jsx` calls `useAssistant({… enabled: true})`, and
`useAssistant.js` runs `probe(false)` from an effect guarded only by `autoTestedRef` (one probe per
distinct `baseUrl`). Merely opening the AI Assistant tab therefore issues a `GET /models` to whatever
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
  pair), and the full displacement clouds (only per-site scalars).

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
  numbers listed in Steps 2–11.

**Local vs cloud.** `PROVIDER_PRESETS` in
[`llm/provider/presets.js`](../../web_app/frontend/src/llm/provider/presets.js) marks Ollama
(`http://localhost:11434/v1`) and LM Studio (`http://localhost:1234/v1`) as `cloud: false`, and OpenAI,
Gemini, and Anthropic as `cloud: true`. `isLocalUrl()` classifies any base URL whose hostname is in
`{localhost, 127.0.0.1, [::1], ::1, 0.0.0.0}` as local. When the active provider is cloud, the settings
drawer renders a "Sends your run data to *provider* — it leaves your device" banner and an API-key
field ([`ConnectionSettings.jsx`](../../web_app/frontend/src/llm/components/ConnectionSettings.jsx)).
Cloud use is opt-in: the user must select the preset and paste a key — but note that once a cloud preset
is active, the row-1 probe above authenticates with that key on every fresh page mount.

**No rmc-toolkits server is ever contacted.** There is no telemetry, no proxy, no relay. On GitHub Pages
(static mode) the page is served from `drthyang.github.io` but the model traffic goes browser →
`localhost` (or browser → the cloud provider you chose) directly. The API key is stored only in
`localStorage` and sent only as an `Authorization: Bearer …` header, and only when non-empty
(`authHeaders()`). For `api.anthropic.com` two extra headers are added
(`anthropic-version: 2023-06-01`, `anthropic-dangerous-direct-browser-access: true`, in
`requestHeaders()`); no other provider-specific handling exists.

**You can inspect the exact payload — but only before your first message.**
[`ChatView.jsx`](../../web_app/frontend/src/llm/components/ChatView.jsx) renders a `<details>` panel
labeled *"Context sent to the model"* containing `contextToJson(context)` — the identical string that
`contextBlock()` embeds in the request. That panel lives inside the **empty-state** branch
(`hasConversation = turns.length > 0 || streaming`), so it disappears as soon as you send the first
message and does not return until the page is reloaded (which also clears the conversation). It also
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
([`llm/useAssistant.js`](../../web_app/frontend/src/llm/useAssistant.js)). `AssistantPage` is mounted only
after the user first visits the AI Assistant tab (`visitedPages.assistant` in
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
[`browserData.js`](../../web_app/frontend/src/browserData.js). **Static mode only** — the Flask
`/api/structure` response carries no `moves` key.

**Parse gate.** `readMovesMetadata()` scans **only the header**: `text.slice(0, …)` up to the index of
the literal `Atoms:` marker, or the **first 4,000 characters** when that marker is absent or at index 0.
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

Ratios are emitted only when the denominator is $>0$; the block is dropped entirely if empty.

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
  `parseRunSettings()` itself caps the dataset list at 8 blocks.

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
- `ladder` — a list of rungs `{sg, holds_A: [from, to], n_ops}`: the space group the average sites
  satisfy over each Cartesian tolerance interval in Å. **Emitted only when the ladder has more than one
  rung** (`symmetry.ladder.length > 1`); a structure whose ladder collapses to a single rung gets no
  `ladder` key at all, and nothing records that it was suppressed rather than not computed.
- `sites` — one entry per Wyckoff orbit:
  `{element, multiplicity, wyckoff?, site_sym?, frac?, mean_disp_A?, max_disp_A?}`, where
  - `multiplicity` = `orbit.size`, the number of basis sites in the orbit;
  - `wyckoff` = multiplicity and letter concatenated (`` `${orbit.size}${orbit.wyckoff}` ``, e.g.
    `"16e"`), present only when a letter was identified;
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
reference site $s$ with $N_c$ copies (atoms), axis $i$, and within-unit-cell fraction
$w_i = \operatorname{frac}(x_i N_i)$, the circular resultant and circular standard deviation are

$$\bar R_i = \frac{1}{N_c}\left\lVert \left(\textstyle\sum \cos 2\pi w_i,\ \sum \sin 2\pi w_i\right)\right\rVert,
\qquad \sigma_i^{\mathrm{frac}} = \frac{\sqrt{-2\ln \bar R_i}}{2\pi},$$

and the site's scalar displacement is the quadrature sum converted to Å with the conventional-cell edge
lengths:

$$u_s = \sqrt{\sum_{i=1}^{3}\left(\sigma_i^{\mathrm{frac}}\, a_i\right)^2}, \qquad a_i = \lVert\mathbf L_i\rVert / N_i .$$

Axes with $\bar R_i \ge 1$ (a single atom, or zero spread) contribute nothing; $\bar R_i$ is floored at $10^{-6}$.
This is a **single-snapshot** spread: it mixes static disorder with thermal motion, and it uses edge
lengths rather than the full metric tensor, so for a strongly non-orthogonal cell the Å conversion is an
approximation.

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
| `anisotropy` | $\sqrt{\lambda_1/\lambda_3}$ = rms₁/rms₃ | — | **no** — assigned with no finite check |
| `non_gaussianity` | mean over the three PCs of the excess kurtosis $m_4/m_2^2 - 3$ | — | only when `Number.isFinite(site.nonGaussianity)` |
| `degenerate` | present only when `true` (see below) | — | conditional |

**Silent field dropping.** `U_iso_A2`, `rms_axes_A` and `anisotropy` are assigned unconditionally; when
the source value is missing or non-finite, `roundSig` returns `undefined` and `JSON.stringify` deletes
the key. There is no `null` and no `*_omitted` marker, so an incomplete PCA row can arrive at the model
as `{ref, element}` alone with nothing explaining why. `non_gaussianity` is the only field with an
explicit finite guard.

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

**Ordering:** descending `non_gaussianity`, ties broken by descending `U_iso_A2` (missing values sort as
$-\infty$), then truncated to `MAX_SITES = 12` with `sites_omitted`. A fixed `note` string travels with
the block defining every symbol, because small models otherwise misread the kurtosis sign convention.

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
5. The entry is kept only if at least one of `avg_structure_d_A` / `gr_peaks_A` is present.

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
CSV has at least 3 columns, using columns 1/2/3 as $(x, y^{\mathrm{obs}}, y^{\mathrm{calc}})$. `pdf_partials`,
`exafs_q`, `exafs_r` (and `r_value`, which is excluded anyway) never have one, so in a real run most
`datasets` entries are `{kind}` or `{kind, title}` with no residual — their absence is not an error.

The metric is computed by `rwp(x, obs, fit)` in `browserData.js`; the Python equivalent is `rwp()` in
[`rmc_toolkits/parsers.py`](../../rmc_toolkits/parsers.py). **Both compute the same unweighted quantity, and
they agree exactly, including the zero-denominator branch:**

$$R = \sqrt{\frac{\sum_i (y^{\mathrm{calc}}_i - y^{\mathrm{obs}}_i)^2}{\sum_i (y^{\mathrm{obs}}_i)^2}}$$

with no weights $w_i$ — despite the name `rwp` and despite the system prompt describing it as "a weighted
profile residual". Treat it as a relative goodness indicator between datasets of the same run, not as a
Rietveld $R_{wp}$. **Zero-denominator caveat:** when $\sum_i (y^{\mathrm{obs}}_i)^2 = 0$ both
implementations return exactly `0` (not `NaN`), so an `rwp: 0` entry means *"the observed column was all
zeros"*, not *"a perfect fit"*.

**Code:** `runContext.js` → `datasetContext()`; `browserData.js` → `plotDataFromText()`, `rwp()`.

---

### Step 11 — Convergence series (`context.convergence`)

**Input:** `rValueFile.plotData.series[0].y`.

**What the values are.** RMCProfile's `<stem>-NN.log` files are read by `readChi()` in
`browserData.js`: skip the first 2 lines, take the **last whitespace-separated number** of every line
with ≥2 fields. The dashboard then stores $y_k = \ln\!\left(\max(\chi_k, 10^{-12})\right)$ — the natural
log, with a hard floor of $10^{-12}$ on the raw $\chi$ before the log, so a $\chi$ of exactly 0 maps to
$\ln(10^{-12}) \approx -27.63$ rather than $-\infty$.

**How the logs are combined.** `combineRValueFiles()` in
[`Dashboard.jsx`](../../web_app/frontend/src/components/Dashboard.jsx) **concatenates** the log-transformed
series of every visible R-value log and re-indexes $x$ as $0,1,2,\dots$; `final_chi_r` is taken from the
*last* parsed file. The order is `comparePlotFiles`: first by `plotOrder` index of `plotKind`, then — for
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

with $m$ the tail least-squares slope defined in Step 14. On output, `n_steps` is the **exact integer
point count and is never rounded**; `first`, `last`, `min`, `max` are 3 s.f.; `recent_slope_per_step`
is 2 s.f.

**Downsampling (`downsampleSeries`).** Uniform-stride resampling to at most `HISTORY_POINTS = 48` points
that always retains the endpoints:

$$\mathrm{step} = \frac{N-1}{P-1}, \qquad \tilde y_i = y_{\,\mathrm{round}(i\cdot\mathrm{step})},\quad i = 0,\dots,P-1$$

for $P = 48$ (or fewer, see Step 12). Series with $N \le P$ pass through at full length. Every retained
value is rounded to 3 s.f. This is *decimation, not averaging* — no local minimum between kept indices
survives, which is why the min/max/slope statistics are computed on the full series first.
(`HISTORY_POINTS` is also the default argument of the exported `downsampleSeries`.)

**Also emitted:** `quantity` — the literal string
`"ln of chi^2 goodness metric (natural log; lower is better)"` — and `final_chi_squared` =
`plotData.metrics.final_chi_r`, the **un-logged** last value, so
$\mathrm{last} = \ln(\text{final\_chi\_squared})$ up to rounding.

> **Naming discrepancies (trust the code).** Three labels attached to this one array disagree with each
> other and with the parser:
>
> 1. **$\chi$ vs $\chi^2$.** The parser and the Python plotter label the quantity $\chi$
>    (`yLabel: 'log(χ)'` in `browserData.js`; `ax.set_ylabel(r"log($\chi$)")` in `plots.py`; the metric
>    key is `final_chi_r`), while the assistant context and the watchdog comments call it $\chi^2$. The
>    numbers are identical either way — only the label the model is told differs. It matters for
>    interpreting relative changes: a 0.02 shift in the stored $\ln$ value is a 2 % change in *the
>    quantity as parsed*, which would be 4 % if that quantity is really $\chi$ and you wanted $\chi^2$.
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

**Input:** the assembled context object and a budget, default `CONTEXT_CHAR_BUDGET = 4500` **characters**
(not tokens). Serialization is `JSON.stringify(context, null, 1)` — pretty-printed with a 1-space indent,
so newlines and indentation count against the budget.

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
never whole pairs). A pathological run can therefore exceed 4,500 characters after all eight steps. The
tests assert the budget holds for a 30-dataset, 50,000-step fixture and that a forced budget of 10
characters triggers every trim step in order, but no test asserts the budget is met unconditionally.

**Code:** `runContext.js` → `contextToJson()`.

---

### Step 13 — Message assembly (what the request body actually contains)

**Chat** (`buildChatMessages()` in
[`llm/prompts/templates.js`](../../web_app/frontend/src/llm/prompts/templates.js)), in order:

1. `system` — the fixed `SYSTEM_PROMPT` from
   [`llm/prompts/system.js`](../../web_app/frontend/src/llm/prompts/system.js). No run data, but it does
   carry hard numeric interpretation rules (below).
2. `user` — `"Context for this RMC modeling run:"` followed by ```` ```json … ``` ```` containing
   `contextToJson(context)`. The fence exists so the model can distinguish data from instructions.
3. `assistant` — a fixed one-line acknowledgement.
4. the last `CHAT_HISTORY_TURNS = 8` prior messages (a flat `slice(-8)` over `{role, content}` — 8
   *messages*, roughly 4 exchanges).
5. `user` — the question you typed.

**What "history" contains.** The `turns` array lives in `ChatView`'s local `useState` — never persisted,
cleared on unmount or reload. The user message is appended **before** the request is sent; the assistant
message is appended **only if `reply.text` is non-empty**. A stopped stream, or a reasoning-only reply
that produced no answer text, is therefore dropped, leaving the user's question in the re-sent history
with no answer beside it. `buildChatMessages()` maps each turn to `{role, content}` only, so
`reasoning` / `reasoningMs` (chain-of-thought) are **never re-sent**.

The full context is **re-sent on every message**, which is why old turns can be dropped cheaply.
Request body: `{model, messages, temperature, stream: true}` with `temperature` from settings
(default **0.2**).

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
`JSON.stringify(stats)` (Step 14), the heuristic label, the previous label if any, and a demand for
exactly one line of the form `` `STATUS: improving|converged|stalled|diverging — <one short sentence
citing a number>` ``. **This payload is not fenced** — unlike the chat context, the stats JSON is
interpolated into the middle of an English sentence, so the data/instruction separation the chat path
argues for is not applied on the second request path. (The payload is five numbers, so the risk is
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

**Watchdog payload (`watchdogStats`).** Exactly five fields:
`{n_steps, first, last, min, recent_window_delta}`. `n_steps` is the exact integer point count (never
rounded); `first`, `last`, `min` are 3 s.f.; `recent_window_delta` is 2 s.f. The full history is *never*
sent to the watchdog. Three implementation details worth knowing:

- It **recomputes** `windowDelta(values)` — a second full `recentSlope` pass over the tail — instead of
  reusing `stats.recentSlopePerStep` from the `seriesStats()` call it already made.
- It calls `Number(x.toPrecision(n))` **directly, not `roundSig`**, so unlike everywhere else in the
  module a non-numeric `first` / `last` / `min` throws a `TypeError` rather than passing through.
- `max` and `recentSlopePerStep` from `seriesStats()` are deliberately dropped.

**Re-ask gate (`significantChange`).** Given the stats of the last LLM call and the current stats:

$$\mathrm{significant} \iff
(\,n_{\mathrm{steps}}^{\mathrm{new}} - n_{\mathrm{steps}}^{\mathrm{old}} \ge 200\,)
\ \vee\
\big|\,y^{\mathrm{new}}_{\mathrm{last}} - y^{\mathrm{old}}_{\mathrm{last}}\big| \ge \ln(1 + 0.02) \approx 0.0198$$

i.e. 200 new log lines, or a ≥2 % move in $\chi$. It returns `true` when there are no previous stats and
`false` when there are no current stats. **Non-finite fallback:** if either `last` is not a finite number
the $\ln(1.02)$ test is skipped entirely and the function returns the strict inequality
`prevStats.last !== nextStats.last` — which for two `NaN`s is `true`, forcing a model call on every poll
that clears the interval gate. Unlike `WINDOW_DELTA_EPSILON` and `MIN_TOTAL_DROP` (module constants), the
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
  `[enabled, nSteps, lastValue, baseUrl, model, apiKey, watchdogIntervalMin]` — not object identity,
  because Live Data produces a fresh object on every poll. **Blind spot:** a history whose *length and
  final value* are unchanged but whose interior changed does not re-classify. During Live Data this is a
  real (if brief) staleness window.
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
  blocked" (indistinguishable from the page) and quotes `window.location.origin` so the
  `OLLAMA_ORIGINS` value is copy-ready; 401/403/404/429 get specific messages. **One exception:** an
  `AbortError` (the caller's own abort signal) is **re-thrown**, not translated. The function's own
  comment says it "never throws"; the code disagrees, and the code is what runs.
- `streamChat()` parses the SSE body manually: split on `\n`, keep only lines starting with `data:`,
  buffer the incomplete tail across reads (chunks can split mid-line), stop at `data: [DONE]`, ignore
  unparseable lines. It yields `{content}` or `{reasoning}` deltas — reasoning models report
  chain-of-thought in `delta.reasoning` (Ollama) or `delta.reasoning_content` (DeepSeek and others),
  surfaced separately in the collapsible "Thinking" panel. The reader is always cancelled in a `finally`.
- `completeChat()` is the non-streaming variant used by the watchdog.
- Rendering of replies is Markdown via `react-markdown` with `rehypeRaw → rehypeSanitize → rehypeKatex`;
  `<img>` is stripped from the sanitize schema so model output cannot trigger an external load
  ([`ChatView.jsx`](../../web_app/frontend/src/llm/components/ChatView.jsx)).

---

### Runtime-mode differences (static vs Flask) — read this before quoting the assistant

The context is far richer in **static mode** (the browser parses the run itself) than in **Flask mode**,
because the Flask endpoints return figures and metadata rather than series and basis data. What survives
in each mode:

| Context block | Static mode (browser parsing) | Flask mode (backend) |
|---|---|---|
| `run`, `live_mode` | yes | yes |
| `structure` (cell, angles, composition) | yes | yes (`/api/structure`) |
| `configuration_optimization` | yes | **no** — `/api/structure` returns no move counters |
| `run_settings` | yes | **no** — the `.dat` is read only from a locally picked folder |
| `symmetry` (+ per-site displacements, `frac`) | yes | **no** — `/api/structure` returns no `basis`, so `describeSymmetry()` returns `null` |
| `pair_correlations` | yes (needs parsed partials) | **no** — plot files carry no `plotData.series` |
| `datasets` | `kind` + `title` + `rwp` (never `stog`) | `kind` only (never `stog`) |
| `convergence` + watchdog badge | yes | **no** — no parsed R-value series, so the badge stays hidden |
| `pca_displacements` | after opening the PCA Ellipsoid page | after opening the PCA Ellipsoid page (`/api/pca/sites`) |

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
| `CONTEXT_CHAR_BUDGET` | 4500 | characters | `runContext.js` | Target size of the serialized context |
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
| `CHAT_HISTORY_TURNS` | 8 | messages | `templates.js` | Prior chat messages re-sent (in-memory only) |
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
- **Cloud providers receive the context.** Selecting OpenAI, Gemini, or Anthropic sends the summarized run
  context (including run/dataset names, composition, and the `frac` triples above) to that provider under
  your API key, and merely opening the page issues an authenticated `GET /models` before you type
  anything. This is opt-in and warned about in the UI, but it is a real disclosure. Local providers keep
  everything on the machine.
- **The context inspector is only available before your first message.** It lives in the chat's empty
  state, so once a conversation starts you cannot re-read what was sent without reloading (which clears
  the conversation). It also re-serializes the *current* context, which under Live Data may differ from
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
- **`mean_disp_A` / `max_disp_A` mix static disorder with thermal motion** and are derived from a
  single-snapshot circular standard deviation, converted to Å with cell *edge lengths* — an approximation
  for non-orthogonal cells.
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
- **`rwp` is unweighted** (Step 10) despite its name and the system prompt's description, exists for only
  five dataset kinds, and returns exactly `0` for an all-zero observed column — which reads as a perfect
  fit but means the opposite.
- **The quantity is labeled $\chi^2$ in the context but $\chi$ in the parser and the Python plots**
  (Step 11). Same numbers, contradictory labels.
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
  - `rwp` — `browserData.js` and `rmc_toolkits/parsers.py` implement the identical unweighted formula
    *and* the identical zero-denominator fallback. **Exact match.**
  - `detect_plot_kind` — the `pdf_partials` (`*PDFpartials*.csv`), `npdf`, `xpdf`, `xray_sq`,
    `neutron_sq`, `bragg`, `exafs_q`, `exafs_r` and `r_value` (`*-NN.log`) branches match exactly, which
    covers everything this section depends on. The **`stog` branch differs**: the JS returns `'stog'` for
    any name matching `/\.(gr|sq|fq)$/i`, while the Python returns `'stog'` only for the three literal
    default names `{scale_ft.gr, scale_ft.sq, scale_ft_rmc.fq}`. (This has no effect on the assistant,
    which never sees `stog` files at all.)
  - `DEGENERATE_RATIO` in the PCA engines (Step 6) also matches exactly between `pcaKde.js` and
    `pca_kde.py`.
