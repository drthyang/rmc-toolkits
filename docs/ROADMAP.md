# Development Roadmap

## Product Vision

Build a **browser-first** analysis app for RMCProfile modeling workflows: atomistic configuration
optimization under experimental constraints. Anyone can open the hosted dashboard, import a run
directory, inspect detected outputs, compare fits, explore atomic structures, generate KDE slices,
monitor R-values, and export publication-ready figures — with no install and no data ever leaving
their device. An optional Flask backend adds server-side file browsing and reference-grade
computation for local or self-hosted use (plus API routes, such as `.rmc6f` → `Frac_coord`
conversion, that no page calls).

The next frontier is lowering the technical barrier to *running* RMCProfile itself: guided run setup
with validation and ready-to-use input files, so researchers spend less time hand-editing command
files and command lines. STOG-style scaling/Fourier-transform work is preprocessing, so it stays out
of the visible app scope for now.

## Phase 1: Foundation

Goal: make the project reliable enough to build on. **Status: largely complete.**

- Package reusable analysis logic under `rmc_toolkits/`. **Done for the active parser, plot, and KDE paths.**
- Keep CLI scripts as thin wrappers over package functions. *(Pending — legacy `src/` scripts are still standalone; see backlog.)*
- Add and maintain tests for: RMC CSV parsing, log chi parsing, Rwp calculation, plot-kind detection,
  `.rmc6f` lattice/atom parsing, and `.rmc6f` → `Frac_coord*.txt` conversion. **Done.**
- Add a small fixture-based smoke test using `data/`. **Done** (sample-backed tests skip cleanly when the dataset is absent, so CI stays green).
- Add formatting and linting commands for Python and frontend code. **Done** (ESLint + CI run lint/build/tests on every push).
- Document app startup, data-root behavior, and expected file patterns. **Done** (README, QuickStart, AGENTS).

## Phase 2: Project Workspace

Goal: move from file-by-file plotting to project-level analysis. **Status: mostly complete.**

- Project scanner that detects RMC CSV outputs; PDF, S(Q), Bragg, partials, logs, and RMCProfile
  EXAFS dataset Q/R outputs; and `.rmc6f`/`Frac*.txt` structure files. **Done.**
- Frontend workspace layout: header data-path controls, a dashboard page for run plots and model
  summary, and a KDE / 3D page for structure exploration. **Done.**
- Optional Live Data monitoring so charts refresh when files change. **Done** — both the Flask
  dashboard (server-side watch) and the hosted browser app (File System Access API, Chromium).
- Loaded-file controls for hiding/showing individual dashboard charts. **Done.**
- Project summary JSON with file roles, available plots, metrics, lattice metadata, element list, and warnings. *(Partial — metadata is surfaced per file; a single `/api/project/scan` summary is still planned.)*
- Run comparison across multiple directories. *(Planned.)*

## Phase 3: Interactive Analysis

Goal: replace static inspection with interactive scientific workflows. **Status: mostly complete.**

- Browser-native SVG plots with hover readouts, legend toggles, and drag-to-zoom. **Done.**
- Plot controls — drag-zoom range and log scale (KDE) **Done**; per-figure export resolution **Done**
  (PNG / 3× / SVG). Residual-panel display is *planned*.
- KDE slice viewer: element selector, z-position slider (plus **drag-the-slab-band** to set the
  slice), slab thickness, bandwidth, colormap selector, and contour toggle. **Done.**
- Server-side SciPy KDE path **Done**; export/reproducibility controls (parameter capture) *planned*.
- Hosted GitHub Pages path with browser-side parsing and Web Worker KDE (WebGPU-accelerated with an
  automatic CPU fallback). **Done.** Live Data also works here in Chromium via the File System Access
  API — it is no longer Flask-only.
- Three.js structure viewer: per-element coloring, unit-cell display, and figure/screenshot export
  (PNG, native or 3×). **Done.** Standalone element visibility toggles and camera presets are *planned*.
- PCA Ellipsoid page: per-site PCA of RMC displacement clouds → thermal ellipsoid + separable 3D
  Gaussian-KDE isosurface with wall projections/contours, non-Gaussianity readout, and a clickable
  unit-cell structure picker. **Done** (method after Maksim Eremenko's PCA_KDE). Isosurface/CSV
  export of the volume is *planned*.

## Phase 4: Background Jobs

Goal: make expensive analysis responsive and reproducible. **Status: planned.**

- Add a job model for KDE, structure transforms, batch plots, and report generation.
- Start with SQLite-backed local jobs.
- Track job status, input paths, parameters, output artifacts, runtime, and error messages.
- Add frontend job status indicators and retry controls.
- Cache expensive computed arrays and generated plots.

## Phase 5: Reporting And Export

Goal: make the app useful at the end of a research session. **Status: in progress.**

- Export individual plots as PNG and SVG, and bundle a whole dashboard into one `.zip`. **Done.**
  CSV (raw series) export is *planned*.
- Export project summary as JSON. *(Partial — the AI Copilot's context builder assembles a
  compact run-summary JSON, viewable in its "context sent to the model" inspector.)*
- Generate a reproducible report containing input directory, detected files, software version, plot
  parameters, Rwp metrics, lattice metadata, and selected figures. *(Planned. The chat-only
  assistant's Markdown run report was removed with it.)*
- Add figure presets for manuscript, talk, and notebook usage. *(Planned.)*

### Experimental: AI Copilot

An experimental track (module `web_app/frontend/src/llm/`, see its README and
[COPILOT_PLAN.md](COPILOT_PLAN.md)) connecting the Workbench to a model straight from the browser:
a local one the user runs (Ollama, LM Studio; run data stays on the machine) or a cloud one with
the user's key.

- The AI Copilot docked beside the pages: it runs the Workbench analyses itself, checks each result
  and its own answer, asks when a choice changes the result (and the Workbench asks itself where the
  run data show one), follows written skills, points at the cards it explains, and keeps a
  convergence watchdog badge (heuristics-first, LLM-narrated). **Done (experimental).**
- Possible next steps: feed KDE-slice findings into the context, run comparison Q&A, and a guided
  "why is my fit bad?" diagnostic flow.

### AI Copilot speed (planned; deferred)

The maintainer found local models "a bit slow" and asked whether the laptop or the agent's workflow
was at fault. On 2026-10-11 they put speed **after** correctness: the tools, the workflow's
scientific accuracy and its robustness come first (the `ask_user` question card and the context
guard shipped then). The diagnosis below comes from Ollama 0.40.2's own server log on the
maintainer's Apple M5 Pro (64 GB, 51.8 GiB usable by the GPU), 2026-10-09 to 10-11, while another
job (Athanor) shared the server with qwen3:32b and gemma4:26b.

**Where the time went:**

- **Waiting for another job, not computing.** In the view-tool checks (2026-10-11, 00:31–01:14),
  each Copilot request on gemma4:26b took 8–22 s of server time: 8,000–10,200 prompt tokens
  processed in 7–13 s, and 20–450 tokens generated in 0.3–9 s. Yet each request took 2.5–6.7 min
  end to end, because it queued behind the other job's requests (one slot) and model swaps:
  gemma4:26b (17 GB, plus 5.4 GB of KV cache at the 262,144-token default context) and qwen3:32b
  (31 GB at 40,960) do not fit together in 51.8 GiB, so the log shows 137 model loads in 39 h.
  The README's 223–847 s per question for that check were mostly this.
- **Reloads between questions.** With the default 5-minute keep-alive, questions more than five
  minutes apart reload the model (17 gemma4 loads on 2026-10-10 between 10:00 and 19:35, 2–9 s
  each) and lose the prompt cache.
- **The prompt cache misses on every new question (gemma4).** gemma4 uses sliding-window attention
  (1,024-token window). llama-server can reuse such a cache only from a saved checkpoint at or
  below the prefix the new prompt shares, and its checkpoints sit near the ends of earlier
  prompts. So each new question re-processes the whole 7,000–9,000-token prompt ("forcing full
  prompt re-processing due to lack of cache data", 2026-10-10 13:08, 6 s at 1,400 tokens/s).
  Within a question, each step extends the last prompt and is reused (8,428 of 8,645 tokens).
- **Dense models are an order of magnitude slower.** Prefill: gemma4:26b (mixture of experts)
  850–1,400 tokens/s against qwen3:32b's 130–180 tokens/s (~60 s for one 8,400-token prompt).
  Generation: 50–63 against ~10 tokens/s. qwen/qwen3.8-27b (LM Studio, a dense reasoning model at
  its default reasoning effort) took 1–5 min per question in the 2026-10-10 spot check, and one
  reasoning step 10.5 min. The Copilot sends no reasoning control.
- **The laptop is not the bottleneck** for gemma4:26b: before the other job started (the
  2026-10-10 spot check, no other model loaded) a question took 8–95 s.

**Plan** (each item measured before and after on an idle server: the six spot-check questions of
the llm README, gemma4:26b, `main` against the branch):

1. **Measure first:** a timing line on each step card and answer (load, prompt processing,
   thinking, writing, tool time; details on expand), from Ollama's native durations
   (`load_duration`, `prompt_eval_*`, `eval_*`), OpenAI-style `usage` with
   `stream_options.include_usage`, and client timestamps otherwise.
2. **Keep the model loaded** while the Copilot is open (`keep_alive`), and optionally warm it when
   the Copilot opens: load it and process the static part of the prompt once, which also leaves a
   sliding-window checkpoint at the end of that prefix so every new question reuses it.
3. **Prompt caching:** keep the static prefix (rules, tool schemas, then the run context)
   byte-identical across questions; the tool schemas are the largest part (~14,800 characters,
   measured on the Demo before `read_skill` and the skills index were added).
   Sending only a page's or question's tool group would change that prefix, so do it only if the
   cache keeps working.
4. **Reasoning control, "Fast / Thorough":** little thinking while choosing tools and normal
   thinking for the answer, through Ollama's `think`, OpenAI-style `reasoning_effort` and LM
   Studio's equivalent, detected per server, with a silent fallback when a server rejects the
   field. Check that toggling it per step does not change the prompt's start (some templates put
   the thinking switch in the system turn), or the cache breaks.
5. **Fewer round trips:** let the model request independent analyses in one step and run them
   concurrently; keep answering from the context when it holds the answer.
6. **Recommend gemma4:26b** for local use in the settings, noting that dense reasoning models are
   much slower; and say in the setup guide that a second model on the same server queues the
   Copilot's requests.
7. **Context size and memory:** the context guard (shipped) asks for a context that fits the
   prompt. A smaller context than the 262,144 default would also shrink gemma4's KV cache from
   5.4 GB to under 1 GB, but a context that differs from another client's reloads the model on
   every alternation; decide with the measurements.

## Phase 6: Lab-Ready App

Goal: make the tool safe and pleasant for broader use. **Status: planned (partially started).**

- Add authentication if served beyond localhost. *(Planned.)*
- Keep data-root restrictions enabled by default. **Done.**
- Add project persistence and recent projects. *(Planned.)*
- Package as a one-command local app. *(Partial — a `Dockerfile` builds and serves the full stack.)*
- Add robust error messages for malformed files. *(Partial.)*
- Add documentation with example workflows and screenshots. **Done** (README + QuickStart).

## Phase 7: Guided RMCProfile Setup

Goal: **minimize the technical barrier to running RMCProfile.** Today the app visualizes finished
outputs; the next step is form-driven run setup, so going from measured/reduced data to a configured
RMC modeling run does not require hand-editing input files or the command line. **Status: next major
focus.**

### Guided RMCProfile run setup

- Form-driven generation of RMCProfile input files (data sets, fit ranges, constraints,
  swap/move/translate moves, supercell) from curated templates, with validation and sensible
  defaults instead of hand-edited text.
- A pre-flight check before a run: missing data files, unit/format mismatches, and density/lattice
  sanity, with plain-language diagnostics.
- Presets for common dataset combinations (neutron / x-ray total scattering, combined runs, and
  runs that include RMCProfile EXAFS datasets).

### Preprocessing — Auto StoG (ACTIVE; shipped)

The former "deferred preprocessing" precondition — a dedicated preprocessing module with a clear
user path — is met: `rmc_toolkits.scaling`/`transforms`/`scattering` (the auto-scaling engine,
validated against three complete classic-Fortran stog runs), the `rmc-autoscale` CLI and the
`/api/scaling/preview|run` endpoints replace the classic stog "try again" loop with an automatic
physics-anchored fit (high-Q level sweep + low-r density limit, with the Faber-Ziman Q→0 limit as
an independent amplitude criterion/cross-check). The **Auto StoG** page runs the Web-Worker port
(`workers/autoScale.js`) in both runtimes, but the tab is hidden in the shipped build
(`SHOW_AUTO_STOG = false` in `App.jsx`); the engine, CLI and API are supported. Full plan and
validation record: [STOG_SCALING_PLAN.md](STOG_SCALING_PLAN.md).

### Lowering the barrier end-to-end

- A guided path **reduced data → configured RMC run**, with plain-language explanations and links
  back to the RMCProfile documentation.
- Stay browser-first wherever possible; offload only genuinely native steps (executing the
  RMCProfile binary) to the optional local backend or to a downloadable, ready-to-run input bundle.

## Phase 8: Remote / HPC Run Monitoring

Goal: **monitor RMCProfile runs executing on an HPC cluster** with the same feature set available
today (dashboard R-values/convergence, Atomic Density, PCA Ellipsoid, AI Copilot, Live Data),
without hand-copying the run down. **Status: planned.** Full design in
[`docs/HPC_MONITORING_PLAN.md`](HPC_MONITORING_PLAN.md) (written for both this project and the
RMCProfile team).

Principles: **security first, performance second**; **listen-only / read-only** (the app monitors and
never writes to or controls the run); **no third parties** (data flows only between the HPC and the
user's own machine); auth via the user's **existing SSH trust** (public key / `ssh-agent`, bastion
and MFA aware) — the app never stores credentials.

- **Phase 8a — MVP (our side only, no RMCProfile changes):** the local Flask backend pulls the run's
  output files read-only over SSH (`rsync`/`sftp`, incremental) into a local cache; every existing
  page reads that cache like a local folder. Requires the local backend (browser-only static mode
  can't open SSH).
- **Phase 8b — status file:** consume an RMCProfile-emitted, atomically-updated JSON status/heartbeat
  file (step, χ²/Rwp, phase, ETA) for cheaper, more robust convergence monitoring.
- **Phase 8c — optional:** opt-in low-latency push over an SSH reverse tunnel; SLURM/PBS scheduler
  status (`squeue`/`sacct`, read-only). See the plan for the RMCProfile-team asks and open questions.

## Architecture Target

- `rmc_toolkits/`: pure Python package for parsing, analysis, plotting, structure transforms, and
  (planned) RMC input-file generation.
- `web_app/backend/`: API server, project scanner, optional run-setup/input generation, jobs, artifact storage.
- `web_app/frontend/`: React app with project workspace, interactive plots, KDE, structure viewer,
  figure export, and (planned) guided run-setup pages.
- `.github/workflows/pages.yml`: builds the static GitHub Pages dashboard from `web_app/frontend`.
- `data/`: small example fixtures.
- `docs/`: development changelog, roadmap, and architecture notes (agent guide at repo-root `AGENTS.md`).

## Suggested Immediate Backlog

1. Add trimmed fixtures of the GaNb₄Se₈ and `stog_tests` runs so their real-data tests run in CI
   (since 0.6.0 the committed demo run backs the parser, KDE-golden and plot-parity tests there).
2. Refactor `src/RMC_plot.py` into a CLI wrapper around `rmc_toolkits.plots`.
3. Refactor `src/RMC_3D.py` to avoid Mayavi import and execution at import time.
4. Add `/api/project/scan` for directory-level summaries.
5. Add project-level warnings for missing expected files and malformed outputs.
6. Export controls for plots — PNG/SVG and dashboard `.zip` **Done**; KDE/3D PNG (native/3×) **Done**; raw-series CSV export remaining.
7. Move heavier static-mode parsing/KDE work toward transferable typed arrays and profile large `.rmc6f` files.
8. Resolve GitHub Pages/static dashboard loading on iPhone Safari by adding a unified run-source
   abstraction, explicit import diagnostics, and a verified mobile folder-access strategy.
9. Add recent-project persistence for local desktop use.
10. Draft RMCProfile input-file templates and a form-to-input generator with pre-flight validation.
11. Prototype Phase 8a remote monitoring: a read-only SSH pull of an HPC run directory into a local
    cache, surfaced through the existing run-source abstraction (see `docs/HPC_MONITORING_PLAN.md`).

## Candidates after 0.6.0

Maintainer decisions the 0.6.0 audit considered and did not take, grouped by engine (moved here
from the 0.6.0 entry of [CHANGELOG.md](CHANGELOG.md)). Until one is taken, the behaviour
documented in [ALGORITHMS.md](ALGORITHMS.md) stands.

- **Structure KDE:**
  - A physical kernel, isotropic in the real plane with a width in Å, which would remove the
    slab-layout kernel artefact. That artefact is documented and flagged on the map.
  - Normalising oblique slices by the section integral.
  - An `[uvw]` input mode, and `(100)`-style labels for the a/b/c presets.
  - A decline instead of a warning for `unresolved` maps.
- **χ² history:**
  - Plotting every χ² column and the total, and having the watchdog classify on the total. Today
    both the chart and the badge follow the last column, and the badge names it.
  - A dedicated "blown-up" watchdog status.
  - A weighted Rwp.
- **Symmetry finder:**
  - A Python port or an spglib cross-check in CI; the finder is browser-only.
  - Moving it into a Web Worker so the 2000-site and 384-operation caps can rise.
  - An origin-shift search, so Wyckoff letters are unique across origins (Ga 4c vs 4d).
  - Cell reduction for boxes declared 1×1×1, and a reduced-cell re-description to name "≥"
    lower bounds fully.
  - Subgroup enumeration instead of the greedy maximal group.
- **Displacement Directions:**
  - ~~Making Auto (or a coarser ν) the default resolution instead of ν = 10.~~ Done (plan item
    2.5, 2026-10-19): the page opens at Auto, the tool's tiling, and reads its readouts at smoothing
    0. Measured during the audit, 1 of 52 real sites exceeded 2σ at ν = 10, against 8–11 on Auto.
  - A direction finer than its cell: on Auto a small site's tiling is coarse (ν = 1 names a cube
    axis only to within 31.7°), and the 2× display smoothing blurs such a map further.
  - Removing or renaming the legacy `significance`/`peakZScore` fields.
  - A variance-matched map test, and an exact antipodal null at coarse ν.
  - A reference-position option, since coherent off-centring is invisible by construction.
  - A physical sign rule for the PCA frame.
- **PCA ellipsoid:**
  - ~~Identical subsamples in both engines above 20 000 copies.~~ Done (0.8.0): one draw.
  - `AXIS_RESOLUTION_SIGMAS` = 2 instead of 3. At 3, PC2 is never resolved on the heavy-tailed
    5 K sites, so the κ column mostly reads "—".
  - U_cif/β_ij beside the Cartesian tensor.
  - A fixed-band shell colour scale.
  - Linking the isosurface and ellipsoid levels.
  - Per-species clouds for mixed sites.
  - Exact analytic crystal-frame marginals.
- **Auto StoG:**
  - A lower default `r_cutoff`, or lowering it automatically for short bonds.
  - Switching to the FZ amplitude automatically when the density limit is degenerate (the Mn₃Sn
    refusals).
  - Refusing an unconverged, an unreliable-FZ or an r₀-above-the-shell fit.
  - Exposing the ρ₀ physical range.
  - An "enforcement not applied" chip and explicit first-peak-window fields on the page.
- **Bond angles:**
  - Pairing A = C distinct-window requests as (w12 bond, w23 bond), so the work and not only the
    angle count is bounded by the budget.
  - A different `APP_MAX_ANGLES`.
  - A CLI `--max-angles`.
  - Showing the exact count before Compute.
  - Snapping the length-histogram bins.
  - Checking the distinct-window convention against an RMCProfile TRIPLETS run.
- **Parsers and API:**
  - Flipping `iter_rmc6f_atoms`'s full-layout-only default in 2.0.
  - Renaming the `xray_sq`/`neutron_sq` kinds.
  - The mid-write guard on the analysis pages.
  - A "configuration changed" banner with Live Data off.
  - Rejecting out-of-range grids and z instead of clamping them.
  - Caching `rmc6f_problem` per file signature (each request re-reads a 64 KiB head per
    candidate, ~0.4 ms).
