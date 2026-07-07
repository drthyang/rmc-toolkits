# Development Roadmap

## Product Vision

Build a **browser-first** analysis app for RMCProfile refinement workflows. Anyone can open the
hosted dashboard, import a run directory, inspect detected outputs, compare fits, explore atomic
structures, generate KDE slices, monitor R-values, and export publication-ready figures — with no
install and no data ever leaving their device. An optional Flask backend adds server-side file
browsing, structure conversion, and reference-grade computation for local or self-hosted use.

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
- Export project summary as JSON. *(Partial — the AI assistant's context builder assembles a
  compact run-summary JSON, viewable in its "context sent to the model" inspector.)*
- Generate a reproducible report containing input directory, detected files, software version, plot
  parameters, Rwp metrics, lattice metadata, and selected figures. *(Partial — the experimental AI
  assistant exports a Markdown run report with model summary, Rwp metrics, and convergence tables,
  plus an optional LLM-written assessment; figures and plot parameters still planned.)*
- Add figure presets for manuscript, talk, and notebook usage. *(Planned.)*

### Experimental: AI assistant (local LLM)

An experimental track (module `web_app/frontend/src/llm/`, see its README) connecting the dashboard
to a user-run local LLM (Ollama / LM Studio) straight from the browser — no server, no API keys,
data stays local:

- Run summary/assessment, chat Q&A over the loaded run, Markdown report generation, and a live
  convergence watchdog (heuristics-first, LLM-narrated). **Done (experimental).**
- Possible next steps: feed KDE/symmetry findings into the context, run comparison Q&A, and a
  guided "why is my fit bad?" diagnostic flow.

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
refinement does not require hand-editing input files or the command line. **Status: next major
focus.**

### Guided RMCProfile run setup

- Form-driven generation of RMCProfile input files (data sets, fit ranges, constraints,
  swap/move/translate moves, supercell) from curated templates, with validation and sensible
  defaults instead of hand-edited text.
- A pre-flight check before a run: missing data files, unit/format mismatches, and density/lattice
  sanity, with plain-language diagnostics.
- Presets for common dataset combinations (neutron / x-ray total scattering, combined runs, and
  runs that include RMCProfile EXAFS datasets).

### Deferred preprocessing

STOG is a preprocessing/reduction tool, not a refinement or fit. Keep STOG-specific scaling and
Fourier-transform workflows out of the current app surface until there is a dedicated preprocessing
module with a clear user path.

### Lowering the barrier end-to-end

- A guided path **reduced data → configured RMC run**, with plain-language explanations and links
  back to the RMCProfile documentation.
- Stay browser-first wherever possible; offload only genuinely native steps (executing the
  RMCProfile binary) to the optional local backend or to a downloadable, ready-to-run input bundle.

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

1. Add a committed standard example run or trimmed fixtures so sample-backed tests run in CI.
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
