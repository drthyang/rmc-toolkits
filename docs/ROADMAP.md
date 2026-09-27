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

### Preprocessing — Auto StoG (ACTIVE; shipped 2026-07-17)

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
today (dashboard R-values/convergence, Atomic Density, PCA Ellipsoid, AI assistant, Live Data),
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
  - Making Auto (or a coarser ν) the default resolution instead of ν = 10. Measured during the
    audit, 1 of 52 real sites exceeded 2σ at ν = 10, against 8–11 on Auto.
  - Removing or renaming the legacy `significance`/`peakZScore` fields.
  - A variance-matched map test, and an exact antipodal null at coarse ν.
  - A reference-position option, since coherent off-centring is invisible by construction.
  - A physical sign rule for the PCA frame.
- **PCA ellipsoid:**
  - Identical subsamples in both engines above 20 000 copies.
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
