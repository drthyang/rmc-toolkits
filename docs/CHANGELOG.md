# Development Log

Chronological record of notable changes, newest first. For current architecture and conventions see
[AGENTS.md](../AGENTS.md); for forward plans see [ROADMAP.md](ROADMAP.md).

## Unreleased

PCA / thermal-ellipsoid KDE computation engine.

- New engine turns per-site RMC displacement clouds into anisotropic displacement tensors (the
  thermal ellipsoids) and a smooth 3D probability density. An atom's offset from the average
  structure is `coords − cellIndices/supercell` folded over the supercell boundary, mean-subtracted
  per reference site, then mapped to Cartesian Å; the covariance of each site's cloud is its ADP,
  and a Gaussian KDE of the cloud is the density the ellipsoid only approximates. This ports Maxim
  Eremenko's PCA_KDE utilities into the toolkit's dual-mode architecture.
- **The KDE is separable, and exact.** Sampling on a grid aligned with the cloud's principal axes
  makes SciPy's `H = factor²·C` bandwidth diagonal, so the 3D Gaussian factorizes into three 1D
  kernels and the volume becomes their tensor product: `N·3·grid` exponentials contracted through
  BLAS instead of `N·grid³` kernel evaluations. Measured **36–51× faster than a naïve
  `scipy.stats.gaussian_kde` sweep of the same volume, with max abs error ~1e-13** (machine
  precision) — it is the same estimator, not an approximation. A 52-site, 52 000-atom configuration
  builds all clouds in ~80 ms and all 52 ellipsoids in ~1 ms; a single-site 48³ volume solves in
  ~7 ms server-side.
- `rmc_toolkits/pca_kde.py` is the source of truth: `load_site_displacements`, `site_ellipsoids`
  (batched, one pass), `pca_kde_volume` / `site_pca_kde` (volume + PC-plane projections + iso
  thresholds by enclosed probability mass and by raw density). Axes are sign- and
  handedness-canonicalized for reproducibility; a floored eigenvalue keeps flat/linear clouds
  non-singular and flags them `degenerate`.
- `web_app/frontend/src/workers/pcaKde.js` is a straight JS port for static mode (a 3×3 Jacobi
  eigensolver stands in for `eigh`), driven by `pcaKdeWorker.js` off the main thread — ~96 ms for a
  1000-point 48³ volume in-browser, no GPU needed. `tests/test_pca_kde.py` and
  `src/workers/__tests__/pcaKde.test.js` both assert the separable volume equals the full
  brute-force estimator to round-off (15 + 14 tests), plus site extraction, supercell-boundary
  folding, mass normalization, and known-anisotropy recovery.
- Flask exposes `/api/pca/sites` (per-site ellipsoid table) and `/api/pca/kde` (one site's or one
  element's volume), both behind a per-(path, mtime) LRU cache.
- **Thermal Ellipsoids page** (new top-level tab). A site picker over all reference sites, an
  ellipsoid summary table, a Three.js scene showing the KDE isosurface (extracted by a self-written
  marching-cubes module — `workers/marchingCubes.js`, since Three's `MarchingCubes` only builds
  metaballs) nested with the p% thermal-ellipsoid wireframe and PC-axis triad, and the three
  PC-plane KDE projections as heatmaps. An isosurface-mass slider sweeps the enclosed-probability
  threshold. Works in both runtimes: Flask endpoints in backend mode, the `pcaKdeWorker` +
  `pcaKde.js` engine off the main thread in static mode (verified on the bundled Demo run).
- **Non-Gaussianity (excess kurtosis)** is reported per site. It quantifies why a KDE isosurface can
  sit inside its harmonic ellipsoid: the covariance is inflated by fat tails while the isosurface
  tracks the peaked core. Verified on `data/5K_try1` — Nb sites read ~5–7 (strongly anharmonic,
  isosurface well inside the ellipsoid), near-isotropic Ga sites read ~1 (isosurface ≈ ellipsoid),
  and a synthetic Gaussian cloud reads ~0. Marching cubes has its own sphere/normals tests; the
  isosurface-vs-ellipsoid scaling was validated to agree with theory (KDE surface ≈ 1.06× the
  ellipsoid for a true Gaussian).

Periodic boundary conditions for the KDE slice.

- The KDE under-counted near cell boundaries: folding atoms into one unit cell drops their
  periodic neighbors, so the density decayed toward faces/edges/corners (an edge blob showed
  roughly half the interior amplitude, a corner roughly a quarter), and a slab centered on the
  z=0 face missed the atoms just below z=1 entirely.
- Both KDE paths now tile periodic images from the 26 neighbor cells within a margin that covers
  the kernel reach and the slab depth (`min(0.5, max(0.1, 2*bw, thickness))`). In
  `rmc_toolkits/kde.py`, `_augment_periodic_images` feeds `oriented_kde_slice`, and `kde_slice`
  reports `slabCount` as unique source atoms and rescales the density (SciPy's `gaussian_kde`
  divides by every fit point, images included). In `localKdeWorker.js` the same augmentation runs
  before `makeSlab`, and the image factor rides in the kernel normalizer, so the WebGPU path
  inherits the fix unchanged. `kde_slice` called directly (no `source_index`) keeps the old
  truncated behavior for generic point clouds.
- Verified on `snao.rmc6f` at a z=0 slab: slab population 3110 → 6034 atoms, corner density
  20.8 → 55.4, edge-column peaks now match interior peaks (~73 vs ~75), opposite edges agree to
  ~5% (fit-subsampling noise). New tests: in-plane wrap symmetry and depth-wrap slab selection in
  `tests/test_kde.py`, plus a mirrored vitest suite for the worker (whose `onmessage` registration
  is now guarded so tests can import the module outside a worker).

Source data file names on figures (user-experience feedback).

- Every dashboard plot card shows its source file name in monospace under the title — the two
  otherwise-identical `EXAFS Q-space` cards are now distinguishable — with the full path on hover.
- The collapsed R-value strip shows its source log name(s); a combined multi-log strip lists every
  parsed log.
- The Model information title cell now displays the source `.rmc6f` file name (previously hidden
  in a tooltip), which also covers the KDE / 3D panels since they all derive from that one file.
- File names sit outside the card `h3`, so "Save all figures" export names are unchanged.

Math rendering in AI Assistant replies.

- Assistant Markdown now typesets LaTeX: integrated `remark-math` + `rehype-katex` so
  dollar-delimited math (e.g. `$R_{wp}$`, `$\chi^2$`, `$\text{\AA}$`) in replies renders as real
  math instead of literal text — including inline math inside GFM table cells, which is how models
  tend to format result summaries. KaTeX runs *after* `rehype-sanitize` (raw → sanitize → katex) so
  model HTML stays fully sanitized while KaTeX emits its own trusted markup, and KaTeX fonts are
  bundled locally so the static GitHub Pages build needs no CDN.

Demo run and refreshed screenshots.

- Added a header **Demo** toggle that loads a bundled GaTa4Se8 250 K example run (under
  `web_app/frontend/public/demo/`) so first-time visitors see a populated dashboard; a second click
  clears it.
- Recaptured the Dashboard and KDE/3D screenshots against the demo run and added an AI Assistant demo
  GIF (a local Ollama model summarizing the run as a table with inline LaTeX math). Renamed the screenshots to
  `assets/rmc-toolkits-dashboard-demo.png` / `assets/rmc-toolkits-kde-demo.png` (new paths bust the
  stale README image cache) and removed the unused legacy screenshots from `assets/`.

Relicensed to AGPLv3.

- Switched the project license from MIT to the **GNU Affero General Public License v3.0** to keep
  the code open while requiring that modified versions — including those run as a hosted network
  service (AGPL §13) — release their source. Updated `LICENSE`, `pyproject.toml` (license field +
  trove classifier), the README badge/notice, and the in-app footers.
- Added an **About & documentation** link in the Dashboard and KDE/3D footers pointing to the repo,
  which also surfaces the source for AGPL §13 network-use compliance.

RMCProfile-first positioning.

- README, QuickStart, roadmap, and agent notes now frame the current app as an RMCProfile
  modeling-output dashboard. STOG is treated as deferred preprocessing/legacy parser support, not
  a current RMC modeling feature.
- The dashboard excludes STOG plot kinds from the visible plot list and assistant context while
  keeping the underlying parser support in place.

RMCProfile wording cleanup.

- Replaced RMCProfile "refinement" language in docs and assistant prompts with RMC modeling /
  atomistic configuration optimization wording, to avoid implying Rietveld-style parameter
  refinement.
- Renamed the assistant's move-counter context from `refinement` to
  `configuration_optimization`.

Gentler AI Assistant startup.

- The automatic on-load connection probe no longer shows a red "Connection failed" before the user
  has set anything up. The status reads a neutral **"Connect a model to start"** (grey dot), and the
  settings drawer stays clean; the full error + setup hint appear only after the user actively
  presses **Test**.

Assistant "Beta" badge and a cleaner chat speaker label.

- The header badge now reads **Beta** (was "Experimental") in **amber** — the complement of the
  app's cobalt accent, and the same palette as the cloud-provider notices — so the beta reminder
  stands out instead of blending into the blue theme.
- Chat messages are labeled **Assistant** rather than the raw model id (which stays visible in the
  header's model switcher), for a cleaner transcript.

Run history and run-control settings in the AI Assistant's context.

- **`configuration_optimization` block** from the `.rmc6f` header: moves generated/tried/accepted
  plus the derived **acceptance ratio** and **accepted moves per atom** — the standard gauge of
  whether the configuration has been sampled long enough — and accumulated running time.
- **`run_settings` block** from the RMCProfile run-control file. The correct `.dat` is the one whose
  basename matches the chosen structure stem (so `chi2.dat`, `optimization.dat`, … are never picked
  up). Extracts title/material/phase/temperature, **minimum distances labeled per element pair**
  (hard closest-approach constraints — the model is told a g(r) peak pinned there may be
  constraint-limited), max move sizes per element, time/save limits, flags, and the fitted-data
  blocks. Static mode only, like the rest of the browser-parsed context.
- System prompt teaches both blocks; 9 new tests (header counters, `.dat` parser, stem-matched
  selection, derived stats, pair labeling).

Chat rendering fixes and a window-filling chat box.

- Markdown now parses **and sanitizes** inline HTML (`rehype-raw` + `rehype-sanitize`), so `<br>`
  line breaks inside table cells render instead of collapsing onto one line; scripts, event
  handlers, and images are still stripped (no injection surface, no external loads from model
  output) — covered by tests.
- The AI Assistant chat box now **fills the window height** (message log scrolls internally, the
  composer stays pinned at the bottom) and the column is a bit wider (880 → 1000 px).

Rendered Markdown in AI Assistant chat replies.

- Assistant answers now render **Markdown** — GitHub-flavored **tables**, lists, code blocks, inline
  code, headings — via `react-markdown` + `remark-gfm`, so a "summarize this as a table" request
  shows a real table instead of raw `| --- |` text. Rendering builds a React tree with raw HTML
  disallowed, so the no-injection property is preserved (covered by a test); user messages stay
  plain. (eslint: enable `ignoreRestSiblings` for the `{ node, ...props }` component overrides.)

App identity — a minimalist pair-distribution "wave" mark.

- Replaced the placeholder "R" brand with a single shared g(r)-wave SVG (flat brand blue, `#2563eb`):
  a base-aware favicon (`public/favicon.svg`), the header brand mark, and the AI Assistant
  avatar / empty-state icon. Removed the default `vite.svg`.

ChatGPT-style "Thinking" indicator for reasoning models.

- Reasoning models (e.g. qwen3 via Ollama) stream their chain-of-thought in a separate field; the
  streaming client now yields structured content/reasoning chunks and the chat shows a collapsible,
  shimmering **"Thinking"** panel (live reasoning) that collapses to **"Thought for Ns"**
  (re-expandable), plus a bouncing-dots indicator for non-reasoning models slow to their first token.
  Reasoning is preserved on each turn.

Local-distortion evidence in the AI Assistant's run context.

- The context JSON gains a **`symmetry` block** — space group, the symmetry-vs-tolerance **ladder**
  (distortion magnitude + character), and the Wyckoff-orbit **sites ranked by rms displacement** —
  plus **`pair_correlations`**: each partial-PDF pair's first g(r) peaks next to the
  nearest-neighbour distance the average structure predicts, so the model can point at the sites
  and pairs participating in short-range correlations.
- Per-site rms displacements (`dispA`) fall out of the circular-mean accumulators
  `structureFromRmc6f` already had — derived per axis from the resultant length, scaled to Å; the
  symmetry orbits now carry `members` (basis indices) so the context can aggregate per orbit.
- New `llm/context/pairCorrelations.js` (average-structure NN distances over periodic images +
  a smoothed local-maxima peak finder); context budget raised to ~4.5k chars with an explicit trim
  order (extra peaks → ladder middle rungs → history → low-displacement sites → datasets, each
  recorded with an `*_omitted` count); system prompt teaches the model how to read the new blocks.
- Tests: 51 passing (peak finder, NN distances, displacement math on synthetic .rmc6f fixtures,
  symmetry-block assembly, trim order).

AI Assistant rework (shipped as PRs #3/#4, catching the log up): moved from a dashboard card to a
dedicated **AI Assistant page** (chat-only — Summary/Report tabs removed, everything happens in the
chat), modern chat UI with a persistent connection bar (status dot, model switcher, settings gear,
auto-connect), redesigned Connection Settings, and optional **cloud providers** (OpenAI, Gemini)
with Bearer API-key auth and explicit data-leaves-your-device warnings; local Ollama/LM Studio
remains the default. Ollama CORS + Safari setup documented in-app and in QuickStart.

Experimental AI Assistant — local LLM in the data-monitor pipeline.

- Added `web_app/frontend/src/llm/`, a self-contained experimental module connecting the dashboard
  to a **local LLM** (Ollama `/v1` or LM Studio, both OpenAI-compatible) **directly from the
  browser** — no server, no API keys; run-derived summaries go only to the model server the user
  runs. Four features: one-click **run summary/assessment** (streamed), **chat Q&A** with the run
  context injected, one-click **Markdown run report** (deterministic metrics tables + clearly
  labeled AI narrative), and a **live convergence watchdog** badge on the R-value card
  (slope heuristics are the source of truth; the LLM only writes the note; piggybacks on the
  existing Live Data poll).
- Pipeline seams are explicit as a learning artifact (see the module README): context builder with
  history downsampling + a ~3k-char budget → prompt templates → hand-rolled SSE streaming client →
  streaming UI with a "context sent to the model" inspector. Connection test translates failures
  into actionable CORS/setup hints.
- Dashboard mounts the collapsible **AI Assistant** card (zero network activity until used) and the
  watchdog badge; strict import boundary (props in, one `figureExport` helper out) keeps the module
  extractable to its own repository.
- Added **vitest** (first frontend unit tests, 43 passing) covering the context builder,
  convergence heuristics, SSE parsing, prompt templates, and report assembly; CI now runs
  `npm test` between lint and build.

## v0.3.0 — 2026-07-02

Symmetry analysis in the browser.

- Added a **Detected SG** card beside *Model information* on both the Dashboard and KDE/3D pages: a
  client-side, table-free **FINDSYM-like space-group finder** (ported from the RMC-phonon-dynamics
  project — no spglib/WASM). From the folded `.rmc6f` basis it detects the space group (H–M symbol +
  number), point group, and operation count, and draws an interactive **tolerance ladder** — each
  brick is the space group that holds over a range of atomic-position tolerance; click a rung to
  select it. New frontend modules `symmetry.js` (the finder) and `symmetryModel.js` (structure →
  conventional cell + basis glue); `browserData.structureFromRmc6f` now also returns a
  circular-mean per-reference-number basis.
- The space-group **number table is complete for every producible symbol** (all point-group ×
  allowed-centering combinations), and the headline always shows the point group even when a number
  is unavailable.
- The **Detected SG tolerance selection persists** when switching between the Dashboard and KDE/3D
  pages (shared via context).
- KDE/3D: **contours and log-scale density are on by default.**

## v0.2.0 — 2026-06-26

Figure export, axis labels, and browser-first positioning.

- Added per-figure **Save** controls. Dashboard charts (inline SVG) export as **PNG** (raster) or
  true-vector **SVG**; the KDE/3D panels (canvas + WebGL) export **PNG** at native or **3×**
  resolution. The earlier raster-PDF attempt was dropped — a PDF of a chart isn't true vector, and
  SVG is the proper vector format (convertible to PDF offline if needed).
- Added **Save all figures** in the dashboard's *Loaded N plot files* header. It bundles every
  visible chart into a single `.zip` (PNG or SVG), avoiding the browser's multi-download blocking.
- Added supporting modules: `figureExport.js` (SVG/canvas rasterization + standalone-SVG
  serialization with inlined computed styles), a shared `SaveMenu` component, and a dependency-free
  store-method ZIP writer (`zipArchive.js`, with CRC-32).
- High-resolution panel export re-renders the same drawing onto a 3× offscreen canvas (KDE slice,
  slab) or re-renders the Three.js scene at a higher pixel ratio (`preserveDrawingBuffer`), so the
  output is genuinely higher-resolution rather than upscaled.
- Gave proper vertical-axis labels to plots that previously showed a generic `data`: PDF/partials →
  `G(r)`, x-ray/neutron S(Q) → `S(Q)`, Bragg → `Intensity`. Applied in both the Flask
  `/api/plot/data` endpoint and the browser static-mode parser.
- Reframed the README and roadmap around the hosted browser app (GitHub Pages) as the primary,
  no-install way to use the dashboard; the Flask backend and Python package are now positioned as
  optional local/advanced paths. Added a Phase 7 plan for interactive STOG reduction and guided
  RMCProfile run setup.
- Capitalized **RMCProfile** consistently in the app title, header, and messages.

## v0.1.0 — 2026-06-21 (first tagged release)

First public, releasable version: the reusable Python package, Flask API, interactive dashboard,
server-side + browser (WebGPU/CPU) KDE, and Three.js structure viewer with a draggable slice band.

RMCProfile EXAFS dataset output parsing and large-structure display (added late in the v0.1.0 cycle):

- Added plot-kind detection for RMCProfile EXAFS dataset output files (`*-EXAFS-*_Q_OUTPUT.csv` and
  `*-EXAFS-*_R_OUTPUT.csv`) in both the Python package and browser static mode.
- Added `read_exafs_csv`, which handles Q-output files with a descriptive title row before the column
  header and R-output files with real/imaginary/modulus transform columns.
- Added backend and frontend chart labels for these dataset outputs:
  - Q-space: `k (Å^{-1})` vs `χ(k) k²`.
  - R-space: `r (Å)` vs `FT[χ(k) k²]`.
- Added dashboard file badges/order and tests for parser, plot metadata, file listing, and
  `/api/plot/data`.
- Raised the structure endpoint and Structure page point limits to 1,000,000 and raised the Slab In
  Cell canvas draw cap to match, so moderate structures such as `data/RMC/snao.rmc6f` render all
  returned atoms instead of every other point.

Release engineering added in this version:

- **License:** added an MIT `LICENSE`.
- **Packaging:** added `pyproject.toml` so `rmc_toolkits` is pip-installable (`pip install -e .`) and
  exposes `rmc_toolkits.__version__`.
- **CI:** added `.github/workflows/tests.yml` — runs the Python test suite plus frontend lint/build
  on every push and PR to `main`.
- **Tests:** sample-data-backed tests now skip cleanly when the gitignored GNSe example dataset is
  absent, so the suite is green on a fresh clone and in CI (17 run, 16 skip without the sample).

## 2026-06-21 — Draggable slice band + structure-view color work

- **Slab In Cell drag-to-move slice.** Users can now grab the highlighted band in the *Slab In Cell*
  panel and drag it along the slice axis to set the slice position (`zCenter`) live; thickness is
  unchanged. Implemented by adding an `invert()` to `makePlaneMapper` (recovers plane coordinates
  from a screen point), publishing the band geometry each render into `slabGeometryRef`, and wiring
  pointer handlers on the slab canvas (`grab`/`grabbing` cursors, pointer capture,
  `touch-action: none`).
- **Bug fix:** the drag listeners initially never attached — the effect used `[]` deps but the slab
  canvas only renders after a structure loads (`{structure && ...}`), so `slabCanvasRef` was null at
  mount. Fixed by depending on `[structure]`. Verified end-to-end against `data/5K_try1`: drag moved
  SLICE 0.39 → 0.709 with the band and in-slab atoms updating.
- Earlier in the day: gave each element a distinct color across the slab and 3D views, added a shared
  atom color legend above the structure views, and relocated the atom badges.

## 2026-06-19 — WebGPU-accelerated browser KDE

- Moved the static-mode density-map hot loop (`O(grid² · samples)`, up to ~400M `exp()` per slice)
  onto the GPU. Each grid cell is independent → one WebGPU compute-shader invocation per cell.
- Added `workers/gpuKde.js`: inline WGSL shader, lazily-cached adapter/device init, `computeDensityGpu`
  (writes buffers, dispatches, reads back via `mapAsync`, reshapes to the same `density[y][x]` grid as
  the CPU path), and a `shouldUseGpu` work-size heuristic.
- Refactored `localKdeWorker.js` into `computeDensityCpu` + a GPU-or-CPU branch; made `computeKde`
  and `onmessage` async. Same `{ id, result }` message shape, so `StructurePage.jsx` needed no
  changes (its request-`id` guard already drops stale replies). Result carries `backend: 'gpu'|'cpu'`.
- GPU used only when `grid·grid·samples >= 2_000_000`; init attempted once per worker and cached.
- **Robustness is the design point:** missing `navigator.gpu`, no adapter, device/shader error, lost
  device, or sub-threshold work all fall back to the JS loop with identical output; no GPU failure
  surfaces through `worker.onerror`.
- Verified (Chromium/Metal): GPU(f32) vs CPU(f64) parity ~1.7e-6 relative across grids 16/120/260;
  density step ~59× faster at grid 120 and ~107× faster at grid 260. Still subsamples to 6000 slab
  points and uses f32 — a visualization path, not a substitute for server-side SciPy KDE.

## 2026-06-18 — Live Data and loaded-file controls

- Added optional Live Data monitoring for the local Flask dashboard: polls the selected folder for
  supported-file changes and refreshes metadata/charts without a manual reload.
- Added file mtime + size to `/api/files` so the frontend can detect changes cheaply.
- Added a collapsed `Loaded N plot files` panel; expanding shows detected files as badges that can
  hide/show their chart.
- Aligned local folder-selector copy with the hosted dashboard (`Run folder` / `Select Folder`) and
  hardened the native picker startup so a missing default folder falls back to the nearest existing
  directory.

## 2026-06-16 — Hosted static dashboard

- Added the GitHub Pages workflow (`.github/workflows/pages.yml`); Pages source set to **GitHub
  Actions** so the built React/Vite app is served instead of a Jekyll README page.
- Added static-mode local file loading (`browserData.js`): open the hosted dashboard, select a local
  run folder, and parse RMCProfile CSV/log/STOG outputs with no Python and no upload.
- Extended static mode to parse uploaded `.rmc6f`, populate the model summary, and render the slab
  projection + Three.js 3D view.
- Added `workers/localKdeWorker.js`: off-thread browser KDE capped at 6000 slab points, deterministic
  pseudo-random sampling, contour segments for the existing overlay.
- Kept Flask `/api/kde/slice` as the reference server-side SciPy path.

## 2026-06-15 — Palette and docs

- Three.js atom palette → Nature-style (Ga deep blue, Nb vermillion, Se teal).
- Refreshed `assets/rmc-toolkits-KDE.png`; updated README + frontend docs for the local preview flow.

## 2026-05-27 — UI refresh and interactive dashboard

- Bright/dark theme system via CSS variables with a persisted header toggle.
- Removed the sidebar-first workflow; header now hosts data-path input, Dashboard/KDE nav, and theme
  toggle. App defaults to the `data` path.
- Added `GET /api/plot/data` (parsed series + normalized scientific labels: `χ`, `Å`, `Q (Å⁻¹)`).
- Replaced PNG cards with `InteractivePlot.jsx` (SVG, hover readouts, legend toggles, integer x
  ticks, drag-to-zoom + reset). Simplified the Dashboard to a three-card grid.
- Reworked KDE/3D into three side-by-side panels (XY slice, slab x-z projection, Three.js model);
  2D panels preserve lattice aspect ratios. Added gray slab-edge outlines in the 3D model.

## 2026-05-27 — Real server-side KDE

- Added `rmc_toolkits/kde.py`: loads unit-cell-folded cartesian (Å) positions from a `.rmc6f`
  (optional element filter) and computes an XY `scipy.stats.gaussian_kde` density for a z-slab
  (ported from `src/RMC_KDE.py`). Returns plain arrays (grid, extent, contour polylines via
  `contourpy`, slab count, vmin/vmax); fit subsampled to 6000 points for slider responsiveness.
- Added `GET /api/kde/slice`. `z`/`dz` are cell-edge fractions converted to Å internally; loaded
  positions cached per (path, mtime, element) with `lru_cache`.
- Made the backend port configurable via `RMC_TOOLKITS_PORT` (default 5000).
- Replaced the browser box-blur "density" with the fetched real KDE grid (colormap LUT + contour
  overlay); fetches debounced with `AbortController`. Added bandwidth/colormap/grid/contour/log-scale
  controls and `colormaps.js`. Default z-slice auto-snaps to the densest band on load.
- Added `.venv/`, `__pycache__/`, `*.pyc` to `.gitignore`.

## 2026-05-27 — Package + tests foundation

- Added `rmc_toolkits/parsers.py` (RMC CSV/log/STOG, `.rmc6f`, `Frac*.txt`) and
  `write_frac_from_rmc6f` conversion.
- Added `rmc_toolkits/plots.py` (plot detection, figures, metrics, PNG serialization); backend now
  calls into the package instead of duplicating logic.
- Added backend endpoints: `/api/health`, `/api/files`, `/api/plot`, `/api/plot/metadata`,
  `/api/convert/frac`, `/api/structure`, `/api/kde/slice`; plus a data-root guard
  (`RMC_TOOLKITS_DATA_ROOT`).
- Frontend uses `VITE_API_BASE_URL` (default `http://localhost:5000`); file-path field no longer
  fires a request per keystroke.
- Added the Dashboard all-plots view and the KDE/3D page (element filter, z-slice controls, density
  KDE canvas, Three.js folded unit cell with OrbitControls + translucent slab overlay, and the
  slab-in-cell x-z projection).
- Fixed structure sampling so the sample renders all 52,000 atoms and preserves all 52 reference
  sites (sampling grouped by reference number, not raw stride).
- Added a `unittest` suite under `tests/` covering CSV/log parsing, Rwp, `.rmc6f` metadata +
  conversion, structure loading, plot detection/metadata/PNG, and KDE loading/slice.
