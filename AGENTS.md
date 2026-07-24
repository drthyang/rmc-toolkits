# Agent Guide

Onboarding for AI agents and new contributors. Human users want [README.md](README.md) /
[QuickStart.md](QuickStart.md). This file is the "pick up where we left off" record: architecture,
key files, conventions, and current state. Full chronological history lives in
[docs/CHANGELOG.md](docs/CHANGELOG.md); forward plans in [docs/ROADMAP.md](docs/ROADMAP.md).

## What this project is

Post-processing for **RMCProfile** modeling outputs, in three layers. RMCProfile performs atomistic
configuration optimization under experimental constraints; avoid calling those runs
"refinements" (Rietveld-style parameter refinement is a different workflow). Since 2026-07-17 the
app also covers the *pre*-processing step: **Auto StoG** (`rmc_toolkits/scaling.py` + the
`rmc-autoscale` CLI + the Auto StoG tab) automatically scales measured total-scattering S(Q) and
writes the classic stog/RMCProfile-ready file family — see
[docs/STOG_SCALING_PLAN.md](docs/STOG_SCALING_PLAN.md) for the verified math and validation record.

1. **`rmc_toolkits/`** — pure-Python package (parsing, plots, KDE). The source of truth; new app
   code should call into this, not the legacy scripts.
2. **`web_app/`** — Flask API (`backend/app.py`) + React/Vite SPA (`frontend/`).
3. **`src/`** — original standalone research scripts, kept for CLI workflows.

The same React app ships in two runtime modes:
- **Flask mode** — backend serves the built SPA and provides server-side file browsing, SciPy KDE,
  conversion, and Live Data.
- **Static mode** (`VITE_STATIC_MODE=true`) — GitHub Pages build. No backend; the browser parses
  files locally and computes KDE in a Web Worker (WebGPU + CPU fallback).

## Architecture map

```
rmc_toolkits/
  parsers.py     RMC CSV/log, STOG parsing (incl. stog.inp / STOG xy / .dat headers), .rmc6f metadata + atom iteration, Frac*.txt conversion, structure loading
  plots.py       plot-kind detection, matplotlib figures, Rwp/chi metrics, PNG serialization
  kde.py         unit-cell position loading + server-side gaussian_kde slice (+ contours)
  pca_kde.py     per-site RMC displacement clouds → PCA/thermal-ellipsoid stats + separable 3D KDE volume (source of truth)
  orientation.py displacement-direction distribution: Goldberg (hex + 12 pentagon) sphere tiling, solid-angle histogram, enhancement/z-score, orientation tensor (source of truth)
  transforms.py  Keen-2001 conversions, sine-FT pair, Lorch, omitted-low-Q correction, Fourier filter, stog low-r enforcement (pure numpy, no I/O)
  scaling.py     Auto StoG engine: level sweep, closed-form (a,b) fit, self-consistent filter loop, FZ amplitude mode, estimate_rho0 (density from amplitude-criteria concordance), diagnostics (no I/O)
  scattering.py  Faber-Ziman coefficients from a chemical formula (Sears table, 89 elements)
  scaling_cli.py rmc-autoscale CLI: stog.inp/--data → engine → classic stog output family + provenance JSON (all scaling file I/O)

web_app/backend/app.py    Flask API; data-root guard; per-(path,mtime,…) LRU caches for KDE, PCA-KDE, and scaling; /api/scaling/preview|run share the CLI writer

web_app/frontend/src/
  browserData.js                 static-mode local file parsing + run assembly
  colormaps.js                   colormap LUTs for the KDE canvas
  orientationSphere.js           pure helpers for the orientation sphere (cell mesh/outline typed arrays, relief radii, colorbar gradient) — unit-tested, no Three.js
  useSiteCloud.js                shared hook: .rmc6f text loading, worker/API request routing, per-site ellipsoid table, selected site (one app-lifetime worker → shared parse cache across the PCA Ellipsoid and Orientation pages)
  api.js                         frontend API base URL config (VITE_API_BASE_URL)
  llm/                           experimental AI assistant — local LLM (Ollama/LM Studio) or cloud (OpenAI/Gemini); see its README
    context/                     dashboard state → compact LLM context JSON (symmetry + per-site displacements, pair correlations, char budget)
    provider/client.js           OpenAI-compatible client (models, SSE streaming incl. reasoning, connection hints)
    prompts/                     shared system prompt + chat/watchdog message builders
    watchdog/                    convergence heuristics (source of truth) + LLM-narrated badge hook
    useAssistant.js              shared hook: settings, connection probe/auto-connect, run context
    components/                  AssistantPage (chat-only) + connection bar, settings drawer, ChatView (Thinking panel), WatchdogBadge
  components/
    App.jsx                      shell, run-folder selection, page nav, Live Data
    AutoStogPage.jsx             Auto StoG tab — pre-processing, fully client-side in BOTH runtimes and independent of the run folder: page-local S(Q) upload (± stog.inp) → grouped params (fieldsets w/ descriptions) → worker auto-scale (+ rho0 self-consistency estimate when rho0 is empty) → readout + S(Q)/GK/D(r) plots → zip export. Does NOT call /api/scaling/* (those remain for API/CLI use)
    Dashboard.jsx                all-plots run dashboard
    InteractivePlot.jsx          browser-native SVG plot renderer (hover, legend, drag-zoom)
    PlotViewer.jsx               PNG plot rendering + metadata
    StructurePage.jsx            KDE slice, Slab In Cell, Three.js 3D view  ← most complex component
    PcaKdePage.jsx               PCA Ellipsoid tab: site picker, ADP table, Three.js isosurface + ellipsoid + wall projections + clickable unit-cell structure
    OrientationPage.jsx          Orientation tab (own workspace page — not a PCA product): site picker + OrientationView
    OrientationView.jsx          hex-binned orientation sphere (Three.js flat-shaded cells, amplitude relief), per-cell hover readout, fixed-angle a/b/c | PC1/2/3 mini views (click to snap), colorbar + asymmetry/significance strip
    FileExplorer.jsx             file navigation
  workers/
    localKdeWorker.js            static-mode KDE worker (GPU-or-CPU density map, contours, slab math)
    gpuKde.js                    WGSL compute-shader density map + shouldUseGpu heuristic + cached device init
    pcaKde.js                    static-mode PCA-KDE engine (JS port of pca_kde.py): 3×3 Jacobi eigensolver, per-site clouds, separable volume + projections
    autoScale.js                 static-mode Auto StoG engine (JS port of scaling.py + transforms.py + stog parsers + Faber-Ziman); parity-tested against Python goldens (autoScale.test.js)
    autoScaleWorker.js           off-thread runner for autoScale.js (transferable buffers)
    orientation.js               static-mode displacement-orientation engine (JS port of orientation.py: Goldberg hex+pentagon sphere tiling, solid-angle histogram)
    pcaKdeWorker.js              static-mode PCA-KDE worker (parses clouds once, answers 'sites'/'kde'/'orientation' requests off-thread)
    marchingCubes.js             isosurface extraction over a scalar field (Lorensen-Cline tables) for the Three.js KDE surface
```

## Key conventions & gotchas

- **`z` / `dz` are cell-edge fractions** at the API/slider boundary, converted to Ångström inside
  `kde.py`. Keep that contract when touching KDE code.
- **KDE fit is subsampled to 6000 slab points** (deterministic pseudo-random, to avoid RMC
  atom-order aliasing). Static-mode KDE is a *visualization* path — the server-side SciPy
  `gaussian_kde` is the reference for publication values.
- **GPU KDE must always degrade gracefully.** Missing `navigator.gpu`, no adapter, device/shader
  error, lost device, or sub-threshold work all fall back to the CPU loop with identical output.
  GPU is used only when `grid*grid*samples >= 2_000_000`.
- **PCA-KDE is separable, not approximate.** `pca_kde.py` samples the 3D KDE on a grid aligned with
  the cloud's principal axes; with SciPy's `H = factor²·C` bandwidth, `C` and `H` are both diagonal
  in that frame, so the Gaussian factorizes into three 1D kernels and the volume is their tensor
  product (`N·3·grid` exponentials, contracted via BLAS, instead of `N·grid³`). The result equals
  `scipy.stats.gaussian_kde` to round-off — `tests/test_pca_kde.py` and the JS
  `pcaKde.test.js` both assert exact agreement against the full estimator. `pcaKde.js` is a straight
  port (a 3×3 Jacobi eigensolver stands in for `eigh`); keep the two in sync when touching either.
- **PCA displacement convention**: an atom's offset is `coords − cellIndices/supercell` folded over
  the *supercell* boundary (only that edge wraps), mean-subtracted per reference site, then mapped to
  Cartesian Å through the full supercell `latticeVectors`. Clouds pooled by element are meaningful
  because each site is already centered on its own average position.
- **Orientation histogram bins are hexes + exactly 12 pentagons, area-normalized.** A sphere cannot
  be tiled by hexagons alone (Euler), so `orientation.py` uses a Goldberg polyhedron (10ν²+2 cells)
  and divides every count by the cell's exact solid angle — never plot raw counts, or the 12
  pentagons print the icosahedron onto the map. `enhancement = 4π·density` is 1 for an isotropic
  site; `zScore` comes from raw counts only (never after smoothing). The map is deliberately
  **never antipodally folded** — the +u/−u imbalance (off-centring, odd anharmonicity) is the
  signal the ellipsoid cannot show; `antipodalAsymmetry` (vs its Poisson `…Null` floor) quantifies
  it. `workers/orientation.js` is a straight port — keep the two in sync, same construction order
  so cell indices agree.
- **`StructurePage.jsx` canvases render conditionally** (`{structure && (...)}`). Effects that attach
  listeners to those canvases must depend on `structure` (or the canvas ref), not `[]` — otherwise
  they run at mount before the canvas exists and never attach. (This was the slab-drag bug, fixed
  2026-06-21.)
- **Slab In Cell drag**: cursor→slice mapping inverts the 2D plane projection in `makePlaneMapper`
  (`invert()`); the band geometry is published each render into `slabGeometryRef` for the pointer
  handlers. Drag updates `zCenter` live.
- **Three.js atom palette** is a Nature-style scheme; each element gets a distinct color shared by
  the slab and 3D views via a legend above them.
- **rho0 self-consistency is criteria concordance, not a new fit.** The density-limit amplitude
  depends on rho0 (C2 rows scale with `-4·pi·rho0·r`) — degenerate with the scale from low-r alone;
  the Q→0 Faber-Ziman amplitude does not (needs only ⟨b²⟩ + the measured level). `estimate_rho0`
  (scaling.py, JS port in `autoScale.js` — keep in sync) root-finds
  `a_fz/a_density(rho0) = 1` by fixed-point iteration. It REQUIRES a composition; cross-engine
  parity of the iterated result is bounded by the rtol stopping rule (~1e-4), not single-pass
  transform precision.
- **Backend data-root guard**: relative paths resolve under `RMC_TOOLKITS_DATA_ROOT` (default repo
  root); absolute paths are rejected unless inside the root or a natively-picked folder.
- **`src/llm/` import boundary**: the AI assistant module receives run data **only as props**
  (`runName`, `plotFiles`, `rValueFile`, `structure`, `symmetry`, `liveData`) and must not import
  from the rest of the app except `figureExport.js` (`downloadBlob`/`sanitizeFilename`). Cell math
  is duplicated from `ModelSummary.jsx` on purpose. This keeps the module extractable — don't
  "clean up" the duplication by adding host imports. The R-value series it receives is **ln(χ²)**
  (browserData applies `Math.log`); the context builder labels it so the model reads it correctly.

## Run & test

```bash
# Backend (venv with numpy scipy flask flask-cors matplotlib)
source .venv/bin/activate
RMC_TOOLKITS_PORT=5050 python web_app/backend/app.py

# Frontend dev
cd web_app/frontend
VITE_API_BASE_URL=http://localhost:5050 npm run dev

# Static-mode dev (no backend)
VITE_STATIC_MODE=true npm run dev

# Tests (sample-data-backed tests skip when data/ GNSe sample is absent)
MPLCONFIGDIR=/tmp/rmc_toolkits_matplotlib python -m unittest discover -s tests

# Lint frontend
cd web_app/frontend && npx eslint src

# Frontend unit tests (vitest — src/llm module)
cd web_app/frontend && npm test
```

CI (`.github/workflows/tests.yml`) runs the Python suite + frontend lint/build on every push/PR to
`main`. `.github/workflows/pages.yml` deploys the static dashboard. The `rmc_toolkits` package is
pip-installable (`pip install -e .`, see `pyproject.toml`) and exposes `__version__`.

The repo's sample data lives in `data/` (a GaNb₄Se₈ run). Point the run folder at a subdirectory
containing a `.rmc6f` (e.g. `data/5K_try1`) to exercise the KDE/3D page.

> The machine's Anaconda Python has a broken numpy and no Flask — always use the dedicated `.venv`.

## Current known issues

- **iPhone Safari static mode** (2026-06-17): unreliable after selecting a local run folder; desktop
  static mode works. Likely in the mobile folder-selection / file-enumeration path. Keep the local
  Flask workflow as the supported path for mobile until a unified run-source abstraction lands.
- **Zero-atom `.rmc6f` parse** (2026-06-18): some datasets load plots but show 0 atoms / no KDE in
  static mode — the browser `.rmc6f` parser assumes one exact atom-line format. Make it tolerant to
  RMCProfile variants and validate against the declared atom count.
- `src/RMC_3D.py` imports Mayavi and runs work at import time. `src/STOG_plot.py` also has
  top-level plotting, but STOG should stay hidden until a dedicated preprocessing workflow returns.
- The GNSe example dataset (the tests' and README's reference sample) is gitignored and not in the
  repo. Sample-backed tests skip without it, so CI exercises only logic/synthetic-fixture tests. A
  committed standard example run is still wanted so the full suite runs in CI.

## Next best steps

1. Commit a small standard example run (or trimmed fixtures) so the sample-backed tests run in CI
   instead of skipping.
2. Refactor `src/RMC_plot.py` and `src/RMC_3D.py` into thin wrappers with no import-time work.
   Defer `src/STOG_plot.py` unless preprocessing becomes a visible workflow again.
3. Make browser `.rmc6f` parsing tolerant + diagnostic (fixes the zero-atom mobile issue).
4. Add the z-distribution histogram and global x-z projection panels to match `src/RMC_KDE.py`.
5. Export controls (plot PNG/SVG/CSV, KDE/3D screenshots) and `/api/project/scan` summaries.
6. ~~Thermal-ellipsoid ("PCA_KDE") view~~ **DONE** — `rmc_toolkits/pca_kde.py` (source of truth) +
   `workers/pcaKde.js` (static mode), `/api/pca/sites` + `/api/pca/kde`, and the **Thermal
   Ellipsoids** tab (`PcaKdePage.jsx` + `marchingCubes.js`): site picker, ADP table with a
   non-Gaussianity (excess-kurtosis) readout, Three.js isosurface + p% ellipsoid + PC-axis triad,
   and three PC-plane projections. Verified in both runtimes. Possible follow-ups: element-pooled
   clouds in the picker (the engine already supports `element=`), a per-site ellipsoid overlay in
   the main structure view, PNG/CSV export of the volume, and letting the user compare two sites
   side by side. The displacement-*orientation* view (hex-binned sphere; `orientation.py` +
   `workers/orientation.js` + `/api/pca/orientation` + worker kind `'orientation'`, UI as the
   **Orientation** tab of the main panel in `OrientationView.jsx`) landed 2026-07-24. Reference: Maksim Eremenko's PCA_KDE utilities at
   <https://github.com/MaximEremenko/Utilities/tree/main/RMCProfileUtilities/PCA_KDE>.

See [docs/ROADMAP.md](docs/ROADMAP.md) for the full phased plan.
