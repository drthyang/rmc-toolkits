# RMC Toolkits Hand-Off Record

Date: 2026-05-27

## Current State

This repository contains research utilities for post-processing RMCProfile and STOG outputs. The original codebase had script-first plotting tools in `src/` and a lightweight Flask/React web viewer in `web_app/`.

This hand-off adds a reusable package layer in `rmc_toolkits/` and wires the web backend to it. The scripts in `src/` are still present for familiar CLI/research workflows, but new app development should depend on `rmc_toolkits` instead of importing script entry points directly.

## What Changed

- Added `rmc_toolkits/parsers.py` for reusable RMC CSV, log, STOG, `.rmc6f`, and `Frac*.txt` parsing.
- Added `.rmc6f` to `Frac_coord*.txt` conversion via `write_frac_from_rmc6f`.
- Added `rmc_toolkits/plots.py` for plot detection, figure generation, metric calculation, and PNG serialization.
- Replaced duplicated backend plotting logic in `web_app/backend/app.py` with calls into `rmc_toolkits`.
- Added backend endpoints:
  - `GET /api/health`
  - `GET /api/files`
  - `GET /api/plot`
  - `GET /api/plot/metadata`
  - `POST /api/convert/frac`
  - `GET /api/structure`
  - `GET /api/kde/slice`
- Added a configured data root guard. By default the app can only browse files under the project root. Override with `RMC_TOOLKITS_DATA_ROOT=/path/to/data`.
- Updated frontend API calls to use `VITE_API_BASE_URL` with `http://localhost:5000` as the default.
- Updated the frontend file explorer so typing in the path field does not trigger a request on every keystroke.
- Added plot metadata display for titles and numeric metrics such as `rwp`.
- Added a browser action for `.rmc6f` files to generate matching `Frac_coord_<stem>.txt` files.
- Added a dashboard workspace view that renders all detected plots in the current run folder.
- Added a KDE / 3D page with element filtering, z-slice controls, a metric-aspect density-style KDE slice canvas, and a Three.js folded-unit-cell structure viewer.
- Added OrbitControls to the Three.js view for drag/pan/zoom interaction, plus a translucent 3D slab overlay tied to z/dz.
- Added a slab-in-cell panel next to the KDE slice showing an x-z projection and the current slab band as z/dz changes.
- Fixed structure sampling so the GNSe sample renders all 52,000 atoms and preserves all 52 reference-number sites. For larger datasets, backend sampling is grouped by reference number rather than by raw stride.

## 2026-05-27 Update: Real Server-Side KDE

- Added `rmc_toolkits/kde.py`: a reusable compute layer that loads unit-cell-folded cartesian (Angstrom) positions from a `.rmc6f` file (optional element filter) and computes an XY `scipy.stats.gaussian_kde` density for a z-slab. Ports the slab math from `src/RMC_KDE.py`. Returns plain arrays so the frontend owns colormap/contour styling: density grid, plot extent, contour polylines (via `contourpy`), slab atom count, and vmin/vmax. The KDE fit is subsampled to 6000 slab points to keep slider interaction responsive.
- Added the `GET /api/kde/slice` backend endpoint. `z` and `dz` are passed as fractions of the cell edge (matching the existing slider semantics) and converted to Angstrom internally. Loaded positions are cached per (rmc6f path, mtime, element) with an `lru_cache`.
- Made the backend port configurable via `RMC_TOOLKITS_PORT` (defaults to 5000) so a test instance can run alongside an existing dev server.
- Replaced the browser-side box-blur "density" in `StructurePage.jsx` with the fetched real KDE grid, rendered to canvas with a colormap lookup table plus a contour overlay. KDE fetches are debounced and use an `AbortController` to cancel in-flight requests on rapid slider changes.
- Added KDE controls: bandwidth, colormap (viridis/magma/seismic/reds/greys), grid resolution, contour toggle, and log-scale toggle. Colormap and contour visibility are pure client-side re-renders; bandwidth, grid, element, z, dz, and log trigger a recompute.
- Added `web_app/frontend/src/colormaps.js` with interpolated colormap LUTs.
- The default z-slice now auto-snaps to the densest z-band on load, because the geometric cell midpoint can fall in an empty gap between atomic layers (as it does for the GNSe sample).
- Added `.venv/`, `__pycache__/`, and `*.pyc` to `.gitignore`.

## 2026-05-27 Update: Package Tests

- Added a standard-library `unittest` suite under `tests/`.
- Covered sample-file parsing for RMC CSV outputs, RMC log chi values, Rwp calculation, `.rmc6f` lattice and atom metadata, `.rmc6f` to `Frac_coord*.txt` conversion, full folded structure loading, plot-kind detection, plot metadata/PNG serialization, and KDE position loading/slice computation.
- Documented the test command in `README.md`.

## 2026-05-27 Update: UI Refresh And Interactive Dashboard

- Added a bright/dark theme system using CSS variables, with a persisted theme toggle in the app header.
- Removed the old sidebar-first workflow from the primary app shell. The app now defaults to the `data` path and exposes the data-path input, Dashboard/KDE navigation, and theme toggle directly in the header.
- Added `GET /api/plot/data`, which returns parsed plot series and normalized scientific labels for browser-native rendering. Labels now use `χ`, `Å`, and `Q (Å^{-1})` with frontend superscript rendering.
- Replaced static dashboard PNG cards with `InteractivePlot.jsx`, a lightweight SVG renderer with hover readouts, legend toggles, integer x-axis ticks, and drag-to-zoom x-range selection with a reset button.
- Simplified the Dashboard by removing the summary tiles and arranging plots in a three-card desktop grid.
- Reworked KDE/3D into three side-by-side panels on wide screens: KDE XY slice, slab-in-cell x-z projection, and Three.js model. The 2D panels preserve lattice-parameter aspect ratios; the 3D model uses lattice-scaled positions.
- Added gray slab-edge outlines in the Three.js model to make the current z/dz slab boundaries visible.

## 2026-06-15 Update: KDE / 3D Palette And Docs

- Updated the Three.js atom colors in `StructurePage.jsx` to a Nature Publishing Group-style
  palette: Ga deep blue, Nb vermillion, and Se teal.
- Refreshed `assets/rmc-toolkits-KDE.png` with the current KDE / slab / 3D layout.
- Updated `README.md` and frontend documentation to describe the current local preview workflow.

## 2026-06-16 Update: Hosted Static Dashboard

- Added a GitHub Pages deployment workflow at `.github/workflows/pages.yml`. The repository Pages
  source should be set to **GitHub Actions** so `https://drthyang.github.io/rmc-toolkits/` serves
  the built React/Vite dashboard instead of the README/Jekyll page.
- Added static-mode local file loading in `web_app/frontend/src/browserData.js`. Users can open the
  hosted dashboard directly in a browser, select a local run folder, and parse supported RMCProfile
  CSV/log/STOG outputs without installing the Python package or uploading data.
- Extended static mode to parse uploaded `.rmc6f` files, populate the model summary, render the
  slab-in-cell projection, and show the Three.js 3D structure view.
- Added `web_app/frontend/src/workers/localKdeWorker.js`: a browser-side Gaussian KDE worker for the
  hosted dashboard. It evaluates KDE off the UI thread, caps the fit population at 6000 slab points,
  uses deterministic pseudo-random sampling to avoid RMC atom-order aliasing, and emits contour
  segments for the existing overlay renderer.
- Kept the Flask app as the reference path for server-side SciPy KDE through `/api/kde/slice`.
  Browser KDE is useful for the hosted tool and local privacy-preserving inspection, but it is still
  slower than the local Flask/SciPy path and should be treated as a visualization workflow rather
  than a replacement for validated publication calculations.
- Added the hosted dashboard link near the top of `README.md` and tightened the app header so the
  Dashboard and KDE/3D tabs sit closer to the RMCprofile Run Monitor logo.

## 2026-06-18 Update: Live Data And Loaded-File Controls

- Added optional Live Data monitoring for the local Flask dashboard. When enabled, the dashboard
  polls the selected folder for supported file changes and refreshes plot metadata/charts without a
  manual reload.
- Added file modification time and size to `/api/files` responses so the frontend can detect changes
  cheaply.
- Added a collapsed `Loaded x plot files` dashboard panel. Expanding it shows detected plot files as
  badges, and each badge can hide or show its chart.
- Updated the local folder selector copy to match the hosted dashboard more closely: `Run folder`
  and `Select Folder`.
- Hardened the native folder picker startup path so a missing default folder falls back to the
  nearest existing directory before opening the OS dialog.

## 2026-06-19 Update: WebGPU-Accelerated Browser KDE

- Addressed slow browser-side KDE on mobile by moving the density-map hot loop
  (`O(grid^2 * samples)`, up to ~400M `exp()` evaluations per slice) onto the GPU. Each grid cell
  is independent, so it maps cleanly to a WebGPU compute shader with one invocation per cell.
- Added `web_app/frontend/src/workers/gpuKde.js`: an inline WGSL compute shader, a lazily-cached
  adapter/device initializer, `computeDensityGpu` (writes sample/param/output buffers, dispatches,
  reads back via `mapAsync`, reshapes to the same `density[y][x]` grid the CPU path returns), and a
  `shouldUseGpu` work-size heuristic.
- Refactored `localKdeWorker.js`: extracted the existing loop into `computeDensityCpu`, added a
  GPU-or-CPU branch (`density = gpuResult ?? computeDensityCpu(args)`), and made `computeKde` and the
  `onmessage` handler async. The worker still posts the same `{ id, result }` message, so
  `StructurePage.jsx` needed no changes; its existing request-`id` guard already drops stale async
  replies during rapid slider drags. The result now carries a `backend: 'gpu' | 'cpu'` diagnostic tag.
- The GPU path is used only when `grid * grid * sampleCount >= 2_000_000` (`GPU_MIN_WORK`); smaller
  slices stay on the CPU to avoid GPU setup/readback overhead. Init is attempted once per worker and
  the result cached, so unsupported devices never re-probe on every message.
- **Robustness is the design point.** Missing `navigator.gpu`, no adapter, a device/shader error, a
  lost device, or sub-threshold work all fall back to the existing JS loop. Output is identical and
  no GPU failure ever surfaces through `worker.onerror`. Devices without WebGPU behave exactly as
  before; modern iOS Safari (26+) and Android Chrome (121+) ship WebGPU on by default and get the win.
- Verified in-browser (Chromium/Metal): numeric parity between GPU (f32) and CPU (f64) is
  ~1.7e-6 relative across grids 16/120/260 — well under the 1e-4 tolerance — and the density step
  measured ~59× faster at grid 120 (1.9 ms vs 112.6 ms) and ~107× faster at grid 260
  (4.3 ms vs 461.8 ms). Confirmed the fallback returns cleanly (cached null, no retry storm) and the
  worker reports `backend: 'cpu'` below the threshold, with no console errors.
- The GPU acceleration speeds up the *density-map computation* only; the worker still subsamples to
  6000 slab points and uses f32 on the GPU, so static-mode KDE remains a visualization workflow, not
  a substitute for the server-side SciPy `gaussian_kde` reference path.

## Important Files

- `rmc_toolkits/parsers.py`: parsing and structure-loading functions.
- `rmc_toolkits/plots.py`: reusable plot builders and plot-kind detection.
- `rmc_toolkits/kde.py`: server-side KDE slice computation (positions loading + `gaussian_kde`).
- `web_app/backend/app.py`: Flask API layer.
- `web_app/frontend/src/colormaps.js`: colormap lookup tables for the KDE canvas.
- `web_app/frontend/src/api.js`: frontend API base URL config.
- `web_app/frontend/src/components/FileExplorer.jsx`: file navigation.
- `web_app/frontend/src/components/Dashboard.jsx`: all-plots run dashboard.
- `web_app/frontend/src/components/InteractivePlot.jsx`: browser-native SVG plot renderer for the dashboard.
- `web_app/frontend/src/components/StructurePage.jsx`: KDE slice and 3D model page.
- `web_app/frontend/src/browserData.js`: static-mode local file parsing and run assembly.
- `web_app/frontend/src/workers/localKdeWorker.js`: browser-side KDE worker for GitHub Pages/static mode (GPU-or-CPU density map, contours, slab math).
- `web_app/frontend/src/workers/gpuKde.js`: WebGPU compute-shader density map plus `shouldUseGpu` heuristic and cached device init; falls back to the worker's CPU loop when WebGPU is unavailable.
- `web_app/frontend/src/components/PlotViewer.jsx`: PNG plot rendering and metadata display.
- `docs/ROADMAP.md`: development roadmap for the larger application.

## Known Limitations

- **Known issue, 2026-06-17:** GitHub Pages/static mode currently does not run reliably on iPhone
  Safari/mobile browsers after selecting a local run folder. Desktop static mode works. The failure
  is likely in the mobile browser folder-selection/file-enumeration path rather than the dashboard
  rendering itself, but this still needs a focused device-side debug pass. Keep the local
  Flask/backend workflow as the supported path for mobile-adjacent use until a unified import/source
  abstraction is implemented and verified on iOS.
- **Known issue, 2026-06-18:** Some datasets load plots on iPhone/static mode but show `0` total
  atoms and no KDE/3D structure. This suggests the `.rmc6f` file is being selected/read, but the
  browser-side atom parser is skipping all atom rows. The current parser assumes an exact atom-line
  format after `Atoms:`: atom number, element, type label, x/y/z, reference number, and cell indices.
  Future work should make `.rmc6f` parsing tolerant to alternate RMCProfile variants, support rows
  with or without a type-label column, validate parsed atom count against the declared `Number of
  atoms`, and surface a clear parse error when zero atoms are parsed.
- The old `src/RMC_3D.py` still imports Mayavi and executes visualization at import time. It should be refactored before being reused by the web app.
- `src/STOG_plot.py` still contains top-level plotting code. The new package has basic STOG single-file plotting, but not the full multi-panel STOG workflow yet.
- The Dashboard renders interactive SVG plots directly from parsed data. The PNG plot endpoint is
  still available for API consumers and future export workflows.
- The dashboard renderer is intentionally lightweight and custom. It supports hover, legend toggles, integer x ticks, and x-range drag zoom, but not full Plotly-style pan/selection/export yet.
- The KDE page uses the real SciPy `gaussian_kde` via `GET /api/kde/slice` (resolved 2026-05-27).
  Remaining gaps vs. `src/RMC_KDE.py`: the desktop tool also shows a z-distribution histogram and
  a global x-z projection panel alongside the slice, and it supports non-orthorhombic limits; the
  web page mirrors only the slab x-z projection so far. The KDE fit is subsampled to 6000 slab
  points, which is fine for visualization but not an exact full-population estimate.
- GitHub Pages/static mode has no Python backend. It uses local browser parsing and a Web Worker KDE
  implementation whose density map runs on the GPU via a WebGPU compute shader when available (with
  an automatic CPU fallback), so the heaviest slices are far faster than the old CPU-only loop. On
  devices without WebGPU it uses the CPU path and can still trail the local Flask/SciPy app on very
  large structures. The worker uses deterministic random subsampling when a slab exceeds 6000 points
  to avoid biased densities from file-order stride sampling.
- Live Data monitoring is only available in the local Flask app. GitHub Pages/static mode cannot
  watch filesystem changes after a folder is selected.
- The native folder picker is intended for local desktop sessions. It depends on OS-level dialog
  availability and is not expected to work when the app is served remotely or inside restricted
  browser/sandbox environments.
- The test suite currently targets the reusable package layer and sample fixtures. Expand it around backend API behavior and edge-case fixture files before broadening app features.
- The backend is Flask. It is acceptable for the current local viewer, but FastAPI would be a better fit for typed analysis APIs and background job status.
- The `.rmc6f` to `Frac_coord*.txt` converter preserves the current observed format exactly for the sample data. Add fixtures for non-cubic or unusual RMCProfile outputs before relying on it for every dataset.

## Local Run Notes

Backend:

```bash
cd web_app/backend
python app.py
```

Frontend:

```bash
cd web_app/frontend
npm install
npm run dev
```

Optional data root and port:

```bash
RMC_TOOLKITS_DATA_ROOT=/absolute/path/to/rmc/data RMC_TOOLKITS_PORT=5050 python web_app/backend/app.py
```

To point the frontend at a non-default backend, set `VITE_API_BASE_URL` (e.g. in `web_app/frontend/.env.local`).

Note: the machine's Anaconda Python has a broken numpy and no flask. Use a dedicated venv (`python3.13 -m venv .venv`) with `numpy scipy flask flask-cors matplotlib` to run the backend.

## Next Best Engineering Step

Structure metadata (`/api/structure`), KDE slice (`/api/kde/slice`), and plot-data (`/api/plot/data`) endpoints now exist, and the reusable package layer has initial sample-backed tests. Next: add backend API tests around the Flask endpoints, then refactor `src/RMC_3D.py` and `src/STOG_plot.py` so no analysis module performs work at import time. Remaining viewer work: continue polishing the dashboard UI, especially the model-information summary and plot-card density; add the z-distribution histogram and global x-z projection panels to match `src/RMC_KDE.py`; add export/screenshot controls; and keep tightening dashboard interactions if researchers need pan or publication export. Batch run summaries are still future work.
