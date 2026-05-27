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
- `web_app/frontend/src/components/PlotViewer.jsx`: PNG plot rendering and metadata display.
- `docs/ROADMAP.md`: development roadmap for the larger application.

## Known Limitations

- The old `src/RMC_3D.py` still imports Mayavi and executes visualization at import time. It should be refactored before being reused by the web app.
- `src/STOG_plot.py` still contains top-level plotting code. The new package has basic STOG single-file plotting, but not the full multi-panel STOG workflow yet.
- The Dashboard now renders interactive SVG plots directly from parsed data. The standalone File view still uses the PNG plot endpoint.
- The dashboard renderer is intentionally lightweight and custom. It supports hover, legend toggles, integer x ticks, and x-range drag zoom, but not full Plotly-style pan/selection/export yet.
- The KDE page now uses the real SciPy `gaussian_kde` via `GET /api/kde/slice` (resolved 2026-05-27). Remaining gaps vs. `src/RMC_KDE.py`: the desktop tool also shows a z-distribution histogram and a global x-z projection panel alongside the slice, and it supports non-orthorhombic limits; the web page mirrors only the slab x-z projection so far. The KDE fit is subsampled to 6000 slab points, which is fine for visualization but not an exact full-population estimate.
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

Structure metadata (`/api/structure`), KDE slice (`/api/kde/slice`), and plot-data (`/api/plot/data`) endpoints now exist, and the reusable package layer has initial sample-backed tests. Next: add backend API tests around the Flask endpoints, then refactor `src/RMC_3D.py` and `src/STOG_plot.py` so no analysis module performs work at import time. Remaining viewer work: add the z-distribution histogram and global x-z projection panels to match `src/RMC_KDE.py`, add export/screenshot controls, and keep tightening dashboard interactions if researchers need pan or publication export. Batch run summaries are still future work.
