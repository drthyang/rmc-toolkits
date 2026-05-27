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

## Important Files

- `rmc_toolkits/parsers.py`: parsing and structure-loading functions.
- `rmc_toolkits/plots.py`: reusable plot builders and plot-kind detection.
- `web_app/backend/app.py`: Flask API layer.
- `web_app/frontend/src/api.js`: frontend API base URL config.
- `web_app/frontend/src/components/FileExplorer.jsx`: file navigation.
- `web_app/frontend/src/components/Dashboard.jsx`: all-plots run dashboard.
- `web_app/frontend/src/components/StructurePage.jsx`: KDE slice and 3D model page.
- `web_app/frontend/src/components/PlotViewer.jsx`: PNG plot rendering and metadata display.
- `docs/ROADMAP.md`: development roadmap for the larger application.

## Known Limitations

- The old `src/RMC_3D.py` still imports Mayavi and executes visualization at import time. It should be refactored before being reused by the web app.
- `src/STOG_plot.py` still contains top-level plotting code. The new package has basic STOG single-file plotting, but not the full multi-panel STOG workflow yet.
- The web app still renders static PNG plots. Interactive Plotly/Three.js visualizations are a future milestone.
- The dashboard currently uses static PNG plot cards. Replacing those cards with interactive Plotly panels is still a future milestone.
- The KDE page currently renders a browser-side density-style slice from sampled, unit-cell-folded structure points. It is not yet the full SciPy `gaussian_kde` implementation from `src/RMC_KDE.py`.
- There are no automated tests yet. Add tests around the package layer before expanding the app.
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

Optional data root:

```bash
RMC_TOOLKITS_DATA_ROOT=/absolute/path/to/rmc/data python web_app/backend/app.py
```

## Next Best Engineering Step

Add package-level tests using the sample files in `data/`, then refactor `src/RMC_3D.py` and `src/STOG_plot.py` so no analysis module performs work at import time. After that, build API endpoints for structure metadata, atom coordinates, KDE slices, and batch run summaries.
