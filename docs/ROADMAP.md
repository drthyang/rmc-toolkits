# Development Roadmap

## Product Vision

Build a sophisticated local-first analysis app for RMCProfile and STOG workflows. The app should help researchers import a run directory, inspect detected outputs, compare fits, explore atomic structures, generate KDE slices, monitor R-values, and export publication-ready figures and reproducible reports.

## Phase 1: Foundation

Goal: make the project reliable enough to build on.

- Package reusable analysis logic under `rmc_toolkits/`. **Done for the active parser, plot, and KDE paths.**
- Keep CLI scripts as thin wrappers over package functions.
- Add and maintain tests for:
  - RMC CSV parsing.
  - Log chi parsing.
  - Rwp calculation.
  - Plot-kind detection.
  - `.rmc6f` lattice and atom parsing.
  - `.rmc6f` to `Frac_coord*.txt` conversion.
- Add a small fixture-based smoke test using `data/`.
- Add formatting and linting commands for Python and frontend code.
- Document app startup, data-root behavior, and expected file patterns.

## Phase 2: Project Workspace

Goal: move from file-by-file plotting to project-level analysis.

- Add a project scanner that detects:
  - RMC CSV outputs.
  - PDF, S(Q), Bragg, partials, and logs.
  - `.rmc6f` and `Frac*.txt` structure files.
  - STOG input/output files.
- Create project summary JSON with file roles, available plots, metrics, lattice metadata, element list, and warnings.
- Add frontend workspace layout:
  - header-level data path controls,
  - dashboard page for run plots and model summary,
  - KDE / 3D page for structure exploration.
- Expand the current dashboard with richer run metadata, warnings, and optional comparison summaries.
- Add run comparison support for multiple directories.

## Phase 3: Interactive Analysis

Goal: replace static inspection with interactive scientific workflows.

- Continue with the current browser-native SVG plots unless Plotly-level pan/export features become
  necessary for researcher workflows.
- Add controls for axes, range, residual display, log scale, and export resolution.
- Extend the KDE slice viewer:
  - element selector,
  - z-position slider,
  - slab thickness control,
  - bandwidth control,
  - colormap selector,
  - contour toggle.
- Maintain the current server-side SciPy KDE path and add export/reproducibility controls.
- Maintain the hosted GitHub Pages path with browser-side parsing and Web Worker KDE for users who
  want to inspect local files without installing the package.
- Extend the Three.js structure viewer:
  - element visibility toggles,
  - configurable atom coloring,
  - unit cell/supercell display,
  - camera presets,
  - screenshot export.

## Phase 4: Background Jobs

Goal: make expensive analysis responsive and reproducible.

- Add a job model for KDE, structure transforms, batch plots, and report generation.
- Start with SQLite-backed local jobs.
- Track:
  - job status,
  - input paths,
  - parameters,
  - output artifacts,
  - runtime,
  - error messages.
- Add frontend job status indicators and retry controls.
- Cache expensive computed arrays and generated plots.

## Phase 5: Reporting And Export

Goal: make the app useful at the end of a research session.

- Export individual plots as PNG, SVG, and CSV.
- Export project summary as JSON.
- Generate a reproducible report containing:
  - input directory,
  - detected files,
  - software version,
  - plot parameters,
  - Rwp metrics,
  - lattice metadata,
  - selected figures.
- Add figure presets for manuscript, talk, and notebook usage.

## Phase 6: Lab-Ready App

Goal: make the tool safe and pleasant for broader use.

- Add authentication if served beyond localhost.
- Keep data-root restrictions enabled by default.
- Add project persistence and recent projects.
- Package as a one-command local app.
- Add robust error messages for malformed files.
- Add documentation with example workflows and screenshots.

## Architecture Target

- `rmc_toolkits/`: pure Python package for parsing, analysis, plotting, and structure transforms.
- `web_app/backend/`: API server, project scanner, jobs, artifact storage.
- `web_app/frontend/`: React app with project workspace, interactive plots, KDE, structure viewer, and exports.
- `.github/workflows/pages.yml`: builds the static GitHub Pages dashboard from `web_app/frontend`.
- `data/`: small example fixtures.
- `docs/`: hand-off records, roadmap, workflow notes, and architecture decisions.

## Suggested Immediate Backlog

1. Add backend API tests around `/api/files`, `/api/plot/data`, `/api/structure`, and `/api/kde/slice`.
2. Refactor `src/RMC_plot.py` into a CLI wrapper around `rmc_toolkits.plots`.
3. Refactor `src/STOG_plot.py` so all top-level plotting becomes callable functions.
4. Refactor `src/RMC_3D.py` to avoid Mayavi import and execution at import time.
5. Add `/api/project/scan` for directory-level summaries.
6. Add project-level warnings for missing expected files, malformed outputs, and unsupported STOG layouts.
7. Add export controls for plot PNG/SVG/CSV and KDE/3D screenshots.
8. Move heavier static-mode parsing/KDE work toward transferable typed arrays and profile large `.rmc6f` files.
9. Add recent-project persistence for local desktop use.
