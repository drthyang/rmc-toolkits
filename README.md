# rmc-toolkits

[![Tests](https://github.com/drthyang/rmc-toolkits/actions/workflows/tests.yml/badge.svg)](https://github.com/drthyang/rmc-toolkits/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Post-processing utilities and a local-first dashboard for **RMCProfile** and **STOG** outputs.

🔒 **Your data stays yours.** Run files are read and rendered entirely in your own browser (or on
your own machine) and are never uploaded to any server — a private, secure way to monitor an RMC
refinement.

🔗 **Hosted dashboard:** [drthyang.github.io/rmc-toolkits](https://drthyang.github.io/rmc-toolkits/)
📖 **New here?** Start with [QuickStart.md](QuickStart.md).

## What's Inside

- **Python package (`rmc_toolkits/`)** — parse RMC/STOG/EXAFS outputs, build plots, convert
  `.rmc6f` files, load folded structures, and compute SciPy KDE slices.
- **Web app (`web_app/`)** — Flask API + React/Vite frontend for loading a run directory,
  inspecting interactive plots, converting structures, and exploring KDE/3D views.
- **Legacy scripts (`src/`)** — original standalone CLI/desktop research scripts.

## Features

- **Run dashboard** — auto-detects RMCProfile outputs in a folder and renders browser-native SVG
  charts with hover readouts, legend toggles, and drag-to-zoom.
- **Figure export** — every chart has a Save badge for **PNG** or true-vector **SVG**, *Save all
  figures* bundles the dashboard into a single `.zip`, and the KDE/3D panels export PNG at native or
  high (3×) resolution.
- **EXAFS outputs** — recognizes `*-EXAFS-*_Q_OUTPUT.csv` and `*-EXAFS-*_R_OUTPUT.csv`, including
  Q-space title rows and R-space real/imaginary/modulus columns, with appropriate `k` and `r` axes.
- **Live Data** — the local Flask app watches the selected folder and refreshes charts when files
  change.
- **KDE / 3D page** — model summary, server-side `scipy.stats.gaussian_kde` density slices, a
  slab-in-cell projection, and a Three.js folded unit-cell view. **Drag the highlighted band in the
  Slab In Cell panel** to move the slice position directly.
- **GPU-accelerated browser KDE** — in static/offline mode a Web Worker computes the density map
  with a WebGPU compute shader when available (with an automatic CPU fallback) — up to ~100× faster
  on the heaviest grids, with numerically identical output.
- **Structure tools** — `.rmc6f` → `Frac_coord_<stem>.txt` conversion and reference-number-preserving
  atom sampling for large structures.

## Screenshots

| Run dashboard | KDE / slab / 3D |
| --- | --- |
| ![dashboard](assets/rmc-toolkits-dashboard.png) | ![KDE](assets/rmc-toolkits-KDE.png) |

## Repository Layout

| Path | Purpose |
| --- | --- |
| `rmc_toolkits/` | Reusable package: parsing, plots, KDE. |
| `web_app/backend/app.py` | Flask API server with data-root guarding. |
| `web_app/frontend/` | React + Vite single-page app. |
| `src/` | Original standalone CLI/desktop scripts. |
| `data/` | Sample run folder for quick testing. |
| `tests/` | `unittest` suite for the package and backend. |
| `docs/` | Changelog, roadmap, architecture notes. |

## Setup

```bash
# Python (from repo root)
python3 -m venv .venv
source .venv/bin/activate
pip install -r web_app/backend/requirements.txt   # web app
pip install -e .                                   # rmc_toolkits package (editable)

# Frontend (Node 20.19+ or 22.12+, for Vite 7)
cd web_app/frontend
npm install
```

## Running The Web App

Build the frontend once, then serve it from Flask:

```bash
cd web_app/frontend && npm run build && cd ../..
source .venv/bin/activate
python web_app/backend/app.py        # http://127.0.0.1:5000/
```

On macOS, port 5000 may be taken by AirPlay Receiver — use another port:

```bash
RMC_TOOLKITS_PORT=5050 python web_app/backend/app.py
```

By default the backend only serves paths under the repo root. To browse another data root:

```bash
RMC_TOOLKITS_DATA_ROOT=/absolute/path/to/data python web_app/backend/app.py
```

In the app, use `Select Folder` to pick a run folder and toggle `Live Data` for auto-refresh.

### Frontend development

```bash
cd web_app/frontend
VITE_API_BASE_URL=http://localhost:5050 npm run dev   # http://localhost:5173/
```

## Hosting The Dashboard

Two supported modes:

- **GitHub Pages (static)** — easiest public link. Users select a local run folder in the browser;
  data stays on their machine, no Python server required. Plot parsing, `.rmc6f` model summaries,
  the slab view, browser-side KDE (WebGPU + CPU fallback), and the 3D view all run client-side.
  Live Data works in Chromium browsers (Chrome, Edge, Arc, Opera) via the File System Access API.
  Deployed automatically from `main` by `.github/workflows/pages.yml`.
- **Flask web service** — full-featured: server-side file browsing, structure sampling, conversion,
  SciPy KDE, and Live Data. The included `Dockerfile` builds the frontend and serves it via
  Flask/Gunicorn:

  ```bash
  docker build -t rmc-toolkits-dashboard .
  docker run --rm -p 5000:5000 rmc-toolkits-dashboard
  ```

  For a public deployment (Render, Fly.io, Railway, a VPS), the container honors the provider's
  `PORT` and falls back to `5000`. Set `RMC_TOOLKITS_DATA_ROOT` to expose your run folders.

To test the static build locally:

```bash
cd web_app/frontend
VITE_STATIC_MODE=true VITE_BASE_PATH=/ npm run build
npm run preview
```

## Python Package Usage

```python
from rmc_toolkits import (
    detect_plot_kind, kde_slice, load_unit_cell_positions,
    make_plot, plot_to_png, read_exafs_csv, read_structure, write_frac_from_rmc6f,
)

frac_path = write_frac_from_rmc6f("data/GNSe.rmc6f", overwrite=True)
structure = read_structure("data")

positions = load_unit_cell_positions("data/GNSe.rmc6f", element="Ga")
payload = kde_slice(
    positions.positions,
    z_center=0.5 * positions.cell_lengths[2],
    dz=0.08 * positions.cell_lengths[2],
    xlim=(0.0, float(positions.cell_lengths[0])),
    ylim=(0.0, float(positions.cell_lengths[1])),
)

png_bytes = plot_to_png(make_plot("data/GNSe_FQ1.csv"))
```

Lower-level parser helpers are also exported: `read_rmc_csv`, `read_exafs_csv`, `read_chi`, `read_stog`,
`read_atom_indices`, `read_cell_vectors`, `iter_rmc6f_atoms`, `frac_lines_from_rmc6f`, `rwp`.

## Backend API

All endpoints are under `/api`. Relative paths resolve under `RMC_TOOLKITS_DATA_ROOT`; absolute
paths are rejected unless inside the configured root or a folder selected via the native picker.

| Method & path | Description |
| --- | --- |
| `GET /api/health` | Service status and active data root. |
| `GET /api/files?dir=data` | Directory listing with detected `plotKind` values. |
| `POST /api/dialog/folder` | Open a native folder picker; register it as an allowed root. |
| `GET /api/plot?path=...` | Render one supported file as a PNG. |
| `GET /api/plot/metadata?path=...` | Plot kind, title, and metrics (`rwp`, `final_chi_r`). |
| `GET /api/plot/data?path=...` | Parsed plot series and normalized axis labels for SVG rendering. |
| `POST /api/convert/frac` | Convert `.rmc6f` → `Frac_coord_<stem>.txt`. |
| `GET /api/structure?dir=...&maxPoints=1000000` | Folded atom positions, counts, lattice. |
| `GET /api/kde/slice?dir=...&element=Ga&z=0.5&dz=0.08` | KDE density grid, contours, slab counts. |

## Supported File Patterns

- Real-space PDF/G(r): `*_FT_XFQ1.csv`, `*PDF*.csv`
- Reciprocal-space S(Q): `*_FQ1.csv`, `*_SQ1.csv`
- EXAFS: `*-EXAFS-*_Q_OUTPUT.csv` (`k` vs `χ(k) k²`) and `*-EXAFS-*_R_OUTPUT.csv` (`r` vs Fourier-transform components)
- Bragg profiles: `*_bragg.csv`
- R-value logs: `*.log`
- Basic STOG outputs: `scale_ft.gr`, `scale_ft.sq`, `scale_ft_rmc.fq`
- Structure files: `*.rmc6f`, `Frac*.txt`

Most RMCProfile CSV parsers expect first-row labels followed by numeric rows. EXAFS Q-output files
may include a descriptive title row before the column header; `read_exafs_csv` handles that layout.

## Legacy CLI Scripts

```bash
pip install numpy matplotlib scipy seaborn
python src/RMC_plot.py --dir data [--save --no-show]
python src/RMC_KDE.py [--el Mn]
python src/RMC_3D.py            # needs mayavi
```

`RMC_KDE.py` and `RMC_3D.py` expect `Frac*.txt` plus `.rmc6f` in the working directory;
`STOG_plot.py` expects `stog_input.dat` and STOG output files.

## Tests

```bash
source .venv/bin/activate
MPLCONFIGDIR=/tmp/rmc_toolkits_matplotlib python -m unittest discover -s tests
```

## Status

The reusable Python package, Flask API, interactive dashboard, EXAFS plotting support, server-side
KDE, and Three.js viewer are in place. Current priorities: richer run summaries, export/report
workflows, and refactoring the legacy scripts into thin wrappers around `rmc_toolkits`. See
[docs/ROADMAP.md](docs/ROADMAP.md) for the full plan and [docs/CHANGELOG.md](docs/CHANGELOG.md) for
history.

## License

Released under the [MIT License](LICENSE) © 2026 Tsung-Han Yang.
