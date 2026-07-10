# rmc-toolkits

[![Tests](https://github.com/drthyang/rmc-toolkits/actions/workflows/tests.yml/badge.svg)](https://github.com/drthyang/rmc-toolkits/actions/workflows/tests.yml)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](LICENSE)

A **browser-first dashboard** for inspecting **RMCProfile** modeling run folders. Open the hosted
app, select a run directory, and review plots, model information, KDE slices, symmetry, and 3D
structure views without installing anything.

### ▶️ Open the app — [drthyang.github.io/rmc-toolkits](https://drthyang.github.io/rmc-toolkits/)

1. Visit the link above.
2. Click **Select Folder** and choose your RMCProfile run directory.
3. Your plots, KDE slices, symmetry summary, and folded 3D structure render in your browser.

🔒 **Your raw run files never leave your device.** They are read and rendered entirely in your
browser and are never uploaded to rmc-toolkits or any project server — a private, secure way to
monitor an RMC modeling run. (Your browser's picker may say “Upload”, but nothing is sent
anywhere.)[^cloud-llm-privacy]

⚡ **Live monitoring** auto-refreshes charts as new files are written, in Chromium browsers (Chrome,
Edge, Arc, Opera).

📖 **New here?** Start with [QuickStart.md](QuickStart.md). Want server-side features or a Python
API? See [Run locally](#run-locally-optional) and [Python package usage](#python-package-usage).

## What's Inside

- **Web app (`web_app/`)** — a React/Vite single-page app for the hosted, no-install dashboard. An
  **optional** Flask backend adds server-side file browsing, `.rmc6f` conversion, and
  reference-grade SciPy KDE for local or self-hosted use.
- **Python package (`rmc_toolkits/`)** — reusable parsing, plotting, `.rmc6f` conversion, folded
  structure loading, and SciPy KDE helpers for RMCProfile outputs, including RMCProfile EXAFS dataset
  CSVs.
- **Legacy scripts (`src/`)** — original standalone CLI/desktop research scripts.

## Features

- **Run dashboard** — auto-detects supported RMCProfile files and renders browser-native SVG
  charts for PDFs/G(r), S(Q), Bragg profiles, R-value logs, partials, and RMCProfile EXAFS dataset
  Q/R CSVs.
- **Interactive plots and export** — hover readouts, legend toggles, drag-to-zoom, per-chart PNG/SVG
  export, and a *Save all figures* `.zip`; KDE, slab, and 3D panels export PNG at native or 3×
  resolution.
- **Live Data** — auto-refreshes charts as your RMC modeling run writes new files. The hosted app
  uses the File System Access API in Chromium browsers; the optional Flask backend watches folders
  server-side.
- **Structure workspace** — reads folded `.rmc6f` structures, summarizes the model, draws KDE density
  slices, shows a slab-in-cell projection, and renders a Three.js unit-cell view. Drag the
  highlighted slab band to move the slice position directly.
- **Browser KDE with fallback** — computes density maps in a Web Worker, using WebGPU when available
  and falling back automatically to CPU. The optional backend keeps the reference SciPy KDE path for
  local use.
- **Symmetry analysis** — a client-side, FINDSYM-like *Detected SG* panel reports space group, point
  group, operation count, and tolerance-dependent changes from the folded structure.
- **AI Assistant (beta)** — optional chat over the loaded run with a local LLM by default (Ollama or
  LM Studio) or an opt-in cloud model (OpenAI, Gemini). It sends compact run context such as Rwp,
  ln(χ²), symmetry, and pair-correlation cues; raw files are not uploaded to
  rmc-toolkits.[^cloud-llm-privacy] See
  [`web_app/frontend/src/llm/README.md`](web_app/frontend/src/llm/README.md) for setup.

[^cloud-llm-privacy]: If you opt into a cloud LLM, the compact summarized run context used for
    assistant responses is sent directly to the cloud LLM server you selected. Raw run files are not
    uploaded to rmc-toolkits.

## Screenshots

| Run dashboard | KDE / slab / 3D |
| --- | --- |
| ![dashboard](assets/rmc-toolkits-dashboard-demo.png) | ![KDE](assets/rmc-toolkits-kde-demo.png) |

### AI Assistant

Ask about the loaded run in plain language. The run's metrics, symmetry, and convergence history
travel with every message, so answers quote the actual numbers — and reasoning models stream their
chain of thought in a collapsible *Thinking* panel. Below, a local model (LM Studio) summarizes the
bundled demo run as a set of tables.

![AI Assistant summarizing a run as tables](assets/rmc-toolkits-assistant.gif)

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

> Only needed to run the Flask backend, use the Python package, or develop locally. **To just use
> the dashboard, open the [hosted app](https://drthyang.github.io/rmc-toolkits/) — nothing to install.**

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

## Run Locally (optional)

Most users don't need this — the [hosted app](https://drthyang.github.io/rmc-toolkits/) covers
monitoring and visualization entirely in the browser. Run the Flask backend only when you want
server-side file browsing, `.rmc6f` conversion, reference-grade SciPy KDE, or to self-host on a
network.

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

Two ways to deploy:

- **GitHub Pages (static)** — how the public app at
  [drthyang.github.io/rmc-toolkits](https://drthyang.github.io/rmc-toolkits/) is served, and the
  recommended way to share it. Users select a local run folder in the browser; data stays on their
  machine, no Python server required. Plot parsing, `.rmc6f` model summaries, the slab view,
  browser-side KDE (WebGPU + CPU fallback), and the 3D view all run client-side. Live Data works in
  Chromium browsers (Chrome, Edge, Arc, Opera) via the File System Access API. Deployed automatically
  from `main` by `.github/workflows/pages.yml`.
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

Lower-level parser helpers are also exported: `read_rmc_csv`, `read_exafs_csv`, `read_chi`,
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
- RMCProfile EXAFS dataset outputs: `*-EXAFS-*_Q_OUTPUT.csv` (`k` vs `χ(k) k²`) and
  `*-EXAFS-*_R_OUTPUT.csv` (`r` vs Fourier-transform components)
- Bragg profiles: `*_bragg.csv`
- R-value logs: `*.log`
- Structure files: `*.rmc6f`, `Frac*.txt`

Most RMCProfile CSV parsers expect first-row labels followed by numeric rows. RMCProfile EXAFS
dataset Q-output files may include a descriptive title row before the column header;
`read_exafs_csv` handles that layout.

## Legacy CLI Scripts

```bash
pip install numpy matplotlib scipy seaborn
python src/RMC_plot.py --dir data [--save --no-show]
python src/RMC_KDE.py [--el Mn]
python src/RMC_3D.py            # needs mayavi
```

`RMC_KDE.py` and `RMC_3D.py` expect `Frac*.txt` plus `.rmc6f` in the working directory.

## Tests

```bash
source .venv/bin/activate
MPLCONFIGDIR=/tmp/rmc_toolkits_matplotlib python -m unittest discover -s tests

# Frontend unit tests (vitest — AI assistant module)
cd web_app/frontend && npm test
```

## Status

The reusable Python package, hosted dashboard, optional Flask API, RMCProfile plot handling,
server-side/browser KDE paths, symmetry panel, AI assistant, and Three.js viewer are in place.
Current priorities: richer run summaries, export/report workflows, and refactoring the legacy
scripts into thin wrappers around `rmc_toolkits`. See [docs/ROADMAP.md](docs/ROADMAP.md) for the
full plan and [docs/CHANGELOG.md](docs/CHANGELOG.md) for history.

## License

Released under the [GNU Affero General Public License v3.0](LICENSE) © 2026 Tsung-Han Yang.

The AGPL is a strong copyleft license: you may use, study, modify, and redistribute this
software, but derivative works must also be released under the AGPLv3. Notably, if you run a
modified version as a **network service**, you must offer its complete source code to the users
of that service (AGPL §13). If you use rmc-toolkits in published research, please cite it.
