# Setup & Reference

Everything beyond [using the hosted app](https://drthyang.github.io/rmc-toolkits/): running the
dashboard locally, self-hosting, the Python package, the backend API, and supported file formats.
For a guided tour of the app itself, start with [QuickStart.md](../QuickStart.md); for the
mathematics behind each page — every formula, default, and approximation, anchored to the code —
see [ALGORITHMS.md](ALGORITHMS.md).

## Repository Layout

| Path | Purpose |
| --- | --- |
| `rmc_toolkits/` | Reusable package: parsing, plots, KDE, PCA ellipsoid, displacement directions. |
| `web_app/backend/app.py` | Flask API server with data-root guarding. |
| `web_app/frontend/` | React + Vite single-page app. |
| `web_app/frontend/public/demo/` | Bundled GaTa₄Se₈ 250 K demo run (the app's **Demo** button). |
| `src/` | Original standalone CLI/desktop scripts. |
| `tests/` | `unittest` suite for the package and backend. |
| `docs/` | Changelog, roadmap, architecture notes. |
| `docs/ALGORITHMS.md` + `docs/algorithms/` | Per-page math reference: every operation, anchored to the code. |

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

## Run Locally

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

demo = "web_app/frontend/public/demo"  # bundled GaTa4Se8 250 K example run

frac_path = write_frac_from_rmc6f(f"{demo}/GTS_250K.rmc6f", overwrite=True)
structure = read_structure(demo)

positions = load_unit_cell_positions(f"{demo}/GTS_250K.rmc6f", element="Ga")
payload = kde_slice(
    positions.positions,
    z_center=0.5 * positions.cell_lengths[2],
    dz=0.08 * positions.cell_lengths[2],
    xlim=(0.0, float(positions.cell_lengths[0])),
    ylim=(0.0, float(positions.cell_lengths[1])),
)

png_bytes = plot_to_png(make_plot(f"{demo}/GTS_250K_FQ1.csv"))
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
python src/RMC_plot.py --dir web_app/frontend/public/demo [--save --no-show]
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
