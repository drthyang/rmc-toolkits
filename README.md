# RMCProfile Workbench

[![Tests](https://github.com/drthyang/rmc-toolkits/actions/workflows/tests.yml/badge.svg)](https://github.com/drthyang/rmc-toolkits/actions/workflows/tests.yml)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](LICENSE)

A **browser-first dashboard** for inspecting **RMCProfile** modeling run folders. Open the hosted
app, select a run directory, and review plots, model information, atomic-density KDE slices, PCA
thermal ellipsoids, displacement-direction maps, symmetry, and 3D structure views without
installing anything.

### ▶️ Open the app — [drthyang.github.io/rmc-toolkits](https://drthyang.github.io/rmc-toolkits/)

1. Visit the link above.
2. Click **Select Folder** and choose your RMCProfile run directory — or press **Demo** for a
   bundled example run.
3. Everything renders in your browser.

🔒 **Your raw run files never leave your device.** They are read and rendered entirely in your
browser and are never uploaded to rmc-toolkits or any project server. (Your browser's picker may
say “Upload”, but nothing is sent anywhere.)[^cloud-llm-privacy]

⚡ **Live monitoring** auto-refreshes charts as new files are written, in Chromium browsers (Chrome,
Edge, Arc, Opera).

📖 **New here?** Start with [QuickStart.md](QuickStart.md). Everything else — local/self-hosted
setup, backend API, file formats — lives in [docs/REFERENCE.md](docs/REFERENCE.md).

## Features

- **Run dashboard** — auto-detects RMCProfile outputs (PDF/G(r), S(Q), Bragg profiles, partials,
  EXAFS Q/R CSVs, R-value logs) and renders interactive charts with hover readouts, drag-to-zoom,
  and PNG/SVG/`.zip` export.
- **Live Data** — charts auto-refresh while your run writes new files: client-side in Chromium
  browsers, or server-side through the optional Flask backend.
- **Atomic Density** — KDE density slices (WebGPU with automatic CPU fallback), a draggable
  slab-in-cell projection, and a Three.js folded unit-cell view of `.rmc6f` structures.
- **PCA Ellipsoid** — per-site thermal ellipsoids from the RMC displacement clouds: anisotropic
  displacement tensor, 3D KDE isosurface with wall projections, and non-Gaussianity readouts.
  Follows Maksim Eremenko's
  [PCA_KDE utilities](https://github.com/MaximEremenko/Utilities/tree/main/RMCProfileUtilities/PCA_KDE)
  (independent reimplementation).
- **Displacement Directions** — the direction-space counterpart to the ellipsoid: displacement
  directions binned in solid angle on a hex-tiled sphere reveal discrete hop directions and ±u
  asymmetry that the U tensor cannot see.
- **Symmetry analysis** — a client-side, FINDSYM-like panel reports the detected space group and
  how it changes with tolerance.
- **AI Assistant (beta)** — chat about the loaded run with a local LLM (Ollama, LM Studio) or an
  opt-in cloud model (OpenAI, Gemini). Only compact run context is sent, never raw
  files.[^cloud-llm-privacy] Setup:
  [`web_app/frontend/src/llm/README.md`](web_app/frontend/src/llm/README.md).
- **Python package (`rmc_toolkits/`)** — the same parsing, plotting, KDE, PCA-ellipsoid, and
  displacement-direction analyses as a reusable library, plus `.rmc6f` conversion helpers.

[^cloud-llm-privacy]: If you opt into a cloud LLM, the compact summarized run context used for
    assistant responses is sent directly to the cloud LLM server you selected. Raw run files are not
    uploaded to rmc-toolkits.

## Screenshots

| Run dashboard | Atomic density (KDE / slab / 3D) |
| --- | --- |
| ![Run dashboard](assets/rmc-toolkits-dashboard-demo.png) | ![Atomic density](assets/rmc-toolkits-kde-demo.png) |

| PCA ellipsoid (thermal ellipsoids) | Displacement directions |
| --- | --- |
| ![PCA ellipsoid](assets/rmc-toolkits-pca-demo.png) | ![Displacement directions](assets/rmc-toolkits-displacement-demo.png) |

### AI Assistant

Ask about the loaded run in plain language. The run's metrics, symmetry, and convergence history
travel with every message, so answers quote the actual numbers — and reasoning models stream their
chain of thought in a collapsible *Thinking* panel. Below, a local model (Ollama) summarizes the
bundled demo run as a table, with LaTeX math such as Rwp and χ² rendered inline.

![AI Assistant summarizing a run as tables](assets/rmc-toolkits-assistant.gif)

## Run It Locally (optional)

The hosted app needs no install. Run the Flask backend when you want server-side file browsing,
`.rmc6f` conversion, reference-grade SciPy KDE, or to self-host on a network:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r web_app/backend/requirements.txt && pip install -e .
(cd web_app/frontend && npm install && npm run build)   # Node 20.19+ or 22.12+
python web_app/backend/app.py                           # http://127.0.0.1:5000/
```

Ports, data roots, dev servers, Docker/GitHub Pages deployment, the backend API, and supported
file patterns are covered in [docs/REFERENCE.md](docs/REFERENCE.md).

## Python Package

```python
from rmc_toolkits import kde_slice, load_unit_cell_positions, make_plot, plot_to_png

demo = "web_app/frontend/public/demo"  # bundled GaTa4Se8 250 K example run

positions = load_unit_cell_positions(f"{demo}/GTS_250K.rmc6f", element="Ga")
density = kde_slice(
    positions.positions,
    z_center=0.5 * positions.cell_lengths[2],
    dz=0.08 * positions.cell_lengths[2],
    xlim=(0.0, float(positions.cell_lengths[0])),
    ylim=(0.0, float(positions.cell_lengths[1])),
)
png_bytes = plot_to_png(make_plot(f"{demo}/GTS_250K_FQ1.csv"))
```

Full usage, parser helpers, and the legacy CLI scripts:
[docs/REFERENCE.md](docs/REFERENCE.md#python-package-usage).

## Documentation

- [QuickStart.md](QuickStart.md) — guided tour of the hosted app, including AI-assistant setup.
- [docs/REFERENCE.md](docs/REFERENCE.md) — repository layout, setup, self-hosting, backend API,
  supported file patterns, package usage, legacy CLI scripts, tests.
- [docs/ROADMAP.md](docs/ROADMAP.md) · [docs/CHANGELOG.md](docs/CHANGELOG.md) — plans and history.
- [AGENTS.md](AGENTS.md) — architecture notes and contributor onboarding.

## License

Released under the [GNU Affero General Public License v3.0](LICENSE) © 2026 Tsung-Han Yang.

The AGPL is a strong copyleft license: you may use, study, modify, and redistribute this
software, but derivative works must also be released under the AGPLv3. Notably, if you run a
modified version as a **network service**, you must offer its complete source code to the users
of that service (AGPL §13). If you use rmc-toolkits in published research, please cite it.

*This project is personal work, developed and maintained in my personal capacity.*
