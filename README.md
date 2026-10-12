# RMCProfile Workbench

[![Tests](https://github.com/drthyang/rmc-toolkits/actions/workflows/tests.yml/badge.svg)](https://github.com/drthyang/rmc-toolkits/actions/workflows/tests.yml)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](LICENSE)

A browser-first workbench for inspecting **RMCProfile** modeling runs:
**[drthyang.github.io/rmc-toolkits](https://drthyang.github.io/rmc-toolkits/)**.

## What it does

- **Run dashboard:** interactive charts of the fits (G(r), F(Q)/S(Q), Bragg, partials, EXAFS) and
  the χ² history, with model information, the detected space group, and PNG/SVG/ZIP export.
- **Four analysis pages:** Atomic Density (KDE slices and a 3D folded cell), Bond Geometry
  (bond-angle distributions), PCA Ellipsoid (per-site thermal ellipsoids, after Maksim Eremenko's
  [PCA_KDE utilities](https://github.com/MaximEremenko/Utilities/tree/main/RMCProfileUtilities/PCA_KDE))
  and Displacement Directions (where atoms move, with significance tests).
- **Private:** run files are read in your browser and never uploaded; a cloud AI model, if you
  choose one, receives summaries only.
- **Live Data:** charts refresh while a run writes new files (Chromium browsers, or the local
  Flask backend).
- **Python package and CLIs:** the same engines as a library (`rmc_toolkits`), bond angles from the
  command line (`rmc-triplets`), automatic S(Q) scaling (`rmc-autoscale`), and an optional Flask
  backend.

## Quick start

Open [the app](https://drthyang.github.io/rmc-toolkits/), then **Select Folder** for your run
directory, or press **Demo** for a bundled example run (GaTa₄Se₈, 250 K). To run it locally:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r web_app/backend/requirements.txt && pip install -e .
(cd web_app/frontend && npm install && npm run build)   # Node 20.19+ or 22.12+
python web_app/backend/app.py                           # http://127.0.0.1:5000/
```

## AI Copilot

A panel beside the pages that answers questions about the open run by running the Workbench's
analyses itself, with deterministic checks on each result and on its answer. It asks when a choice
changes the result, and points at the cards it explains. It works with a local model (Ollama,
LM Studio) or, with your key, a cloud one. See the
[Copilot guide](web_app/frontend/src/llm/README.md) and [QuickStart](QuickStart.md#4-ai-copilot-optional).

## Status

Beta: version 0.7.0, with unreleased changes on main ([changelog](docs/CHANGELOG.md)). The AI
Copilot is experimental: trust its step cards' checks over its own verdict. An answer the checks
contradict is marked **Not verified**, with the check that contradicts it.

## Documentation

- [QuickStart.md](QuickStart.md): a tour of the app, Live Data and Copilot setup, screenshots.
- [docs/REFERENCE.md](docs/REFERENCE.md): local and self-hosted setup, the backend API, file
  patterns, the Python package and command-line tools, tests.
- [docs/ALGORITHMS.md](docs/ALGORITHMS.md): the math of every page, anchored to the code.
- [AI Copilot guide](web_app/frontend/src/llm/README.md): how it works, its tools, privacy, setup.
- [docs/ROADMAP.md](docs/ROADMAP.md) · [docs/CHANGELOG.md](docs/CHANGELOG.md) · [AGENTS.md](AGENTS.md)
  (contributor guide).

## License

[AGPL-3.0](LICENSE) © 2026 Tsung-Han Yang ([what that means, and citing](docs/REFERENCE.md#license-and-citation)).

*This project is personal work, developed and maintained in my personal capacity.*
