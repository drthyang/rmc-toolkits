# rmc-toolkits frontend

React/Vite frontend for the RMCProfile Run Monitor dashboard.

The production build can run in two modes:

- Static GitHub Pages mode: users choose a local folder in the browser.
- Local Flask mode: the backend serves the built frontend and enables server-side file browsing,
  SciPy KDE, and Live Data folder monitoring.

## Setup

```bash
npm install
```

Use Node.js compatible with Vite 7, such as Node `20.19+` or `22.12+`.

## Development

For a browser-only preview, start Vite directly:

```bash
npm run dev
```

Open the URL printed by Vite, usually `http://localhost:5173/`. This mode reads a selected local
folder in the browser and does not require Flask.

For the full Flask-backed development path, start the backend from the repository root. Use a
non-default port if needed:

```bash
source .venv/bin/activate
RMC_TOOLKITS_PORT=5050 python web_app/backend/app.py
```

Start Vite:

```bash
VITE_API_BASE_URL=http://localhost:5050 npm run dev
```

Open the URL printed by Vite, usually `http://localhost:5173/`.

## Local Production Build

```bash
npm run build
cd ../..
python web_app/backend/app.py
```

Open `http://127.0.0.1:5000/`, click `Select Folder`, and turn on `Live Data` if charts should
refresh when files in the selected folder change.

## Preview Build

```bash
VITE_API_BASE_URL=http://localhost:5050 npm run build
npm run preview -- --host 127.0.0.1 --port 5174
```

Open `http://127.0.0.1:5174/`.

## Notes

- `Select Folder` uses the backend native folder picker in local Flask mode, so it is intended for
  local desktop sessions.
- GitHub Pages/static mode does not support Live Data monitoring because the browser cannot watch
  local filesystem changes.
- The Dashboard page renders parsed plot data as browser-native SVG.
- The KDE / 3D page uses backend SciPy KDE data and a Three.js structure viewer. In static/offline
  mode it instead computes the KDE in a Web Worker, using a WebGPU compute shader when the browser
  supports it and falling back to a CPU loop otherwise.
- On the KDE / 3D page, the `Slab In Cell` band is draggable: grab it to move the slice position
  (`zCenter`) live without using the slider.
- Working on this codebase? See [`AGENTS.md`](../../AGENTS.md) for architecture and conventions, and
  [`docs/CHANGELOG.md`](../../docs/CHANGELOG.md) for history.
