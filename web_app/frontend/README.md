# rmc-toolkits frontend

React/Vite frontend for the local RMCProfile dashboard.

The frontend expects the Flask backend API to be running. In development, Vite reads
`VITE_API_BASE_URL`; if it is not set, development mode falls back to `http://localhost:5000`.
For production builds served by Flask or Docker, the API base defaults to the same origin.

## Setup

```bash
npm install
```

Use Node.js compatible with Vite 7, such as Node `20.19+` or `22.12+`.

## Development

Start the backend from the repository root:

```bash
source .venv/bin/activate
RMC_TOOLKITS_PORT=5050 python web_app/backend/app.py
```

Start Vite:

```bash
VITE_API_BASE_URL=http://localhost:5050 npm run dev
```

Open the URL printed by Vite, usually `http://localhost:5173/`.

## Preview Build

```bash
VITE_API_BASE_URL=http://localhost:5050 npm run build
npm run preview -- --host 127.0.0.1 --port 5174
```

Open `http://127.0.0.1:5174/`.

## Notes

- The app defaults to the repository `data/` sample folder.
- The `Browse` button uses the backend native folder picker, so it is intended for local desktop
  sessions.
- The Dashboard page renders parsed plot data as browser-native SVG.
- The KDE / 3D page uses backend SciPy KDE data and a Three.js structure viewer.
