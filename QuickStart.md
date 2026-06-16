# QuickStart

Run `rmc-toolkits` locally to inspect RMCProfile output files, plots, structures, KDE slices, and
the 3D model viewer.

## 1. Install

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r web_app/backend/requirements.txt

cd web_app/frontend
npm install
```

Use Node.js compatible with Vite 7, such as Node `20.19+` or `22.12+`.

## 2. Start The Backend

From the repository root:

```bash
source .venv/bin/activate
RMC_TOOLKITS_PORT=5050 python web_app/backend/app.py
```

The API runs at `http://localhost:5050`.

## 3. Start The Frontend

In a second terminal:

```bash
cd web_app/frontend
VITE_API_BASE_URL=http://localhost:5050 npm run dev
```

Open the Vite URL, usually `http://localhost:5173/`.

## 4. Read The Sample Dataset

The app opens the included sample folder by default:

```text
data/
```

You should see:

- Dashboard plots for R-value, Bragg, S(Q), xPDF, and PDF partials.
- Model information from `data/GNSe.rmc6f`.
- KDE / 3D controls on the `KDE / 3D` page.

## 5. Load Your Own Dataset

Use a folder containing RMCProfile-style outputs, for example:

```text
your-run/
  sample.rmc6f
  sample.log
  sample_FQ1.csv
  sample_FT_XFQ1.csv
  sample_bragg.csv
```

Either type the folder path in the app's `Data path` field and click `Load`, or start the backend
with a custom data root:

```bash
RMC_TOOLKITS_DATA_ROOT=/absolute/path/to/your/runs RMC_TOOLKITS_PORT=5050 python web_app/backend/app.py
```

Then load a folder path relative to that data root.

## Common File Patterns

- `*.rmc6f` for structure and KDE / 3D views.
- `*.log` for R-value history.
- `*_FQ1.csv` or `*_SQ1.csv` for S(Q).
- `*_FT_XFQ1.csv` or `*PDF*.csv` for real-space PDF/G(r).
- `*_bragg.csv` for Bragg profiles.
- `*partials*.csv` for partial PDFs.

## Quick Check

Backend health:

```bash
curl http://localhost:5050/api/health
```

Expected response includes:

```json
{"status": "ok"}
```
