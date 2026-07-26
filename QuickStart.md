# QuickStart

Use `rmc-toolkits` to inspect RMCProfile run folders with interactive plots, model information,
KDE slices, thermal ellipsoids, displacement-direction maps, and a 3D structure view.

## 1. Open The Hosted Dashboard

Go to:

[https://drthyang.github.io/rmc-toolkits/](https://drthyang.github.io/rmc-toolkits/)

Click `Select Folder` and choose a local RMCProfile run folder. Your files stay on your machine.

## 2. Use A Run Folder

A typical folder contains files like:

```text
sample.rmc6f
sample-01.log
sample_FQ1.csv
sample_SQ1.csv
sample_FT_XFQ1.csv
sample_bragg.csv
Nb-EXAFS-1_Q_OUTPUT.csv
Nb-EXAFS-1_R_OUTPUT.csv
```

After loading, use:

- `Dashboard` for plots, loaded-file badges, and hide/show chart toggles.
- `Atomic Density` for model information, KDE slices, and the folded unit-cell view. Drag the
  highlighted band in the `Slab In Cell` panel to move the slice position directly.
- `PCA Ellipsoid` for per-site thermal ellipsoids: pick a site (from the list or by clicking an atom
  in the unit-cell view) to see its PCA displacement ellipsoid, the 3D KDE isosurface with density
  projected on the box walls, and a non-Gaussianity readout.
- `Displacement Directions` for where atoms move rather than how far: the site's displacement
  directions binned on a hex-tiled sphere, highlighting preferred hop directions and ±u asymmetry
  the ellipsoid can't show.

## 3. Live Data (Auto-Refresh)

`Live Data` watches your selected folder and refreshes charts as files are updated. On the hosted
dashboard it works in Chromium browsers (Chrome, Edge, Arc, Opera). In Safari or Firefox, run the
local Flask app for the same feature:

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r web_app/backend/requirements.txt

cd web_app/frontend
npm install
npm run build

cd ../..
python web_app/backend/app.py
```

Open:

[http://127.0.0.1:5000/](http://127.0.0.1:5000/)

Click `Select Folder`, choose your run folder, then turn on `Live Data`.

## 4. AI Assistant (Optional)

The `AI Assistant` tab lets you chat about the loaded run using an LLM. It runs against a model on
your own machine (private) — [Ollama](https://ollama.com) or [LM Studio](https://lmstudio.ai) — or a
cloud provider (OpenAI, Gemini) with your own API key.

Because the hosted dashboard is served over HTTPS, a **local** model server must allow this page's
origin (`https://drthyang.github.io`) via CORS. For Ollama:

```bash
# start Ollama allowing the hosted app:
OLLAMA_ORIGINS="https://drthyang.github.io" ollama serve

# or, for the macOS menu-bar app, set it once and reopen Ollama (resets on logout):
launchctl setenv OLLAMA_ORIGINS "https://drthyang.github.io"
```

For LM Studio, enable CORS in the server settings. Then open `AI Assistant` → pick the provider →
`Test` → choose a model.

**Safari note:** Safari blocks HTTPS pages from calling `http://localhost`. Use Chrome, Edge, or
Firefox, or run the app locally over http so the page and the model share the `localhost` origin:

```bash
cd web_app/frontend && VITE_STATIC_MODE=true npm run dev
```

**Cloud note:** OpenAI/Gemini need no server, but your (summarized) run data is sent to that provider
— the app labels those options and warns before you use them. API keys are stored in your browser only.

## 5. If Port 5000 Is Busy

Use another port:

```bash
RMC_TOOLKITS_PORT=5050 python web_app/backend/app.py
```

Then open:

[http://127.0.0.1:5050/](http://127.0.0.1:5050/)

## Supported Files

- `*.rmc6f` for model information, KDE, and 3D structure views.
- `*.log` for R-value history.
- `*_FQ1.csv` and `*_SQ1.csv` for S(Q).
- `*_FT_XFQ1.csv` and `*PDF*.csv` for PDF/G(r).
- `*-EXAFS-*_Q_OUTPUT.csv` and `*-EXAFS-*_R_OUTPUT.csv` for RMCProfile EXAFS dataset Q-space and
  R-space outputs.
- `*_bragg.csv` for Bragg profiles.

## Want the math?

Every number on every page — how a density map is smoothed, how the space group is decided, how
S(Q) is put on absolute scale, what the AI assistant sends and where — is written up, formula by
formula and anchored to the function that implements it, in
[docs/ALGORITHMS.md](docs/ALGORITHMS.md). The approximations and limitations are listed there too.
