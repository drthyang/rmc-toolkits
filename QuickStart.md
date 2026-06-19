# QuickStart

Use `rmc-toolkits` to inspect RMCProfile run folders with interactive plots, model information,
KDE slices, and a 3D structure view.

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
```

After loading, use:

- `Dashboard` for plots, loaded-file badges, and hide/show chart toggles.
- `KDE / 3D` for model information, KDE slices, and the folded unit-cell view.

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

## 4. If Port 5000 Is Busy

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
- `*_bragg.csv` for Bragg profiles.
- `scale_ft.gr`, `scale_ft.sq`, and `scale_ft_rmc.fq` for basic STOG outputs.
