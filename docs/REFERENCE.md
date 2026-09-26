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
server-side file browsing, `.rmc6f` conversion or reference-grade SciPy KDE on your own machine.
To self-host on a network, use Gunicorn or the Docker image ([Hosting The
Dashboard](#hosting-the-dashboard)), not the development server below.

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

`python web_app/backend/app.py` is the Flask development server. It listens on `127.0.0.1`
(this machine only) with debug mode off, and its first log line states the bind address. The
environment variables that change that:

| Variable | Default | Effect |
|---|---|---|
| `PORT`, else `RMC_TOOLKITS_PORT` | `5000` | Port (1–65535). `PORT` wins when both are set. |
| `RMC_TOOLKITS_HOST` | `127.0.0.1` | Bind address. `0.0.0.0` listens on every interface. |
| `RMC_TOOLKITS_DEBUG` | off | `1`/`true`/`yes`/`on` enables Flask debug mode (auto-reload + interactive debugger). |
| `RMC_TOOLKITS_DATA_ROOT` | repo root | The folder the API may read (see above). |

A malformed value stops the server with an error naming the variable. Anyone who can reach the
server can read every file under the data root through the API, so bind a network address only
on a network you trust; for network hosting prefer Gunicorn or the Docker image (below).
**Never set `RMC_TOOLKITS_DEBUG` on a network-reachable server**: the Werkzeug debugger runs
arbitrary Python for whoever reaches it. It is for local development only.

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
  docker run --rm -p 5000:5000 -v /absolute/path/to/runs:/app/data rmc-toolkits-dashboard
  ```

  The image builds from a clean clone and contains no run data: the local `data/` folder is
  excluded by `.dockerignore`, and the image has an empty `/app/data` for you to mount run folders
  on (writable, so Frac conversion and Auto StoG outputs can be written next to the runs). The
  data root is `/app` (`RMC_TOOLKITS_DATA_ROOT`), so mounted runs appear under `data/` in the file
  browser. For a public deployment (Render, Fly.io, Railway, a VPS), the container honors the
  provider's `PORT` and falls back to `5000`; set `RMC_TOOLKITS_DATA_ROOT` to point elsewhere.

  Without Docker, run the same production server from the repository root (Gunicorn is in
  `web_app/backend/requirements.txt`; Flask debug mode is never on under it):

  ```bash
  RMC_TOOLKITS_DATA_ROOT=/absolute/path/to/runs \
    gunicorn --bind 0.0.0.0:5000 --timeout 120 web_app.backend.app:app
  ```

  Either way the whole data root is readable by anyone who reaches the port: put it behind your
  own authentication or firewall on anything but a trusted network.

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

# With no destination the Frac file is written beside the input (into the tracked demo folder);
# point it at scratch space. read_structure pairs it with GTS_250K.rmc6f by stem.
frac_path = write_frac_from_rmc6f(
    f"{demo}/GTS_250K.rmc6f", "/tmp/Frac_coord_GTS_250K.txt", overwrite=True
)
structure = read_structure(demo, frac_path=frac_path)

positions = load_unit_cell_positions(f"{demo}/GTS_250K.rmc6f", element="Se")
payload = kde_slice(
    positions.positions,
    z_center=0.12 * positions.cell_lengths[2],
    dz=0.08 * positions.cell_lengths[2],
    xlim=(0.0, float(positions.cell_lengths[0])),
    ylim=(0.0, float(positions.cell_lengths[1])),
)

png_bytes = plot_to_png(make_plot(f"{demo}/GTS_250K_FQ1.csv"))
```

Lower-level parser helpers are also exported: `read_rmc_csv`, `read_exafs_csv`, `read_chi`,
`read_chi_log`, `read_atom_indices`, `read_cell_vectors`, `iter_rmc6f_atoms` /
`parse_rmc6f_atoms` / `classify_rmc6f_atom_line` (the `.rmc6f` atom-line grammar shared with the
browser, with its `Rmc6fParseReport`), `find_run_configuration` (which `.rmc6f` of a run folder the
app and the CLIs analyse), `frac_lines_from_rmc6f`, `rwp` / `fit_rwp` / `rwp_columns`. The engines
are exported too: `pca_kde_volume`, `site_ellipsoids`, `site_orientation_histogram`,
`bond_angle_summary` / `bond_angle_summary_from_file`, `autoscale`, `estimate_rho0` and the
`transforms` functions; `rmc_toolkits.__all__` lists them all.

## Backend API

All endpoints are under `/api`. Relative paths resolve under `RMC_TOOLKITS_DATA_ROOT`; absolute
paths are rejected unless inside the configured root or a folder selected via the native picker.

**Parameter validation.** Every numeric query-string or JSON-body parameter is parsed by
`_number()` in `web_app/backend/app.py`: it must be a finite number (text, lists, objects,
booleans, `NaN` and `±Infinity` are rejected), integer parameters must be integral, and each value
must lie in its documented range. A violation is **HTTP 400** with an `error` message naming the
parameter; a missing or blank parameter takes its default. Grid sizes are the exception: they are
clamped to the engine's limits instead of rejected (and `z` is clamped to [0, 1] and echoed as
`center`). The KDE-slice, PCA-KDE and orientation routes also refuse to serialize a result that
came out `NaN`/`Infinity` for finite but extreme values (e.g. a PCA-KDE `extent` or `bwScale` of
`1e300`): that is a 400 too, never a 200 whose body is invalid JSON or whose map is all `null`.
That check runs the strict stdlib encoder itself (`_strict_result_response()`), so it holds
whatever JSON provider the app installs: writing a non-finite value as `null` is right for a masked
*data series* (a gap in a chart), never for a computed density. The scaling routes apply the same
rule to their fitted `a`/`b` and output series (a manual `a = b = 1e308` is a 400, and nothing is
written). A KDE-slice bandwidth so small that the kernel cannot be evaluated (σ below $10^{-10}$
in-plane fractional units, e.g. `bw=1e-200`) is not an error: both runtimes return the finite zero
map with its `kernel`, flagged `subgrid` + `unresolved`. A bandwidth so large that forming
`bw²·C` overflows (e.g. `bw=1e300`) is a 400 naming `bw`. Error statuses: 400 bad parameter or unusable input (including a
plot file that cannot be parsed or decoded, a JSON body that is not an object, and a non-string
`path`/`kind`/`mode`/`formula`/`outDir`… field), 403 path
outside the data roots, 404 missing file/folder, 409 output exists (`/api/scaling/run` without
`force`) or source file still being written (see below), 500 unexpected failure.

**The static-mode workers apply the same rules** (`workers/requestGuards.js`: `requestNumber()`
mirrors `_number()`, `assertFiniteResult()` mirrors `_strict_result_response()` with its exact
message). The PCA worker requires an integer orientation `frequency` and an integer `smoothing`
in [0, 64], refuses a KDE or orientation result holding `NaN`/`Infinity`, and rejects an unknown
request `kind`; the Auto StoG worker requires a finite, non-zero `a` and a finite `b` in manual
mode and refuses a non-finite fit (as `_require_finite_scaling()`); the structure worker reads
`maxPoints` as an integer clamped to [100, 10⁶], as `/api/structure` does. Every worker answers a
null or malformed message with an error, so the page's request always settles.

**Caching and freshness.** The KDE-slice, PCA (`sites`/`kde`/`orientation`), triplets and scaling
routes keep small in-process LRU caches of parsed files (`_FileCache` in `app.py`). Every cache key
is the file's signature `(st_mtime_ns, st_ctime_ns, st_size, st_ino)` from `_file_signature()`,
never `st_mtime` alone, so a file rewritten within the same whole-second mtime (sshfs/SFTP mounts,
`scp -p`, rsync from a coarse filesystem) is re-parsed. The signature is taken again after each
parse: a parse of a file that changed while it was being read is never cached — it is re-read once
under the new signature, and if the file changes again (a writer still busy) the request fails
with **409** and a "changed while it was being read" message instead of returning a torn result.

The table lists all 15 routes. A `dir` parameter names a run folder (the backend picks its
`.rmc6f` by output-stem match, else the first alphabetically); `path` names one file. Booleans in
query strings are true for `1`/`true`/`yes`. Defaults are in parentheses.

| Method & path | Parameters (default; accepted range) | Returns |
| --- | --- | --- |
| `GET /api/health` | — | `status`, active `dataRoot`. |
| `GET /api/files` | `dir` (`.`) | Sub-folders plus files matching the supported patterns: `name`, `path`, `type`, `plotKind`, `modified` (`st_mtime`, seconds), `size` (bytes). Live Data polls this. |
| `POST /api/dialog/folder` | JSON `dir` (`.`): where the native picker opens | `path`, `name` of the chosen folder, which becomes an allowed data root; 400 when cancelled. |
| `GET /api/plot` | `path` | One supported file rendered as a PNG. |
| `GET /api/plot/metadata` | `path` | `kind`, `title`, `metrics` (`rwp`, `final_chi_r`). |
| `GET /api/plot/data` | `path` | Metadata plus `xLabel`, `yLabel` and `series` (`label`, `x`, `y`) for the SVG plots; 400 for an unsupported file. |
| `POST /api/convert/frac` | JSON `path` (a `.rmc6f`), `outputPath` (next to the source), `overwrite` (false) | `path`, `name` of the written `Frac_coord_<stem>.txt`; 409 if it exists and `overwrite` is false. |
| `GET /api/structure` | `dir`; `maxPoints` (1 000 000; integer, clamped to [100, 1 000 000]) | Atoms folded into one unit cell (sampled per site above `maxPoints`), `totalAtoms`, `elementCounts`, `atomIndices`, `supercell`, `latticeVectors`, and the `.rmc6f` move counters `moves` (`generated`, `tried`, `accepted`, `accumulatedTimeS`) when the header has them. |
| `GET /api/kde/slice` | `dir`; `element` (all; case-insensitive; an element the file lacks is a 400 naming those it has, and a file with no parseable atom is a 400); `orientation` `a`/`b`/`c` (`c`), any other value = custom normal `nx`, `ny`, `nz` (0, 0, 1; finite, not all zero) with an optional in-plane frame `ux`, `uy`, `uz`, `vx`, `vy`, `vz` (all six or none; each axis non-zero and orthogonal to the normal and to the other, within 10⁻⁶; ignored for the presets; the page always sends it, and without it the library's own frame is used, which for e.g. (1 1 0) is the page's rotated by 90°); `z` (0.5; finite, clamped to [0, 1]); `dz` (0.08; 0 < dz ≤ 1); `bw` (0.03; > 0, SciPy `gaussian_kde` scalar factor); `grid` (120; integer, clamped to [16, 400]); `levels` (8; integer in [0, 64]); `log` (false) | Density grid, `extent`, contour polylines, `slabCount`, `fitCount`, slab/plane geometry, `kernel` (`covariance`, `sigmaMinor`, `sigmaMajor`; in-plane fractional units; `H = bw²·C`, `C` the covariance of the slab's atoms), `message` (why nothing was drawn — the same texts as the browser worker, plus `engine` when the server's SciPy cannot evaluate the kernel) and `warnings` (`subgrid`, `unresolved`). **`z` and `dz` are fractions of the unit cube's projection range along the slice normal** (equal to cell-edge fractions only for the `a`/`b`/`c` presets) and are echoed as given; `depth`/`depthThickness` are in depth-projection units. The slab, bandwidth and grid stay in fractional coordinates, with no conversion to Å. The KDE fit uses at most 6000 rows — slab atoms plus their periodic images within the kernel margin, so at larger `bw` a few thousand atoms already reach the cap; above it the server and the browser sum different subsamples (maps differ by several percent of the peak). |
| `GET /api/pca/sites` | `dir`; `probability` (0.5; 0 < p < 1) | Per-site displacement tensor and thermal ellipsoid table (`sites`: per site also `nonGaussianity` — Mardia's multivariate excess kurtosis — `excessKurtosis` with `axisResolved`, `zeroSpread`, `elementCounts` and `mixed`), `referenceNumbers`, `elements`, `totalAtoms`, `latticeVectors`, `supercell`, `parseWarning` (atom lines the `.rmc6f` grammar skipped — non-finite coordinates, unparsed lines, a header count mismatch — or `null`). |
| `GET /api/pca/kde` | `dir`; `referenceNumber` (integer) or `element` (pools that element's sites; neither pools every atom); `bw` (`scott`; `scott`, `silverman` or a number > 0); `bwScale` (1.0; > 0); `grid` (48; integer, clamped to [8, 128]); `extent` (3.0; > 0, box half-width in kernel-broadened σ); `cubicBox` (false: only sizes the display box `boxHalfWidths`; the volume is always sampled on the per-axis box `halfWidths`); `probability` (0.5; 0 < p < 1); `projections` (true) | PCA frame, ellipsoid, and the separable 3D KDE volume (+ three wall projections), with the captured `mass`, `nonGaussianity`, `axisResolved` and `boxHalfWidths`. A zero-spread site (λ₁ < 10⁻⁸ Å², e.g. an average configuration) is a 400. A volume that captures less than 10⁻⁶ of the density (every kernel between the grid nodes: a tiny `bw`/`bwScale` or a huge `extent`) is a 400, not an all-zero volume. |
| `GET /api/pca/orientation` | `dir`; `referenceNumber` or `element` as above; `frequency` (auto: the recommended value; integer in [1, 64]); `weight` (`count`; `count`, `amplitude`, `amplitude2`); `minAmplitude` (0 Å; ≥ 0); `minAmplitudeQuantile` (0; in [0, 1)); `smoothing` (0; integer in [0, 64]); `frame` (`cartesian`; or `pca`); `geometry` (true: include cell polygons) | Goldberg-cell histogram of displacement directions: `enhancement`, `zScore` (local, uncorrected), peak direction with `peakSignificance` (exact Poisson tail, Šidák over all cells), `mapSignificance` / `mapPValue` (Pearson's X² against its exact-moment isotropic reference; `null` below 0.1 expected coincident pairs) with `mapNullSd`, `mapNullSkewness`, `mapExpectedPairs`, `antipodalAsymmetry` against its exact null (`…Null`, `…NullSd`, `…Z`, `…Significant` at z > 3), `orientationAnisotropySignificance` (Bingham); the legacy `significance` (an RMS local z) and `peakZScore` are kept but are not significances. `parseWarning` as in `/api/pca/sites`. |
| `GET /api/triplets` | `dir`; `end1`, `apex`, `end2` (elements, case-insensitive; `apex` is the central atom); `r12Min` + `r12Max` (required, Å, inclusive; `r12Max` ≤ 15); `r23Min` + `r23Max` (optional pair, default = the 1-2 window; `r23Max` ≤ 15); `binWidth` (1.0°; ≥ 0.05) | Bond-angle histogram, bond-length statistics (`count` = B-centred bond vectors, `uniqueBonds` = physical bonds), coordination histogram (payload of `triplets.bond_angle_summary`), `parseWarning`. A spec that would form more than `APP_MAX_ANGLES` = 5×10⁷ angles is a 400 naming the count, before any angle is formed. The 15 Å, 0.05° and 5×10⁷ caps are API limits; the engine itself is unrestricted. |
| `POST /api/scaling/preview` | JSON `path` (a `stog.inp` or an S(Q) data file); `kind` (`auto`; `inp`, `data`); `inspect` (false: only parse the source). Numeric overrides, each a finite number: `qmin`, `qmax`, `rho0`, `bAvgSq`, `bSqAvg`, `massDensity`, `rCutoff` (≥ 0; data mode 1.0), `rmax` (data mode 50), `nr` (data mode 5000; integer), `r0`, `rFitMin`, `rFitMax`. Text: `formula`, `c1Mode` (`sweep`), `amplitude` (`density`). Booleans (JSON `true`/`false`; the strings `1`/`true`/`yes` are true): `lorch` (data mode false), `lowQCorrection` (true), `robust` (true), `despike` (false), `useSigma` (true). `mode` (`auto`; `manual` takes `a`, a finite non-zero number, and `b` (0; finite), falling back to the `stog.inp` values); low-r enforcement controls `enforce` (boolean), `enforceCutoff` (finite number), `peakWindow` (`[rmin, rmax]`, two finite numbers) — defaults in [auto-stog.md](algorithms/auto-stog.md) | Fitted `a`, `b`, convergence history, diagnostics, provenance, plot guides, `warnings` (the coefficient warnings the CLI prints, e.g. a `formula`'s ⟨b²⟩ left unused because its ⟨b⟩² disagrees with `bAvgSq`; `[]` when none), and the S(Q)/G_K(r)/D(r) series (+ low-r–enforced versions). A `stog.inp` supplies every value not overridden; data mode requires `qmin`, `qmax`, a density (`rho0`, `massDensity` + `formula`, or a `NUMBER_DENSITY ::` header) and ⟨b⟩² (`bAvgSq` or `formula`). |
| `POST /api/scaling/run` | As `preview`, plus `outDir` (`<source folder>/autoscale`), `outStem`, `force` (false) | Writes the classic stog file family + `stog_provenance.json` (the `rmc-autoscale` CLI writer); returns `a`, `b`, `outputs`, `outDir`, diagnostics, `warnings` (as `preview`). 409 when an output exists and `force` is false; 400 (even with `force`) when an output would overwrite the input data file or the stog.inp, when two outputs name one file (names compared case-insensitively), when an output path is a directory, or when an output folder cannot be created (a path component is a file). These checks run before any computation, and the family is written through temporary files renamed into place only after every write succeeded, so a failed run leaves no partial family. |

## Supported File Patterns

- Real-space PDF/G(r): `*_FT_XFQ1.csv`, `*PDF*.csv`
- Reciprocal-space S(Q): `*_FQ1.csv`, `*_SQ1.csv`
- RMCProfile EXAFS dataset outputs: `*-EXAFS-*_Q_OUTPUT.csv` (`k` vs `χ(k) k²`) and
  `*-EXAFS-*_R_OUTPUT.csv` (`r` vs Fourier-transform components)
- Bragg profiles: `*_bragg.csv`
- χ² logs: `*.log` (the run's `<stem>-NN.log` files; the last column is charted)
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

# Frontend unit tests (vitest: engine ports and their Python goldens, components, AI assistant)
cd web_app/frontend && npm test
```
