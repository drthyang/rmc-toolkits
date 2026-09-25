# Run Dashboard — algorithm reference

> Part of the [Algorithms and Math Reference](../ALGORITHMS.md). Every step is anchored to the source; if this document and the code disagree, **the code wins**.

How a run folder becomes plots and numbers: file detection and parsing, the fit-quality metrics, the plot renderer and export path, and the model summary with the client-side symmetry finder.

## Contents

- [Run Dashboard — file detection, parsing, and fit-quality metrics](#run-dashboard--file-detection-parsing-and-fit-quality-metrics)
  - [What this section covers — parsing and metrics](#what-this-section-covers--parsing-and-metrics)
  - [Step 1 — Enumerate the run folder](#step-1--enumerate-the-run-folder)
  - [Step 2 — Classify each file: `detect_plot_kind()`](#step-2--classify-each-file-detect_plot_kind)
  - [Step 3 — Pick the structure file and the run-control `.dat`](#step-3--pick-the-structure-file-and-the-run-control-dat)
  - [Step 4 — Parse the file into series](#step-4--parse-the-file-into-series)
  - [Step 5 — Compute the numbers](#step-5--compute-the-numbers)
  - [Step 6 — Assemble titles, axis labels, and series for drawing](#step-6--assemble-titles-axis-labels-and-series-for-drawing)
  - [Step 7 — Build the model summary from the `.rmc6f`](#step-7--build-the-model-summary-from-the-rmc6f)
  - [Step 8 — Order, group, and render](#step-8--order-group-and-render)
  - [Step 9 — Live Data: change detection and refresh](#step-9--live-data-change-detection-and-refresh)
  - [Parameters and defaults — parsing and metrics](#parameters-and-defaults--parsing-and-metrics)
  - [Caveats / what this is not](#caveats--what-this-is-not)
- [Plot rendering, interaction, and figure export](#plot-rendering-interaction-and-figure-export)
  - [What this section covers — plot rendering and export](#what-this-section-covers--plot-rendering-and-export)
  - [Step 1 — Where the series arrays come from](#step-1--where-the-series-arrays-come-from)
  - [Step 2 — Series ordering, roles, colours and marker style](#step-2--series-ordering-roles-colours-and-marker-style)
  - [Step 3 — Legend visibility filter (and the label-identity assumption)](#step-3--legend-visibility-filter-and-the-label-identity-assumption)
  - [Step 4 — Axis domain selection (autoscale and padding)](#step-4--axis-domain-selection-autoscale-and-padding)
  - [Step 5 — Viewport geometry and the data → pixel affine map](#step-5--viewport-geometry-and-the-data--pixel-affine-map)
  - [Step 6 — Tick generation ("nice numbers"), tick labels and tick geometry](#step-6--tick-generation-nice-numbers-tick-labels-and-tick-geometry)
  - [Step 7 — Path construction: clipping, NaN handling, markers, no decimation](#step-7--path-construction-clipping-nan-handling-markers-no-decimation)
  - [Step 8 — Pointer → user-space coordinates](#step-8--pointer--user-space-coordinates)
  - [Step 9 — Nearest-point hover search](#step-9--nearest-point-hover-search)
  - [Step 10 — Drag-to-zoom: rectangle → data range inversion](#step-10--drag-to-zoom-rectangle--data-range-inversion)
  - [Step 11 — Wheel zoom (x only)](#step-11--wheel-zoom-x-only)
  - [Step 12 — Reset, double-click, and view-state lifetime](#step-12--reset-double-click-and-view-state-lifetime)
  - [Step 13 — Export path: SVG serialisation and PNG rasterisation](#step-13--export-path-svg-serialisation-and-png-rasterisation)
  - [Step 14 — "Save all figures" → one ZIP](#step-14--save-all-figures--one-zip)
  - [Step 15 — Canvas figures: device pixel ratio, the 3× export, and colormaps](#step-15--canvas-figures-device-pixel-ratio-the-3-export-and-colormaps)
  - [Step 16 — The matplotlib path (`plots.py`) and the metrics](#step-16--the-matplotlib-path-plotspy-and-the-metrics)
  - [Parameters and defaults — plot rendering and export](#parameters-and-defaults--plot-rendering-and-export)
  - [Screen vs. export: what actually changes](#screen-vs-export-what-actually-changes)
  - [Caveats / what this is not](#caveats--what-this-is-not-1)
- [Model summary and the Detected SG symmetry finder](#model-summary-and-the-detected-sg-symmetry-finder)
  - [What this section covers — model summary and symmetry](#what-this-section-covers--model-summary-and-symmetry)
  - [Part A — Model information](#part-a--model-information)
  - [Part B — The Detected SG symmetry finder](#part-b--the-detected-sg-symmetry-finder)
  - [Parameters and defaults — model summary and symmetry](#parameters-and-defaults--model-summary-and-symmetry)
  - [Caveats / what this is not](#caveats--what-this-is-not-2)

---

## Run Dashboard — file detection, parsing, and fit-quality metrics

### What this section covers — parsing and metrics

The **Run Dashboard** turns an RMCProfile run folder into a page of charts without asking the user
to name anything. It enumerates the folder, classifies each file by *filename pattern only*, parses
the recognized ones into (x, y…) series, computes a small set of numbers per file (an R-factor, a
PDF index, a final χ value), and draws every series as an interactive SVG chart. It also shows a
collapsible **R-value** strip built from the run's `.log` files, a model summary from the `.rmc6f`,
and — when Live Data is on — refreshes all of that on a fixed poll while the run is still going.

Two independent implementations exist and both are documented here:

| Runtime | File listing | Parsing + metrics | Rendering |
|---|---|---|---|
| **Flask mode** | `GET /api/files` ([web_app/backend/app.py](../../web_app/backend/app.py) → `list_files()`) | Python, [rmc_toolkits/parsers.py](../../rmc_toolkits/parsers.py) + [rmc_toolkits/plots.py](../../rmc_toolkits/plots.py), served by `/api/plot/metadata` and `/api/plot/data` | browser SVG ([InteractivePlot.jsx](../../web_app/frontend/src/components/InteractivePlot.jsx)) |
| **Static mode** (GitHub Pages) | browser directory pick, or the bundled demo | JavaScript, [web_app/frontend/src/browserData.js](../../web_app/frontend/src/browserData.js) | same SVG renderer |

**Which implementation runs is decided at load time by a four-branch rule**, not by a single
environment variable (`browserData.js` → `isStaticMode()`):

1. `VITE_STATIC_MODE === 'true'` → static;
2. `VITE_STATIC_MODE === 'false'` → Flask;
3. otherwise `import.meta.env.DEV && !VITE_API_BASE_URL` → static;
4. otherwise `window.location.hostname.endsWith('github.io')`.

Branch 3 matters in practice: a plain `vite dev` session with no API base URL silently uses the
**JavaScript** parsers even though a Flask server may be running. Everywhere the two
implementations disagree — and they do, in the ways catalogued below — the deployment decides which
answer the scientist sees.

The matplotlib figures in `plots.py` (`make_plot()` → `plot_to_png()`) are *not* what the dashboard
draws. They are reachable only through `GET /api/plot`, whose only consumer,
`PlotViewer.jsx`, is not mounted anywhere in the current app. In Flask mode `make_plot()` is still
called by `/api/plot/metadata` and `/api/plot/data` — a full matplotlib `Figure` is constructed and
immediately closed — purely to obtain the `metrics` dictionary and the title. That detour has two
observable consequences (a hard column-count precondition, Step 6, and a repeated-parse cost,
Step 9).

#### Notation and units

$Q$ = momentum transfer (Å⁻¹); $r$ = radial distance (Å); $k$ = photoelectron
wavenumber (Å⁻¹); ToF = neutron time of flight (µs); $\chi_r$ = RMCProfile's goodness metric as parsed from the last `.log` column
(dimensionless; described as `chi^2` elsewhere in the codebase — see 5d); $N$ = number of data rows in a file.

**Reproducing the examples.** Every filename-level example below is taken from the run bundled with
the repository at [web_app/frontend/public/demo/](../../web_app/frontend/public/demo/) — a GaTa₄Se₈
250 K run, 10 tracked files: `GTS_250K.rmc6f`, `GTS_250K.dat`, `GTS_250K-00.log`,
`GTS_250K-01.log`, `GTS_250K-02.log`, `GTS_250K_FQ1.csv`, `GTS_250K_FQ1partials.csv`,
`GTS_250K_FT_XFQ1.csv`, `GTS_250K_PDFpartials.csv`, `GTS_250K_XFQ1.csv`. Where an example needs a
file the demo does not contain (notably `scale_ft_rmc.fq`), that is stated explicitly: `data/` is
gitignored ([.gitignore](../../.gitignore) line 2) and no run folder under it ships with the repo.

---

### Step 1 — Enumerate the run folder

**Inputs:** one of three things — a folder path (Flask), a browser directory selection (static), or
the bundled demo run (either mode).

**Flask mode.** `list_files()` resolves the requested directory through `_resolve_inside_root()`
(rejects anything outside `RMC_TOOLKITS_DATA_ROOT`, default = repo root, plus any folder the user
picked with the native picker), lists immediate subdirectories whose name does not start with `.`,
then globs the directory for these patterns (`SUPPORTED_PATTERNS`, app.py):

```
*.csv  *.log  *.rmc6f  Frac*.txt  scale_ft.*  stog_input.dat  *.inp  *.sq  *.dat
```

The glob is **non-recursive** (one directory level). Hits are collected into a dict **keyed by the
resolved `Path`**, so a file matched by several patterns (e.g. `scale_ft.sq` matches both
`scale_ft.*` and `*.sq`) appears exactly once. Each hit becomes
`{name, path, type, plotKind, modified, size}` via `_file_payload()`, where `modified` is
`st_mtime` in **seconds since the epoch (float)** and `size` is bytes. The returned list is sorted
by `(item['type'] != 'directory', item['name'].lower())` — **directories first, then
case-insensitive name**.

`_file_payload()` swallows `OSError` from `stat()` and returns `modified = size = None`. Those
nulls become empty strings in the change signature (Step 9), so **a file whose `stat()` fails is
invisible to change detection** — it will never trigger a Live Data refresh.

**Static mode.** `makeRunFromEntries()` in `browserData.js` filters `{path, file}` pairs through
`isSupportedFile()`, which accepts a name if it ends with `.csv`, `.log`, `.rmc6f`, `.inp`, or
`/\.(gr|sq|fq)$/i`, **or** starts with `Frac`, **or** is one of
`{scale_ft.gr, scale_ft.sq, scale_ft_rmc.fq, stog_input.dat}`. `modified` here is
`File.lastModified`, **milliseconds** since the epoch.

**Demo run (both modes).** The header's **Demo** button calls `loadDemoRun()`, which `fetch`es a
hard-coded list of the 10 files above from `${BASE_URL}demo/`, wraps each blob as
`new File([blob], name, { lastModified: Date.now() })` under the path prefix `Demo/`, and feeds the
result through the *same* `makeRunFromEntries()` pipeline as a folder pick. Two consequences worth
knowing: the demo is always parsed by the **JavaScript** readers even when a Flask backend is
present, and because `lastModified` is the moment of loading, **the demo's file timestamps carry no
information about the run** and every demo load produces a different change signature (Step 9).

> **The two file lists differ.** Flask accepts every `*.dat` and every `*.sq` but **not** `*.gr` or
> `*.fq` (`scale_ft.*` matches only names with the literal prefix `scale_ft.`, so a file named
> `scale_ft_rmc.fq` is *not* listed by Flask — an example from a local, gitignored run folder, not
> from the bundled demo). The browser accepts every `.gr`/`.sq`/`.fq` but no generic `.dat`.
> A folder therefore yields slightly different file sets in the two runtimes.

Directory recursion also differs: the File System Access path walks subdirectories recursively
(`collectHandleEntries()`), the `<input webkitdirectory>` path receives whatever the browser
enumerates, and Flask never recurses.

**Code:** `app.py` → `list_files()`, `_file_payload()`; `browserData.js` → `isSupportedFile()`,
`makeRunFromEntries()`, `buildLocalRun()`, `buildLocalRunFromHandle()`, `collectHandleEntries()`,
`loadDemoRun()`.

---

### Step 2 — Classify each file: `detect_plot_kind()`

**Inputs:** the **basename only**. No file content, no header sniffing, no size check.

The Python `detect_plot_kind()` ([plots.py](../../rmc_toolkits/plots.py)) and the JavaScript
`detectPlotKind()` ([browserData.js](../../web_app/frontend/src/browserData.js)) test the same patterns
**in this order** and return on the first match:

| # | Pattern (regex or test) | Kind | Meaning |
|---|---|---|---|
| 1 | `-EXAFS-.+_Q_OUTPUT\.csv$` | `exafs_q` | EXAFS k-space output |
| 2 | `-EXAFS-.+_R_OUTPUT\.csv$` | `exafs_r` | EXAFS R-space (Fourier-transformed) output |
| 3 | `_FT_XFQ\d+\.csv$` | `xpdf` | X-ray PDF obtained by Fourier-transforming F(Q) |
| 4 | name contains `PDF` **and** ends `.csv`, and contains `PDFpartials` | `pdf_partials` | per-pair partial PDFs |
| 5 | name contains `PDF` **and** ends `.csv` (otherwise) | `npdf` | neutron PDF |
| 6 | ends `_FQ1.csv` | `xray_sq` | titled "S(Q) (x-ray)" |
| 7 | ends `_SQ1.csv` | `neutron_sq` | titled "S(Q) (neutron)" |
| 8 | `_bragg(?:_.+)?\.csv$` | `bragg` | Bragg profile |
| 9 | `-\d{2,}\.log$` (an **inline** regex in both files, not `R_VALUE_LOG_RE`) | `r_value` | RMCProfile run log (χ history) |
| 10 | **Python:** name ∈ `{scale_ft.gr, scale_ft.sq, scale_ft_rmc.fq}`<br>**JS:** `\.(gr\|sq\|fq)$` (case-insensitive) | `stog` | STOG/preprocessing data file |
| — | anything else | `null` | not plotted |

Ordering is load-bearing where the patterns genuinely overlap: **rules 4–5 fire before rules 6–7**,
so a name containing both `PDF` and `_FQ1.csv` is classified as a PDF, never `xray_sq` — as
`pdf_partials` if it also contains `PDFpartials` (rule 4), otherwise as `npdf` (rule 5). Both are
one `if "PDF" in name and name.endswith(".csv")` branch in the source, split into two rows here
because they return different kinds. The rule-3-before-rules-4–5 ordering only bites for an xPDF
file whose own stem contains the literal substring `PDF` (e.g. `myPDFrun_FT_XFQ1.csv`); the demo's
`GTS_250K_FT_XFQ1.csv` contains no `PDF`
substring and would be `xpdf` under either ordering.

**Known consequences of pattern-only classification** (all reproducible against the bundled demo
folder, [web_app/frontend/public/demo/](../../web_app/frontend/public/demo/)):

* `GTS_250K_XFQ1.csv` (the x-ray F(Q)) matches **nothing** — `_FQ1.csv` requires an underscore
  immediately before `FQ1`, and `XFQ1.csv` supplies an `X`. That file is indexed but never charted.
* `GTS_250K_FQ1partials.csv` matches nothing (it does not end in `_FQ1.csv`), so partial F(Q) files
  are silently skipped. Partial *PDF* files are not, because rules 4–5 key on the substring `PDF`
  anywhere in the name (`GTS_250K_PDFpartials.csv` → rule 4).
* Rule 6 asserts that `_FQ1.csv` is **x-ray** and rule 7 that `_SQ1.csv` is **neutron**. This is a
  naming convention, not something read from the data; the title chip will say "S(Q) (x-ray)" for
  any `*_FQ1.csv` regardless of radiation.
* Rule 9 needs **two or more** digits: `run-01.log` is an R-value log, `run-1.log` and
  `derivative.log` are not (`tests/test_plots.py::test_detect_plot_kind_for_supported_outputs`
  pins `GNSe.log → None`, `run-info.log → None`, `GNSe-123.log → r_value`).
* Rule 8 requires the character after `bragg` to be `_` or `.`; the same test pins
  `GNSe_braggish.csv → None`.

**Python/JS agreement:** rules 1–9 are character-for-character equivalent. Rule 10 **differs**: the
browser treats *any* `.gr`/`.sq`/`.fq` as `stog`, the Python only three fixed filenames. This is
deliberate (`browserData.js` comment: "runs often use descriptive data-file names") and pinned by
`web_app/frontend/src/__tests__/browserData.test.js` ("detects any .gr / .sq / .fq file as a STOG
plot").

**`stog` files never reach the Dashboard.** `Dashboard.jsx` filters with
`isDashboardPlotFile = (file) => file.plotKind && file.plotKind !== 'stog'`, applied *before* any
metadata or plot-data call, so STOG files are detected and counted in the run diagnostics but are
neither charted nor listed on this page.

---

### Step 3 — Pick the structure file and the run-control `.dat`

**Inputs:** the classified file list.

A run folder usually holds several `.rmc6f` files (e.g. the working configuration and an
`…AVERAGE.rmc6f`). The chosen one is the `.rmc6f` whose **stem matches an output file's stem**,
searched in this priority order. The stem is always **capture group 1** of the matching pattern; the
last column names the kind of output file the pattern matches, not the value extracted:

| Priority | Output-name pattern | Matches (stem = group 1) |
|---|---|---|
| 0 | `^(.+)-\d{2,}\.log$` | run log (`GTS_250K-01.log` → `GTS_250K`) |
| 1 | `^(.+)-EXAFS-.+_[QR]_OUTPUT\.csv$` | EXAFS |
| 1 | `^(.+)_FT_XFQ\d+\.csv$` | xPDF |
| 1 | `^(.+)_[FS]Q\d+\.csv$` | S(Q) |
| 1 | `^(.+)_bragg(?:_.+)?\.csv$` | Bragg |
| 1 | `^(.+)_PDF(?:partials\|\d+)?\.csv$` | PDF |
| 2 | `^Frac_coord_(.+)\.txt$` | fractional-coordinate export |

Candidates are sorted by (priority, lowercase filename) and the first whose stem has a matching
`.rmc6f` wins; otherwise the **first** `.rmc6f` in the list is used. Python sorts the directory
alphabetically before this scan, so "first" is alphabetical; the browser uses the enumeration order
of the directory pick, so the two can disagree on the fallback. The browser additionally keys the
stem map by `dirname/stem`, so matching is per-subfolder; the Python only ever looks at one
directory.

**Code:** `app.py` → `_run_stem_from_output_name()`, `_find_rmc6f()`; `browserData.js` →
`runStemFromOutputName()`, `chooseStructureFile()`.

**Run-control file (static mode only).** `chooseSettingsEntry()` computes
`wanted = structureFile.path.replace(/\.rmc6f$/, '.dat')` and returns the raw entry whose `path`
equals it exactly. **The exact path/stem match is what distinguishes the run-control file** from the
auxiliary `.dat` files a run folder carries (`chi2.dat`, `optimization.dat`, `adv_opt.dat`, …). The
lookup runs over the **raw** entries rather than the filtered list only because `isSupportedFile()`
would already have discarded it — generic `.dat` is not a supported browser file type
(`SUPPORTED_NAMES` admits only `stog_input.dat`).

`pairFitTypes()` then reads the **first 131072 bytes (128 KiB)** of up to **6** candidate `.dat`
files (`runControlCandidates()` — the stem-matched entry first, then any other `.dat`) and runs
`parseRunSettings()` on that head. Its control flow has three properties that are easy to miss:

* It **returns after the first candidate that yields a non-empty fit-type map.** A wrong-but-parseable
  `.dat` earlier in the list therefore wins outright, and the correct one is never read.
* Every read/parse failure is swallowed by a bare `catch {}`. A truncated, binary or unreadable
  `.dat` produces no warning; labels silently fall back to the extension.
* The map from `fitTypeByFilename()` is keyed by `basename(dataset.file).toLowerCase()`. A
  `FILENAME` declared with a directory path still matches on basename alone, and two datasets whose
  files share a basename collide (last one wins). `FIT_TYPE` wins over `DATA_TYPE`; at most **8**
  `*_DATA` blocks are kept.
* `pairFitTypes` then iterates the **whole supported-file list**, not just plot files, so any
  supported file whose lowercase basename appears in the map gets a `fitType` property — including
  a `.rmc6f` or `Frac*.txt` if a dataset block happened to name it.

**What `fitType` actually does on this page: nothing.** It is consumed only inside the
`kind === 'stog'` branches of `plotMetadataFromFile()` and `plotDataFromText()`, and
`isDashboardPlotFile` drops every `stog` file before either function is called. The `D(r)`-instead-of-
`G(r)` label override is therefore **unreachable from the Run Dashboard**; it lives in
`browserData.js` and is exercised only by `browserData.test.js`. **No conversion is applied to the
numbers** in any case — the same column is plotted either way.

The run-control `.dat`'s one live effect on this page is on the **AI-assistant run context**:
`localRun.settingsFile` is re-read and passed to `parseRunSettings()` in a separate effect that is
gated on `wantAssistantData`, i.e. it stays idle until the user opens the AI Assistant page.

---

### Step 4 — Parse the file into series

**Inputs:** file text. **Outputs:** `{labels, data}` where `data` is column-major (transposed).

#### 4a. Standard RMC CSV — `read_rmc_csv()` / `readRmcCsv()`

Used for `xpdf`, `npdf`, `pdf_partials`, `xray_sq`, `neutron_sq`, `bragg`.

* Line 1 is the header; labels = comma-split, whitespace-stripped, **empties preserved**.
* Every subsequent line is comma-split, stripped, and **empty fields are dropped** (this is what
  lets RMCProfile's trailing-comma data rows, e.g. `GTS_250K_PDFpartials.csv`, parse against a
  header that has no trailing comma).
* A row whose surviving field count ≠ the label count raises
  `"<path> line <n> has <k> values; expected <m>"`.
* Result is `np.asarray(rows, float).T` — column 0 is x, columns 1… are the y series.

**Blank-line handling differs, and it moves the header.** The JavaScript reader pre-filters *all*
blank/whitespace-only lines (`text.split(/\r?\n/).filter((line) => line.trim())`) and then takes
element 0 as the header; Python takes `lines[0]` of the raw file. A file with a leading blank line
therefore gets its real header in the browser and a **blank header (label count 1)** in Flask. The
same pre-filter shifts the error message's line number: Python enumerates true file lines
(`enumerate(lines[1:], start=2)`) while JavaScript reports `index + 2` over the blank-filtered
array. The message text is identical; **the line number is the true file line in Python and the
index among non-blank lines in JavaScript**, so the two disagree on any file containing blank lines.

**Python/JS discrepancy on non-numeric tokens:** Python calls `float(value)` and lets a `ValueError`
propagate, so a stray non-numeric row **fails the whole file**. JavaScript calls `Number(value)`,
which yields `NaN` silently; the SVG renderer then drops non-finite points at draw time
(`Number.isFinite` guard in `InteractivePlot.jsx` → `seriesShapes`). The same file can therefore
plot in static mode and error in Flask mode — and, as 5a shows, produce a *misleading metric* rather
than an error.

**For this reader only,** both runtimes accept the literal token `NaN` as a number and neither masks
it before computing metrics (see 5a). That is not true of the other three readers: see 4b, 4c
and 4d.

#### 4b. EXAFS CSV — `read_exafs_csv()` / `readExafsCsv()`

RMCProfile's `_Q_OUTPUT` files carry a descriptive title row *above* the column header; `_R_OUTPUT`
files do not. Both readers therefore scan for the **first fully numeric comma-separated row**, take
the line immediately above it as the header, and skip every later non-numeric line. A file whose
first line is already numeric (`data_start == 0`) is rejected — there would be no header. Pinned by
`tests/test_parsers.py::test_read_exafs_csv_skips_q_output_title_row` and
`…accepts_r_output_header`.

**"Fully numeric" is defined differently in the two runtimes.** Python's `_numeric_csv_values()`
uses `float()`, which accepts `NaN` and `Inf`; JavaScript's `numericCsvValues()` requires
`parsed.every(Number.isFinite)`, which rejects them. A row of `NaN`s is therefore numeric to Python
and non-numeric to JavaScript, so the same EXAFS file can pick a **different line as the header**
and keep a **different number of rows** in the two runtimes.

Typical columns: Q-space `k, calculated, experiment`; R-space
`r, Re_Calc, Im_Calc, Mod_Calc, Re_Ex, Im_Ex, Mod_Ex`. Every column after the first is drawn as its
own series.

#### 4c. STOG data file — `read_stog()` / `readStog()`

Both skip **exactly the first two lines** (the classic STOG layout is `count`, then a title line),
whitespace-split the rest, transpose, and plot only columns 0 and 1. Beyond that they are not the
same function:

* **Python `read_stog()`** ([parsers.py](../../rmc_toolkits/parsers.py)) has **no** `try`/`except`. Every
  remaining non-empty line goes through `[float(value) for value in parts]`, so a single
  non-numeric token — a trailing text line, a units row — raises `ValueError` and **fails the whole
  file**. `np.asarray(rows).T` is then built with no column-count check, so ragged rows raise as
  well. `float('NaN')` succeeds, so `NaN` tokens are **retained**.
* **JavaScript `readStog()`** delegates to `parseNumberRows(lines, 2)`, which keeps a row only when
  `parsed.every(Number.isFinite)` and **silently drops** every other row — text rows and rows
  containing `NaN` alike. `transpose()` then uses the first surviving row's length as the column
  template, so a ragged row yields `undefined` entries rather than an error.

(The tolerant "keep the rows that parse" behaviour that *does* exist in Python lives in
`read_stog_xy()`, a different function the dashboard path never calls.)

#### 4d. R-value log — `read_chi()` / `readChi()`

**Skips exactly the first two lines** of each `.log` (RMCProfile writes a column-name row and a
`WEIGHT PARAMETERS` row), whitespace-splits the rest, and for every line with ≥ 2 tokens takes:

$$\chi_Q \leftarrow \texttt{parts[-2]} , \qquad \chi_r \leftarrow \texttt{parts[-1]}$$

i.e. **the last and second-to-last whitespace-separated columns, by position**. Neither
implementation reads the header to identify which dataset those columns belong to. In the bundled
demo run the log columns are
`Time, moves_acc, moves_gen, F(Q)_1, Curvature ×6, X_ray_(R)1`, so `parts[-1]` is the X-ray R
column and `parts[-2]` is a (identically zero) curvature-constraint column. **For a run with a
different dataset/constraint ordering, a different quantity is plotted.** This is the single most
fragile heuristic on the page.

The skipped second line is where the run's own per-dataset weights live — in the demo run it reads
`h/m/s/.th  WEIGHT PARAMETERS  0.100E+01 …`. Those weights are **discarded here and never applied
anywhere in the app** (see 5a).

The JavaScript reads only the last column (it never builds `chi_q`) and keeps a value only
`if (Number.isFinite(value))`. The Python builds both inside a `try`/`except ValueError`, which lets
`NaN`/`Inf` through — so a log containing such a token yields a **longer series in Flask than in the
browser**. The Python also has an ordering quirk worth knowing about: it appends to `chi_q` *before*
parsing `parts[-1]`, so a line whose last token is non-numeric leaves `chi_q` one element longer
than `chi_r`. Only `chi_r` is plotted, so the dashboard is unaffected.

**Multiple logs are concatenated — but by different code in each mode.** Python
`related_r_value_logs()` re-scans the log's parent directory for every file matching
`^(.+)-(\d{2,})\.log$` (`R_VALUE_LOG_RE`) with the **same stem**, and `sort_r_value_logs()` orders
them by (lowercase stem, integer sequence, lowercase name) — so `run-01.log, run-02.log, run-10.log`,
*numerically*, not lexicographically
(`tests/test_parsers.py::test_related_r_value_logs_use_numeric_suffix_order`). This runs **server-side
inside `/api/plot/data`**, on every request.

The browser's `combineRValueFiles()` in `Dashboard.jsx` concatenates the already-parsed y arrays of
the r-value files in `comparePlotFiles()` order (same stem/sequence rule, via a fourth copy of the
pattern, `rValueLogParts`). It **short-circuits and returns `rValueFiles[0]` unchanged** when any of
these hold:

1. there is only one r-value file;
2. no file carries `sourceFile`, `plotData` or `parseError` — which is exactly the **Flask** case,
   where files come from `/api/files` and carry none of those;
3. any browser-backed log is still being parsed (`sourceFile && !plotData && !parseError`).

So the concatenation shown here happens **client-side only in static mode**. Condition 3 means that
during a Live Data re-parse the strip transiently shows only the first log's curve. And the
"hiding a log removes its segment from the combined curve" behaviour (Step 8) is **static-mode
only**: in Flask mode `related_r_value_logs()` always reads every same-stem log on disk regardless
of what the user hid.

---

### Step 5 — Compute the numbers

Three metric keys exist across the whole page: `rwp`, `pdf_index`, `final_chi_r`.

#### 5a. The R-factor reported as "Rwp"

**Inputs:** the parsed CSV columns. **Applies to:** `xpdf`, `npdf`, `xray_sq`, `neutron_sq`,
`bragg` — and only when the file has **≥ 3 columns**. Never for `pdf_partials`, `exafs_*`, `stog`,
or `r_value`.

`rwp()` in [parsers.py](../../rmc_toolkits/parsers.py) (and the near-identical `rwp()` in
`browserData.js`):

```python
paired = isfinite(observed) & isfinite(fitted)   # NaN on either side drops the point
if not paired.any(): return None                 # nothing to compare
denom = sum(observed[paired]**2)
if denom == 0: return None                       # no scale to normalize against
residual = fitted[paired] - observed[paired]
return sqrt(sum(residual**2) / denom)
```

Called as `rwp(series.data[0], series.data[1], series.data[2])` — the x column is passed and
**ignored** (`void x` in the JS). So the reported number is

$$R \;=\; \sqrt{\dfrac{\sum_{i=1}^{N}\bigl(y^{(3)}_i - y^{(2)}_i\bigr)^2}{\sum_{i=1}^{N}\bigl(y^{(2)}_i\bigr)^2}}$$

where $y^{(2)}$ is the file's **second** column and $y^{(3)}$ its **third**.

Four things must be stated plainly:

1. **It is not weighted.** There is no $w_i = 1/\sigma_i^2$ anywhere. Despite the "Rwp" chip label
   and the parameter names, this is the unweighted R-factor $R = \|\Delta\|_2 / \|y\|_2$. The run's
   own per-dataset weights *are* physically present, on line 2 of the `.log`
   (`WEIGHT PARAMETERS 0.100E+01 …`), and are skipped by `read_chi()`/`readChi()` and used nowhere.
2. **The denominator is the calculated curve, not the data.** RMCProfile writes these CSVs in the
   order `(x, calculated, experimental)` — verified in the demo run:
   `Q, F(Q)_RMC, F(Q)_Expt` and `r(A), X_ray-calc, X_ray_exp_renorm`. The parameter named
   `observed` therefore receives the **RMC-calculated** column and the one named `fitted` receives
   the **experiment**. The numerator is symmetric so the residual is right, but the normalization
   is $\sum y_\mathrm{calc}^2$, not the conventional $\sum y_\mathrm{obs}^2$. The two agree only
   insofar as $\sum y_\mathrm{calc}^2 \approx \sum y_\mathrm{obs}^2$ — close for a good fit,
   not identical, and **not** the standard crystallographic $R_\mathrm{wp}$.
3. **Only columns 2 and 3 ever enter it.** For a file with four or more numeric columns (multi-bank
   outputs, partial-inclusive exports) every column beyond the third is **drawn but contributes
   nothing to the number**, and which two curves get compared is purely positional.
4. **It is strictly per-dataset.** One value per file, over **every row in the file** — no Q or r
   window, no exclusion region, no point weighting by Δx (so a non-uniform grid is summed as a
   plain point sum). There is no combined/global R anywhere in the app, and the app never
   reproduces RMCProfile's own weighted total χ².

`tests/test_parsers.py::test_rwp_uses_observed_series_as_denominator` pins the formula:
`rwp([0,1,2], [2,4,4], [1,5,6]) == sqrt(6/36)`.

**Python vs JavaScript:** the formula and inputs are the same, and on clean numeric data the two
agree to floating-point round-off (relative ~1e-15, since NumPy sums pairwise and the JS `reduce`
sums sequentially) rather than bit-for-bit. **On degenerate input they now agree as well:**

| Input | Python `rwp()` | JavaScript `rwp()` |
|---|---|---|
| `denom == 0` (flat-zero 2nd column) | `None` | `null` |
| `NaN` throughout the **2nd** (calculated) column | `None` | `null` |
| `NaN` throughout the **3rd** column | `None` | `null` |
| `NaN` in *some* rows of either column | value over the finite rows | value over the finite rows |

This matters because it is reachable in static mode: `readRmcCsv` turns unparseable cells into
`NaN` instead of raising (4a). Until 2026-07-25 the guards were `if denom == 0` (parsers.py) versus
`if (!denom) return 0` (browserData.js), and `!NaN === true` — so a file whose calculated column
was entirely non-numeric rendered a chart with a confident-looking **"Rwp 0.000"** chip, while
Python returned `NaN` for the same file and `json.dumps` emitted it as the bare token `NaN`, which
is not valid JSON, leaving Flask-mode callers with no usable metrics object at all. Both now return
the unavailable sentinel, which serializes as JSON `null` and reaches both runtimes intact. The
degenerate cases are pinned by `tests/test_parsers.py::test_rwp_is_none_when_*` /
`::test_rwp_uses_only_points_finite_in_both_series` and the matching `describe('Rwp metric')` block
in `browserData.test.js`.

Displayed as a chip, `Number(metrics.rwp).toPrecision(4)` — **4 significant digits** — or `—` when
the metric is the null sentinel (`Dashboard.jsx` → `renderRwpChip()`, shared by the plot card and
the R-value panel).

#### 5b. `pdf_index`

For `npdf` only: `pdf_index()` extracts the integer from `PDF(\d+)\.csv$`, or **0** when the name
has no number. It is metadata (which of several neutron PDF banks this is), not a fit quality
number, and it is not currently surfaced in the dashboard UI.

#### 5c. `final_chi_r`

For `r_value`: the **last** element of the raw `chi_r` array — **not** log-transformed —
after concatenating all related logs. `tests/test_plots.py::test_log_plot_combines_related_logs_in_numeric_order`
pins that three logs `run-01/02/10` yield the value from `run-10.log`. In static mode
`combineRValueFiles()` takes `final_chi_r` from the **last parsed** log file, matching the Python;
in Flask mode `combineRValueFiles()` short-circuits (4d) and the value is whatever the server
computed from its own concatenation.

#### 5d. The logarithm — where and why

The R-value series is plotted as the **natural logarithm** of the parsed column:

$$y_i \;=\; \ln\!\bigl(\max(v_i,\, 10^{-12})\bigr)$$

where $v_i$ is the raw value from the last whitespace column of log row $i$.

* `browserData.js` → `plotDataFromText()` (kind `r_value`): `Math.log(Math.max(value, 1e-12))`.
* `app.py` → `plot_data()` (kind `r_value`): `math.log(max(value, 1e-12))`.
* `plots.py` → `_chi_plot()` (matplotlib PNG path): `np.log(chi_r)` — **no clamp**. A zero or
  negative entry produces `-inf`/`nan` and a NumPy warning here, where the other two paths floor at
  $\ln(10^{-12}) \approx -27.63$.

The reason for the log is dynamic range: the metric falls by orders of magnitude over a run, and
additive shifts in ln-space are *relative* changes in the metric, which is exactly what the
convergence watchdog thresholds on.

> **What is actually being plotted — flagged.** The axis label is `log(χ)` in all three paths
> (`plots.py` `r"log($\chi$)"`, `app.py` `"log(χ)"`, `browserData.js` `'log(χ)'`), while
> `src/llm/context/runContext.js` describes the same array as
> `"ln of chi^2 goodness metric (natural log; lower is better)"` and names the metric key
> `final_chi_r`. **Neither label is verified by the code**, which reads `parts[-1]` positionally and
> never looks at the header. In the bundled demo run that column is headed `X_ray_(R)1` — an
> R-factor for the X-ray dataset, not a χ² and not a χ. The only statements this document can stand
> behind are: the transform is a **natural** log with a $10^{-12}$ floor, and the quantity is
> *whatever the last whitespace column of the run's `.log` happens to hold*. Treat both `log(χ)` and
> "ln(χ²)" as conventions the code does not check.

The x-axis is the **row index** (`0, 1, 2, …`) labelled `"Time steps"`. The `.log`'s actual first
column is a wall-clock stamp (`hhmmss.sss`) and is discarded, so one "time step" is one log row,
i.e. one RMCProfile print/save period — for the demo run, `SAVE_PERIOD :: 10.00 MINUTES`, not a
Monte-Carlo step.

**Downstream use (dashboard badge) — opt-in, and usually invisible.** `WatchdogBadge` classifies the
same ln-series via `src/llm/watchdog/heuristics.js`, but `useWatchdog()` returns
`{status: 'off', …}` unless **both** `settings.watchdogEnabled` is true (an AI-assistant setting)
**and** the combined series has ≥ 2 points, and `WatchdogBadge` renders `null` on `status === 'off'`.
By default nothing is drawn.

When it does run, over the ln-series $y_0…y_{n-1}$:

* window length $L = \max\bigl(\min(n,5),\, \lceil 0.2n \rceil\bigr)$ points (the tail of the series);
* `recentSlope()` is the ordinary least-squares slope of that tail against its **local index**
  $i = 0…L-1$,
  $$m \;=\; \frac{\sum_i (i-\bar i)(y_i-\bar y)}{\sum_i (i-\bar i)^2},$$
  returning `0` when the denominator is 0. **Units: ln(metric) per log row**, not per unit time.
* $\texttt{windowDelta} = m\,(L - 1)$;
* `classifyConvergence()` returns `'unknown'` for `values.length < 2` (badge label "Watching"),
  else `diverging` if `windowDelta > 0.02`, `improving` if `< −0.02`, else `stalled` when the total
  drop $y_0 - y_{n-1} < 0.1$, else `converged`.

Because the values are in ln units, 0.02 ≈ a 2 % change in the metric and 0.1 ≈ a 10 % improvement,
independent of magnitude. These are pure heuristics; the LLM only narrates them, and the narration
gate (`useWatchdog.js`) is stricter than the heuristic itself. A call is made only when a provider
`baseUrl` **and** `model` are configured, no call is already in flight, and either

* the classified `status` **changed** since the last LLM call, **or**
* **both** at least `watchdogIntervalMin` minutes (default `5`) have elapsed since the last call
  **and** `significantChange()` fires — ≥ 200 new steps, or a last-value move of at least
  `ln(1.02)`.

`significantChange()` alone is therefore *not* sufficient: outside a status flip it is ANDed with the
interval. The attempt timestamp is recorded in the `finally` block, so a failing provider is retried
on the same throttle schedule rather than on every poll.

---

### Step 6 — Assemble titles, axis labels, and series for drawing

**Inputs:** kind + parsed columns. **Outputs:** the object the SVG renderer consumes:
`{kind, title, metrics, xLabel, yLabel, series[]}`.

**A Python-only precondition applies first.** `_series_plot()` raises
`ValueError("<path> needs at least two numeric columns")` whenever the parsed file has fewer than
2 columns. Because `/api/plot/metadata` *and* `/api/plot/data` both call `make_plot()` before doing
anything else, a single-column file makes **both endpoints return HTTP 500 in Flask mode** — no
title, no metrics, no chart. The JavaScript has no such guard: `plotDataFromText()` returns
`series: []` and the card renders an **empty chart**. Note that this check lives inside the
matplotlib builder, so it applies even though the dashboard never uses the PNG.

| Kind | Title | x label | y label | Rwp? |
|---|---|---|---|---|
| `exafs_q` | `EXAFS Q-space` | `k (Å⁻¹)` | `χ(k) k²` | no |
| `exafs_r` | `EXAFS R-space` | `r (Å)` | `FT[χ(k) k²]` | no |
| `xpdf` | `xPDF` | `r (Å)` | `G(r)` | yes |
| `npdf` | last `_`-segment of the stem (e.g. `PDF1`) | `r (Å)` | `G(r)` | yes |
| `pdf_partials` | last `_`-segment of the stem | `r (Å)` | `G(r)` | no |
| `xray_sq` | `S(Q) (x-ray)` | `Q (Å⁻¹)` | `S(Q)` | yes |
| `neutron_sq` | `S(Q) (neutron)` | `Q (Å⁻¹)` | `S(Q)` | yes |
| `bragg` | `BRAGG` | `ToF (µs)` **or** `Q (Å⁻¹)` | `Intensity` | yes |
| `r_value` | `R-value` | `Time steps` | `log(χ)` (see 5d) | no |
| `stog` † | see below | `r (Å)` if `.gr`, else `Q (Å⁻¹)` | see below | no |

† **The `stog` row is unreachable from the Run Dashboard** (`isDashboardPlotFile` drops it, Step 2)
and the three implementations do not agree on it, so it is recorded here only for completeness:

* browser **metadata** (`plotMetadataFromFile`): `title = fitType` if known, else `G(r)`/`S(Q)` by
  extension;
* browser **plot data** (`plotDataFromText`): `title = file.name`; the fit-function label lands in
  `yLabel` only;
* **Flask**: no `fitType` concept at all — `/api/plot/data` returns
  `yLabel = "G(r)" if name.endswith(".gr") else "S(Q)"`, and the title comes from `make_plot()` →
  `_stog_plot()`, which sets `title = path.name`.

**Bragg axis selection** — `bragg_is_tof(header)` (plots.py) / `braggAxis(header)`
(browserData.js) is a case-insensitive regex on the **first column's header text**:

```
/tof|flight|time/  →  ToF (µs)      anything else  →  Q (Å⁻¹)
```

The raw x values are used **as written** — no ToF↔Q conversion, no unit rescaling. RMCProfile's
time-of-flight Bragg CSVs label the column `Flight time (us)`; the app takes microseconds as given.
`tests/test_plots.py::BraggAxisTests` pins `"Flight time (us)"`, `"TOF,ms"`, `"Time"` → ToF, and
`"Q or theta"`, `"2-theta, deg"`, `""`, `None` → not-ToF. Note the last two: a `"TOF,ms"` header is
still labelled **µs**, and a 2θ-in-degrees Bragg dataset is labelled **Q (Å⁻¹)**. Both are
mislabels the code makes knowingly.

**`_clean_axis_label()` / `cleanAxisLabel()` never run on this page.** Both sit in the final `else`
of a branch chain that has already handled every kind `detect_plot_kind()` can return, so the branch
is dead for the dashboard. For the record, the substitution table is
`Q` → `Q (Å^{-1})`, `r`/`R` → `r (Å)`, `(A^-1)`/`(A^{-1})` → `(Å^{-1})`, `(A)` → `(Å)` — a cleanup
step the axis labels **do not** get.

**Series construction.** Every column after the first becomes one series with the file's own header
text as its label; an empty label becomes `Series <n>` where **n is the 1-based index among the
non-first labels** (so column 2 of the file is `Series 1`). A label that is declared but has no
matching column is skipped (`if idx < len(series.data)` / `.filter((series) => series.y)`).

**Series identity is the label string.** The legend uses `key={series.label}` and the hidden-series
set stores labels, so two columns sharing a header — common in partials exports — collapse into one
legend button and are hidden/shown together. The renderer also honours per-series `role: 'guide'`,
`defaultHidden` and a chart-level `initialYDomain`; **none of these are ever produced by the
dashboard's parsers** (they exist for the Auto StoG page) and can be ignored when reading this
section.

**Transformations applied to the numbers before drawing: none.** No offset, no normalization, no
scaling, no smoothing, no rebinning, no subsampling, and **no difference/residual curve is ever
computed or drawn** — `grep` for `residual` in the dashboard path finds only the two lines inside
`rwp()`. What you see is the file's columns, in file order, plotted against the file's first column.
The one exception is the R-value panel's natural log (Step 5d).

The SVG renderer does apply presentation-only adjustments. These are documented in full in the
**Plot rendering, interaction, and figure export** section; the ones that change what a dashboard
user sees or reads off a chart are:

* **`/exp/i` markers are conditional.** Series whose label matches `/exp/i` are treated as measured
  data and drawn first as hollow markers (radius 2.8) with the calculated curve on top **only when
  the chart contains both at least one `/exp/i` series and at least one non-`/exp/i` series**
  (`paired = experimental.length > 0 && calculated.length > 0`). Otherwise file order is kept and
  every series draws as a plain polyline.
* **Autoscale in y follows the x window.** The y domain is computed only from `y` values whose
  matching `x` lies inside the *current* x domain (falling back to the full y range when that filter
  is empty), so zooming in x silently rescales y. Both domains are then padded by 5 % of their
  range; if min = max the domain is first widened by ±1, and if no finite values remain it is
  `[0, 1]`. Ticks are snapped to 1/2/5 × 10ⁿ, 7 x × 6 y for grid cards and 11 x × 4 y for the wide
  R-value strip.
* **Clipping is done in data space with no edge interpolation.** `seriesShapes` drops every point
  whose x lies outside the current x domain before building the path, so a zoomed line ends at the
  last in-window data point rather than at the axis. Non-finite points are dropped too.
* **Hover is a nearest-neighbour lookup, not an interpolation.** `nearestHover()` scans every point
  of every visible non-guide series and picks $\arg\min_i |x_i - x_\mathrm{cursor}|$ **independently
  per series**, so the tooltip can report points at different x values for different series. The
  crosshair is anchored to the *first* visible series' chosen point, not to the cursor. The scan is
  O(N) per series per pointer move.
* **The tooltip numbers are rounded.** `formatNumber()` returns `value.toExponential(2)` when
  $|v| \ge 10^6$ or $0 < |v| < 0.01$, `Math.round(value).toLocaleString()` when $|v| \ge 1000$ (so
  any value above 1000 is shown as an **integer**), and `value.toPrecision(4)` otherwise. Tick
  labels are compacted separately (`12.5k`, `3.2M`; exponential below 1e-4 or above 1e9). Hover
  readouts are therefore 4-significant-digit renderings of the parsed value, not the parsed value.

**Code:** `plots.py` → `make_plot()`, `_series_plot()`, `_chi_plot()`, `_stog_plot()`;
`app.py` → `plot_data()`, `plot_metadata()`;
`browserData.js` → `plotDataFromText()`, `plotMetadataFromFile()`;
`InteractivePlot.jsx` → `orderedSeries`, `domains`, `seriesShapes`, `nearestHover`, `formatNumber`.

> **PNG-path discrepancies (matplotlib only).** The figures built by `_series_plot()` differ from
> the interactive charts in two places in their **axis labels** (the styling differences — `lw=1.0,
> alpha=0.65` strokes, and a legend and title baked into the figure — are covered in Step 16b and in
> item 1 of "Screen vs. export"):
> * **y label** — `_series_plot()`'s default `ylabel="data"` is used for
>   `xpdf`/`npdf`/`pdf_partials`/`xray_sq`/`neutron_sq`/`bragg`, instead of `G(r)`/`S(Q)`/`Intensity`.
> * **x label** — only `xray_sq` and `neutron_sq` pass `labels[0]`, the file's **raw first-column
>   header**, instead of `Q (Å⁻¹)`. `xpdf`, `npdf` and `pdf_partials` pass the literal
>   `r ($\mathrm{\AA}$)` and `bragg` passes the same ToF/Q label the interactive path derives from
>   `bragg_is_tof()`, so those match (in LaTeX form).
>
> `_stog_plot()` additionally draws a **dashed black horizontal reference line** at $y = 1$ (or
> $y = 0$ when the name ends `.fq`) spanning `data[0][0]` to `data[0][-1]`, colours its single
> series red, and uses a taller figure (6.75 × 4.725 in). That line is the **only** reference/guide
> line anywhere in this page's code, and it is on the PNG path the UI never calls.

---

### Step 7 — Build the model summary from the `.rmc6f`

**Inputs:** the `.rmc6f` chosen in Step 3. **Outputs:** the object `ModelSummary.jsx` renders
(unconditionally, above the charts) and the point cloud the Structure pages reuse. This runs
**concurrently** with Steps 4–6.

**Static mode.** `Dashboard.jsx` spawns `workers/localStructureWorker.js`, which reads the file text
and calls `browserData.js` → `structureFromRmc6f(file, maxPoints)` with **`maxPoints = 100`**.

**Flask mode.** `loadServerDashboard()` requests `GET /api/structure?dir=…&maxPoints=100`. The query
argument is clamped to $[100,\ \texttt{MAX\_STRUCTURE\_POINTS} = 1{,}000{,}000]$.

The two produce **different subsamples and different fields**:

* **Flask — site-stratified quota.** `_sample_atoms_by_site()` returns all atoms unchanged when
  $N_\mathrm{atoms} \le$ `max_points`. Otherwise it groups atoms by `reference_number`, sets
  $\mathrm{quota} = \max(1, \lfloor \text{max\_points}/n_\mathrm{sites}\rfloor)$, and from each group
  (ascending reference number) takes `group[::max(1, len(group)//quota)][:quota]`, finally
  truncating to `max_points`. With the dashboard's `maxPoints = 100` and, say, 52 reference sites,
  that is **1 atom per site**. The `sampleStride` reported back to the UI is a *different* number,
  $\max\bigl(1,\ \lfloor N_\mathrm{atoms}/\text{max\_points}\rfloor\bigr)$ — not the stride actually
  applied to any group. The `max(1, …)` is defensive — this branch is reached only when
  $N_\mathrm{atoms} >$ `max_points`, so the floor is already $\ge 1$; the
  $N_\mathrm{atoms} \le$ `max_points` early return reports `1` outright.
* **Static — uniform global stride.** `structureFromRmc6f()` uses
  $s = \lceil N_\mathrm{atoms}/\mathrm{maxPoints}\rceil$, keeps every $s$-th atom and truncates with
  `.slice(0, maxPoints)`. No site stratification, so rare sites can vanish from the sample entirely.

Unit-cell folding also differs in form (both are equivalent up to the removal of an integer):
Python computes `reduced = coords − cell_indices/supercell` then `(reduced * supercell) % 1.0`;
the browser computes `((coords * supercell) % 1 + 1) % 1` directly, which also works for old files
that carry no per-atom cell index.

**Counts are over all atoms; points are the subsample.** `elementCounts`, `totalAtoms` and
`atomIndices` are accumulated over **every** atom in the file in both runtimes; only `points` is
subsampled (`sampledAtoms`). The element rows and totals a user reads off the card are therefore
exact even though the plotted cloud is a sample.

**Fields only static mode produces.** The Flask `/api/structure` response contains **no `basis` and
no `moves`**:

* `basis` — one representative site per `reference_number`, built from per-axis **circular means** of
  the within-cell fraction, with a per-site rms displacement `dispA` in Å derived from the circular
  resultant. The full derivation (circular mean, the $\sigma = \sqrt{-2\ln R}/2\pi$ wrapped-normal
  relation, the resultant floor of `1e-6`, and the cell-edge normalization
  $a_i = |\text{lattice row } i| / \max(\mathrm{supercell}_i, 1)$) is documented in the
  **Model summary and the Detected SG symmetry finder** section; it is not repeated here.
* `moves` — run-history counters scraped from the `.rmc6f` header by `readMovesMetadata()` with four
  regexes (`Number of moves generated/tried/accepted:`, `Accumulated time (s)…:`). It reads only the
  text before the `Atoms:` marker, or the **first 4000 characters** if that marker is not found in
  the leading text. These feed the AI-assistant run context only.

**Element names are normalized differently.** Python `iter_rmc6f_atoms()` stores
`parts[1].capitalize()`; the browser's `parseAtomLine()` keeps the raw token. A `.rmc6f` written with
`SE` or `se` yields one merged `Se` row in Flask mode and one or two raw `SE`/`se` rows in static
mode — **the element table of the model summary differs between runtimes for such a file**. The
browser parser also accepts a 5–6 field coordinates-only line (`referenceNumber = null`, excluded
from the site basis); the Python requires ≥ 9 fields and drops those lines entirely.

**What `ModelSummary.jsx` computes on top of this** — box lengths $|\text{lattice row}|$,
conventional cell lengths $|\text{lattice row}_i| / \max(\mathrm{supercell}_i, 1)$, the three cell
angles from `acos` of the normalized dot products (with a `1e-12` denominator floor and a $[-1,1]$
clamp), the display rounding (3 fractional digits for lengths, 1 for angles, locale-formatted), and
the client-side symmetry finder at a **default tolerance of 0.2 Å** with a tolerance ladder capped
at **1.0 Å** — is documented in the **Model summary and the Detected SG symmetry finder** section.
It is listed here because the card is part of this page and its numbers come from the data path
above.

**Code:** `Dashboard.jsx` → the `localRun` effect (worker spawn, `maxPoints: 100`) and
`loadServerDashboard()`; `workers/localStructureWorker.js`; `browserData.js` →
`structureFromRmc6f()`, `readCellVectors()`, `readMovesMetadata()`; `rmc6f.js` → `parseAtomLine()`;
`app.py` → `structure()`, `_sample_atoms_by_site()`; `parsers.py` → `iter_rmc6f_atoms()`,
`read_cell_vectors()`, `read_atom_indices()`; `ModelSummary.jsx`.

---

### Step 8 — Order, group, and render

`Dashboard.jsx` sorts the **chartable** plot files (`allPlotFiles = files.filter(isDashboardPlotFile)`)
with `comparePlotFiles()`:

1. by kind, using the fixed order
   `['r_value', 'bragg', 'xray_sq', 'neutron_sq', 'exafs_q', 'exafs_r', 'xpdf', 'npdf', 'pdf_partials']`
   (a kind absent from the list sorts *first*, because `indexOf` returns −1 — reachable only if a
   new kind is added without updating this array);
2. for two r-value logs, by lowercase stem then integer sequence (`rValueLogParts`, a browser-side
   copy of `R_VALUE_LOG_RE`);
3. otherwise by filename with `localeCompare(..., {numeric: true, sensitivity: 'base'})`.

`r_value` files are then pulled out into the collapsible **R-value** strip (collapsed by default,
`showRValue = false`, rendered in the `wide` 1440×320 viewport); everything else goes into the plot
grid (720×450 cards).

**Only one R-value card is ever rendered**, for the combined (or, per 4d, first) r-value file. The
remaining logs still appear in the loaded-files list with an `r_value` badge but get no chart of
their own.

The **"Loaded N plot files"** panel lists every *chartable* plot file — i.e. those with a non-null
`plotKind` other than `stog` — with its kind badge, and lets the user hide individual charts.
Unrecognized files (e.g. `GTS_250K_XFQ1.csv`, `GTS_250K_FQ1partials.csv`) and detected `stog` files
appear **neither in the list nor in the count N**; there is no "unrecognized files" report anywhere
on the page. Hidden paths are removed from `plotFiles`, which — in static mode only (4d) — also
removes them from the combined R-value curve.

Per-file parse failures are non-fatal: the file gets a `parseError` string, its card renders no
chart, and a dismissible alert shows the message.

---

### Step 9 — Live Data: change detection and refresh

**The change signal** is a string fingerprint, `fileSignature()` in `browserData.js`:

```js
items.map(f => `${f.path}:${f.modified ?? ''}:${f.size ?? ''}:${f.plotKind ?? ''}`).sort().join('|')
```

so any added, removed, resized, or re-timestamped file changes it. `WATCH_INTERVAL_MS = 3000`
(3 seconds) is the single polling constant, shared by both runtimes.

**Flask mode** (`Dashboard.jsx`, the `watchFiles` effect): every 3 s, `GET /api/files`; if the new
signature differs from `signatureRef.current`, call `loadServerDashboard({silent: true})`, which
re-fetches `/api/plot/metadata` for every plot file and `/api/structure`. `pollInFlightRef` blocks
overlapping polls. `silent: true` suppresses the loading spinner and preserves the user's
hidden-chart selection. Freshness is `st_mtime` (seconds, float) plus byte size.

**Static mode** (`App.jsx`, the `fsAccess && watchFiles` effect): every 3 s,
`buildLocalRunFromHandle(dirHandleRef.current)` re-walks the picked directory handle recursively.
`FileSystemFileHandle.getFile()` returns a **fresh `File`** each call, so `size` and `lastModified`
reflect the current on-disk state; only metadata is read at this stage, not contents. If the
signature changed, `setLocalRun({...nextRun, runId: runIdRef.current})` — reusing the **same
`runId`** is what tells the Dashboard "same run, new files".

**In-place refresh (both modes).** When `localRun.runId` matches the previous one, `Dashboard.jsx`
does not tear down the view. Instead:

* Per plot file: if `prev.plotData` exists **and** `prev.modified === file.modified` **and**
  `prev.size === file.size`, the already-parsed data is reused. Otherwise the file is re-read
  (`readAndParseLocalPlotFile()` → `sourceFile.text()`, the **whole** file) and re-parsed.
* The `.rmc6f` is re-parsed (in `workers/localStructureWorker.js`) only if its mtime or size
  changed; the previous model summary stays on screen until the new one arrives, so there is no
  flash to empty.
* `InteractivePlot` receives `refreshKey = "${file.modified}:${file.size}"`, which forces the
  Flask-mode `/api/plot/data` fetch to re-run when a file changes.

**Cost of a refresh — every parse is a full parse.** Nothing here is incremental: an appending
`.log` is never tail-read, it is re-read and re-parsed from byte 0.

* In **Flask mode** a single silent poll re-issues `/api/plot/metadata` for **every** plot file, not
  only the changed ones, and each of those calls `make_plot()` — which parses the file with
  matplotlib and immediately throws the figure away. For `xray_sq`, `neutron_sq` and `bragg`,
  `make_plot()` additionally calls `read_rmc_csv(path)` a *second* time purely to read the header
  row for the axis label. `/api/plot/data` then parses the file again for the changed cards. Budget
  **2–3 full parses per file per refresh**. `related_r_value_logs()` also re-globs the log directory
  on every r-value request.
* In **static mode** a changed file is re-read in full via `sourceFile.text()` and re-parsed from
  scratch in the main thread (the `.rmc6f` goes to a worker; the plot files do not).

**Availability.** Static-mode Live Data requires the **File System Access API**
(`window.showDirectoryPicker`), i.e. Chromium browsers (Chrome, Edge, Arc, Opera).
`supportsFileSystemAccess()` additionally returns `false` in the Vite dev server when no
`VITE_API_BASE_URL` is set. Firefox/Safari fall back to `<input webkitdirectory>`, which yields a
**one-shot snapshot** with no watching — re-pick the folder to refresh. The demo run is likewise a
one-shot snapshot. In Flask mode Live Data works in any browser because the polling is server-side
file listing.

**Limitations of the change signal.** It has the resolution of the filesystem's mtime as exposed to
the runtime (milliseconds in the browser, `st_mtime` float seconds from `os.stat`) — a file rewritten
within the same tick to the same byte length would not be noticed. A file whose `stat()` raised
`OSError` carries `modified = size = null` and can never trigger a refresh (Step 1). And a file that
is being written when the poll lands may be read half-complete; the resulting parse error surfaces
as a per-card alert and is corrected on the next poll.

---

### Parameters and defaults — parsing and metrics

| Name | Value | Where | Meaning / units |
|---|---|---|---|
| `WATCH_INTERVAL_MS` | `3000` | `browserData.js` | Live Data poll period, ms (both runtimes) |
| `SUPPORTED_PATTERNS` | `*.csv, *.log, *.rmc6f, Frac*.txt, scale_ft.*, stog_input.dat, *.inp, *.sq, *.dat` | `app.py` | Flask file-listing globs (non-recursive) |
| `SUPPORTED_NAMES` | `scale_ft.gr, scale_ft.sq, scale_ft_rmc.fq, stog_input.dat` | `browserData.js` | extra exact names the browser accepts |
| demo file list | 10 files under `<base>/demo/` | `browserData.js` → `DEMO_FILES` | bundled GaTa₄Se₈ 250 K run; `lastModified = Date.now()` |
| `maxPoints` (dashboard) | `100` | `Dashboard.jsx` | atoms sampled for the model summary |
| `MAX_STRUCTURE_POINTS` | `1_000_000` | `app.py` | clamp on `/api/structure?maxPoints` (min 100) |
| structure sampling | site quota $\lfloor 100/n_\mathrm{sites}\rfloor$ (Flask) vs stride $\lceil N/100 \rceil$ (browser) | `app.py`, `browserData.js` | atoms; counts still use all atoms |
| circular-resultant floor | `1e-6` | `browserData.js` → `structureFromRmc6f` | guards $\sqrt{-2\ln R}$ |
| `.rmc6f` header window | `4000` characters | `browserData.js` → `readMovesMetadata` | fallback when `Atoms:` is not found early |
| symmetry tolerance (default) | `0.2` Å | `ModelSummary.jsx` | model-summary space-group search |
| symmetry ladder cap | `1.0` Å | `ModelSummary.jsx` → `toleranceLadder` | upper bound of the tolerance sweep |
| χ² clamp | `1e-12` | `browserData.js`, `app.py` | floor before `ln`; ⇒ y ≥ −27.63. **Absent** in `plots.py` |
| R-value log-header skip | first **2** lines | `read_chi()` / `readChi()` | fixed, not sniffed (line 2 holds the discarded `WEIGHT PARAMETERS`) |
| STOG header skip | first **2** lines | `read_stog()` / `readStog()` | fixed, not sniffed |
| R-value classification | inline `-\d{2,}\.log$` | `plots.py`, `browserData.js` | ≥ 2 digits required |
| R-value grouping/sorting | `R_VALUE_LOG_RE = ^(.+)-(\d{2,})\.log$` | `parsers.py`; mirrored by `rValueLogParts` in `Dashboard.jsx` | anchored stem + integer sequence |
| run-control head read | `131072` bytes | `pairFitTypes()` | 128 KiB per candidate `.dat` |
| run-control candidates | `6` | `runControlCandidates()` | max `.dat` files tried; stops at first non-empty map |
| datasets parsed | `8` | `parseRunSettings()` | max `*_DATA` blocks kept |
| Rwp chip precision | `4` significant digits | `Dashboard.jsx` | display only |
| min columns for Rwp | `3` | `plots.py`, `browserData.js` | x + 2 y columns |
| min columns to plot at all | `2` | `plots.py` → `_series_plot` | **Flask only**; fewer ⇒ HTTP 500 |
| PNG dpi | `150` (tests use 72) | `plot_to_png()` | `/api/plot` only |
| matplotlib figure size | `6.75 × 4.05 in` (`6.75 × 4.725` for STOG) | `plots.py` | `/api/plot` only |
| STOG PNG reference line | dashed at `y = 1` (`y = 0` for `.fq`) | `plots.py` → `_stog_plot` | `/api/plot` only |
| SVG viewport | `720 × 450` grid, `1440 × 320` wide | `InteractivePlot.jsx` | display only |
| axis padding | 5 % of range each side; ±1 if min = max; `[0,1]` if no finite values | `InteractivePlot.jsx` → `niceDomain` | display only |
| tick counts | 7 x × 6 y (grid), 11 x × 4 y (wide) | `InteractivePlot.jsx` | display only |
| tick multiplier rule | ν < 1.5 → 1; < 3 → 2; < 7 → 5; else 10 | `InteractivePlot.jsx` → `niceTicks` | display only |
| wheel zoom factors | `1.22` out / `0.82` in, clamped to the base x domain, min span `1e-9` | `InteractivePlot.jsx` → `zoom` | display only |
| drag-zoom threshold | `8` px per axis | `InteractivePlot.jsx` → `finishDrag` | display only |
| marker radius / palette | `2.8`; 8 colours cycled mod 8 | `InteractivePlot.jsx` | display only |
| tooltip precision | 4 significant digits; `toExponential(2)` outside `[0.01, 1e6)`; integer above 1000 | `InteractivePlot.jsx` → `formatNumber` | display only |
| watchdog enable | `settings.watchdogEnabled` (default off) + ≥ 2 points | `useWatchdog.js` | badge renders nothing otherwise |
| watchdog window | `max(min(N,5), ⌈0.2N⌉)` points | `heuristics.js` | recent-slope fit length, per log row |
| watchdog thresholds | `0.02` (ln units), total drop `0.1` | `heuristics.js` | ≈ 2 % and 10 % in the metric |
| watchdog LLM gate (data) | `≥ 200` new steps, or `\|Δ last\| ≥ ln(1.02)` | `heuristics.js` → `significantChange` | necessary, not sufficient — ANDed with the interval below |
| `watchdogIntervalMin` | `5` | `useWatchdog.js` | minutes; minimum spacing between LLM calls when the status has *not* flipped |
| watchdog LLM gate (overall) | `statusChanged \|\| (intervalElapsed && significantChange)`, and `baseUrl` + `model` set | `useWatchdog.js` | throttles narration only |
| watchdog stat rounding | 3 sig. digits (values), 2 (window delta) | `heuristics.js` → `watchdogStats` | prompt payload only |

---

### Caveats / what this is not

* **Classification is by filename only.** Nothing on this page inspects file content to decide what
  a file is. Rename a file and the app will plot it as something else; use a naming convention the
  patterns in Step 2 do not cover and the file is silently ignored. There is no "unrecognized files"
  report — the "Loaded N plot files" panel counts only *chartable* files.
* **"Rwp" is neither weighted nor conventionally normalized.** It is
  $\sqrt{\sum(\mathrm{col3}-\mathrm{col2})^2 / \sum \mathrm{col2}^2}$ with unit weights, and for
  RMCProfile CSVs column 2 is the *calculated* curve. Columns 4 and beyond are drawn but never
  enter it. The run's own `WEIGHT PARAMETERS` line is parsed past and discarded. Do not quote this
  number as $R_\mathrm{wp}$ in a paper without recomputing from the columns yourself. It is
  per-file; there is no combined R across datasets.
* **A degenerate R-factor is reported as unavailable, not as a number.** Both implementations sum
  only the points where the calculated *and* experimental values are finite, and return
  `None`/`null` — displayed as the chip **"Rwp —"** — when no such point exists or when the
  denominator is zero. Neither case is a fit quality, so no number is offered for one. A *partly*
  NaN column still produces a value, computed over the finite points only and therefore over
  silently fewer rows than the file has; the chip does not say how many (5a).
* **The R-value curve is "the last column of the log file".** No header is consulted, so neither
  `log(χ)` (the axis text) nor `ln(χ²)` (the codebase's description) is verified by the code. In the
  bundled demo run that column is headed `X_ray_(R)1`, an R-factor. The only verified statements are
  the positional column choice and the transform `ln(max(v, 1e-12))`. The second-to-last column is
  parsed into `chi_q` (Python only) and never displayed.
* **"Time steps" are log rows**, one per RMCProfile print/save period — not Monte-Carlo steps, and
  not uniform wall-clock time (the actual timestamps in column 1 are discarded).
* **The convergence badge is off by default.** It renders only when the user enables the watchdog in
  the AI-assistant settings *and* the series has ≥ 2 points; below 2 points the classification is
  `'unknown'`.
* **No residual/difference curve exists** anywhere on this page, and no data is offset, normalized,
  scaled, smoothed, rebinned, or subsampled before drawing. Every parsed point is handed to the
  renderer. The only reference line in the codebase is on the unused matplotlib STOG path.
* **What the tooltip shows is rounded.** 4 significant digits, or an integer above 1000. Read exact
  values from the file, not the hover readout.
* **The model summary's point cloud is a 100-atom subsample**, drawn by two different algorithms in
  the two runtimes (site-stratified quota in Flask, uniform stride in the browser). The element
  counts and totals shown on the card are over *all* atoms, so counts and cloud disagree by
  construction — that is intended, not a bug.
* **`fitType` (`FIT_TYPE` from the run-control `.dat`) has no effect on this page.** The
  `D(r)`-instead-of-`G(r)` relabel exists only in the `stog` branches of `plotMetadataFromFile`/
  `plotDataFromText`, and `stog` files are filtered out before those run. The `.dat`'s live effect
  here is on the AI-assistant context only, and only after the assistant page has been opened. In
  no case is any conversion applied to the numbers — the app never converts G(r) ↔ D(r) ↔ T(r).
* **`stog` files are detected but never charted here**, and the Flask file listing does not even
  return `.gr`/`.fq` files (only `*.sq` and literal `scale_ft.*`).
* **Static and Flask parsing are separate implementations, and which one runs depends on the
  deployment** (`isStaticMode()`, four branches — a dev server without `VITE_API_BASE_URL` uses the
  JavaScript parsers). They agree on detection rules 1–9, on the Rwp formula for clean numeric data
  (to floating-point round-off) *and* on its degenerate cases (5a), on the χ² clamp, and on the
  log-combination order. They differ on:
  the `stog` rule (any `.gr/.sq/.fq` vs three fixed names); the file-listing patterns; non-numeric
  and `NaN` handling in **all four** readers (raise vs `NaN` vs silently-dropped rows), including
  which line becomes the EXAFS header; blank-line handling and error line numbers in the CSV reader;
  the `.rmc6f` fallback pick; whether hidden logs are excluded
  from the combined R-value curve; the atom-sampling strategy; the presence of `basis`/`moves` in
  the structure payload; and element-name normalization (`parts[1].capitalize()` in Python vs the
  raw token in the browser, so `SE`/`se`/`Se` merge in Flask and may not in the browser).
* **The repository's reference run folders are not in the repository.** `data/` is gitignored, so
  examples such as `scale_ft_rmc.fq` cannot be reproduced from a clean clone; the reproducible run
  is [web_app/frontend/public/demo/](../../web_app/frontend/public/demo/).
* **Sample-backed tests skip in CI.** The GNSe reference dataset is likewise gitignored, so the
  assertions that pin real-file shapes, `Rwp > 0`, and `final_chi_r ≈ 0.00405`
  (`tests/test_plots.py`) do not run on CI — only the synthetic-fixture and pure-logic tests do
  (`AGENTS.md`, "Current known issues").
* **The matplotlib rendering path is effectively dead UI.** `PlotViewer.jsx` and `FileExplorer.jsx`
  are not mounted; `GET /api/plot` still works as an API. Its axis labels, its STOG reference line,
  and its unclamped `np.log(chi_r)` differ from what the dashboard draws — but its
  two-numeric-column precondition still gates the live Flask endpoints, because they build the
  figure to get the metrics.


## Plot rendering, interaction, and figure export

### What this section covers — plot rendering and export

The **Run Dashboard** ([`Dashboard.jsx`](../../web_app/frontend/src/components/Dashboard.jsx)) renders
every recognised RMCProfile output file in a run folder as its own chart card: Bragg profiles,
S(Q) (neutron and x-ray), EXAFS in Q- and R-space, PDF/xPDF fits and partials, plus a collapsible
full-width R-value (convergence) strip. Each chart is drawn by
[`InteractivePlot.jsx`](../../web_app/frontend/src/components/InteractivePlot.jsx), which is a
hand-written SVG renderer — **no plotting library is involved in the browser**. Every pixel
coordinate in the figure comes from the ~30 lines of arithmetic documented below.

Two things matter for a scientist reading this section:

1. **In the SVG chart path, nothing is smoothed, binned, decimated, or interpolated.** Every finite
   data point inside the current x-window becomes one vertex of the drawn polyline. The only lossy
   step is rounding the *pixel* coordinates to two decimals. This invariant covers the dashboard and
   Auto StoG charts only — the `<canvas>` KDE/3D panels of Step 15 *do* interpolate (bilinear
   heatmap smoothing) and *do* approximate (5-anchor colormaps).
2. **The exported figure is not a screenshot.** It is a re-serialisation of the same SVG DOM at a
   fixed nominal size, with the HTML chrome (legend, title, tooltip, $R_\mathrm{wp}$ chip) left
   behind. See Step 13 and the "Screen vs. export" list.

A second, matplotlib-based rendering path exists in
[`rmc_toolkits/plots.py`](../../rmc_toolkits/plots.py). It is the **package-API** figure path
(`rmc_toolkits.make_plot` / `rmc_toolkits.plot_to_png`, re-exported from
[`rmc_toolkits/__init__.py`](../../rmc_toolkits/__init__.py)), exposed over HTTP as `GET /api/plot`, and
described in Step 16. Two honest qualifications:

- The repo's only console script is `rmc-autoscale = rmc_toolkits.scaling_cli:main`
  ([`pyproject.toml`](../../pyproject.toml)), and `scaling_cli.py` imports no matplotlib and no
  `plots.py`. So "CLI figure path" would be wrong: `make_plot`/`plot_to_png` are used by
  [`app.py`](../../web_app/backend/app.py), by the tests, and by library callers, nothing else.
- `GET /api/plot` (the PNG endpoint) is **unreachable from the shipped UI**. Its only frontend
  consumer is [`PlotViewer.jsx`](../../web_app/frontend/src/components/PlotViewer.jsx), which is
  imported by no component — it is dead code in the current app. What the dashboard shows is always
  the SVG renderer.

#### Notation and units

| Symbol | Meaning |
| --- | --- |
| $x_i, y_i$ | the $i$-th data pair of one series, in the file's own units (Å⁻¹ for Q/k, Å for r, µs for ToF, arbitrary for intensity/S(Q)/G(r)) |
| $[x_0, x_1]$, $[y_0, y_1]$ | the current x- and y-**domains** (data units) |
| $p_x, p_y$ | coordinates in the SVG **user space** = viewBox units (dimensionless; *not* CSS pixels) |
| $L, R, T, B$ | left/right/top/bottom margins of the viewBox (user units) |
| $W, H$ | viewBox width/height (user units) |
| $W_p = W - L - R$, $H_p = H - T - B$ | plot-area width/height (user units) |
| $s$ | screen scale factor, CSS px per user unit |
| $W_c, H_c$ | canvas panel width/height in CSS px (Step 15) |

---

### Step 1 — Where the series arrays come from

**Inputs:** one output file (CSV, `.log`, or STOG text).
**Outputs:** a JSON-shaped `plotData` object

```
{ kind, title, metrics, xLabel, yLabel, series: [{ label, x: [...], y: [...] }] }
```

There are two producers, one per runtime mode:

- **Flask mode** — `GET /api/plot/data`
  ([`web_app/backend/app.py`](../../web_app/backend/app.py) → `plot_data()`), which calls
  `rmc_toolkits.parsers.read_rmc_csv` / `read_exafs_csv` / `read_chi` / `read_stog`. Column 0 of the
  CSV is the shared x array; every remaining column becomes one series, labelled by its header cell
  (an *empty* header cell is replaced by `Series {idx}`, 1-based over the data columns; a
  *duplicate* header is left as-is — see Step 3).
  `InteractivePlot` fetches this itself when no `plotData` prop is supplied.
- **Static mode** — [`browserData.js`](../../web_app/frontend/src/browserData.js) →
  `plotDataFromText()`, a JS re-implementation of the same parsing and labelling, run in the browser
  on a locally-picked file. `Dashboard.jsx` passes the result down as the `plotData` prop, and the
  card heading/metrics come from `plotMetadataFromFile()`.

#### 1a — Axis labels are hard-coded per plot kind, not read from the file

For every kind except the fallback branch, the axis strings are constants chosen by `kind` — i.e.
**an assumption about the file's units**, not a measurement of them. Both producers use the same
table (`app.py::plot_data`, `browserData.js::plotDataFromText`):

| kind | xLabel | yLabel |
| --- | --- | --- |
| `exafs_q` | `k (Å^{-1})` | `χ(k) k²` |
| `exafs_r` | `r (Å)` | `FT[χ(k) k²]` |
| `xpdf`, `npdf`, `pdf_partials` | `r (Å)` | `G(r)` |
| `xray_sq`, `neutron_sq` | `Q (Å^{-1})` | `S(Q)` |
| `bragg` | `ToF (µs)` or `Q (Å^{-1})` (see 1b) | `Intensity` |
| `r_value` | `Time steps` | `log(χ)` |
| `stog` | `r (Å)` if `.gr`, else `Q (Å^{-1})` | `G(r)`/`S(Q)` (Python) — see the third bullet in 1e |
| anything else | `cleanAxisLabel(header[0])` | `data` |

Only the fallback branch reads the file's own first-column header, through
`app.py::_clean_axis_label` / `browserData.js::cleanAxisLabel`, which apply four rewrites:
a header of exactly `Q` → `Q (Å^{-1})`; exactly `r` or `R` → `r (Å)`; the substring `(A^-1)` or
`(A^{-1})` → `(Å^{-1})`; the substring `(A)` → `(Å)`. Nothing else is touched. The `^{…}`
syntax in these strings is what `AxisLabel` renders as a superscript (Step 6).

#### 1b — The Bragg ToF-vs-Q axis is a substring guess on the column header

`plots.py::bragg_is_tof(header)` (Python) and `browserData.js::braggAxis(header)` (JS) both test the
**first column header** of the CSV against the case-insensitive regex `/tof|flight|time/`. A match
labels the x-axis `ToF (µs)`; anything else labels it `Q (Å^{-1})`. The two implementations agree
exactly. Two honest limitations:

- **No unit conversion is ever applied.** The raw column-0 values are plotted as-is under whichever
  label the regex picked. If the label is wrong, the axis is mislabelled and nothing else changes.
- The test is a plain substring search, so a header containing `lifetime`, `timestamp` or
  `time constant` would also be read as time-of-flight.

#### 1c — The R-value series is a concatenation of several log files

R-value ("convergence") charts are not one file. In **Flask mode**, `plot_data()` calls
`read_chi(related_r_value_logs(path))`; `parsers.related_r_value_logs` globs the *parent directory*
for every sibling matching `R_VALUE_LOG_RE = ^(.+)-(\d{2,})\.log$` with the same stem, and
`sort_r_value_logs` orders them by `(stem.lower(), sequence, name.lower())`. All of their chi values
are concatenated into one array. In **static mode**, `Dashboard.jsx::combineRValueFiles` does the
equivalent client-side: it `flatMap`s `plotData.series[0].y` over every parsed R-value file and
re-indexes `x = 0 … N-1`.

Consequences worth stating:

- The x-axis ("Time steps") is a **running index over the concatenation**, with no marker at the
  stage boundaries. A discontinuity in the curve is a run restart, not a physical jump.
- `final_chi_r` is taken from the **last** file only (`_chi_plot` returns `float(chi_r[-1])`;
  `combineRValueFiles` copies `lastParsed.plotData.metrics.final_chi_r`). It is the **raw, un-logged**
  last chi, even though the plotted quantity is its logarithm.
- `combineRValueFiles` bails out to `rValueFiles[0]` when there is only one file, when none of them
  is a locally-parsed file (Flask mode — the server does the concatenation instead), or while any of
  them is still parsing.

**Which column is chi.** `parsers.read_chi` skips the first **2 lines** of every log, splits each
remaining line on whitespace, requires ≥ 2 tokens, and takes `parts[-2]` as $\chi_Q$ and `parts[-1]`
as $\chi_r$. **Only $\chi_r$ is ever plotted**; $\chi_Q$ is parsed and discarded. Lines that fail
`float()` are skipped (Python) or dropped by `Number.isFinite` (JS `browserData.js::readChi`, which
reads only the last token). So the plotted quantity is *the last whitespace-separated field of each
post-header line*, and the sample index is an index into the **surviving** lines, not into the
file's lines.

Both interactive producers plot $\ln\!\big(\max(\chi_r, 10^{-12})\big)$ against that index.

#### 1d — Metrics

`metrics` carries at most `rwp`, `pdf_index` and `final_chi_r`. Their definitions, applicability
conditions and display rounding are given in Step 16, since the same numbers are produced by the
matplotlib path.

#### 1e — Known differences between the Python and JS producers

The two agree on structure, series ordering, and (with the exceptions below) axis-label strings and
metrics.

- **Non-numeric CSV cells.** `read_rmc_csv` (Python) calls `float(value)` and raises, so the request
  fails with an error message. `readRmcCsv` (JS) uses `values.map(Number)`, which yields `NaN`
  silently; those points are then dropped at draw time (Step 7) and the polyline bridges the gap.
- **R-value log clamp.** Both interactive producers clamp with $\max(\chi,10^{-12})$; the matplotlib
  path (`_chi_plot()`) uses `np.log(chi_r)` with **no** clamp, so a zero/negative entry gives
  $-\infty$/NaN there but $\ln 10^{-12} = -27.63$ in the interactive chart.
- **STOG y-label and card title.** The browser prefers the fit-function form declared in the
  run-control `.dat` file (`file.fitType`, e.g. `D(r)`, harvested by `browserData.js::pairFitTypes`
  → `fitTypeByFilename`) for *both* the y-label (`plotDataFromText`) and the card heading
  (`plotMetadataFromFile`). The Flask path always uses the extension default for the y-label
  (`"G(r)" if path.name.endswith(".gr") else "S(Q)"`) and the bare file name for the title (from
  `_stog_plot`). The same file therefore shows y-label `D(r)` / heading `D(r)` in static mode but
  y-label `G(r)` / heading `scale_ft.gr` in Flask mode. Python's `.gr` test is **case-sensitive**;
  the JS one lower-cases first.
- **CSV line numbering in error messages.** `readRmcCsv` (JS) filters blank lines *before* numbering
  rows, so its reported "line N" counts non-blank lines; `read_rmc_csv` (Python) numbers against the
  raw file. The EXAFS readers agree (both number against the raw line list).
- **EXAFS data-row detection.** `readExafsCsv` (JS) locates the first data row by requiring
  `Number.isFinite` on every token; `read_exafs_csv` (Python) uses `float()` in a `try/except`.
  Tokens Python accepts as non-finite floats (`inf`, `nan`) make a row "numeric" for Python but not
  for JS, which can shift the detected header line.
- **Strict column count (both).** Every data row must have exactly `len(labels)` values or the read
  raises — a hard failure, not a skipped row. In static mode this surfaces as the card's parse
  error; in Flask mode as a 500 from `/api/plot/data`.

Whichever producer ran, `InteractivePlot` treats `x` and `y` as opaque numbers in the file's own
units and performs **no unit conversion whatsoever**.

**Code:** `app.py::plot_data`, `app.py::_clean_axis_label`, `plots.py::bragg_is_tof`,
`parsers.py::read_chi`, `parsers.py::related_r_value_logs`, `parsers.py::sort_r_value_logs`,
`browserData.js::plotDataFromText`, `browserData.js::plotMetadataFromFile`, `browserData.js::readChi`,
`browserData.js::braggAxis`, `browserData.js::cleanAxisLabel`,
`Dashboard.jsx::combineRValueFiles`, `InteractivePlot.jsx` (the `useEffect` that GETs
`/api/plot/data`).

---

### Step 2 — Series ordering, roles, colours and marker style

**Inputs:** `plotData.series`. **Outputs:** `orderedSeries`, each with a resolved `color`,
`marker` and `guide` flag.

`InteractivePlot.jsx` → `orderedSeries` (a `useMemo`) applies, in order:

1. **Split off guides.** Any series with `role: 'guide'` is removed from the data set. Guides are
   theory/reference lines (e.g. the `S → 1` line and the $-\langle b\rangle^2$ level on the Auto
   StoG page, Step 4a). They are appended *after* the data series, so they draw on top; they default
   to the colour `#98a2b3` (`GUIDE_STROKE`), **never consume a palette slot**, and are **excluded
   from the hover search** (Step 9).
   Their CSS class is `series-path series-path--guide`. In
   [`InteractivePlot.css`](../../web_app/frontend/src/components/InteractivePlot.css)
   `.series-path--guide` declares `stroke-dasharray: 7 5`, `stroke-width: 1.3` and `opacity: 0.9`,
   **but the `stroke-width: 1.3` never takes effect**: `.series-path` is declared *later* in the same
   file with `stroke-width: 2` and has identical specificity (0,0,1,0), so the later rule wins.
   Guides therefore draw at **2 units, non-scaling**, dashed `7 5`, at opacity 0.9. This propagates
   into every export, because `figureExport.js::inlineComputedStyles` freezes `getComputedStyle`
   values — the exported SVG/PNG carries `stroke-width:2` on guides too.
2. **Classify measured vs. calculated.** `isExperimental(label)` is the regular expression
   `/exp/i` applied to the series label. It matches `Experiment`, `F(Q)_Expt`,
   `X_ray_exp_renorm`, … and would also match any other label containing the letters "exp".
3. **Reorder.** If *both* an experimental and a non-experimental series exist (`paired`), the
   experimental ones are moved first so they are drawn underneath; otherwise the original order is
   kept.
4. **Colour.** `entry.color || palette[index % 8]` with
   `palette = ['#1f6fd6', '#e8590c', '#099268', '#d6336c', '#6741d9', '#66a80f', '#0c8599', '#e67700']`
   and `index` counted over the reordered *data* series only.
5. **Marker flag.** `marker = paired && isExperimental(label)`. Marker series are drawn as hollow
   circles instead of a line (Step 7).

The memo's only dependency is `effectivePlot`, so the assignment is stable while the user toggles
the legend but is recomputed whenever the series list itself changes — see the colour caveat at the
end of this section.

**Code:** `InteractivePlot.jsx` → `orderedSeries`, `isExperimental`, `palette`, `GUIDE_STROKE`;
`InteractivePlot.css` → `.series-path`, `.series-path--guide`.

---

### Step 3 — Legend visibility filter (and the label-identity assumption)

**Inputs:** `orderedSeries`, the `hidden` set of labels. **Outputs:** `visibleSeries`.

Each legend chip is an HTML `<button>`; clicking it toggles the series label in a `Set`
(`setHidden`). `visibleSeries = orderedSeries.filter(s => !hidden.has(s.label))`.

Series may start hidden: when the caller supplies `plotData` directly, any entry flagged
`defaultHidden: true` is seeded into `hidden` (used by the Auto StoG page to park the "Measured
(unscaled)" and "Pre-enforcement fit" diagnostics one click away).

**Series identity is the label string, not the column index.** `hidden` is a `Set` of labels, the
legend buttons use `key={series.label}`, and the emitted `<path>` elements use `key={series.label}`
(Step 7). Two series that share a header — common in partials CSVs, or any file with a repeated
column name — therefore toggle together from a single legend chip, collide as React reconciliation
keys, and cannot be told apart in the tooltip. Nothing deduplicates labels: `plot_data()` substitutes
`Series {idx}` only for an *empty* header cell, never for a duplicate one.

**This filter feeds autoscale.** Everything in Step 4 is computed from `visibleSeries` only, so
hiding the measured curve rescales both axes to the calculated curve alone, and *vice versa*.
Hiding a series also changes the wheel-zoom clamp (`baseX`, Step 11). Guide series **do** count
toward the domains, so a guide line that extends past the data widens the axes.

**Code:** `InteractivePlot.jsx` → `visibleSeries`, the `.plot-legend` buttons.

---

### Step 4 — Axis domain selection (autoscale and padding)

**Inputs:** `visibleSeries`, the user zoom state `xDomain` / `yDomain`, an optional
`plotData.initialYDomain`, and the `fullExtent` toggle.
**Outputs:** `domains = { x, y, baseX, baseY }` in data units.

The padding helper is

```js
niceDomain(values):  m = min(finite), M = max(finite)
                     if m === M: m -= 1, M += 1
                     pad = (M - m) * 0.05
                     return [m - pad, M + pad]
```

i.e. a symmetric **5 % of the data range** added on each side:

$$[\,\min_i v_i - 0.05\,\Delta,\; \max_i v_i + 0.05\,\Delta\,],\qquad \Delta = \max_i v_i - \min_i v_i .$$

Degenerate cases: an empty/all-non-finite set returns $[0,1]$; a constant series is first widened to
$[v-1, v+1]$ and *then* padded by 5 % of the widened range, giving $[v-1.1, v+1.1]$.

The x- and y-domains are then chosen as:

- $\mathrm{baseX} = \texttt{niceDomain}$ over the concatenated `x` of all visible series.
- $[x_0, x_1] = \texttt{xDomain}$ if the user has zoomed, else $\mathrm{baseX}$.
- $\mathrm{baseY} = \texttt{niceDomain}$ over the `y` values **of points whose `x` lies inside the
  current x-window**, i.e. over $\{\,y_i : x_0 \le x_i \le x_1\,\}$. If that set is empty the
  fallback is all `y` of all visible series.
- $[y_0, y_1] = \texttt{yDomain}$ (user y-zoom) **else** `plotData.initialYDomain` (a caller-supplied
  default window, used for the Auto StoG low-r zoom; suppressed while `fullExtent` is true) **else**
  $\mathrm{baseY}$.

The in-window filter is implemented as
`series.y.filter((_, index) => series.x[index] >= currentX[0] && series.x[index] <= currentX[1])` —
an **index-based pairing** that assumes `x` and `y` have the same length. Two silent exclusions
follow: for caller-supplied `plotData` (Auto StoG guides, hand-built series) a length mismatch makes
the comparison `undefined >= number`, which is false, so those `y` values drop out of the
y-autoscale; and a NaN `x` fails both comparisons, so its `y` is excluded too. `niceDomain` then
independently drops any non-finite `y` via `Number.isFinite`.

Consequence worth stating plainly: **the y-axis re-autoscales on every x-zoom** unless the user has
explicitly y-zoomed. Panning/zooming in x therefore changes the apparent amplitude of features.
`initialYDomain` is used verbatim — it is *not* padded and *not* rounded to nice numbers.

#### 4a — The Auto StoG guide constants

[`AutoStogPage.jsx`](../../web_app/frontend/src/components/AutoStogPage.jsx) is the only caller that
supplies guides, `defaultHidden` and `initialYDomain`, so its numbers belong here:

| chart | guide series (`role: 'guide'`) | x span | colour |
| --- | --- | --- | --- |
| S(Q) (`variant="wide"`) | `S → 1` at $y = 1$ | `q[0]` → `q[-1]` | default `#98a2b3` |
| S(Q) | `Level <fmt(level,5)>` at $y = $ `guides.level` | `guides.levelWindow` | `#4c7df0` |
| S(Q) | `S(0) FZ target <fmt(s0Target,3)>` | `q[0]` → `min(q[0] + 1.5, q[-1])` (Å⁻¹) | `#0c8599` |
| $G_K(r)$ | `−⟨b⟩² theory`, flat at `level = guides.gkLowR` | `0` → `min(RMAX_DISPLAY = 8, r[-1])` Å | default |
| $D(r)$ | `−4πρ₀⟨b⟩²r theory`, sloped $0 \to$ `drSlope · guideEnd` | `0` → `min(8, r[-1])` Å | default |

The $G_K(r)$ chart is the only one with `initialYDomain`, set to
$[\,2.1\,\ell,\; -3.2\,\ell\,]$ with $\ell = $ `guides.gkLowR` (negative), i.e. an asymmetric window
around the $-\langle b\rangle^2$ level. `defaultHidden: true` is set on `Measured (unscaled)` (S(Q))
and on `Pre-enforcement fit` (both real-space charts).

Two consequences: the guide **labels embed display-rounded numbers** (`fmt(level, 5)`,
`fmt(s0Target, 3)`), so the legend text is itself a rounded quantity, not the working value; and
because guides count toward the domains (Step 3) and the real-space guides start at $r = 0$, the
x-axis of the $G_K(r)$ and $D(r)$ charts is widened to include the origin even when the data start
later.

**Code:** `InteractivePlot.jsx` → `niceDomain`, `domains` (`useMemo`); `AutoStogPage.jsx` →
`sqPlot`, `gkPlot`, `drPlot`, `RMAX_DISPLAY`.

---

### Step 5 — Viewport geometry and the data → pixel affine map

**Inputs:** `domains`, the `variant` prop. **Outputs:** `xScale`, `yScale`, `xInvert`, `yInvert`.

Two fixed viewBoxes, chosen by `variant`:

| variant | $W \times H$ | $L$ | $R$ | $T$ | $B$ | plot area $W_p \times H_p$ |
| --- | --- | --- | --- | --- | --- | --- |
| grid card (default) | 720 × 450 (8:5) | 60 | 18 | 16 | 58 | 642 × 376 |
| `wide` (R-value strip, Auto StoG S(Q)) | 1440 × 320 | 64 | 20 | 18 | 58 | 1356 × 244 |

The `<svg>` is `width: 100%; height: auto` with a CSS `aspect-ratio` that exactly matches the
viewBox (8/5 and 1440/320), so `preserveAspectRatio` never letterboxes and the user-space → CSS-px
map is a single uniform scale $s = W_{\mathrm{CSS}}/W$. The element also carries
`touch-action: none`, so touch drags zoom the chart instead of scrolling the page.

The data → user-space map is affine and axis-independent:

$$p_x(x) = L + \frac{x - x_0}{x_1 - x_0}\,W_p, \qquad
p_y(y) = T + H_p - \frac{y - y_0}{y_1 - y_0}\,H_p .$$

The $y$ form is written with the subtraction because SVG $y$ grows downward: $y_0$ maps to the
bottom edge $T + H_p$ and $y_1$ to the top edge $T$. In the code the denominators are written
`(domains.x[1] - domains.x[0] || 1)`, i.e. a zero (or NaN) span is replaced by 1, collapsing the
series onto the left/bottom edge instead of producing NaN geometry.

The inverses, used by every pointer interaction, are

$$x(p_x) = x_0 + \frac{p_x - L}{W_p}\,(x_1 - x_0), \qquad
y(p_y) = y_0 + \frac{T + H_p - p_y}{H_p}\,(y_1 - y_0).$$

Note the inverses carry **no `|| 1` guard** on the span, so with a degenerate domain they return the
domain minimum rather than NaN.

No logarithmic axes exist anywhere in the chart renderer. The R-value plot is "log" only because the
*data* were log-transformed upstream (Step 1c); the axis itself is linear in $\ln\chi$.

**Code:** `InteractivePlot.jsx` → `view`, `plotWidth`, `plotHeight`, `xScale`, `yScale`,
`xInvert`, `yInvert`; `InteractivePlot.css` → `.interactive-plot svg`.

---

### Step 6 — Tick generation ("nice numbers"), tick labels and tick geometry

**Inputs:** a domain $[a,b]$ and a target tick count $n$. **Outputs:** `{ ticks, step }`.

`niceTicks(domain, count = 5)` is the classic 1–2–5 decade algorithm (the default `count = 5` is
never used by the chart — both axes pass an explicit target):

$$\mathrm{raw} = \frac{b-a}{\max(1,n)},\qquad
m = 10^{\lfloor \log_{10}\mathrm{raw}\rfloor},\qquad
\nu = \frac{\mathrm{raw}}{m} \in [1,10)$$

$$\mu = \begin{cases}
1 & \nu < 1.5\\
2 & 1.5 \le \nu < 3\\
5 & 3 \le \nu < 7\\
10 & \nu \ge 7
\end{cases}
\qquad \mathrm{step} = \mu\, m$$

Ticks start at the first multiple of `step` inside the domain,
$t_0 = \lceil a/\mathrm{step}\rceil \cdot \mathrm{step}$,
and are generated by repeated addition while $t \le b + 10^{-6}\,\mathrm{step}$
(the $10^{-6}$ slop admits a final tick that floating-point accumulation would otherwise push just
past $b$). Any tick with $|t| < 10^{-9}\,\mathrm{step}$ is snapped to exactly `0`. If the span is
non-finite or $\le 0$, the result is the single tick $[a]$ with `step = 1`.

Because ticks are produced by accumulation rather than $t_0 + k\cdot\mathrm{step}$, the last tick of a
long axis carries $O(N\varepsilon)$ float drift; this affects only the printed label digits, never
the drawn data.

Target counts (these are *targets*; the realised count is whatever the 1–2–5 step yields):

| axis | grid card | wide |
| --- | --- | --- |
| x | 7 | 11 |
| y | 6 | 4 |

**Tick label text** — `formatTick(value, step, axisMax)`, where
$\mathrm{axisMax} = \max_k |t_k|$ over that axis's ticks so a whole axis shares one unit suffix.
Branch selection uses $\mathrm{magnitude} = \max(|\mathrm{axisMax}|, |v|)$, but the number printed is the
raw $v$, so a small tick on a large axis is printed with the axis's suffix (e.g. a $5\times10^5$
tick on a $2\times10^6$ axis prints `0.5M`):

- $|v| < 10^{-12}$ → `"0"` (this short-circuit runs *before* any `axisMax` logic).
- magnitude $\ge 10^{9}$ or $< 10^{-4}$ → `toExponential(1)` with `e+` collapsed to `e`
  (e.g. `1.2e9` for the large branch, `3.0e-5` for the small one).
- magnitude $\ge 10^{6}$ → $v/10^6$ with an `M` suffix, at `decimalsForStep(step / 1e6)` decimals.
- magnitude $\ge 10^{4}$ → $v/10^3$ with a `k` suffix, at `decimalsForStep(step / 1e3)` decimals.
- otherwise → fixed point at `decimalsForStep(step)` decimals.

**The decimal count is per branch, not global**: the `M` and `k` branches feed `decimalsForStep` the
*rescaled* step, so the suffixed labels keep consecutive ticks distinct; the exponential branch is a
flat one digit regardless of `step`.

`decimalsForStep(step)`: with $e = \lfloor\log_{10}\mathrm{step}\rfloor$, zero decimals when
$e \ge 0$, else $\min(6, -e)$ decimals. It falls back to **2 decimals** when `step` is non-finite or
$\le 0$ — a branch unreachable from `niceTicks`, which always returns a positive step.

**Grid and tick geometry are asymmetric by design.** For each y-tick: a full-width horizontal grid
line (`.plot-grid-line`) from $x=L$ to $x=W-R$, plus a label at $(L-10,\; p_y(t)+4.5)$ with
`text-anchor="end"`. For each x-tick: only a **5-unit** tick stub from $y = T+H_p$ to
$y = T+H_p+5$ (`.plot-tick-mark`), plus a label at $(p_x(t),\; H-36)$ anchored middle. There is no
vertical grid.

**The grid is not clipped.** The `.plot-bg` rect, the grid lines, the tick marks, the tick labels,
the `.plot-frame` border and the axis labels are all emitted *outside* the `<g clipPath=…>` group;
only the series paths and the hover marks are inside it (Step 7). The final tick admitted by the
`step × 1e-6` slop therefore draws its grid line a hair beyond the plot rectangle.

**Axis labels** are rendered by `AxisLabel`, which recognises the pattern `…^{…}…` and emits a
`<tspan baselineShift="super" fontSize="70%">` for the exponent (so `Q (Å^{-1})` renders as
Q (Å⁻¹)). Only one superscript group per label is handled, and the greedy `(.*)` prefix means a
label with two groups renders only the last one as a superscript. The x-label sits at
$(L + W_p/2,\; H-10)$; the y-label at $(18,\; T + H_p/2)$ rotated −90°. Tooltip text uses a
different converter, `labelToText`, which substitutes Unicode superscript glyphs (⁻¹ etc.).

**Code:** `InteractivePlot.jsx` → `niceTicks`, `formatTick`, `decimalsForStep`, `AxisLabel`,
`labelToText`, the `yTicks`/`xTicks` render blocks.

---

### Step 7 — Path construction: clipping, NaN handling, markers, no decimation

**Inputs:** `visibleSeries`, `domains`. **Outputs:** one SVG `path` `d` string per series
(`seriesShapes`, memoised on `[visibleSeries, domains, view…]` so hover re-renders do not rebuild
them).

For each series, in file order:

1. **X-window filter.** A point is skipped if `x < domains.x[0] || x > domains.x[1]`.
2. **Finiteness filter.** A point is skipped if `x` or `y` is not finite (NaN comparisons are false,
   so NaN x survives step 1 and is caught here).
3. **Map** through $p_x, p_y$ of Step 5 and round each coordinate with `toFixed(2)`.
4. **Emit.** Line series become `M x y L x y L x y …`. Marker series become one closed circle per
   point, built from two 180° arcs:
   `M (c_x - r) c_y  a r r 0 1 0 (2r) 0  a r r 0 1 0 (-2r) 0` with $r = 2.8$ user units
   (`MARKER_RADIUS`), filled with `var(--plot-bg)` and stroked in the series colour (a hollow
   circle) at `stroke-width: 1`, non-scaling.

**There is no decimation, no LTTB, no min/max binning, no averaging, and no point budget.** A 2649-
point S(Q) (the size asserted in [`tests/test_backend_api.py`](../../tests/test_backend_api.py)
`test_plot_metadata_and_data_endpoints`) produces a 2649-vertex path. The consequences to be honest
about:

- **Pixel-coordinate rounding is the only quantisation.** Two decimals in user space is
  $0.01/W_p$ of the axis: for a 0–30 Å⁻¹ Q-axis on a grid card that is
  $0.01/642 \times 30 \approx 4.7\times10^{-4}$ Å⁻¹ in x, and $0.01/376$ of the y-range in y.
  Two data points closer than that render at the same coordinate; nothing is *dropped*, the polyline
  just has a zero-length segment. **No feature can be hidden by decimation, because there is none**
  — but a feature narrower than one rendered device pixel can still be invisible on screen, exactly
  as with any un-decimated plot.
- **Non-finite points are bridged, not gapped.** A NaN in the middle of a series does not break the
  polyline: the neighbouring points are joined by a straight segment, so a hole in the data looks
  like a straight interpolation. Nothing in the UI marks this.
- **Edge clipping is by point, not by segment.** Points outside $[x_0,x_1]$ are discarded whole, so
  the segment that would cross the domain boundary is never drawn. After a zoom, the curve therefore
  starts at the first *in-range* sample rather than at the axis, leaving a sub-sample-spacing gap at
  each edge. A `clipPath` (unique per chart instance via `useId()`) additionally crops all series
  and hover marks to the plot rectangle so a y-zoom cannot draw over the axes.
- **Dense marker series merge.** With 2649 experimental points across 642 user units the spacing is
  0.24 units while the marker radius is 2.8 — the hollow circles overlap ~20-fold and read as a
  solid band. This is cosmetic, not numerical.

Line style comes from CSS, not from the data: `.series-path` is `fill: none; stroke-width: 2;
stroke-linejoin/linecap: round; vector-effect: non-scaling-stroke`. `non-scaling-stroke` means the
2-unit stroke is realised as **2 CSS px regardless of how wide the chart is displayed**, so relative
line weight changes with the panel size (and again on export — Step 13). Guides inherit the same
2 px (Step 2), differing only in dash pattern and opacity.

**Code:** `InteractivePlot.jsx` → `seriesShapes`, `MARKER_RADIUS`, the `<g clipPath=…>` block;
`InteractivePlot.css` → `.series-path`, `.series-markers`.

---

### Step 8 — Pointer → user-space coordinates

**Inputs:** a `PointerEvent` (`clientX`, `clientY` in CSS px).
**Outputs:** $(p_x, p_y)$ in viewBox user units.

`pointerToView` uses the browser's own matrix rather than manual arithmetic:

```js
const m = svgRef.current.getScreenCTM();
const p = new DOMPoint(event.clientX, event.clientY).matrixTransform(m.inverse());
```

This is exact for any page zoom, CSS transform, or scroll position; it returns `{x: L, y: T}` if the
ref or the CTM is unavailable. Two clamps keep interactions inside the axes:
`clampPlotX` to $[L,\;W-R]$ and `clampPlotY` to $[T,\;H-B]$.

**Every pixel threshold in the interaction code (the 8-unit drag threshold below) is expressed in
user units, not CSS pixels.** At display scale $s$ a threshold of 8 user units is $8s$ CSS px — about
5.3 CSS px for a grid card rendered 480 px wide, 8 px at 720 px wide, 16 px at 1440 px wide.

**Code:** `InteractivePlot.jsx` → `pointerToView`, `pointerToViewX`, `clampPlotX`, `clampPlotY`.

---

### Step 9 — Nearest-point hover search

**Inputs:** pointer position; `visibleSeries` minus guides. **Outputs:** `hover = { x, px, values[] }`.

The metric is **x-distance in data units only** — not Euclidean, not screen-space, and with no y
term at all:

$$i^\star_{(s)} = \arg\min_{i}\; \big| x^{(s)}_i - \hat{x} \big|, \qquad
\hat{x} = x(p_x)\ \text{from Step 5's inverse.}$$

Because $x \mapsto p_x$ is affine and strictly monotonic, minimising the data-x distance is
equivalent to minimising the screen-x distance; the choice of units does not change the winner.

Details that matter:

- The search is an **independent linear scan per series** (`series.x.forEach`), $O(N)$ per series per
  `pointermove`. No binary search is used even though the x arrays are monotonic in practice, and no
  spatial index is built. For a handful of few-thousand-point series this is imperceptible; for very
  large files it is the dominant hover cost.
- The scan runs over the **full** `series.x` array, including points outside the current x-window.
  With a zoom that excludes all data, the reported point can lie outside the visible axes; its dot is
  clipped away but the tooltip still prints its value.
- **There is no proximity cut-off.** The nearest point is always found, however far the cursor is;
  hovering an empty region snaps to the closest endpoint.
- **Ties go to the lowest index.** `best` is initialised to `0` and the comparison is strictly
  `distance < bestDistance`, so an exact tie keeps the earlier sample.
- **A series with no finite x still reports index 0.** If every `x` is NaN (all comparisons false) or
  the array is empty, `best` stays `0`, `series.x[0]` may be `undefined`, and the row's `cx`/`cy`
  become NaN — a NaN dot and a NaN tooltip entry, with no error.
- **Series with different x-grids report different x.** Each series answers with *its own* nearest
  sample, but the crosshair position and the tooltip header x are taken from `values[0]` — the first
  visible non-guide series. If two series are on different grids (e.g. an experimental file and a
  finer calculated grid), the tooltip header shows series 0's x while the other rows show y at their
  own, slightly different x. The tooltip does not display each row's x. The fallback
  `values[0]?.x ?? dataX` (the raw cursor data-x) applies **only** when there is no first series at
  all.
- Guides (`role: 'guide'`) are excluded from the search entirely.
- Hover is cleared when the pointer leaves the SVG (`onPointerLeave`) or moves outside
  $[L,\,W-R]$ horizontally; there is no vertical bound, so hovering over the tick-label strip still
  produces a crosshair.

Rendering: a dashed vertical line at $p_x$ of the series-0 match, one filled circle of radius 3.6
user units per series, and an **HTML** tooltip positioned at
$\mathrm{left} = 100\,p_x/W\,\% + 14\mathrm{px}$ (or mirrored to `right` when $p_x > W/2$) so it sits on
the emptier side. Tooltip numbers use `formatNumber`: `toExponential(2)` when $|v|\ge10^6$ or
$0<|v|<0.01$; rounded and locale-grouped when $|v|\ge1000$; otherwise `toPrecision(4)`.

**Code:** `InteractivePlot.jsx` → `nearestHover`, `formatNumber`, the `.plot-tooltip` block.

---

### Step 10 — Drag-to-zoom: rectangle → data range inversion

**Inputs:** pointer-down/move/up in user space. **Outputs:** new `xDomain` and/or `yDomain`.

- `startDrag` rejects a press outside $[L,\;W-R]$ horizontally, tries `setPointerCapture` (failure is
  tolerated), stores `{x0, y0, x1, y1}` all clamped into the plot rectangle, and clears hover.
- `moveDrag` updates `(x1, y1)`; the preview rectangle is drawn at
  $(\min(x_0,x_1), \min(y_0,y_1))$ with size $(|x_1-x_0|, |y_1-y_0|)$.
- `finishDrag` computes the clamped extremes
  $x_{lo}, x_{hi}, y_{lo}, y_{hi}$ and applies **per-axis thresholds of 8 user units**:

$$\mathrm{zoomX} \Leftrightarrow x_{hi}-x_{lo} > 8, \qquad
\mathrm{zoomY} \Leftrightarrow y_{hi}-y_{lo} > 8 .$$

  A genuine box zooms both axes; a thin horizontal swipe zooms x only; a thin vertical swipe zooms y
  only; a click (both below 8) zooms nothing and simply clears the drag.

The inversion uses Step 5's inverses evaluated with the **pre-zoom** domain:

$$\texttt{xDomain} \leftarrow [\,x(x_{lo}),\; x(x_{hi})\,], \qquad
\texttt{yDomain} \leftarrow [\,y(y_{hi}),\; y(y_{lo})\,].$$

Note the y swap: screen y grows downward, so the *larger* pixel value $y_{hi}$ is the *smaller* data
value and must become the domain minimum. The resulting domains are exactly the dragged rectangle —
**no 5 % padding and no nice-number rounding is applied to a user zoom**, so the axis limits after a
drag are arbitrary reals and the ticks are re-derived from them by Step 6.

**Interaction-ordering hazard.** The press-time coordinates `drag.x0`/`drag.y0` are captured and
clamped in *user space* at pointerdown, but they are inverted at release through whichever
`xInvert`/`yInvert` closure is current then. A wheel notch during an active drag (Step 11) replaces
`xDomain`, so the release inverts old pixels against a new domain and the resulting window does not
correspond to the rectangle the user saw. Nothing guards against this.

There is no limit on how far in one may zoom, and no pan gesture.

**Code:** `InteractivePlot.jsx` → `startDrag`, `moveDrag`, `finishDrag`, the `.zoom-selection` rect.

---

### Step 11 — Wheel zoom (x only)

**Inputs:** a `wheel` event. **Outputs:** a new `xDomain`.

`zoom(event)` calls `preventDefault()`, converts the pointer to a data coordinate
$c = x(\mathrm{clampPlotX}(p_x))$, and rescales the **x-domain only** about that point:

$$\Delta' = (x_1 - x_0)\cdot f,\qquad
f = \begin{cases} 1.22 & \Delta_{\mathrm{wheel}} > 0 \ (\text{zoom out})\\ 0.82 & \mathrm{otherwise}\end{cases}$$

$$\texttt{xDomain} \leftarrow \big[\max(\mathrm{baseX}_0,\; c - \Delta'/2),\;
\min(\mathrm{baseX}_1,\; c + \Delta'/2)\big]$$

applied only if the resulting span exceeds $10^{-9}$.

Honest details:

- The factors are **not exact inverses**: $1.22 \times 0.82 = 1.0004$, so a zoom-out followed by a
  zoom-in does not return to the identical window (0.04 % drift per pair).
- The wheel **delta magnitude is ignored** — only its sign. A high-resolution trackpad emits many
  small events, each applying a full 22 %/18 % step, so trackpad zooming is much faster than a
  notched mouse wheel.
- The clamp is applied to each endpoint independently against `baseX` (the padded full extent of the
  *visible* series). Zooming out near an axis end therefore produces an asymmetric window rather than
  re-centring, and the zoom-out limit is the un-zoomed extent — you cannot wheel out past the data.
- The wheel never touches y. Because `yDomain` stays null until the user drags a y-box, the y-axis
  re-autoscales to the new x-window on every wheel notch (Step 4).
- `zoom` clears **neither** an in-progress `drag` nor the `hover` state; see the hazard in Step 10.

**Code:** `InteractivePlot.jsx` → `zoom`.

---

### Step 12 — Reset, double-click, and view-state lifetime

- The **Reset zoom** chip appears whenever `xDomain`, `yDomain`, or `fullExtent` is set, and clears
  all three.
- **Double-click** always clears `xDomain` and `yDomain`. Additionally, when the chart was supplied
  with an `initialYDomain` *and* was already at its default view, it toggles `fullExtent`, i.e. it
  switches between the caller's preferred y-window and the full 5 %-padded data extent. Without an
  `initialYDomain` a double-click is a plain zoom reset.
- When a *new* `plotData` object arrives (static mode / Auto StoG), the component resets `hidden`
  (re-seeded from `defaultHidden`), `xDomain`, `yDomain`, `hover`, `drag`, and `fullExtent` during
  render, using React's "adjust state when a prop changes" pattern rather than an effect.
- **In Flask mode the reset is asymmetric.** The corresponding branch runs after a successful fetch
  when `file.path` changed, and it resets `hidden`, `xDomain`, `hover` and `drag` but **omits
  `setYDomain(null)` and `setFullExtent(false)`**. A y-zoom therefore survives a switch to a
  different file, and the new file is drawn into the previous file's y-window. This is a real
  inconsistency between the two modes, not a documented design choice.
- `refreshKey` is `"<modified>:<size>"` of the file; it is a dependency of the fetch effect, so Live
  Data re-polls the series when the file changes on disk. There is no client- or server-side cache:
  every `refreshKey` change re-reads and re-parses the whole file (Step 16).

**Code:** `InteractivePlot.jsx` → the `onDoubleClick` handler, the `plotData !== lastPlotData`
render-time reset, the fetch `useEffect`; `Dashboard.jsx` → `renderPlotBody` (passes `refreshKey`).

---

### Step 13 — Export path: SVG serialisation and PNG rasterisation

All figure export lives in
[`figureExport.js`](../../web_app/frontend/src/figureExport.js). Charts offer two formats
(`CHART_SAVE_OPTIONS = png | svg`) through the shared
[`SaveMenu.jsx`](../../web_app/frontend/src/components/SaveMenu.jsx) badge (which saves immediately when
only one format is offered, otherwise opens a small menu that closes on outside-pointerdown or
Escape).

**13a — Clone and freeze styles.** `standaloneSvg(svgElement)`:

1. Reads the nominal size from `svgElement.viewBox.baseVal` → 720 × 450 (or 1440 × 320), falling back
   to `clientWidth/clientHeight` and finally to 720 × 450.
2. `cloneNode(true)`.
3. `inlineComputedStyles(source, clone)` walks both trees in lockstep over `.children` and, for each
   node **including the root `<svg>` itself**, copies `window.getComputedStyle(source)` values for
   exactly these 16 properties onto the clone's `style` attribute: `fill, fill-opacity, stroke,
   stroke-width, stroke-dasharray, stroke-linejoin, stroke-linecap, opacity, font-family, font-size,
   font-weight, font-style, letter-spacing, text-anchor, font-variant-numeric, vector-effect`.
   This is what makes the file self-contained: the charts are styled through external CSS classes and
   CSS custom properties (`var(--plot-bg)`, `var(--accent)`), and the computed value resolves those to
   literal colours. **Any presentation property not on that list is silently dropped** (e.g.
   `stroke-opacity`, `paint-order`, `dominant-baseline`, filters). The recursion walks `.children`
   only — it assumes the clone is structurally identical (true for `cloneNode(true)`) and never
   touches text nodes.
4. Stamps `xmlns`, `width`, `height` attributes. **`width`/`height` are the viewBox numbers**, so the
   exported SVG's intrinsic size is 720 × 450 (or 1440 × 320) CSS px regardless of how large the
   chart was on screen. This is the fact the line-weight argument in item 5 of "Screen vs. export"
   rests on.

Two further properties of the clone: it keeps its `class` attributes (harmless, since no stylesheet
travels with it and every visual property is already inlined), and it keeps its
`clip-path="url(#plot-clip-…)"` reference. That reference resolves because `<defs>` is cloned with
it, but the id is only unique **within one page** (`useId()`), so two exported figures pasted into a
single host document can collide on the id.

**13b — SVG file.** `svgToSvgString` runs `XMLSerializer().serializeToString(clone)` and prefixes
`<?xml version="1.0" encoding="UTF-8"?>`. Saved as `<sanitised title>.svg` with MIME
`image/svg+xml;charset=utf-8`. This is **true vector output** — the same coordinates, to the same
two decimals, as on screen. No background rectangle is added beyond the chart's own `.plot-bg` rect,
so the margins are transparent.

**13c — PNG file.** `svgToPngBlob(svgElement, { scale = 2, background = '#ffffff' })`:

1. Serialise the styled clone and wrap it in a
   `data:image/svg+xml;charset=utf-8,<encodeURIComponent(...)>` URL. `encodeURIComponent` typically
   inflates the payload 2–3×, and this is the practical size ceiling of the PNG path: a chart with
   very long path data can exceed the browser's data-URI limit, and the failure surfaces as the
   `Could not rasterize the figure` rejection below.
2. Load it into an `Image`; a load error rejects with `Could not rasterize the figure`.
3. Create a canvas of $\lfloor W\cdot\mathrm{scale}\rceil \times \lfloor H\cdot\mathrm{scale}\rceil$ —
   **1440 × 900** for a grid card, **2880 × 640** for a wide strip, at the default `scale = 2`.
4. Fill it with opaque `#ffffff` **first**, then `drawImage(image, 0, 0, canvas.width, canvas.height)`.
   Because the fill precedes the draw, the chart's own `.plot-bg` rect composites *over* white rather
   than replacing it. The four-argument `drawImage` stretches the image to the whole canvas box; this
   is distortion-free only because the canvas is exactly `viewBox × scale` in both dimensions.
5. `canvas.toBlob(…, 'image/png')`; a null blob rejects with `Could not encode the figure`.

Both rejection messages are caught by `InteractivePlot.jsx::saveFigure` and shown in the chart's own
error state.

**The chart PNG resolution is fixed by the viewBox and `scale`, not by the browser window and not by
`window.devicePixelRatio`.** `devicePixelRatio` is never consulted on the chart path; the exported
raster is identical on a 1× and a 3× display. `scale` is a parameter but no caller overrides it, so
**2× is the effective constant** for charts. (The 3× export in this codebase belongs to the
canvas-based KDE/3D panels — Step 15.)

**13d — Filenames.** `sanitizeFilename(name)` strips `^{…}` braces down to the exponent text,
replaces every run of characters outside `[A-Za-z0-9_.-]` with `_`, trims leading/trailing
underscores, and falls back to `figure`. Because JS `\w` is ASCII-only, `Å`, `χ` and `µ` are
replaced by `_`, so filenames are always ASCII-safe. The name passed in is
`effectivePlot.title || file.name`.

**13e — Download.** `downloadBlob` creates an object URL, clicks a synthetic `<a download>`, removes
it, and revokes the URL on the next tick.

**Code:** `figureExport.js` → `standaloneSvg`, `inlineComputedStyles`, `STYLE_PROPS`,
`svgToSvgString`, `svgToPngBlob`, `loadImage`, `saveSvgAsPng`, `saveSvgAsSvg`, `saveSvgFigure`,
`sanitizeFilename`, `downloadBlob`; `InteractivePlot.jsx` → `saveFigure`.

---

### Step 14 — "Save all figures" → one ZIP

**Inputs:** the live DOM of the dashboard page. **Outputs:** `figures-png.zip` or `figures-svg.zip`.

`Dashboard.jsx` → `handleSaveAllFigures(format)` queries the page root for every `.plot-card`, takes
its `.interactive-plot svg` node and the `.plot-card-header h3` text as the name, and hands the list
to `saveSvgFiguresAsZip`. A `savingAll` flag guards re-entry. Selection and naming rules:

- **A card with no `.interactive-plot svg` is skipped entirely** — a still-parsing card, or one
  showing a parse error, contributes nothing to the archive.
- The fallback name is `figure-<n>`, where `n` counts only the cards that actually yielded an SVG,
  not the card's position on the page.
- Because it reads the *live* nodes, each chart is exported **in its current zoom and legend state**;
  the R-value strip is included only when the panel is expanded; and charts hidden via the
  "Loaded N plot files" panel are excluded (they are not rendered).

`saveSvgFiguresAsZip` sanitises each name and disambiguates repeats by appending `_2`, `_3`, …
from a `Map` of use counts — which can itself collide with a genuine figure whose sanitised title
already ends in `_2`. It then encodes each figure (`TextEncoder` for SVG, `svgToPngBlob` at the
default 2× for PNG) and calls `buildZip`.

**Performance and memory, honestly:** the PNG entries are produced **sequentially** with `await`
inside the loop — N full rasterisations one after another, each allocating a complete
`viewBox × 2` canvas — and `buildZip` allocates **one `ArrayBuffer` sized to the entire archive**, so
the whole ZIP is materialised in memory before the download starts. That is the practical limit on
"save all" for a dashboard carrying many large charts.

[`zipArchive.js`](../../web_app/frontend/src/zipArchive.js) → `buildZip(entries)` is a dependency-free
writer: local file headers (signature `0x04034b50`, version-needed 20, general-purpose flags 0),
**compression method 0 = store (no compression)**, CRC-32 over each entry (standard reflected
polynomial `0xEDB88320`, table built once), a central directory, and an EOCD record. Its limits:

- Last-mod time and date are both written as **0**, which is not a legal DOS date; extractors render
  it inconsistently (1979-11-30 is common, not 1980-01-01).
- **General-purpose flag bit 11 (the UTF-8 filename flag) is not set**, although the names are
  encoded with `TextEncoder` (UTF-8). This is harmless only because `sanitizeFilename` guarantees
  ASCII names — that dependency is load-bearing.
- The EOCD entry counts and the filename-length fields are `uint16`, so the archive is capped at
  **65535 entries**, and sizes/offsets are `uint32` (not ZIP64), so it tops out at 4 GiB.
- It writes no directory entries.

**Code:** `Dashboard.jsx::handleSaveAllFigures`, `figureExport.js::saveSvgFiguresAsZip`,
`zipArchive.js::buildZip`, `zipArchive.js::crc32`.

---

### Step 15 — Canvas figures: device pixel ratio, the 3× export, and colormaps

Two pages render to `<canvas>` rather than SVG, so they export as PNG only. They are **not** the same
machinery and are described separately. Their science is documented in the KDE/3D sections; what
belongs here is the geometry of drawing and saving them.

#### 15a — `StructurePage.jsx` (KDE slice, Slab-In-Cell, folded-unit-cell 3D)

[`StructurePage.jsx`](../../web_app/frontend/src/components/StructurePage.jsx) owns
`PANEL_SAVE_OPTIONS = [{ id: 'png', hint: '1×' }, { id: 'png3x', hint: '3×' }]`, the two 2D panels,
`save2dPanel`, `makePlaneMapper` and the LUT heatmap. Everything in this subsection is
`StructurePage.jsx` only.

**Live canvas sizing.** Each 2D panel measures its CSS box with `getBoundingClientRect()`, floors it,
clamps it to a minimum (**320 × 260** CSS px for the KDE slice, **220 × 260** for the slab), reads
`dpr = window.devicePixelRatio || 1`, sets the backing store to $W_c\!\cdot\!dpr \times H_c\!\cdot\!dpr$,
and applies `ctx.setTransform(dpr,0,0,dpr,0,0)` so all drawing code works in CSS units. A
`ResizeObserver` bumps a `sizeTick` so the panel is redrawn crisp once layout settles.

**Native PNG export** calls `canvas.toBlob` on that live element, so its pixel dimensions are
$W_c\!\cdot\!dpr \times H_c\!\cdot\!dpr$ — i.e. **the native PNG size depends on the user's display**.

**3× PNG export** (`save2dPanel(..., format='png3x')`) does not upscale pixels: it re-measures the
element, allocates an offscreen canvas of $3W_c \times 3H_c$, sets `ctx.setTransform(3,0,0,3,0,0)`,
and **re-runs the same draw function**, so vector content (outlines, contours, text) is genuinely
re-rasterised at 3×.

**Canvas plane mapping.** `makeProjectedPlane(u⃗, v⃗, uvPoints)` + `makePlaneMapper(plane, W_c, H_c,
padding = 18)` is a **two-stage** map, not a single axis-aligned fit. Given plane fractional
coordinates $(u_f, v_f)$ and the plane's Cartesian basis vectors $\vec u, \vec v$ (Å, obtained from
the fractional slice directions through the unit-cell vectors):

$$\text{(1) projection:}\qquad x = u_x u_f + v_x v_f,\qquad y = u_y u_f + v_y v_f$$

$$\text{(2) padded uniform fit:}\qquad
k = \min\!\left(\frac{W_c - 2\cdot 18}{\Delta x},\; \frac{H_c - 2\cdot 18}{\Delta y}\right)$$

$$p_x = \frac{W_c - k\,\Delta x}{2} + (x - x_{\min})\,k,\qquad
p_y = \frac{H_c - k\,\Delta y}{2} + (y_{\max} - y)\,k$$

where $\Delta x = x_{\max}-x_{\min}$ and $\Delta y = y_{\max}-y_{\min}$ are the spans of the
**projected polygon bounds** computed by `makeProjectedPlane` over `planePolygon` — *not* the spans
of $u_f$ and $v_f$. Stage (1) is what makes the map a general linear projection: for a skewed
(non-orthogonal) cell the plane axes are not perpendicular on screen, and the drawn cell outline is a
parallelogram. Guards: `makePlaneMapper`'s own `spanX`/`spanY` fall back to `1` when zero (`|| 1`).
The mapper computes no aspect ratio; separately, `makePlane` (line 139) and `makeProjectedPlane`
(line 159) each store an `aspect = Δx / max(Δy, 1e-9)`, and it is `makeProjectedPlane`'s that
`slicePanelGeometry` (line 531) reads back as `planeAspect` (line 535) and `sideAspect` (line 540)
and hands to the panels' CSS `--panel-aspect` custom property — `planeAspect` for the KDE-slice panel
(line 1321), `sideAspect` for the "Slab In Cell" panel (line 1351), and `Math.max(planeAspect, 1)`
for the "Folded Unit Cell" 3D panel (line 1364), which clamps a taller-than-wide plane back to at
least square. A layout quantity, not part of the point map.

`invert(p_x, p_y)` undoes both stages: first the algebraic inverse of the offset/scale
($x = (p_x - \mathrm{offset}_x)/k + x_{\min}$, $y = y_{\max} - (p_y - \mathrm{offset}_y)/k$), then a 2×2
determinant solve back to $(u_f, v_f)$ with $\det = u_x v_y - v_x u_y$, returning $(0,0)$ when
$|\det| < 10^{-12}$. This is what the slab drag handler uses.

**Heatmap drawing** (`drawKdeSlice`). The heatmap branch runs **only** when
`density && grid > 0 && kde.vmax > kde.vmin`; otherwise no heatmap is drawn at all and the panel
shows the text "Computing KDE…" or "No atoms in this slab". A perfectly flat density field therefore
renders **empty, not uniform**. When it does run:

1. `span = vmax - vmin || 1`; each cell is normalised, quantised to the 256 LUT levels with
   `Math.max(0, Math.min(255, Math.round(normalized * 255)))`, and written into a `grid × grid`
   `ImageData` with **alpha forced to 255** for every cell.
2. `ctx.imageSmoothingEnabled = true`, then the projected cell outline is used as a `ctx.clip()`
   path.
3. The unit-square placement is an explicit matrix built from `mapper.map` of the extent corners:
   `ctx.transform(uEdge.x-origin.x, uEdge.y-origin.y, vEdge.x-origin.x, vEdge.y-origin.y, origin.x,
   origin.y)` with `origin = map(xMin,yMin)`, `uEdge = map(xMax,yMin)`, `vEdge = map(xMin,yMax)`,
   followed by `drawImage(off, 0,0,grid,grid, 0,0,1,1)`.
4. Contour polylines (when enabled) are stroked **over** the smoothed image at `lineWidth = 1` in the
   theme contour colour, then the cell outline is stroked in the theme border colour.

**So the on-screen heatmap is bilinearly interpolated from the KDE grid, not shown at native grid
resolution** — the one place in this section where the renderer smooths data.

**3D view.** The live three.js renderer uses `renderer.setPixelRatio(Math.min(window.devicePixelRatio,
2))` — **capped at 2** — and is created with `preserveDrawingBuffer: true` so the frame stays
readable. `captureModelBlob(3)` sets `renderer.setPixelRatio(3)`, calls `setSize(w, h, false)`,
re-renders, captures with `toBlob`, then restores the previous ratio, re-sizes and re-renders.

#### 15b — `PcaKdePage.jsx` (WebGL only)

[`PcaKdePage.jsx`](../../web_app/frontend/src/components/PcaKdePage.jsx) has **no 2D slice or slab
canvas, no `save2dPanel`, and no `PANEL_SAVE_OPTIONS`**. Its own constant is
`SAVE_OPTIONS = [{ id: 'png', label: 'Standard PNG', hint: '1×' }, { id: 'png3x', label: 'High
quality PNG', hint: '3×' }]` — the same two ids, different labels. `saveMainView(format)` operates on
the three.js renderer canvas: the standard path re-renders and calls `saveCanvasAsPng` on
`renderer.domElement` (again relying on `preserveDrawingBuffer`); the `png3x` path reads the current
size, sets `setPixelRatio(3)`, `setSize(size.x, size.y, false)`, re-renders, `toBlob`s, then restores
the previous ratio and re-renders. There is no offscreen 2D re-draw anywhere on this page.

**Its live pixel ratio is uncapped**: `renderer.setPixelRatio(window.devicePixelRatio || 1)`.

#### 15c — What "3×" actually means

The **3× export is absolute** (3 × the CSS size) on both pages, while the **native export is
display-dependent and page-dependent**:

| page | live pixel ratio | native PNG | "3×" PNG | ratio of the two |
| --- | --- | --- | --- | --- |
| `StructurePage.jsx` 2D panels | `dpr` | $W_c\,dpr \times H_c\,dpr$ | $3W_c \times 3H_c$ | $3/dpr$ |
| `StructurePage.jsx` 3D | $\min(dpr, 2)$ | $W\min(dpr,2)$ | $3W$ | $3/\min(dpr,2)$ — never below 1.5 |
| `PcaKdePage.jsx` 3D | $dpr$ (uncapped) | $W\,dpr$ | $3W$ | $3/dpr$ |

So on a 2× display the "3×" file is only 1.5× the native one everywhere; on a 3× display it equals
the native file on `PcaKdePage` but is still 1.5× on `StructurePage`'s 3D view because of the
$\min(dpr,2)$ cap.

#### 15d — Colormaps

[`colormaps.js`](../../web_app/frontend/src/colormaps.js) holds **5-anchor piecewise-linear RGB ramps**
named `viridis`, `magma`, `seismic`, `reds`, `greys`, baked once into a 256-entry
`Uint8ClampedArray` LUT and cached per name. For LUT index $i$: $t = i/255$, $\sigma = 4t$, lower
anchor $\lfloor\sigma\rfloor$, upper anchor $\min(4, \lfloor\sigma\rfloor + 1)$, and each channel is
linearly interpolated by $\sigma - \lfloor\sigma\rfloor$.

Honest qualifications:

- `viridis`, `magma`, `seismic` and (loosely) `reds` are 5-anchor samples of the matplotlib maps of
  the same name — **not** the published 256-entry tables, so a figure using them is not
  byte-identical to a matplotlib render.
- **`greys` is not matplotlib `Greys`.** The anchors run
  `[12,13,15] → [70,74,80] → [130,136,144] → [195,200,207] → [245,247,250]`, i.e. **dark → light**
  (matplotlib `Greys` runs white → black, so this is closer to `Greys_r`) and slightly **blue-tinted**
  rather than neutral. A reader who assumed a matplotlib-equivalent ramp would read the density scale
  backwards.
- An **unrecognised colormap name silently falls back to `viridis`** (`ANCHORS[name] ||
  ANCHORS.viridis`) and is cached under the unrecognised name.
- The `upper` anchor is clamped to the last anchor, so $t = 1$ interpolates the last anchor against
  itself with `frac = 0`.
- The LUT is a `Uint8ClampedArray`, so every interpolated channel is **rounded to an integer**
  (ties-to-even) — a second quantisation on top of the 256 levels.
- Two consumer entry points, with different responsibilities:
  `getLut(name)` returns the raw table (used by `StructurePage.jsx::drawKdeSlice`,
  `PcaKdePage.jsx` and `orientationSphere.js::colorbarGradient`), and the caller does its own
  normalisation and clamping; `sampleColormap(name, v)` maps one value by
  `clamp(round(v · 255), 0, 255)` — it **clamps but does not normalise**, so out-of-range inputs
  saturate at the ends and callers must map their data into $[0,1]$ themselves (`PcaKdePage.jsx`
  does exactly that with `t < 0 ? 0 : t > 1 ? 1 : t`).

**Code:** `StructurePage.jsx` → `PANEL_SAVE_OPTIONS`, `save2dPanel`, `captureModelBlob`,
`drawKdeSlice`, `drawSlab`, `makeProjectedPlane`, `makePlaneMapper`;
`PcaKdePage.jsx` → `SAVE_OPTIONS`, `saveMainView`; `colormaps.js` → `ANCHORS`, `buildLut`, `getLut`,
`sampleColormap`; `figureExport.js` → `canvasToPngBlob`, `saveCanvasAsPng`.

---

### Step 16 — The matplotlib path (`plots.py`) and the metrics

[`rmc_toolkits/plots.py`](../../rmc_toolkits/plots.py) is the **package-API** figure path —
`make_plot()` / `plot_to_png()`, re-exported from `rmc_toolkits/__init__.py` — exposed over HTTP as
`GET /api/plot` (PNG) and `GET /api/plot/metadata` (JSON). It is a *separate* renderer from the
interactive chart and agrees with it only in what data it reads. As noted in the intro, no console
script uses it (`rmc-autoscale` is the only one, and it never imports matplotlib), and its PNG
endpoint has no live consumer in the shipped UI.

#### 16a — `detect_plot_kind(path)` and its precedence

Classification is by **file name**, tested in this exact order (first match wins). It is not all
regex — two branches are substring/`endswith` tests and the last is a set membership test:

| # | test on `Path(path).name` | kind |
| --- | --- | --- |
| 1 | `re.search(r"-EXAFS-.+_Q_OUTPUT\.csv$")` | `exafs_q` |
| 2 | `re.search(r"-EXAFS-.+_R_OUTPUT\.csv$")` | `exafs_r` |
| 3 | `re.search(r"_FT_XFQ\d+\.csv$")` | `xpdf` |
| 4 | `"PDF" in name and name.endswith(".csv")` | `pdf_partials` if `"PDFpartials" in name` else `npdf` |
| 5 | `name.endswith("_FQ1.csv")` | `xray_sq` |
| 6 | `name.endswith("_SQ1.csv")` | `neutron_sq` |
| 7 | `re.search(r"_bragg(?:_.+)?\.csv$")` | `bragg` |
| 8 | `re.search(r"-\d{2,}\.log$")` | `r_value` |
| 9 | `name in {"scale_ft.gr", "scale_ft.sq", "scale_ft_rmc.fq"}` | `stog` |
| — | otherwise | `None` → `/api/plot/data` answers **400** |

Note rule 4's precedence: it fires **before** the `_FQ1`/`_SQ1` tests, so a file whose name contains
both `PDF` and `_SQ1.csv` is classified `npdf`, not `neutron_sq`.

The JS counterpart is `browserData.js::detectPlotKind`, and the Python side is pinned by
[`tests/test_plots.py`](../../tests/test_plots.py) `test_detect_plot_kind_for_supported_outputs`.
**The two detectors agree on every rule except `stog`:** Python accepts only the three legacy exact
names above, while the browser accepts *any* file with a `.gr`, `.sq` or `.fq` extension,
case-insensitively (`/\.(gr|sq|fq)$/i`). The browser therefore *classifies* as `stog` a set of
`.gr`/`.sq`/`.fq` files for which Flask returns **400 `Unsupported plot file type: <name>`**
(`app.py` → `plot_data()`). **Neither runtime charts them**: `Dashboard.jsx`'s
`isDashboardPlotFile` drops every `stog` file before any metadata or parse call (Step 2), and
`readAndParseLocalPlotFile` / `plotDataFromText` / `plotMetadataFromFile` have no other production
caller — a `grep` across `web_app/frontend/src` finds only `Dashboard.jsx` and
[`__tests__/browserData.test.js`](../../web_app/frontend/src/__tests__/browserData.test.js). The
classification difference is observable only in the run diagnostics — `plotFileCount` counts every
file with a truthy `plotKind`, `stog` included — never in a chart.

#### 16b — The figure builders

- `_series_plot()` **raises `ValueError` when the file has fewer than 2 numeric columns**. Otherwise
  it builds `plt.figure(figsize=(6.75, 4.05))` (inches) with one `ax.plot` per CSV column index
  `1 ≤ idx < len(series.data)` against column 0, all at `lw=1.0, alpha=0.65` — i.e. **every curve is
  65 % opaque**, so overlapping curves blend, unlike the opaque browser strokes. A header row with
  more labels than data columns therefore silently drops the extra labels (the JSON API mirrors this
  and substitutes `Series {idx}` for a *blank* header cell). Legend `loc=1` (upper right),
  `fontsize=9`, `frameon=False`; axis labels at 11 pt; `fig.suptitle(title, fontsize=14)`.
  $R_\mathrm{wp}$ is computed here only when `calculate_rwp` is set *and* the file has ≥ 3 columns.
- `_chi_plot()` plots `np.log(chi_r)` against the implicit index — **unclamped**, unlike the
  interactive path's $\max(\chi,10^{-12})$ guard — over the concatenated sibling logs of Step 1c, and
  returns the metric `final_chi_r = float(chi_r[-1])`, the **raw, un-logged** last chi. The dashboard
  and the static path report the same quantity.
- `_stog_plot()` uses `figsize=(6.75, 4.725)` and a single **opaque** red curve (`alpha=1.0,
  color="r"`). Its dashed black horizontal reference is `ax.hlines(y, data[0][0], data[0][-1],
  ls="--", lw=0.5, color="black")` with $y = 0$ for `.fq` and $y = 1$ otherwise — spanning the
  **first and last x samples**, not the axis limits. It emits **no `fig.suptitle`**: the only text
  identifying the file is the legend label.
- `npdf` additionally reports `pdf_index` from `parsers.pdf_index`, the filename regex
  `PDF(\d+)\.csv$` (default 0, matching the JS `pdfIndex`), and takes its title from
  `path.stem.split("_")[-1]`; `pdf_partials` uses the same title rule.
- The `xray_sq`, `neutron_sq` and `bragg` branches pass the **raw first CSV header** (or the
  ToF/Q string) as the matplotlib x-label and leave the y-label at `_series_plot`'s default
  `"data"` — so the matplotlib axis text differs from the interactive chart's hard-coded strings of
  Step 1a.
- `plot_to_png(result, dpi=150)` writes with `bbox_inches="tight"`, so the raster is
  $6.75\times150 = 1012$ px wide *before* the tight crop trims whitespace; the final pixel size is
  therefore layout-dependent and not a fixed number.

#### 16c — Cost of the JSON endpoint

`/api/plot/data` calls `make_plot(path)` purely to obtain `title` and `metrics`, then immediately
`close_plot()`s it. **The interactive endpoint therefore builds and discards a complete matplotlib
figure on every request**, and there is no server-side cache: each Live Data `refreshKey` change
re-reads and re-parses the whole file (and re-globs the sibling logs, for R-value).

#### 16d — Where it differs from the browser chart

| | browser SVG chart | matplotlib PNG |
| --- | --- | --- |
| axis limits | explicit 5 % padded `niceDomain`, or the user's zoom | matplotlib autoscale (the code sets no `set_xlim`/`set_ylim`), padded by the `axes.xmargin`/`axes.ymargin` rcParams — 0.05 in a default install; `plots.py` sets **no** rcParams |
| ticks | 1–2–5 `niceTicks`, target 7/6 (or 11/4) | matplotlib's own `MaxNLocator` |
| curve opacity | opaque, including guides (Step 2) | `alpha=0.65` on every series **except** `_stog_plot`, which is opaque (`alpha=1.0`) |
| markers | hollow circles for `*exp*` series when paired | none; all series are lines |
| R-value log | $\ln\max(\chi,10^{-12})$ | $\ln\chi$, unclamped |
| title | HTML card header only (not in the figure) | `fig.suptitle` for `_series_plot`/`_chi_plot`; **`_stog_plot` adds none** (its file name appears only in the legend) |
| axis label text | hard-coded per kind (Step 1a) | same for EXAFS/PDF; raw CSV header for `_FQ1`/`_SQ1`/bragg x, and y-label `"data"` |
| legend | HTML chips outside the SVG (not exported) | inside the figure, upper right |
| interactivity | zoom/hover/legend toggles | none (static PNG) |
| output size | fixed 1440 × 900 px (PNG) or true vector (SVG) | ≈1012 px wide at dpi 150, then `bbox_inches='tight'` |

#### 16e — The $R_\mathrm{wp}$ metric

The formula is computed identically in both languages — `parsers.rwp(x, observed, fitted)` and
`browserData.js::rwp`:

$$R_\mathrm{wp} = \sqrt{\frac{\sum_i (f_i - o_i)^2}{\sum_i o_i^2}}$$

with $o$ = CSV column 1 and $f$ = CSV column 2 **by position**, so $R_\mathrm{wp}$ is only meaningful when
the file's column order really is (x, observed, calculated). The `x` argument is accepted and
ignored by both. The two runtimes agree to floating-point round-off. The conditions and the display
rounding matter as much as the formula:

- It is computed **only** for kinds `xpdf`, `npdf`, `xray_sq`, `neutron_sq` and `bragg` — never for
  `exafs_q`, `exafs_r`, `pdf_partials`, `stog` or `r_value` — and **only when the CSV has ≥ 3
  columns**. Otherwise the card shows no chip.
- When the denominator $\sum_i o_i^2$ is exactly zero, both implementations return **`0.0` by
  convention** — a silent "perfect fit" reading rather than an error or a blank.
- The dashboard chip prints `Number(rwp).toPrecision(4)` (4 significant figures). The unmounted
  `PlotViewer.jsx` metric strip prints every metric at `toPrecision(5)`.

**Code:** `plots.py` → `detect_plot_kind`, `bragg_is_tof`, `_series_plot`, `_stog_plot`, `_chi_plot`,
`make_plot`, `plot_to_png`, `close_plot`; `parsers.py` → `rwp`, `pdf_index`;
`app.py` → `plot_file`, `plot_metadata`, `plot_data`; `browserData.js` → `detectPlotKind`, `rwp`,
`pdfIndex`.

---

### Parameters and defaults — plot rendering and export

**Geometry and scales**

| Parameter | Value | Units | Where |
| --- | --- | --- | --- |
| viewBox, grid card | 720 × 450 | user units | `InteractivePlot.jsx` → `view` |
| viewBox, `wide` variant | 1440 × 320 | user units | `InteractivePlot.jsx` → `view` |
| margins L/R/T/B, grid | 60 / 18 / 16 / 58 | user units | `InteractivePlot.jsx` → `view` |
| margins L/R/T/B, wide | 64 / 20 / 18 / 58 | user units | `InteractivePlot.jsx` → `view` |
| domain padding | 5 % of data range, each side | — | `niceDomain` |
| degenerate-series widening | ±1 before padding | data units | `niceDomain` |
| touch behaviour | `touch-action: none` on the `<svg>` | — | `InteractivePlot.css` |

**Ticks and labels**

| Parameter | Value | Units | Where |
| --- | --- | --- | --- |
| tick step family | 1, 2, 5, 10 × 10ⁿ | — | `niceTicks` |
| tick multiplier thresholds | ν < 1.5 → 1; < 3 → 2; < 7 → 5; else 10 | — | `niceTicks` |
| `niceTicks` default count | 5 (never used — both axes pass a target) | — | `niceTicks` |
| target tick count, x | 7 (grid) / 11 (wide) | — | `InteractivePlot.jsx` |
| target tick count, y | 6 (grid) / 4 (wide) | — | `InteractivePlot.jsx` |
| tick loop slop | `step × 1e-6` | — | `niceTicks` |
| zero-snap threshold | `step × 1e-9` | — | `niceTicks` |
| tick-label decimals, fixed-point | 0 if ⌊log₁₀step⌋ ≥ 0, else min(6, −⌊log₁₀step⌋) | — | `decimalsForStep` |
| tick-label decimals, `k` / `M` branch | same rule on `step/1e3` / `step/1e6` | — | `formatTick` |
| tick-label decimals, exponential branch | fixed 1 digit | — | `formatTick` |
| `decimalsForStep` fallback | 2 (non-finite or ≤ 0 step; unreachable from `niceTicks`) | — | `decimalsForStep` |
| k/M suffix thresholds | 1e4 → k, 1e6 → M, ≥1e9 or <1e-4 → exponential | — | `formatTick` |
| zero short-circuit | \|v\| < 1e-12 → `"0"` | — | `formatTick` |
| y-tick label offset | x − 10, y + 4.5; `text-anchor: end` | user units | `InteractivePlot.jsx` |
| x-tick label position | x = `xScale(t)`, y = H − 36; anchor middle | user units | `InteractivePlot.jsx` |
| x tick stub length | 5 | user units | `InteractivePlot.jsx` |
| x-axis label position | (L + W_p/2, H − 10) | user units | `AxisLabel` |
| y-axis label position | (18, T + H_p/2), rotated −90° | user units | `AxisLabel` |
| tick text | `#475467`, 14.5 px, weight 450, `tabular-nums` | CSS px | `.plot-tick` |
| axis-label text | `#1d2939`, 16 px, weight 600, letter-spacing −0.01em | CSS px | `.axis-label` |

**Series, hover and interaction**

| Parameter | Value | Units | Where |
| --- | --- | --- | --- |
| palette | 8 colours, cycled by position | — | `palette` |
| marker radius | 2.8 | user units | `MARKER_RADIUS` |
| marker stroke | 1, non-scaling; fill `var(--plot-bg)` | CSS px | `.series-markers` |
| series stroke width | **2 for all line series** (the guide rule's declared 1.3 is overridden by the later, equal-specificity `.series-path`), non-scaling | CSS px | `InteractivePlot.css` |
| guide stroke | `#98a2b3` (`GUIDE_STROKE`), dash `7 5`, opacity 0.9 | — | `GUIDE_STROKE`, `.series-path--guide` |
| path coordinate precision | 2 decimals | user units | `seriesShapes` |
| grid line | `rgba(16,24,40,0.06)`, width 1 | — | `.plot-grid-line` |
| tick mark | `rgba(16,24,40,0.35)`, width 1 | — | `.plot-tick-mark` |
| plot frame | `rgba(16,24,40,0.18)`, width 1, `fill: none` | — | `.plot-frame` |
| hover line | `rgba(16,24,40,0.3)`, width 1, dash `3 4` | — | `.hover-line` |
| hover dot | radius 3.6 user units, stroke `#fff` at 1.6 | — | `.hover-dot` |
| zoom-selection rect | fill `var(--accent-soft)`, stroke `var(--accent)` at 1.2 | — | `.zoom-selection` |
| tooltip offset / placement | 14 px from the crosshair, `top: 0.7rem`, `max-width: min(17rem, 100% − 1.5rem)` | CSS px | `.plot-tooltip` |
| legend strip | `max-height: 4.2rem`, overflow auto | — | `.plot-legend` |
| tooltip number format | exp(2) if ≥1e6 or <0.01; rounded if ≥1000; else 4 sig figs | — | `formatNumber` |
| drag-zoom threshold | 8, per axis | user units | `finishDrag` |
| wheel zoom-out / zoom-in factor | 1.22 / 0.82 | — | `zoom` |
| minimum zoom span | 1e-9 | data units | `zoom` |
| plot-area background | `#f8fafc` (default/dark theme) / `#ffffff` (light) | — | `index.css` → `--plot-bg` |
| $R_\mathrm{wp}$ chip precision | 4 significant figures | — | `Dashboard.jsx` |
| `PlotViewer` metric precision | 5 significant figures (dead code) | — | `PlotViewer.jsx` |

**Auto StoG chart constants**

| Parameter | Value | Units | Where |
| --- | --- | --- | --- |
| `RMAX_DISPLAY` (real-space guide end) | 8 | Å | `AutoStogPage.jsx` |
| $G_K(r)$ `initialYDomain` | `[2.1·level, −3.2·level]`, `level = guides.gkLowR` | data units | `AutoStogPage.jsx` |
| `S(0) FZ target` guide span | `q[0]` → `min(q[0] + 1.5, q[-1])` | Å⁻¹ | `AutoStogPage.jsx` |
| guide colour overrides | `#4c7df0` (Level), `#0c8599` (S(0) FZ target) | — | `AutoStogPage.jsx` |
| guide label rounding | `fmt(level, 5)`, `fmt(s0Target, 3)` | sig figs | `AutoStogPage.jsx` |

**Export**

| Parameter | Value | Units | Where |
| --- | --- | --- | --- |
| chart PNG scale | 2 (never overridden) | × viewBox | `svgToPngBlob` |
| chart PNG background | `#ffffff`, opaque, filled *before* `drawImage` | — | `svgToPngBlob` |
| chart PNG size | 1440 × 900 (grid), 2880 × 640 (wide) | px | derived |
| exported SVG intrinsic size | 720 × 450 or 1440 × 320 (the viewBox numbers) | CSS px | `standaloneSvg` |
| inlined style properties | 16, listed in `STYLE_PROPS` | — | `figureExport.js` |
| canvas panel export | native (`W_c·dpr`) or 3× re-render | px | `save2dPanel`, `captureModelBlob` |
| KDE slice canvas minimum | 320 × 260 | CSS px | `StructurePage.jsx` |
| slab canvas minimum | 220 × 260 | CSS px | `StructurePage.jsx` |
| canvas plane padding | 18 | CSS px | `makePlaneMapper` |
| live pixel ratio, `StructurePage` 3D | `min(devicePixelRatio, 2)` | — | `StructurePage.jsx` |
| live pixel ratio, `PcaKdePage` 3D | `devicePixelRatio` (uncapped) | — | `PcaKdePage.jsx` |
| colormap LUT | 256 entries from 5 anchors, `Uint8ClampedArray` | — | `colormaps.js` |
| ZIP compression | method 0 (store) | — | `zipArchive.js` |
| ZIP entry cap | 65535 (uint16 counts); 4 GiB (uint32 sizes, no ZIP64) | — | `zipArchive.js` |
| matplotlib figure size | 6.75 × 4.05 in (6.75 × 4.725 for STOG) | inches | `plots.py` |
| matplotlib dpi | 150 (default arg), `bbox_inches='tight'` | dpi | `plot_to_png` |
| matplotlib line style | `lw=1.0, alpha=0.65` (`alpha=1.0` for STOG) | — | `_series_plot`, `_stog_plot` |
| matplotlib STOG reference line | `ls='--', lw=0.5, color='black'`, spanning `data[0][0] … data[0][-1]` | — | `_stog_plot` |

---

### Screen vs. export: what actually changes

1. **The legend, chart title and $R_\mathrm{wp}$ chip are not in the exported figure.** They are HTML
   siblings of the `<svg>`; only the SVG is serialised. An exported PNG/SVG therefore has axes,
   grid, curves and axis labels — **but no legend and no title**. The title survives only as the
   filename.
2. **The hover tooltip is HTML and never exported.** The crosshair line and hover dots *are* inside
   the SVG, but moving the pointer to the Save button fires `onPointerLeave`, which clears them, so
   in practice exports are clean.
3. **Margins turn white in PNG; the plot rectangle was already near-white.** On screen the area
   outside the plot rectangle is transparent and shows the (dark or light) panel background; the PNG
   is composited over opaque `#ffffff`. The exported SVG leaves those margins transparent, so an SVG
   viewed on a dark background shows the chart's hard-coded dark inks against it. The reason those
   inks are legible in the app's dark theme is that `--plot-bg` is `#f8fafc` — near-white — in the
   **default/dark** theme too, not only in the light theme (`#ffffff`). That is what makes hard-coded
   dark ink and a white PNG composite mutually consistent. Every stroke colour in the chart is
   theme-independent and hard-coded: ticks `#475467`, axis labels `#1d2939`, grid
   `rgba(16,24,40,0.06)`, tick marks `rgba(16,24,40,0.35)`, frame `rgba(16,24,40,0.18)`, hover line
   `rgba(16,24,40,0.3)`, hover-dot outline `#fff`. Only the series colours, `--plot-bg` and the
   zoom-selection rect come from variables.
4. **Fonts may differ.** The UI loads Inter from Google Fonts; the exported SVG inlines
   `font-family: Inter, ui-sans-serif, …` but carries no `@font-face`, and an SVG loaded through an
   `<img>` (the PNG rasterisation path) cannot fetch external resources. Text therefore falls back to
   a locally installed font. Glyph shapes and widths change; anchor positions do not, since every
   text element carries explicit `x`, `y` and `text-anchor`, and the font sizes are absolute px
   (14.5 / 16), hence identical user-unit sizes on screen and in the export.
5. **Relative line weight changes.** `vector-effect: non-scaling-stroke` fixes the stroke at 2 px of
   the current viewport. On screen that is 2 CSS px whatever the panel width; in the exported
   standalone SVG the viewport is exactly 720 (or 1440) units wide (item 4 of Step 13a), so the same
   declaration produces a different line weight relative to the data than a chart displayed at, say,
   480 px.
6. **The zoom, the legend state and the y-autoscale are baked in.** The export renders exactly the
   current `domains` — including a y-range that was auto-derived from whatever x-window happened to
   be showing.
7. **Numerically, nothing changes.** Both export formats use the same `seriesShapes` coordinates,
   rounded to the same two decimals, from the same domains. The PNG additionally quantises to its
   1440 × 900 raster.

---

### Caveats / what this is not

- **No decimation exists — and the failure mode is a crash, not slow rendering.** Nothing is hidden
  by subsampling, but nothing protects you either. `domains` builds `visibleSeries.flatMap(s => s.x)`
  (and a second flatMap for y) on **every** zoom step and legend toggle, and `niceDomain` then calls
  `Math.min(...finite)` / `Math.max(...finite)` by argument spread. Argument-count limits (~10⁵ in
  practice, engine-dependent) make that **throw `RangeError: Maximum call stack size exceeded`** once
  the total point count across visible series passes the limit. Nothing catches it inside the memo
  and there is no error boundary on that path, so the chart unmounts into a React error rather than
  degrading. There is no point budget, no progressive rendering, and no warning. The same spread
  pattern appears in `yAxisMax`/`xAxisMax` and in `makeProjectedPlane`, but those operate on tick
  lists and polygon corners and are harmless.
- **The hover search is x-only and unbounded.** It ignores y entirely, so with several overlapping
  curves the tooltip reports every visible series' value at (its own) nearest x, not the curve you
  are pointing at. There is no "snap radius", so a cursor far from any data still produces a reading.
  Ties resolve to the lowest index, and a series whose x values are all NaN still reports index 0 —
  yielding a NaN dot and a NaN tooltip row.
- **The tooltip header x belongs to the first visible non-guide series.** On mismatched x-grids the
  other rows' y values are sampled at slightly different x than the header states.
- **Series identity is the label string.** Duplicate column headers merge into one legend chip, share
  one visibility flag, and collide as React keys. Only *empty* headers are renamed (`Series {idx}`).
- **NaN/Inf are bridged silently.** Gaps in a data file appear as straight-line interpolation with
  no visual cue.
- **Clipping is per point, and the grid is not clipped at all.** A zoomed view omits the partial
  segments that would cross the axis boundary, so curves can appear to start slightly inside the
  frame; conversely grid lines, tick marks, the frame and the axis labels are drawn outside the
  clip group.
- **Zoom is not exactly invertible by wheel.** 1.22 × 0.82 = 1.0004; use "Reset zoom" or
  double-click for an exact return. A wheel notch *during* a drag makes the released rectangle
  disagree with the preview (Step 10).
- **A y-zoom survives a file switch in Flask mode.** The fetch-path reset omits `yDomain` and
  `fullExtent`; the static/`plotData` path resets both (Step 12).
- **`isExperimental` is a substring test.** Marker styling and the experimental/calculated ordering
  are driven by the regex `/exp/i` on the series label, not by any file metadata. The Bragg
  ToF-vs-Q axis label is likewise a substring guess (`/tof|flight|time/i`) on the column header,
  with no unit conversion behind it.
- **Series colours are positional, and stable only against legend toggles.** They come from an
  8-colour palette indexed by position in the reordered list. Hiding a series does **not** recolour
  the others, because `orderedSeries` depends only on `plotData`. But a *changed series list* re-runs
  the memo and can recolour everything: the reordering step moves all `/exp/i` series to the front
  only when `paired` is true, so adding the first experimental series flips `paired`, reorders the
  list, shifts every palette index, *and* switches the measured series to hollow markers. Two
  different files can also use the same colour for different quantities. Colours are not
  accessible-contrast validated and repeat after 8 series.
- **These are inspection figures.** The chart renderer applies no error bars (error columns in the
  source CSVs are plotted as ordinary series if present as columns, and otherwise ignored), no
  weighting, no fit overlays beyond what the file contains, and no uncertainty propagation.
- **The colormaps are 5-anchor approximations** of `viridis`/`magma`/`seismic`/`reds`, not the
  published 256-entry tables — and `greys` is not an approximation of matplotlib `Greys` at all: it
  is a custom dark→light, slightly blue-tinted ramp (nearer `Greys_r`). An unknown name silently
  falls back to `viridis`. The KDE heatmap is additionally bilinearly smoothed on screen, and renders
  **empty** (not uniform) when `vmax == vmin`.
- **The ZIP writer is minimal**: stored (uncompressed), invalid zero timestamps, no ZIP64, no UTF-8
  flag (safe only because filenames are sanitised to ASCII), no directory entries, 65535-entry cap,
  and the whole archive is built in memory. Fine for a handful of figures; not a general archiver.
- **`PlotViewer.jsx` and the `GET /api/plot` matplotlib PNG path are dead code in the shipped UI.**
  They are exercised only by direct API calls, library callers and the tests. Nothing in this section
  describing the SVG chart depends on them.
- **Nothing on this page re-analyses data.** Every number drawn or reported comes from the parsers;
  the renderer's only arithmetic is the affine map, the tick algorithm and the display rounding
  described above.


## Model summary and the Detected SG symmetry finder

### What this section covers — model summary and symmetry

The Run Dashboard renders two cards side by side, both produced by
[`web_app/frontend/src/components/ModelSummary.jsx`](../../web_app/frontend/src/components/ModelSummary.jsx)
(the same component is mounted again on the KDE/3D page,
[`StructurePage.jsx`](../../web_app/frontend/src/components/StructurePage.jsx)):

- **Model information** — the unit-cell edge lengths and angles, the supercell tiling, the atom
  count per element, and the number of distinct *reference sites* per element, all read from the
  run's `.rmc6f` configuration file.
- **Detected SG** — a client-side, FINDSYM-like space-group finder run on the
  **average (reference-site) structure** of the RMC configuration, plus an interactive
  *symmetry-vs-tolerance ladder*.

**The Detected SG card is computed entirely in the browser**, in JavaScript, with no backend call,
no `spglib`, and no WASM. It consults **no external space-group database**; classification uses
in-repo tables — the 230 space groups with their standard short symbols and crystal classes
(`spaceGroupTable.js`), the Bravais centering vectors and the per-system symmetry directions
(`spaceGroupSymbol.js`), point-group orders (`symmetry.js`), and the Wyckoff positions of all 230
groups in their ITA standard setting (`wyckoffTable.js`) — all documented in Steps 10–14. There is
**no Python counterpart** anywhere in `rmc_toolkits/` or `web_app/backend/` — a repo-wide search for
`spglib`, `space_group`, or `point_group` in the Python tree returns nothing.

The **Model information** card is browser-computed too in static / local-folder mode. In
server-directory (Flask) mode its raw inputs come from Python: `Dashboard.jsx` and `StructurePage.jsx`
issue `GET /api/structure`, and the endpoint
([`web_app/backend/app.py`](../../web_app/backend/app.py) → `structure()`) returns
`latticeVectors`, `supercell`, `elementCounts`, `atomIndices`, `totalAtoms` and sampled `points`
produced by the Python parsers; only the *derived* quantities (edge lengths, angles) are computed in
JS on that path. The response carries **no `basis` field**, so in Flask mode the *Detected SG* card
is suppressed entirely (`describeSymmetry` returns `null` without a basis). The symmetry analysis
exists only on the browser-parsed path (`browserData.structureFromRmc6f`), which is used in static
mode and whenever a local run folder is selected.

#### Notation and units

| Symbol | Meaning | Units |
| --- | --- | --- |
| $L$ (rows $\mathbf L_i$) | $3\times3$ matrix whose **rows** are the RMC simulation-box (supercell) lattice vectors | Å |
| $N=(N_1,N_2,N_3)$ | supercell tiling declared in the `.rmc6f` header | dimensionless (integers) |
| $A$ (= $A_\mathrm{conv}$) | conventional-cell matrix, rows $\mathbf a_1,\mathbf a_2,\mathbf a_3$ | Å |
| $\mathbf x$ | fractional coordinate in the **conventional** cell, column vector, $x_i\in[0,1)$ | dimensionless |
| $G$ | metric tensor $G=AA^{\mathsf T}$ | Å² |
| $R$ | integer rotation part of a symmetry operation, acting on fractional columns | dimensionless |
| $\mathbf t$ | fractional translation part | dimensionless |
| $\tau$ (`tol`, `symTol`) | Cartesian atomic-position tolerance for symmetry acceptance | Å |
| $\varrho$ (`residual`) | worst-site Cartesian mapping error of one operation | Å |
| $\tau_L$ (`latticeTol`) | Cartesian tolerance on the lattice strain of a point operation (default $\tau$) | Å |

---

### Part A — Model information

`ModelSummary.jsx` returns `null` outright — **no card at all, not even Model information** — when
`structure.latticeVectors` or `structure.supercell` is missing (`if (!summary) return null`). The
card header shows the file *basename* (`structure.source.split('/').pop()`, fallback
`'Structure file'`) with the full path only as the hover `title`; `structure.source` is the parsed
file's `path`.

#### Step 1. Read the lattice and supercell from the `.rmc6f` header

**Input**: the raw text of the selected `.rmc6f` file. The choice of *which* `.rmc6f` the whole page
analyses is made by `chooseStructureFile()` in
[`browserData.js`](../../web_app/frontend/src/browserData.js), which stem-matches the model against the
run's output files:

1. `runStemFromOutputName(name)` maps each file name to a `{priority, stem}` via the first matching
   pattern: priority 0 — `^(.+)-\d{2,}\.log$`; priority 1 — `^(.+)-EXAFS-.+_[QR]_OUTPUT\.csv$`,
   `^(.+)_FT_XFQ\d+\.csv$`, `^(.+)_[FS]Q\d+\.csv$`, `^(.+)_bragg(?:_.+)?\.csv$`,
   `^(.+)_PDF(?:partials|\d+)?\.csv$`; priority 2 — `^Frac_coord_(.+)\.txt$`.
2. Those candidates are sorted by `(priority, lowercase file name)`.
3. Each is looked up in a map keyed `` `${dirname(path)}/${rmc6f stem}` `` — the match is
   **directory-scoped**, so a stem only matches an `.rmc6f` sitting in the same folder. The first hit
   wins.
4. If nothing matches, the fallback is `rmc6fFiles[0]` — the first `.rmc6f` in the **unsorted** input
   file list, i.e. directory-enumeration order, not alphabetical order.

`readCellVectors()` scans every line and takes:

- the line whose first token is `Supercell` → $N$ = the **last three** whitespace tokens, parsed as
  numbers (`parts.slice(-3).map(Number)`);
- the line whose first token is `Lattice` → the **next three lines**, each split on whitespace and
  parsed as numbers, become the rows of $L$ (Å).

If either is missing the parse throws `Missing lattice or supercell metadata` and no summary is
shown. Note that the last `Supercell`/`Lattice` occurrence in the file wins (the loop overwrites).

**No numeric validation.** `readCellVectors` checks only that the two *markers* exist. The three
lattice rows are read as `row.trim().split(/\s+/).map(Number)` with no `filter(Boolean)`, no length
check and no finiteness check, and the supercell tokens are `parts.slice(-3).map(Number)` with no
positive-integer check. A blank or short line after `Lattice` therefore yields `NaN` (or `undefined`)
entries that propagate silently into $G$. Because `Math.abs(NaN - NaN) > eps` is `false` and the
`|| 1` guard turns a `NaN` trace into a scale of 1, `latticePointOps` then accepts **all 6960**
unimodular $\{-1,0,1\}$ patterns instead of throwing (verified by running the code). Downstream the
damage is contained rather than silent-but-wrong: `cartDist` returns `NaN` for every pair, so
`mappingResidual` rejects every candidate and the card degrades to `P1` / No. 1 / **0 operations**
with an empty ladder, while the Model information card prints `NaN` cell edges. Nothing is raised.

**Code**: `browserData.js` → `readCellVectors()`; Python equivalent
`rmc_toolkits/parsers.py` → `read_cell_vectors()` uses the identical rule (`parts[-3:]` and the
three following lines) and the two agree exactly.

#### Step 2. Per-atom parsing, element counts, and reference sites

`structureFromRmc6f(file, maxPoints = 100)` walks the lines after the `Atoms:` marker and hands each
whitespace-split line to `parseAtomLine()` in
[`rmc6f.js`](../../web_app/frontend/src/rmc6f.js). That function indexes **from the end** of the line so
that any number of label columns between the element and the coordinates is tolerated:

- $\ge 9$ fields ("full" format): the last four fields are the reference number and the three cell
  indices $(c_1,c_2,c_3)$; the three before them are the fractional box coordinates.
- 5–6 fields ("coords-only", oldest format): the last three fields are the coordinates;
  `referenceNumber` and `cellIndices` come back `null`.
- 7–8 fields, or non-finite numbers → `null` (line skipped).

From this:

$$\texttt{totalAtoms} = \#\{\text{parsed atom lines}\},\qquad
\texttt{elementCounts}[e] = \#\{\text{atoms with element } e\}$$

$$\texttt{atomIndices}[e] = \{\, \text{distinct reference numbers of element } e \,\}\ \text{(sorted ascending)}$$

`elementCounts` counts atoms in the **whole simulation box**; `atomIndices[e].length` is the number
of distinct **reference sites**, i.e. sites of the conventional cell — that is the "N sites"
sub-label under each element in the card. Note that `atomIndices[element].add(referenceNumber)` is
keyed by *the atom's own element token*, so a reference number carried by atoms of more than one
species (a chemically mixed / swapped site) lands in more than one element's set; see Step 5 for the
consequences. The card's "Total atoms" line likewise shows `totalAtoms` with
$\sum_e |\texttt{atomIndices}[e]|$ as its sub-label (`ModelSummary.jsx` → `summary.referenceSites`).
The per-element and total "N sites" sub-labels render only when the count is $>0$.

**Display order.** `structure.elements` is `Object.keys(counts).sort()` (plain lexicographic, by code
point); the card's element rows are sorted independently by `a.localeCompare(b)` over the element
token (locale-aware). The two orders can differ under a non-default locale.

The `points` array (the 3D/slab preview sample) is a stride subsample followed by a hard cap:

```
stride = max(1, ceil(nAtoms / maxPoints))
points = atoms.filter((_, i) => i % stride === 0).slice(0, maxPoints)
```

The `.slice()` is defensive: since `ceil(n / ceil(n/m)) ≤ m` for all $n$, the stride filter can never
exceed `maxPoints` and the slice never actually truncates (checked by enumeration for
$n \le 5000$, $m = 100$). Each sampled atom's box coordinate is folded straight into one unit cell as
$((x_i N_i)\bmod 1 + 1)\bmod 1$; subtracting the atom's cell index first would be equivalent (it only
removes an integer before the modulo), which is why this fold **also works for the old coords-only
format** that carries no cell index. **The subsample never touches the summary numbers or the
symmetry finder** — those use all atoms / all reference sites.

`maxPoints` is **not** a single app-wide constant:

- **Run Dashboard** — `100`, both for the Flask query (`Dashboard.jsx`, `params: { dir, maxPoints: 100 }`)
  and for the local worker (`Dashboard.jsx` posts `maxPoints: 100` to
  `workers/localStructureWorker.js`).
- **KDE/3D page** — `STRUCTURE_MAX_POINTS = 1000000` (`StructurePage.jsx`), used for both the worker
  post and the Flask query, so `points` there is effectively unsampled (stride 1 for any realistic
  configuration).

`ModelSummary` behaves identically on both pages because the summary and symmetry paths ignore
`points` entirely — only `sampleStride` / `sampledAtoms` / `points` differ. The Flask endpoint
additionally clamps the request to $[100, 10^6]$ (`app.py`, `MAX_STRUCTURE_POINTS = 1_000_000`) and
samples *per reference site* (`_sample_atoms_by_site()`) rather than by a flat stride.

The `.rmc6f` header's declared `Number of atoms:` is parsed for nothing and is **not** validated
against the number of atom lines actually read (listed as a known issue in
[AGENTS.md](../../AGENTS.md)).

**Code**: `browserData.js` → `structureFromRmc6f()`; `rmc6f.js` → `parseAtomLine()`.

**Python counterpart** (`rmc_toolkits/parsers.py`), used only in Flask mode — it differs from the JS
parser in two ways, not one:

- `iter_rmc6f_atoms()` uses the same index-from-the-end rule **for the ≥ 9-field full format only**:
  it hard-rejects shorter lines (`n = len(parts); if n < 9: continue`). The 5–6-field "coords-only"
  form that `parseAtomLine` tolerates yields **zero atoms** in Flask mode, so an old file that gives
  a full Model information card in browser mode gives an empty one through `/api/structure`.
- It **capitalizes** the element token (`parts[1].capitalize()`, which also lowercases the tail:
  `SE → Se`) while the JavaScript parser keeps it verbatim. Element identity is compared by exact
  string equality downstream, so a file mixing `SE` and `Se` is two species in browser mode and one
  in Flask mode. Worse, the sibling function `read_atom_indices()` — which produces the `atomIndices`
  in the same response — does **not** capitalize and accepts `len(parts) >= 5` with `int(parts[-4])`.
  A file written with upper-case tokens therefore returns `elementCounts` keyed `Se` and
  `atomIndices` keyed `SE`, and the card shows the element row with **no** "N sites" sub-label while
  the "Total atoms" sub-label still counts those sites.

#### Step 2b. Run counters from the header (`readMovesMetadata`)

`structureFromRmc6f` also calls `readMovesMetadata(file.text)` and returns the result as
`structure.moves`. It slices the header —
`text.slice(0, text.indexOf('Atoms:') > 0 ? text.indexOf('Atoms:') : 4000)`, i.e. a 4000-character
fallback when the marker is absent or at index 0 — and applies four regexes:
`Number of moves generated:`, `… tried:`, `… accepted:` (each `([\d.]+)`) and
`Accumulated time \(s\)[^:]*:\s*([\d.]+)`. Each field is `Number(match[1])` or `null`, and the whole
object collapses to `null` unless at least one value is finite.

It is **not rendered on either card**. It feeds the AI-assistant run context
(`runContext.js`: acceptance ratio $=$ accepted/tried, accepted moves per atom, accumulated time in
hours). There is no Python equivalent on the `/api/structure` path, so `structure.moves` is absent in
Flask mode.

**Code**: `browserData.js` → `readMovesMetadata()`, called from `structureFromRmc6f()`.

#### Step 3. Conventional-cell edge lengths

The `.rmc6f` lattice vectors describe the *box*; the conventional cell is obtained by dividing each
row by its supercell multiplicity:

$$\mathbf a_i = \frac{\mathbf L_i}{\max(N_i,1)},\qquad
a_i = |\mathbf a_i| = \frac{\sqrt{\sum_k L_{ik}^2}}{\max(N_i,1)} \ \ [\text{Å}]$$

Displayed as `a × b × c` with **at most** 3 fractional digits — `formatNumber()` is
`Number(value).toLocaleString(undefined, { maximumFractionDigits: digits })`, a *maximum*, not a
fixed precision, so trailing zeros are dropped (a 5 Å edge renders as `5`, not `5.000`). The same
helper renders every figure on the card: `digits = 3` for cell edges, `1` for angles, `0` for counts
and supercell multiplicities (counts are *rounded* for display, never validated as integers). Passing
`undefined` as the locale means the **viewer's** locale decides the decimal separator and digit
grouping — a French locale renders 10.532 Å as `10,532`.

**Code**: `ModelSummary.jsx` → `vectorLength()` and the `summary` memo; the same expression is
duplicated verbatim in
[`llm/context/runContext.js`](../../web_app/frontend/src/llm/context/runContext.js) →
`structureContext()` (deliberate duplication — the `src/llm/` module is not allowed to import from
the host app, per [AGENTS.md](../../AGENTS.md)), and again in `browserData.js` as `cellEdgeA`.

#### Step 4. Cell angles

$$\alpha = \angle(\mathbf A_{\mathrm{box},2},\mathbf A_{\mathrm{box},3}),\quad
\beta  = \angle(\mathbf A_{\mathrm{box},1},\mathbf A_{\mathrm{box},3}),\quad
\gamma = \angle(\mathbf A_{\mathrm{box},1},\mathbf A_{\mathrm{box},2})$$

with

$$\angle(\mathbf u,\mathbf v) = \frac{180}{\pi}\,
\arccos\!\left(\mathrm{clamp}\!\left[\frac{\mathbf u\cdot\mathbf v}{\max(|\mathbf u||\mathbf v|,\,10^{-12})},\,-1,\,1\right]\right)$$

The angles are computed from the **box** vectors, not the conventional-cell vectors. This is exact,
not an approximation: dividing each row by a positive scalar $N_i$ leaves all pairwise angles
unchanged. The $10^{-12}$ floor on the denominator and the $[-1,1]$ clamp guard against a
degenerate/zero row and against `acos` domain error from round-off. Displayed in degrees with **at
most** 1 fractional digit (`maximumFractionDigits: 1`), so a right angle renders as `90°`, not
`90.0°`.

**Code**: `ModelSummary.jsx` → `angleBetween()`; duplicated in `runContext.js` → `angleBetween()`
(identical, except `runContext` rounds to 4 significant digits for the LLM context).

#### Step 5. The average-structure basis (what the symmetry finder actually sees)

This is the step that turns a disordered RMC snapshot into an "average structure". Every atom
carries a **reference number** identifying which site of the conventional cell it is a copy of. The
loop **returns early for any atom with `referenceNumber === null`** — before `atomIndices` and before
the `rnAcc` accumulator — so the old coords-only format contributes to `totalAtoms` and
`elementCounts` and to `points`, but produces an *empty basis* and hence no symmetry analysis at all.

For each remaining atom, its position *within one conventional cell* is the box fraction taken
modulo 1 after multiplying by the supercell multiplicity:

$$w_i = \big((x_i \, N_i) \bmod 1 + 1\big) \bmod 1,\qquad i=1,2,3$$

where $x_i$ is the fractional box coordinate (the doubled modulo is the JavaScript idiom for a true
modulus, since `%` is a sign-following remainder).

**Supercell guard inconsistency.** This fold uses the **raw** `supercell[i]`, whereas every other use
of the multiplicity divides by `Math.max(supercell[i], 1)` — the card's cell lengths
(`ModelSummary.jsx`), `cellEdgeA` in `browserData.js`, and `conventionalCell()` in
`symmetryModel.js`. Since `readCellVectors` never validates that the three `Supercell` tokens are
positive integers, a header declaring `0` (or a non-integer) yields a *guarded*, finite conventional
edge on that axis while collapsing every atom's $w_i$ to 0 — a one-site basis and a spurious
high-symmetry answer, with no error raised.

The per-site representative is the **circular mean** of $w$ over all box copies of that reference
number, computed independently on each axis:

$$\bar w_i = \frac{1}{2\pi}\,\operatorname{atan2}\!\Big(\textstyle\sum_c \sin 2\pi w_i^{(c)},\ \sum_c \cos 2\pi w_i^{(c)}\Big)
\ \ \text{wrapped into } [0,1)$$

A circular (not arithmetic) mean is required so that a site straddling the cell boundary
($w\approx0 \equiv w\approx1$) averages to the true position rather than to $0.5$. This is asserted
by `structureFromRmc6f site displacement (dispA) › handles a boundary-wrapping site (mean at 0 ≡ 1)`
in [`__tests__/browserData.test.js`](../../web_app/frontend/src/__tests__/browserData.test.js).

The same accumulators give the per-site spread for free. With resultant length
$\bar R_i = \big|\sum_c(\cos,\sin)\big| / N_c$ over the $N_c$ copies, the circular standard
deviation in cell fractions is $\sigma_i^{\mathrm{frac}} = \sqrt{-2\ln \bar R_i}\,/\,2\pi$ (with $\bar R_i$ floored
at $10^{-6}$, and $\sigma_i^{\mathrm{frac}}$ taken as 0 when $\bar R_i\ge1$, i.e. a single copy or zero spread), and

$$u_s \;\equiv\; \texttt{dispA} = \sqrt{\sum_{i=1}^{3}\big(\sigma_i^{\mathrm{frac}} \, a_i\big)^2}\ \ [\text{Å}]$$

This rms displacement is **not shown on the card**; it is consumed by the AI-assistant context
(`runContext.js` → `symmetryContext()` aggregates `mean_disp_A` / `max_disp_A` per Wyckoff orbit).
Two approximations are worth naming: (i) $\sigma_i^{\mathrm{frac}} a_i$ multiplies a fractional spread by an edge
*length*, which ignores the metric cross-terms and is therefore exact only for orthogonal axes;
(ii) the circular-std formula is the von-Mises/wrapped-normal relation, exact only for a wrapped
Gaussian. The unit test pins the value for a two-copy $\pm0.02$-fraction case on a 10 Å edge at
$0.2003$ Å.

Sites are keyed by reference number and emitted **sorted by reference number**, as
`{ el, referenceNumber, frac, dispA }`. Per the early return above, an oldest-format file yields an
empty basis, hence the Model information card and no Detected SG card.

**A site's species is fixed by the first atom that carries its reference number.** The accumulator is
created as `if (!acc) { acc = { element, n: 0, sc, ss }; }` and `acc.element` is never revisited, so
every later copy of that reference number is accumulated under whatever species happened to appear
first in the file, regardless of its own element token. Combined with `atomIndices` being keyed by the
atom's own element (Step 2), a configuration in which one reference number carries more than one
species — e.g. after RMCProfile swap moves, which change an atom's species while leaving the
reference-number column alone — has two visible consequences:

1. $\sum_e |\texttt{atomIndices}[e]|$, the card's "sites" sub-label, **counts that reference number
   once per species** and can therefore exceed `basis.length`, the number of sites the symmetry
   finder actually uses.
2. The finder is handed a **fictitious ordered decoration**: a chemically mixed site is presented as
   a single pure species, and the element-equality matching of Step 9 (property 2) treats it as such.

**Code**: `browserData.js` → `structureFromRmc6f()`, the `rnAcc` accumulator and the `basis` map.

#### What Model information does *not* compute

There is **no cell volume and no number density** on this page. `ModelSummary.jsx` computes only
lengths, angles, counts and site counts; `runContext.js` mirrors the same set.

The only determinant in this subsystem is `det3()` in
[`symmetry.js`](../../web_app/frontend/src/symmetry.js), and it has exactly two call sites — neither of
which is a volume:

- `latticePointOps()` — keep only unimodular candidates, $|\det R| = 1$ (Step 7);
- `classifyRotation()` — split proper ($\det = +1$) from improper / roto-inversion ($\det = -1$)
  operations, which is what drives the rotation-type table of Step 10(d).

The metric tensor $G$ is likewise built only for the automorphism test; $V = \sqrt{\det G}$ is never
evaluated.

Number density $\rho_0$ (Å⁻³) is computed elsewhere in the app, for a different purpose — the Auto
StoG scaling pipeline — and never from the `.rmc6f` model. The conversion from a chemical formula and
a mass density lives in [`rmc_toolkits/scattering.py`](../../rmc_toolkits/scattering.py) →
`number_density_from_mass_density()`; the `NUMBER_DENSITY` header of a STOG `.dat` input is parsed in
[`rmc_toolkits/parsers.py`](../../rmc_toolkits/parsers.py) → `read_dat_header()` (→ `number_density`).
The scaling engines ([`rmc_toolkits/scaling.py`](../../rmc_toolkits/scaling.py),
[`workers/autoScale.js`](../../web_app/frontend/src/workers/autoScale.js)) only *consume* `config.rho0`
(and `estimate_rho0()` iterates it). That path is documented in its own section.

---

### Part B — The Detected SG symmetry finder

Part B is pure JavaScript (no React, no I/O) in six modules under `web_app/frontend/src/`:

| Module | Role |
| --- | --- |
| [`symmetry.js`](../../web_app/frontend/src/symmetry.js) | lattice rotations, operation search and refinement, closure walk, ladder, point group, orbits (Steps 7–13), the operation estimate |
| [`spaceGroupSymbol.js`](../../web_app/frontend/src/spaceGroupSymbol.js) | screw/glide analysis, centering, the standard-setting search, the H–M symbol, the lower-bound check (Step 10) |
| [`spaceGroupTable.js`](../../web_app/frontend/src/spaceGroupTable.js) | the 230 groups: number, standard short symbol, crystal class, pre-2002 `e`-glide spellings |
| [`wyckoff.js`](../../web_app/frontend/src/wyckoff.js) + [`wyckoffTable.js`](../../web_app/frontend/src/wyckoffTable.js) | Wyckoff positions of the 230 groups and the letter assignment (Step 14) |
| [`symmetryModel.js`](../../web_app/frontend/src/symmetryModel.js) | structure → finder glue: the cell, the basis-size cap and operation budget, orbits and letters in the naming cell |

Tests: `web_app/frontend/src/__tests__/symmetry*.test.js` and `wyckoff.test.js` (the 230-group
fixtures in `__tests__/fixtures/spaceGroups.js`, well-known structures in
`__tests__/fixtures/symmetryStructures.js`).

**Convention**: fractional coordinates are **column vectors** and operations act as
$\mathbf x' = R\,\mathbf x + \mathbf t$. The lattice matrix $A$ stores lattice vectors as **rows**,
so a Cartesian position is $\mathbf r = A^{\mathsf T}\mathbf x$ (Å) and
$|\mathbf r|^2 = \mathbf x^{\mathsf T} G\, \mathbf x$ with $G = AA^{\mathsf T}$.

#### Step 6. Build the conventional cell

$$A_{ij} = \frac{L_{ij}}{\max(N_i,1)}$$

i.e. row $i$ of the box lattice divided by the $i$-th supercell multiplicity. This is a *declared*
conventional cell — it is whatever cell the `.rmc6f` header says the box tiles. If the RMC box was
built on a primitive cell, or on a non-standard setting, the finder analyses that cell.

**Code**: `symmetryModel.js` → `conventionalCell()`.

#### Step 7. Lattice point operations: integer automorphisms of the metric

**Input**: $A$ (Å), Cartesian lattice tolerance $\tau_L$ (`latticeTol`, Å; default: the pass
tolerance $\tau$).

The metric tensor is built explicitly,

$$G_{ij} = \mathbf a_i\cdot\mathbf a_j = \sum_{k}A_{ik}A_{jk} \qquad (G = AA^{\mathsf T},\ \text{Å}^2)$$

A linear map $\mathbf x \mapsto R\mathbf x$ preserves all lattice distances iff
$(R\mathbf x)^{\mathsf T}G(R\mathbf x) = \mathbf x^{\mathsf T}G\mathbf x$ for all $\mathbf x$, i.e.

$$\boxed{\,R^{\mathsf T} G R = G\,}$$

and it maps the lattice onto itself iff $R$ is an integer matrix with $|\det R| = 1$. The code
enumerates candidates **exhaustively, in a reduced basis** of the lattice:

- **Reduce** (`reduceBasis()`): an integer matrix $M$ ($\det M=+1$, rows = new basis vectors in the
  given basis) such that $A_r = MA$ is a reduced basis of the same lattice — each vector is
  repeatedly shortened by integer multiples of the others (pairwise size reduction, then
  $\mathbf b_k \pm \mathbf b_i \pm \mathbf b_j$), shortest first, until nothing shortens (the
  Minkowski conditions in three dimensions). A conventional cell is already reduced, so this only
  matters for an oblique cell.
- Every $3\times3$ matrix with entries drawn from $\{-1,0,1\}$: $3^9 = 19\,683$ patterns, iterated
  as a base-3 counter over `code = 0 … 19682`, taken as $R_r$ in the reduced basis.
- Keep those with $\det R_r \in \{+1,-1\}$ — **6960** of the 19 683 (verified by direct enumeration).
- Keep those whose **Cartesian lattice strain** (on $A_r$) is at most $\tau_L$, and carry each back to
  the given basis, $R = M^{\mathsf T}R_rM^{-\mathsf T}$ (an integer matrix, since $M$ is unimodular;
  fractional columns transform as $\mathbf x_r = M^{-\mathsf T}\mathbf x$).

**Lattice strain** (`latticeStrain()`). Used as if it were an isometry, $R$ acts on Cartesian
vectors as $M = A^{\mathsf T}RA^{-\mathsf T}$, whose Green strain is

$$E=\tfrac12\big(M^{\mathsf T}M-I\big)=\tfrac12\,A^{-1}DA^{-\mathsf T},\qquad D = R^{\mathsf T}GR-G\ [\text{Å}^2].$$

The cell edge $\mathbf a_i = A^{\mathsf T}\mathbf e_i$ is displaced by $|E\mathbf a_i| = \tfrac12|A^{-1}D\mathbf e_i|$,
and the strain of $R$ is the largest of the three (evaluated with $A = A_r$, the reduced cell's
edges — the same three edges as the given cell's for a conventional cell),

$$\varrho_L(R) = \max_i\ \tfrac12\,\big|A^{-1}D\,\mathbf e_i\big|\quad[\text{Å}].$$

It is 0 for an exact lattice symmetry and is on the **same Å scale as an atomic-position residual**
(Step 9). It is not an absolute threshold on $D$ scaled by the mean squared edge, so a long $c$ axis
does not loosen the test on $a$ and $b$ (a 2.8 % $a/b$ splitting in a cell with $c = 20$ Å is a
0.14 Å strain whatever $c$ is). A `NaN` lattice gives `NaN` strain and every $R$ is rejected.

Each surviving $R$ carries its strain into the operation's residual: $\varrho = \max(\varrho_\text{atoms},
\varrho_L)$ (Step 9). So a homogeneous strain of a higher-symmetry cell whose atoms still sit at
parent-compatible *fractional* coordinates is no longer invisible: a perovskite with $c/a = 1.004$
(0.016 Å) reads `P4/mmm` below 0.016 Å and `Pm-3m` above it, and the cubic rung of the ladder starts
at ≈ 0.016 Å. Because the accepted set of lattice operations need not compose (two operations that
each strain the cell by $\le\tau_L$ can compose to one that strains it by more), closure is not
assumed here: it is enforced on the final operation set (Step 11).

Verified by enumeration: a cubic metric ($a=10$ Å) yields exactly 48 operations; a primitive
hexagonal metric ($a=b=5$ Å, $c=8$ Å, $\gamma=120°$) yields exactly 24. The full holohedry order
(48, 24, 16, 12, 8, 4) also comes out for cubic P/F/I, tetragonal P/I, hexagonal, rhombohedral,
orthorhombic P/C/I/F and monoclinic P lattices each re-described by 59 random unimodular basis
changes (products of six random shears; $\tau_L = 10^{-6}$ Å), and the test suite pins it on four
strongly oblique cells (`symmetryObliqueCells.test.js`).

**Why $\{-1,0,1\}$ suffices — in a reduced basis.** In a Minkowski-reduced basis every lattice
automorphism has entries in $\{-1,0,1\}$ (the enumeration spglib uses after Delaunay reduction); a
conventional cubic, tetragonal, orthorhombic, monoclinic, hexagonal (six-fold
$\big[\begin{smallmatrix}1&-1&0\\ 1&0&0\\ 0&0&1\end{smallmatrix}\big]$) or rhombohedral cell is already such a
basis. In an **oblique** description of the same lattice — a tetragonal crystal on a 60° cell
($\mathbf a' = \mathbf a+\mathbf b$), a sheared cell — some automorphisms need entries of magnitude 2 or
more in the given basis. Before 1.0 the scan ran in the given basis, missed them, and named the
subgroup that was left with its own number (rutile on a 60° cell: `Cmmm` No. 65; a Pnma perovskite
on a sheared cell: `P-1` No. 2). Scanning the reduced basis and carrying the matrices back finds
every lattice rotation of the given cell, whatever its shape.

**What it still excludes.** Any symmetry of the structure that does not map the **declared cell's
lattice** $\mathbb Z^3$ onto itself: the finder works in the declared cell, and $A$ is divided by $N$
before the search. When the declared cell is a supercell of the crystal's own translation lattice
$\mathbb Z^3+T$ (Step 10c) — a perovskite on a $2\times2\times1$ or $\sqrt2\times\sqrt2\times2$ cell,
rocksalt on a $1\times1\times2$ cell — rotations of $\mathbb Z^3+T$ that do not preserve $\mathbb Z^3$
(the cubic 3-folds there) are never tested. The naming step detects this and reports the result
as a lower bound (`≥ <symbol>`, Step 10f).

**Output**: an array of integer $3\times3$ matrices (row-major), typically 2–48 entries.

**Code**: `symmetry.js` → `metricTensor()`, `det3()`, `inv3()`, `conjugate()`, `latticeStrain()`,
`latticeCandidates()`, `latticePointOps()` (the exported wrapper, default $\tau_L = 0.01$ Å);
`spaceGroupSymbol.js` → `reduceBasis()` (re-exported by `symmetry.js`; the naming step uses it too). The
finder passes $\tau_L = \min(\texttt{latticeTol}, \tau)$ with `latticeTol` defaulting to $\tau$.

#### Step 8. Candidate translations

**Input**: the point-operation list, the basis $\{(e_s, \mathbf x_s)\}$ (element + fractional
position), the Cartesian tolerance $\tau$ Å.

The basis is bucketed by element (`byEl`). A reference atom $a_0$ is chosen from the **rarest
element** (fewest sites in the cell) — a speed choice, since it minimises the number of candidate
partners — and both choices are made from the sites themselves, not from their order:

```
refEl   = the element with fewest sites; ties → the smaller element name
refAtom = the refEl site with the lexicographically smallest (x, y, z), each wrapped to [0, 1)
```

The reference atom fixes the *seed* translations, and a seed is refined by least squares (Step 9)
to the nearest local optimum: at a tight pairing radius every seed of one operation reaches the
same one, but at the ladder's 2 Å (Step 9) a seed can settle elsewhere. So $a_0$ must not depend on
the order of the basis. Before 1.0 it was the first site of the rarest element in basis order (ties
by first appearance), and a shuffled basis could move the ladder's bricks: a noisy $P6_3/mmc$
(σ = 0.03 Å) had its $P6_3/mmc 	o P6/mmm$ boundary at 0.794 Å in one site order and 0.767 Å in
another. With the order-free choice the operations, residuals and bricks are the same for any order
of the basis (up to floating-point summation order in the mean offset).

For each rotation $R$ and each same-element candidate partner $\mathbf x_b \in$ `byEl[refEl]`:

$$\mathbf t = \operatorname{frac}\!\big(\mathbf x_b - R\,\mathbf x_{a_0}\big),\qquad
\operatorname{frac}(u)=u-\lfloor u\rfloor \in [0,1)$$

This is exhaustive over the possible images of $a_0$: any genuine operation $\{R|\mathbf t\}$ must
send $a_0$ to *some* same-element site, so its translation part is one of these candidates (up to
the mapping error at that site). The candidate is only a **seed**: it carries the displacement noise
of $a_0$ and of its partner, so Step 9 refines it by least squares over every site before the
operation is accepted or its residual recorded. Seeds of the same operation from nearby
starting points refine to the same translation; which atom is $a_0$ is fixed independently of the
basis order (above), so even a seed that settles in another local optimum does so in every order.

**Deduplication**: every seed is refined (Step 9). Refined translations of one $R$ closer than
$\tau$ Å to each other — $\tau$ being the **pass** tolerance, 1 Å in the card's ladder and headline
(Step 12) — (minimum-image Cartesian distance) are one operation at this resolution, and
the **best-fitting** one is kept: they are taken in order of residual (ties by translation) and
each is dropped if it lies within $\tau$ of one already kept. Before 1.0 the first one found, in
seed order, was kept and seeds near it were skipped; at the ladder's loose 1 Å a poor
near-duplicate could then shadow the true operation for good (a noisy $P4_322$: a 4-fold at 0.8 Å
residual kept, 0.7 Å from the $4_3$ at 0.2 Å, so the ladder read `P222_1` up to 0.96 Å and then a
set that was not a group, while the headline read `P4_322`).

**Code**: `symmetry.js` → `detectOperations()`; `applyR()`.

#### Step 9. Acceptance test and residual

For a candidate $\{R|\mathbf t\}$, every basis site $s$ is mapped,

$$\mathbf y_s = R\,\mathbf x_s + \mathbf t ,$$

and paired with its **nearest same-element** site $o(s)$ under a minimum-image Cartesian distance:

$$\boldsymbol\delta_s = (\mathbf x_{o} - \mathbf y_s) - \operatorname{round}(\mathbf x_{o} - \mathbf y_s),\qquad
d_s = \big\| A^{\mathsf T}\boldsymbol\delta_s \big\|_2$$

(the component-wise nearest-integer wrap; result in Å).

**Least-squares translation** (`refineOperation()`). With the pairing fixed, the translation that
minimises $\sum_s d_s^2$ is

$$\mathbf t \leftarrow \mathbf t + \frac1N\sum_s \boldsymbol\delta_s ,$$

because every pair shares the metric $G$. The sites are re-paired at the new $\mathbf t$ and the
step repeated until the mean offset vanishes (at most four passes). Pairing is within $2\tau$ of
the pass (the seed can be off by the noise of the two atoms that defined it; 2 Å in the card's
1 Å pass), so the pass tolerance, not only the acceptance threshold, decides which local optimum a
seed refines to; the refined operation is **accepted** iff every site then has a partner and

$$\varrho(R,\mathbf t) \;=\; \max_{s}\ d_s \;\le\; \tau \quad [\text{Å}]$$

i.e. the *worst-site* nearest-image error at the least-squares translation — an
$L_\infty$-over-sites, $L_2$-in-space measure. It is independent of the basis order (with $a_0$
chosen as above) and never
larger than the error read off a single reference pair (on the bundled demo the full-group residual
drops from 0.038 Å to 0.031 Å). `!(ϱ ≤ τ)` also rejects `NaN`.

Three properties matter for interpretation:

1. The test is a **one-way covering test**, not a bijection test: two distinct sites may map onto the
   same partner and the operation is still accepted. No permutation/injectivity check is performed.
2. Occupancy/element matching is by **exact string equality** of the element token.
3. The minimum image uses component-wise rounding of the fractional difference. For orthogonal cells
   this is the true minimum image; for oblique cells it can return a *larger* distance than the true
   minimum, so residuals are, if anything, over-estimated and symmetry under-detected there.

**Strict vs. non-strict comparisons.** The code is not uniform about the boundary, and the
difference is only visible at exactly the tolerance, so it is worth stating per test:

| Test | Code | Boundary |
| --- | --- | --- |
| operation acceptance (Step 9) | `if (!m \|\| !(m.worst <= tol)) return null` | **non-strict**: $\varrho = \tau$ is accepted |
| translation dedup (Step 8) | `cartDist(u, t, A) < tol` | strict |
| orbit union + stabiliser (Step 14) | `bestD = tol; if (d < bestD)`, `cartDist(...) < tol` | strict |
| Wyckoff coordinate form (Step 14) | `Math.abs(fitted - p[i]) > tol` rejects | **non-strict**: a deviation of exactly `tolFrac` fits |
| translation snap (Step 10c) | `<= tol` with `tol = 0.02` | non-strict |
| residual threshold filters (Steps 12–13) | `o.residual <= r + 1e-9` | non-strict, with $10^{-9}$ slack |

`findSpaceGroupOps` returns `{ ops: [{R, t, residual}], order = ops.length, maxResidual }` merged
with the classification of Step 10. (`maxResidual` here is the max over *accepted* ops.)

**Code**: `symmetry.js` → `matchImages()`, `refineOperation()`, `offset()`, `findSpaceGroupOps()`.

#### Step 10. Classification: centering, point group, Hermann–Mauguin symbol

`classifyOperations(ops, tolFrac)` reduces an operation set to labels.

**(a) Distinct rotation parts.** Operations are keyed by the flattened matrix string;
`nPoint = ` number of distinct $R$. This set *is* the point group of the detected space group.

**(b) Pure translations.** Operations whose rotation is the identity are collected. Their count,
`nTrans`, uses the key `Math.round(wrap01(t_i) * 1000) % 1000` per component — a fixed $10^{-3}$
fractional grid, folded mod 1 so that a translation recovered as 0.9997 and an exact 0 count once.

`tolFrac` $=\tau/\overline a$ with $\overline a = \tfrac13(|\mathbf a_1|+|\mathbf a_2|+|\mathbf a_3|)$
(`meanEdge()`) is used only by the closure check of a set handed in directly (Step 11).

**(c) Centering letter** (`centeringOfOps()` → `bravaisCentering()`, `spaceGroupSymbol.js`). The
pure translations are snapped to exact fractions as a **group** (`pureTranslations()`): the $n$
pure translations (zero included) form a finite group, so each has an order $d$ dividing $n$
($d\boldsymbol\tau$ is a lattice vector), and each is snapped to $\operatorname{round}(d\boldsymbol\tau)/d$ with
the smallest such $d$ for which every component lies within 0.02 of the $1/d$ grid. $d$ is capped
at $1/(2\cdot0.02) = 25$, where neighbouring grid points are $2\cdot0.02$ apart and a snap would be a
guess. This covers every centering and any supercell fraction up to $1/25$ (fifths, sevenths, …)
with the loose tolerance a refined translation needs at a loose $\tau$ with few sites (0.0075 off
on the noisy R-3m fixture at $\tau = 0.6$ Å). A set that does not snap, or in which two
translations snap together, gives no letter. Before 1.0 each component was snapped to a fixed
$1/24$ grid, which holds no fifths: CsCl in a $5\times5\times5$ cell went unnamed. The **whole
set**, zero included, must equal one Bravais centering exactly:

| Letter | Translations besides 0 (fractional) |
| --- | --- |
| `P` | none |
| `A` / `B` / `C` | $(0,\tfrac12,\tfrac12)$ / $(\tfrac12,0,\tfrac12)$ / $(\tfrac12,\tfrac12,0)$ |
| `I` | $(\tfrac12,\tfrac12,\tfrac12)$ |
| `F` | $(0,\tfrac12,\tfrac12)$, $(\tfrac12,0,\tfrac12)$, $(\tfrac12,\tfrac12,0)$ |
| `R` | obverse $(\tfrac23,\tfrac13,\tfrac13)$, $(\tfrac13,\tfrac23,\tfrac23)$; reverse $(\tfrac13,\tfrac23,\tfrac13)$, $(\tfrac23,\tfrac13,\tfrac23)$ |

A translation set that is not exactly one of these — the finer lattice of a **supercell** of the
true cell (a perovskite in a $2\times2\times2$ cell has all eight $(i/2, j/2, k/2)$, which contain the F
vectors but are not an F lattice) — gives `centering = null`. Before 1.0 the letter came from the
mere presence of the vectors, so such cells were read as F, I or C. The letter is the centering of
the **given** cell; the symbol is named in whatever cell Step 10(f) finds, with that cell's letter.

**(d) Rotation type from determinant and trace.** Both are similarity invariants, so this is
basis-independent. With $t = \operatorname{tr}R$:

| $\det R$ | $t=3$ | $t=2$ | $t=1$ | $t=0$ | $t=-1$ | $t=-2$ | $t=-3$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| $+1$ | `1` | `6` | `4` | `3` | `2` | (falls to `2`) | (falls to `2`) |
| $-1$ | (falls to `-6`) | (falls to `-6`) | `m` | `-3` | `-4` | `-6` | `-1` |

This is the standard $t = 1+2\cos\theta$ relation for proper rotations
($\theta = 0, 60°, 90°, 120°, 180°$ → $t = 3, 2, 1, 0, -1$) and its roto-inversion counterpart
($\bar1, m, \bar4, \bar3, \bar6$ → $t = -3, 1, -1, 0, -2$). Both branches are chains of ternaries
whose final `else` is a catch-all (`2` and `-6` respectively), so an unexpected trace would be
silently mislabelled; for a genuine isometry this cannot occur.

**(e) Point-group symbol.** `pointGroupOf(rotations)` builds a histogram $h$ over the ten types and
uses `order = rotations.length`, `inv = h['-1']>0`, `nm = h['m']`:

```
if h['3'] >= 8            → cubic : order 12 → 23 ; order 48 → m-3m ;
                                    else inv ? m-3 : (h['4']>0 ? 432 : -43m)
elif h['6']>0 or h['-6']>0 → hex   : order 24 → 6/mmm ; order 12 → inv ? 6/m :
                                    (h['6']>0 ? (nm>0 ? 6mm : 622) : -6m2) ;
                                    else h['6']>0 ? 6 : -6
elif h['4']>0 or h['-4']>0 → tet   : order 16 → 4/mmm ; order 8 → inv ? 4/m :
                                    (h['4']>0 ? (nm>0 ? 4mm : 422) : -42m) ;
                                    else h['4']>0 ? 4 : -4
elif h['3']>0 or h['-3']>0 → trig  : order 12 → -3m ; order 6 → inv ? -3 :
                                    (nm>0 ? 3m : 32) ; else 3
elif h['2']>0 or nm>0      → o/m   : order 8 → mmm ; order 4 → inv ? 2/m :
                                    (h['2']>=3 ? 222 : mm2) ; else h['2']>0 ? 2 : m
else                       → inv ? -1 : 1
```

The crystal **class is derived from the rotation content itself**, not from the lattice metric —
deliberate, so that a structure whose symmetry is a proper subgroup of its lattice's holohedry (the
generic case partway up the tolerance ladder) is classified correctly.

**(f) Space-group symbol** (`spaceGroupHM()` → `spaceGroupSymbol.js`). Triclinic groups are `P1`
(No. 1) or `P-1` (No. 2) whatever cell describes them (unless item 4 marks them as a lower bound).
Every other group is named in a **standard setting**, which the RMC cell need not be:

1. **Symmetry elements** (`classifyElement()`): each operation's characteristic direction (rotation
   axis, or mirror normal), order, and — from the **intrinsic translation**
   $\mathbf t_\text{int} = \tfrac1n\sum_{k<n}R^k\mathbf t$ (independent of the origin) — whether a
   rotation is a screw $n_m$ ($m = \mathrm{round}(n\,\mathbf t_\text{int}\!\cdot\mathbf d/\mathbf d\!\cdot\!\mathbf d) \bmod n$
   along the axis $\mathbf d$, from the unreduced $\mathbf t_\text{int}$: reducing each component
   mod 1 is not a lattice translation along an axis with components of both signs) and which
   glide ($a,b,c,n,d$, or $e$) a reflection is ($\mathbf t_\text{int}$ reduced mod 1 and snapped to
   quarters). An operation stands for its coset $\{R\,|\,\mathbf t+\boldsymbol\ell\}$, and where
   the lattice projects onto the axis or plane in a fraction of its period ($\tfrac1n\sum_k R^k\boldsymbol\ell$
   not a lattice vector) the representatives differ in kind: the $[100]$ 2-fold of a hexagonal cell
   and the $2_1$ half a cell away, a cubic $\langle111\rangle$ 3-fold and its $3_1$, $3_2$, a mirror on
   a tetragonal or cubic diagonal and an $n$-glide. Which representative the finder holds is an
   accident of wrapping (noise turns an exact 0 into 0.9995), so the symbol reads **every** element
   of each coset (`cosetElements()`, $\boldsymbol\ell\in\{-1,0,1\}^3$) and the short-symbol rules
   pick the highest-priority one per position (rotation before screw, $m$ before glides). Before 1.0
   the held representative decided: a noisy $\{E, 2_{[100]}\}$ subgroup of $P6_3/mmc$ read `P2_1`
   No. 4 (it is `C2` No. 5 — item 2: its conventional cell is the orthohexagonal C cell).
2. **Candidate cells** (`hmSymbolInStandardSetting()`). The group is re-expressed (`applySetting()`:
   $R' = Q^{-1}RQ$, $\mathbf t' = Q^{-1}\mathbf t$, and the new cell's translations $Q^{-1}(\mathbb Z^3 + T)$
   mod 1) in a sequence of cells $Q$ (columns = new basis vectors in the old fractional basis), and
   named in the first that is conventional:
   - the given cell and its five other axis orders (`PERMUTATIONS`, all right-handed so a screw's
     handedness survives), for every crystal system;
   - cells built from the symmetry elements (`derivedBases()`), each basis vector the **shortest
     lattice vector** of the full translation lattice $\mathbb Z^3 + T$ along its direction (the
     candidates, `shortLatticeVectors()`, are $\mathbf n + \boldsymbol\tau$ with $\mathbf n$ over a
     $\pm3$ box of a **reduced** basis of $\mathbb Z^3$ — `reduceBasis()`, Step 7 — so the axes of a
     crystal on a strongly oblique cell are in reach; $\pm3$ covers the $c$ axis of an R lattice on
     a reduced rhombohedral basis):
     cubic — $a,b,c$ on the three 4-fold ($\bar4$ for $\bar43m$, 2-fold for $23$, $m\bar3$) axes;
     tetragonal — $c$ on the 4 / $\bar4$ axis, $a$ the shortest lattice vector $\perp c$ (and its
     45° diagonal), $b = 4\cdot a$; trigonal/hexagonal — $c$ on the 3-fold, $a = \pm$ the shortest
     lattice vector $\perp c$, $b = 3\cdot a$; orthorhombic — the three 2-fold axes / mirror normals
     in all six orders; monoclinic — $b$ on the unique axis and $(a, c)$ every unimodular pair from
     the two shortest lattice vectors $\perp b$ (covers the cell choices, e.g. I2/a → C2/c, and
     P2₁/n → P2₁/c), tried shortest $|a|^2+|c|^2$ first and then non-acute β — the reduced cell
     choice, so where several cell choices spell the symbol the letters are those of the reduced
     cell (the cell choices $a' = a + kc$ of P2₁/c swap which inversion centres are 2b and 2d).
     "Perpendicular" and "along" are decided by the rotations themselves
     ($R\mathbf v = \pm\mathbf v$, $\sum_k R^k\mathbf v = 0$), exactly, not by the metric.

   The new cell's translation set is $\mathbb Z^3 + T$ modulo $Q\mathbb Z^3$, generated mod 1 from
   $Q^{-1}\mathbf e_i$ and $Q^{-1}\boldsymbol\tau$ by closure (no search box), and must have
   $|T|\det Q$ members. A cell is kept only if its basis vectors are lattice vectors, it is
   right-handed, every rotation is an integer matrix in it, its translation set is exactly a
   Bravais centering (Step 10c) that a
   standard setting of the system uses (monoclinic P, C; orthorhombic P, A, C, I, F; tetragonal P,
   I; trigonal P, R obverse; hexagonal P; cubic P, I, F), every rotation has the **block form** of a
   conventional cell (`elementsFitSetting()` → `conventionalForm()`: the monoclinic $b$, or the
   tetragonal / trigonal / hexagonal $c$, is perpendicular to the other two basis vectors, so $R$
   maps it onto $\pm$ itself and the other two into their own plane — $R_{01}=R_{10}=R_{12}=R_{21}=0$,
   resp. $R_{02}=R_{12}=R_{20}=R_{21}=0$; orthorhombic rotations diagonal, cubic ones signed
   permutations), and every element lies along a direction family that may carry its type
   (`elementsFitSetting()`: a tetragonal 4-fold only on [001], cubic 3-folds only on
   $\langle111\rangle$, …). The element directions alone do not show that a cell is conventional:
   on the primitive cell $((\mathbf a+\mathbf b)/2, \mathbf b, \mathbf c)$ of a C-monoclinic lattice
   the 2-fold still runs along the second basis vector, and on
   $((\mathbf a+\mathbf b+\mathbf c)/2, (-\mathbf a+\mathbf b+\mathbf c)/2, \mathbf c)$ of an I-tetragonal
   one the 4-fold along the third, but the other vectors lean on the axis and the cell's pure
   translations are primitive where the conventional cell's are centred. Before the block-form
   check (and on 0.x) C2, Cm, Cc, C2/m, C2/c read `P2` No. 3 … `P2/c` No. 13, I4, I-4, I4/m
   `P4` No. 75, `P-4` No. 81, `P4/m` No. 83, and R3, R-3 `P3` No. 143, `P-3` No. 147 on such cells,
   with the P group's Wyckoff letters; and every monoclinic rung of a hexagonal or trigonal ladder
   was a P group — a 2-fold along $a_\text{hex}$ (or a mirror normal to it) of a P-hexagonal lattice
   has the orthohexagonal C lattice, so noisy wurtzite's `Pm` is `Cm` No. 8 and bismuth's
   `P2`, `P2/m` are `C2` No. 5, `C2/m` No. 12. Because derived cells use the full translation lattice,
   a supercell of the true cell is named in the true cell (CsCl in a $2\times2\times2$ cell: `Pm-3m`), and
   a subgroup that keeps its parent's centering is named in its own conventional cell (Ga shifted
   along [111] in the F-cubic lacunar spinel: `R3m` on hexagonal axes; tetragonally strained
   rocksalt in its F cell: `I4/mmm`).
3. **Symbol** (`hmSymbolCandidates()`): the system's ordered symmetry directions
   (`SYSTEM_DIRECTIONS`) are the positions of the Hermann–Mauguin symbol; each is filled with the
   elements along it by the short-symbol rules, and all defensible spellings are listed best-first.
   A candidate is **accepted** only if it is a tabulated standard symbol (230 groups plus the
   pre-2002 `e`-glide spellings, `spaceGroupTable.js`) **of the detected crystal class**
   (`pointGroupOfSymbol()`) and starts with the cell's centering letter. The number is looked up
   only for an accepted symbol; `Cmca`-style spellings are shown in the current form.
4. **Lower bound** (`allLatticeOpsTested()`, `latticeFullyTested()`). The finder only tries
   rotations that map the **given** cell's lattice $\mathbb Z^3$ onto itself (Step 7). When the
   given cell is a supercell of the crystal's own translation lattice $\mathbb Z^3+T$, rotations of
   that lattice which do not preserve $\mathbb Z^3$ were never tested — a perovskite in a
   $\sqrt2\times\sqrt2\times2$ or $2\times2\times1$ cell, or rocksalt in a $1\times1\times2$ cell, cannot test the
   cubic 3-folds. The check takes a primitive basis $P$ of $\mathbb Z^3+T$ (`primitiveBasis()`: its
   successive minima, which in three dimensions always form a basis), lists that lattice's
   rotations with the Step 7 scan at a strain of the group's worst residual (at least
   $10^{-3}$ Å), and requires every one, carried to the given cell ($PRP^{-1}$), to be an integer
   matrix. If one is not, the result is only a lower bound and is shown as `≥ P4/mmm`
   (`lowerBoundLabel()`), with no number and no Wyckoff letters — whatever the branch: a named
   symbol, `P1`/`P-1`, or a crystal class (`≥ 4/mmm class`). Before 1.0 only the naming cell's own
   rotations were checked, so rocksalt in a $1\times1\times2$ cell read `I4/mmm` No. 139.
5. **Crystal class.** If no cell gives an accepted symbol, the card shows `"<point group> class"`
   (e.g. `4/mmm class`, `classLabel()`) with no number and no letters. Positional assembly on a cell
   that is not conventional can spell another group's symbol (rocksalt on its primitive cell gives
   `Pmmm`); before 1.0 such spellings, and the symmorphic `centering + point group` fallback, were
   shown with that group's ITA number.

Checked on all 230 groups (the test fixtures) in all six axis orders, and on four oblique cells
(a 60° cell, a sheared cell and two strongly oblique unimodular cells: 839 group–cell pairs, exact
and with 0.005 Å noise): every group is named correctly (the two location-degenerate pairs by
Step 10g). Before 1.0, 753 of those 839 oblique descriptions got another group's number (rutile on
a 60° cell: `Cmmm` No. 65), because the Step 7 scan missed the rotations that need entries of
magnitude 2 there. The centred groups above on their primitive axis-keeping cells get the
conventional description's number and Wyckoff labels, and no P-named rung of the noisy
trigonal/hexagonal fixture ladders hides a centring (an independent check: along the principal
axis, the axial part $\tfrac1n\sum_k R^k\mathbf g$ of every lattice vector $\mathbf g$ is itself a
lattice vector; `symmetryHiddenCentring.test.js`).

**Code**: `symmetry.js` → `classifyOperations()`, `classifyRotation()`, `pointGroupOf()`,
`spaceGroupHM()`, `classLabel()`, `lowerBoundLabel()`, `POINT_GROUP_ORDER`; `spaceGroupSymbol.js` →
`intrinsicTranslation()`, `classifyElement()`, `cosetElements()`, `centeringOfOps()`, `bravaisCentering()`,
`applySetting()`, `derivedBases()`, `elementsFitSetting()`, `hmSymbolCandidates()`,
`hmSymbolInStandardSetting()`, `reduceBasis()`, `shortLatticeVectors()`, `primitiveBasis()`,
`allLatticeOpsTested()`, `latticeFullyTested()`, `pureTranslations()`, `twoFoldsMeet()`, `transformOps()`;
`spaceGroupTable.js` →
`SPACE_GROUPS`, `spaceGroupNumber()`, `pointGroupOfSymbol()`, `canonicalSymbol()`.

**(g) Location-degenerate pairs** (`twoFoldsMeet()`, `resolveLocationPair()`). I222/I2₁2₁2₁
(Nos. 23/24) and I23/I2₁3 (Nos. 197/199) contain the same element types along the same directions
(the I centering turns every 2-fold into a 2₁ half a cell away and vice versa), so the symbol
candidates cannot separate them; they differ in where the axes sit. In the standard setting a pure
2-fold along axis $i$, $\{R|\mathbf t\}$ with $t_i\equiv0$, fixes the line $2p_j \equiv t_j$ ($j\ne i$); the
2-folds along $a$, $b$, $c$ have a common point iff, for some choice of one pure 2-fold per axis,
their translations agree mod 1 on the shared components. They meet in I222 and I23 and never in
I2₁2₁2₁ and I2₁3, whatever the origin. Before 1.0 both members of each pair were reported as the
symmorphic one (with its number).

#### Step 11. Group closure: the largest closed group at each threshold

With noise, every operation of the true group has its own residual, so the operations with
$\varrho \le r$ are an arbitrary **subset** of the group — and a subset of the right *size* is not a
group. (Before 1.0 only the size was checked, and the ladder of the bundled demo showed sets such as
`{E, 2, m, m'}` with a missing product as "P2 (No. 3)".) Nothing is classified unless it is
**closed under composition** modulo lattice translations.

**Products.** $\{R_a|\mathbf t_a\}\{R_b|\mathbf t_b\} = \{R_aR_b\,|\,R_a\mathbf t_b+\mathbf t_a\}$ (`composeOps()`),
matched to the detected operation $k$ with the same rotation and the nearest translation, accepted
when the Cartesian translation mismatch is at most $\varrho_a+\varrho_b+\varrho_k+10^{-5}$ Å
(`productTable()`, evaluated lazily and memoised). The bound follows from the residuals: each
operation maps every site to within its residual of a same-element site, so the product maps every
site to within $\varrho_a+\varrho_b$, and two same-rotation operations sending a site to the same
partner differ by at most the sum of all three. Translations of one rotation are at least $\tau$
apart (Step 8 dedup), so the nearest one is the only candidate.

**The walk** (`groupsByThreshold()`). The distinct residuals are visited tight → loose; at each
threshold $r$ the newly admitted operations are offered, best residual first, to the current group
$H$ (which starts as the identity):

1. **Growth** (`growingGroup().extend()`): the group generated by $H$ and the candidate $x$ is built
   from generators — add $H\cdot x$, then right-multiply every new element by every generator (a
   finite set holding the identity and closed under right multiplication by the generators *is* the
   group they generate). It replaces $H$ only if every element of it has $\varrho \le r$. A candidate
   that fails waits for the operation that blocked it to be admitted; a product that does not exist
   among the detected operations blocks it for good.
2. **Elimination** (`closeUnder()`, only when there are at most 256 candidate operations): starting
   from every operation with $\varrho\le r$, repeatedly drop the one with the **worst** residual among
   those taking part in a missing product, until the set is closed. If that group is larger than
   $H$, it replaces $H$ (growth from a smaller group cannot reach a group that does not contain it).

A group closed at a tighter threshold is still a closed subset at a looser one, so the reported
group never shrinks as the tolerance loosens. Both procedures are greedy: the result is always
closed, and maximal in practice, but maximality is not proven.

**Cardinality pre-filter.** `isValidGroup(cls)` additionally requires

$$n_\mathrm{space} \;=\; |P|\times \max(n_\mathrm{trans},1),$$

which any closed set satisfies ($|P|$ = tabulated order of the point group, `POINT_GROUP_ORDER`;
$n_\mathrm{trans}$ = distinct pure translations, keyed on a $10^{-3}$ grid **after folding mod 1**, so
0.9997 and 0 count once). `classifyOperations(ops, tolFrac, { closed })` names a set only if it
passes this and is closed — by construction when the walk built it (`closed: true`), otherwise by an
all-pairs check (`isClosedSet()`, products matched within $3\cdot$`tolFrac` per fractional
component). A set that fails is labelled `not a group` with no number. The walk's products are
matched within the residuals, so this is not excluded by construction at loose thresholds. Before
1.0 it was reached through the Step 8 near-duplicates (a noisy $P4_322$ ladder ended in `not a
group`). It can still be reached: a 4-site noisy Pm structure read `not a group` (6–7 operations
with two pure translations) at 0.90–0.98 Å in a pass at that τ. In the ladder's 1 Å pass —
the one the card uses for its headline too (Step 12) — that structure reads `Pm` and `≥ Amm2`. The
ladder itself can still end in one: in two sweeps of noisy fixture groups at physical density in
random cells (2065 rungs), one rung was `not a group` (a two-orbit Ccc2 from 0.98 Å to the 1 Å
end); the earlier sweeps of 1380 noisy fixture ladders to 1 Å in their conventional cells and 790
well-known structures in eight cells found none. The card shows it as such, with no number.

**Code**: `symmetry.js` → `composeOps()`, `productTable()`, `growingGroup()`, `closeUnder()`,
`groupsByThreshold()`, `isClosedSet()`, `isValidGroup()`.

#### Step 12. The reported space group at the selected tolerance

`spaceGroupAtTolerance(A, basis, τ, passTol)` produces the headline of the Detected SG card:

1. Run one detection pass (Steps 7–9) at $\tau' = \max(\text{passTol}, \tau, 10^{-3})$ Å — its
   lattice strain bound, pairing radius $2\tau'$ and near-duplicate radius $\tau'$.
   `describeSymmetry` passes `passTol = max(τ, 1 Å)`, the ladder's own pass (Step 13); the default
   `passTol = τ` is a cheaper pass at τ alone.
2. Walk its residual thresholds $\le\tau$ (Step 11) and take the group holding at the last one.
3. Classify it (Step 10) with `tolFrac = max(τ, 10⁻⁶)/meanEdge(A)`.

`maxResidual` is the **worst residual of the operations returned** (0 only when that group is the
identity alone, whose residual is exactly 0) — shown in the card's tooltip as `fits to <maxResidual>
Å`. `nSpace` is the number of operations of that group ("Operations"). On an empty basis, or when no
operation survives — not even the identity: a lattice with a non-finite entry, or a singular one,
rejects every candidate in Step 7 — the result is `undetermined` (`UNDETERMINED`: no number, point
group `—`, 0 operations, `maxResidual` `NaN`), `describeSymmetry` returns no orbits, and the ladder
is empty. Before 1.0 it was `P1` / No. 1 — a space-group number for a structure never analysed.

With the ladder's pass the headline is **exactly** the group the ladder shows at $\tau$: the same
operations with the same residuals (the refinement of a seed does not depend on the acceptance
threshold, and the near-duplicate merge keeps the best fit first, so the operations kept below τ
are those the ladder keeps below τ), and the same walk up to τ over the same product table. A pass
at τ alone need not agree: its smaller pairing radius drops seeds the 1 Å pass refines, and its
smaller merge radius keeps near-duplicates the 1 Å pass merges. Before 1.0 the card's headline ran
such a pass, and in sweeps of noisy fixture groups (physical density, random cells) about 1 brick
midpoint in 400 disagreed with its brick: a noisy P-6m2 whose ladder read `P3m1` at 0.15 Å had the
headline `Cm` there, and a 4-site Pm read `not a group` (Step 11) inside a `Pm` brick
(`symmetryHeadlineLadder.test.js`). The price is that every headline runs the ladder's 1 Å pass —
tens of milliseconds for the 52-site GaNb₄Se₈ and GaTa₄Se₈ bases.

**Code**: `symmetry.js` → `spaceGroupAtTolerance()`, `UNDETERMINED`; glue in `symmetryModel.js` → `describeSymmetry()`.

#### Step 13. The tolerance ladder

`symmetryLadder(A, basis, tolMax, latticeTol)` builds the coloured brick strip from **one** detection
pass at the loosest tolerance (`tolMax`; the app passes **1.0 Å**):

1. Walk every distinct residual threshold $r_i$ (Step 11). The group at $r_i$ holds over
   $[r_i, r_{i+1})$; the last one extends to `tolMax`.
2. Classify each distinct group once (`tolFrac = tolMax/meanEdge`).
3. Consecutive thresholds with the same symbol are **merged**; the merged brick keeps its first
   `from` and takes the latest `to` and `nSpace`.

The first threshold is always 0 — the identity maps every site onto itself with residual exactly 0
— so the first brick starts at 0 without being forced there. Operation counts are non-decreasing
left → right (Step 11). Each brick carries `{ from, to, spaceGroup, spaceGroupNumber, pointGroup, nSpace }`.

**One-pass caveat.** The ladder's translation dedup radius (Step 8) is `tolMax` for every rung
(the card's headline shares it, Step 12): distinct translations of one rotation closer than 1 Å are
one operation at every tolerance.

**Code**: `symmetry.js` → `symmetryLadder()`; `symmetryModel.js` → `toleranceLadder()`.

#### Step 14. Site orbits (Wyckoff structure) and site symmetry

`siteOrbits(A, basis, ops, τ)` partitions the basis into symmetry orbits with a union–find
(disjoint-set with path halving):

- For every operation of the reported group and every site $i$: map $\mathbf x_i \to \mathbf y$, find the
  **nearest same-element** site $j$ with $d(\mathbf y,\mathbf x_j) < \tau$, and `union(i, j)`.
- Group members by root. The orbit's `size` is its multiplicity **in the given cell**; the
  representative is its **first member in basis order** (i.e. lowest reference number).
- The **site symmetry** is the stabiliser of the representative: the distinct rotation parts $R$ of
  operations with $d(R\mathbf x_\mathrm{rep}+\mathbf t,\ \mathbf x_\mathrm{rep}) < \tau$, fed through the same
  `pointGroupOf()`. It is derived from the operations, not looked up.
- Orbits are returned largest-first.

`describeSymmetry` calls `siteOrbits` with the user tolerance $\tau$, which is at least the worst
residual of the reported group (`maxResidual`), so the orbit and stabiliser tests are at least as
generous as the group itself.

**Wyckoff letters** (`symmetryModel.js` → `lettersInSetting()`, `wyckoff.js` →
`assignWyckoffLetters()`). The table (`wyckoffTable.js`) lists the Wyckoff positions of all 230
groups in their ITA standard setting (all 1731 positions; for the groups with two origin choices, the
one the table was built from — origin choice 2 for the centrosymmetric ones, e.g. Fd-3m). Each row
was checked against the group's own operations (its coordinate expanded to exactly its
multiplicity), and its site symmetry is the stabiliser computed the same way as above, so the two
are comparable.

1. Letters are read in the **standard cell the group was named in** (Step 10f). Each orbit's members
   are carried there, $\mathbf x' = Q^{-1}\mathbf x$, together with every translation of that cell, and
   its multiplicity is scaled by the cell-volume ratio (the four Ga of a lacunar spinel in its
   F-cubic cell are one 3a orbit of R3m on hexagonal axes). With no such cell — a crystal class,
   a `≥` lower bound, or P1/P-1 in a non-primitive cell — every letter is withheld.
2. The candidates are the rows with the orbit's multiplicity **and** site symmetry. One candidate
   → that letter. Several (Pm-3m 3c and 3d are both 3 × 4/mmm) → the tie is broken by the
   **coordinate form**: a letter is kept if some member of the orbit (any lattice-equivalent
   representative) fits its form, e.g. $(x,x,0)$, within `tolFrac` $= \tau/\overline{a'}$ per
   component ($\overline{a'}$ = mean edge of the naming cell; `fitsForm()` solves for the free
   parameters by ridge-regularised least squares and checks the residual). Exactly one fit → that
   letter; otherwise none.
3. A letter is returned with the multiplicity **of the cell it is read in**,
   `wyckoffMultiplicity`, beside the given-cell `size`: a Wyckoff label pairs the two, so the Ga
   of the lacunar spinel is `3a` (never `4a`, which R3m does not have) and rocksalt on its
   primitive cell is `4a`/`4b` (never `1a`/`1b`). `orbitLabel()` prints
   `"<wyckoffMultiplicity><letter>"`.
4. No letter is ever guessed: with none, `wyckoffMultiplicity` is `null` and `orbitLabel()` shows
   `"<given-cell multiplicity> (<site symmetry>)"`.

**Origin.** The standard cell is found by a change of basis only (Step 10f), so letters that need
the coordinate form assume the structure's origin is the table's origin or an equivalent one.
Letters that multiplicity and site symmetry fix alone do not depend on the origin; a tie at a
shifted origin usually fits no form and gets no letter (diamond described with its atoms at 0 and
¼ — origin choice 1 — gets no letter for its 8-fold site, because the table uses origin choice 2).

Orbits are not rendered on the Run Dashboard card itself; they flow into the AI-assistant context
(`runContext.js` → `symmetryContext()` → `symmetry.sites`, ranked by mean displacement, capped at 12
sites), where `multiplicity` is the orbit's size in the given cell and `wyckoff` should be the
naming-cell pair `wyckoffMultiplicity` + letter. **Open:** `symmetryContext()` still builds the
label as `` `${orbit.size}${orbit.wyckoff}` ``, so wherever the naming cell is not the given cell
the assistant is handed a label the named group does not have: `16i` for the Nb of the 5 K
GaNb₄Se₈ run at τ = 0.02 Å (P-4n2 is named in a cell half the F-cubic one; the card's pair is
`8i`), `4a`/`12b` for the R3m lacunar spinel in its F cell (`3a`/`9b`), `1a`/`1b` for rocksalt on
its primitive cell (`4a`/`4b`). The fix is one line in `runContext.js`,
`` `${orbit.wyckoffMultiplicity ?? orbit.size}${orbit.wyckoff}` ``; it reads a field of the
`symmetry` prop, so the `llm/` import boundary holds (it needs no `orbitLabel()` import). It is
outside the symmetry finder's files and handed to the `llm/` module's owner, with a test.

**Code**: `symmetry.js` → `siteOrbits()`; `symmetryModel.js` → `describeSymmetry()`,
`lettersInSetting()`, `orbitLabel()`; `wyckoff.js` → `wyckoffPositions()`, `parseCoordinateForm()`,
`fitsForm()`, `assignWyckoffLetters()`; data in `wyckoffTable.js`.

#### Step 15. What the card renders

- **Space group**: `symmetry.spaceGroup`, with `No. <n> · <pointGroup>` beneath, and a tooltip
  `Point group <pg> · fits to <maxResidual.toFixed(3)> Å`. The `No.` part is omitted whenever
  `spaceGroupNumber` is `null` — a crystal class, a `≥` lower bound, `undetermined` or
  `not analysed`. For the last two `maxResidual` is `NaN` and the tooltip currently prints
  `fits to NaN Å`, and the `reason` of a skipped structure is not shown anywhere on the card;
  the card's info badge also still says the symbol is reported in the given cell "up to an axis
  permutation", which Step 10f superseded. Those are `ModelSummary.jsx` wording issues, not
  finder results.
- **Operations**: `symmetry.nSpace`.
- **Space group vs. tolerance**: the ladder bricks. Brick **fill** is
  `color-mix(in srgb, var(--accent) P%, var(--panel-raised))` with

  $$\lambda = \begin{cases} \dfrac{\ln n_\mathrm{space}}{\ln n_\mathrm{space}^{\max}} & n_\mathrm{space}^{\max} > 1\\[6pt] 0 & \mathrm{otherwise}\end{cases}
  \qquad P = \operatorname{round}\!\left(12 + 74\,\lambda\right)$$

  where $n_\mathrm{space}^{\max}$ is the largest operation count over all bricks. The guard matters: a
  degenerate ladder with $n_\mathrm{space}^{\max} \le 1$ short-circuits to $\lambda = 0$, hence
  $P = 12$, instead of the indeterminate $\ln 1/\ln 1$ the bare formula implies. (Label text turns
  white above $P=50$.) Brick **width is deliberately not proportional to the tolerance range**: with
  $\ell$ bricks, the widest rung is pinned at 34 % and the remaining $\ell-1$ rungs share 66 %
  equally, because the full-symmetry rung otherwise occupies almost the whole axis. A single rung
  takes 100 %.
- Each brick's tooltip reads
  `<sg> (No. <n>) · holds <from.toFixed(2)>–<to.toFixed(2)> Å · <nSpace> ops`. Note the **two
  tooltips use different precisions**: the ladder prints 2 decimals, the headline
  `maxResidual.toFixed(3)` prints 3 — a brick labelled "holds 0.00–0.05 Å" and a headline "fits to
  0.048 Å" can be the same number at different rounding.
- Clicking a brick sets $\tau \leftarrow (\texttt{from}+\texttt{to})/2$; the brick is marked active
  when $\texttt{from} \le \tau < \texttt{to}$ — **a half-open interval, closed at `from` and open at
  `to`**. In practice this covers
  every reachable $\tau$: `bricks[0].from` is 0 (the identity's residual) and the last brick's `to` is `tolMax`, the
  bricks are contiguous, and the only setter is the midpoint click (which always yields
  $\tau < \texttt{tolMax}$), so exactly one brick highlights for any $\tau\in[0,1)$. A $\tau$ of
  exactly 1.0 Å would highlight nothing, but no code path produces it.

The tolerance itself is held in `SymTolContext`
([`symTolContext.js`](../../web_app/frontend/src/symTolContext.js)), a `[value, setValue]` pair
provided by [`App.jsx`](../../web_app/frontend/src/App.jsx) (`useState(0.2)`), so the selection persists
when switching between the Dashboard and the KDE/3D page. `ModelSummary.jsx` falls back to its own
local `useState(0.2)` if no provider is present.

**Code**: `ModelSummary.jsx` → `brickStyle()`, `brickWidth()`, the JSX for `.sym-ladder`.

#### Where and how often this runs

`ModelSummary.jsx` computes `describeSymmetry(structure, symTol)` and `toleranceLadder(structure, 1.0)`
in two `useMemo`s — that is **two full detection passes** per structure (the ladder one re-runs only
when the structure changes; the headline one re-runs on every tolerance change). A **third** pass
plus ladder runs in [`Dashboard.jsx`](../../web_app/frontend/src/components/Dashboard.jsx) to feed the
AI assistant, gated on `wantAssistantData` so it stays idle until the assistant page is opened.

The two parts run in **different places**. Parsing and the average-structure basis (Steps 1–5,
`structureFromRmc6f`) run inside [`workers/localStructureWorker.js`](../../web_app/frontend/src/workers/localStructureWorker.js),
a Web Worker, in browser-local mode — so the $O(N_\mathrm{atoms})$ pass is off the main thread (in
Flask mode the equivalent work happens server-side in Python). The **symmetry finder** (Steps 6–14,
`describeSymmetry` / `toleranceLadder`) is the part that runs **synchronously on the main thread**
inside `useMemo`, unlike the KDE and PCA-KDE paths, which use Web Workers.

**Cost and the basis cap.** The $3^9$ lattice scan (Step 7) is repeated on every pass. The
partner search of Steps 9 and 14 uses a **cell list** (`partnerIndex()`): same-element sites are
binned on a fractional grid whose bins are at least $r\,|\mathbf b_i|$ wide ($r$ = the matching
radius, $\mathbf b_i$ = reciprocal vectors), so an image only visits its own and the neighbouring
bins; each seed costs about $N$ times a few neighbours instead of $N\cdot\bar N_e$. Measured on a
random two-species basis (Node, one pass each): 2000 sites — `describeSymmetry` ≈ 0.1 s and the
ladder ≈ 0.1 s; 4000 sites ≈ 0.4 s + 0.2 s (before 1.0: 1.8 s + 5.3 s at 2000 sites, quadratic).
The closure walk (Step 11) adds $O(n\cdot\text{generators})$ product look-ups per group growth for
$n$ operations, plus the quadratic elimination only up to 256 operations; a $2\times2\times2$ supercell
of rocksalt (1536 operations) with noise takes ≈ 0.8 s per pass.

Because all of it is synchronous on the main thread, `symmetryModel.js` refuses a basis of more than
**`MAX_SYMMETRY_SITES` = 2000** sites: `describeSymmetry` returns
`{ skipped: true, spaceGroup: 'not analysed', pointGroup: '<N> sites > 2000 limit', nSpace: '—',
maxResidual: NaN, orbits: [], reason }` (the card shows the first two as headline and subtitle;
the residual is `NaN`, never `0`, because nothing was fitted and `0` would read as an exact fit — the
LLM context omits a non-finite residual), and `toleranceLadder` returns no bricks. A box whose `.rmc6f` declares a $1\times1\times1$ supercell with
one reference number per atom (a glass, an imported P1 configuration) reaches this; it is not a
unit-cell configuration, and a symmetry search on it would only ever return `P1`.

The site cap alone does not bound the time. If such a box is **crystalline**, every lattice rotation
holds with every pure translation of its repeat unit, so the candidate operations grow with the
square of the box: a $4\times4\times4$ rocksalt box (512 sites) has $48\times256 = 12\,288$, and the
ladder of a noisy $3\times3\times3$ one (216 sites, 5184) took over a minute (Node). Before either pass, `symmetryModel.js` therefore asks
`operationEstimate(A, basis, tol, limit)` (`symmetry.js`) for the lattice rotations admitted at the
tolerance (Step 7) times the pure translations $\{I\,|\,\mathbf t\}$ of the structure there
(Steps 8–9 for the identity alone, stopped as soon as the product passes `limit`). It is judged at
the ladder's loosest tolerance, $\max(\tau, 1.0)$ Å, for the headline as well, so the card either
analyses a structure at every tolerance or says why not. Above **`MAX_SYMMETRY_OPS` = 384** — twice
the 48 × 4 = 192 a correctly declared cell can reach (F centring), enough for a $2\times2\times2$
supercell of a primitive cubic cell — `describeSymmetry` returns the same `skipped` shape with
`pointGroup: '≥ <t> translations per cell'` and a reason naming the supercell, and
`toleranceLadder` returns no bricks.

### Parameters and defaults — model summary and symmetry

| Name | Where | Default | Units | Meaning |
| --- | --- | --- | --- | --- |
| `symTol` ($\tau$) | `App.jsx` `useState`, `ModelSummary.jsx` fallback, `Dashboard.jsx` fallback | `0.2` | Å (Cartesian) | atomic-position tolerance for accepting an operation |
| `tol` | `spaceGroupAtTolerance()` | `0.2` | Å | same, function-level default |
| `tol` | `findSpaceGroupOps()` | `0.1` | Å | never used from the app (callers always pass) |
| `tolMax` | `toleranceLadder()` (called with `1.0`) | `1.0` | Å | loosest tolerance the ladder explores |
| `tolMax` | `symmetryLadder()` signature | `1.5` | Å | never used (callers pass `1.0`) |
| `latticeTol` ($\tau_L$) | `findSpaceGroupOps`, `symmetryLadder`, `spaceGroupAtTolerance` | $\tau$ (the pass tolerance) | Å | largest Cartesian lattice strain $\varrho_L(R)$ of an admitted point operation (Step 7); the strain is also a floor on that operation's residual. Not forwarded by `symmetryModel.js`, so it follows $\tau$ |
| `tol` | `latticePointOps()` signature | `0.01` | Å | exported wrapper only; the finder calls `latticeCandidates()` with $\tau_L$ |
| `tol` | `siteOrbits()` signature | `0.1` | Å | never used (`describeSymmetry` always passes `symTol`) |
| `tolFrac` | `classifyOperations()` signature | `0.02` | cell fractions | never used (all callers pass `tol / meanEdge(A)`) |
| translation snap | `pureTranslations()` (`spaceGroupSymbol.js`) | grid $1/d$, $d$ = the smallest order dividing the group order $n$ with every component within `0.02`; $d\le25$ | cell fractions | pure translations are snapped as a group before the exact Bravais match (Step 10c) |
| lattice-completeness strain | `spaceGroupHM()` → `allLatticeOpsTested()` | the group's worst residual, at least `1e-3` | Å | strain up to which a rotation of the translation lattice counts in the lower-bound check (Step 10f) |
| translation-key granularity | `classifyOperations()` | `1e-3` | cell fractions | rounding used to count distinct pure translations (no fold of 1000 → 0) |
| Wyckoff coordinate tolerance | `assignWyckoffLetters()` (from `describeSymmetry`) | $\tau/\overline{a'}$ | cell fractions | per-component fit to a tabulated coordinate form, in the naming cell |
| threshold epsilon | `symmetryLadder`, `spaceGroupAtTolerance` | `1e-9` | Å | float-safety slack on `residual ≤ r` |
| $\tau$ floor (detection) | `spaceGroupAtTolerance()` | `1e-3` | Å | `Math.max(tol, 1e-3)` for the full detection pass |
| $\tau$ floor (`tolFrac`) | `spaceGroupAtTolerance()` | `1e-6` | Å | `Math.max(tol, 1e-6)` before dividing by `meanEdge(A)` for the classification tolerance — **distinct** from the `1e-3` detection floor |
| resultant floor | `structureFromRmc6f()` | `1e-6` | dimensionless | floor on $\bar R$ inside $\sqrt{-2\ln\bar R}$ |
| `maxPoints` | `Dashboard.jsx` (Flask query and `localStructureWorker` post) | `100` | atoms | 3D-preview subsample only; **not** used by the summary or symmetry |
| `maxPoints` = `STRUCTURE_MAX_POINTS` | `StructurePage.jsx` (Flask query and worker post) | `1000000` | atoms | KDE/3D page — effectively stride 1, no subsampling |
| `MAX_STRUCTURE_POINTS` | `web_app/backend/app.py` | `1_000_000` | atoms | server-side clamp: `max(100, min(maxPoints, 1e6))` |
| brick fill span | `ModelSummary.jsx` `brickStyle()` | `12 %`–`86 %` accent | — | $P = 12 + 74\lambda$ with $\lambda=\ln n/\ln n_{\max}$, and $\lambda = 0$ when $n_{\max}\le1$ |
| widest-brick width | `ModelSummary.jsx` `brickWidth()` | `34 %` | — | remaining `66 %` split evenly |
| search space | `latticePointOps()` | $3^9=19\,683$ | — | 6960 have $\lvert\det R\rvert=1$ |
| `MAX_SYMMETRY_SITES` | `symmetryModel.js` | `2000` | sites | larger bases are not analysed (`skipped: true`), and the ladder is empty |
| `MAX_SYMMETRY_OPS` | `symmetryModel.js` | `384` | operations | lattice rotations × pure translations (`operationEstimate`, at $\max(\tau,1.0)$ Å) above which a structure is not analysed (`skipped: true`) and the ladder is empty |
| `ELIMINATION_MAX_OPS` | `symmetry.js` | `256` | operations | the quadratic elimination of Step 11 runs only up to this many candidate operations |

### Caveats / what this is not

- **This is not spglib, and not FINDSYM.** It is a bounded approximation written for interactive use
  in a browser, with **no external space-group database** (no spglib tables, no WASM); it rests on
  the in-repo tables listed at the top of Part B (`SPACE_GROUPS`, `POINT_GROUP_SYSTEM`,
  `SYSTEM_DIRECTIONS`, `BRAVAIS`, `STANDARD_CENTERING`, `POINT_GROUP_ORDER`, the Wyckoff table). Use it
  to see *how* symmetry changes with tolerance on a disordered configuration, not to produce a
  published space-group assignment.
- **A group that cannot be named is shown as its crystal class.** The symbol is built positionally
  and accepted only when it is a tabulated symbol of the detected class and centering (Step 10f);
  otherwise the card shows e.g. `4/mmm class` with no number and no Wyckoff letters.
- **Supercells of the true cell give a lower bound.** Rotations are found on any basis of the given
  cell's lattice (Step 7, in a reduced basis), but only those that map that lattice onto itself. A
  symmetry of the crystal's own translation lattice that does not — the cubic 3-folds of a
  structure modelled in a $\sqrt2\times\sqrt2\times2$, $2\times2\times1$ or $1\times1\times2$ cell — is never
  tested, and the result is marked `≥ <symbol>` with no number (Step 10f). The cell handed to the
  finder is the `.rmc6f` box divided by its declared supercell, so this happens only when that
  declared cell is itself a supercell.
- **Lattice strain is measured on the atom scale.** A point operation is admitted when the Cartesian
  displacement it implies for the cell edges, $\varrho_L$ (Step 7), is within $\tau$, and that
  strain is a floor on the operation's residual. A strained cell therefore reads as the lower
  symmetry below its strain and as the higher one above it, and the ladder shows where. The strain
  is measured on the three cell edges; an atom far from the origin of a large cell can be displaced
  more by the same strain than an edge-length figure suggests. There is still no separate knob: the
  only user-adjustable parameter on this page is $\tau$, set by clicking a ladder brick.
- **The answer is tolerance-dependent by design.** An RMC configuration is disordered; the basis is
  a circular mean over supercell copies, and residual displacements of 0.05–0.5 Å are normal. There
  is no "correct" $\tau$ — the ladder exists precisely because the reported group is a function of
  $\tau$, and the default 0.2 Å is a UI convenience, not a physical constant.
- **The acceptance test is a covering test, not a bijection.** Distinct sites may map onto the same
  partner without being rejected. There is no check that an accepted operation permutes the basis.
- **Closure is enforced, maximality is heuristic.** Every reported operation set is closed under
  composition (Step 11), but the largest closed group at a threshold is found greedily (growth, plus
  elimination for up to 256 operations), not by enumerating subgroups.
- **Minimum image is component-wise.** `cartDist()` rounds each fractional component independently,
  which is exact for orthogonal cells and can over-estimate distances for strongly oblique ones.
- **Translation dedup radius equals the detection tolerance.** At the ladder's `tolMax = 1.0` Å,
  distinct translations less than 1 Å apart are merged for *every* rung. On a small cell this can
  suppress real centering vectors.
- **The headline pays for the ladder's pass.** The card's headline is walked from the ladder's
  own 1 Å detection pass (Step 12), so it is exactly the ladder's group at $\tau$, and every τ change
  reruns that pass. A direct `spaceGroupAtTolerance(A, basis, τ)` call without `passTol` runs a
  cheaper pass at $\tau$ and can differ from the ladder at a brick midpoint.
- **A merged brick shows its loosest rung's operation count.** Merging keeps `from` and takes the
  latest `to` and `nSpace` (Step 13).
- **Wyckoff letters assume the table's origin.** Letters are read in the standard cell the group is
  named in, for all 230 groups, but no origin shift is searched: a tie between positions of equal
  multiplicity and site symmetry is broken only if the structure's origin is the table's (or an
  equivalent one), and otherwise left without a letter (Step 14).
- **No origin shift.** The standard cell is found by a change of basis only; the origin stays where
  the `.rmc6f` puts it. Space-group names do not depend on the origin, Wyckoff letters do (Step 14).
- **Header input is unvalidated.** Neither the `Lattice` numbers nor the `Supercell` multiplicities
  are checked. A non-finite or singular lattice gives `NaN` strains, every candidate operation is
  rejected (Step 7), and the card shows `undetermined` with 0 operations, no ladder and `NaN` cell
  edges (Step 12). A zero or non-integer
  supercell entry is guarded in the *cell* division but not in the *fold into one cell*, giving a
  collapsed one-site basis and a spurious high-symmetry answer (Steps 1 and 5).
- **A reference site's species is whatever appeared first.** `acc.element` is set once and never
  re-checked, while `atomIndices` is keyed per (element, reference number). A chemically mixed or
  swapped site is presented to the finder as a single pure species, and the card's "sites"
  sub-labels can exceed the number of basis sites actually analysed (Step 5).
- **Orbits and site symmetries are evaluated at a looser radius than the group.** `siteOrbits` is
  called with the raw $\tau$, not with the `maxResidual` at which the group was accepted (Step 14).
- **No cell volume, no number density** is reported by this page (see the note after Step 5).
- **Flask/server-directory mode has no Detected SG card** — `/api/structure` returns no basis. It is
  also a different parser: the Python `iter_rmc6f_atoms()` skips every line with fewer than 9 fields,
  so coords-only `.rmc6f` files yield zero atoms there, and it capitalizes element tokens while
  `read_atom_indices()` (same response) does not (Step 2).
- **What the tests pin.** `symmetry.test.js` recovers all 230 fixture groups (built from ITA
  generators) from their atoms, names rocksalt, perovskite, diamond, hcp and an I4/mcm perovskite,
  and checks the point-group and space-group tables; `symmetrySettings.test.js` names the 230
  groups in all six axis orders and well-known structures in centred, primitive, rhombohedral and
  supercell cells, with the principal axis off $c$, and as subgroups of a centred parent;
  `symmetryObliqueCells.test.js` pins the lattice-rotation count on oblique bases, naming on oblique
  and strongly oblique cells, the lower bound and the translation snap; `symmetryFinder.test.js`
  covers basis-order independence, lattice strain, closure of every rung (the bundled GTS_250K demo,
  a noisy lacunar spinel) and class-consistent symbols; `symmetryRepresentatives.test.js` checks
  that a symbol does not depend on the lattice representative of a translation, and that every
  noisy trigonal, hexagonal and cubic subgroup is named after a group with the same element types
  per coset (rotation or screw, mirror or glide); `symmetryHiddenCentring.test.js` names the
  centred groups C2 … C2/c, I4 … I4/m, R3 and R-3 on primitive cells that keep the axis as a basis
  vector (with the conventional description's Wyckoff labels) and checks that no P-named rung of a
  noisy trigonal or hexagonal ladder hides a centring; `symmetryHeadlineLadder.test.js` that the
  card's headline is the ladder's group at every brick midpoint and that the bricks do not move
  with the order of the basis; `symmetryDuplicates.test.js` checks that of two
  near-duplicate operations the better-fitting one is kept, so the ladder reaches the group the
  headline finds whatever the site order;
  `symmetryWyckoff.test.js` and `wyckoff.test.js` check every Wyckoff row against its group's
  operations, the demo's letters and letters in permuted settings, and
  `symmetryWyckoffLabels.test.js` that a letter is paired with the multiplicity of the cell it is
  read in; `symmetryLimits.test.js`, `symmetryPerformance.test.js` and
  `symmetryUndetermined.test.js` cover the basis cap, the operation budget and an unanalysable
  lattice. The two GaNb₄Se₈ runs are gitignored, so they are not
  in the suite; every rung of their ladders was checked by hand to be closed and to have the
  element types of the group it is named after.
