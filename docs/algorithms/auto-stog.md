# Auto StoG — algorithm reference

> Part of the [Algorithms and Math Reference](../ALGORITHMS.md). Every step is anchored to the source; if this document and the code disagree, **the code wins**.

Pre-processing: putting a measured total-scattering $S(Q)$ on absolute scale and writing the classic stog / RMCProfile-ready file family. Four sections: the inputs and composition constants, the transform layer, the auto-scaling estimator, and the page workflow with its outputs.

## Contents

- [Auto StoG — Step 0: inputs, composition constants, and data conditioning](#auto-stog--step-0-inputs-composition-constants-and-data-conditioning)
  - [What this page shows](#what-this-page-shows)
  - [Step 1 — Accepting files (browser page only)](#step-1--accepting-files-browser-page-only)
  - [Step 2 — Reading the S(Q) columns](#step-2--reading-the-sq-columns)
  - [Step 3 — The `.dat` metadata header (`KEY :: value`)](#step-3--the-dat-metadata-header-key--value)
  - [Step 4 — The classic `stog.inp` (optional), and the divide-vs-multiply footgun](#step-4--the-classic-stoginp-optional-and-the-divide-vs-multiply-footgun)
  - [Step 5 — Chemical formula → atomic fractions](#step-5--chemical-formula--atomic-fractions)
  - [Step 6 — Sears/NIST scattering lengths → the Faber-Ziman coefficients](#step-6--searsnist-scattering-lengths--the-faber-ziman-coefficients)
  - [Step 7 — The composition-derived $S(0)$ target](#step-7--the-composition-derived-s0-target)
  - [Step 8 — Number density $\rho_0$](#step-8--number-density-rho_0)
  - [Step 9 — Assembling the configuration: what is required, what is derived](#step-9--assembling-the-configuration-what-is-required-what-is-derived)
  - [Step 10 — Cropping $Q$ to $(0,\; Q_\mathrm{min}\ldots Q_\mathrm{max}]$](#step-10--cropping-q-to-0-q_mathrmminldots-q_mathrmmax)
  - [Step 11 — Despiking (optional, OFF by default)](#step-11--despiking-optional-off-by-default)
  - [Step 12 — The σ column and the x-ray / non-neutron path](#step-12--the-σ-column-and-the-x-ray--non-neutron-path)
  - [Step 13 — What is handed to the fit](#step-13--what-is-handed-to-the-fit)
  - [Parameters and defaults](#parameters-and-defaults)
  - [Python ↔ JavaScript parity for Step 0](#python--javascript-parity-for-step-0)
  - [Caveats / what this is not](#caveats--what-this-is-not)
- [Auto StoG — the transform layer: Fourier pair, Lorch, filter, Keen conversions](#auto-stog--the-transform-layer-fourier-pair-lorch-filter-keen-conversions)
  - [What this page shows](#what-this-page-shows-1)
  - [Step 1 — Algebraic conversions (Keen 2001 definitions)](#step-1--algebraic-conversions-keen-2001-definitions)
  - [Step 2 — The grids](#step-2--the-grids)
  - [Step 3 — The discrete sine transform pair](#step-3--the-discrete-sine-transform-pair)
  - [Step 4 — The Lorch modification window](#step-4--the-lorch-modification-window)
  - [Step 5 — The analytic omitted-low-$Q$ correction](#step-5--the-analytic-omitted-low-q-correction)
  - [Step 6 — The assembled forward transform](#step-6--the-assembled-forward-transform)
  - [Step 7 — The Fourier filter](#step-7--the-fourier-filter)
  - [Step 8 — Conversion of the filtered result to the RMCProfile functions](#step-8--conversion-of-the-filtered-result-to-the-rmcprofile-functions)
  - [Step 9 — Low-$r$ enforcement (classic `stog` parity)](#step-9--low-r-enforcement-classic-stog-parity)
  - [Parameters and defaults](#parameters-and-defaults-1)
  - [Python ↔ JavaScript parity](#python--javascript-parity)
  - [What the tests assert against real data](#what-the-tests-assert-against-real-data)
  - [Caveats / what this is not](#caveats--what-this-is-not-1)
- [Auto StoG — the auto-scaling engine: level sweep, closed-form fit, self-consistency](#auto-stog--the-auto-scaling-engine-level-sweep-closed-form-fit-self-consistency)
  - [What this page shows](#what-this-page-shows-2)
  - [Step 1 — From files to arrays: parse, resolve, crop, (optionally) despike](#step-1--from-files-to-arrays-parse-resolve-crop-optionally-despike)
  - [Step 2 — The model and its physics targets](#step-2--the-model-and-its-physics-targets)
  - [Step 3 — The high-$Q$ level sweep](#step-3--the-high-q-level-sweep)
  - [Step 4 — Anchoring the offset: $b = 1 - aL$](#step-4--anchoring-the-offset-b--1---al)
  - [Step 5 — The closed-form affine solve](#step-5--the-closed-form-affine-solve)
  - [Step 6 — Huber IRLS robust re-weighting](#step-6--huber-irls-robust-re-weighting)
  - [Step 7 — The self-consistent loop with the Fourier filter](#step-7--the-self-consistent-loop-with-the-fourier-filter)
  - [Step 8 — $r_0$ detection from the first-shell $|g|$ flank, and the refinement pass](#step-8--r_0-detection-from-the-first-shell-g-flank-and-the-refinement-pass)
  - [Step 9 — The alternative amplitude: the $Q\to 0$ Faber-Ziman criterion](#step-9--the-alternative-amplitude-the-qto-0-faber-ziman-criterion)
  - [Step 10 — `estimate_rho0`: the density from criteria concordance](#step-10--estimate_rho0-the-density-from-criteria-concordance)
  - [Step 11 — Final pipeline and outputs](#step-11--final-pipeline-and-outputs)
  - [Step 12 — Every diagnostic returned, and how to read it](#step-12--every-diagnostic-returned-and-how-to-read-it)
  - [Python ↔ JavaScript parity](#python--javascript-parity-1)
  - [Parameters and defaults](#parameters-and-defaults-2)
  - [Caveats / what this is not](#caveats--what-this-is-not-2)
- [Auto StoG — page workflow, outputs, and the written file family](#auto-stog--page-workflow-outputs-and-the-written-file-family)
  - [What this page shows](#what-this-page-shows-3)
  - [Step 0 — Where each piece of the computation runs](#step-0--where-each-piece-of-the-computation-runs)
  - [Step 1 — Page-local upload and source classification](#step-1--page-local-upload-and-source-classification)
  - [Step 2 — Parsing the selected source and prefilling the form](#step-2--parsing-the-selected-source-and-prefilling-the-form)
  - [Step 3 — Grouped parameter fieldsets: what each knob changes in the engine](#step-3--grouped-parameter-fieldsets-what-each-knob-changes-in-the-engine)
  - [Step 4 — Config resolution and the $\rho_0$ self-consistency trigger](#step-4--config-resolution-and-the-rho_0-self-consistency-trigger)
  - [Step 5 — The worker call protocol](#step-5--the-worker-call-protocol)
  - [Step 6 — Post-engine work done inside the worker](#step-6--post-engine-work-done-inside-the-worker)
  - [Step 7 — The numeric readout cards](#step-7--the-numeric-readout-cards)
  - [Step 8 — The three plots](#step-8--the-three-plots)
  - [Step 9 — Export: the written file family](#step-9--export-the-written-file-family)
  - [Step 10 — The provenance JSON](#step-10--the-provenance-json)
  - [Parameters and defaults](#parameters-and-defaults-3)
  - [Cross-engine agreement (browser vs Python)](#cross-engine-agreement-browser-vs-python)
  - [Caveats / what this is not](#caveats--what-this-is-not-3)

---

## Auto StoG — Step 0: inputs, composition constants, and data conditioning

### What this page shows

**Auto StoG** is the app's *pre*-processing tab. It takes a measured total-scattering
structure factor $S(Q)$ — the rebinned output of a data-reduction pipeline — and produces the
classic stog / RMCProfile-ready file family ($S(Q)$, $g(r)-1$, filtered pair, $F_K(Q)$,
$G_K(r)$, $D(r)$, `ft.dat`) with the scale and offset determined from physics instead of the
classic stog interactive "try again" loop.

This section documents everything that happens **before any scaling math**: how files are read
and interpreted, how a chemical formula becomes scattering constants, how the number density is
resolved, and how the $Q$ array is cropped and (optionally) despiked. The fit itself
(level sweep, closed-form $(a,b)$ solve, Fourier-filter loop, $r_0$ detection, diagnostics) is
covered in the following sections.

Two facts frame everything below:

1. **The page is fully client-side in both runtimes.** `AutoStogPage.jsx` uploads files into
   page memory, parses them with the JS port, and runs the engine in a Web Worker. It does
   **not** call `/api/scaling/*` (those endpoints remain for CLI/API consumers) and it is
   independent of the run folder used by the other tabs. Files never leave the browser.
   Reference: [`web_app/frontend/src/components/AutoStogPage.jsx`](../../web_app/frontend/src/components/AutoStogPage.jsx),
   [`web_app/frontend/src/workers/autoScaleWorker.js`](../../web_app/frontend/src/workers/autoScaleWorker.js).
2. **Two engines, one math.** The reference implementation is Python
   ([`rmc_toolkits/scaling.py`](../../rmc_toolkits/scaling.py),
   [`rmc_toolkits/scattering.py`](../../rmc_toolkits/scattering.py),
   [`rmc_toolkits/parsers.py`](../../rmc_toolkits/parsers.py),
   [`rmc_toolkits/transforms.py`](../../rmc_toolkits/transforms.py)), driven by the `rmc-autoscale`
   CLI ([`rmc_toolkits/scaling_cli.py`](../../rmc_toolkits/scaling_cli.py)) and the Flask API.
   The browser runs a straight port,
   [`web_app/frontend/src/workers/autoScale.js`](../../web_app/frontend/src/workers/autoScale.js),
   parity-tested against Python-generated goldens
   ([`tests/generate_autoscale_fixture.py`](../../tests/generate_autoscale_fixture.py) →
   [`src/__tests__/autoScale.test.js`](../../web_app/frontend/src/__tests__/autoScale.test.js)).
   Per-step agreement is stated explicitly below; a summary table is at the end.

#### Notation and units

| Symbol | Meaning | Units |
| --- | --- | --- |
| $Q$ | momentum transfer | Å⁻¹ |
| $S(Q)$ | normalized total-scattering structure factor (Keen 2001 Eq. 19), $\to 1$ at high $Q$ | dimensionless |
| $\sigma(Q)$ | per-point uncertainty of $S(Q)$ (optional third data column) | same as $S$ |
| $r$ | real-space distance | Å |
| $r_0$ | closest interatomic approach (first-shell onset) | Å |
| $r_\mathrm{cut}$ | Fourier-filter cutoff radius | Å |
| $c_i$ | atomic fraction of element $i$, $\sum_i c_i = 1$ | — |
| $b_i$ | bound coherent neutron scattering length of element $i$ | fm |
| $\langle b\rangle^2 = (\sum_i c_i b_i)^2$ | classic stog "Faber-Ziman coefficient" | barn (also fm² internally) |
| $\langle b^2\rangle = \sum_i c_i b_i^2$ | Laue/self term, sets the $Q\to0$ limit | barn (also fm²) |
| $\rho_0$ | atomic number density | Å⁻³ |
| $\rho_m$ | mass density | g/cm³ |
| $M$ | molar mass of one formula unit | g/mol |
| $n$ | atoms per formula unit | — |
| $a, b$ | affine correction $S_\mathrm{corr}(Q) = a\,S_\mathrm{meas}(Q) + b$ | — |
| **C1** | the *high-$Q$ residual block* of the fit: $S_\mathrm{corr}\to1$ on $Q \ge Q_\mathrm{tail,min}$ | — |
| **C2** | the *low-$r$ residual block*: $g_\mathrm{corr}\to0$ on $[r_\mathrm{fit,lo}, r_\mathrm{fit,hi}]$ | — |
| $f_\mathrm{tail}$ | tail fraction of the $Q$ window that defines C1 (0.15) | — |
| $Q_\mathrm{tail,min}$ | start of the C1 window | Å⁻¹ |
| $r_\mathrm{fit,lo}$, $r_\mathrm{fit,hi}$ | bounds of the C2 window | Å |
| $a_\mathrm{density}$ | amplitude implied by the low-$r$ density limit (the default criterion) | — |
| $a_\mathrm{fz}$ | amplitude implied by the $Q\to0$ Faber-Ziman limit (independent criterion) | — |
| concordance | $a_\mathrm{fz}/a_\mathrm{density}$ — equals 1 when the two criteria agree | — |
| $\delta S(Q)$ | the Fourier-filter correction term (`delta_sq`), held fixed inside each fit iteration | — |

C1 and C2 are the two stacked least-squares blocks that the whole configuration exists to
define; the fit itself is documented in the following section — here they matter only because
every window, weight, and threshold below feeds one of them.

Unit convention throughout: $1\ \mathrm{barn} = 100\ \mathrm{fm}^2$. Scattering lengths are stored
in fm; every coefficient that reaches the engine is in **barns**.

---

### Step 1 — Accepting files (browser page only)

**Inputs:** files chosen through the file picker or dropped on the drop zone.
**Operation:** each `File` is read as text with `file.text()` and kept in React state as
`{name, text}`. Only two name patterns are accepted (`AutoStogPage.jsx`):

- **stog input candidate** — `isInpCandidate`: name ends in `.inp`, or is exactly
  `stog_input.dat`. **Both tests are case-sensitive** (`name.endsWith('.inp')`,
  `name === 'stog_input.dat'`).
- **data candidate** — `isDataCandidate`: name matches `/\.(sq|fq|dat)$/i` (case-**in**sensitive)
  **and** does not start with `scale.`/`scale_`/`ft.dat` **and** does not end in
  `_rmc.fq|gr|dr`. The exclusions exist so that dropping a whole finished stog run does not
  offer that run's *outputs* as new inputs.

Two consequences of the literal matching are worth knowing, because neither produces a
diagnostic: a file named `SAMPLE.INP` or `RUN.INP` matches **neither** predicate and is silently
discarded (the page then reports "No usable files" for a file it fully supports), and
`stog_input.dat` matches **both** — `selectSource` tests `isInpCandidate` first, so a plain data
file that happens to carry that name is always parsed as a classic stog input and fails with a
stog-input error rather than being read as data.

Everything else is silently ignored; if nothing is accepted the page reports either "drop
individual files, not a folder" (when a read failed, e.g. a directory entry) or a hint listing
the accepted formats.

**Merging and selection (`ingestFiles`).** Each drop/pick is *merged* into the existing upload
set, not appended: an incoming file replaces any already-held entry with the same name (so
re-dropping a corrected file silently swaps its contents), and the merged list is sorted with
`.inp` candidates first, then by `localeCompare`. Auto-selection ranges over the **newly
accepted** files only — `accepted.find(isInpCandidate) || accepted[0]`, i.e. drop order, not the
sorted dropdown order. Selecting a `.inp` requires its referenced data file to have been
uploaded too (`stog input references '<name>' — upload that file too`).

**Session persistence (browser only).** Every change to the parameter form and to the Advanced
panel's open state is written to `sessionStorage['autostog-session']` and restored on mount
(`STORAGE_KEY`, the two `useEffect`s near the top of the component). **Uploads are not
persisted; every physics parameter is** — composition, $\rho_0$, mass density, $\langle b\rangle^2$/$\langle b^2\rangle$ overrides,
the $Q$ window, $r_\mathrm{cut}$, $r_0$, the fit-window overrides, all toggles, and the manual
$(a,b)$. A reload therefore keeps the *previous* sample's physics and loses only its data file.
This is the mechanism behind the stale-override warning described in Step 12. **Reset params**
clears the form to `EMPTY_FORM` and then re-runs `selectSource`, which re-prefills from the
currently selected file (see Step 8).

**Outputs:** `dataRef.current = { q, sq, sigma, inp, header, name }` — the arrays and metadata
that the rest of the pipeline consumes.

**Code:** `AutoStogPage.jsx` → `isInpCandidate`, `isDataCandidate`, `ingestFiles`,
`selectSource`. In CLI mode the equivalent is `scaling_cli.main()`, which takes exactly one of
a positional `stog.inp` or `--data FILE`.

---

### Step 2 — Reading the S(Q) columns

**Inputs:** the raw text of the data file.
**Operation (both engines):** a tolerant whitespace-column reader.

1. Split each line on whitespace; skip lines with fewer than 2 tokens (this drops the classic
   **count header**, which is a single integer such as `        3401`, and any stray scalar
   line).
2. Attempt to convert *every* token on the line to a float. If any token fails, the whole line
   is skipped (this drops the title line and text banners). `NaN` tokens **succeed** — the
   NaN padding rows that rebinned files carry below the data's own $Q_\mathrm{min}$ are
   deliberately retained at this stage and masked later (Step 9).
3. Group the surviving rows by their token count; **keep the group with the most rows**. This
   is what prevents a stray two-number header line (e.g. `count qmin` printed on one line) from
   becoming the column template and discarding every real data row.
4. Return the kept rows transposed: column 0 is $Q$, column 1 is $S(Q)$, column 2 (if present)
   is $\sigma$.

**More than three columns.** The modal-group rule can select a 4-, 5- or 6-column group, and
every consumer then takes **index 2 as $\sigma$ regardless of what it actually holds** —
`AutoStogPage.jsx` (`columns.length >= 3 ? columns[2] : null`), `scaling_cli.py::_load_dataset`
and `app.py::_cached_scaling` (`columns.shape[0] >= 3`). Columns 3 and beyond are dropped
without a warning, and nothing checks that column 2 looks like an uncertainty. On a >3-column
file, turn the **σ column** toggle off unless you know the third column is the error bar.

**Outputs:** `q`, `sq`, and optionally `sigma`, all in file order.

**Code:** [`rmc_toolkits/parsers.py`](../../rmc_toolkits/parsers.py) → `read_stog_xy()`;
JS port `autoScale.js` → `readStogXy()`. Writer counterpart: `write_stog_xy()` /
`writeStogXy()` emit `count` line, title line, then `%.16E` columns.

**Python ↔ JS divergences (real, small):**

- Python uses `float(token)`; JS uses a regex gate
  `^[+-]?((\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?|nan|inf(inity)?)$` followed by a `D→e` rewrite.
  A **Fortran `D` exponent** (`1.234D+00`) therefore parses in the browser and is *skipped* by
  the Python reader (that row is dropped; if all rows use `D`, Python raises "does not contain
  STOG numeric rows").
- Tie-breaking when two column-count groups have the same number of rows follows insertion
  order in both engines, i.e. the group whose first row appeared earliest wins.
- Neither engine sorts, deduplicates, or checks that $Q$ is monotonically increasing. Ascending
  $Q$ is assumed by every downstream trapezoid transform.

**Test coverage:** `autoScale.test.js` → *"readStogXy keeps NaN padding and picks the modal
column count"* and the `writeStogXy` → `readStogXy` round-trip.

---

### Step 3 — The `.dat` metadata header (`KEY :: value`)

**Inputs:** the same data-file text.
**Operation:** every line containing `::` is split at the first occurrence; the key is
upper-cased and stripped, the value stripped. Three fields are then interpreted:

| Key | Parsed as | Used for |
| --- | --- | --- |
| `TITLE` | string | parsed, then surfaced **only** by the Flask API response (`app.py::_header_payload`); neither the browser page nor the CLI ever reads it |
| `NUMBER_DENSITY` | first token that parses as a float (e.g. `0.057329 Angstrom^(-3)` → `0.057329`) | $\rho_0$ in Å⁻³ (Step 8) |
| `MINIMUM_DISTANCES` | **minimum** over all float tokens (e.g. `2.4 2.2` → `2.2`) | $r_0$, the closest approach (Step 9) |

No unit conversion is attempted — the numeric token is taken at face value, so
`NUMBER_DENSITY` **must already be in atoms/Å³** and `MINIMUM_DISTANCES` in Å.

**Outputs:** `{raw, title?, number_density?, min_distance?}`.

**Code:** `parsers.py` → `read_dat_header()`; JS `autoScale.js` → `readDatHeader()`. The split
itself matches (the JS version slices at `indexOf('::')`, the Python version uses
`str.partition`), and the two agree on all realistic headers — but the numeric gate differs and
there are three edge divergences:

- **Non-finite tokens.** Python's test is `float(token)`, which **accepts** `nan`, `inf`,
  `-inf`. JS's is `Number.isFinite`, which skips such a token and tries the next one. So
  `NUMBER_DENSITY :: nan` yields `number_density = nan` in Python — later rejected by
  `ScalingConfig` with the unrelated-sounding *"rho0 must be finite and positive"* — while the
  browser silently ignores it and falls through to the next precedence level. A `nan` among
  `MINIMUM_DISTANCES` values likewise poisons Python's `min()` and is dropped by JS.
- **Literal forms.** Python's `float()` accepts underscore separators (`1_0` → 10.0) and rejects
  hex; JS's `Number()` accepts `0x10` (→ 16) and rejects `1_0`.
- **Empty values and encoding.** Python sets `title` whenever the key is present (`if "TITLE" in
  raw`), including for an empty value; JS uses a truthiness test (`if (raw.TITLE)`) and omits
  the key entirely for `TITLE ::` with nothing after it. Python opens the file with
  `errors="replace"`, the browser decodes with `File.text()` (UTF-8), so a latin-1 byte inside a
  key can produce a different key string in the two engines.

Test: `autoScale.test.js` → *"readDatHeader parses :: metadata"*;
`tests/test_scaling_cli.py::DataModeTests::test_header_metadata_and_formula` asserts the header
values reach the config (`rho0 = 0.05`, `r0 = 2.65`).

---

### Step 4 — The classic `stog.inp` (optional), and the divide-vs-multiply footgun

**Inputs:** the text of a classic single-dataset, filter-enabled stog input.
**Operation:** blank lines are stripped, then the remaining lines are read **positionally**.
The first column below is the 1-based line number of the surviving (non-blank) lines — the
convention used everywhere else in this document; the second is the 0-based `lines[i]` index the
code actually writes:

| Line (1-based) | `lines[i]` | Field | Notes |
| --- | --- | --- | --- |
| 1 | 0 | number of data files | must be `1`, else rejected |
| 2 | 1 | data file name | resolved relative to the `.inp` (CLI) or matched against the uploaded set (browser) |
| 3 | 2 | `qmin qmax` | Å⁻¹ |
| 4 | 3 | `yoffset yscale` | **Fortran divide convention** (see below); `yscale` must be finite and non-zero, **and `yoffset` must be finite** |
| 5 | 4 | Q offset | must be `0`, else rejected |
| 6, 7 | 5, 6 | output names: scaled $S(Q)$, unfiltered $g-1$ | |
| 8 | 7 | `rmax` | Å |
| 9 | 8 | `nr` | number of $r$ points |
| 10 | 9 | Lorch flag | `Y…`/anything-else |
| 11 | 10 | $\rho_0$ | Å⁻³ |
| 12 | 11 | second y offset | must be `0`, else `NotImplementedError`. What the Fortran does with a non-zero value is **not modelled anywhere in this codebase** — the reader only refuses it |
| 13 | 12 | "try again" flag | must be **not** `Y`, else rejected (that loop is what Auto StoG replaces) |
| 14 | 13 | Fourier-filter flag | must be `Y`, else rejected |
| 15 | 14 | filter $r$ cutoff | Å |
| 16, 17 | 15, 16 | filtered output names | |
| 18 | 17 | "Faber-Ziman coefficient" $\langle b\rangle^2$ | **barns** |
| 19–21 | 18–20 | RMC output names (`.fq`, `.gr`, `.dr`) | |
| 22 | 21 | `peak_cutoff peak_rmin peak_rmax` | first-peak ripple-cleanup window, Å |

At least 22 non-empty lines are required. Boolean fields are decided by
`token.strip().upper().startswith("Y")`.

**The scale convention.** Classic Fortran `stog_new3` *divides*:

$$S_\mathrm{scaled}(Q) = \frac{S_\mathrm{raw}(Q)}{\texttt{yscale}} + \texttt{yoffset}$$

while pystog and this app's engine *multiply*:

$$S_\mathrm{corr}(Q) = a\,S_\mathrm{meas}(Q) + b .$$

The reader converts once, at the boundary:

$$a = \frac{1}{\texttt{yscale}},\qquad b = \texttt{yoffset}$$

so the validation run's `-9 0.1` becomes $a = 10$, $b = -9$, i.e. $10\,S - 9 = 10(S-1)+1$ — an
amplitude change that leaves the high-$Q$ asymptote at 1. Every API, plot label, provenance
field, and UI control in this app uses the **multiply** form; `yscale`/`yoffset` appear only in
the parsed `StogInput` and in the provenance `stog_inp_reference` block.

**Rejections are hard, by design.** Non-canonical layouts (multiple datasets, non-zero $Q$
offset, non-zero second $y$ offset, the interactive rescale loop, filter disabled) raise
`NotImplementedError` (Python) / `Error` (JS) rather than being silently misparsed. **A non-zero
$Q$ offset is never applied anywhere in the codebase** — it is parsed, checked, and refused.

**Code:** `parsers.py` → `read_stog_inp()`, `StogInput` (with properties `a = 1/yscale`,
`b = yoffset`); JS `autoScale.js` → `readStogInp()` (same line indices, same five rejections,
also exposes `a`, `b`). Test: `autoScale.test.js` → *"readStogInp decodes the classic layout"*
(asserts $a=10$, $b=-9$, $\rho_0=0.063049$, and that a `Y` try-again flag throws).

**Python ↔ JS divergence — numeric coercion.** The line indices and the rejection set match, the
number parsing does not. Python uses `float()`/`int()`, which **raise** on any field that does
not parse; the JS port uses `Number()`/`parseInt()`, which **coerce silently**. Concretely:
`nr = 5000.0` raises `ValueError` in Python but yields `5000` in the browser (`parseInt` drops
the fractional part, as it does for `n_files`); a garbage `rmax` line raises in Python but
produces `NaN` in JS, deferred until `makeConfig` reports *"rmax must be finite and positive"*.
The browser therefore accepts a few malformed `.inp` files that the reference implementation
refuses, and reports the failure later and less specifically when it does refuse.

**Doc/code note:** the `StogInput` docstring calls it "the 23-line `stog.inp`", while the code
requires ≥ 22 non-empty lines and indexes 0…21. Because empty lines are stripped *before*
indexing, a file with a genuinely blank field (e.g. an empty output-name line) will shift the
whole layout and misparse without warning.

---

### Step 5 — Chemical formula → atomic fractions

**Inputs:** a formula string typed into the **Composition** field (browser) or passed as
`--formula` (CLI). Also accepted programmatically: a `{element: count}` dict.
**Operation:** a small recursive-descent parser supporting decimal stoichiometries and nested
parentheses:

- token: `([A-Z][a-z]?)(\d*\.?\d*)` — element symbol plus optional (possibly fractional) count,
  default 1;
- `(` opens a group, `)` closes it and consumes an optional multiplier applied to every count
  inside;
- repeated elements accumulate (`CH3COOH` → C 2, H 4, O 2);
- unbalanced parentheses, an empty formula, or an unparseable remainder raise `ValueError`.

Counts $N_i$ are converted to fractions by

$$c_i = \frac{N_i}{\sum_j N_j}.$$

Examples pinned by tests: `GaTa4Se8` → {Ga 1, Ta 4, Se 8}; `Sr0.5Ba0.5TiO3` → {Sr 0.5, Ba 0.5,
Ti 1, O 3}; `Al2(SO4)3` → {Al 2, S 3, O 12}.

**Code:** `scattering.py` → `parse_formula()`; JS `autoScale.js` → `parseFormula()`.
Tests: `tests/test_scattering.py::ParseFormulaTests`; `autoScale.test.js` (`Sr0.5Ba0.5TiO3`,
`Al2(SO4)3`). The two parsers agree on every tested form, with one untested divergence: the
count regex `\d*\.?\d*` can match a bare `.`, and the engines then diverge — Python's
`float(".")` raises `ValueError` and the formula is rejected with a clear message, while JS's
`Number(".")` returns `NaN` and stores a NaN count, which propagates through the total, every
fraction, $\langle b\rangle^2$ and $\langle b^2\rangle$, and finally surfaces as `makeConfig`'s unrelated-sounding
*"bAvgSq must be finite and positive"*. This is the same `Number` vs `float` tolerance gap as in
Steps 3 and 4.

*(Implementation nit: a module-level regex `_FORMULA_TOKEN` in `scattering.py` is unused — the
actual parsing is done by the inline `re.match` calls in `parse_formula`.)*

---

### Step 6 — Sears/NIST scattering lengths → the Faber-Ziman coefficients

**Inputs:** the fractions $c_i$ from Step 5, plus optional per-element overrides
(`b_overrides_fm`, e.g. for isotopically enriched samples — Python and the JS function accept
this; the browser UI does not expose it).

**The table.** `COHERENT_B_FM` holds **bound coherent neutron scattering lengths in fm** for 89
natural elements, from the NIST NCNR table (after Sears, *Neutron News* **3**, 26 (1992)). For
the strong absorbers with complex $b$ — B, Cd, Dy, Eu, Gd, In, Sm (listed in
`COMPLEX_B_ELEMENTS`) — **only the real part is stored**, which is the standard choice for
Faber-Ziman weighting but means absorption/resonance effects are not represented. An element
absent from the table raises rather than defaulting to anything.

**Math.** With $b_i$ in fm:

$$\langle b\rangle^2 = \Bigl(\sum_i c_i b_i\Bigr)^2 \quad [\mathrm{fm}^2],
\qquad
\langle b^2\rangle = \sum_i c_i b_i^2 \quad [\mathrm{fm}^2],$$

$$\langle b\rangle^2_{\mathrm{barn}} = \frac{\langle b\rangle^2_{\mathrm{fm}^2}}{100},
\qquad
\langle b^2\rangle_{\mathrm{barn}} = \frac{\langle b^2\rangle_{\mathrm{fm}^2}}{100}
\qquad (1\ \mathrm{barn} = 100\ \mathrm{fm}^2).$$

Both numbers are returned in **both** unit systems on purpose, because the ecosystem is
inconsistent (pystog's bundled argon config quotes `<b_coh>^2 = 3.644` in fm²; a classic
`stog.inp` quotes 0.015407 in barns).

The Python `FaberZiman` dataclass additionally returns the pairwise Faber-Ziman weights

$$w_{ij} = \frac{c_i c_j b_i b_j}{\langle b\rangle^2}, \qquad \sum_{ij} w_{ij} = 1,$$

used for diagnostics and partial-weight reasoning. **The JS port does not compute
`fractions`, `b_coh_fm`, or `weights`** — it returns only the four scalars
(`bAvgSqFm2`, `bSqAvgFm2`, `bAvgSqBarn`, `bSqAvgBarn`). Nothing in the browser pipeline needs
the rest.

**Guard.** A near-null-matrix composition is refused:

$$\langle b\rangle^2 < 10^{-4}\,\langle b^2\rangle \;\Rightarrow\; \mathrm{ValueError}$$

because Faber-Ziman normalization by $\langle b\rangle^2$ is then undefined (e.g. a Ti/Zr null
matrix). Both engines apply the identical test.

**Validation pins (tests/test_scattering.py, test_scaling.py, autoScale.test.js):**

- Mn₃Sn → $\langle b\rangle^2 = 0.015407$ barn — the value the POWGEN 59438 validation run was
  actually scaled with. Note the run's `stog.inp` is **not** in the repository
  (`data/stog_tests/stog_59438/` holds only the data and output files);
  [`tests/test_scaling.py`](../../tests/test_scaling.py) reconstructs the whole parameter set from
  the run's own outputs — the flat $-\langle b\rangle^2$ stretch of `scale_ft_rmc.gr` is what
  pins 0.015407. A comment in `tests/test_scattering.py` describes the number as the expert's
  hand-computed coefficient; that attribution is not verifiable from anything in the repo.
  $\langle b^2\rangle = 0.201223$ barn; ratio $\langle b^2\rangle/\langle b\rangle^2 = 13.06$.
- Ar → $\langle b\rangle^2 = 3.644$ fm² = 0.03644 barn (pystog cross-check).
- Monatomic (Ni) → $\langle b^2\rangle = \langle b\rangle^2$ exactly (Laue term vanishes).
- Polyatomic → $\langle b^2\rangle > \langle b\rangle^2$ (Cauchy–Schwarz), weights sum to 1.

**Code:** `scattering.py` → `COHERENT_B_FM`, `faber_ziman()`, `FaberZiman`; JS `autoScale.js` →
`COHERENT_B_FM`, `faberZiman()`. The two tables are byte-identical value lists and must be kept
in sync manually.

---

### Step 7 — The composition-derived $S(0)$ target

**Inputs:** $\langle b\rangle^2$, $\langle b^2\rangle$ (both barns), or an explicit override.
**Math.** Keen 2001 Eq. 21 gives the $Q\to0$ limit of the normalized structure factor for a
dense solid (the compressibility term is dropped):

$$S(0) = 1 - \frac{\langle b^2\rangle}{\langle b\rangle^2}.$$

**Resolution rule** (`ScalingConfig.effective_s0_target`, JS `effectiveS0Target`):

1. explicit `s0_target` / `s0Target` if set — a field of **both** engines' config objects
   (`ScalingConfig.s0_target`, `defaultConfig.s0Target`), but wired to no CLI flag, no Flask API
   field, and no UI control, so only a direct caller of the Python or JS library can set it (the
   same is true of `fit_offset`/`fitOffset` and `c2_weight`/`c2Weight` in the browser, which
   `AutoStogPage.jsx` never passes to `makeConfig`);
2. else $1 - \langle b^2\rangle / \langle b\rangle^2$ when $\langle b^2\rangle$ is available;
3. else **0.0** — the pystog "solid-state limit" special case.

This target is *not* a fitted quantity at this stage; it enters the analytic correction for the
unmeasured $[0, Q_\mathrm{min}]$ range (`transforms.low_q_correction_basis`, documented in the
transform section), where it changes **only the constant term** of the affine correction basis:
$\mathrm{const}' = (1-S_0)\,\mathrm{const} + S_0\,\mathrm{coef}$ (asserted exactly in
`tests/test_scaling.py::DetectionAndDensityTests::test_s0_target_changes_only_the_constant_term`).
The same $S(0)$ is what the `fz` amplitude criterion pins the extrapolated $Q\to0$ value onto.

**Why it matters:** for negative-$b$ compositions the number is large and negative — Mn₃Sn gives
$S(0) = 1 - 13.06 = -12.06$. Assuming pystog's $S(0)=0$ there introduces an $O(1)$ bias in the
low-$r$ transform; using the composition-aware target cut the low-$r$ residual by ~40 % on the
PG3 runs (recorded in [`docs/SCALING_PROCEDURE.md`](../SCALING_PROCEDURE.md) §3).

The browser shows the value live in the **coefficients chip** ("⟨b⟩² … · ⟨b²⟩ … · S(0) …") and
draws it as a guide line over the first 1.5 Å⁻¹ of the $S(Q)$ plot.

---

### Step 8 — Number density $\rho_0$

$\rho_0$ (atoms/Å³) is load-bearing: the low-$r$ density limit that pins the amplitude is the
line $D(r) = -4\pi\rho_0\langle b\rangle^2 r$, so a wrong $\rho_0$ yields a confidently wrong
scale.

**Resolution chain (browser, `AutoStogPage.jsx` → `resolveConfig`):**

1. the **ρ₀ field** if the user typed one;
2. else the **`stog.inp`** value (line 11), when an `.inp` was loaded;
3. else the data file's **`NUMBER_DENSITY ::`** header;
4. else **mass density + composition**, converted as below (needs both);
5. else, **if $\langle b^2\rangle$ is available**, seed $\rho_0 = 0.05$ Å⁻³ and set
   `wantEstimate = true` — the worker then runs `estimateRho0` first and adopts the
   self-consistent result for the fit (it also refuses to fit if that estimate does not
   converge);
6. else throw: *"number density unknown: set ρ₀ or mass density — or give a composition and
   Auto StoG estimates ρ₀ self-consistently"*.

**The CLI/API chain is not the same list — it branches on mode.** `scaling_cli.py::_build_config`
(and `app.py::_resolve_scaling_config`, which mirrors it) has two disjoint arms:

- **with a `stog.inp`:** `--rho0` → `inp.rho0` (line 11), **and nothing else is consulted**. The
  `NUMBER_DENSITY ::` header and `--mass-density` / `massDensity` live entirely inside the
  `else:` (data-mode) branch and are unreachable here. Passing `--mass-density` alongside a
  `stog.inp` produces no error and no effect — it is silently ignored.
- **in `--data` mode:** `--rho0` → `NUMBER_DENSITY ::` → `--mass-density` + `--formula` (the
  formula is required, else *"--mass-density needs --formula"*) → error.

Neither front-end seeds; the self-consistency route is the explicit `--estimate-rho0` flag
(`scaling_cli.py` → `main`), which is a hard error on non-convergence. The Flask API has **no**
self-consistency path at all.

**Pre-fill is not uniformly non-destructive (browser).** `selectSource` runs the moment a file is
selected and writes straight into the form. It is *not* a "fill the empty fields" pass:

- $Q_\mathrm{min}$/$Q_\mathrm{max}$ **preserve** a value you already typed
  (`current.qmin || extent…`);
- **$\rho_0$ and $r_0$ do not.** The code is `rho0: inp ? String(inp.rho0) : (header?.numberDensity !=
  null ? String(header.numberDensity) : current.rho0)` and `r0: header?.minDistance != null ?
  String(header.minDistance) : current.r0` — no `current.… ||` guard. Whenever the selected data
  file carries `NUMBER_DENSITY ::` / `MINIMUM_DISTANCES ::`, those values **overwrite** whatever
  you had entered;
- every `stog.inp`-derived field — $\langle b\rangle^2$, $r_\mathrm{cut}$, `rmax`, `nr`, Lorch, the enforcement
  cutoff, and the manual $(a, b)$ — is overwritten whenever an `.inp` is loaded, and loading an
  `.inp` forces the **Enforce low-r** toggle on regardless of the user's setting.

This fires on selecting a second file, on re-selecting the same file from the dropdown, and on
**Reset params** (which clears the form and then calls `selectSource` again). So the user does
see these as editable defaults — but a hand-entered $\rho_0$ or $r_0$ can disappear without notice.

**Mass-density conversion (ADDIE convention).**

$$\rho_0 \;[\text{Å}^{-3}] \;=\; \rho_m \;[\text{g/cm}^3] \times \frac{N_A}{10^{24}} \times \frac{n}{M},$$

with $N_A = 6.02214076\times10^{23}$ mol⁻¹ (so the constant `_AVOGADRO_PER_A3 = 0.602214076`),
$M$ the molar mass of **one formula unit** in g/mol and $n$ its atom count. The choice of
formula unit cancels, since both $M$ and $n$ scale with it. The inverse
`mass_density_from_number_density()` is provided and is used by the docs/tests for sanity
statements. Atomic weights come from `ATOMIC_MASS_U` — standard atomic weights (u), CIAAW
values generated from the `periodictable` package for exactly the `COHERENT_B_FM` element set.

Round-trip pinned by `tests/test_scaling.py` and `autoScale.test.js`: Mn₃Sn at
$\rho_0 = 0.063049$ Å⁻³ ↔ $\rho_m = 7.4209$ g/cm³, and back to 9 decimal places.

**Code:** `scattering.py` → `molar_mass()`, `number_density_from_mass_density()`,
`mass_density_from_number_density()`; JS `autoScale.js` → `molarMass()`,
`numberDensityFromMassDensity()`, `massDensityFromNumberDensity()`. Numerically identical
(same constant, same formula); the JS versions throw on non-positive input exactly like Python.

**Self-consistent estimate (referenced here, detailed in the fit section).**
`estimate_rho0` / `estimateRho0` root-find $a_\mathrm{fz}/a_\mathrm{density}(\rho_0) = 1$ by the
fixed-point update $\rho \leftarrow \rho\cdot\mathrm{concordance}$, with `rtol = 1e-3`,
`max_iter = 8`, and clipping to $[10^{-4}, 1]$ Å⁻³. It **requires** $\langle b^2\rangle$, and
sets an `extrapolated` flag when $Q_\mathrm{min} > 1$ Å⁻¹ (the $Q\to0$ extrapolation then owns
the estimate — a starting point, not a measurement). Cross-engine agreement of the *iterated*
result is bounded by the stopping rule (~1e-4 relative), not by single-pass transform precision.

Four behaviours of the estimator are load-bearing and easy to miss:

1. **It overrides your fit architecture.** The first thing it does is
   `replace(config, amplitude_criterion="density", c1_mode="sweep")` (JS: the same spread). An
   `amplitude = 'fz'` or `c1_mode = 'joint'` selection therefore does **not** apply during
   estimation, only to the final fit that follows.
2. **It can stop without converging.** Besides `|concordance − 1| ≤ rtol`, the loop breaks when
   `result.a <= 0 or concordance <= 0` (no density can reconcile the two criteria) and when the
   clipped update lands on the same value (`rho_next == rho`, i.e. pinned at $10^{-4}$ or 1
   Å⁻³). Both exits return `converged: False`.
3. **The returned `rho0` is `history[-1][0]`** — the density the *last fit was run at*, not the
   last update. On a non-converged return you get the last density tried, not a refined one.
4. **Cost.** Each iteration is a full two-pass `autoscale`, so the estimate costs up to
   8 × 2 = 16 pipeline passes before the run's own fit even starts.

The three front-ends treat non-convergence differently: the Web Worker **throws** (the page
refuses to fit on a garbage density), the CLI raises a `CliError`, and the Python
`estimate_rho0` simply returns `converged: False` for the caller to check.

**The standalone "Estimate ρ₀" button (browser only)** is a separate entry point with its own
rules (`AutoStogPage.jsx` → `runEstimate`, gate `canEstimate`). It is enabled whenever a file is
loaded **and** either the composition field is non-empty or $\langle b^2\rangle$ is set — it does **not** require
$\rho_0$ to be unresolved, so it will happily run (and overwrite) an $\rho_0$ you supplied. It calls
`resolveConfig` with the default `mode = 'auto'`, so it can itself fall through to the 0.05 Å⁻³
seed, and posts a `kind: 'estimateRho0'` job.

**Reproducibility caveat.** On success — both from this button and after an auto run that used
the seeding path — the estimated density is written back into the ρ₀ **form field rounded to 5
significant digits** (`String(Number(x.toPrecision(5)))`), while the run that produced it used
the *unrounded* value, which is what lands in the exported provenance (`rho0Used` → the
provenance `config.rho0`). Re-running therefore uses a marginally different density than the one
recorded in the previous run's outputs.

---

### Step 9 — Assembling the configuration: what is required, what is derived

**Inputs:** the form state (or CLI flags), the parsed `stog.inp`, the `.dat` header, and the
composition.
**Operation:** one precedence pass, then eager validation.

**Precedence — there is no single global ordering.** `AutoStogPage.jsx::resolveConfig`,
`scaling_cli.py::_build_config`, and `app.py::_resolve_scaling_config` agree *quantity by
quantity* (table below), but the often-quoted chain "explicit → `stog.inp` → header →
composition → default" is wrong in three ways:

- **$r_0$ takes the data header *before* the `stog.inp`**, the reverse of that rule
  (`_default_r0`: `--r0` → `header['min_distance']` → the `.inp` peak line; `resolveConfig` does
  the same).
- **$\langle b^2\rangle$ never consults the `.inp` or the header** — it comes from the
  form/payload or the composition, nothing else.
- **The CLI and the API branch on mode; the browser does not.** With a `stog.inp` loaded, the
  `.inp` fields are the *only* fallback in the Python front-ends: the `NUMBER_DENSITY ::`
  header, `--mass-density`/`massDensity`, and the hard-coded defaults for
  $r_\mathrm{cut}$/`rmax`/`nr`/Lorch are all unreachable (only $r_0$ still reads
  `MINIMUM_DISTANCES ::`). The browser's `resolveConfig` runs one unbranched chain that consults
  the header and the mass density regardless, so it is slightly more permissive than the
  reference implementation.

The path the header is read *from* also differs: the browser matches `inp.dataFile` against the
uploaded set, the CLI resolves `inp_path.parent / inp.data_file`, and the API additionally
rejects a data file that resolves outside the configured data roots.

Per-quantity (browser column; read the CLI/API as "the `.inp` arm only, when an `.inp` is
given"):

| Quantity | Precedence | Hard default |
| --- | --- | --- |
| $Q_\mathrm{min}$, $Q_\mathrm{max}$ | form → `stog.inp` (browser also pre-fills from the data's finite extent) | **none — required** |
| $\langle b\rangle^2$ (barn) | form → `stog.inp` line 18 → composition | **none — required** |
| $\langle b^2\rangle$ (barn) | form → composition (never the `.inp`, never the header) | `None` (then $S(0)$ target = 0, no FZ amplitude, no $\rho_0$ estimate, no $Q\to0$ diagnostic) |
| $\rho_0$ | see Step 8 | seed 0.05 + self-consistent estimate (browser only) |
| $r_\mathrm{cut}$ | form → `stog.inp` line 15 | 1.0 Å |
| $r_0$ | form → **header `MINIMUM_DISTANCES`** → `stog.inp` $\max(\text{peak\_cutoff}, \text{peak\_rmin})$, but only if $r_0 - 0.25 > r_\mathrm{cut} + 0.2$ | `None` → detected from the data after the first pass |
| `rmax`, `nr` | form → `stog.inp` | 50.0 Å, 5000 points |
| Lorch | form → `stog.inp` flag | off |
| enforcement `(cutoff, peak_rmin, peak_rmax)` | resolved *outside* the config — see below | `'auto'` → the detected $r_0$ |

**Enforcement is a triple, and it does not travel in the config.** `ScalingConfig.enforce_cutoff`
exists and drives `transforms.enforce_low_r` (a hard $G_K(r) = -\langle b\rangle^2$ for all
$r \le$ cutoff), but **no shipping path ever sets it**: `_build_config` does not, the API's
`_resolve_scaling_config` does not, and the JS `defaultConfig` has no `enforceCutoff` key at all.
It is a library-only knob. All three front-ends instead resolve a separate three-value tuple
alongside the config and apply `transforms.first_peak_zero` *after* the fit — in the CLI writer
(`_write_outputs`) and in the Web Worker. `first_peak_zero` zeroes $g(r)$ where
$r \le \mathrm{cutoff}$ **and** ($r \ge \text{peak\_rmax}$ **or** $r \le \text{peak\_rmin}$), i.e.
it can preserve a first-peak window that starts inside the cutoff; `enforce_low_r` is the special
case in which the peak window lies entirely outside $[0, \mathrm{cutoff}]$.

How the triple is resolved (`AutoStogPage.jsx::resolveEnforcement`,
`scaling_cli.py::_resolve_enforcement`, `app.py::_resolve_scaling_enforcement`):

- `cutoff`: form **Cutoff** / `--enforce-cutoff` → `stog.inp` line 22 `peak_cutoff` → `'auto'`,
  which the worker (and the CLI's post-run block) resolves to the **detected $r_0$**.
- `peak_rmin`, `peak_rmax`: taken from the `.inp` **only** when the user supplied no explicit
  cutoff. A user-typed cutoff, and the `'auto'` path, both collapse the window to
  `peak_rmin = peak_rmax = cutoff` — a flat replacement of everything below the cutoff.
- The CLI exposes `--peak-window RMIN RMAX` (same semantics) and rejects
  `--no-enforce` combined with `--enforce-cutoff`/`--peak-window`. The browser has no
  peak-window control.
- **Browser wrinkle:** because `selectSource` pre-fills the Cutoff field with `inp.peakCutoff`
  (Step 8), `numberOr(form.enforceCutoff)` is *defined* whenever an `.inp` is loaded, so
  `usingInpWindow` is false and the `.inp`'s `peak_rmin`/`peak_rmax` are **not** used — the
  browser degenerates to a flat replacement where the CLI would keep the peak window. This is
  invisible on the validation runs (where `peak_rmin > cutoff`, so both forms zero exactly
  $r \le$ cutoff) and only bites when `peak_rmin < cutoff`. Clear the Cutoff field to get the
  `.inp` window back.

**So the genuinely required inputs are:** an $S(Q)$ file, the $[Q_\mathrm{min}, Q_\mathrm{max}]$
window, $\langle b\rangle^2$, and $\rho_0$ — where a **chemical composition alone supplies both
$\langle b\rangle^2$ and $\langle b^2\rangle$**, and $\langle b^2\rangle$ in turn makes $\rho_0$
recoverable from the data. That is the sense in which
[`docs/STOG_SCALING_PLAN.md`](../STOG_SCALING_PLAN.md) §5c says the ρ₀ resolution chain ends in
the self-consistent estimate, "making composition + Q window the only required inputs"
([`docs/SCALING_PROCEDURE.md`](../SCALING_PROCEDURE.md) §1 says the same thing as a table, with
composition and $[Q_\mathrm{min}, Q_\mathrm{max}]$ the only two "Required? **yes**" rows).
Strictly, the composition is *not* required if you supply
$\langle b\rangle^2$ yourself (from a `stog.inp` or the Advanced panel) — that is exactly the
x-ray path in Step 12.

**Q-window pre-fill (browser only).** When no `stog.inp` is present and the Q fields are empty,
`selectSource` fills them from the **finite** extent of the uploaded data (`dataExtent` walks
the rows and records the first and last index where both $Q$ and $S$ are finite — this is what
makes NaN padding invisible to the user), rounded to **4 significant digits**
(`Number(x.toPrecision(4))`). Rounding can move the bound by up to half a unit in the 4th digit,
which may include or exclude one edge point relative to the raw file extent. These are
suggestions; the user is expected to set the physics window (see the Caveats).

**Malformed numbers are silently discarded, not rejected (browser).** Every numeric form field
goes through `numberOr`, which returns `undefined` **both** for an empty field **and** for
anything failing `Number.isFinite`. An unparseable entry — `0.05x`, `1,5`, `2 .4`, or a stray
unit like `0.063 A^-3` — is therefore treated as *not supplied* and falls through to the next
precedence level with no "could not parse" message anywhere in the UI. A typo in $\rho_0$ quietly
drops you onto the `stog.inp` value, the header value, the mass-density route, or the 0.05 seed
plus self-consistent estimate; a typo in a $Q$ bound falls back to the `.inp` window or raises
the generic *"set the Q window (Qmin and Qmax)"*. If a run behaves as though a field were empty,
check the resolved values in the exported provenance JSON.

**A fixed-$(a,b)$ run silently forces the amplitude criterion.** `resolveConfig` takes a `mode`
argument and sets `amplitudeCriterion: mode === 'manual' ? 'density' : form.amplitude`, so a
leftover `fz` selection cannot reject an otherwise valid manual run (the `fz` validation needs
$\langle b^2\rangle$ and `c1_mode='sweep'`). The CLI implements the same intent *loudly*: it errors with
*"--amplitude selects the auto-fit criterion; it cannot be combined with
--manual/--scale/--offset"*, and likewise for `--estimate-rho0` + `--manual`.

**Validation (eager, before any compute).** `ScalingConfig.__post_init__` / `makeConfig` reject:
`c1_mode ∉ {sweep, joint}`; `amplitude_criterion ∉ {density, fz}`; `amplitude_criterion == 'fz'`
without $\langle b^2\rangle$ or with `c1_mode != 'sweep'`; non-finite or non-positive $\rho_0$;
non-finite or non-positive $\langle b\rangle^2$; $Q_\mathrm{max} \le Q_\mathrm{min}$;
non-integer or non-positive `nr`; non-finite or non-positive `rmax`. The CLI, the API, and
`makeConfig` additionally evaluate the low-$r$ fit window immediately, so an empty window
$[\,r_\mathrm{fit,lo}, r_\mathrm{fit,hi}\,]$ fails fast with a clean message instead of surfacing
mid-fit.

Two validation asymmetries:

- **NaN $Q$ bounds.** Python tests `if self.qmax <= self.qmin`, which is `False` when either is
  NaN — so a config with NaN $Q$ bounds passes `__post_init__` and fails later. JS tests
  `if (!(config.qmax > config.qmin))`, which throws on NaN.
- **Fractional `nr`.** The browser never lets one reach validation: `resolveConfig` does
  `nr: Math.round(pick(form.nr, …))`, so a typed `4999.6` is silently snapped to 5000 and
  $\Delta r = r_\mathrm{max}/n_r$ changes accordingly. `ScalingConfig` would raise
  (*"nr must be a positive integer"*) and argparse rejects `--nr 4999.6` outright.

**Derived grids and windows** (properties, not stored arrays):

$$r_k = k\,\frac{r_\mathrm{max}}{n_r},\quad k = 1\ldots n_r
\quad\Rightarrow\quad \Delta r = 0.01\ \text{Å for the defaults (50 Å / 5000)},$$

$$Q_\mathrm{tail,min} = Q_\mathrm{max} - f_\mathrm{tail}\,(Q_\mathrm{max} - Q_\mathrm{min}),
\qquad f_\mathrm{tail} = 0.15,$$

$$[r_\mathrm{fit,lo},\, r_\mathrm{fit,hi}] =
\bigl[\,r_\mathrm{cut} + 0.2,\; r_0 - 0.25\,\bigr]
\quad\text{(or } r_\mathrm{fit,lo} + 1.0 \text{ when } r_0 \text{ is unknown)} .$$

Note $r = 0$ is **excluded** from the grid (the $g(r)$ conversions divide by $r$). Both fit
windows must contain at least 2 points or `_fit_windows` / `fitWindows` raises.

**The C1 window is anchored on the *configured* $Q_\mathrm{max}$, not on the last surviving data
point.** `q_tail_min` is a pure property of the config; the tail mask is then `q >= q_tail_min`
over whatever the crop left. So setting $Q_\mathrm{max}$ beyond the end of the measured data
silently **shrinks** the C1 block — in the limit below 2 points, which surfaces as the generic
*"fit windows contain fewer than 2 points"* with no hint that the $Q$ window is the cause. (The
4-significant-digit pre-fill makes overshooting easy, and the field is freely editable.)
Conversely a $Q_\mathrm{max}$ well inside the data leaves C1 correctly sized but discards
measurements. The level sweep behaves differently — it operates on the surviving data, not the
configured window — so the two can disagree about where "high $Q$" is. Keep $Q_\mathrm{max}$
inside the measured range.

**Code:** `scaling.py` → `ScalingConfig` (+ `r_grid`, `q_tail_min`, `r_fit_window`,
`effective_s0_target`), `_fit_windows()`; JS `autoScale.js` → `defaultConfig`, `makeConfig`,
`rGrid`, `qTailMin`, `rFitWindow`, `fitWindows`. Test:
`autoScale.test.js` → *"rejects invalid configurations like the Python engine"*.

---

### Step 10 — Cropping $Q$ to $(0,\; Q_\mathrm{min}\ldots Q_\mathrm{max}]$

**Inputs:** the raw `q`, `sq`, optional `sigma`, and the config.
**Operation:** a single boolean keep-mask, applied to all columns together:

$$\mathrm{keep}_i = \operatorname{isfinite}(Q_i)\;\wedge\;\operatorname{isfinite}(S_i)
\;\wedge\; Q_i > 0
\;\wedge\; Q_i \ge Q_\mathrm{min} - 10^{-12}
\;\wedge\; Q_i \le Q_\mathrm{max} + 10^{-12}.$$

This is where the **NaN padding** of rebinned files is finally removed, together with any
$Q \le 0$ row. The $Q>0$ test is unconditional and deliberate: every $S(Q)\leftrightarrow F(Q)$
conversion divides by $Q$, and the analytic low-$Q$ correction already models the omitted
$[0, Q_\mathrm{min}]$ range, so a $Q=0$ point carries no usable information (and
`fourier_filter` refuses a grid whose first point is $\le 0$).

If fewer than **16** points survive, both engines raise *"fewer than 16 usable S(Q) points after
cropping"*.

**The 16-point floor is checked *before* despiking and never re-checked.** `crop_sq`/`cropSq`
raise on `keep.sum() < 16`, then apply the despike mask unconditionally and return whatever is
left (Step 11). Despiking can therefore push the dataset below 16 points, and the failure then
surfaces downstream with no mention of despiking — as *"level_sweep needs at least 32 finite
points"* or *"fit windows contain fewer than 2 points"*.

**The stricter gate is 32, not 16.** With the default `c1_mode = 'sweep'`, `level_sweep` /
`levelSweep` require **at least 32 finite points** (and each candidate window needs ≥ 24 points
and a width ≥ 3 Å⁻¹). 16 is only the binding limit in `joint` mode or in a fixed-$(a,b)$ run.

The $10^{-12}$ Å⁻¹ tolerance on both bounds means a user-typed bound that matches a grid point
to within float noise still includes that point.

**Everything low-$Q$ is anchored on `q[0]` — the first *surviving* point, not on your
$Q_\mathrm{min}$.** After the crop, `transforms.low_q_correction_basis` builds the analytic
$[0, q_0]$ correction from `q0 = q[0]`; `_solve_affine` uses `sq[0]`/`delta_sq[0]` as the
$S(Q_\mathrm{min})$ term of the affine basis; `fourier_filter` refuses a grid with
`q[0] <= 0`; and `amplitude_from_fz_limit` fits the head `q <= q[0] + 1.0` Å⁻¹. A one-point
change at the low-$Q$ edge — from the $10^{-12}$ tolerance, from the 4-significant-digit
pre-fill rounding, or from a NaN padding row — moves that anchor and shifts both the analytic
correction and the $Q\to0$ Faber-Ziman extrapolation.

**Outputs:** cropped `(q, sq, sigma)` in original order.
**Code:** `scaling.py` → `crop_sq()`; JS `autoScale.js` → `cropSq()`. Behaviourally identical.
`sigma` is carried along by the same mask but is **not** itself tested for finiteness here.

---

### Step 11 — Despiking (optional, OFF by default)

**Inputs:** the cropped `sq`.
**Operation:** a rolling-median outlier rejection applied *after* the crop and *before* any
transform.

1. Rolling median $m_i$ over a window of `despike_window = 7` points (`pad = window // 2 = 3`),
   with **edge replication** at the array ends.
2. Residual $d_i = S_i - m_i$.
3. Robust scale $\mathrm{MAD} = 1.4826 \cdot \operatorname{median}_i |d_i|$ — note this is a
   median of absolute residuals about **zero**, not about the residual median.
4. Keep-mask $|d_i| \le n_\sigma \max(\mathrm{MAD}, 10^{-12})$ with `despike_nsigma = 6.0`.

**Python despikes twice in the auto path; JS despikes once — the engines diverge whenever
`despike=True`.** Despiking lives inside `crop_sq`, and the Python auto path calls `crop_sq`
twice: `_autoscale_pass` rebinds `q, sq, sigma = crop_sq(q, sq, config, sigma)` and then hands
those **already cropped-and-despiked** arrays to `scale_pipeline`, which calls `crop_sq` again.
The second pass recomputes the rolling median and the MAD on cleaned data, so the MAD is
smaller, the $6\sigma$ bound is tighter, and additional points are dropped. The JS
`autoscalePass` instead calls `scalePipeline(qIn, sqIn, …)` with the **original** arrays, so the
crop+despike is recomputed identically and exactly one despike pass is applied. Two consequences:

- with `despike=True` the two engines return **different point sets**, so cross-engine parity
  holds only for `despike=False` (the default, and the only setting the parity fixtures cover);
- in Python the arrays written to the output files (`q`, `sq_scaled`, `fk`, …) are a strict
  *subset* of the arrays the fit actually used.

**`n_despiked` / `nDespiked` is not the total in either engine.** Python computes it inside
`scale_pipeline` as (points before *that* despike) − (points after it), so in the **auto** path
it counts only the extra points removed by the redundant second pass and therefore
**under-reports**; it is the true count only in the manual path (`scale_pipeline` called
directly). In the browser, `cropSq` computes `nDespiked` and `autoScaleWorker.js` returns it —
but `AutoStogPage.jsx` never reads it: `setPreview` copies only
`{a, b, converged, iterations, lowRRms, c1TailMean, history, c1ModeEffective}`, nothing renders
it, and `writeFiles` omits it from the exported provenance JSON. **There is currently no way for
a browser user to see how many points despiking removed**; the count is only visible from the
CLI/API provenance, with the auto-path caveat above.

**Why it is off by default:** despiking measurably restores clean scale recovery under
detector-glitch spikes — which otherwise ring through the sine transform into the low-$r$ fit
window, a channel the Huber IRLS re-weighting cannot reject — but it **also flags real Bragg
maxima on crystalline data**: the `ScalingConfig` docstring (and a test comment) record 12 % of
points on the 59438 benchmark, a figure no test or artifact in this repo reproduces. Enable it
only for glitch-type contamination.
Test: `tests/test_scaling.py::test_despike_restores_recovery_under_tail_glitches` asserts only
that the despiked fit recovers $a$ to 2 %, beats the spiked fit by ≥ 5×, and reports
`n_despiked >= 8`.

**Code:** `scaling.py` → `_despike_mask()` (uses `numpy.lib.stride_tricks.sliding_window_view`),
called from `crop_sq()`; JS `autoScale.js` → `despikeKeepMask()` (explicit clamped-index buffer).
The two mask *functions* produce the same mask for odd windows — it is the surrounding pipelines
that differ. **An even `despike_window` is not supported**: the NumPy path would raise on a shape
mismatch while the JS path would silently use an asymmetric window. The UI exposes only the
on/off toggle, so window and $n_\sigma$ are always 7 and 6.0 in the browser.

---

### Step 12 — The σ column and the x-ray / non-neutron path

#### σ (third column)

If the data file has a third column it is taken as a per-point uncertainty and used to
$1/\sigma$-weight the high-$Q$ (C1) rows of the fit, normalized to unit mean over the tail block:

$$w_i = \frac{1/\max(\sigma_i, 10^{-12})}{\overline{\left(1/\max(\sigma, 10^{-12})\right)}}.$$

Toggle: **σ column** in the Advanced panel (default on) / `--sigma/--no-sigma` (default on) /
`useSigma` in the API payload. σ never affects the low-$r$ (C2) rows.

Note $\max(\sigma, 10^{-12})$ is a **lower clamp, not a rejection**. A zero or negative σ — common
filler in reduced files — is clamped to $10^{-12}$ and produces a weight of order $10^{12}$. That
single row then dominates the entire C1 block and effectively pins the high-$Q$ level to it: no
NaN, no error, just a silently distorted fit.

**Divergence — the browser skips the CLI's validity gate.** `scaling_cli.py::_load_dataset` and
`app.py::_cached_scaling` both drop the σ column entirely if any σ on a usable row is
non-finite or $\le 0$ ("a broken uncertainty column must not poison the fit"), so the Python
front-ends proceed unweighted rather than distorted. `AutoStogPage.jsx` takes `columns[2]` as-is
when it exists, and `cropSq` masks σ only on the finiteness of $Q$ and $S$. In the browser,
therefore, a σ column with `NaN` on rows where $S$ is finite produces `NaN` weights and a
failed/garbage solve, and a σ column with zeros or negatives produces the $10^{12}$-weight
distortion above — neither is reported.

#### X-ray / electron data

**There is no form-factor path.** No $f(Q)$, $Z$, Compton, or electron scattering-factor table
exists anywhere in `rmc_toolkits/` or in `autoScale.js`; the Sears table is neutron-only and the
correction model $S_\mathrm{corr} = aS + b$ is $Q$-independent by construction. The supported
x-ray workflow is entirely manual, documented in the UI tooltips and
[`docs/SCALING_PROCEDURE.md`](../SCALING_PROCEDURE.md) §5:

- set $\langle b\rangle^2$ explicitly — **1** for a normalized x-ray $S(Q)$ — in
  **Advanced → Coefficients** (or `--b-avg-sq`);
- set $\langle b^2\rangle = \langle Z^2\rangle/\langle Z\rangle^2$ explicitly (using
  $f(0) = Z$), which is what supplies the $S(0)$ target and the FZ amplitude;
- leave the composition field empty, or accept the warning chip.

The page guards the classic stale-session trap (Step 1's `sessionStorage` persistence is what
creates it): if a composition *and* an Advanced override are both present and disagree by more
than **2 %**, the coefficients chip turns into a warning ("⚠ … Advanced overrides"). The CLI
prints an analogous warning to stderr and keeps the configured value — but **the two tests are
not equivalent**:

| | Browser (`AutoStogPage.jsx`, `coefficients` memo) | CLI (`scaling_cli.py::_build_config`) |
| --- | --- | --- |
| Quantities checked | **both** $\langle b\rangle^2$ and $\langle b^2\rangle$ | **$\langle b\rangle^2$ only** — a contradictory $\langle b^2\rangle$ override is never flagged |
| Test | $\lvert x_\mathrm{override} - x_\mathrm{fz}\rvert > 0.02\,\lvert x_\mathrm{fz}\rvert$ (relative to the **composition** value) | $\lvert \langle b\rangle^2_\mathrm{fz} - \langle b\rangle^2_\mathrm{config}\rvert > 0.02\,\lvert \langle b\rangle^2_\mathrm{config}\rvert$ (relative to the **configured** value) |
| When it can fire | whenever a parseable composition and an override coexist | only when `--formula` is combined with an already-resolved $\langle b\rangle^2$ — i.e. `--b-avg-sq` **or** a `stog.inp` line 18 |

Because the denominators differ, the two disagree near the boundary. The CLI's version also
fires in `stog.inp` mode, comparing `--formula` against the `.inp`'s hand-entered line-17 value.
The Flask API performs neither check.

The only nod to x-ray/Placzek drift in the engine is the Python-only, experimental
`c1_slope_nuisance` flag (an extra $m\,(Q - \bar Q)$ column on the C1 rows). It is **not** in the
JS engine, **not** exposed by the CLI or the UI, and its own docstring notes that a drift
spanning the whole $Q$ range also enters through the transform, which it does not correct.
`docs/STOG_SCALING_PLAN.md` §5 still states "neutron / constant-b only in v1" and labels x-ray
use out of scope; the FeCoSn x-ray runs have since been used as validation data with
hand-supplied coefficients, so read that scope note as "no automatic x-ray corrections", not
"x-ray files are rejected".

---

### Step 13 — What is handed to the fit

At the end of Step 0 the engine holds:

- `q`, `sq` (float64, cropped, optionally despiked, ascending by assumption), plus `sigma`
  or `None`;
- the $r$ grid $r_k = k\,r_\mathrm{max}/n_r$;
- the tail mask $Q \ge Q_\mathrm{tail,min}$ and the low-$r$ window mask
  $r \in [r_\mathrm{cut}+0.2,\; r_0-0.25]$;
- the scalar physics constants $\langle b\rangle^2$, $\langle b^2\rangle$, $S(0)$, $\rho_0$;
- the flags (Lorch, low-$Q$ correction, robust, `c1_mode`, `amplitude_criterion`, enforcement).

Most of this is recorded in the provenance JSON that ships with the outputs — but **there are two
differently-shaped provenance documents, and they are not field-comparable.**

**Python.** `scale_pipeline()` builds `result.provenance` =
`{model, mode, a, b, config, q_tail_window, r_fit_window, n_q_points, n_despiked}`, plus
`c1_mode_effective` and `level_sweep` (and `r0_detected`/`window_refined` when detection ran) in
auto mode. `provenance.config` is a **fixed 26-key tuple** of `ScalingConfig` attributes; note it
omits `r_fit_min`/`r_fit_max`, so an explicit fit-window override is *not* recorded there — only
the resolved `r_fit_window` is. `scaling_cli.main()` then wraps that dict in a larger payload
carrying `tool`, `rmc_toolkits_version`, `argv`, `stog_inp`, `data_file`, `stog_inp_reference`,
`enforcement`, `outputs`, `rho0_estimate`, `diagnostics`, and `provenance` — so
`stog_inp_reference`, `rho0_estimate` and `diagnostics` are **already present in the CLI JSON**,
as siblings of `provenance`.

**Browser.** `AutoStogPage.jsx::writeFiles` builds a different, flatter document:
`{tool: "rmc-autoscale (browser engine)", source, a, b, mode, c1ModeEffective, stogInpReference,
history, enforcement, rho0Estimate, config, diagnostics}`. Its `config` is the **camelCase JS
config object** (the 25 keys of `defaultConfig`) — it *includes* `rFitMin`/`rFitMax`, which
Python's tuple omits, and *lacks* `c2Bins`, `c1SlopeNuisance` and `enforceCutoff`, which have no
JS counterpart. It carries **no** `q_tail_window`, **no** `n_q_points`, and **no** `n_despiked`
(see Step 11), records `r0_detected`/`window_refined` only inside `diagnostics`, and has no
version token — the literal string `"rmc-autoscale (browser engine)"` stands in for one.

Only the fit `history` is genuinely browser-only: Python's `scale_pipeline` provenance dict has
no `history` key (the `ScalingResult` carries it as a field, but it is not serialized into the
provenance JSON).

---

### Parameters and defaults

| Parameter | Default | Where set | Units / notes |
| --- | --- | --- | --- |
| `qmin`, `qmax` | **required** | form / `--qmin`,`--qmax` / `stog.inp` line 3 | Å⁻¹ |
| `b_avg_sq` $\langle b\rangle^2$ | **required** (composition or explicit) | form / `--b-avg-sq` / `--formula` / `stog.inp` line 18 | barn |
| `b_sq_avg` $\langle b^2\rangle$ | `None` | form / `--b-sq-avg` / `--formula` | barn; enables FZ amplitude, $\rho_0$ estimate, S(0) target, Q→0 diagnostic |
| `rho0` | **required** (or estimated) | form / `--rho0` / `--mass-density` / header / `stog.inp` line 11 | Å⁻³ |
| `s0_target` / `s0Target` | `None` → $1-\langle b^2\rangle/\langle b\rangle^2$ → `0.0` | Python/JS config only (not exposed) | dimensionless |
| `r_cutoff` | `1.0` | form / `--r-cutoff` / `stog.inp` line 15 | Å |
| `r0` | `None` → header → `stog.inp` → detected | form / `--r0` / `MINIMUM_DISTANCES ::` | Å |
| `r_fit_min` / `r_fit_max` | `r_cutoff + 0.2` / `r0 − 0.25` (or `lo + 1.0`) | form / `--r-fit-min`,`--r-fit-max` | Å |
| `rmax` | `50.0` | form / `--rmax` / `stog.inp` line 8 | Å |
| `nr` | `5000` | form / `--nr` / `stog.inp` line 9 | → Δr = 0.01 Å at defaults |
| `q_tail_frac` | `0.15` | Python/JS config only (not exposed) | fraction of the Q window |
| `lorch` | `False` | form toggle / `--lorch` / `stog.inp` line 10 | — |
| `low_q_correction` | `True` | form toggle / `--low-q-correction` | disable only for strict classic-stog parity |
| `robust` (Huber IRLS) | `True` | form toggle / `--robust` | Huber $c = 1.345$, 3 IRLS passes |
| `despike` | `False` | form toggle / `--despike` | — |
| `despike_window` | `7` | Python/JS config only | points; must be odd |
| `despike_nsigma` | `6.0` | Python/JS config only | × MAD |
| `c1_mode` | `"sweep"` | form select / `--c1-mode` | `sweep` \| `joint` |
| `amplitude_criterion` | `"density"` | form select / `--amplitude` | `density` \| `fz` |
| `fit_offset` | `True` | `--fit-offset` (joint mode only); Python/JS config only in the browser | superseded in sweep mode |
| `c2_weight` | `1.0` | Python/JS config only | — |
| `c2_bins` | `0` (pointwise) | **Python only** | — |
| `c1_slope_nuisance` | `False` | **Python only**, experimental | — |
| `max_iter`, `tol` | `50`, `1e-6` | Python/JS config only | self-consistent loop |
| enforcement `(cutoff, peak_rmin, peak_rmax)` | browser default "Enforce low-r" ON → `'auto'` (detected $r_0$) | form Cutoff / `--enforce-cutoff` + `--peak-window` / `stog.inp` line 22 | Å; resolved *outside* the config, applied with `first_peak_zero` |
| `enforce_cutoff` (config field) | `None` | **library only** — never set by the CLI, the API or the browser | Å; drives `enforce_low_r` |
| min. usable points after crop | `16` (checked *before* despiking) | hard-coded | both engines |
| min. finite points for the level sweep | `32` (window ≥ 24 pts, ≥ 3 Å⁻¹ wide) | hard-coded | the binding limit in the default `sweep` mode |
| $\rho_0$ estimate `rtol`, `max_iter`, bounds | `1e-3`, `8`, `[1e-4, 1]` | hard-coded | Å⁻³; each iteration is a full two-pass `autoscale` |
| $\rho_0$ seed when unresolved | `0.05` | browser only | Å⁻³ |
| $\rho_0$ write-back rounding | 5 significant digits | browser only | run uses the unrounded value |
| formula-vs-override warning threshold | 2 % (of the *composition* value in the browser, of the *configured* value in the CLI) | browser chip ($\langle b\rangle^2$ and $\langle b^2\rangle$) + CLI stderr ($\langle b\rangle^2$ only) | — |

### Python ↔ JavaScript parity for Step 0

| Component | Python | JS | Agreement |
| --- | --- | --- | --- |
| Formula parsing | `parse_formula` | `parseFormula` | exact on all tested forms |
| Sears table / $\langle b\rangle^2$, $\langle b^2\rangle$ | `faber_ziman` | `faberZiman` | exact (same table, same arithmetic); JS omits `fractions`/`weights`/`b_coh_fm` |
| Mass ↔ number density | `number_density_from_mass_density` | `numberDensityFromMassDensity` | exact (same $N_A$ constant) |
| `S(0)` target | `effective_s0_target` | `effectiveS0Target` | exact |
| Column reader | `read_stog_xy` | `readStogXy` | same rules **except** Fortran `D` exponents (JS accepts, Python skips the row) |
| `stog.inp` reader | `read_stog_inp` | `readStogInp` | same line indices, same five rejections; Python raises on any field that is not a valid float/int, JS coerces (`Number`/`parseInt`) and defers NaN to `makeConfig` — and accepts non-integer `nr`/`n_files` |
| `.dat` header | `read_dat_header` | `readDatHeader` | matches on realistic headers; diverges on `nan`/`inf` tokens (Python accepts, JS skips), on underscore vs hex literals, on an empty `TITLE ::`, and on file decoding — see Step 3 |
| Q crop | `crop_sq` | `cropSq` | exact |
| Despike mask | `_despike_mask` | `despikeKeepMask` | exact for odd windows; even windows differ (NumPy raises) |
| Despike **pipeline** | applied **twice** in the auto path | applied **once** | **divergent** whenever `despike=True` — see Step 11 |
| `n_despiked` reporting | honest only in manual mode | computed and returned, then dropped by the page | neither is trustworthy end-to-end |
| σ validity gate | CLI/API drop a broken column | **absent** in the browser page | divergent — see Step 12 |
| Config validation | `ScalingConfig.__post_init__` | `makeConfig` | same checks, same messages in spirit, except that JS's `!(qmax > qmin)` also rejects NaN $Q$ bounds that Python's `qmax <= qmin` lets through; both evaluate the low-$r$ fit window eagerly in the shipping paths |
| Enforcement | `first_peak_zero` in `_write_outputs` | `firstPeakZero` in the worker | same function; the browser cannot reach the `.inp` peak window (Step 9) |

Overall engine parity is asserted numerically: $(a,b)$ to 1e-6 relative, level to 1e-9, sampled
$G_K(r)$/filtered $S(Q)$ to 9 decimals, $\rho_0$ estimate to 1e-4 relative
([`src/__tests__/autoScale.test.js`](../../web_app/frontend/src/__tests__/autoScale.test.js)).

### Caveats / what this is not

- **The Q window is a physics choice, not a formatting detail.** Auto StoG crops to
  $(0, Q_\mathrm{min}\ldots Q_\mathrm{max}]$ and everything downstream — the analytic low-$Q$
  correction, the level sweep, the FZ extrapolation — is defined relative to that window. The
  browser's pre-filled values are the data's finite extent rounded to 4 significant digits, not
  a recommendation. Cut $Q_\mathrm{max}$ before detector rolloff; push $Q_\mathrm{min}$ as low
  as the reduction allows.
- **No monotonicity, uniqueness, or unit checking of $Q$.** Rows are kept in file order. A
  non-monotonic or duplicated $Q$ column will silently corrupt every trapezoid transform.
  Likewise `NUMBER_DENSITY ::` and `MINIMUM_DISTANCES ::` are read as bare numbers with the
  units in the line ignored.
- **The scattering-length table is neutron, natural-abundance, real-part-only.** Isotopic
  samples need `b_overrides_fm`, the second argument of `faber_ziman()` / `faberZiman()` — a
  **library-only** knob with no CLI flag (`build_parser` exposes only `--formula`,
  `--b-avg-sq`, `--b-sq-avg`), no Flask API field, and no browser control; the
  seven complex-$b$ absorbers keep only $\mathrm{Re}(b)$; an element outside the 89-element
  table raises rather than guessing.
- **$\langle b\rangle^2$ and $\langle b^2\rangle$ are different numbers** and are easy to swap. $\langle b\rangle^2$ is the classic stog
  "Faber-Ziman coefficient" and the low-$r$ level $-\langle b\rangle^2$; $\langle b^2\rangle$ only sets the
  $Q\to0$ limits. Both are quoted in barns internally, but pystog configs quote fm² in places —
  a factor-100 trap.
- **Without $\langle b^2\rangle$ the pipeline quietly assumes $S(0)=0$.** That is right for monatomic/dense
  solids and badly wrong for negative-$b$ compositions. Supply a composition (or an explicit
  $\langle b^2\rangle$) whenever the contrast is not simple.
- **A wrong $\rho_0$ produces a confidently wrong scale**, because the C2 target line scales with it.
  The self-consistent estimate is not a substitute for a measured density when
  $Q_\mathrm{min} \gtrsim 1$ Å⁻¹ (it is flagged `extrapolated` there), and it refuses to run
  without $\langle b^2\rangle$.
- **Despiking is not a general cleanup, and its accounting is unreliable.** It removes narrow
  rolling-median outliers, which on crystalline data includes genuine Bragg peaks. It is off by
  default. When on: the Python auto path despikes twice and reports only the second pass's
  count, the browser reports nothing at all, the two engines no longer produce the same point
  set, and the 16-point floor is not re-checked afterwards (Steps 10–11). Treat `n_despiked` as
  a lower bound, and prefer the CLI in manual mode if you need the real number.
- **The browser's σ handling is less defensive than the CLI's** (Step 12). A partially invalid
  uncertainty column gives NaN weights (from NaN σ) or a single row with a $10^{12}$ weight
  (from a zero or negative σ) instead of being dropped. Either clean the column or turn the
  **σ column** toggle off.
- **Typos in the form are silently ignored, not rejected** (Step 9). An unparseable number is
  indistinguishable from an empty field, so a mistyped $\rho_0$ or $Q$ bound quietly falls through to
  the next fallback. Verify the resolved values in the provenance JSON.
- **Parameters persist across reloads; uploads do not** (Step 1). A `sessionStorage` restore can
  carry a previous sample's composition, $\rho_0$ and coefficient overrides into a new file. Watch for
  the ⚠ coefficients chip, and use **Reset params** when switching samples — but note that
  selecting a file also *overwrites* a hand-entered $\rho_0$ or $r_0$ from the file's header (Step 8).
- **`stog.inp` support is deliberately narrow.** Single dataset, zero $Q$ offset, zero second
  $y$ offset, no interactive rescale loop, filter enabled. Anything else is refused rather than
  guessed at, and blank lines inside the file will shift the positional layout.
- **No x-ray/electron form factors exist.** Any use on non-neutron data requires you to supply
  $\langle b\rangle^2$ and $\langle b^2\rangle$ yourself; the constant-$(a,b)$ model cannot represent $Q$-dependent $f(Q)$ or
  Compton contributions.
- **Nothing in Step 0 certifies the data.** The conditioning here removes padding, out-of-range
  points, and (optionally) glitches. It does not detect a bad normalization, a missing low-$Q$
  region, or a wrong composition — those show up, if at all, in the fit diagnostics documented
  in the following sections.


## Auto StoG — the transform layer: Fourier pair, Lorch, filter, Keen conversions

### What this page shows

The **Auto StoG** tab turns a measured total-scattering structure factor $S(Q)$ into the
RMCProfile-ready file family that the classic Fortran `stog` program produces — scaled $S(Q)$,
$g(r)$, Keen's $G(r)$, $D(r)$, $F(Q)$, and the Fourier-filter correction. Everything the page
draws and exports passes through one small, pure-numerical layer: a set of algebraic
Keen-convention conversions plus a discrete sine Fourier-transform pair, a Lorch window, an
analytic correction for the unmeasured $[0, Q_\mathrm{min}]$ range, and a Fourier filter.

That layer is [`rmc_toolkits/transforms.py`](../../rmc_toolkits/transforms.py) (source of truth,
pure NumPy, no file I/O) and its line-by-line JavaScript port at
[`web_app/frontend/src/workers/autoScale.js`](../../web_app/frontend/src/workers/autoScale.js)
(the top ~200 lines, before the level-sweep/fitting code). **The Auto StoG page runs the
JavaScript engine in *both* runtimes** — Flask mode and static GitHub-Pages mode — inside a Web
Worker ([`autoScaleWorker.js`](../../web_app/frontend/src/workers/autoScaleWorker.js)); it does not
call the backend (`AGENTS.md` architecture map; the Flask `/api/scaling/*` endpoints exist for
API/CLI consumers and run the Python engine). So the Python code is the *reference* and the
*CLI/API* path; the *page* runs the port. Parity is asserted by test, not by construction, and
the tests do not cover everything: **one real divergence is documented in Step 7** (the JS
auto-scale loop calls the Fourier filter without the $S(0)$ target that Python passes). See
[Python ↔ JavaScript parity](#python--javascript-parity) below.

Golden reference for every convention and equation number cited here: **D. A. Keen,
*J. Appl. Cryst.* **34**, 172–177 (2001)**, doi:10.1107/S0021889800019993. Equation numbers
below are the ones the code itself cites in its docstrings; the paper is not in the repository,
so the numbering is reported as the code states it, not independently verified. Two warnings
that the Caveats section expands on: the code's own citations are **not self-consistent** (the
`transforms.py` module docstring labels $F_K$ "Eq. 9" and $G_K$ "Eq. 10" while `sq_to_fk` /
`g_to_gk` label the *same two relations* "Eq. 19 inverted" and "Eq. 16 inverted"), and where
this section derives something the code does not label with an equation number, it says so
explicitly rather than inventing a citation.

The scale/offset fitting, level sweep, and $\rho_0$ self-consistency that *drive* this layer are
documented in the sibling Auto StoG engine section; this section covers only what the transforms
do to the numbers.

---

#### Notation and units

| Symbol | Meaning | Units |
| --- | --- | --- |
| $Q$ | momentum transfer, $4\pi\sin\theta/\lambda$ | Å⁻¹ |
| $r$ | interatomic separation | Å |
| $\rho_0$ | atomic number density | Å⁻³ (atoms per Å³) |
| $b_i$ | coherent neutron scattering length of species $i$ | fm (tables) → barn when squared |
| $c_i$ | atomic fraction of species $i$ | dimensionless |
| $\langle b\rangle^2 = \left(\sum_i c_i b_i\right)^2$ | "Faber-Ziman coefficient" of the classic `stog.inp` | barn |
| $\langle b^2\rangle = \sum_i c_i b_i^2$ | **different number** whenever there is more than one element | barn |
| $S(Q)$ | normalized total-scattering structure factor, $\to 1$ at high $Q$ | dimensionless |
| $F(Q) = Q\,[S(Q)-1]$ | the pystog / PDF-community $F(Q)$ | Å⁻¹ |
| $F_K(Q) = \langle b\rangle^2 [S(Q)-1]$ | **Keen's** $F(Q)$ | barn |
| $g(r)$ | pair distribution function; $0$ below the closest approach $r_0$, $\to 1$ | dimensionless |
| $G_\mathrm{PDF}(r) = 4\pi\rho_0 r\,[g(r)-1]$ | PDFFIT-style $G(r)$ | Å⁻² |
| $G_K(r) = \langle b\rangle^2 [g(r)-1]$ | **Keen's** $G(r)$; flat at $-\langle b\rangle^2$ below $r_0$ | barn |
| $D(r) = 4\pi\rho_0 r\,G_K(r)$ | Keen's $D(r)$ | barn·Å⁻² |
| $r_0$ | closest interatomic approach | Å |
| $a, b$ | the fitted affine correction $S_\mathrm{corr} = a\,S_\mathrm{meas} + b$ ($b$ dimensionless too) | dimensionless |
| $\Delta Q$, $\Delta r$ | grid spacings (`q` spacing as supplied; $r_\mathrm{max}/n_r$) | Å⁻¹, Å |
| $Q_0$ | **first supplied** $Q$ point, `q[0]` — the lower limit of the numerical integral | Å⁻¹ |
| $Q_\mathrm{max}$ | **last supplied** $Q$ point, `q[-1]` — truncation edge and Lorch width | Å⁻¹ |
| $M(Q)$ | the Lorch modification window (Step 4), code `lorch_window` / `lorchWindow` | dimensionless |
| $a_L = \pi/Q_\mathrm{max}$ | Lorch length: real-space resolution scale of the window; local `a` inside `low_q_correction_basis` | Å |
| $s_0 = F(Q_0)/Q_0 + 1$ | the measured $S(Q_0)$, code `s0` | dimensionless |
| $s_0^\mathrm{target}$ | assumed $S(0)$ for the low-$Q$ extrapolation, code `s0_target` / `effective_s0_target` | dimensionless |
| $f_1, f_2$ | the two closed-form moment integrals of the low-$Q$ correction, code `f1`, `f2` | Å⁻³, Å⁻² (non-Lorch) |
| $v = Q_0 r$; $v_\mp = Q_0(r \mp a_L)$; $v_a = 2a_LQ_0$ | dimensionless arguments of $f_1, f_2$, code `v`, `vm`/`vp`, `vpa` | dimensionless |
| $K(r)$ | the numerical "offset kernel", the transform of $F \equiv Q$, code `g_one` / `gOne` | Å⁻² |

**Symbol-collision warning.** The fitted scale $a$ and the Lorch length are *different
quantities* that the source code both calls `a` (the fitted `a` returned by `_solve_affine`;
the local `a = np.pi / q[-1]` inside `low_q_correction_basis`). This section writes the Lorch
length as $a_L$ throughout to keep them apart; when reading the source, disambiguate by scope.

The `b`-in-barn convention comes from the classic `stog.inp` line 18 (`0.015407` for Mn₃Sn) and
is produced by `faberZiman()` / `rmc_toolkits/scattering.py` as `bAvgSqBarn = bAvgSqFm2 / 100`
(1 barn = 100 fm²).

> **Do not conflate $\langle b\rangle^2$ and $\langle b^2\rangle$.** $\langle b\rangle^2$ sets the
> normalization and the low-$r$ level; $\langle b^2\rangle$ sets the $Q\to0$ limit. For Mn₃Sn
> ($b_\mathrm{Mn} = -3.73$ fm, $b_\mathrm{Sn} = 6.225$ fm) they differ by a factor
> $\langle b^2\rangle/\langle b\rangle^2 \approx 13$.

---

### Step 1 — Algebraic conversions (Keen 2001 definitions)

**Inputs** → arrays on a common grid plus the two scalars $\rho_0$ and $\langle b\rangle^2$.
**Outputs** → the same arrays in a different convention. These are pure, grid-preserving,
element-wise maps: no interpolation, no smoothing, no grid change.

$$F(Q) = Q\,[S(Q)-1] \qquad\Longleftrightarrow\qquad S(Q) = \frac{F(Q)}{Q} + 1$$

$$F_K(Q) = \langle b\rangle^2\,[S(Q)-1] \qquad\Longleftrightarrow\qquad S(Q) = \frac{F_K(Q)}{\langle b\rangle^2} + 1$$

$$G_\mathrm{PDF}(r) = 4\pi\rho_0 r\,[g(r)-1] \qquad\Longleftrightarrow\qquad g(r) = \frac{G_\mathrm{PDF}(r)}{4\pi\rho_0 r} + 1$$

$$G_K(r) = \langle b\rangle^2\,[g(r)-1] \qquad\Longleftrightarrow\qquad g(r) = \frac{G_K(r)}{\langle b\rangle^2} + 1$$

$$D(r) = 4\pi\rho_0 r\,G_K(r) \qquad\text{(Keen Eq. 29)}$$

$$D_\mathrm{low-}r(r) = -4\pi\rho_0 \langle b\rangle^2 r \qquad\text{(Keen Eqs. 29 + 15)}$$

The last line is the *theoretical* straight line that $D(r)$ must follow below $r_0$, because
$g(r)=0$ there makes $G_K = -\langle b\rangle^2$ (Keen Eq. 15). It is the C2 target of the
auto-scaler and the guide line drawn on the page.

**Two of these relations carry two different equation numbers in the source.** The
`transforms.py` module docstring labels $F_K$ "Eq. 9" and $G_K$ "Eq. 10"; the function
docstrings of `sq_to_fk` and `g_to_gk` label the identical formulas "Keen Eq. 19 inverted" and
"Keen Eq. 16 inverted". $G_\mathrm{PDF}$ is attributed to a "Keen Eq. 43/44 lineage" (module
docstring only). Nothing in the repo resolves the conflict — see the Caveats.

**Code** — [`rmc_toolkits/transforms.py`](../../rmc_toolkits/transforms.py):
`sq_to_fq()`, `fq_to_sq()`, `sq_to_fk()`, `fk_to_sq()`, `g_to_gpdf()`, `gpdf_to_g()`,
`g_to_gk()`, `gk_to_g()`, `gk_to_dr()`, `density_line()`.
The inversions divide by $Q$ and by $4\pi\rho_0 r$ respectively, so both require strictly
positive grids — enforced upstream: `crop_sq()` in [`scaling.py`](../../rmc_toolkits/scaling.py)
drops $Q \le 0$ unconditionally, and the $r$ grid starts at $\Delta r$, never $0$ (Step 2).

In JavaScript these conversions are **not** separate exported functions; they are inlined at
their call sites in [`autoScale.js`](../../web_app/frontend/src/workers/autoScale.js)
(`fourierFilter()`, `scalePipeline()`) as e.g. `fq[i] = q[i] * (sq[i] - 1)` and
`gk[i] = config.bAvgSq * (gFiltered[i] - 1)`. `density_line()` has no JS port; the equivalent
number is computed inline in `diagnosticsSummary()` as `d_r_low_r_slope_theory = -4 * Math.PI *
config.rho0 * config.bAvgSq`.

`tests/test_transforms.py::ConversionTests` asserts every round trip to `atol=1e-12`, the two
Keen forms ($F_K = \langle b\rangle^2[S-1]$, $D = 4\pi\rho_0 r G_K$) to `atol=1e-15`, and that
`density_line()` equals `gk_to_dr()` of a constant $-\langle b\rangle^2$ array.

---

### Step 2 — The grids

#### The $Q$ grid

The transforms **never resample $Q$**. They integrate on whatever grid the data arrive on,
after cropping (`crop_sq()` / `cropSq()`): finite rows only, $Q > 0$, and
$Q \in [q_\mathrm{min} - 10^{-12},\; q_\mathrm{max} + 10^{-12}]$. Non-uniform spacing is
supported by the quadrature rule (Step 3). Classic `stog` data files are typically already
rebinned to a uniform $\Delta Q$ (0.01 Å⁻¹ in the Mn₃Sn validation run) with `NaN` padding below
the true $Q_\mathrm{min}$; the `NaN` rows are dropped by the crop, so the working
$Q_\mathrm{min}$ is the **first finite measured point**, which can differ from the configured
`qmin`.

Two derived scalars matter downstream and both are taken from the *supplied* array, not from the
config:

- $Q_0 \equiv$ `q[0]` — the lower integration limit, used as $Q_\mathrm{min}$ by the omitted-low-$Q$ correction.
- $Q_\mathrm{max} \equiv$ `q[-1]` — the upper limit, used as the Lorch window width.

#### Two things the crop does besides cropping

**A hard minimum-point guard.** If fewer than **16** rows survive the finite/positive/window
test, the whole run aborts with `ValueError("fewer than 16 usable S(Q) points after cropping")`
(JS throws the same message). The guard is evaluated *before* despiking, so despiking can leave
fewer than 16 points without tripping it.

**Optional despiking — a destructive statistical filter, OFF by default.** With
`despike=True` the crop additionally deletes rows from the array that is then integrated:

1. a rolling median of width `despike_window = 7` (edge-padded, so the output length matches);
2. the residual $\varepsilon_i = S_i - \mathrm{median}_i$ and a robust scale
   $\mathrm{MAD} = 1.4826\,\mathrm{median}(|\varepsilon|)$;
3. keep only $|\varepsilon_i| \le$ `despike_nsigma` $\times \max(\mathrm{MAD}, 10^{-12})$, with
   `despike_nsigma = 6.0`.

Rejected points are **removed**, not interpolated or flagged, so the $Q$ grid handed to the
trapezoid rule acquires gaps that the quadrature then bridges with one wide panel. On crystalline
data this removes genuine Bragg maxima — the `ScalingConfig` docstring records **12% of points
flagged on the 59438 benchmark** — so it is intended for detector-glitch contamination only, and
the number actually dropped is reported as `n_despiked`. **Code** — `_despike_mask()` and
`crop_sq()` in [`scaling.py`](../../rmc_toolkits/scaling.py); `despikeKeepMask()` and `cropSq()` in
[`autoScale.js`](../../web_app/frontend/src/workers/autoScale.js). Same rule, same constants, same
order.

#### The $r$ grid

$$r_i = i\,\frac{r_\mathrm{max}}{n_r}, \qquad i = 1,2,\dots,n_r$$

so the grid **excludes $r=0$**, includes $r = r_\mathrm{max}$ exactly, and has uniform spacing
$\Delta r = r_\mathrm{max}/n_r$. Defaults $r_\mathrm{max} = 50$ Å, $n_r = 5000$ give
$\Delta r = 0.01$ Å — the classic `stog` grid.

**Code** — `ScalingConfig.r_grid` (property) in [`scaling.py`](../../rmc_toolkits/scaling.py):
`np.arange(1, nr + 1) * (rmax / nr)`; JS `rGrid(config)` in
[`autoScale.js`](../../web_app/frontend/src/workers/autoScale.js): `r[i] = (i + 1) * step`.
Identical construction.

**The $r$ grid is chosen independently of the data.** `r_grid` is a function of `rmax` and `nr`
alone — it never looks at `q` — and `fq_to_gpdf` accepts any `r` array with no consistency check
against the $Q$ range. Two consequences the plotted curves do not advertise, expanded in the
Caveats:

- With $Q_\mathrm{max} = 28$ Å⁻¹ the transform's intrinsic real-space resolution is
  $\pi/Q_\mathrm{max} \approx 0.11$ Å while $\Delta r = 0.01$ Å, so the output is **~11×
  oversampled** and neighbouring $r$ points are strongly correlated. Lorch (Step 4) makes the
  true resolution coarser still.
- The discrete transform's aliasing period is $2\pi/\Delta Q$ — 628 Å at $\Delta Q = 0.01$ Å⁻¹,
  comfortably beyond $r_\mathrm{max} = 50$ Å, but only 63 Å at $\Delta Q = 0.1$ Å⁻¹. Nothing in
  the code compares the two.

---

### Step 3 — The discrete sine transform pair

#### The continuum pair

**What the code asserts.** The `transforms.py` module docstring states the pair the code
implements, and **attaches no equation number to it**:

$$\boxed{\,G_\mathrm{PDF}(r) = \frac{2}{\pi}\int_0^\infty F(Q)\,\sin(Qr)\,dQ\,}\qquad
\boxed{\,F(Q) = \int_0^\infty G_\mathrm{PDF}(r)\,\sin(Qr)\,dr\,}$$

with $2/\pi$ forward and a bare $1$ backward. That, verbatim, is the whole of the code's claim
about the transform pair; the docstring's only justification is empirical ("validated against
pystog and a classic Fortran stog run").

**Reconstruction, flagged as such.** The standard Keen barns-scale forward/inverse pair is
written for $F_K(Q)$ and $G_K(r)$ with a $4\pi r^2 \sin(Qr)/(Qr)$ kernel and a
$1/[(2\pi)^3\rho_0]$ inverse prefactor. Substituting $G_K = \langle b\rangle^2 (g-1)$ and
$F_K = \langle b\rangle^2 (S-1)$ gives the intermediate forms
$\langle b\rangle^2 (S-1) = \frac{4\pi\rho_0}{Q}\int r\,\langle b\rangle^2 (g-1)\sin(Qr)\,dr$ and
$(g-1) = \frac{1}{2\pi^2\rho_0 r}\int Q\,(S-1)\sin(Qr)\,dQ$, which rearrange into the boxed pair
once both sides of the latter are multiplied by $4\pi\rho_0 r$. So the $1/[(2\pi)^3\rho_0]$
prefactor of the Keen form and the $2/\pi$ prefactor in the code are the *same normalization*
expressed in two conventions; the $\rho_0$ cancels because $G_\mathrm{PDF}$ carries a
$4\pi\rho_0 r$ that $G_K$ does not.

> **This derivation is this document's own reconstruction, not a code citation.** The repository
> cites **no equation number anywhere for the transform pair** — a grep for `Eq.` across
> `rmc_toolkits/` and `web_app/frontend/src/` returns only Eqs. 9, 10, 14, 15, 16, 19, 21, 29 and
> the "43/44 lineage". Keen 2001 is not vendored, so the kernel, the $1/[(2\pi)^3\rho_0]$
> prefactor and the intermediate $1/(2\pi^2\rho_0 r)$ form above are unverified here. **The code
> never implements the Keen-form prefactor directly** — it implements the boxed pair.

#### The quadrature

Both directions go through one primitive:

$$\texttt{sine\_transform}(x, y, w) \;=\; \int y(x)\,\sin(x\,w)\,dx \;\approx\; \sum_{i=1}^{N-1} \frac{x_i - x_{i-1}}{2}\,\Big[\,y_i \sin(x_i w) + y_{i-1}\sin(x_{i-1} w)\Big]$$

for every output point $w$. That is the **trapezoid rule on the supplied grid** — explicitly not
an FFT. The two sources in the repo justify that choice differently, and they are not the same
justification:

- `docs/STOG_SCALING_PLAN.md` §1.4 gives the design reason: "direct trapezoidal sine transform
  (`np.trapezoid`), not FFT — grids are short and the r-grid is arbitrary".
- The `transforms.py` **module docstring** instead justifies it by discretization compatibility:
  the trapezoid rule "matches the pystog `Transformer` discretization (and the Fortran stog to
  ~6e-4 rms)". `sine_transform`'s own docstring says nothing about either — only "Chunked over
  the output grid to bound the kernel-matrix memory."

Consequences worth being explicit about:

- The integral runs from $Q_0$ to $Q_\mathrm{max}$ (or $\Delta r$ to $r_\mathrm{max}$), **not**
  from 0 to $\infty$. The missing $[0, Q_0]$ panel is what Step 5 corrects analytically; the
  missing $(Q_\mathrm{max}, \infty)$ tail is a hard truncation and produces termination ripple
  unless the Lorch window (Step 4) is switched on.
- Cost is $O(N_Q \times N_r)$ sine evaluations. With defaults (~2700 $Q$ points — the 59438 run's
  1.0–28.0 Å⁻¹ at $\Delta Q = 0.01$ Å⁻¹ is exactly 2701 — and 5000 $r$ points) that is
  ~1.4×10⁷ sines *per transform*.
- One Fourier-filter call performs three sine transforms, but of **two different sizes**: the two
  forward transforms run on the full $r$ grid ($N_Q \times n_r$, ~1.4×10⁷ sine evaluations each
  at defaults), while the backward transform of the removed section runs only over
  $\lfloor r_\mathrm{cut}/\Delta r\rfloor \approx 100$ points at defaults — ~2.8×10⁵ sines, about
  50× cheaper. Budget a filter call as ~2, not 3, full transforms.
- **Neither engine guards a degenerate filter section**: $r_\mathrm{cut} < \Delta r$ yields an
  empty integration range, and $r_\mathrm{cut} \ge r_\mathrm{max}$ back-transforms the entire
  grid.
- Python bounds the kernel-matrix memory by chunking the **output** grid at
  `_SINE_CHUNK = 512` points (transforms.py line 40) and forming a `(512, N_Q)` matrix per
  chunk, reduced with `np.trapezoid` (`np.trapz` on NumPy < 2.0). JavaScript uses a plain
  sequential double loop. Same rule, different summation order — see parity below.

**Code** — `sine_transform()`, `fq_to_gpdf()`, `gpdf_to_fq()` in
[`transforms.py`](../../rmc_toolkits/transforms.py); `sineTransform()`, `fqToGpdf()` in
[`autoScale.js`](../../web_app/frontend/src/workers/autoScale.js). Note the JS port has **no
`gpdfToFq`** — the backward transform is called as bare `sineTransform(rSection, ySection, q)`
inside `fourierFilter()`, which is the identical operation (the Python `gpdf_to_fq()` is a
one-line alias for `sine_transform`).

**Tests** — `tests/test_transforms.py::FourierTransformTests`:
`test_round_trip_recovers_g` builds a synthetic $g(r)$ (tanh onset at 2.65 Å + Gaussian first
shell at 2.8 Å), forward-transforms it to $S(Q)$ on 0.01–35 Å⁻¹, transforms back, and requires
`atol=2e-2` against the model over $r\in[0.5, 8]$ Å.
`test_high_q_asymptote_and_low_r_zero` requires the synthetic $S(Q)$ tail mean above 25 Å⁻¹ to
sit within 5×10⁻³ of 1 (Keen Eq. 21) and the back-transformed $g$ to have rms < 3×10⁻² over
1.0–2.2 Å (below the closest approach).

---

### Step 4 — The Lorch modification window

**Input** → $F(Q)$ on the measured grid. **Math** →

$$M(Q) = \frac{\sin(\pi Q / Q_\mathrm{max})}{\pi Q / Q_\mathrm{max}}, \qquad M(0) \equiv 1,\quad M(Q_\mathrm{max}) = 0$$

**Output** → $F(Q)\,M(Q)$, which is then transformed exactly as in Step 3.

$Q_\mathrm{max}$ is `q[-1]`, the **last point of the supplied array** — not `config.qmax`. The
$Q=0$ singularity is removed explicitly (Python `np.divide(..., out=np.ones_like(x), where=x != 0)`;
JS `x === 0 ? 1 : Math.sin(x) / x`), which matters only for grids that actually contain $Q=0$.

Effect, stated honestly: multiplying $F(Q)$ by $M(Q)$ tapers the truncation edge to zero, which
suppresses the $\sin(Q_\mathrm{max} r)$ termination ripple; in exchange the real-space result is
convolved with the transform of the window, broadening every feature. The characteristic length
of that broadening is the **Lorch length** $a_L = \pi/Q_\mathrm{max}$ (0.112 Å at
$Q_\mathrm{max} = 28$ Å⁻¹) — the same quantity that appears as a removable singularity in the
Lorch branch of the low-$Q$ correction (Step 5), and the code's local variable `a` inside
`low_q_correction_basis` (not the fitted scale $a$). **The code computes no
resolution function and applies no deconvolution**; the page simply shows the windowed result.
The `ScalingConfig` docstring notes the practical trade: Lorch "damps termination ripples and
speeds loop convergence at the cost of real-space resolution".

The window is applied **only in the forward direction** (`fq_to_gpdf`). The backward transform
inside the Fourier filter is un-windowed.

**Code** — `lorch_window()` / `lorchWindow()`; applied in `fq_to_gpdf()` / `fqToGpdf()`.
**Default** — `lorch=False` in `ScalingConfig`, `defaultConfig`, and the page form
(`AutoStogPage.jsx`); a parsed `stog.inp` line 10 (`N`/`Y`) overrides it.
**Test** — `test_lorch_window_limits` pins $M \to 1$ as $Q\to0$,
$M(Q_\mathrm{max}/2) = \sin(\pi/2)/(\pi/2)$, and $M(Q_\mathrm{max}) = 0$.

---

### Step 5 — The analytic omitted-low-$Q$ correction

**Input** → the measured $F(Q)$, its first point, and the target $S(0)$.
**Output** → an additive term on the $r$ grid, added to $G_\mathrm{PDF}$ after the numerical
transform.

#### The model

Measured data start at $Q_0 > 0$, so the trapezoid integral silently drops
$\int_0^{Q_0} F(Q)\sin(Qr)\,dQ$. The correction models $S(Q)$ on the unmeasured interval as a
**straight line** from an assumed $S(0) = s_0^\mathrm{target}$ to the first measured value
$S(Q_0) \equiv s_0 = F(Q_0)/Q_0 + 1$:

$$S(Q) = s_0^\mathrm{target} + \left(s_0 - s_0^\mathrm{target}\right)\frac{Q}{Q_0}, \qquad 0 \le Q \le Q_0$$

so that $F(Q) = Q[S(Q)-1] = (s_0^\mathrm{target}-1)\,Q + \dfrac{s_0 - s_0^\mathrm{target}}{Q_0}Q^2$ and

$$\Delta G_\mathrm{PDF}(r) = \frac{2}{\pi}\left[\frac{s_0 - s_0^\mathrm{target}}{Q_0}\,f_1(r) \;+\; \left(s_0^\mathrm{target}-1\right) f_2(r)\right]$$

with the two closed-form moment integrals

$$f_1(r) = \int_0^{Q_0}\! Q^2 \sin(Qr)\,dQ = \frac{2v\sin v - (v^2-2)\cos v - 2}{r^3},\qquad
f_2(r) = \int_0^{Q_0}\! Q \sin(Qr)\,dQ = \frac{\sin v - v\cos v}{r^2},\qquad v \equiv Q_0 r$$

#### How the code stores it (why it is an affine *basis*)

`low_q_correction_basis()` returns the pair

$$\mathrm{coef}(r) = \frac{2}{\pi}\frac{f_1(r)}{Q_0}, \qquad \mathrm{const}(r) = \frac{2}{\pi}f_2(r)$$

and then, **only when $s_0^\mathrm{target} \neq 0$**, retargets the constant term in place:

$$\mathrm{const}'(r) = \left(1 - s_0^\mathrm{target}\right)\mathrm{const}(r) + s_0^\mathrm{target}\,\mathrm{coef}(r)$$

so that the correction is always evaluated as

$$\Delta G_\mathrm{PDF}(r) = \mathrm{coef}(r)\,s_0 - \mathrm{const}'(r)$$

Expanding $\mathrm{const}'$ reproduces the boxed expression above exactly. The point of the
split is that **the correction stays affine in $s_0$**, hence affine in the fitted $(a,b)$ — the
auto-scaler's closed-form least-squares solve depends on that, and `_solve_affine()` in
[`scaling.py`](../../rmc_toolkits/scaling.py) folds `coef`/`const` into its design matrix by hand
rather than letting `fq_to_gpdf` add the correction (a code comment states the reason:
"computing it inside each basis transform would multiply-count its constant term").

With $s_0^\mathrm{target}=0$ the formula reduces to `pystog`'s `Transformer._low_x_correction`
(original author Jack Carpenter) — the "solid-state limit" special case, per the
`low_q_correction_basis` docstring.

#### The $S(0)$ target

$$
s_0^\mathrm{target} = \begin{cases}
\texttt{s0\_target} & \text{if the user pinned one}\\[2pt]
1 - \dfrac{\langle b^2\rangle}{\langle b\rangle^2} & \text{if a composition supplied } \langle b^2\rangle \quad\text{(Keen Eq. 21)}\\[6pt]
0 & \text{otherwise (pystog solid-state limit)}
\end{cases}
$$

**Code** — `ScalingConfig.effective_s0_target` (property) / `effectiveS0Target(config)` in JS.
For a strongly negative-$b$ composition such as Mn₃Sn the Keen limit is
$1 - 13.06 = -12.06$, i.e. $O(-10)$, and the docstring records that using it instead of 0
"removes an O(1) bias in the low-r transform".

#### The Lorch branch

When `lorch=True` the correction must include the same window over the extrapolated segment, so
the integrand carries $M(Q) = \sin(a_LQ)/(a_LQ)$ with $a_L = \pi/Q_\mathrm{max}$. Using
$\sin(a_LQ)\sin(rQ) = \tfrac12[\cos((r-a_L)Q) - \cos((r+a_L)Q)]$ the two moments become, with
$v_\mp = Q_0(r \mp a_L)$:

$$f_1(r) = \frac{1}{2a_L}\left[\frac{v_-\sin v_- + \cos v_- - 1}{(r-a_L)^2} - \frac{v_+\sin v_+ + \cos v_+ - 1}{(r+a_L)^2}\right],\qquad
f_2(r) = \frac{1}{2a_L}\left[\frac{\sin v_-}{r-a_L} - \frac{\sin v_+}{r+a_L}\right]$$

These have a **removable singularity at $r = a_L = \pi/Q_\mathrm{max}$**, patched with the
analytic limits ($\lim_{k\to0}(kc\sin kc + \cos kc - 1)/k^2 = c^2/2$ and
$\lim_{k\to 0}\sin(kc)/k = c$, with $c = Q_0$):

$$f_1(a_L) = \frac{1}{2a_L}\left[\frac{Q_0^2}{2} - \frac{v_a \sin v_a + \cos v_a - 1}{(2a_L)^2}\right],\qquad
f_2(a_L) = \frac{1}{2a_L}\left[Q_0 - \frac{\sin v_a}{2a_L}\right], \qquad v_a = 2a_LQ_0$$

applied wherever $|r - a_L| \le 10^{-9}\max(1, a_L)$ (Python `np.isclose(..., rtol=0.0,
atol=1e-9*max(1.0, a))`; JS the same explicit bound).

#### Edge cases and defaults

- If `q[0] == 0` the correction is identically zero (`coef = const = 0`, and
  `omitted_low_q_correction` returns zeros): data that start at $Q=0$ omit nothing, and the
  $[0, q_1]$ panel is already inside the trapezoid sum. This is stated as pystog parity
  ($F_1 = F_2 = 0$).
- $f_1, f_2$ are forced to 0 at $r = 0$ **in the non-Lorch branch only** (both are $O(r)$ there
  anyway): Python `np.where(r == 0, 0.0, ...)` after an `errstate` guard, JS `continue` on
  `ri === 0`. **The Lorch branch has no $r = 0$ guard in either engine** — its only patch is the
  removable singularity at $r = a_L$. This is harmless on the shipped grid, which starts at
  $\Delta r$, but `low_q_correction_basis` / `lowQCorrectionBasis` are public functions that
  accept an arbitrary $r$ array and are unprotected there.
- **Two different default policies, one layer apart.** `transforms.py`'s own signatures default
  `lorch=False`, `low_q_correction=False`, `s0_target=0.0` — every entry point
  (`fq_to_gpdf`, `low_q_correction_basis`, `fourier_filter`) is **uncorrected and un-windowed by
  default**. The ON default lives one layer up, in `ScalingConfig.low_q_correction = True` /
  `defaultConfig.lowQCorrection = true` (page checkbox on), which `_pipeline` / `scalePipeline`
  pass through. This matters concretely: a direct caller of the documented reference API — this
  repo's own `StogRunParityTests` included — runs the **uncorrected** transform and gets a
  different number from the page on identical data. The `ScalingConfig` docstring records the
  measured motivation for the ON policy: on the synthetic benchmark with
  $Q_\mathrm{min} = 0.6$ Å⁻¹ the correction cuts the recovered-scale bias from 8% to 0.3%. Turn
  it off only for strict classic-`stog` parity.

**Code** — `low_q_correction_basis()`, `omitted_low_q_correction()` in
[`transforms.py`](../../rmc_toolkits/transforms.py); `lowQCorrectionBasis()` and the inline
correction in `fqToGpdf()` in [`autoScale.js`](../../web_app/frontend/src/workers/autoScale.js).

**Tests** —
`test_low_q_correction_basis_consistency` pins
$\texttt{omitted\_low\_q\_correction} = \mathrm{coef}\cdot s_0 - \mathrm{const}$ to `atol=1e-14`;
`test_low_q_correction_zero_when_data_start_at_zero` requires *exact* zeros (`atol=0.0`) when
`q[0] == 0`;
`test_low_q_correction_lorch_finite_at_singular_r` requires finiteness at $r = \pi/Q_\mathrm{max}$
and that the patched value is bracketed by its two neighbours;
`test_low_q_correction_improves_truncated_transform` truncates the synthetic data at 0.9 Å⁻¹ and
requires the corrected low-$r$ rms to beat the uncorrected one;
`tests/test_scaling.py::DetectionAndDensityTests::test_s0_target_changes_only_the_constant_term`
pins the retarget algebra ($\mathrm{coef}$ unchanged to `atol=1e-15`, $\mathrm{const}$ to the
formula above at `atol=1e-12`).

---

### Step 6 — The assembled forward transform

Putting Steps 3–5 together,
`fq_to_gpdf(q, fq, r, *, lorch=False, low_q_correction=False, s0_target=0.0)` — note again that
those are the *function's* defaults, all inert; the engine supplies the live values — computes

$$G_\mathrm{PDF}(r_k) \;=\; \frac{2}{\pi}\!\!\sum_{i=1}^{N_Q-1}\!\frac{Q_i - Q_{i-1}}{2}\Big[w_i F_i \sin(Q_i r_k) + w_{i-1}F_{i-1}\sin(Q_{i-1} r_k)\Big] \;+\; \underbrace{\Big[\mathrm{coef}(r_k)\,s_0 - \mathrm{const}'(r_k)\Big]}_{\text{only if } \texttt{low\_q\_correction}}$$

with $w_i = M(Q_i)$ when `lorch=True` and $w_i \equiv 1$ otherwise, and
$s_0 = F_0/Q_0 + 1$ computed from the **un-windowed** $F(Q)$ (the window enters the correction
analytically instead). $g(r)$ then follows by `gpdf_to_g()`.

Because the transform is linear in the data, $G_\mathrm{PDF}[a S + b]$ is affine in $(a, b)$;
that is what lets the auto-scaler precompute three basis transforms
(`fq_to_gpdf(q, sq_to_fq(q, sq), rw)`, `fq_to_gpdf(q, q, rw)` — the transform of $F \equiv Q$,
i.e. the numerical "offset kernel" $K(r)$ — and the filter term) once per solve and finish in
closed form. Crucially, **the offset kernel is computed numerically on the same $Q$ window, grid,
and window function as the data**, not analytically, so truncation effects cancel
self-consistently (`docs/STOG_SCALING_PLAN.md` §1.4). Details in the engine section.

---

### Step 7 — The Fourier filter

**Input** → $Q$, $S(Q)$ (already scaled), the $r$ grid, $\rho_0$, a cutoff $r_\mathrm{cut}$, and
the Lorch / low-$Q$ flags. **Output** → `(sq_filtered, sq_ft, g_filtered)`.

The physical statement: no atom pairs exist below the closest approach, so any structure in
$g(r)$ below $r_\mathrm{cut}$ is an artifact (background, normalization, truncation). The filter
back-transforms exactly that artifact into $Q$ space and subtracts it.

**Procedure** (verbatim order of `fourier_filter()`):

1. Reject the call if `q[0] <= 0` (`ValueError` / JS `Error`) — the $S \leftrightarrow F$
   conversions divide by $Q$.
2. $F(Q) = Q[S(Q)-1]$.
3. $G_\mathrm{PDF}(r)$ by the Step-6 forward transform (Lorch, low-$Q$ correction and
   $s_0^\mathrm{target}$ as passed in — see the divergence note below for the one call site where
   the JS engine passes a different $s_0^\mathrm{target}$ than Python), then
   $g(r) = G_\mathrm{PDF}/(4\pi\rho_0 r) + 1$.
4. Select the section $\{r : r \le r_\mathrm{cut}\}$ and back-transform the **removed** content,

   $$F_\mathrm{ft}(Q) = \int_{\Delta r}^{r_\mathrm{cut}} 4\pi\rho_0\, r\, g(r)\,\sin(Qr)\,dr$$

   which is precisely the change in $G_\mathrm{PDF}$ caused by setting $g \to 0$ on the section:
   $4\pi\rho_0 r(g-1) \to -4\pi\rho_0 r$, a difference of $4\pi\rho_0 r g$. (The code comment
   describes the same thing as "pystog shifts the section by +1 and re-derives $G_\mathrm{PDF}$".)
   This back-transform uses **no Lorch window and no low-$Q$ correction**, and its lower limit is
   $\Delta r$, so the $[0, \Delta r]$ panel is dropped — negligible because the integrand
   $\propto r\,g(r)$ vanishes at the origin, but it is a genuine (tiny) discretization choice.
5. $S_\mathrm{ft}(Q) = F_\mathrm{ft}(Q)/Q + 1$ — **this array is the classic `ft.dat` file**
   (written by `scaling_cli.py::_write_outputs` as `targets["ft_correction"]`).
6. $F_\mathrm{filt} = F - F_\mathrm{ft}$, $\;S_\mathrm{filt} = F_\mathrm{filt}/Q + 1$.
7. Re-run the Step-6 forward transform on $F_\mathrm{filt}$ (Lorch + low-$Q$ correction again,
   the latter now using $S_\mathrm{filt}(Q_0)$) and convert to $g_\mathrm{filt}(r)$.

#### The `ft.dat` oracle relation

Steps 5–6 imply, exactly and by construction,

$$S_\mathrm{filtered}(Q) = S(Q) - \big[S_\mathrm{ft}(Q) - 1\big]$$

which is the relation `scale_ft.sq = scale.fq - (ft.dat - 1)` decoded from the real Fortran run
(`docs/STOG_SCALING_PLAN.md` §1.3). Because `ft.dat` is present in the validation fixture, it is
a **direct oracle for the filter's discretization, independent of any documentation**.
`tests/test_transforms.py::test_fourier_filter_identity_and_cleanup` asserts the identity to
`atol=1e-10` on synthetic data, and also that the filter *reduces* the sub-cutoff rms of $g$
relative to the unfiltered transform.

#### Honest limits of the filter

The filter is **not an exact projection**. The subtracted term is band-limited to
$[Q_0, Q_\mathrm{max}]$ and re-transformed with the same truncation, so sub-cutoff content is
reduced, not annihilated; the test asserts only "less than before". It also *injects an
assumption* — that $g \equiv 0$ below $r_\mathrm{cut}$ — into the $S(Q)$ that is subsequently
exported, so the filtered $S(Q)$ is no longer a raw measurement. Finally, the classic workflow's
filter performs part of the normalization: the plan records that the hand-scaled Mn₃Sn example's
*unfiltered* tail sits at ~0.94 and the filter's correction lifts it to 1.000, which is why the
auto-scaler always evaluates its high-$Q$ criterion on filtered data.

**Code** — `fourier_filter()` in [`transforms.py`](../../rmc_toolkits/transforms.py);
`fourierFilter()` in [`autoScale.js`](../../web_app/frontend/src/workers/autoScale.js). One
cosmetic implementation difference: Python selects the section with the boolean mask
`r <= cutoff`, JS with a monotone prefix scan (`while (r[sectionEnd] <= cutoff) sectionEnd += 1`).
These agree for any ascending $r$ grid, which is all `r_grid`/`rGrid` ever produces.

> #### Known divergence: the JS auto loop drops `s0Target`
>
> The two engines do **not** call the filter identically inside the auto-scale iteration.
> Python's `_pipeline()` always passes `s0_target=config.effective_s0_target`, and it is the only
> Python caller of `fourier_filter`. On the JavaScript side `scalePipeline()` does pass
> `s0Target: effectiveS0Target(config)` — but the iteration inside `autoscale()` calls
> `fourierFilter(q, sqScaled, r, { rho0, cutoff, lorch, lowQCorrection })` **with no `s0Target`
> key**, so the destructuring default `s0Target = 0` applies.
>
> The consequence is not cosmetic. That call's `sqFt` becomes the fed-back
> $\Delta S = S_\mathrm{ft} - 1$ that seeds the next `solveAffine`, so **whenever
> $s_0^\mathrm{target} \neq 0$ — i.e. whenever a composition supplies $\langle b^2\rangle$, which
> is the composition-first workflow the page advertises — the fitted $(a,b)$ differ between the
> two engines.** For Mn₃Sn the dropped target is $\approx -12$. The *final* `scalePipeline` filter
> does pass the target in both engines, so the exported arrays agree **given the same $(a,b)$**;
> it is the $(a,b)$ themselves that diverge.
>
> The golden-fixture parity tests do not catch it: the `auto` / `manual` / detection fixture
> configs set no `b_sq_avg` (so $s_0^\mathrm{target} = 0$ in both engines), the one `fz` case that
> does set it takes `autoscale`'s early-return branch and never enters the loop. The
> `estimate_rho0` case *does* exercise the loop with a nonzero target — the fixture's
> $\langle b\rangle^2 = 0.02$, $\langle b^2\rangle = 0.0347972$ barn give
> $s_0^\mathrm{target} = -0.74$ — but only under the loose `1e-4` relative tolerance on $\rho_0$,
> which the measured 2.8·10⁻⁵ shift sits inside (engine section
> [parity](#python--javascript-parity-1), difference 4). This is an
> unreported divergence in the shipped code, not a rounding-level difference — it is recorded
> here because a reader has no other way to know the page's fitted scale comes from a
> differently-corrected filter term than the CLI's.

**Default** — `r_cutoff = 1.0` Å in both engines (`ScalingConfig.r_cutoff`,
`defaultConfig.rCutoff`); a parsed `stog.inp` line 15 overrides it.

---

### Step 8 — Conversion of the filtered result to the RMCProfile functions

From $g_\mathrm{filt}$ and $S_\mathrm{filt}$, the exported functions are

$$G_K(r) = \langle b\rangle^2\left[g_\mathrm{filt}(r) - 1\right],\qquad
D(r) = 4\pi\rho_0\, r\, G_K(r),\qquad
F_K(Q) = \langle b\rangle^2\left[S_\mathrm{filt}(Q) - 1\right]$$

**Code** — `scale_pipeline()` in [`scaling.py`](../../rmc_toolkits/scaling.py) (`g_to_gk`,
`gk_to_dr`, `sq_to_fk`); `scalePipeline()` in
[`autoScale.js`](../../web_app/frontend/src/workers/autoScale.js) (inlined).

The classic file family maps onto these as follows (`docs/STOG_SCALING_PLAN.md` §1.3, pinned
empirically against the Fortran run, and `scaling_cli.py::_write_outputs`):

| Classic file | Contents written by this repo |
| --- | --- |
| `scale.fq` | the scaled $S(Q)$ (S-convention despite the name) |
| `scale.gr` | $g(r) - 1$ from the **unfiltered** transform, recomputed with the same discretization |
| `ft.dat` | $S_\mathrm{ft}(Q)$, the filter correction (Step 7.5) |
| `scale_ft.sq` | $S_\mathrm{filt}(Q)$ |
| `scale_ft.gr` | column 1 = $r$, column 2 = $g_\mathrm{filt}-1$, column 3 = $4\pi\rho_0 r (g_\mathrm{filt}-1) = G_\mathrm{PDF}$ |
| `scale_ft_rmc.fq` | $F_K(Q)$ (barn) |
| `scale_ft_rmc.gr` | $G_K(r)$ (barn), **after enforcement if enabled** |
| `scale_ft_rmc.dr` | $D(r)$ (barn·Å⁻²), **after enforcement if enabled** |

Every one of these is written by `write_stog_xy(path, x, y, extra=…)`, which emits $x$ first, so
**column 1 is always $r$ (or $Q$)**; the table numbers physical columns. Note that
`docs/STOG_SCALING_PLAN.md` §1.3 numbers only the *value* columns and labels `scale_ft.gr`'s
second value column "D(r)" — but the code writes
$4\pi\rho_0 r (g-1) = G_\mathrm{PDF} = D(r)/\langle b\rangle^2$, which for the 59438 run
differs from $D(r)$ by a factor $\langle b\rangle^2 = 0.015407$ barn. The code is what ships;
a reader cross-checking against
§1.3 will hit this contradiction, so it is flagged here rather than papered over.

`scale_ft.gr` always holds **pre-enforcement** values; only the three `_rmc` files carry the
enforced arrays (and only `.gr`/`.dr` — `_rmc.fq` is a $Q$-space function and is untouched).

**The browser export path re-derives these columns on the main thread, not in the worker.**
`AutoStogPage.jsx::writeFiles` recovers $g-1$ from the already-converted Keen function as
`gm1 = series.gk.map(v => v / config.bAvgSq)` — an exact float round-trip of the Step-8
conversion — and builds the `_ft.gr` third column as `4πρ₀ r · gm1` there. The CLI does not round
trip: `scaling_cli.py::_write_outputs` uses `result.g_filtered - 1.0` directly. Two further
browser-only details:

- The unfiltered `scale.gr` uses the worker's separately recomputed `gm1Unfiltered`
  (`autoScaleWorker.js` re-runs `fqToGpdf` on the scaled $F(Q)$ with the full config, exactly as
  the CLI does, because it is not part of the engine result) — **falling back to the *filtered*
  `gm1` if that array is missing** (`series.gm1Unfiltered || gm1`), which would silently
  substitute filtered for unfiltered content.
- `_rmc.gr` / `_rmc.dr` likewise fall back to the un-enforced `series.gk` / `series.dr` when
  `gkEnforced` / `drEnforced` are absent (see Step 9 for when that happens).

---

### Step 9 — Low-$r$ enforcement (classic `stog` parity)

The final Fortran `stog` stage **hard-replaces** the low-$r$ region of the RMC outputs with the
exact theoretical values. Two functions implement it:

`first_peak_zero(r, g, cutoff, peak_rmin, peak_rmax)` — the exact Fortran semantics
(`stog_new3.f90` L398 per the plan's risk table): set $g(r) = 0$ wherever

$$r \le \texttt{cutoff} \quad\textbf{and}\quad \big(r \ge \texttt{peak\_rmax}\ \textbf{or}\ r \le \texttt{peak\_rmin}\big)$$

i.e. remove sub-peak and inter-peak ripples while keeping a first peak that starts below the
cutoff. A zeroed $g$ becomes $-\langle b\rangle^2$ in the written $G_K$.

`enforce_low_r(r, gk, cutoff, b_avg_sq)` — the degenerate case, applied directly in $G_K$ space:
$G_K(r \le \texttt{cutoff}) = -\langle b\rangle^2$. This is what `first_peak_zero` reduces to
when the peak window lies entirely above the cutoff, as in the validation run
(`2.48 2.65 3.1`: `rmcmin = 2.65 > rmccut = 2.48`).

#### Which function actually runs depends on the entry point

There are **three enforcement paths with two different semantics**, and the choice is not a
configuration knob — it is a property of how you called the code:

| Entry point | Function used | Semantics |
| --- | --- | --- |
| Python engine (`scale_pipeline()`, hence `ScalingResult.gk_enforced` / `d_r_enforced`) | `enforce_low_r`, when `config.enforce_cutoff` is not `None` | **always** the flat $-\langle b\rangle^2$ replacement; `first_peak_zero` is never called |
| CLI writer (`scaling_cli.py::_write_outputs`) | `first_peak_zero` with the full `(cutoff, peak_rmin, peak_rmax)` triple | true Fortran ripple removal when the window is genuine |
| Browser worker (`autoScaleWorker.js`) | `firstPeakZero` with whatever window the page resolved | same as the CLI |

So a Python API caller with a *genuine* first-peak window (`peak_rmin < cutoff`) silently gets
the flat replacement instead of the Fortran ripple removal, and
`ScalingResult.gk_enforced` is **not** the same array as the CLI's `scale_ft_rmc.gr` for the same
run. Only the CLI and the browser reach the general semantics.

#### How the enforcement boundary is chosen in `'auto'` mode

`'auto'` resolves the cutoff to a **data-derived** closest approach from
`detect_first_peak_onset()` ([`scaling.py`](../../rmc_toolkits/scaling.py); ported verbatim as
`detectFirstPeakOnset()` in [`autoScale.js`](../../web_app/frontend/src/workers/autoScale.js)).
Because Step 9's headline caveat is that the Keen limits become true *by construction* below that
boundary, the boundary's provenance is part of the transparency story. The procedure, with every
constant:

1. Restrict to the search window $[\,$`search_min`$,\,$`search_max`$\,]$. Defaults are
   $[1.0, 6.0]$ Å, but **every caller passes `search_min = r_cutoff + 0.3`** (= 1.3 Å at the
   default $r_\mathrm{cut}$) — `autoscale()`, `scaling_cli.py`, and `autoScaleWorker.js` alike.
   Fewer than 3 points in the window → `None`.
2. Take the maximum of $|g|$ in that window. **Absolute value**, because Faber-Ziman totals of
   negative-$b$ compositions (Mn₃Sn) can have an *inverted* first shell.
3. Reject if that peak $<$ `floor` $= 0.5$ → `None`.
4. Threshold `level` $= \max(0.5,\ 0.35 \times \mathrm{peak})$ (`fraction = 0.35`). Peak-relative
   because sub-$r_0$ truncation ripples scale with the fitted amplitude — the docstring notes
   they can reach $O(\mathrm{peak}/3)$ on missing-low-$Q$ data, so no fixed threshold separates
   them.
5. Walk left from the peak index while $|g| >$ `level`; return `r[index + 1]`, the first grid
   point at or above the crossing. If $|g|$ never drops below `level` inside the window →
   `None`.

The `qmax` parameter is accepted but unused (`# noqa: ARG001`, kept for signature stability); the
JS port does not take it at all. **All five constants are heuristics, and no test pins their
sensitivity** — there is no test that perturbs `fraction`, `floor`, or the search bounds and
checks how $r_0$ moves.

#### Parity in the browser, and two undocumented branches

`autoScale.js` ports only `firstPeakZero`. `autoScaleWorker.js` obtains the flat-replacement
behaviour by passing a degenerate window, `{cutoff: r0, peakRmin: r0, peakRmax: r0}`, when
enforcement is set to `'auto'` — the boolean then collapses to $r \le r_0$, identical to
`enforce_low_r`. The same convention is used by the CLI (`scaling_cli.py`,
`enforcement = (float(r0_detected),) * 3`). Two control-flow branches on the page deserve
stating:

- **`'auto'` enforcement is not guaranteed to happen.** The worker first tries to recover `r0`
  when the engine did not detect one (manual runs skip `autoscale`'s two-pass detection, so the
  worker re-runs `detectFirstPeakOnset` with `searchMin: config.rCutoff + 0.3`, exactly like the
  CLI's post-run detection). But if detection still returns `null`, `effectiveEnforcement` is set
  to `null` and **no enforcement is applied at all**: `gkEnforced` / `drEnforced` come back
  `null`, the exported `_rmc.gr` / `_rmc.dr` fall back to the un-enforced arrays, and a checked
  "Enforce low-r" box has silently become a no-op. The only signal is that the page reports
  `r0_detected` as absent. The worker's own comment flags this as the failure mode the recovery
  guards against — the guard cannot manufacture an `r0` when detection genuinely fails.
- **The window shape depends on where the cutoff came from.** In
  `AutoStogPage.jsx::resolveEnforcement`, a cutoff **typed into the form** always produces a
  degenerate window (`peakRmin = peakRmax = cutoff`), i.e. the flat $-\langle b\rangle^2$
  replacement. Only an **empty cutoff field with a parsed `stog.inp` loaded** uses that file's
  `(peakCutoff, peakRmin, peakRmax)` — which is the one path on the page that reaches the true
  Fortran inter-peak ripple removal.

> **This is the single most important transparency point on the page.** After enforcement, Keen
> Eq. 15 ($G_K = -\langle b\rangle^2$ below $r_0$) and Eq. 29+15 (the $D(r)$ straight line) hold
> **by construction, not as evidence of a good fit**. `enforce_low_r`'s docstring says so
> explicitly — "The published files therefore satisfy Keen Eq. 15 by construction — judge fit
> quality on the *pre*-enforcement residual"; `first_peak_zero`'s docstring states only the
> Fortran zeroing rule and that a zeroed $g$ "becomes `-<b>^2` in the written Keen G(r)", with no
> by-construction warning. Both engines nevertheless always report the *pre-enforcement* low-$r$
> rms (`low_r_rms_pre_enforcement` in `diagnostics_summary()` / `diagnosticsSummary()`) as the
> honest fit-quality metric.

**Display-only choices on the page** (`AutoStogPage.jsx`; none of them alter the exported
numbers): the enforced curve is plotted as primary and the pre-enforcement curve as a
hidden-by-default series; the $-\langle b\rangle^2$ and $-4\pi\rho_0\langle b\rangle^2 r$ theory
guide lines are drawn only out to `RMAX_DISPLAY = 8` Å even though the curves themselves run to
$r_\mathrm{max} = 50$ Å; and the $G_K(r)$ panel opens on a hard-coded y-window
$[2.1\,(-\langle b\rangle^2),\ -3.2\,(-\langle b\rangle^2)]$ scaled off the theory level, so the
initially visible amplitude range is a design choice, not the data range.

**Tests** — `test_enforce_low_r` (flat below, untouched above, `atol=1e-15`);
`test_first_peak_zero_fortran_semantics` (both the general window and the degenerate 2.48/2.65/3.1
layout); `autoScale.test.js` asserts `firstPeakZero` produces exact `0` for every
$r \le 2.48$.

---

### Parameters and defaults

**Read the two default columns separately.** The three flags that control the transform have
*different* defaults in `transforms.py` (the pure reference API) and in
`ScalingConfig`/`defaultConfig` (the engine policy the page and CLI apply). Calling the reference
API directly gives you the uncorrected, un-windowed transform.

| Parameter | `transforms.py` default | `ScalingConfig` / `defaultConfig` default | Units | Where else set | Effect on this layer |
| --- | --- | --- | --- | --- | --- |
| `lorch` | `False` | `False` | — | `stog.inp` line 10 | applies $M(Q)$ to $F(Q)$ before the forward transform |
| `low_q_correction` | **`False`** | **`True`** | — | page checkbox | adds the analytic $[0, Q_0]$ term |
| `s0_target` | **`0.0`** | `None` → `effective_s0_target` | — | library only — never set by the page, CLI or API | resolves to $1-\langle b^2\rangle/\langle b\rangle^2$ if $\langle b^2\rangle$ known, else 0 |

Everything else:

| Parameter | Default | Units | Where set | Effect on this layer |
| --- | --- | --- | --- | --- |
| `qmin`, `qmax` | required | Å⁻¹ | `ScalingConfig` / `defaultConfig` / `stog.inp` line 3 | crop window; **not** the $Q_0$/$Q_\mathrm{max}$ the transforms use |
| `rmax` | 50.0 | Å | `ScalingConfig` / `defaultConfig` / `stog.inp` line 8 | upper end of the $r$ grid |
| `nr` | 5000 | — | `ScalingConfig` / `defaultConfig` / **`stog.inp` line 9** | $r$ grid size; $\Delta r = r_\mathrm{max}/n_r = 0.01$ Å |
| `r_cutoff` | 1.0 | Å | `ScalingConfig` / `defaultConfig` / `stog.inp` line 15 | Fourier-filter section boundary |
| `rho0` | required | Å⁻³ | user / `NUMBER_DENSITY ::` header / `stog.inp` line 11 / self-consistent estimate | $G_\mathrm{PDF}\leftrightarrow g$ and $D(r)$ |
| `b_avg_sq` | required | barn | user / formula / `stog.inp` line 18 | $G_K$, $F_K$, the $-\langle b\rangle^2$ level |
| `b_sq_avg` | `None` | barn | formula (`faber_ziman`) | sets the default $S(0)$ target and the Keen Eq. 14 diagnostic |
| `despike` | `False` | — | `ScalingConfig` / `defaultConfig` / page checkbox | **deletes** rolling-median outliers from $S(Q)$ before any transform |
| `despike_window` | 7 | points | same | rolling-median width (edge-padded) |
| `despike_nsigma` | 6.0 | MAD units | same | keep threshold $\lvert\varepsilon\rvert \le n\sigma\max(\mathrm{MAD}, 10^{-12})$ |
| minimum surviving points | 16 | points | `crop_sq` / `cropSq` (hard-coded) | fewer → `ValueError`; checked **before** despiking |
| `enforce_cutoff` | `None` (Python engine) | Å | `ScalingConfig`; CLI/page default to `stog.inp` `peak_cutoff` or the detected $r_0$ | low-$r$ hard replacement (Python engine: always `enforce_low_r`) |
| `search_min` (peak onset) | 1.0 default, **1.3** in practice ($r_\mathrm{cut} + 0.3$) | Å | `detect_first_peak_onset`; every caller overrides | lower bound of the $r_0$ search |
| `search_max` (peak onset) | 6.0 | Å | same | upper bound of the $r_0$ search |
| `fraction` (peak onset) | 0.35 | — | same | flank threshold $=\max(\mathrm{floor}, \mathrm{fraction}\times\mathrm{peak})$ |
| `floor` (peak onset) | 0.5 | — | same | minimum $\lvert g\rvert$ for a feature to count as the first shell |
| `_SINE_CHUNK` | 512 | output points | `transforms.py` line 40 | memory bound only; no numerical effect |
| Lorch singularity tolerance | $10^{-9}\max(1, a_L)$ | Å | `low_q_correction_basis` | switches to the analytic limit at $r = a_L = \pi/Q_\mathrm{max}$ |
| $Q_\mathrm{max}$ for Lorch | `q[-1]` | Å⁻¹ | implicit | last *supplied* point, not `config.qmax` |
| $Q_0$ for the low-$Q$ correction | `q[0]` | Å⁻¹ | implicit | first *finite measured* point, not `config.qmin` |
| `RMAX_DISPLAY` | 8 | Å | `AutoStogPage.jsx` | length of the plotted theory guide lines only; no effect on exports |

All `stog.inp` line numbers above are **1-based indices into the file's non-empty lines** —
`readStogInp` / `StogInput` strip blanks before indexing (`lines[7]` is `rmax`, `lines[8]` is
`nr`, and so on).

---

### Python ↔ JavaScript parity

| Item | Python | JavaScript | Agreement |
| --- | --- | --- | --- |
| Quadrature rule | `np.trapezoid` on a chunked kernel matrix | sequential trapezoid loop | Same rule; **summation order differs** (NumPy pairwise vs JS sequential), so results differ at float-rounding level only |
| Lorch window | `lorch_window` | `lorchWindow` | Identical formula and $Q=0$ guard |
| Low-$Q$ basis (both branches, singular patch) | `low_q_correction_basis` | `lowQCorrectionBasis` | Identical formulas, identical tolerance |
| Forward transform | `fq_to_gpdf` | `fqToGpdf` | Identical |
| Backward transform | `gpdf_to_fq` (alias of `sine_transform`) | no alias; `sineTransform` called directly | Identical operation |
| Fourier filter (the function) | `fourier_filter` | `fourierFilter` | Identical formulas; section selected by mask vs prefix scan (equivalent for ascending $r$) |
| Filter **call inside the auto loop** | `_pipeline` passes `s0_target=config.effective_s0_target` | `autoscale` omits `s0Target` → default `0` | **Differs** whenever $s_0^\mathrm{target}\neq0$ (i.e. whenever a composition supplies $\langle b^2\rangle$); the fed-back $\Delta S$ and hence the fitted $(a,b)$ diverge. See Step 7 |
| Filter call in the final pipeline | `_pipeline` (same) | `scalePipeline` **does** pass `s0Target` | Identical |
| Crop + despike | `crop_sq` / `_despike_mask` | `cropSq` / `despikeKeepMask` | Identical rule, constants and ordering |
| Peak-onset detection | `detect_first_peak_onset` | `detectFirstPeakOnset` | Identical heuristics and all five constants; JS drops the unused `qmax` argument |
| Algebraic conversions | named functions | inlined at call sites | Identical arithmetic |
| `enforce_low_r` | present, and the **only** function `scale_pipeline` calls | **absent** | Equivalent behaviour via `firstPeakZero` with a degenerate peak window — but see Step 9: the general window is reachable only from the CLI and the browser |
| `density_line` | present | **absent** | Value recomputed inline in `diagnosticsSummary` |
| `first_peak_zero` | present; used by `scaling_cli.py`, never by `scale_pipeline` | `firstPeakZero`; used by `autoScaleWorker.js` | Identical function; different callers |
| Config knobs `c2_bins`, `c1_slope_nuisance` | present (defaults 0 / `False`) | **absent** | Results agree only because the defaults are inert |

Parity is enforced by golden-fixture tests: `tests/generate_autoscale_fixture.py` runs the Python
engine on a shared synthetic model and writes
`web_app/frontend/src/__tests__/fixtures/autoscale_fixture.json`;
`web_app/frontend/src/__tests__/autoScale.test.js` (vitest) then asserts against it. The
transform-layer-relevant tolerances:

- fitted $a$, $b$: relative error < `1e-6`
- pre-enforcement low-$r$ rms: relative error < `1e-5` (auto) / `1e-6` (manual)
- high-$Q$ tail mean: relative error < `1e-8`
- sampled $G_K(r)$ and $S_\mathrm{filt}(Q)$ values: `toBeCloseTo(..., 9)` — 9 decimal places
- the $\rho_0$ fixed-point estimate: relative error < `1e-4`, explicitly documented in the test
  as "bounded by the `rtol=1e-3` stopping rule, not by single-pass transform precision"

The test file's own header states the tolerance philosophy: "loose enough for summation-order
float noise (numpy pairwise vs JS sequential) and tight enough that any real math drift fails."

**Which fixture case sees the `s0Target` divergence.** Only one. The `auto`, `manual` and
`autoDetect` configs set no `b_sq_avg`, so $s_0^\mathrm{target} = 0$ in both engines and the Step-7
`s0Target` divergence cannot show up in the `1e-6`/`1e-8`/9-decimal tolerances above. The `fz` case
does set `b_sq_avg`, but `amplitude_criterion='fz'` returns from `autoscale` before the iteration,
so it never enters the loop. The `rho0Estimate` case **does** exercise it: it passes
`bSqAvg = fixture.fzBSqAvg` and seeds `rho0 = 0.02` away from the model's 0.05, so the fixed-point
loop runs with a non-zero $s_0^\mathrm{target}$ — which is precisely why its tolerance is relative
error < `1e-4` on $\rho_0$ and $\lvert\text{concordance} - 1\rvert < 1.5\times10^{-3}$, two orders
looser than the `1e-6` asserted for $a$ and $b$, while `converged` and `iterations` must still match
exactly (`autoScale.test.js` → *"rho0 self-consistency matches Python and recovers the density"*).

So the composition-first path is covered **only** through `estimateRho0`. No fixture case runs the
main `autoscale` auto path with a non-zero $s_0^\mathrm{target}$, which is the configuration the Auto
StoG page actually uses once a formula is entered.

### What the tests assert against real data

`tests/test_transforms.py::StogRunParityTests` reproduces the intermediates of a complete Fortran
`stog` run (POWGEN run 59438, Mn₃Sn; `data/stog_tests/stog_59438/`). **The data directory is
gitignored, so these tests skip in CI** and only run on a machine that has the fixture:

| Comparison | Tolerance asserted |
| --- | --- |
| Scaled $S(Q)$ vs `scale.fq` | `atol=1e-9` |
| Forward transform vs `scale.gr` | $\mathrm{rms}/\max\lvert\mathrm{ref}\rvert < 5\times10^{-3}$ |
| Filter correction vs `ft.dat` | rms < 2×10⁻³ |
| Filtered $S(Q)$ vs `scale_ft.sq` | rms < 2×10⁻³ |
| Enforced $G_K$, $D$ below the cutoff vs `scale_ft_rmc.{gr,dr}` | `atol=1e-9`; and $-\langle b\rangle^2$ / the density line to `atol=1e-12` |

Two qualifications that the tolerances alone do not convey:

- **These are post-interpolation numbers.** The repo's outputs and the Fortran files are not on
  identical grids, so every comparison first resamples one side onto the other with linear
  `np.interp` (`np.interp(ref[0], self.q, self.sq_scaled)`,
  `np.interp(self.r, ref[0], ref[1])`, and so on). Linear resampling smooths, so sub-grid-scale
  disagreement can be hidden; the enforcement comparison additionally interpolates a
  *piecewise-constant* reference. The `ft.dat` comparison also restricts to
  $[q_\mathrm{min}, q_\mathrm{max}]$ before computing the rms, discarding the sub-$Q_\mathrm{min}$
  stub rows the Fortran writes and this repo does not (`scaling_cli.py`'s `_FT_NAME` comment:
  "the Fortran's extra sub-Qmin stub carries no data").
- **These tolerances are for the *uncorrected* transform.** `StogRunParityTests` calls the
  `transforms.py` functions bare — `fq_to_gpdf(self.q, sq_to_fq(...), self.r)` and
  `fourier_filter(self.q, self.sq_scaled, self.r, rho0=…, cutoff=…)` — so the function defaults
  apply: **no Lorch window, no low-$Q$ correction, $s_0^\mathrm{target} = 0$**. That is the right
  choice for classic-`stog` parity (the Fortran does none of them either), but it means the
  5×10⁻³ / 2×10⁻³ figures do **not** validate the correction the page applies by default.

The independently-executed `pystog` 0.6.7 cross-validation recorded in
`docs/STOG_SCALING_PLAN.md` §3.1 (run *before* this code existed) reported rms 6.1×10⁻⁴ for the
filter section vs `ft.dat` and 1.4×10⁻¹³ for the scaling — that 6×10⁻⁴ figure is the one quoted
in the `transforms.py` module docstring, and it comes from the cross-run, not from this repo's
test suite (whose committed thresholds are the looser 2×10⁻³ above).

---

### Caveats / what this is not

- **Trapezoid quadrature, not an FFT.** Accuracy is set by the data's own $\Delta Q$ and by
  $\Delta r$; there is no zero-padding, no interpolation, and no spectral-accuracy claim. Cost is
  quadratic in grid sizes.
- **Hard truncation at $Q_\mathrm{max}$.** Nothing corrects the missing
  $(Q_\mathrm{max}, \infty)$ tail. Termination ripple is real and visible in $g(r)$; the only
  mitigation offered is
  the Lorch window, which trades ripple for resolution and is **off by default**.
- **The low-$Q$ correction is a model, not data.** It assumes $S(Q)$ is *linear* on
  $[0, Q_\mathrm{min}]$ and that $S(0)$ equals a target you supply (or that the composition
  implies). Real $S(Q)$ below $Q_\mathrm{min}$ can contain genuine structure — for the Mn₃Sn run
  (`Qmin = 1.0` Å⁻¹, crystalline) `docs/STOG_SCALING_PLAN.md` §3.1b records that the dataset is
  **missing $O(1)$ structure below $Q_\mathrm{min}$, so no affine $(a,b)$ can satisfy the
  absolute normalization criteria**; the analytic correction models that interval as a straight
  line and cannot recover it. (Do not confuse that $O(1)$ deficiency with the Keen Eq. 21 *target*
  $S(0) = 1 - \langle b^2\rangle/\langle b\rangle^2 \approx -12$ for the same composition, quoted
  in Step 5 — the two are different quantities.) A wrong composition (wrong
  $\langle b^2\rangle$) moves the target and biases the low-$r$ region directly.
- **The $r$ grid is not tied to the data.** `r_grid` depends only on `rmax`/`nr` and is never
  checked against $\Delta Q$ or $Q_\mathrm{max}$. With the defaults $\Delta r = 0.01$ Å while the
  transform's resolution is $\pi/Q_\mathrm{max} \approx 0.11$ Å at $Q_\mathrm{max} = 28$ Å⁻¹, so
  the plotted curves are ~11× oversampled and neighbouring points are **not independent** —
  structure at the grid scale is interpolation, not information. Likewise $r_\mathrm{max} = 50$ Å
  is only safe while the aliasing period $2\pi/\Delta Q$ exceeds it (628 Å at
  $\Delta Q = 0.01$ Å⁻¹, but only 63 Å at $\Delta Q = 0.1$ Å⁻¹); coarse-$\Delta Q$ data can wrap.
- **Despiking, if enabled, deletes measured points.** It is OFF by default, but when on it
  removes rows from the array that is then integrated — 12% of points on the crystalline 59438
  benchmark, i.e. real Bragg maxima — and leaves gaps the trapezoid rule bridges with one wide
  panel. Check the reported `n_despiked` before trusting a despiked run.
- **`transforms.py`'s defaults are not the app's defaults.** The reference API is uncorrected and
  un-windowed (`low_q_correction=False`, `s0_target=0.0`, `lorch=False`); the ON policy is a
  `ScalingConfig`/`defaultConfig` decision. Numbers reproduced by calling the reference API
  directly will not match the page unless you pass the flags yourself.
- **The Python and JavaScript engines are known to disagree in one place.** Inside the auto-scale
  iteration the JS filter call omits `s0Target` (Step 7), so the fitted $(a,b)$ differ between the
  CLI/API and the page whenever a composition supplies $\langle b^2\rangle$. The golden fixture
  does not cover that case.
- **The Fourier filter injects an assumption.** After filtering, $S(Q)$ carries the constraint
  "$g \equiv 0$ below $r_\mathrm{cut}$". It is a data-conditioning step, not a measurement, and
  it does part of the normalization work in the classic workflow.
- **Enforcement makes the Keen limits true by construction.** Any assessment of fit quality must
  use the reported pre-enforcement residual. In `'auto'` mode the boundary itself comes from a
  five-constant heuristic (`detect_first_peak_onset`, Step 9) with no test pinning its
  sensitivity — and if that heuristic returns `None`, the browser applies **no enforcement at
  all** while the checkbox still reads as checked.
- **No uncertainty propagation through the transforms.** Per-point $\sigma$ from the data file is
  used to weight the high-$Q$ fit rows only; it is never transformed into real space, and no
  error bars are produced on $g(r)$, $G_K(r)$, or $D(r)$.
- **Constant, $Q$-independent scattering lengths.** The whole layer assumes neutron data with
  real constant $b$. X-ray data (with $Q$-dependent $f(Q)$ and Compton scattering) is explicitly
  labelled out of scope in `docs/STOG_SCALING_PLAN.md` §5; nothing in the transform layer models
  it. No Placzek/inelastic, absorption, or multiple-scattering correction is applied here — the
  input is assumed already reduced.
- **$\rho_0$ is load-bearing and unvalidated by the transforms.** It enters
  $g \leftrightarrow G_\mathrm{PDF}$ and $D(r)$ linearly; a wrong `NUMBER_DENSITY ::` yields a
  confidently wrong low-$r$ target line and hence a wrong scale.
- **The equation numbers are the code's own citations, and the code is not self-consistent.**
  Keen 2001 is not vendored in the repo, so nothing here is verified against the paper.
  `transforms.py` cites **Eq. 9, 10, 15, 16, 19, 21, 29 and a "Eq. 43/44 lineage"** for
  $G_\mathrm{PDF}$ — the very quantity the whole transform layer works in; `scaling.py` and
  `scattering.py` additionally cite **Eq. 14** (the $Q\to0$ Faber-Ziman diagnostic). The
  labelling conflicts with itself: the module docstring calls $F_K = \langle b\rangle^2[S-1]$
  "Eq. 9" and $G_K = \langle b\rangle^2[g-1]$ "Eq. 10", while `sq_to_fk` and `g_to_gk` label the
  **same two relations** "Keen Eq. 19 inverted" and "Keen Eq. 16 inverted". Step 1 prints those
  two formulas; which number they carry depends on which docstring you read. **Eq. 11 and Eq. 12
  appear nowhere in the repository** — the Fourier-pair derivation in Step 3 is this document's
  own reconstruction and is labelled as such. Treat every equation number in this section as a
  pointer into the source, not as a verified citation.


## Auto StoG — the auto-scaling engine: level sweep, closed-form fit, self-consistency

### What this page shows

The **Auto StoG** tab replaces the classic `stog` program's interactive *"try again"* loop — where a
human repeatedly guesses a `yscale`/`yoffset` pair, transforms, looks at the low-$r$ region, and
guesses again — with a determined estimator. You upload a measured total-scattering $S(Q)$ (plus,
optionally, a classic `stog.inp`), give a composition and a $[Q_\mathrm{min}, Q_\mathrm{max}]$
window, and press **Auto-scale**. The engine returns the affine correction

$$S_\mathrm{corr}(Q) \;=\; a\,S_\mathrm{meas}(Q) \;+\; b$$

together with the Fourier-filtered $S(Q)$, the RMCProfile-ready $F_K(Q)$, $G_K(r)$, $D(r)$, and a
set of diagnostics that say *how much the data actually determined the answer*.

Two runtimes compute the same thing. The Python package
[`rmc_toolkits/scaling.py`](../../rmc_toolkits/scaling.py) is the source of truth (used by the
`rmc-autoscale` CLI and the Flask `/api/scaling/*` endpoints); the browser worker
[`web_app/frontend/src/workers/autoScale.js`](../../web_app/frontend/src/workers/autoScale.js) is a
hand-written port that the Auto StoG page uses in **both** Flask and static mode, so files never
leave your machine. Section [Python ↔ JavaScript parity](#python--javascript-parity-1) states exactly
where they agree and where they do not.

#### Notation and units

| Symbol | Meaning | Units |
| --- | --- | --- |
| $Q$ | momentum transfer | Å⁻¹ |
| $Q_0$ | the **first** point of the cropped $Q$ grid (`q[0]`) | Å⁻¹ |
| $Q_{N-1}$ | the **last** point of the cropped $Q$ grid (`q[-1]`) — *not* `config.qmax` | Å⁻¹ |
| $r$ | pair separation | Å |
| $S_\mathrm{meas}(Q)$ | the measured (mis-scaled) structure factor as read from file | dimensionless |
| $a, b$ | multiplicative scale and additive offset of the correction | dimensionless |
| $L$ | the measured high-$Q$ level (asymptote) of $S_\mathrm{meas}$ | dimensionless |
| $\rho_0$ | atomic number density | Å⁻³ |
| $\langle b\rangle^2 = (\sum_i c_i b_i)^2$ | classic stog "Faber-Ziman coefficient" | barn (1 barn = 100 fm²) |
| $\langle b^2\rangle = \sum_i c_i b_i^2$ | a *different* number; only equal to $\langle b\rangle^2$ for one element | barn |
| $F(Q) = Q[S(Q)-1]$ | PDF-community $F(Q)$ used inside the transforms | Å⁻¹ |
| $F_K(Q) = \langle b\rangle^2[S(Q)-1]$ | Keen Eq. 9/19 | barn |
| $G_\mathrm{PDF}(r) = 4\pi\rho_0 r[g(r)-1]$ | PDFFIT-style intermediate | Å⁻² |
| $G_K(r) = \langle b\rangle^2[g(r)-1]$ | Keen Eq. 10/16; flat $-\langle b\rangle^2$ below $r_0$ | barn |
| $D(r) = 4\pi\rho_0 r\,G_K(r)$ | Keen Eq. 29 | barn Å⁻² |
| $r_0$ | closest interatomic approach | Å |
| $\delta(Q) = S_\mathrm{ft}(Q) - 1$ | the Fourier filter's subtraction term (classic `ft.dat`, shifted) | dimensionless |

$\langle b\rangle^2$ and $\langle b^2\rangle$ must be given in the same unit system as each other;
whichever you choose fixes the units of $F_K$, $G_K$, $D$. For normalized x-ray $S(Q)$ the
convention is $\langle b\rangle^2 = 1$ and
$\langle b^2\rangle = \langle Z^2\rangle/\langle Z\rangle^2$.

**Reference:** D. A. Keen, *J. Appl. Cryst.* **34**, 172 (2001). Equation numbers below are that
paper's, *as cited by the code* — see the note on Eq. 21 in Step 2.

#### The Fourier pair, written out

Every transform in this section is one of exactly two operations, and the pair is **asymmetric** —
the $r\to Q$ direction carries no $2/\pi$ prefactor:

$$G_\mathrm{PDF}(r) \;=\; \frac{2}{\pi}\int F(Q)\,\sin(Qr)\,dQ, \qquad
F(Q) \;=\; \int G_\mathrm{PDF}(r)\,\sin(Qr)\,dr .$$

Both are evaluated by the **trapezoid rule on the supplied grids** — no resampling, no analytic
quadrature — which matches the pystog `Transformer` discretization (and a classic Fortran stog run
to ≈6·10⁻⁴ rms). Write $T[\cdot]$ for the forward (first) operation.

**Code:** [`rmc_toolkits/transforms.py`](../../rmc_toolkits/transforms.py) → `sine_transform()`,
`fq_to_gpdf()` (forward, applies the $2/\pi$), `gpdf_to_fq()` (backward, returns
`sine_transform(r, gpdf, q)` **unscaled**). JS: `sineTransform()`, `fqToGpdf()`; the JS engine has
no separate `gpdfToFq` — the backward transform is `sineTransform` called directly inside
`fourierFilter()`.

**The Lorch window** (`lorch=False` by default) multiplies $F(Q)$ before the forward transform only:

$$M(Q) \;=\; \frac{\sin(\pi Q/Q_{N-1})}{\pi Q/Q_{N-1}}, \qquad M(0)\equiv 1 .$$

Note the $Q$ in the denominator: `fq_to_gpdf()` calls `lorch_window(q, q[-1])`, i.e. the window is
built from the **last cropped data point**, not from `config.qmax`. If you crop at
$Q_\mathrm{max} = 30$ Å⁻¹ but your data end at 29.40 Å⁻¹ (as on the parity fixture) the window is
the 29.40 Å⁻¹ one. The same $Q_{N-1}$ sets $a = \pi/Q_{N-1}$ in the Lorch low-$Q$ correction basis
(Step 5b).

#### Composition → $\langle b\rangle^2$, $\langle b^2\rangle$, and $\rho_0$

Composition-derived values come from [`rmc_toolkits/scattering.py`](../../rmc_toolkits/scattering.py) →
`faber_ziman()`, ported to JS as `faberZiman()` in
[`autoScale.js`](../../web_app/frontend/src/workers/autoScale.js). `parse_formula()` /`parseFormula()`
turn a formula string into element counts $n_i$ (decimal stoichiometries and parentheses supported:
`Sr0.5Ba0.5TiO3`, `Al2(SO4)3`). Then, with mole fractions $c_i$ and bound coherent lengths $b_i$
**in fm**:

$$c_i = \frac{n_i}{\sum_j n_j},\qquad
\langle b\rangle = \sum_i c_i b_i,\qquad
\langle b\rangle^2 = \Big(\sum_i c_i b_i\Big)^2,\qquad
\langle b^2\rangle = \sum_i c_i b_i^2 .$$

Both are returned in fm² **and** in barn, the barn value being the fm² value divided by 100
(`FaberZiman.b_avg_sq_barn`, JS `bAvgSqBarn`/`bSqAvgBarn`). **The page always uses the barn form.**

- **Table:** 89 elements, bound coherent lengths in fm from the NIST NCNR table (after Sears,
  *Neutron News* **3**, 26 (1992)). For the seven strong absorbers with complex $b$ — **B, Cd, Dy,
  Eu, Gd, In, Sm** — the table stores the **real part** only (`COMPLEX_B_ELEMENTS`), the standard
  Faber-Ziman choice. `b_overrides_fm` / `bOverridesFm` substitutes per-element values (isotopic
  enrichment).
- **Hard stop:** if $\langle b\rangle^2 < 10^{-4}\,\langle b^2\rangle$ the call **raises**
  ("near-null-matrix composition") — this is how Ti/Zr-style null matrices are rejected rather than
  silently producing an undefined normalization.

**Mass density instead of $\rho_0$.** The page accepts a mass density and converts it
(`number_density_from_mass_density()`, JS `numberDensityFromMassDensity()`), the ADDIE convention:

$$\rho_0\,[\text{Å}^{-3}] \;=\; \rho_m\,[\text{g cm}^{-3}]\times\frac{N_A}{10^{24}}\times
\frac{n_\mathrm{atoms}}{M},\qquad N_A = 6.02214076\times10^{23},$$

with $M$ the molar mass (g mol⁻¹) of one formula unit from the CIAAW `ATOMIC_MASS_U` table (the same
89 elements) and $n_\mathrm{atoms}$ that unit's atom count — the choice of formula unit cancels.
A non-positive mass density raises. The inverse `mass_density_from_number_density()` /
`massDensityFromNumberDensity()` is also exposed (the page uses it nowhere; the vitest suite
round-trips Mn₃Sn at 0.063049 Å⁻³ ↔ 7.4209 g cm⁻³).

---

### Step 1 — From files to arrays: parse, resolve, crop, (optionally) despike

**Inputs:** the uploaded text files. **Outputs:** cropped $(Q, S_\mathrm{meas}, \sigma)$ arrays, a
validated `ScalingConfig`, `n_despiked`.

#### 1a. Reading the files

Three readers turn uploads into numbers. They exist in both runtimes
(`rmc_toolkits/parsers.py` → `read_stog_xy()`, `read_stog_inp()`, `read_dat_header()`; JS
`readStogXy()`, `readStogInp()`, `readDatHeader()` in `autoScale.js`).

**`read_stog_xy` — the data file.** This is a *filtering* operation on your data, not a plain load:

- each line is split on whitespace; lines with fewer than 2 tokens are skipped;
- a line where **any** token fails to parse as a number is skipped — this is how titles and the
  classic leading point-count line are dropped;
- surviving rows are grouped by **column count** and only the **largest group** is kept, so a stray
  numeric header line (e.g. `count qmin` on one line) cannot become the column template and silently
  discard every real row;
- `nan` / `inf` tokens are *preserved*: rebin padding survives to `crop_sq`, which drops it there;
- column 0 = $Q$, column 1 = $S_\mathrm{meas}$, column 2 (when present) = $\sigma$.

> Divergence: the JS reader additionally rewrites Fortran `D` exponents (`1.5D+02` → `1.5e+02`) via
> `tokenToFloat`. The Python reader calls bare `float()`, which rejects `D` — so a `D`-exponent row
> is **skipped** by Python and **kept** by JS. Verified: `float("1.5D+02")` raises.

**`read_stog_inp` — the classic input.** Read by **fixed line index**, requires ≥22 non-empty lines,
and hard-rejects the variants whose layout it cannot know: `n_files ≠ 1`, a nonzero $Q$ offset
(line 5), a nonzero second $y$ offset (line 12), an enabled interactive "try again" flag (line 13), a
disabled filter (line 14), and a zero/non-finite `yscale`. What it supplies to the run:
`qmin`/`qmax` (line 3), `yoffset`/`yscale` (line 4) → $(a, b)$, `rmax` (8), `nr` (9), `lorch` (10),
$\rho_0$ (11), `r_cutoff` (15), $\langle b\rangle^2$ (18), and the enforcement triple
`peak_cutoff`/`peak_rmin`/`peak_rmax` (22). Python and JS validate identically.

**`read_dat_header` — `KEY :: value` metadata.** Lines containing `::` are split into an
upper-cased key and a value. `TITLE` is kept as-is; `NUMBER_DENSITY` contributes its **first finite
token** as $\rho_0$; `MINIMUM_DISTANCES` contributes the **minimum** of its finite values as $r_0$.
Both of these are silent auto-fills of physics parameters — see 1b.

#### 1b. Where each parameter actually comes from

The Auto StoG page resolves the run's inputs from several sources. Your data can be scaled with a
$\rho_0$ or $r_0$ you never typed. The precedence (`AutoStogPage.jsx` → `resolveConfig()`, mirroring
the CLI):

| Quantity | Resolution order (first hit wins) |
| --- | --- |
| $\langle b\rangle^2$ | form override → `stog.inp` line 18 → composition (barn). Absent ⇒ hard error. |
| $\langle b^2\rangle$ | form override → composition (barn). Absent ⇒ FZ mode, `estimate_rho0` and the C3 diagnostic are unavailable. |
| $\rho_0$ | form → `stog.inp` → the data file's `NUMBER_DENSITY ::` header → mass density + composition → **seed 0.05 Å⁻³ with `estimateRho0` switched on**. If even $\langle b^2\rangle$ is missing at that last step the run refuses with "number density unknown". |
| $r_0$ | form → the **minimum** of the header's `MINIMUM_DISTANCES ::` list → $\max(\texttt{peak\_cutoff}, \texttt{peak\_rmin})$ from `stog.inp`, but only when that candidate $-\,0.25 > r_\mathrm{cutoff}+0.2$ (otherwise it would produce an empty fit window). |
| $Q_\mathrm{min}$, $Q_\mathrm{max}$ | the form field → `stog.inp`. On upload, `selectSource()` pre-fills the form from the `stog.inp`, or from the data's finite extent rounded to **4 significant digits**. |
| `r_cutoff`, `rmax`, `nr`, `lorch`, enforcement cutoff, manual $(a,b)$ | form → `stog.inp` → package default. |

Two consequences worth stating plainly: a run with an empty $\rho_0$ field and a `stog.inp` or a
`NUMBER_DENSITY ::` header present does **not** trigger the self-consistent estimate; and the
`0.05 Å⁻³` seed is a placeholder that only ever survives if the estimate is then adopted (Step 10).
The UI warns (`⚠` chip) when a typed $\langle b\rangle^2$/$\langle b^2\rangle$ override differs from
the composition-derived value by more than **2 %** — the stale-session trap.

#### 1c. Crop

`crop_sq()` keeps rows satisfying all of

$$\mathrm{finite}(Q)\;\wedge\;\mathrm{finite}(S)\;\wedge\;Q > 0\;\wedge\;Q \ge Q_\mathrm{min}-10^{-12}\;\wedge\;Q \le Q_\mathrm{max}+10^{-12}.$$

$Q\le 0$ rows are dropped **unconditionally** — every $S\leftrightarrow F$ conversion divides by
$Q$, and the analytic low-$Q$ correction (Step 5b) already models the omitted $[0, Q_\mathrm{min}]$
range, so a $Q=0$ point carries no usable information. Fewer than 16 surviving points is a hard
error. NaN padding below the data's true $Q_\mathrm{min}$ (normal in rebinned STOG files) is removed
here.

> The ≥16 check runs on the cropped array **before** despiking and is never re-applied. With
> `despike=True` a run can therefore proceed on fewer than 16 points; what actually catches an
> over-aggressive despike downstream is the level sweep's ≥32-point requirement or the ≥2-points-per
> -fit-window check, with a different error message.

**Grid assumptions (both engines, never checked).** The input $Q$ column is used **in file order** —
no sorting, deduplication, interpolation or rebinning is performed anywhere. The window sweep's
prefix sums, both transforms' trapezoid rule, the low-$Q$ correction hook ($Q_0 =$ first cropped
point), the Lorch $Q_{N-1}$ (last cropped point) and the FZ head selection all assume a strictly
ascending grid. On a non-monotonic file the results are undefined, and the two engines even
disagree: the sweep's failure path uses `np.median(sq[q >= q.max() - min_width])` in Python and
`q[n-1]` in JS.

The $r$ grid is built independently of the data: `ScalingConfig.r_grid` is
$r_i = i\cdot (r_\mathrm{max}/n_r)$ for $i = 1\ldots n_r$ (defaults $r_\mathrm{max}=50$ Å,
$n_r=5000$, so
$\Delta r = 0.01$ Å). **$r=0$ is deliberately excluded** — the $g = G_\mathrm{PDF}/(4\pi\rho_0 r)+1$
conversion divides by $r$.

#### 1d. Optional despiking

`despike=False` by default. `_despike_mask()` computes a rolling median over `despike_window = 7`
points (edge-padded), forms $\varepsilon_i = S_i - \mathrm{med}_i$, a global scale
$s = 1.4826\,\mathrm{median}(|\varepsilon|)$ (note: *not* median-centred), and keeps points with
$|\varepsilon_i| \le 6.0\,\max(s, 10^{-12})$ (`despike_nsigma = 6.0`). This exists because
detector-glitch spikes ring through the sine transform into the low-$r$ window — a channel Huber
re-weighting cannot reject. It is **off by default** because it also flags genuine Bragg maxima: the
class docstring records 12 % of points removed on the crystalline POWGEN 59438 benchmark.

`despike_window` must be **odd**, and *nothing validates it*. On an even value Python's rolling
median returns $n+1$ values against $n$ points and `_despike_mask` dies with a numpy broadcast error;
JS clamps indices instead and silently uses an off-centre window. (Verified on both.)

**What `n_despiked` actually counts — a defect worth knowing.** In `autoscale()` the data are
despiked **twice**. `_autoscale_pass()` calls `crop_sq()` (which despikes) and hands those
*already-despiked* arrays to `scale_pipeline()`, which crops and despikes them **again** and computes
`n_despiked` from that second pass. Measured on the repo's own despike fixture (synthetic + 8
injected spikes, `despike_window=7`, `despike_nsigma=6.0`): 2941 points survive the crop, the first
pass removes **107**, the second removes a further **75**, and `result.provenance["n_despiked"]`
reports **75** — neither the total removed nor the number of injected spikes. A second consequence:
the fit runs on singly-despiked data while the published arrays come from doubly-despiked data, so
they can sit on slightly different $Q$ grids. The JS engine does not do this (`autoscalePass()`
passes the *raw* arrays to `scalePipeline()`), so it despikes once and its `nDespiked` is the true
count — see [parity](#python--javascript-parity-1) item 6.

`tests/test_scaling.py` → `test_despike_restores_recovery_under_tail_glitches` asserts that **at
least 8 points are dropped** (`n_despiked >= 8` — a loose lower bound on the total, not an assertion
that the 8 injected spikes were identified; the actual reported value on that fixture is 75), that
the recovered scale returns to within 2 % of truth, and that the residual scale error is ≤ 20 % of
the un-despiked run's.

**Code:** `rmc_toolkits/scaling.py` → `crop_sq()`, `_despike_mask()`; JS `cropSq()`,
`despikeKeepMask()`.

---

### Step 2 — The model and its physics targets

The correction is affine, and the convention is **multiply**:

$$S_\mathrm{corr}(Q) = a\,S_\mathrm{meas}(Q) + b .$$

Classic `stog_new3` *divides* (`S/yscale + yoffset`), pystog multiplies. The `stog.inp` reader
converts on the way in: $a = 1/\texttt{yscale}$, $b = \texttt{yoffset}$ (the validation run's
`-9 0.1` line ⇒ $a=10$, $b=-9$). Everything in the API and the page uses the multiply form.

Two independent physics targets pin $(a, b)$, and a third is reported:

- **C1 — high-$Q$ asymptote:** $S_\mathrm{corr}(Q\to\infty)\to 1$. Evaluated on a tail
  window $Q \ge Q_\mathrm{max} - f_\mathrm{tail}(Q_\mathrm{max}-Q_\mathrm{min})$ with
  $f_\mathrm{tail} = $ `q_tail_frac` $= 0.15$. **The window is built from the *configured* `qmin`
  and `qmax`, not from where the data actually stop** — a $Q_\mathrm{max}$ set above the end of the
  data shifts it and can empty it, which raises "fit windows contain fewer than 2 points".
- **C2 — low-$r$ density limit (Keen Eqs. 15/29 in $g$-space):** $g_\mathrm{corr}(r) \equiv 0$ for
  $r < r_0$, equivalently $G_K \to -\langle b\rangle^2$ and
  $D(r) \to -4\pi\rho_0\langle b\rangle^2 r$. Evaluated on `r_fit_window`
  $= [\,\texttt{r\_cutoff}+0.2,\ r_0 - 0.25\,]$ Å by default (if
  $r_0$ is unknown and un-detected, the upper edge falls back to lower + 1.0 Å).
- **C3 — the $Q\to 0$ limit,** $F_K(0)\to-\langle b^2\rangle$, is **reported, never fitted** in
  density mode; it becomes the *fitting* criterion only in `amplitude_criterion="fz"` (Step 9).

**C3 and the FZ criterion are the same constraint in different units.** Since
$F_K = \langle b\rangle^2(S-1)$,

$$F_K(0) = -\langle b^2\rangle \iff S(0) = 1 - \frac{\langle b^2\rangle}{\langle b\rangle^2},$$

which is why C3 is a passive diagnostic in density mode and *becomes* the amplitude equation in FZ
mode.

> Citation discrepancy to know about (trust the code, not the docstrings): `scaling.py`'s module
> docstring cites **Keen Eq. 21** for the C1 high-$Q$ asymptote, while `amplitude_from_fz_limit()`'s
> docstring cites **Keen Eq. 21** for the $Q\to0$ limit
> $S(0) = 1 - \langle b^2\rangle/\langle b\rangle^2$, and `diagnostics_summary()` cites
> **Keen Eq. 14** for the $F_K(0)\to-\langle b^2\rangle$ form of that same $Q\to0$ statement.
> One equation number is being used for two
> different physics targets. This document does not guess which attribution is correct; check
> against the paper before quoting an equation number from here.

The load-bearing structural fact: because the sine transform is linear in the data, $g$ computed
from $a S_\mathrm{meas} + b$ is **exactly affine in $(a,b)$** — no linearization, no approximation.
Both residual blocks are therefore rows of one linear least-squares system, solved in closed form.
The only iteration in the whole engine comes from the Fourier filter's dependence on the current
scale (Step 7) and from the optional $r_0$ refinement (Step 8).

#### The other path: a fixed-$(a, b)$ "manual" run

The page has a **second, fully separate code path**, reached from *Advanced → Fixed scaling → Run
fixed (a, b)*. It takes $a$ (and optionally $b$) from the form or, failing that, from the loaded
`stog.inp`'s $1/\texttt{yscale}$ and `yoffset` ($b$ defaults to 0 in data mode), and calls
`scale_pipeline()` / `scalePipeline()` **directly**. Steps 3–10 are skipped entirely:
`iterations = 0`, `converged = True` by construction, `sweep = None`, `a_fz = None`,
`c1ModeEffective = "manual"`, `provenance["mode"] = "manual"`. `resolveConfig()` also **overrides
`amplitudeCriterion` to `"density"`** regardless of the UI selection, so a leftover `fz` choice
(which would require $\langle b^2\rangle$ + sweep mode) cannot spuriously reject a run that never
uses an amplitude criterion. Because the manual path never runs `autoscale()`'s detection, the
worker re-runs `detectFirstPeakOnset()` afterwards so that "auto" low-$r$ enforcement is not a
silent no-op (Step 11).

**Code:** module docstring and `ScalingConfig` in `rmc_toolkits/scaling.py`;
`AutoStogPage.jsx` → `resolveConfig()`, `runScaling()`; `autoScaleWorker.js`;
[`docs/STOG_SCALING_PLAN.md`](../STOG_SCALING_PLAN.md) §1.4.

---

### Step 3 — The high-$Q$ level sweep

**Inputs:** cropped $(Q, S_\mathrm{meas})$ — *measured*, unscaled, unfiltered. Per-point $\sigma$ is
**not** used here.

**Question it answers:** "what $Q$ is high enough to read the asymptote?", without a hand-set
flatness tolerance.

**Enumeration.** Candidate window edges are
`edges = np.unique(np.linspace(0, N-1, n_grid).astype(int))` with `n_grid = 80` (truncation, not
rounding). Every ordered pair $(i, j)$ of edges defines a window $[Q_i, Q_j]$. A window is
*considered* only if

$$Q_j - Q_i \;\ge\; \texttt{min\_width} = 3.0\ \text{Å}^{-1}, \qquad n = j-i+1 \;\ge\; 24, \qquad \det > 0 .$$

That is at most $\binom{80}{2} = 3160$ windows, each evaluated in $O(1)$.

**The $O(1)$ line fit.** Six prefix sums over the whole array ($\sum 1$, $\sum Q$, $\sum Q^2$,
$\sum y$, $\sum Qy$, $\sum y^2$ with $y \equiv S_\mathrm{meas}$) make every window's ordinary
least-squares fit of $y = c + m Q$ a handful of subtractions:

$$\det = n\Sigma_{Q^2} - \Sigma_Q^2,\qquad
m = \frac{n\Sigma_{Qy} - \Sigma_Q\Sigma_y}{\det},\qquad
c = \frac{\Sigma_{Q^2}\Sigma_y - \Sigma_Q\Sigma_{Qy}}{\det}$$

$$\mathrm{RSS} = \Sigma_{y^2} - c\,\Sigma_y - m\,\Sigma_{Qy},\qquad
\hat\sigma^2 = \max\!\left(\frac{\mathrm{RSS}}{n-2},\,0\right),\qquad
\mathrm{Var}(m) = \frac{\hat\sigma^2 n}{\det}.$$

The window's **level** is the fitted value at the window centroid $\bar Q = \Sigma_Q/n$, which is
the minimum-variance point of a line fit:

$$L_\mathrm{win} = c + m\bar Q, \qquad \mathrm{Var}(L_\mathrm{win}) = \frac{\hat\sigma^2}{n}.$$

**Admissibility (the criterion that replaces a tolerance).** A window is admissible iff its slope is
statistically indistinguishable from zero *given its own residual scatter*:

$$|m| \;<\; \texttt{slope\_nsigma}\cdot\sqrt{\mathrm{Var}(m)}, \qquad \texttt{slope\_nsigma} = 2.0 .$$

The data — not the user — define "flat". End-of-range artifacts (detector rolloff, dead tails)
exclude themselves: any window touching them acquires a significant slope.
`tests/test_scaling.py` → `test_excludes_dead_tail` injects an exponential rolloff above 24 Å⁻¹ and
asserts the selected window ends at or below 24.5 Å⁻¹.

**Selection and spread.** Among admissible windows the engine keeps the one with the smallest
$\mathrm{Var}(L_\mathrm{win}) = \hat\sigma^2/n$ — the longest, quietest stretch. The reported
uncertainty is the **population standard deviation of $L_\mathrm{win}$ over *all* admissible
windows** (`np.std(admissible_levels)`, `ddof=0`), not a formal standard error. Those windows
overlap heavily, so this is a spread-of-estimates statistic, not an independent-sample one; read it
as "how much the answer moves if you had chosen a different reasonable window".

**Failure path.** If no window is admissible, `LevelSweepResult` is returned with
`asymptote_found = False`, `level_uncertainty = NaN`, `n_admissible = 0`, and

$$L = \mathrm{median}\big(S_\mathrm{meas}\ \text{over the last } 3.0\ \text{Å}^{-1}\big),\qquad
q_\mathrm{lo} = Q_\mathrm{last} - 3.0,\qquad q_\mathrm{hi} = Q_\mathrm{last}.$$

Those `q_lo`/`q_hi` are **fabricated, not searched** — they are simply the last `min_width` of data.
Always read `level_window` together with `asymptote_found`. `autoscale()` then does **not** anchor
the offset and silently falls back to the joint 2-dof fit
(`provenance["c1_mode_effective"] = "joint"`). `tests/test_scaling.py` →
`test_reports_no_asymptote_on_sloped_data` pins this on uniformly sloped data.

**Outputs:** `LevelSweepResult(level, level_uncertainty, q_lo, q_hi, slope, slope_sigma,
asymptote_found, n_admissible)`.

**Code:** `rmc_toolkits/scaling.py` → `level_sweep()`; JS `levelSweep()`. `min_width`, `n_grid` and
`slope_nsigma` are **hard-wired call defaults** — neither `ScalingConfig` nor the UI exposes them.

> Honest caveat: the admissibility test scales with the data's own noise. On a very noisy tail
> $\sqrt{\mathrm{Var}(m)}$ is large and a physically sloped window can pass (it then usually loses
> the minimum-variance contest); on very precise data $\sqrt{\mathrm{Var}(m)}$ shrinks like
> $n^{-3/2}$ and a physically negligible drift can disqualify a long window. On the repo's synthetic
> benchmark the winning window ran $Q\in[23.91, 29.40]$ Å⁻¹ with $|m|/\sigma_m = 1.996$ — right at
> the boundary — out of 797 admissible windows.

---

### Step 4 — Anchoring the offset: $b = 1 - aL$

With a level $L$ in hand, the offset stops being a free parameter. Requiring
$S_\mathrm{corr}\to 1$ where $S_\mathrm{meas}\to L$ gives

$$aL + b = 1 \quad\Longleftrightarrow\quad b = 1 - aL \quad\Longleftrightarrow\quad
S_\mathrm{corr} = a\,(S_\mathrm{meas} - L) + 1 .$$

This is the "shift by the measured level, then scale" decomposition; it leaves exactly **one**
degree of freedom, the amplitude $a$. (`c1_mode="sweep"` is the default; `c1_mode="joint"` restores
the original 2-dof fit, and `fit_offset` — which freezes $b=0$ — only has meaning in joint mode.
`tests/test_scaling.py` → `test_recovers_scale_with_offset_disabled` documents that.)
`test_sweep_mode_recovery_and_provenance` asserts $b = 1 - a L$ to 10 decimal places on the fit
output.

---

### Step 5 — The closed-form affine solve

**Inputs:** cropped $(Q, S_\mathrm{meas})$, the current filter term $\delta(Q)$ (zero on the first
pass), the tail mask, the low-$r$ window mask, optional $\sigma$, optional $L$.

Write the *effective* corrected function that the residuals are built from — the filter's
subtraction term is held **fixed** during the solve, which is what keeps the system linear:

$$S_\mathrm{eff}(Q) \;=\; a\,S_\mathrm{meas}(Q) + b - \delta(Q).$$

#### 5a. C1 rows (one per tail point)

$$S_\mathrm{eff}(Q) - 1 = 0 \;\Longrightarrow\; \underbrace{S_\mathrm{meas}(Q)}_{\text{col }a}\,a
\;+\; \underbrace{1}_{\text{col }b}\,b \;=\; \underbrace{1 + \delta(Q)}_{\mathrm{rhs}} .$$

If per-point $\sigma$ is supplied, both columns and the rhs of this block — **three** columns when
the slope-nuisance column below is enabled, since `col_m_c1` is scaled too — are multiplied by
$w_\sigma = 1/\max(\sigma, 10^{-12})$ **renormalized to unit mean over the block**, so the relative
weight of the C1 block against C2 does not change with the absolute $\sigma$ scale.

$\sigma$ is **not** an opt-in you have to think about: the Auto StoG page ships `useSigma: true`, so
whenever the uploaded data file has a third column it is used automatically (toggle under
*Advanced → Amplitude & offset → σ column*). It is cropped and despiked alongside $Q$ and $S$, and
it re-weights **only the C1 rows** — the level sweep (Step 3), the C2 block (5b) and the FZ head fit
(Step 9) all ignore it.

If `c1_slope_nuisance=True` (experimental, Python only) an extra column $Q - \bar Q_\mathrm{tail}$
is appended, non-zero on C1 rows only. It absorbs a linear tail drift (x-ray $f(Q)$ residuals,
Placzek leftovers) so the drift cannot bias the level. Its fitted coefficient is discarded. A drift
spanning the whole $Q$ range also enters through the transform, which this does **not** correct.

#### 5b. C2 rows (one per $r$ point in the low-$r$ window)

Substituting the affine model into $F = Q(S-1)$:

$$F_\mathrm{eff}(Q) = Q\left(aS_\mathrm{meas} + b - \delta - 1\right)
= a\big[\underbrace{Q(S_\mathrm{meas}-1)}_{F_\mathrm{data}} + Q\big] + bQ - Q\delta - Q .$$

With $T[\cdot]$ the forward transform defined in [Notation](#the-fourier-pair-written-out), three
basis transforms are computed **numerically on the same $Q$ window, grid and window function as the
data** — so truncation effects cancel self-consistently rather than being modelled analytically:

$$g_\mathrm{data} = T[F_\mathrm{data}],\qquad g_\mathrm{one} = T[\,Q\cdot 1\,],\qquad g_\Delta = T[\,Q\,\delta\,],$$

giving
$T[F_\mathrm{eff}] = a(g_\mathrm{data} + g_\mathrm{one}) + b\,g_\mathrm{one} - g_\Delta - g_\mathrm{one}$.
All three are evaluated **only at the $r$ points inside the low-$r$ fit window**
(`rw = r[window]`), not on the full $r$ grid — a deliberate cost optimization, and the reason the
refinement pass of Step 8 has to recompute them from scratch when the window moves.

**The omitted-low-$Q$ correction** (`low_q_correction=True` by default) is folded in
*analytically*, not by re-running the basis transforms with the flag on — doing the latter would
multiply-count its constant term. (Note the contrast with the Fourier filter in Step 7, whose own
forward transforms *do* run with the correction on.) It models the unmeasured $[0, Q_\mathrm{min}]$
panel as a straight line from $S(0) = s_0$ to the first measured point $S(Q_0)$, whose exact
contribution is

$$\frac{2}{\pi}\int_0^{Q_0}\! Q\big[S(Q)-1\big]\sin(Qr)\,dQ
= \mathrm{coef}(r)\big[S(Q_0) - s_0\big] + (s_0-1)\,\mathrm{const}_0(r),$$

$$\mathrm{coef}(r) = \frac{2}{\pi}\frac{f_1(r)}{Q_0},\quad
\mathrm{const}_0(r) = \frac{2}{\pi}f_2(r),$$

with $Q_0 \equiv$ the first cropped $Q$ point. Without the Lorch window,

$$f_1 = \frac{2v\sin v - (v^2-2)\cos v - 2}{r^3},\qquad
f_2 = \frac{\sin v - v\cos v}{r^2},\qquad v = Q_0 r .$$

With `lorch=True`, writing $a_L = \pi/Q_{N-1}$ (the *last cropped data point*, not `qmax`) and
$v_\mp = Q_0(r \mp a_L)$:

$$f_1 = \frac{1}{2a_L}\left[\frac{v_-\sin v_- + \cos v_- - 1}{(r-a_L)^2}
- \frac{v_+\sin v_+ + \cos v_+ - 1}{(r+a_L)^2}\right],\qquad
f_2 = \frac{1}{2a_L}\left[\frac{\sin v_-}{r-a_L} - \frac{\sin v_+}{r+a_L}\right].$$

At the removable singularity $r = a_L = \pi/Q_{N-1}$ the code substitutes the analytic limits, with
$v_a = 2a_LQ_0$:

$$f_1^{\lim} = \frac{1}{2a_L}\left[\frac{Q_0^2}{2} - \frac{v_a\sin v_a + \cos v_a - 1}{(2a_L)^2}\right],
\qquad
f_2^{\lim} = \frac{1}{2a_L}\left[Q_0 - \frac{\sin v_a}{2a_L}\right],$$

selected wherever $|r - a_L| \le 10^{-9}\max(1, a_L)$ (`atol` in both engines).

The code stores the correction as $\mathrm{coef}\cdot S(Q_0) - \mathrm{const}$ with
$\mathrm{const} = (1-s_0)\mathrm{const}_0 + s_0\,\mathrm{coef}$ — algebraically identical, and
*affine in $S(Q_0)$*, which is what lets it join the linear system. Since
$S_\mathrm{eff}(Q_0) = a S_\mathrm{meas}[0] + b - \delta[0]$, the correction splits cleanly across
the $a$, $b$ and constant blocks.

> Dead branch, for completeness: `low_q_correction_basis()` returns zero coefficient **and** zero
> constant when `q[0] == 0` (pystog parity — data starting at $Q=0$ omit nothing). Step 1c's
> unconditional $Q\le 0$ drop makes that branch unreachable from this pipeline; it matters only for
> direct calls into `transforms.py`.

Assembling, the residual $g_\mathrm{eff}(r) = G_\mathrm{PDF,eff}(r)/(4\pi\rho_0 r) + 1 = 0$ becomes

$$\boxed{\;
a\,\frac{g_\mathrm{data} + g_\mathrm{one} + \mathrm{coef}\cdot S_\mathrm{meas}[0]}{4\pi\rho_0 r}
\;+\;
b\,\frac{g_\mathrm{one} + \mathrm{coef}}{4\pi\rho_0 r}
\;=\;
\frac{g_\Delta + g_\mathrm{one} + \mathrm{coef}\cdot\delta[0] + \mathrm{const}}{4\pi\rho_0 r} - 1 \;}$$

each row additionally multiplied by $w_2 = \sqrt{\texttt{c2\_weight}}$ (default 1, i.e. C1 and C2
rows carry equal unit weight **per row** — see the block-balance warning in 5c).

Two consequences of the $1/(4\pi\rho_0 r)$ normalization worth knowing:

1. **It is a weighting choice.** Residuals are expressed in $g$-units, so within the window each $r$
   point is weighted by $1/r$ — small-$r$ rows count more than large-$r$ rows.
2. **This is where $\rho_0$ enters the amplitude.** Multiplying through by $4\pi\rho_0 r$ shows the
   C2 block is literally "match the transform to the density line":
   $$a\,A(r) + b\,B(r) - C(r) \;=\; -4\pi\rho_0 r ,$$
   i.e. $D_\mathrm{corr}(r) \to -4\pi\rho_0\langle b\rangle^2 r$. The right-hand side scales
   linearly with $\rho_0$, so the density-limit amplitude responds to a $\rho_0$ error roughly
   1:1 — the degeneracy exploited in Step 10.

**Optional C2 binning** (`c2_bins > 0`, Python only, default 0 = off) replaces the pointwise rows by
`c2_bins` contiguous bin means, each scaled by $\sqrt{\text{bin size}}$ (`_bin_means()`). This is the
classic "level matching" criterion. It is **not** the default: on ripple-heavy real data the mean
level is degenerate with the low-$Q$-hole artifact and the solve runs away
([STOG_SCALING_PLAN §3.1b](../STOG_SCALING_PLAN.md)). Pointwise residuals are the stable choice.

#### 5c. Reduction to the design matrix

Stack the two blocks into columns $\mathbf{A}$ (the $a$ column) and $\mathbf{B}$ (the $b$ column)
and a right-hand side $\mathbf{y}$. Then:

- **Sweep-anchored** ($L$ known): substituting $b = 1 - aL$ collapses the system to a **single
  column**
  $$\big(\mathbf{A} - L\,\mathbf{B}\big)\,a \;=\; \mathbf{y} - \mathbf{B}.$$
  On the C1 rows this reads $a\,(S_\mathrm{meas}-L) = \delta$; on the C2 rows it is the density
  limit. Note that the C1 rows are **not dropped** in sweep mode, contrary to what a quick read of
  the docstring ("leaving the density limit a single amplitude dof") suggests — they simply
  contribute weakly. Measured on the repo's synthetic benchmark: drowning the C1 block by setting
  `c2_weight = 10⁶` changes $a$ by $<10^{-6}$ relative, so the amplitude is in practice set by C2.
- **Joint**: columns $[\mathbf{A}]$, plus $[\mathbf{B}]$ when `fit_offset=True`.
- Plus the slope-nuisance column when enabled (in either mode).

> **The real C1:C2 balance is set by the row counts, and one of them depends on the $r$ grid.**
> The C2 block contributes **one row per $r$-grid point inside
> $[r_\mathrm{fit,min}, r_\mathrm{fit,max}]$** — about 100 rows for a 1 Å window at the
> defaults ($\Delta r = 0.01$ Å) —
> while the C1 block contributes one row per $Q$ point in the top 15 % of the $Q$ window. Changing
> `nr` or `rmax` therefore silently rescales the density-limit block against the tail block and
> changes the fitted $a$, **even at `c2_weight = 1`**. `c2_weight` is a per-row multiplier on top of
> that; it is not the whole story.

The system is solved with `np.linalg.lstsq(..., rcond=None)` (SVD least squares). Then

$$
a = \mathrm{solution}[0], \qquad
b = \begin{cases} 1 - aL & \text{sweep-anchored}\\ \mathrm{solution}[1] & \text{joint, fit\_offset}\\ 0 & \text{joint, no offset.}\end{cases}
$$

**Code:** `rmc_toolkits/scaling.py` → `_solve_affine()`, `_bin_means()`, and
[`rmc_toolkits/transforms.py`](../../rmc_toolkits/transforms.py) → `fq_to_gpdf()`, `sine_transform()`,
`low_q_correction_basis()`. JS: `solveAffine()`, `fqToGpdf()`, `sineTransform()`,
`lowQCorrectionBasis()` — the JS solve uses pivoted normal equations instead of an SVD (see parity
section).

---

### Step 6 — Huber IRLS robust re-weighting

**Why:** a handful of residual Bragg spikes in the tail, or a ripple burst in the low-$r$ window,
would otherwise drag a closed-form least-squares solution.

**Operation** (`robust=True` by default). After the unweighted solve, exactly **three** re-weighting
passes run — a fixed count, with no convergence test and no early exit:

1. Residuals $\mathbf{e} = \mathrm{design}\cdot\mathbf{x} - \mathbf{y}$.
2. Weights computed **separately per block** (C1 rows $[0, n_1)$, C2 rows $[n_1, \cdot)$ — the
   latter only if the C2 block has $\ge 4$ rows, otherwise its weights stay 1), each block using its
   own MAD scale:
   $$s = 1.4826\,\mathrm{median}\big(|e - \mathrm{median}(e)|\big),\qquad
   z_i = \frac{|e_i|}{s},\qquad
   w_i = \min\!\left(1,\ \frac{c}{\max(z_i, 10^{-14})}\right),\qquad c = 1.345 .$$
   $c=1.345$ is the standard Huber tuning constant (95 % efficiency at the Gaussian). If
   $s \le 10^{-14}$ all weights are set to 1.
3. Re-solve with rows and rhs multiplied by $w_i$.

Because rows are multiplied by $w_i$ rather than $\sqrt{w_i}$, the effective quadratic objective is
$\sum_i w_i^2 e_i^2$ — a *more* aggressive down-weighting than textbook Huber IRLS. Both engines do
the same thing, so this is a documented characteristic, not a cross-engine divergence.

**Limitation, stated in the code:** IRLS cannot reject a detector glitch that has already been
transformed — a spike in $S(Q)$ rings across the whole low-$r$ window, so its C2 residuals look like
signal, not outliers. That is what `despike` (Step 1d) is for.

**Code:** `_huber_weights()`, the `if config.robust:` loop in `_solve_affine()`; JS `huberWeights()`.

---

### Step 7 — The self-consistent loop with the Fourier filter

The classic workflow's Fourier filter performs part of the normalization: on the validation run, the
hand-scaled *unfiltered* tail sits at ≈0.94 and the filter's correction lifts it to 1.000. So C1
must be evaluated on **filtered** data — but the filter depends on the scale. Hence a fixed-point
loop:

```
delta = 0
repeat (max_iter = 50):
    (a, b) = solve_affine(..., delta)          # closed form + IRLS, delta frozen
    S_scaled = a*S_meas + b
    (_, S_ft, g_filtered) = fourier_filter(S_scaled)
    delta = S_ft - 1
    history.append((a, b, low_r_rms))
until |a - a_prev| <= tol*max(1,|a|) and |b - b_prev| <= tol*max(1,|b|)
```

with `tol = 1e-6` and `max_iter = 50`. The first pass compares against $\pm\infty$, so convergence
takes a minimum of two iterations; typical observed counts are **3–7** (the parity fixture converges
in 3, the joint-mode variant in 6). If the loop exhausts `max_iter`, `converged=False` is returned
along with the last iterate — the page shows a ✗ and the result is still produced, so *always check
the convergence flag*.

**The filter itself** (`transforms.fourier_filter()`, the pystog/classic-stog recipe), precisely:

1. Raises immediately if $Q_0 \le 0$ ("requires a strictly positive Q grid").
2. Forward-transforms $F = Q(S_\mathrm{scaled}-1)$ to $G_\mathrm{PDF}$ and forms
   $g = G_\mathrm{PDF}/(4\pi\rho_0 r) + 1$ on the **full** $r$ grid. This forward transform is run
   **with** `low_q_correction` and `s0_target` on — in contrast to the three basis transforms of
   Step 5b, which are run with the correction *off* and get it added analytically.
3. Selects the section $r \le \texttt{r\_cutoff}$ (default 1.0 Å). Because the $r$ grid starts at
   $\Delta r$, the "section" is $r\in(0, r_\mathrm{cutoff}]$ *as gridded*, never $r=0$.
4. Back-transforms **only those section points** —
   $F_\mathrm{ft}(Q) = \int 4\pi\rho_0 r\,g(r) \sin(Qr)\,dr$, a trapezoid over the section's
   own short grid — and converts to
   $S_\mathrm{ft} = F_\mathrm{ft}/Q + 1$.
5. Subtracts, $F_\mathrm{filtered} = F - F_\mathrm{ft}$, and re-transforms (again with
   `low_q_correction`/`s0_target` on) to get `g_filtered`.

`sq_ft` is exactly the classic `ft.dat` in $S$-convention, satisfying
`sq_filtered = sq - (sq_ft - 1)`. It was validated against a real Fortran stog run to 6·10⁻⁴ rms.

`history` records $(a, b, \mathrm{low-}r\ \mathrm{rms})$ per iteration. The page's trajectory readout
("a: 9.97 → 9.97 → …") shows **only `history.slice(-6)`**, **only** the $a$ column, formatted to
**5 significant digits**, and is suppressed entirely when the history has fewer than 2 rows — so a
short arrow chain does not mean few iterations. The full history goes into the exported provenance
JSON.

**Code:** `rmc_toolkits/scaling.py` → `_autoscale_pass()`, `_pipeline()`;
`rmc_toolkits/transforms.py` → `fourier_filter()`. JS `autoscalePass()`, `fourierFilter()`;
`AutoStogPage.jsx` `trajectory` memo.

> Cost note: the three basis transforms of Step 5b are recomputed every iteration even though only
> $g_\Delta$ changes. Each transform is $O(N_Q \times N_\mathrm{window})$ and the filter is
> $O(N_Q\times n_r)$ twice per iteration — that dominates the runtime, which is why the browser
> version runs in a Web Worker.

---

### Step 8 — $r_0$ detection from the first-shell $|g|$ flank, and the refinement pass

**Why:** the low-$r$ fit window needs an upper edge below the first coordination shell. Requiring the
user to know $r_0$ would defeat "composition + $Q$ window are the only required inputs".

**Operation** (`detect_first_peak_onset()`): within a search range
$[\texttt{search\_min}, \texttt{search\_max}]$ — called from `autoscale()` with
`search_min = r_cutoff + 0.3` (default
1.3 Å) and `search_max = 6.0` Å — **return `None` immediately if fewer than 3 $r$-grid points fall
inside that range**; otherwise find the index of the maximum of $|g_\mathrm{filtered}|$. If that
peak height is below `floor = 0.5` return `None`. Otherwise set the flank level

$$\ell = \max\big(\texttt{floor},\ \texttt{fraction}\times\mathrm{peak}\big),\qquad
\texttt{fraction} = 0.35,\ \texttt{floor} = 0.5,$$

walk left from the peak while $|g| > \ell$, and return the $r$ of the **last point still above**
$\ell$ (`r[index + 1]` after the walk). If the walk hits the left edge of the search range without
dropping below $\ell$, return `None` ("feature not separable from the ripple field").

Two deliberate design choices:

- **Peak-relative, not absolute.** Both the physical peak and the sub-$r_0$ truncation ripples scale
  with the fitted amplitude; on missing-low-$Q$ data ripples can reach $\sim$ peak/3, so no fixed
  threshold separates them, while the dominant shell still towers above.
- **$|g|$, not $g$.** Faber-Ziman totals of negative-$b$ compositions (Mn₃Sn) have an *inverted*
  first shell.

**The refinement pass.** `autoscale()` runs `_autoscale_pass()` once, detects the onset, records it
as `provenance["r0_detected"]`, and re-runs the entire fit with `r0 = onset` **only if all** of:
`config.r0 is None`, `config.r_fit_max is None`, $(\mathrm{onset}-0.25) > r_\mathrm{fit,min}$, and
$|(\mathrm{onset}-0.25) - r_\mathrm{fit,max}^\mathrm{current}| > 0.05$ Å. The refined result carries
`provenance["window_refined"] = True`. So a full auto-scale can cost **two** complete
self-consistent loops. Note this happens for **both** amplitude criteria — the FZ branch is not
exempt (Step 9).

Measured detections quoted in [SCALING_PROCEDURE.md](../SCALING_PROCEDURE.md): 2.73–2.77 Å across the
Mn₃Sn runs and 2.53 Å for FeCoSn, against hand-chosen classic cutoffs of 2.40–2.68 Å.
`tests/test_scaling.py` → `test_detects_first_shell_and_refines_window` (synthetic onset 2.65, peak
2.80) and `test_autoscale_composition_only_detects_first_shell` (real Mn₃Sn) pin the behaviour.

**Code:** `detect_first_peak_onset()`, `autoscale()`; JS `detectFirstPeakOnset()`, `autoscale()`.
The Python signature keeps an unused `qmax` positional argument (`# noqa: ARG001`); the JS port
drops it.

---

### Step 9 — The alternative amplitude: the $Q\to 0$ Faber-Ziman criterion

`amplitude_criterion="fz"` replaces the low-$r$ density limit with the $Q\to 0$ value of $S$.
With the level-anchored model $S_\mathrm{corr} = a(S_\mathrm{meas} - L) + 1$ and
$S_\mathrm{corr}(0) = 1 - \langle b^2\rangle/\langle b\rangle^2$:

$$\boxed{\;a_\mathrm{fz} = \frac{s_0^{\mathrm{target}} - 1}{S_\mathrm{meas}(0) - L},\qquad
s_0^{\mathrm{target}} = 1 - \frac{\langle b^2\rangle}{\langle b\rangle^2}\;}$$

$S_\mathrm{meas}(0)$ is a **robust linear extrapolation** of the data head: take
$Q \le Q_0 + \texttt{fit\_width}$ with `fit_width = 1.0` Å⁻¹ (needs $\ge 8$ points, else `None`),
fit $y = c_0 + c_1(Q - \bar Q)$ by **four** Huber-weighted least-squares passes (weights recomputed
after each solve; the first solve is unweighted and the final weight update is discarded), and
evaluate at $Q=0$: $S_\mathrm{meas}(0) = c_0 - c_1\bar Q$. A denominator with
$|S_\mathrm{meas}(0)-L| < 10^{-9}$ returns `None`.

Properties that matter:

- **Closed form, no loop.** The criterion never touches the Fourier filter, so
  `_autoscale_pass()` returns immediately with `iterations = 0` (asserted by
  `test_fz_amplitude_mode_round_trip`).
- **Requires** `b_sq_avg` **and** `c1_mode="sweep"` — both enforced in `ScalingConfig.__post_init__`
  (JS `makeConfig`). Without a statistically flat level there is nothing to anchor, and the run
  raises.
- **The *amplitude formula* is independent of $\rho_0$ and of the low-$r$ window — the *run* is
  not.** An FZ-mode pass still calls `_fit_windows()` **before** the FZ branch, so a too-narrow or
  empty low-$r$ window still raises "fit windows contain fewer than 2 points"; `autoscale()` still
  runs `detect_first_peak_onset()` and can execute a second refinement pass in FZ mode; and every
  reported real-space quantity (`low_r_rms`, `gk`, `d_r`, `density_limit_satisfied`,
  `d_r_low_r_slope_theory`) still uses $\rho_0$ and the window.

**Why the two criteria are not degenerate in the same way.** This is the core of the engine's
epistemics:

| | density limit ($a$) | Faber-Ziman ($a_\mathrm{fz}$) |
| --- | --- | --- |
| Uses | $\rho_0$, the $r$-window, the filter, the whole $[Q_\mathrm{min},Q_\mathrm{max}]$ transform | composition ($\langle b^2\rangle/\langle b\rangle^2$), the level $L$, the first ~1 Å⁻¹ of data |
| Degenerate with | $\rho_0$ (≈1:1) | the length of the $Q\to 0$ extrapolation; a wrong composition |
| Blind to | a smooth low-$Q$ deficiency (absorbed into a biased $a$ with clean residuals) | the shape of the data between $Q_0+1$ Å⁻¹ and the level window — it uses the head extrapolation and the high-$Q$ level $L$, nothing in between |
| Fails on | negative-$b$/near-null-matrix data with $Q_\mathrm{min}\gtrsim 0.8$ Å⁻¹ | Bragg-contaminated or absent low-$Q$ head |

(The "blind to" row matters: $a_\mathrm{fz}$ depends on $L$ through its denominator, and $L$ is
measured on a window at the *top* of the $Q$ range — $Q\in[23.91, 29.40]$ Å⁻¹ on the parity fixture.
It is not a low-$Q$-only quantity.)

They share little but the data, so their **ratio is informative**: agreement is evidence the
absolute scale is real, disagreement quantifies what the data cannot decide — and, because only one
of them moves with $\rho_0$, the disagreement is also a $\rho_0$ measurement (Step 10). On the
Mn₃Sn runs, where $\langle b^2\rangle/\langle b\rangle^2 = 13.06$ and $S(0) = -12.06$, the density
limit is degenerate on every run (`density_limit_satisfied = False`) and the historical hand
scalings disagree with each other by 5× (×2.5, ×2.05, ×10 for the same material); the FZ criterion
lands consistently at $a \in (5, 25)$ — `test_fz_amplitude_uses_the_composition`.

**Code:** `amplitude_from_fz_limit()`, the `if config.amplitude_criterion == "fz":` branch of
`_autoscale_pass()`; JS `amplitudeFromFzLimit()`.

> Discrepancy to know about: `amplitude_from_fz_limit()` computes its own
> $s_0^{\mathrm{target}} = 1 - \langle b^2\rangle/\langle b\rangle^2$ and **ignores an explicitly
> pinned** `ScalingConfig.s0_target`. So if you override `s0_target`, the omitted-low-$Q$ correction
> (Step 5b) uses your value while the FZ amplitude and the concordance diagnostic use the
> composition-derived one. Same in both engines.

---

### Step 10 — `estimate_rho0`: the density from criteria concordance

**The idea.** $a_\mathrm{density}(\rho_0)$ grows ≈linearly with $\rho_0$ (Step 5b); $a_\mathrm{fz}$
does not depend on $\rho_0$ at all. So the true density is the root of

$$\mathrm{concordance}(\rho_0) \;\equiv\; \frac{a_\mathrm{fz}}{a_\mathrm{density}(\rho_0)} \;=\; 1 .$$

**Operation.** `estimate_rho0()` forces `amplitude_criterion="density"`, `c1_mode="sweep"`, seeds
$\rho \leftarrow \mathrm{clip}(\texttt{config.rho0},\ 10^{-4},\ 1.0)$ Å⁻³, and iterates (max 8
passes, each a full `autoscale()` including its two-pass $r_0$ refinement):

$$\rho \;\leftarrow\; \mathrm{clip}\big(\rho\cdot\mathrm{concordance}(\rho),\ \rho_\mathrm{min},\ \rho_\mathrm{max}\big),\qquad
\rho_\mathrm{min}=10^{-4},\ \rho_\mathrm{max}=1.0\ \text{Å}^{-3},$$

stopping when $|\mathrm{concordance}-1| \le \texttt{rtol} = 10^{-3}$. Because
$a_\mathrm{density} \propto \rho_0$ the update is Newton-like; observed 2–4 passes from seeds
spanning a 10× range.

Two consequences of "each pass is a full `autoscale()`" that the shape of the update hides:

- **The objective is not held fixed.** $r_0$ detection and the window refinement re-run at every
  density iterate, so the C2 window can *move* between passes; $\mathrm{concordance}(\rho_0)$ is
  therefore not a smooth function of $\rho_0$ alone.
- **The level sweep is recomputed on every pass** even though $L$ cannot change with $\rho_0$. A
  single estimate can therefore cost up to $8\times 2 = 16$ self-consistency loops and 16 level
  sweeps.

**Honest failure modes, all implemented as explicit stops:**

- Missing `b_sq_avg` → `ValueError` ("without a composition the density and the amplitude are
  degenerate").
- No usable $a_\mathrm{fz}$ (no flat level, degenerate head) → `ValueError`.
- $a_\mathrm{density} \le 0$ or concordance $\le 0$ → **break with `converged=False`**: the two
  criteria are irreconcilable at *any* density (typically missing low-$Q$ structure) and iterating
  would produce garbage. `test_estimate_rho0_fails_honestly` pins this on the Mn₃Sn 500 K run.
- $\rho$ pinned at a clip bound (no progress) → break, `converged=False`.

Callers refuse to adopt a non-converged estimate. The CLI raises with a message pointing at
`--rho0` / `--mass-density` / the data header and suggesting `--amplitude fz`; the worker throws the
same physics in browser wording — "Set ρ₀ explicitly (value, data header, or mass density) and
consider the Faber-Ziman Q→0 amplitude criterion for the scale" — since there is no CLI flag in the
browser path.

**Returned dict:** `rho0`, `converged`, `iterations`, `concordance`, `a_density`, `a_fz`,
`extrapolated`, `history` (rows `[rho0, a_density, a_fz, concordance]`). `rho0` is the density at
which the **last** pass ran, so on success it is the value that produced the accepted concordance.
`extrapolated` is simply `config.qmin > 1.0` Å⁻¹ — a flag meaning the $Q\to 0$ extrapolation is
longer than the data it rests on, so the estimate is *a starting point, not a measurement*.

**Validation:** synthetic truth $\rho_0 = 0.05$ Å⁻³ recovered from seeds 0.02 and 0.2 to within 5 %
(`test_estimate_rho0_recovers_density`); FeCoSn x-ray recovered as 0.0600 Å⁻³ against the expert's
hand value 0.057329 Å⁻³ (4.7 %, tolerance 10 %) from a seed of 0.03
(`test_estimate_rho0_near_hand_value`).

**The ~10⁻⁴ parity caveat.** The loop stops the *first* time $|\mathrm{concordance}-1| \le 10^{-3}$.
Two engines whose single-pass arithmetic differs at round-off can therefore accept **different
iterates** inside that band, and since concordance varies roughly as $1/\rho$, the accepted $\rho_0$
carries an ambiguity of order `rtol`. Cross-engine agreement on the iterated density is bounded by
the stopping rule — the vitest parity test allows **10⁻⁴ relative cross-engine on $\rho_0$**, and
*separately* checks that the JS concordance lands within **1.5·10⁻³ of 1** (the physical target);
the golden concordance itself is never compared — **not** by transform precision, unlike a single
`autoscale()` call where the tolerance is 10⁻⁶ on $(a,b)$.

**Code:** `rmc_toolkits/scaling.py` → `estimate_rho0()`; JS `estimateRho0()`.
`web_app/frontend/src/workers/autoScaleWorker.js` runs it first and adopts the result when $\rho_0$
resolves from **no source at all** — not merely when the form field is blank, but when the field,
the `stog.inp`, the `NUMBER_DENSITY ::` data header and mass-density-plus-composition all fail
(Step 1b) — *and* a composition supplies $\langle b^2\rangle$; without $\langle b^2\rangle$ the page
refuses to run. The page also exposes it directly as the **Estimate ρ₀** button
(`kind: 'estimateRho0'`), which writes the result back into the $\rho_0$ field at 5 significant
digits.

---

### Step 11 — Final pipeline and outputs

With $(a, b)$ fixed, `scale_pipeline()` runs the full chain once more on the cropped data:

1. $S_\mathrm{scaled} = a S_\mathrm{meas} + b$.
2. `fourier_filter()` → `sq_filtered`, `sq_ft` (the classic `ft.dat`), `g_filtered`.
3. Keen conversions: $G_K = \langle b\rangle^2 (g_\mathrm{filtered} - 1)$,
   $D(r) = 4\pi\rho_0 r\,G_K(r)$, $F_K(Q) = \langle b\rangle^2 (S_\mathrm{filtered} - 1)$.
4. **Optional low-$r$ enforcement** — see below.
5. **A fifth series the engine does not produce.** The *unfiltered* $g(r)-1$ (the classic
   `scale.gr`) is not part of `ScalingResult`. Both the CLI (`scaling_cli.py` → `_write_outputs()`)
   and the browser worker (`autoScaleWorker.js`) recompute it separately: forward-transform
   $F(Q) = Q(S_\mathrm{scaled}-1)$ of the **scaled but unfiltered** data with the same
   `lorch` / `low_q_correction` / `s0_target` settings, then divide by $4\pi\rho_0 r$ (no $+1$),
   giving `gm1Unfiltered`. It is what the exported `<stem>.gr` contains, and it is plotted nowhere —
   it exists only in the file family.

#### Low-$r$ enforcement, precisely

`enforce_low_r()` hard-replaces $G_K(r) = -\langle b\rangle^2$ for every $r \le$ cutoff — exactly
what the Fortran stog does to its RMC outputs (verified to 16 digits on the validation run). It is
the special case of the general Fortran semantics in `first_peak_zero()`: zero $g$ where
$r\le\texttt{cutoff}$ **and** $r$ lies outside $[\texttt{peak\_rmin}, \texttt{peak\_rmax}]$, i.e.
`(r <= cutoff) & ((r >= peak_rmax) | (r <= peak_rmin))`. `scale_pipeline`'s own `enforce_cutoff`
path (Python only) uses `enforce_low_r`; the CLI and the browser worker use `first_peak_zero`.

Which triple they use is **not** simply "the detected $r_0$" — an explicit cutoff always wins:

| Situation | `(cutoff, peak_rmin, peak_rmax)` |
| --- | --- |
| CLI with a `stog.inp` and no `--enforce-cutoff` | the inp's own `(peak_cutoff, peak_rmin, peak_rmax)` — genuinely the general Fortran semantics |
| CLI `--enforce-cutoff C` (optionally `--peak-window lo hi`) | `(C, lo, hi)`, or `(C, C, C)` without a peak window |
| CLI `--data` mode, no explicit cutoff, enforcement not disabled | the detected $r_0$, as $(r_0, r_0, r_0)$ from the post-run detection |
| Page: cutoff field filled | `(C, C, C)` |
| Page: cutoff blank, `stog.inp` loaded | the inp's `(peak_cutoff, peak_rmin, peak_rmax)` |
| Page: cutoff blank, no inp (`enforcement: 'auto'`) | the detected $r_0$, as $(r_0, r_0, r_0)$ |

In the degenerate `(r_0, r_0, r_0)` case the predicate `(r >= peak_rmax) | (r <= peak_rmin)` is true
for **every** $r \le$ cutoff, so `first_peak_zero` **collapses exactly to `enforce_low_r`**:
$g \equiv 0$, hence $G_K \equiv -\langle b\rangle^2$, for all $r \le r_0$. The general two-sided form
only does something different when a real peak window is supplied. Note also that the page ships
**Enforce low-r checked by default** (`EMPTY_FORM.enforce = true`), and that a manual run reaches
this code with `r0Detected` recovered by the worker's own `detectFirstPeakOnset()` call (Step 2).

**Enforcement makes the published files satisfy Keen Eq. 15 by construction, not as evidence of a
good fit.** That is why both the result object and every diagnostic report the **pre-enforcement**
`low_r_rms`. Judge quality only on that number.

#### The exported file family

The page's deliverable is a zip of nine entries, written by `AutoStogPage.jsx` → `writeFiles()`
(the CLI writes the same family into an `autoscale/` directory beside the input file, under the
`stog.inp`'s declared output names or a `--out-stem`):

| File | Contents |
| --- | --- |
| `<stem>.sq` | $S_\mathrm{scaled} = aS_\mathrm{meas}+b$, unfiltered |
| `<stem>.gr` | the **unfiltered** $g(r)-1$ (`gm1Unfiltered`, item 5 above) |
| `<stem>_ft.sq` | $S_\mathrm{filtered}$ |
| `<stem>_ft.gr` | $g_\mathrm{filtered}-1 = G_K/\langle b\rangle^2$, **plus a third column** $4\pi\rho_0 r\,(g-1) = G_\mathrm{PDF}(r)$ — a *different function*, not an error bar |
| `<stem>_rmc.fq` | $F_K(Q)$, barn (RMCProfile input) |
| `<stem>_rmc.gr` | $G_K(r)$, barn — the **enforced** array when enforcement is on |
| `<stem>_rmc.dr` | $D(r)$ — likewise enforced |
| `ft.dat` | $S_\mathrm{ft}$, the classic fixed-name filter correction |
| `<stem>_provenance.json` | $a$, $b$, mode, `c1ModeEffective`, the loaded inp's hand scaling, `history`, `enforcement`, `rho0Estimate`, the config, and the diagnostics summary |

Numeric layout (`writeStogXy` / `write_stog_xy`): a point-count line right-aligned in 12 characters
(`%12d`), a title
line — `rmc-autoscale (browser): a=… b=…` with the values at 8 significant digits in the browser,
`rmc-autoscale <version>: …` in the CLI — then Fortran-style `%.16E` columns; the JS writer pads a
single-digit exponent to two digits so the two writers agree byte-for-byte on format.

**Code:** `scale_pipeline()`, `transforms.enforce_low_r()`, `transforms.first_peak_zero()`,
`scaling_cli.py` → `_resolve_enforcement()`, `_write_outputs()`; JS `scalePipeline()`,
`firstPeakZero()`, `writeStogXy()`, the enforcement block in `autoScaleWorker.js`, and
`AutoStogPage.jsx` → `resolveEnforcement()`, `writeFiles()`.

---

### Step 12 — Every diagnostic returned, and how to read it

`diagnostics_summary(result, config)` (JS: `diagnosticsSummary`) is where the numbers on the page
come from. **The two implementations read their inputs from different places:**

- **Python** reads the coefficients, $\rho_0$, the amplitude criterion *and* the $r$ window from
  `result.provenance["config"]` / `result.provenance["r_fit_window"]`, so a window-refined result is
  judged against the window and parameters it was actually fitted with.
- **JS** results carry no provenance dict. `diagnosticsSummary` takes the window from
  `result.rFitWindowUsed` (equivalent) but reads $\langle b\rangle^2$, $\rho_0$, $\langle b^2\rangle$
  and `amplitudeCriterion` from the **caller's `config` object**. In the worker these coincide,
  because the worker overwrites `config.rho0` with the adopted estimate before calling it — but a
  direct JS caller can get theory lines computed from a different $\rho_0$ than the fit used.

| Key | Definition | How to read it |
| --- | --- | --- |
| `a`, `b` | the fitted correction | In sweep mode $b$ is not free: $b = 1 - aL$ exactly. |
| `converged`, `iterations` | fixed-point loop status | `converged = False` means the loop hit `max_iter = 50`; treat $(a,b)$ as unfinished. FZ mode and manual mode always report `iterations = 0`. |
| `c1_tail_mean` | mean of the **filtered** $S(Q)$ over the C1 tail window $Q \ge q_\mathrm{max} - 0.15\,(q_\mathrm{max}-q_\mathrm{min})$, from the **configured** `qmin`/`qmax` | Should sit on 1. Tests accept 0.02–0.05 depending on dataset. Because the window is configured, not data-derived, a `qmax` above the end of your data shifts it (and can empty it — "fit windows contain fewer than 2 points"). |
| `low_r_rms_pre_enforcement` | $\sqrt{\langle g_\mathrm{filtered}^2\rangle}$ over the low-$r$ window, **before** any enforcement | The honest fit-quality number. Compare it against a hand-scaled run — `test_autoscale_beats_hand_tuning` requires the auto fit to be ≤ the expert's (×1.001). |
| `g_window_mean` | signed mean of $g_\mathrm{filtered}$ over the same window (target 0) | A large positive mean is the signature of missing structure below $Q_\mathrm{min}$. |
| `r_fit_window` | the window actually used, Å | Check it moved to where you expect after $r_0$ detection. |
| `gk_low_r_theory` | $-\langle b\rangle^2$ | The flat level the plotted $G_K$ should reach. |
| `d_r_low_r_slope_theory` | $-4\pi\rho_0\langle b\rangle^2$ | The straight line $D(r)$ should follow below $r_0$. |
| `density_limit_satisfied` | $\lvert$`g_window_mean`$\rvert < 0.1$ | **ONE-SIDED.** False *proves* no affine $(a,b)$ can satisfy the density limit — the absolute scale is not recoverable from self-consistency on this data. True only means the fit reached its target; a smooth low-$Q$ deficiency is generically absorbed into a biased scale with all residuals clean. True does **not** certify the absolute scale. |
| `level`, `level_uncertainty`, `level_window`, `asymptote_found` | the Step-3 sweep result | `asymptote_found = False` ⇒ the fit silently ran in joint 2-dof mode, `level_uncertainty` is `NaN`, **and `level_window` is the fabricated last-3 Å⁻¹ span, not a searched window**. Otherwise `level_uncertainty` is a spread over overlapping admissible windows. |
| `r0_detected`, `window_refined` | Step-8 outputs | If `window_refined` is False, the detected $r_0$ was reported but did not change the fit (because you pinned `r0`/`r_fit_max`, or the change was < 0.05 Å). |
| `a_fz` | the independent $Q\to 0$ amplitude | Present only when $\langle b^2\rangle$ is available **and** the sweep found a flat level (`a_fz` is computed inside `if level is not None`). Otherwise it — and the concordance row — are absent. In FZ mode it *is* `a`. |
| `amplitude_concordance` | $a_\mathrm{fz}/a$ — **omitted in FZ mode** (it would be 1 by construction) | The absolute-scale trust metric. FeCoSn agrees to 4–6 %. |
| `amplitudes_concordant` | $\lvert a_\mathrm{fz}/a - 1\rvert < 0.1$ | Discord ⇒ suspect $\rho_0$ (moves only $a$), or missing low-$Q$ (moves them apart). The page's chip says "check ρ₀ / low-Q, or use the Faber-Ziman Q→0 amplitude". |
| `fk_qmin`, `fk_q0_theory` | the **filtered** $F_K$ at the first cropped $Q$ point (`result.fk[0]`, computed from `sq_filtered`) vs $-\langle b^2\rangle$ | The C3 diagnostic. Data rarely reach $Q\approx 0$, so this is a *trend* check, not an equality. |

Additional fields on the result object / provenance: `history` (per-iteration
$(a, b, \mathrm{low-}r\ \mathrm{rms})$), `sweep` (full `LevelSweepResult` including `slope`,
`slope_sigma`,
`n_admissible`), `n_despiked` (with the double-count caveat of Step 1d), `q_tail_window`,
`n_q_points`, `mode` (`"auto"`/`"manual"`), `c1_mode_effective` (`"sweep"`, `"joint"` or, in JS
manual runs, `"manual"` — the *actual* architecture used), and a **config echo**.

> The config echo is a **fixed 26-key subset** of `ScalingConfig`'s 28 fields: `qmin`, `qmax`,
> `rho0`, `b_avg_sq`, `b_sq_avg`, `r_cutoff`, `r0`, `fit_offset`, `q_tail_frac`, `rmax`, `nr`,
> `lorch`, `low_q_correction`, `c2_weight`, `c2_bins`, `robust`, `c1_mode`, `amplitude_criterion`,
> `s0_target`, `c1_slope_nuisance`, `despike`, `despike_window`, `despike_nsigma`, `max_iter`,
> `tol`, `enforce_cutoff`. **`r_fit_min` and `r_fit_max` are not echoed** — exactly the two fields
> that define the C2 window. Reading a provenance JSON, do not conclude the window edges were
> defaults; the window actually used is recorded separately as `r_fit_window`.

---

### Python ↔ JavaScript parity

Both engines do float64 arithmetic with the same trapezoid sine-transform discretization on the same
grids. `tests/generate_autoscale_fixture.py` produces golden numbers from the Python engine that
`web_app/frontend/src/__tests__/autoScale.test.js` checks the JS engine against.

| Quantity | Asserted cross-engine agreement |
| --- | --- |
| `level` | 10⁻⁹ relative |
| `qLo`, `qHi` | `toBeCloseTo(…, 9)` — **absolute** (≈5·10⁻¹⁰), not relative |
| `nAdmissible` | **exact equality** (797 on the fixture) |
| `levelUncertainty` | 10⁻⁶ relative |
| $a$, $b$ (sweep + density), $a$, $b$ (FZ mode), $a$ (detection pass) | 10⁻⁶ relative; iteration counts equal exactly |
| `lowRRms` | 10⁻⁵–10⁻⁶ relative; `c1TailMean` 10⁻⁸ |
| sampled `gk`, `sqFiltered` | 9 decimal places |
| `estimateRho0.rho0` | 10⁻⁴ relative — bounded by the `rtol` stopping rule, not by transform precision. (The companion assertion `abs(concordance − 1) < 1.5·10⁻³` is a check against the *physical target 1*, **not** a comparison with the Python golden concordance, which the test never asserts on.) |

Genuine implementation differences:

1. **Linear algebra.** Python uses `np.linalg.lstsq` (SVD); JS forms the normal equations
   $A^\top A x = A^\top y$ and solves by Gaussian elimination with partial pivoting (throwing on a
   pivot below 10⁻³⁰⁰). For the ≤3-column, well-conditioned systems here the answers agree to
   round-off, but the JS path is the numerically weaker one if a design ever becomes ill-conditioned.
2. **Summation order.** numpy's pairwise summation vs JS sequential accumulation — the source of the
   tolerances above.
3. **Python-only features.** `c2_bins` (binned C2 levels), `c1_slope_nuisance` (tail-drift column)
   and `enforce_cutoff` inside `scale_pipeline` have **no JS equivalent**. The browser enforces the
   low-$r$ level in `autoScaleWorker.js` via `firstPeakZero()` instead. The Auto StoG page does not
   expose `c2_bins`, `c1_slope_nuisance`, `c2_weight`, `q_tail_frac`, `max_iter` or `tol` in either
   runtime.
4. **A real divergence in the loop's filter call.** Python's `_pipeline()` passes
   `s0_target=config.effective_s0_target` to `fourier_filter`; the JS `autoscalePass()` inner-loop
   call **omits `s0Target`, so it defaults to 0**, while the JS final `scalePipeline()` passes it
   correctly. The two engines therefore compute slightly different $\delta(Q)$ whenever the
   composition-aware $S(0)$ target is non-zero. Measured on the parity fixture ($S(0) = -0.74$): the
   estimated $\rho_0$ shifts by 2.8·10⁻⁵ relative — inside the test tolerance. Bounded on a
   Mn₃Sn-like synthetic ($S(0) = -12.06$): $a$ shifts by ≈1.8·10⁻³ relative. Small, but it is a code
   difference, not float noise.
5. **`np.linspace(...).astype(int)` emulation.** JS reproduces the sweep's edge grid with
   `Math.trunc(k*(n-1)/(nGrid-1))`. This matches numpy's truncation on the tested data; a rare
   floating-point tie could in principle shift one edge index by 1.
6. **Despike scope — the engines fit and report on different datasets when `despike=True`.**
   Python's `_autoscale_pass()` hands `scale_pipeline()` the *already cropped and despiked* arrays,
   so the rolling-median filter runs **twice** and `provenance["n_despiked"]` counts only the second
   application (Step 1d: 107 then 75 on the repo's fixture, 75 reported). The JS `autoscalePass()`
   calls `scalePipeline(qIn, sqIn, …)` with the **raw** input arrays, so JS despikes once and
   `nDespiked` is the true count. Published `sq_filtered`/`gk`/`low_r_rms`/`c1_tail_mean` and the
   despike counts are therefore **not comparable between engines** with despiking on, and the "the
   number actually removed" reading of `n_despiked` holds only for JS. The parity fixture does not
   exercise this — `despike` is off in it.
7. **Diagnostics inputs.** Python's `diagnostics_summary` reads `result.provenance["config"]`; JS's
   reads the caller's `config` for $\langle b\rangle^2$, $\rho_0$, $\langle b^2\rangle$ and
   `amplitudeCriterion` (window from `result.rFitWindowUsed`). See Step 12.
8. **Key casing.** JS `estimateRho0` returns `aDensity`/`aFz` where Python returns
   `a_density`/`a_fz`, while the `summary` object's keys stay snake_case to match Python. Reading an
   exported browser provenance JSON against this document, look for `aDensity`, not `a_density`.
9. **Data-file parsing.** The JS `readStogXy` rewrites Fortran `D` exponents to `e`; Python's
   `read_stog_xy` uses bare `float()`, which rejects them and skips the row.

---

### Parameters and defaults

| Parameter | Default | Units | Where |
| --- | --- | --- | --- |
| `qmin`, `qmax` | *required* | Å⁻¹ | fit + transform window (also define the C1 tail window) |
| `rho0` | *required*, or resolved / self-consistently estimated (Step 1b, Step 10) | Å⁻³ | C2 density line |
| `b_avg_sq` $=\langle b\rangle^2$ | *required*, or from composition | barn | Keen normalization |
| `b_sq_avg` $=\langle b^2\rangle$ | `None`, or from composition | barn | enables FZ amplitude, C3 diagnostic, composition-aware $S(0)$, `estimate_rho0` |
| `r_cutoff` | 1.0 | Å | Fourier-filter cutoff |
| `r0` | `None` → detected | Å | closest approach |
| `r_fit_min` | `None` → `r_cutoff + 0.2` | Å | C2 window lower edge (**not** echoed in provenance) |
| `r_fit_max` | `None` → `r0 - 0.25`, else `r_fit_min + 1.0` | Å | C2 window upper edge (**not** echoed in provenance) |
| `q_tail_frac` | 0.15 | — | C1 window = top 15 % of the *configured* $[Q_\mathrm{min},Q_\mathrm{max}]$ |
| `rmax`, `nr` | 50.0, 5000 | Å, count | $r$ grid ($\Delta r = 0.01$ Å), $r=0$ excluded; also sets the C2 row count |
| `lorch` | `False` | — | Lorch window $M(Q)=\sin(x)/x$, $x=\pi Q/Q_{N-1}$, applied to $F(Q)$ before the forward transform only — built from the **last cropped data point**, not `qmax` |
| `low_q_correction` | `True` | — | analytic $[0,Q_\mathrm{min}]$ correction |
| `s0_target` | `None` → $1-\langle b^2\rangle/\langle b\rangle^2$ if available, else 0 | — | target of the low-$Q$ extrapolation (ignored by `amplitude_from_fz_limit`) |
| `c2_weight` | 1.0 | — | $\sqrt{\cdot}$ applied to C2 rows (per row — block balance is set by row counts) |
| `c2_bins` | 0 (off) | count | Python only |
| `robust` | `True` | — | Huber IRLS, 3 passes, $c = 1.345$ |
| `c1_mode` | `"sweep"` | — | `"joint"` = original 2-dof fit |
| `amplitude_criterion` | `"density"` | — | `"fz"` = $Q\to0$ Faber-Ziman; forced to `"density"` on manual runs |
| `fit_offset` | `True` | — | joint mode only |
| `c1_slope_nuisance` | `False` | — | experimental, Python only |
| `despike`, `despike_window`, `despike_nsigma` | `False`, 7 (**must be odd, unvalidated**), 6.0 | —, points, MAD-σ | glitch removal |
| `max_iter`, `tol` | 50, 1e-6 | count, relative | self-consistency loop |
| `enforce_cutoff` | `None` | Å | classic low-$r$ enforcement (page default: **on**, cutoff auto) |
| `use_sigma` (page only) | `true` | — | σ column used automatically when present; weights C1 rows only |
| `level_sweep`: `min_width`, `n_grid`, `slope_nsigma` | 3.0, 80, 2.0 | Å⁻¹, count, σ | **not configurable** |
| `detect_first_peak_onset`: `search_min`, `search_max`, `fraction`, `floor` | `r_cutoff+0.3`, 6.0, 0.35, 0.5 | Å, Å, —, $\lvert g\rvert$ | **not configurable** |
| `amplitude_from_fz_limit`: `fit_width` | 1.0 | Å⁻¹ | head extrapolation span, ≥8 points |
| `estimate_rho0`: `rtol`, `max_iter`, `rho_min`, `rho_max` | 1e-3, 8, 1e-4, 1.0 | —, count, Å⁻³, Å⁻³ | fixed-point root-find |
| diagnostic thresholds | $\lvert$`g_window_mean`$\rvert<0.1$; $\lvert a_\mathrm{fz}/a - 1\rvert<0.1$; UI coefficient-shadowing warning at 2 % | — | one-sided / concordance verdicts |

**Hard minimums.** ≥16 points after cropping (checked *before* despiking, never re-checked); ≥32
finite points for the level sweep; ≥24 points and ≥3.0 Å⁻¹ width per candidate window; ≥2 points in
each fit window; ≥3 $r$-grid points inside the $r_0$ search range (else `detect_first_peak_onset`
returns `None`); ≥8 points in the FZ head; ≥4 C2 rows for C2 IRLS re-weighting.

**Configuration validation** (all raise *before* any math, in `ScalingConfig.__post_init__` /
`makeConfig`): `c1_mode ∈ {sweep, joint}`; `amplitude_criterion ∈ {density, fz}`; `fz` requires
`b_sq_avg` **and** `c1_mode="sweep"`; $\rho_0$ finite and $>0$; $\langle b\rangle^2$ finite and
$>0$; `qmax > qmin`; `nr` a positive integer; `rmax` finite and $>0$. `r_fit_window` raises
"empty low-r fit window" whenever the upper edge $\le$ the lower edge — the JS `makeConfig`
evaluates this eagerly at construction, the Python property only when first accessed (the CLI
touches `config.r_fit_window` deliberately so the error renders as a CLI error).
**Not validated in either engine:** `q_tail_frac`, `c2_weight`, `tol`, `max_iter`, `despike_window`.

---

### Caveats / what this is not

- **The absolute scale is not certifiable from inside the data.** The engine minimizes unphysical
  low-$r$ content and matches the high-$Q$ asymptote; a *smooth* low-$Q$ deficiency (data missing
  structure below $Q_\mathrm{min}$) is generically absorbed into a biased $a$ with **all residual
  diagnostics clean**. `density_limit_satisfied` is explicitly one-sided. Validate externally
  (known composition/density, or a fixed-$a$ `scale_pipeline` run) whenever $Q_\mathrm{min}$ is not
  small.
- **A wrong $\rho_0$ yields a confidently wrong scale.** The C2 target line is
  $-4\pi\rho_0\langle b\rangle^2 r$, so $a$ tracks $\rho_0$ roughly 1:1. The concordance
  diagnostic is the only internal
  alarm; heed it.
- **Several physics parameters can be filled in silently.** $\rho_0$ from a `NUMBER_DENSITY ::`
  header or a `stog.inp`, $r_0$ from `MINIMUM_DISTANCES ::` or a peak window, the $Q$ window from
  the data's extent at 4 significant digits, and — in the last resort — a 0.05 Å⁻³ seed replaced by
  the self-consistent estimate. Read Step 1b and the provenance JSON before believing a number you
  did not type.
- **Neutron / constant-$b$ is the design case.** X-ray data carry $Q$-dependent $f(Q)$ and Compton
  contributions that a constant $(a,b)$ cannot model. The engine is *used* on x-ray data (the
  FeCoSn validation runs, with $\langle b\rangle^2 = 1$ and
  $\langle b^2\rangle = \langle Z^2\rangle/\langle Z\rangle^2$) and the `c1_slope_nuisance`
  column exists to absorb tail drift,
  but the Sears table is neutron and the plan document still labels general x-ray use as out of
  scope for v1.
- **The input grid is trusted.** No sorting, deduplication, interpolation or rebinning is ever
  performed; a non-monotonic $Q$ column gives undefined results and makes the two engines diverge.
- **The level sweep operates inside the configured $[Q_\mathrm{min}, Q_\mathrm{max}]$**, so it
  cannot warn you about rolloff you already cropped out, and it cannot yet inform the $Q_\mathrm{max}$
  choice itself (flagged as the next refinement in
  [STOG_SCALING_PLAN §5b](../STOG_SCALING_PLAN.md)).
- **`level_uncertainty` is a window-choice spread, not a standard error.** Admissible windows
  overlap heavily. On a failed sweep `level_window` is fabricated, not searched.
- **IRLS is capped at 3 passes with no convergence test**, and re-weights rows by $w$ rather than
  $\sqrt{w}$ (effective objective $\sum w_i^2 e_i^2$).
- **`nr`/`rmax` are not just a display grid.** They set the number of C2 rows and therefore the
  implicit C1:C2 balance; changing the $r$ grid changes the fitted $a$.
- **The RMCProfile-ready files satisfy Keen's limits by construction** when enforcement is on. Only
  `low_r_rms_pre_enforcement` measures the fit.
- **Despiking is a blunt instrument on crystalline data** (12 % of points flagged on the 59438
  benchmark) — use it for detector glitches only, and read `n_despiked` knowing the Python engine
  applies the filter twice and reports only the second pass's count.
- **Some of what you see on screen is display-only** and does not affect the exported files: the
  real-space theory guide lines are drawn only out to $\min(8\ \text{Å}, r_\mathrm{max})$
  (`RMAX_DISPLAY = 8`); the $G_K$ plot opens with a fixed $y$-domain of
  $[2.1\,(-\langle b\rangle^2),\ 3.2\,\langle b\rangle^2]$; the $S(0)$ Faber-Ziman target guide spans
  only $[Q_0, \min(Q_0+1.5, Q_{N-1})]$; and when enforcement is on the **enforced** curve is the
  primary $G_K$/$D(r)$ series while the honest pre-enforcement fit is one legend click away.
- **The engine is I/O-free.** File reading/writing lives in
  [`rmc_toolkits/scaling_cli.py`](../../rmc_toolkits/scaling_cli.py) (Python) and the page/worker
  (browser); nothing in `scaling.py` touches disk.


## Auto StoG — page workflow, outputs, and the written file family

### What this page shows

**Auto StoG** is the only *pre*-processing tab in the app: it takes a measured, rebinned
total-scattering structure factor $S_\mathrm{meas}(Q)$ and produces the classic
stog / RMCProfile-ready file family, choosing the scale and offset of the affine correction

$$S_\mathrm{corr}(Q) \;=\; a\,S_\mathrm{meas}(Q) + b$$

from physics rather than from the operator's eye. The page renders three curves (the scaled
and filtered $S(Q)$, the Keen $G_K(r)$, and $D(r)$), a row of numeric verdict cards, and an
export button that writes a nine-file zip.

> **Reachability as of this commit.** The Auto StoG tab is behind a compile-time flag and is
> **not present in the shipped build**. [`App.jsx`](../../web_app/frontend/src/App.jsx) defines
> `const SHOW_AUTO_STOG = false;` (line 28, comment: “Auto StoG is still under development;
> flip to true to expose the tab again”), and **both** the nav button (`App.jsx:295`,
> `{SHOW_AUTO_STOG && (…)}`) and the page mount (`App.jsx:439`,
> `{SHOW_AUTO_STOG && visitedPages.autostog && (…)}`) are gated on it. The page and its engine
> are complete and parity-tested; flipping the constant to `true` exposes the tab. Everything
> below describes the code as written, not a currently reachable UI.

Two properties matter for reading everything below:

- **It is independent of the run folder** used by the Dashboard / Atomic Density / PCA
  Ellipsoid tabs. The S(Q) is uploaded *into the page* (file picker or drag-and-drop).
  Source: [`AutoStogPage.jsx`](../../web_app/frontend/src/components/AutoStogPage.jsx) header
  comment and `ingestFiles()`. When the flag is on, the tab button is always present and the
  page mounts on first visit (`visitedPages.autostog`) with no props and no run-folder
  precondition — it receives nothing from the rest of the app.
- **It never talks to the backend.** `AutoStogPage.jsx` imports no API helper
  (`../api.js` is absent from its import list) and issues no `fetch`. The Flask endpoints
  `/api/scaling/preview` and `/api/scaling/run` in
  [`web_app/backend/app.py`](../../web_app/backend/app.py) exist and are tested
  (`tests/test_backend_api.py::ScalingApiTests`) but are reachable only from the API/CLI,
  not from this page.

#### Notation and units

| Symbol | Meaning | Units |
| --- | --- | --- |
| $Q$ | scattering-vector magnitude | Å⁻¹ |
| $r$ | interatomic separation | Å |
| $S(Q)$ | Faber-Ziman normalized total structure factor (Keen Eq. 19), $\to 1$ at high $Q$ | dimensionless |
| $a,\ b$ | multiplicative scale / additive offset of the correction | dimensionless |
| $L$ | measured high-$Q$ level of $S_\mathrm{meas}$ (level sweep) | dimensionless |
| $\rho_0$ | atomic number density | Å⁻³ |
| $\langle b\rangle^2 = (\sum_i c_i b_i)^2$ | classic stog "Faber-Ziman coefficient" | barns |
| $\langle b^2\rangle = \sum_i c_i b_i^2$ | Keen Eq. 14 $Q\to0$ coefficient (a *different* number) | barns |
| $g(r)$ | pair distribution function, $\equiv 0$ below closest approach, $\to 1$ | dimensionless |
| $G_\mathrm{PDF}(r) = 4\pi\rho_0 r\,[g(r)-1]$ | PDFFIT-style function | Å⁻² |
| $G_K(r) = \langle b\rangle^2[g(r)-1]$ | Keen Eq. 10/16 | barns |
| $D(r) = 4\pi\rho_0 r\,G_K(r)$ | Keen Eq. 29 | barns·Å⁻² |
| $F_K(Q) = \langle b\rangle^2[S(Q)-1]$ | Keen Eq. 9 | barns |
| $F(Q) = Q[S(Q)-1]$ | PDF-community $F(Q)$, transform argument | Å⁻¹ |
| $r_0$ | closest interatomic approach | Å |

1 barn = 100 fm². Scattering lengths are stored in fm and converted to barns by
dividing $\langle b\rangle^2$, $\langle b^2\rangle$ (fm²) by 100 —
[`autoScale.js`](../../web_app/frontend/src/workers/autoScale.js) → `faberZiman()`,
[`scattering.py`](../../rmc_toolkits/scattering.py) → `faber_ziman()`.

---

### Step 0 — Where each piece of the computation runs

| Stage | Flask mode | Static mode (`VITE_STATIC_MODE=true`) |
| --- | --- | --- |
| File reading / parsing (`.inp`, xy data, `::` header) | browser, `autoScale.js` ports | identical |
| Faber-Ziman coefficients from the composition | browser, `autoScale.js` → `faberZiman()` | identical |
| Level sweep, affine fit, Fourier-filter loop, $\rho_0$ estimate | browser Web Worker | identical |
| Plot rendering | browser SVG (`InteractivePlot.jsx`) | identical |
| Output files + zip | browser (`writeStogXy` + `buildZip`) | identical |
| Anything on the server | **nothing** | **nothing (no server exists)** |

The Python side of the same math is reachable only outside this page:
`rmc-autoscale` (console script → [`scaling_cli.py`](../../rmc_toolkits/scaling_cli.py) `main()`),
and `POST /api/scaling/preview` / `/api/scaling/run`, which share the CLI's writer
(`app.py` imports `_write_outputs` and `_resolve_targets` from `scaling_cli`, so the **eight
data files** are byte-identical between API and CLI). The ninth output is not: the API builds
its own `provenance_payload` with `tool: "rmc-autoscale (web API)"`, no `rmc_toolkits_version`
and no `argv`, a `source` key holding `str(inp_path or data_path)`, and a truncated
`stog_inp_reference` of only `{a, b}` (no `yscale`/`yoffset`). The API also exposes **no $\rho_0$
self-consistency** (`--estimate-rho0` has no HTTP counterpart — `estimate_rho0` is never
imported by `app.py`) and memoizes engine results in an `@lru_cache(maxsize=8)` keyed on
`(data path, mtime, config, mode, a, b, use_sigma)` (`app.py` → `_cached_scaling()`).
The browser engine
[`workers/autoScale.js`](../../web_app/frontend/src/workers/autoScale.js) is a straight port of
[`scaling.py`](../../rmc_toolkits/scaling.py) + [`transforms.py`](../../rmc_toolkits/transforms.py)
+ the stog parsers + `scattering.py`; agreement is pinned by
[`autoScale.test.js`](../../web_app/frontend/src/__tests__/autoScale.test.js) against Python
goldens (`tests/generate_autoscale_fixture.py`) — see “Cross-engine agreement” below for the
measured tolerances.

---

### Step 1 — Page-local upload and source classification

**Inputs:** one or more files from the `<input type="file" multiple accept=".sq,.fq,.dat,.inp">`
or from a drop on the dropzone.

**Operation** (`ingestFiles()`): each `File` is accepted only if its *name* matches one of two
predicates (top of `AutoStogPage.jsx`):

```js
isInpCandidate  = name.endsWith('.inp') || name === 'stog_input.dat'
isDataCandidate = /\.(sq|fq|dat)$/i.test(name)
                  && !/^(scale[._]|ft\.dat)/i.test(name)
                  && !/_rmc\.(fq|gr|dr)$/i.test(name)
```

The two negative patterns exist so that dropping a whole classic-stog run folder does not
offer that run's *outputs* (`scale.fq`, `scale_ft.sq`, `ft.dat`, `*_rmc.gr`, …) as scaling
*inputs*. Accepted files are read with `File.text()` (whole file into a JS string); anything
that throws (e.g. a dropped directory entry) is collected in an `unreadable` list. That list
is **only surfaced when the batch yielded no readable file at all** (`if (!accepted.length)`
→ the “Could not read …  — drop individual files, not a folder.” error); if at least one file
was read successfully the failures are dropped silently.

Files are then merged into the page's `sources` list: existing entries whose name matches a
newly accepted one are filtered out and the accepted list is appended, so **re-uploading a
file with the same name replaces its stored text**. The merged list is sorted `.inp`
candidates first, then alphabetically (`leftInp - rightInp || left.name.localeCompare(right.name)`).
Auto-selection is **per upload**: the first `.inp` *in the current batch* — else that batch's
first accepted file — is selected (`accepted.find(isInpCandidate) || accepted[0]`), not the
first `.inp` in the merged list. A later data-only upload therefore replaces an already-selected
`.inp` as the active source.

Selecting *any* source (auto or from the dropdown) resets the page state: `selectSource()`
clears `preview`, `exportResult`, `error` and `rho0Info` and nulls `dataRef.current`, so the
plots and readout cards disappear until **Auto-scale** is pressed again.

**Outputs:** `sources: [{name, text}]`, a selected source name.

**Code:** `AutoStogPage.jsx` → `isInpCandidate`, `isDataCandidate`, `ingestFiles()`, `onDrop()`,
`selectSource()`.

> Note: `isInpCandidate` is case-sensitive (`.INP` is rejected by both predicates), and
> `stog_input.dat` matches *both* predicates — `selectSource()` tests `isInpCandidate` first,
> so it is always treated as a stog input.

---

### Step 2 — Parsing the selected source and prefilling the form

**Inputs:** the selected file's text, plus the other uploaded files (to resolve a `stog.inp`'s
declared data file).

**Operation** (`selectSource()`):

1. If the selection is an `.inp`, `readStogInp(text)` parses the classic **22-line**
   single-dataset layout — both engines require `>= 22` non-empty lines and raise otherwise,
   and read indices `[0]`–`[21]` (line indices after stripping blank lines):
   `[0]` file count (must be 1), `[1]` data file name, `[2]` `qmin qmax`,
   `[3]` `yoffset yscale`, `[4]` Q offset (must be 0), `[5]`/`[6]` output S(Q)/g(r) names,
   `[7]` `rmax`, `[8]` `nr`, `[9]` Lorch flag, `[10]` $\rho_0$, `[11]` second y-offset
   (must be 0), `[12]` "try again" flag (must be N), `[13]` filter flag (must be Y),
   `[14]` filter r-cutoff, `[15]`/`[16]` filtered output names, `[17]` $\langle b\rangle^2$,
   `[18]`–`[20]` RMC output names, `[21]` `peak_cutoff peak_rmin peak_rmax`.
   The Fortran scaling convention $S_\mathrm{scaled} = S_\mathrm{raw}/\texttt{yscale} + \texttt{yoffset}$
   is converted to the multiply convention used everywhere else:
   $$a = 1/\texttt{yscale},\qquad b = \texttt{yoffset}.$$
   Unsupported variants raise instead of mis-parsing. The declared data file must also have
   been uploaded, otherwise the page errors with
   `stog input references '<file>' — upload that file too`.
   *Code:* `autoScale.js` → `readStogInp()`; Python twin
   [`parsers.py`](../../rmc_toolkits/parsers.py) → `read_stog_inp()` / `StogInput`.
   (A stale comment in `tests/test_scaling_cli.py` calls its `INP_TEMPLATE` “23-line”; the
   template itself is 22 lines — trust the parsers.)
2. `readStogXy(dataText)` tokenizes every line; a line is a data row only if it has ≥2 tokens
   and **every** token matches the numeric regex (decimal/exponent, `nan`, `inf`, Fortran
   `D` exponents). Rows are grouped by column count and the **largest group wins**, so a
   stray numeric header line (`3401` + `0.956998`) cannot become the column template. `NaN`
   padding rows are kept. Column 0 → $Q$, column 1 → $S_\mathrm{meas}$, column 2 (if present)
   → $\sigma$. *Code:* `autoScale.js` → `readStogXy()`; Python twin `parsers.read_stog_xy()`.
3. `readDatHeader(dataText)` collects `KEY :: value` lines, and extracts
   `NUMBER_DENSITY ::` (first numeric token → $\rho_0$ in Å⁻³) and `MINIMUM_DISTANCES ::`
   (the **minimum** of the listed values → $r_0$ in Å). *Code:* `autoScale.js` →
   `readDatHeader()`; Python twin `parsers.read_dat_header()`.
4. `dataExtent()` scans for rows where both $Q$ and $S$ are finite and reports
   `{qlo, qhi, count, hasSigma}` — displayed as the file chip
   `name: N pts · Q lo–hi Å⁻¹ · σ`. **`qlo` is the $Q$ of the *first* finite row and `qhi`
   that of the *last*** — they are **not** `min`/`max`. The code assumes $Q$ ascends; a
   non-monotonic file (concatenated banks, a reversed block) therefore produces a wrong Q
   prefill, which can silently exclude most of the data at crop time.
5. The form is prefilled by `selectSource()`. There is **no single precedence chain** — each
   field has its own rule, and two of them overwrite values the user already typed:

   | Field(s) | Rule in `selectSource()` | Precedence |
   | --- | --- | --- |
   | `qmin`, `qmax` | `inp ? inp.q… : (current.q… \|\| Number(extent.q….toPrecision(4)))` | inp → kept form value → data extent (header never consulted) |
   | `rho0` | `inp ? inp.rho0 : (header.numberDensity ?? current.rho0)` | inp → **header** → kept form value |
   | `r0` | `header.minDistance ?? current.r0` | **header** → kept form value (inp not consulted here) |
   | `bAvgSq`, `rCutoff`, `rmax`, `nr`, `lorch`, `enforceCutoff`, `manualA`, `manualB` | `inp ? … : current.…` | inp → kept form value |
   | `enforce` | `inp ? true : current.enforce` | selecting an inp **force-ticks** “Enforce low-r” |

   Two overwrite traps follow directly: selecting a data file whose header carries
   `NUMBER_DENSITY ::` or `MINIMUM_DISTANCES ::` **overwrites the $\rho_0$ and $r_0$ you typed**, and
   selecting a `stog.inp` re-enables low-r enforcement even if you had unticked it. The Q
   prefill is rounded to 4 significant digits (`Number(q.toPrecision(4))`), so a rounded
   $Q_\mathrm{max}$ can land just *below* the last measured point — the crop tolerance is only
   $10^{-12}$ Å⁻¹, so that point is then dropped.

**Outputs:** `dataRef.current = {q, sq, sigma, inp, header, name}` (raw arrays kept
untouched for every later run) and an `inspect` descriptor driving the chips.

---

### Step 3 — Grouped parameter fieldsets: what each knob changes in the engine

The always-visible controls are `DATA`, `SAMPLE` (composition, $\rho_0$, Estimate $\rho_0$, mass
density), `Q WINDOW` (Qmin, Qmax) and the actions (Auto-scale / Advanced / Reset params).
**Advanced** opens five `<fieldset>` groups. *Most* fields map to exactly one key of the
engine config built by `makeConfig()`; the exceptions — `useSigma`, `enforce`/`enforceCutoff`
and `manualA`/`manualB` — are **page-level** and travel in the worker message instead
(`useSigma` decides whether the σ buffer is packed at all, `enforce`/`enforceCutoff` go
through `resolveEnforcement()` into the message's `enforcement` field, and `manualA`/`manualB`
become the message's `a`/`b`; see Step 5).

Every numeric field passes through `numberOr()`, which returns `undefined` for an empty
string **and for anything `Number()` cannot turn into a finite number**. Unparseable text is
therefore indistinguishable from an empty field: `1.0.2` in rmax, `abc` in $\rho_0$, or a unit
suffix like `0.05 A^-3` silently falls through to the next source in that field's chain (the
stog.inp value, the header value, or the built-in default) with no error and no visual cue.
Because a mistyped $\rho_0$ can fall through to the header value or to the 0.05 estimation seed,
this can change the physics of a run — check the provenance JSON's `config` block to see what
was actually used. *Code:* `AutoStogPage.jsx` → `numberOr()`.

**Amplitude & offset**
- *High-Q* (`c1Mode`, default `sweep`) — `sweep` runs `levelSweep()` first and anchors the
  offset to the measured level, $b = 1 - aL$, leaving a single amplitude degree of freedom;
  `joint` restores the original 2-dof least-squares fit of $(a,b)$. If no statistically flat
  window exists, `sweep` degrades to the joint fit and the effective mode is reported as
  `joint` (`c1ModeEffective`).
- *Amplitude* (`amplitudeCriterion`, default `density`) — `density` pins $a$ with the low-$r$
  density limit inside the self-consistent Fourier-filter loop; `fz` pins it in closed form
  from the $Q\to0$ Faber-Ziman limit $S(0) = 1 - \langle b^2\rangle/\langle b\rangle^2$
  (no loop, `iterations = 0`). `fz` requires $\langle b^2\rangle$ **and** `c1Mode='sweep'`;
  `makeConfig()` throws otherwise. The page's `fzMissing` guard pre-blocks the Auto-scale
  button for the **missing-$\langle b^2\rangle$ half only** —
  `form.amplitude === 'fz' && numberOr(form.bSqAvg) === undefined && !form.formula.trim()`,
  which never inspects `c1Mode`. Amplitude = “Faber-Ziman Q→0” with High-Q = “Joint 2-dof”
  therefore leaves the button enabled and fails at run time inside `makeConfig()` with
  `amplitudeCriterion='fz' requires c1Mode='sweep'`, surfaced as the red error banner.
- *Robust* (`robust`, default **on**) — 3 passes of Huber IRLS re-weighting (MAD scale,
  $c = 1.345$) applied per residual block inside `solveAffine()`.
- *σ column* (`useSigma`, default **on**) — page-level, not a config key: it decides whether
  the third data column is transferred to the worker at all. When present, the high-$Q$ (C1)
  rows are weighted by $1/\max(\sigma_i,10^{-12})$, normalized to unit mean over the tail
  block.
- *Despike* (`despike`, default **off**) — drops rolling-median outliers before any
  transform: window 7 points, threshold $6\times$ MAD ($1.4826\times$ median absolute
  residual). Honest warning carried in the tooltip and the engine docstring: it also flags
  real Bragg maxima on crystalline data.

**Coefficients** — `bAvgSq` ($\langle b\rangle^2$, barns) and `bSqAvg` ($\langle b^2\rangle$,
barns). These *override* the composition-derived neutron values and are the required route
for x-ray data ($\langle b\rangle^2 = 1$ for normalized data,
$\langle b^2\rangle = \langle Z^2\rangle/\langle Z\rangle^2$). A `coefficients` memo displays
the values actually
in effect and raises a **⚠ shadowed** chip when an override differs from the typed
composition's Sears value by more than 2 % — the stale-session trap. The browser predicate is
$|\mathrm{override} - \mathrm{fz}| > 0.02\,|\mathrm{fz}|$, evaluated **independently for
$\langle b\rangle^2$ and $\langle b^2\rangle$** (`AutoStogPage.jsx` → the `coefficients` memo);
the CLI's check is different in three ways — see cross-engine difference 7 below.

**Transform** — `rCutoff` (Fourier-filter radius, default 1.0 Å), `rmax` (default 50 Å),
`nr` (default 5000, rounded to an integer), `lorch` (default off; window
$\sin(\pi Q/Q_\mathrm{max})/(\pi Q/Q_\mathrm{max})$ applied to $F(Q)$ before transforming),
`lowQCorrection` (default **on**; the analytic correction for the unmeasured $[0,Q_\mathrm{min}]$
range, extrapolating $S$ linearly from $S(Q_\mathrm{min})$ to `effectiveS0Target(config)`).
That target resolves as: explicit `s0Target` → $1-\langle b^2\rangle/\langle b\rangle^2$ →
**0 when $\langle b^2\rangle$ is unknown**. The page never sets `s0Target` (it stays `null`),
so with no composition and no $\langle b^2\rangle$ override — e.g. $\langle b\rangle^2$ taken
from a stog.inp alone — the correction extrapolates to $S(0)=0$, *not* to a composition-derived
value, and the teal `S(0) FZ target` guide is not drawn at all.
[Step 7 — The composition-derived $S(0)$ target](#step-7--the-composition-derived-s0-target)
covers the resolution rule in full.

**Low-r region** — `r0` (closest approach; empty ⇒ header value, else the stog.inp peak line,
else detected from the data), `rFitMin` (default $r_\mathrm{cutoff} + 0.2$ Å),
`rFitMax` (default $r_0 - 0.25$ Å, or $r_\mathrm{fit,min} + 1.0$ Å when $r_0$ is unknown),
`enforce` (default **on**) and `enforceCutoff`.

**Fixed scaling** — `manualA` / `manualB` and the **Run fixed (a, b)** button, which calls
`runScaling('manual')`. The two values resolve differently:
$a \leftarrow$ typed value → `inp.a` (and the run is refused up front with
*fixed scaling needs a (and optionally b)* only when `manualA` is empty **and** no stog.inp is
loaded); $b \leftarrow$ typed value → `inp.b` → `0`. Consequence worth knowing: with a
stog.inp loaded and **both fields blank**, **Run fixed (a, b)** reproduces the inp's hand
scaling exactly.

The whole form (plus the Advanced open/closed state) is mirrored into
`sessionStorage['autostog-session']` on every change and restored on mount; uploads are not
persisted. **Reset params** restores `EMPTY_FORM` and immediately re-prefills from the
selected file.

---

### Step 4 — Config resolution and the $\rho_0$ self-consistency trigger

**Inputs:** the form, the parsed `stog.inp` (may be null), the `.dat` header (may be empty),
and the run mode (`'auto'` / `'manual'`).

**Operation** (`resolveConfig(form, inp, header, mode)`), in this exact order:

1. **Coefficients.** $\langle b\rangle^2 \leftarrow$ form → `inp.bAvgSq`;
   $\langle b^2\rangle \leftarrow$ form. If a composition string is present,
   `faberZiman(formula)` fills whichever is still undefined:
   $$\langle b\rangle^2 = \Big(\sum_i c_i b_i\Big)^2,\qquad
     \langle b^2\rangle = \sum_i c_i b_i^2,$$
   with $c_i$ the atom fractions from `parseFormula()` (supports decimals and parentheses,
   e.g. `Sr0.5Ba0.5TiO3`, `Al2(SO4)3`) and $b_i$ the bound coherent neutron scattering
   lengths from the 89-element NIST/Sears table `COHERENT_B_FM` (fm; real part for the
   complex-$b$ absorbers B, Cd, Dy, Eu, Gd, In, Sm). Both are divided by 100 to reach barns.
   `faberZiman` refuses near-null-matrix compositions
   ($\langle b\rangle^2 < 10^{-4}\langle b^2\rangle$). Missing $\langle b\rangle^2$ is a hard
   error.
2. **`rCutoff`** ← form → `inp.rCutoff` → 1.0.
3. **`r0`** ← form → header `MINIMUM_DISTANCES` → the stog.inp peak line
   $\max(\texttt{peakCutoff}, \texttt{peakRmin})$, but only if it leaves a non-empty default
   fit window, i.e. only when $\max(\cdot) - 0.25 > r_\mathrm{cutoff} + 0.2$.
4. **Q window** ← form → `inp`. Missing either is a hard error
   (`set the Q window (Qmin and Qmax)`).
5. **$\rho_0$ chain** ← form → `inp.rho0` → header `NUMBER_DENSITY` → mass density + composition:
   $$\rho_0\ [\text{Å}^{-3}] \;=\; \rho_m\ [\mathrm{g\,cm^{-3}}]\times\frac{N_A}{10^{24}}\times\frac{n}{M},$$
   with $N_A = 6.02214076\times10^{23}$, $n$ the atoms per formula unit and $M$ the molar
   mass in u (ADDIE convention) — `numberDensityFromMassDensity()` / `molarMass()` using the
   `ATOMIC_MASS_U` (CIAAW) table for the same 89 elements.
6. **The self-consistency trigger.** If $\rho_0$ is *still* unresolved:
   - and $\langle b^2\rangle$ is known → `rho0 = 0.05` is used as a **seed only** and
     `wantEstimate = true` is returned, so the worker replaces it with the self-consistent
     estimate before fitting;
   - otherwise → hard error: `number density unknown: set ρ₀ or mass density — or give a
     composition and Auto StoG estimates ρ₀ self-consistently`.
7. `makeConfig({...})` validates eagerly: $\rho_0>0$ finite, $\langle b\rangle^2>0$ finite,
   $Q_\mathrm{max}>Q_\mathrm{min}$, `nr` a positive integer, `rmax > 0`, `c1Mode` ∈
   {sweep, joint}, `amplitudeCriterion` ∈ {density, fz} with its two preconditions, and a
   non-empty low-$r$ fit window.

One deliberate mode-dependent override: for `mode === 'manual'` the amplitude criterion is
forced to `'density'`, because a fixed-$(a,b)$ run never uses the criterion and a leftover
`'fz'` selection would otherwise fail `makeConfig`'s validation for no reason.

`resolveEnforcement(form, inp)` runs alongside and returns one of three things:
`null` (checkbox off), the string `'auto'` (checkbox on, no cutoff resolvable — the worker
enforces at the detected $r_0$), or
`{cutoff, peakRmin, peakRmax}`. When the cutoff came from the stog.inp *and the user has not
typed one*, `peakRmin`/`peakRmax` are the inp's first-peak window; otherwise both collapse to
the cutoff (a flat replacement).

**Code:** `AutoStogPage.jsx` → `resolveConfig()`, `resolveEnforcement()`;
`autoScale.js` → `makeConfig()`, `faberZiman()`, `parseFormula()`,
`numberDensityFromMassDensity()`. Python twin: `scaling_cli.py` → `_build_config()`,
`_default_r0()`, `_resolve_enforcement()`.

---

### Step 5 — The worker call protocol

**Inputs:** the resolved config, mode, optional fixed $(a,b)$, the enforcement descriptor, and
the data arrays.

**Operation** (`postJob()` + `packedData()`): a single module Worker
(`new Worker(new URL('../workers/autoScaleWorker.js', import.meta.url), {type:'module'})`) is
created lazily and reused; it is terminated on unmount. Each job gets a monotonically
increasing integer `id`; replies are matched on `id` and the listener is removed, so
overlapping jobs cannot cross-talk. `packedData()` copies $Q$, $S$ and (when `useSigma` and a
third column exist) $\sigma$ into fresh `Float64Array`s and posts their `ArrayBuffer`s as
**transferables** — the copies are neutered in the page, but `dataRef.current` keeps the
originals, so repeated runs are safe.

Message shape:

```
→ { id, kind?, config, mode, estimateRho0, a, b, enforcement, q, sq, sigma }
← { id, ok: true, result: {...} } | { id, ok: false, error: string }
```

Two job kinds:

- `kind: 'estimateRho0'` — runs **only** `estimateRho0(q, sq, makeConfig(config), sigma)` and
  returns the small estimate dict. This is what the **Estimate $\rho_0$** button
  (`runEstimate()`) sends. The button is enabled only when a file is loaded and either a
  composition or an explicit $\langle b^2\rangle$ is available (`canEstimate`).
  The estimator **silently overrides two of your Advanced choices for the duration of the
  estimate**: `let work = { ...config, amplitudeCriterion: 'density', c1Mode: 'sweep' }`. The
  whole method rests on comparing a $\rho_0$-*dependent* (density-limit) amplitude against a
  $\rho_0$-*independent* (Faber-Ziman) one, which is degenerate in `fz` mode ($a \equiv a_\mathrm{fz}$,
  concordance identically 1), so the estimate is always obtained with the density criterion
  and the level sweep no matter what the Amplitude / High-Q selectors say. Your selections
  still govern the scaling run that follows.
- no `kind` (a scaling job) — if `estimateRho0: true` was requested, the worker runs the
  estimate **first**; a non-converged estimate throws (the page shows the physics message
  rather than fitting with a garbage density), otherwise `config.rho0` is replaced by the
  estimate and everything downstream uses it. Then it calls
  `scalePipeline(..., a, b, {mode:'manual'})` for `mode === 'manual'`, else
  `autoscale(q, sq, config, sigma)`.

**Outputs of a scaling job** (`autoScaleWorker.js`): the scalars
`a, b, converged, iterations, history, lowRRms, c1TailMean, mode, c1ModeEffective, nDespiked,
sweep, aFz, r0Detected, windowRefined, rFitWindowUsed, enforcement, rho0Estimate, rho0Used`,
the `summary` dict from `diagnosticsSummary()`, and ten (twelve with enforcement)
`Float64Array` buffers
(`q, sqRaw, sqScaled, sqFiltered, sqFt, r, gk, dr, fk, gm1Unfiltered` plus
`gkEnforced, drEnforced` when enforcement is active) — all transferred back.

**What a fixed-$(a,b)$ run skips.** `mode === 'manual'` takes the `scalePipeline()` path, which
is structurally different from `autoscale()` in ways the cards and plots do not announce:

- **The σ column is never used.** The worker calls `scalePipeline(qArr, sqArr, config, a, b)`
  with no σ argument, and `scalePipeline` itself calls `cropSq(qIn, sqIn, config)` without
  sigma. The **σ column** toggle is therefore completely inert in manual mode.
- **No level sweep runs**, so `sweep = null` and `aFz = null`. `diagnosticsSummary()` then
  emits no `level`, `level_uncertainty`, `level_window`, `asymptote_found`, `a_fz`,
  `amplitude_concordance` or `amplitudes_concordant` — the **High-Q level** card, the
  **Concordance** card and the blue `Level L` guide in Plot 1 all disappear.
- `c1ModeEffective` is the literal string `'manual'`.
- **No two-pass fit-window refinement** (`windowRefined` is always false). The first-shell
  *detection* can still run — but only via the worker's manual-mode recovery (Step 6, item 1),
  i.e. only when the enforcement descriptor is the string `'auto'`. With any other descriptor
  there is no `r0_detected`, so the **First shell $r_0$** card is absent too.
- `iterations = 0`, `converged = true`, `history = []` by construction (`scalePipeline`'s
  `extras` defaults), which is what the **Mode** card shows as `fixed (a, b)`.

The diagnostics that stay meaningful in manual mode are `low_r_rms_pre_enforcement`,
`c1_tail_mean`, `g_window_mean`, `density_limit_satisfied`, and — when $\langle b^2\rangle$ is
known — `fk_qmin` / `fk_q0_theory`.

The fit math itself (level sweep, closed-form affine solve, Fourier-filter self-consistency,
`amplitudeFromFzLimit`, `estimateRho0`, `detectFirstPeakOnset`) is documented in the sibling
sections on the scaling engine and the $\rho_0$ estimator; this section only records *which*
entry point the page calls and *what it does with the answer*.

---

### Step 6 — Post-engine work done inside the worker

Three things are computed in `autoScaleWorker.js` rather than in the engine, because the CLI
computes them outside the engine too and the outputs must match — items 2 and 3 in
`scaling_cli._write_outputs()`, item 1 in `scaling_cli.main()`'s post-run block:

1. **Manual-mode $r_0$ recovery.** `autoscale()` runs the two-pass first-shell detection
   internally; `scalePipeline()` (manual) does not. So when `enforcement === 'auto'` and
   `result.r0Detected` is null, the worker calls
   `detectFirstPeakOnset(r, gFiltered, {searchMin: rCutoff + 0.3})` itself — otherwise a
   ticked "Enforce low-r" would silently become a no-op. This mirrors `scaling_cli.main()`'s
   post-run block (`if enforcement is None and args.enforce is not False`).

   > **Enforcement can still silently become a no-op.** If the descriptor is `'auto'` and
   > `detectFirstPeakOnset()` returns `null` — the dominant $|g|$ feature in
   > $[r_\mathrm{cutoff}+0.3,\ 6.0]$ Å peaks below the `floor = 0.5`, or its left flank never
   > drops below $\max(0.5,\ 0.35\,|g|_\mathrm{peak})$ — the worker sets
   > `effectiveEnforcement = null`. The `*_rmc.gr` / `*_rmc.dr` entries then fall back to the
   > **un-enforced** `gk` / `dr`, the provenance JSON records `enforcement: null`, and nothing
   > on the page says the tick was ignored: the only place enforcement is reported is the
   > **First shell $r_0$** card, which is itself gated on `diagnostics.r0_detected != null` and so
   > is not rendered either. If you rely on enforcement, check `enforcement` in the provenance
   > JSON. *Code:* `autoScaleWorker.js` (the `effectiveEnforcement` block),
   > `autoScale.js` → `detectFirstPeakOnset()`, `AutoStogPage.jsx` (the $r_0$ card guard).
2. **The unfiltered $g(r)-1$** (the classic `scale.gr`), which is not part of the engine
   result:
   $$F(Q) = Q\,[S_\mathrm{corr}(Q)-1],\qquad
     G_\mathrm{PDF}^{\,\mathrm{unfilt}}(r) = \frac{2}{\pi}\!\int F(Q)\sin(Qr)\,\mathrm{d}Q
       \;+\; \mathrm{low-}Q\text{ correction},$$
   $$\big[g(r)-1\big]_\mathrm{unfilt} = \frac{G_\mathrm{PDF}^{\,\mathrm{unfilt}}(r)}{4\pi\rho_0 r},$$
   evaluated with `fqToGpdf(..., {lorch, lowQCorrection, s0Target: effectiveS0Target(config)})`
   — i.e. the *same* discretization (trapezoid sine transform on the data grids) the filter
   used internally. Python twin: `scaling_cli._write_outputs()` →
   `fq_to_gpdf(result.q, sq_to_fq(result.q, result.sq_scaled), result.r, ...)` then
   `gpdf_to_g()`.
3. **Low-r enforcement**, when an enforcement descriptor survives:
   $$g_\mathrm{final}(r) = 0 \quad\text{where } r\le r_\mathrm{cut}\ \mathrm{and}\
     \big(r \ge r_\mathrm{peak,max}\ \mathrm{or}\ r\le r_\mathrm{peak,min}\big),$$
   then $G_K^\mathrm{enf} = \langle b\rangle^2[g_\mathrm{final}-1]$ (so exactly $-\langle b\rangle^2$
   in the zeroed region) and $D^\mathrm{enf} = 4\pi\rho_0 r\,G_K^\mathrm{enf}$. This is the exact
   Fortran `stog_new3.f90` final ripple removal; `firstPeakZero()` in `autoScale.js` is the
   twin of `transforms.first_peak_zero()`.

---

### Step 7 — The numeric readout cards

All values are formatted with `Number(value).toPrecision(digits)` (`fmt()`); `—` means
non-finite/absent. Cards, in DOM order:

| Card | Value shown | Source field |
| --- | --- | --- |
| **Correction** | `a` (5 s.f.), `b` (5 s.f.); sub-line = the stog.inp hand values when an inp was loaded, else `fixed by you` / `auto-fit` | `result.a`, `result.b`, `inp.a`, `inp.b` |
| **Convergence** (or **Mode** when manual) | `✓/✗ N iterations`, or `fixed (a, b)`; sub-line = the last ≤6 fitted `a` values as a `→` trajectory, else the swept level | `converged`, `iterations`, `history[i][0]` |
| **$\rho_0$ self-consistency** (only when the estimate ran) | adopted $\rho_0$ in Å⁻³ (5 s.f.); sub-line = concordance and either `long Q→0 extrapolation` or the pass count | `rho0Estimate.{rho0,concordance,iterations,extrapolated}` |
| **High-Q level** (only when a trajectory exists) | `L` (5 s.f.) ± its spread (2 s.f.), and the window `Q ∈ [q_lo, q_hi]` | `summary.level`, `level_uncertainty`, `level_window` |
| **Fit quality** | `low-r rms` (3 s.f.); sub-line `C1 tail mean` (5 s.f.) | `summary.low_r_rms_pre_enforcement`, `summary.c1_tail_mean` |
| **Density limit** | `satisfied` (green) / `NOT satisfiable` (red) | `summary.density_limit_satisfied` |
| **First shell $r_0$** | detected $r_0$ in Å (4 s.f.); sub-line says whether the fit window was refined and where enforcement was applied | `summary.r0_detected`, `window_refined`, `preview.enforcement.cutoff` |
| **Concordance** | `a_fz / a` (3 s.f.), green when concordant | `summary.amplitude_concordance`, `amplitudes_concordant` |

Definitions of the three verdicts (all from `diagnosticsSummary()` in `autoScale.js`, twin of
`scaling.diagnostics_summary()`):

- **low-r rms** — the *pre-enforcement* residual of the filtered $g$ against 0 over the
  low-$r$ fit window $[r_\mathrm{fit,min}, r_\mathrm{fit,max}]$:
  $\mathrm{rms} = \sqrt{\langle g_\mathrm{filtered}^2\rangle_\mathrm{window}}$. Reported before
  enforcement on purpose — enforcement makes the written files satisfy Keen Eq. 15 by
  construction, so judging quality on them would be circular.
- **C1 tail mean** — $\langle S_\mathrm{filtered}\rangle$ over the high-$Q$ tail
  $Q \ge Q_\mathrm{max} - 0.15\,(Q_\mathrm{max}-Q_\mathrm{min})$; target 1.
- **density limit satisfied** — `abs(mean(g_filtered over the fit window)) < 0.1`.
  **One-sided**: `False` proves the absolute scale is not recoverable from this data;
  `True` is necessary, not sufficient (a smooth low-$Q$ deficiency is absorbed into a biased
  scale with all residuals clean). The card sub-line says exactly this.
- **concordance** — $a_\mathrm{fz}/a$, the ratio of the $\rho_0$-independent $Q\to0$ Faber-Ziman
  amplitude to the fitted density-limit amplitude; `amplitudes_concordant` is
  $|a_\mathrm{fz}/a - 1| < 0.1$. Only present when a level and $\langle b^2\rangle$ exist and
  the criterion is *not* `fz` (in `fz` mode $a \equiv a_\mathrm{fz}$).

The same function emits five further keys that never reach a card but do reach the provenance
JSON, all in the units declared in the symbol table:

| Key | Definition | Units |
| --- | --- | --- |
| `g_window_mean` | $\langle g_\mathrm{filtered}\rangle$ over $[r_\mathrm{fit,min}, r_\mathrm{fit,max}]$ — the quantity the 0.1 density-limit threshold tests | dimensionless |
| `gk_low_r_theory` | $-\langle b\rangle^2$, the exact $G_K$ below $r_0$ | barns |
| `d_r_low_r_slope_theory` | $-4\pi\rho_0\langle b\rangle^2$, the exact $\mathrm{d}D/\mathrm{d}r$ below $r_0$ | barns·Å⁻² per Å |
| `fk_qmin` | $F_K(Q_\mathrm{min})$ = `result.fk[0]`, the first point of the written `_rmc.fq` | barns |
| `fk_q0_theory` | $-\langle b^2\rangle$, the Keen Eq. 14 limit $F_K(0)$ | barns |

The last two are emitted **only when $\langle b^2\rangle$ is known**; their difference is the
second independent check on the absolute scale (the first being the concordance).

When the **Estimate $\rho_0$** button or an auto-triggered estimate returns, the page also writes
the estimate back into the $\rho_0$ **form field**, rounded to 5 significant digits
(`String(Number(rho0.toPrecision(5)))`). Three consequences worth knowing:

- The *next* Auto-scale press reads that number as a user-supplied $\rho_0$ and does **not**
  re-estimate.
- The run currently on screen used the **unrounded** `rho0Estimate.rho0`, while the form now
  holds the 5-digit value. A second Auto-scale press is therefore **not** a re-run of what you
  are looking at — it is a new run at a slightly different density, producing slightly
  different $a$, $b$, $G_K$, $D$ and written files, while the cards read as if nothing changed.
- `estimateRho0` seeds its fixed-point iteration from `config.rho0`
  (`Math.min(Math.max(work.rho0, 1e-4), 1.0)` Å⁻³), so pressing **Estimate $\rho_0$** twice starts
  the second search from the first answer rather than from a fixed seed.

The estimator's exit conditions (`autoScale.js` → `estimateRho0()`) and its `history` shape
are worth stating because both surface in the provenance JSON:

- $|{\rm concordance}-1|\le$ `rtol` $=10^{-3}$ → `converged = true`;
- `result.a <= 0 || concordance <= 0` (non-physical density-limit amplitude) → break with
  `converged = false`;
- the clamped update `rhoNext === rho` (pinned at the $[10^{-4}, 1.0]$ Å⁻³ bound, no progress
  possible) → break with `converged = false`;
- `result.aFz` null / non-finite / ≤ 0 → **throws** (“no usable Faber-Ziman amplitude …”), a
  different failure mode from non-convergence;
- otherwise the loop runs to `maxIter = 8` passes with `converged = false`.

`rho0Estimate.history` rows are `[$\rho_0$, a_density, a_fz, concordance]` per pass — **not** the
scaling `history` shape `[a, b, low_r_rms]` used elsewhere in this section — and the reported
`rho0` / `concordance` / `aDensity` / `aFz` are the **last** row, i.e. the density at which
the last fit was run.

---

### Step 8 — The three plots

All three are `InteractivePlot` (browser-native SVG: hover, legend toggling, drag-zoom,
double-click reset). Series with `role: 'guide'` are drawn dashed/muted and are excluded from
the **palette rotation** and from **hover snapping** — but they are **not** excluded from the
auto-scaled data extent: `orderedSeries` concatenates the guides onto the data series,
`visibleSeries` filters only by the legend's `hidden` set, and the `domains` memo builds
`allX`/`allY` from `visibleSeries.flatMap(...)`, guides included. That is exactly why the
$-4\pi\rho_0\langle b\rangle^2 r$ line reaching to 8 Å and the `initialYDomain` clamp on Plot 2
matter. Series with `defaultHidden: true` start hidden behind a legend click (and, being
hidden, contribute nothing to the extent until toggled on).

Two rendering behaviours affect what a reader believes they are seeing. Each curve is **one
SVG path with no decimation** — all $n_r$ (default 5000) r-points and all cropped Q-points
become path commands. Points are skipped when they fall outside the current x-domain **or**
are non-finite, and skipping never breaks the `d` string: a `NaN` gap is therefore drawn as a
**straight segment bridging the gap** rather than as a hole, and a zoomed curve **terminates
at the last in-domain sample** instead of being clipped at the axis.
*Code:* [`InteractivePlot.jsx`](../../web_app/frontend/src/components/InteractivePlot.jsx) →
`orderedSeries`, `visibleSeries`, `domains`, `seriesShapes`.

**Plot 1 — `S(Q) — scaled and filtered`** (`sqPlot`), x = Q (Å⁻¹), y = S(Q), full width:

| Series | Curve | Pipeline stage |
| --- | --- | --- |
| `Auto-scaled a·S + b` | `sqScaled` | after the affine correction, **before** the Fourier filter |
| `Filtered S(Q)` | `sqFiltered` | after the filter subtraction, `sq − (sq_ft − 1)` |
| `Measured (unscaled)` (hidden) | `sqRaw` | the cropped, despiked input $S_\mathrm{meas}$ |
| `S → 1` (guide) | constant 1 over the whole Q range | Keen Eq. 21 high-Q limit |
| `Level L` (guide, blue `#4c7df0`) | constant $L$, drawn only over the swept window $[q_\mathrm{lo}, q_\mathrm{hi}]$ | the level sweep's measured asymptote of the *unscaled* data |
| `S(0) FZ target` (guide, teal `#0c8599`) | constant $1-\langle b^2\rangle/\langle b\rangle^2$, drawn over the first 1.5 Å⁻¹ of data | Keen Eq. 21 $Q\to0$ limit |

Keen Eq. 21 is the Faber-Ziman definition of $S(Q)$, so **both** limits in that table descend
from it: the codebase cites it for the high-$Q$ normalization
(`scaling.py` module docstring: “C1 — high-Q asymptote: `S_corr -> 1` over a tail window
(Eq. 21)”) and for the $Q\to0$ value ($S(0) = 1-\langle b^2\rangle/\langle b\rangle^2$,
`scaling.py` → `ScalingConfig` docstring and `amplitude_from_fz_limit()`). The teal guide is
drawn only when `guides.s0Target != null`, i.e. only when $\langle b^2\rangle$ is known.

**Plot 2 — `G_K(r) — full range`** (`gkPlot`), x = r (Å), y = $G_K(r)$ in barns.
When enforcement is active the **primary** series is `G_K(r) output (RMC file)` =
`gkEnforced` — i.e. exactly what is written to `*_rmc.gr` — and the un-enforced
`Pre-enforcement fit` is available as a hidden series (the honesty diagnostic). Without
enforcement the single series is `G_K(r) filtered` = `gk`. A guide line at
$-\langle b\rangle^2$ spans $r\in[0,\min(8, r_\mathrm{max})]$. The initial y-domain is
clamped to $[-2.1\langle b\rangle^2,\ 3.2\langle b\rangle^2]$ so the theory level stays
readable; double-click expands to the full data amplitude. The card header prints the low-$r$
fit window actually used.

**Plot 3 — `D(r) — full range`** (`drPlot`), x = r (Å), y = $D(r)$ in barns·Å⁻². Same
enforced/pre-enforcement structure, with the guide line
$D(r) = -4\pi\rho_0\langle b\rangle^2 r$ drawn from the origin to $\min(8, r_\mathrm{max})$ Å.

Convention notes that matter when comparing with other software: the plotted $G_K(r)$ is
**Keen's** $G(r)$ (barns, $\to 0$ at large $r$, flat $-\langle b\rangle^2$ below $r_0$), *not*
the PDFFIT $G(r) = 4\pi\rho_0 r[g-1]$; and $D(r)$ is Keen Eq. 29, *not* the differential
correlation function of some other lineages. The $r$ grid is
$r_i = (i+1)\,r_\mathrm{max}/n_r$, $i = 0\ldots n_r-1$ — it starts at one step, never at 0.
The guides stop at 8 Å (`RMAX_DISPLAY`) because they only mean something below/around the
first shell; the data curves span the full `rmax`.

---

### Step 9 — Export: the written file family

**Trigger:** the **Download .zip (RMCProfile files)** button → `writeFiles()`.

**Stem:** `sanitizeFilename(exportStem || <data file name minus extension> || 'autoscale')`.
`sanitizeFilename` (from [`figureExport.js`](../../web_app/frontend/src/figureExport.js)) is four
steps: strip `^{…}` markup (keeping the exponent text) → replace every run of
non-`[A-Za-z0-9_.-]` characters with `_` → **strip leading and trailing runs of `_`** →
**fall back to the literal `'figure'` if nothing survives**. The `'autoscale'` default applies
only when there is no typed stem *and* no data file name; a typed stem that sanitizes away
(e.g. `***`, which becomes `_` and is then stripped to the empty string) yields `figure.sq`,
`figure_ft.gr`, … and the archive `figure_autoscale.zip`, **not** `autoscale_autoscale.zip`.
Note that `-` and `.` are *inside* the kept character class, so a stem of `--` survives intact
and exports as `--.sq`, `--_ft.gr`, `--_autoscale.zip`.

**Title line** written into every data file:

```
rmc-autoscale (browser): a=<a, 8 s.f.> b=<b, 8 s.f.>
```

The CLI writes the same line with the package version token instead of `(browser)`
(`scaling_cli._write_outputs()`: `f"rmc-autoscale {__version__}: a={a:.8g} b={b:.8g}"`); the
marker difference is deliberate so a browser-produced file is identifiable. The **number
formatters differ** and agree only over the ordinary range: the browser uses
`String(Number(value.toPrecision(8)))` — JS `Number→String`, which switches to exponential
only below $10^{-6}$ and at/above $10^{21}$ and writes single-digit exponents unpadded
(`1e-7`) — while the CLI uses `f"{a:.8g}"`, which switches to exponential below $10^{-4}$ and
at/above $10^{8}$ and zero-pads exponents to two digits (`1e-07`). For a typical
$a \approx 10$, $b \approx -9$ the strings are identical; for very small or very large
corrections the title lines diverge in a way a byte-comparison flags.

**Column format** (`writeStogXy()` in `autoScale.js`, twin of `parsers.write_stog_xy()`):
line 1 is the point count right-aligned in 12 characters, line 2 is the title, then one row
per point, each value formatted `%.16E`-style (16 digits after the decimal point, uppercase
`E`, exponent zero-padded to at least two digits), two leading spaces before each column and
two spaces between columns. Non-finite values are the **one place the two writers are not
byte-identical**: JS `fortranE` falls back to `String(value).toUpperCase()`, giving
`NAN` / `INFINITY` / `-INFINITY`, whereas Python's `f"{v:.16E}"` gives
`NAN` / `INF` / `-INF`.

**Which grid the columns are on.** The four Q-column files (`<stem>.sq`, `<stem>_ft.sq`,
`<stem>_rmc.fq`, `ft.dat`) are written on the **cropped, optionally despiked** grid returned
by `cropSq()`, *not* on the uploaded file's grid: rows with non-finite $Q$ or $S$, rows with
$Q \le 0$, and rows outside $[Q_\mathrm{min}-10^{-12},\ Q_\mathrm{max}+10^{-12}]$ Å⁻¹ are
removed, as are despiked points when Despike is on (fewer than 16 survivors raises *fewer than
16 usable S(Q) points after cropping*). The 12-character count on line 1 is therefore the
number of **surviving** points, generally not the input file's row count, and the written
$S(Q)$ does not extend below $Q_\mathrm{min}$ or above $Q_\mathrm{max}$. The three r-column
files use the full `rGrid(config)` ($n_r$ points) regardless.
[Step 10 — Cropping $Q$](#step-10--cropping-q-to-0-q_mathrmminldots-q_mathrmmax) and
[Step 11 — Despiking](#step-11--despiking-optional-off-by-default) cover those rules.

**The nine zip entries**, in write order:

| # | File | x column | y column | 3rd column | Units / convention |
| --- | --- | --- | --- | --- | --- |
| 1 | `<stem>.sq` | Q | `sqScaled` = $aS_\mathrm{meas}+b$ | — | dimensionless S(Q), unfiltered |
| 2 | `<stem>.gr` | r | `gm1Unfiltered` = $g_\mathrm{unfilt}(r)-1$ (see fallback below) | — | dimensionless; unfiltered transform |
| 3 | `<stem>_ft.sq` | Q | `sqFiltered` | — | dimensionless S(Q), Fourier-filtered |
| 4 | `<stem>_ft.gr` | r | `gk / $\langle b\rangle^2$` = $g_\mathrm{filtered}(r)-1$, **pre-enforcement** | $4\pi\rho_0 r[g-1]$ | y dimensionless; 3rd column Å⁻² |
| 5 | `<stem>_rmc.fq` | Q | `fk` = $\langle b\rangle^2[S_\mathrm{filtered}-1]$ | — | **barns** — RMCProfile $F_K(Q)$ |
| 6 | `<stem>_rmc.gr` | r | `gkEnforced` (else `gk`) | — | **barns** — RMCProfile Keen $G_K(r)$ |
| 7 | `<stem>_rmc.dr` | r | `drEnforced` (else `dr`) | — | **barns·Å⁻²** — RMCProfile $D(r)$ |
| 8 | `ft.dat` (fixed name) | Q | `sqFt` | — | dimensionless; the Fourier-filter correction section |
| 9 | `<stem>_provenance.json` | — | — | — | JSON, see Step 10 |

Two things about these entries that the table cannot carry:

- **Entry 4 is reconstructed, and is never enforced.** The page does not receive
  $g_\mathrm{filtered}-1$ from the worker; it rebuilds it by dividing the transferred $G_K$ back
  by $\langle b\rangle^2$ — `const gm1 = series.gk.map((value) => value / config.bAvgSq);` —
  so the y column (and the third column $4\pi\rho_0 r\,[g(r)-1]$ derived from it) differ from the
  CLI's direct `result.g_filtered - 1.0` by one multiply/divide float round-trip (~1 ulp).
  And `series.gk` is **always the pre-enforcement curve**: only entries 6 and 7 carry the
  enforced curve, so with "Enforce low-r" on, `_ft.gr` keeps the sub-$r_0$ ripples while
  `_rmc.gr` is flat. That disagreement below the cutoff is deliberate.
- **Entry 2 has an undocumented writer fallback.** The call is
  `writeStogXy(series.r, series.gm1Unfiltered || gm1, …)`. If `gm1Unfiltered` were ever absent
  (the worker returns `null` for its buffer), the file would silently receive the **filtered**
  `gk/$\langle b\rangle^2$` instead — making `<stem>.gr` and `<stem>_ft.gr` identical in y with no marker,
  since the title line is the same. In the current worker `gm1Unfiltered` is always produced,
  so the fallback is unreachable in practice; it is recorded here because nothing in the file
  would reveal it if it fired.

Relations the files satisfy exactly, by construction:

$$S_\mathrm{filtered}(Q) = S_\mathrm{scaled}(Q) - \big[S_\mathrm{ft}(Q) - 1\big]$$
$$G_K(r) = \langle b\rangle^2\big[g(r)-1\big],\qquad
  D(r) = 4\pi\rho_0 r\,G_K(r),\qquad
  F_K(Q) = \langle b\rangle^2\big[S_\mathrm{filtered}(Q)-1\big]$$

and, wherever enforcement zeroed $g$: $G_K = -\langle b\rangle^2$ and
$D(r) = -4\pi\rho_0\langle b\rangle^2 r$ to machine precision. The first identity is asserted
on CLI output in `tests/test_scaling_cli.py::StogInpModeTests::test_auto_end_to_end`
(`atol=1e-10`), together with the enforced-region checks (`atol=1e-12`) and the third-column
relation of `_ft.gr` (`rtol=1e-10`).

**Naming vs the CLI.** The browser always uses the *stem* naming
(`<stem>.sq/.gr/_ft.sq/_ft.gr/_rmc.fq/_rmc.gr/_rmc.dr`), identical to
`scaling_cli._OUTPUTS` in `--out-stem` / `--data` mode. The CLI in `stog.inp` mode without
`--out-stem` instead writes the **names declared inside the stog.inp** (`scale.fq`,
`scale.gr`, `scale_ft.sq`, `scale_ft.gr`, `scale_ft_rmc.*`) — so a browser export of a
stog.inp session will not reproduce those file names. `ft.dat` is a fixed name in both.
There is no clobber risk in the browser (a zip is produced, nothing is written in place); the
CLI refuses to overwrite existing outputs without `--force` and defaults into an
`autoscale/` subdirectory.

**Zip container:** [`zipArchive.js`](../../web_app/frontend/src/zipArchive.js) → `buildZip()`
writes a dependency-free ZIP with the **store** method (no compression), CRC-32 per entry,
zeroed timestamps, and a standard end-of-central-directory record; the archive name is
`<stem>_autoscale.zip` and it is handed to `downloadBlob()` (object URL + synthetic anchor
click, URL revoked on the next tick). Files never leave the device.

---

### Step 10 — The provenance JSON

Entry 9 of the zip, `JSON.stringify(..., null, 2)`:

| Key | Content |
| --- | --- |
| `tool` | `"rmc-autoscale (browser engine)"` |
| `source` | the uploaded **data file name** (a bare name — the browser has no path) |
| `a`, `b` | the applied correction |
| `mode` | `"auto"` or `"manual"` |
| `c1ModeEffective` | `"sweep"`, `"joint"` or `"manual"` — what the fit *actually* did |
| `stogInpReference` | `{a, b, yscale, yoffset}` from the loaded stog.inp, else `null` — so an auto run's zip still records the expert's hand values |
| `history` | the iteration trajectory, rows `[a, b, low_r_rms]` |
| `enforcement` | `{cutoff, peakRmin, peakRmax}` or `null` |
| `rho0Estimate` | `{rho0, converged, iterations, concordance, aDensity, aFz, extrapolated, history}` or `null` |
| `config` | the **effective** engine config (camelCase keys), with `rho0` replaced by the value actually used |
| `diagnostics` | the full `diagnosticsSummary()` dict (snake_case keys: `a`, `b`, `converged`, `iterations`, `c1_tail_mean`, `low_r_rms_pre_enforcement`, `g_window_mean`, `r_fit_window`, `gk_low_r_theory`, `d_r_low_r_slope_theory`, `density_limit_satisfied`, plus `r0_detected`/`window_refined`, `level`/`level_uncertainty`/`level_window`/`asymptote_found`, `a_fz`/`amplitude_concordance`/`amplitudes_concordant`, `fk_qmin`/`fk_q0_theory` when available) |

Differences from the CLI's provenance JSON (`scaling_cli.main()` `payload`), which a reader
comparing the two should expect:

- **Truly CLI-only:** `rmc_toolkits_version`, the literal `argv`, absolute `stog_inp` /
  `data_file` paths, and the `outputs` path map.
- **Relocated, not extra.** The CLI nests a `provenance` block copied from the engine result.
  That block is **not** limited to `model` / `q_tail_window` / `r_fit_window` / `n_q_points` /
  `n_despiked` / `level_sweep`: `scaling.py` also puts `mode`, `a`, `b`, and a full snake_case
  `config` dump (26 keys, including the Python-only `c2_bins` and `c1_slope_nuisance`) into it,
  plus `c1_mode_effective` and — when detection fired — `r0_detected` / `window_refined`. So
  the CLI *does* record the effective config, the applied correction and the detected onset;
  it simply nests them under `provenance` where the browser JSON has them at top level as
  `a` / `b` / `mode` / `c1ModeEffective` / `config` and inside `diagnostics`.
- **Both writers record the stog.inp hand scaling.** `stogInpReference` (browser) and
  `stog_inp_reference` (CLI) are both `{a, b, yscale, yoffset}` whenever an inp was loaded,
  in **auto and manual runs alike** — the CLI's condition is only `inp is not None`, and
  `tests/test_scaling_cli.py::StogInpModeTests::test_auto_end_to_end` asserts
  `payload["stog_inp_reference"]["a"]` on an *auto* run. Only the key casing differs.
  (The web API's copy is truncated to `{a, b}` — see Step 0.)
- **`diagnostics.c1_mode` is CLI-only.** After `diagnostics_summary()` returns, the CLI adds
  `summary["c1_mode"] = result.provenance.get("c1_mode_effective", "manual")`, appending
  `", FZ-limit amplitude"` in `fz` mode. The browser's `diagnostics` has no such key — the
  information lives in the top-level `c1ModeEffective`.
- **The enforcement block uses different key casing.** CLI / `_write_outputs` payload:
  `{cutoff, peak_rmin, peak_rmax}`; browser: `{cutoff, peakRmin, peakRmax}`. A script reading
  both provenance files must handle both spellings.
- **`history` really is browser-only** — the CLI's payload has no `history` key (the
  trajectory is not copied out of the engine result).

Two known gaps in the browser JSON, both harmless but worth stating:

- `config.r0` is the *input* value. If the two-pass refinement fired, the effective closest
  approach is in `diagnostics.r0_detected` and the effective window in
  `diagnostics.r_fit_window` — not in `config`.
- The despiked-point count (`nDespiked`) is returned by the worker but is **not** written to
  the JSON (the CLI records it as `provenance.n_despiked`).

---

### Parameters and defaults

| Control | Config key | Default | Units | Effect |
| --- | --- | --- | --- | --- |
| Composition | (derives `bAvgSq`, `bSqAvg`) | none | — | Sears neutron lengths → coefficients + S(0) target |
| $\rho_0$ | `rho0` | resolved (chain) or estimated from a 0.05 seed | Å⁻³ | density line in every real-space conversion |
| mass density | (derives `rho0`) | none | g cm⁻³ | ADDIE conversion via N_A |
| Qmin / Qmax | `qmin` / `qmax` | from stog.inp or the data extent | Å⁻¹ | crop window for fit and transform |
| High-Q | `c1Mode` | `sweep` | — | level-anchored 1-dof vs joint 2-dof |
| Amplitude | `amplitudeCriterion` | `density` | — | low-r density limit vs Q→0 FZ limit |
| Robust | `robust` | on | — | 3 Huber IRLS passes, c = 1.345, MAD scale |
| σ column | (page-level `useSigma`) | on | — | 1/σ weighting of the C1 rows; **inert in fixed-(a, b) mode** |
| Despike | `despike` | off | — | rolling median, window 7, 6·MAD |
| $\langle b\rangle^2$ | `bAvgSq` | from composition | barns | absolute scale of F_K, G_K, D |
| $\langle b^2\rangle$ | `bSqAvg` | from composition | barns | S(0) target, FZ amplitude, $\rho_0$ estimate |
| Filter r-cut | `rCutoff` | 1.0 | Å | Fourier-filter radius |
| rmax | `rmax` | 50.0 | Å | r-grid extent |
| r points | `nr` | 5000 | — | r-grid size; step = rmax/nr |
| Lorch | `lorch` | off | — | sinc window on F(Q) |
| Low-Q corr. | `lowQCorrection` | on | — | analytic [0, Qmin] correction to S(0) target |
| $r_0$ approach | `r0` | header → stog.inp → detected | Å | sets the low-r fit window top |
| Fit win min | `rFitMin` | rCutoff + 0.2 | Å | low-r window bottom |
| Fit win max | `rFitMax` | $r_0$ − 0.25 (or rFitMin + 1.0) | Å | low-r window top |
| Enforce low-r | (page-level) | on | — | classic stog ripple removal on the RMC files |
| Cutoff | (page-level `enforceCutoff`) | stog.inp peak cutoff, else detected $r_0$ | Å | enforcement radius |
| a / b | (manual mode) | stog.inp hand values when loaded | — | fixed correction, skips the fit |
| — | `qTailFrac` | 0.15 | — | C1 tail = top 15 % of the Q window (no UI) |
| — | `c2Weight` | 1.0 | — | relative weight of the low-r block (no UI) |
| — | `fitOffset` | true | — | joint mode only (no UI) |
| — | `despikeWindow` | 7 | points | rolling-median window (library only — no UI, no CLI flag, no API field; `--despike` / the `despike` payload key toggle the filter but cannot retune it) |
| — | `despikeNsigma` | 6.0 | — | despike threshold in MAD units (library only — no UI, no CLI flag, no API field) |
| — | `s0Target` | `null` | — | explicit low-Q target; the page never sets it, so it always resolves through `effectiveS0Target()` (no UI) |
| — | `maxIter` / `tol` | 50 / 1e-6 | — | self-consistency loop stopping rule (no UI) |
| — | level sweep | minWidth 3.0 Å⁻¹, 80 grid edges, ≥24 pts, 2σ slope test | — | not exposed |
| — | $r_0$ detection | search 1.0–6.0 Å (from rCutoff+0.3), 35 % of peak, floor 0.5 | — | not exposed |
| — | $\rho_0$ estimate | rtol 1e-3, ≤8 passes, ρ clamped to [1e-4, 1.0] Å⁻³; also exits on a ≤ 0 / concordance ≤ 0, or a clamp-pinned update; throws when no usable a_fz | — | not exposed |

---

### Cross-engine agreement (browser vs Python)

Everything this page runs also exists in Python. The parity contract is pinned by
`web_app/frontend/src/__tests__/autoScale.test.js` against goldens from
`tests/generate_autoscale_fixture.py`:

- level sweep level: relative error < 1e-9; window edges to 9 decimals; admissible-window
  count exactly equal;
- auto-scale $(a,b)$: relative error < 1e-6; `iterations` exactly equal; `lowRRms` < 1e-5;
  `c1TailMean` < 1e-8;
- FZ-amplitude mode $(a,b)$: < 1e-6, with `iterations === 0`;
- manual pipeline sampled $G_K$ / $S_\mathrm{filtered}$ values: 9 decimals;
- two-pass $r_0$ detection: 9 decimals, same `windowRefined` flag;
- **$\rho_0$ self-consistency: only ~1e-4 relative.** The test and `AGENTS.md` attribute this to the
  fixed-point iteration compounding summation-order float noise against the `rtol = 1e-3`
  stopping rule. That explanation is incomplete: there is also a genuine algorithmic
  divergence in exactly this path (behavioural difference 6 below), which is the more likely
  dominant term. Treat ~1e-4 as the honest bound on the $\rho_0$ estimate's cross-engine agreement,
  not as pure float noise.

Genuine behavioural differences between this page and the `rmc-autoscale` CLI (not
floating-point noise):

1. **σ-column validation.** The CLI (`_load_dataset()`) and the API (`_cached_scaling()`)
   *discard* the whole σ column if any σ is non-finite or ≤ 0 on the usable rows. The page
   does not: `packedData()` forwards column 3 as-is whenever "σ column" is ticked, and a
   `NaN` there propagates into the weighted C1 rows. Untick "σ column" if the third column of
   your file is not a clean uncertainty.
2. **Estimate + manual.** The CLI refuses `--estimate-rho0` together with
   `--manual/--scale/--offset`. The page permits it: a fixed-$(a,b)$ run with an empty $\rho_0$ and
   a composition will seed 0.05, run the estimator (which internally does density-limit auto
   fits), and then apply your $(a,b)$ with the estimated $\rho_0$.
3. **Fixed-b fallback.** With an inp loaded, a typed `a` and a cleared `b`, the page falls
   back to the inp's `b`; the CLI with `--scale` and no `--offset` uses `b = 0`.
4. **Enforcement window from a stog.inp.** Because the page prefills the Cutoff field from
   `inp.peak_cutoff`, `resolveEnforcement()` treats the cutoff as user-supplied and collapses
   the first-peak window to `peakRmin = peakRmax = cutoff`, whereas the CLI keeps the inp's
   `peak_rmin`/`peak_rmax`. The zeroed set is *identical* whenever `peak_rmin ≥ cutoff` (the
   case in every validation run, e.g. cutoff 2.48 with window 2.65–3.1); it differs only for
   inputs whose first peak starts *below* the cutoff, where the browser zeroes a band the
   Fortran/CLI would keep.
5. **Python-only knobs.** `c2_bins` (binned C2 rows) and `c1_slope_nuisance` (linear tail
   drift term) exist in `ScalingConfig` but have no counterpart in `autoScale.js`. Their
   defaults are 0 / `False`, so default runs agree; a Python run that sets them cannot be
   reproduced in the browser.
6. **Low-Q S(0) target inside the fit loop.** `autoScale.js` → `autoscalePass()` calls
   `fourierFilter(q, sqScaled, r, {rho0, cutoff, lorch, lowQCorrection})` with **no**
   `s0Target`, so the per-iteration filter falls back to the destructuring default
   `s0Target = 0` and extrapolates to $S(0)=0$; Python's `scaling.py` → `_pipeline()` passes
   `s0_target=config.effective_s0_target` on **every** loop iteration. Only the final
   `scalePipeline()` call passes `effectiveS0Target(config)` in both engines (and both use it
   in the affine solve's low-Q basis). Consequence: with `lowQCorrection` on **and**
   $\langle b^2\rangle$ known, the browser's converged $(a,b)$ is the fixed point of a slightly
   different map than Python's. The parity fixture does not cover this — the auto / detect /
   manual goldens use a `base` config with no `b_sq_avg` (so `effective_s0_target == 0` there
   too) and the `fz` golden skips the loop entirely (`iterations = 0`). The one golden that
   *does* exercise the loop with a nonzero $S(0)$ target is the $\rho_0$ estimate — precisely the
   comparison that agrees only to ~1e-4.
7. **The "shadowed coefficients" check.** The browser flags
   $|\mathrm{override} - \mathrm{fz}| > 0.02\,|\mathrm{fz}|$ — relative to the **Sears** value —
   independently for $\langle b\rangle^2$ **and** $\langle b^2\rangle$, as a persistent ⚠ chip.
   The CLI (`scaling_cli.py`) tests
   `abs(coefficients.b_avg_sq_barn - b_avg_sq) > 0.02 * abs(b_avg_sq)` — relative to the
   **configured** value — on $\langle b\rangle^2$ **only**, and emits a one-line `stderr`
   warning. The same numbers can therefore be flagged in one engine and not the other, and a
   shadowed $\langle b^2\rangle$ is never reported by the CLI at all.
8. **Mass density with no composition.** The page's `resolveConfig()` converts only when
   *both* are present (`if (massDensity !== undefined && formula)`); otherwise it silently
   ignores the typed mass density and falls through to the 0.05 seed + self-consistent
   estimate, or to the hard `number density unknown` error. The CLI refuses explicitly with
   `--mass-density needs --formula to convert to rho0`. A user who types a mass density with
   x-ray coefficients (no composition string) therefore gets a silently different run in the
   browser.

---

### Caveats / what this is not

- **The tab is not in the shipped build.** `App.jsx`'s `SHOW_AUTO_STOG = false` removes both
  the nav button and the page mount. Everything above documents the code as written; flipping
  the constant to `true` is what exposes it.
- **A ticked "Enforce low-r" can be silently ignored.** When the cutoff is left to `auto` and
  first-shell detection fails, the RMC files are written un-enforced, `enforcement` is `null`
  in the provenance JSON, and no card says so. Verify in the JSON if it matters.
- **Nothing here is a measurement of the absolute scale on its own.** `density_limit_satisfied`
  is one-sided; a smooth low-$Q$ deficiency is absorbed into a biased $a$ with every residual
  clean. The concordance card ($a_\mathrm{fz}/a$) is the trust metric, and it exists only when
  a composition (or explicit $\langle b^2\rangle$) is provided.
- **The written RMC files satisfy the Keen low-$r$ limits by construction** when "Enforce
  low-r" is on (default). Judge the fit only by the pre-enforcement numbers reported on the
  page (`low-r rms`, `C1 tail mean`, concordance) — never by the flatness of `*_rmc.gr` below
  the cutoff.
- **The Sears table is neutron.** For x-ray data the composition route is wrong; set
  $\langle b\rangle^2$ (usually 1 for normalized data) and
  $\langle b^2\rangle = \langle Z^2\rangle/\langle Z\rangle^2$ under Advanced → Coefficients.
  The page warns when
  overrides shadow a typed composition by > 2 %, but it cannot know which is correct.
- **$\rho_0$ estimation needs a composition and a flat high-$Q$ level**, and is flagged
  `extrapolated` whenever $Q_\mathrm{min} > 1.0$ Å⁻¹ (the $Q\to0$ extrapolation is then longer
  than the ~1 Å⁻¹ of data it rests on) — treat such an estimate as a starting point, not a
  measurement. A non-converged estimate aborts the run with the physics message rather than
  fitting with a garbage density. The estimator also forces `amplitudeCriterion = 'density'`
  and `c1Mode = 'sweep'` internally, so the number it reports is not conditioned on your
  Amplitude / High-Q selections; and the value written back into the form is rounded to 5
  significant digits, so re-pressing **Auto-scale** is a *new* run, not a re-run of the results
  on screen.
- **Despiking also flags real Bragg maxima** on crystalline data (~12 % of points on the
  reference PG3 dataset). Use it only for detector glitches.
- **Only the canonical single-dataset, filter-enabled `stog.inp` layout is parsed.** Multiple
  datasets, a nonzero Q offset, a nonzero second y-offset, the interactive "try again" loop,
  or filter-off inputs raise rather than mis-parse.
- **The engine is a *pre*-processor.** It does not read the run folder, does not know about
  `.rmc6f` configurations, and none of the other tabs consume its output automatically — you
  download the zip and feed the files to RMCProfile yourself.
- **Uploads are not persisted.** Only the parameter form survives a reload
  (`sessionStorage`); re-select the file after a refresh.
