# Structure page — algorithm reference

> Part of the [Algorithms and Math Reference](../ALGORITHMS.md). Every step is anchored to the source; if this document and the code disagree, **the code wins**.

Where the configuration puts atoms inside one unit cell: the 2-D kernel-density slice, the slab-in-cell projection, and the folded 3-D unit-cell view.

## Contents

- [Structure page — KDE density slices](#structure-page--kde-density-slices)
  - [What this page shows](#what-this-page-shows)
  - [Step 1 — Load atom positions and fold the supercell into one unit cell](#step-1--load-atom-positions-and-fold-the-supercell-into-one-unit-cell)
  - [Step 2 — Choose the slice plane: normal, in-plane axes, depth coordinate](#step-2--choose-the-slice-plane-normal-in-plane-axes-depth-coordinate)
  - [Step 3 — Restore periodicity by tiling neighbour images](#step-3--restore-periodicity-by-tiling-neighbour-images)
  - [Step 4 — Slab selection in depth, and the fraction-vs-Ångström contract](#step-4--slab-selection-in-depth-and-the-fraction-vs-ångström-contract)
  - [Step 5 — Deterministic pseudo-random subsampling to 6000 fit points](#step-5--deterministic-pseudo-random-subsampling-to-6000-fit-points)
  - [Step 6 — The Gaussian kernel: bandwidth matrix and normalization](#step-6--the-gaussian-kernel-bandwidth-matrix-and-normalization)
  - [Step 7 — Evaluate on the grid (CPU loop and WGSL compute shader)](#step-7--evaluate-on-the-grid-cpu-loop-and-wgsl-compute-shader)
  - [Step 8 — Optional $\log_{10}$ compression](#step-8--optional-log_10-compression)
  - [Step 9 — Contour extraction](#step-9--contour-extraction)
  - [Step 10 — Density → colour, and the affine map onto the real cell](#step-10--density--colour-and-the-affine-map-onto-the-real-cell)
  - [Step 11 — Companion panels (briefly)](#step-11--companion-panels-briefly)
  - [Request lifecycle, caching and determinism](#request-lifecycle-caching-and-determinism)
  - [What the test suite actually checks (and what it does not)](#what-the-test-suite-actually-checks-and-what-it-does-not)
  - [Parameters and defaults](#parameters-and-defaults)
  - [Python vs JavaScript: exact parity table](#python-vs-javascript-exact-parity-table)
  - [Caveats / what this is not](#caveats--what-this-is-not)
- [Structure page — Slab In Cell projection and the 3D unit-cell view](#structure-page--slab-in-cell-projection-and-the-3d-unit-cell-view)
  - [What this page shows](#what-this-page-shows-1)
  - [Step 1 — Parse the `.rmc6f` configuration and fold it into one cell](#step-1--parse-the-rmc6f-configuration-and-fold-it-into-one-cell)
  - [Step 2 — Unit-cell basis and its normalization](#step-2--unit-cell-basis-and-its-normalization)
  - [Step 3 — Defining the slice: normal, in-plane basis, depth range](#step-3--defining-the-slice-normal-in-plane-basis-depth-range)
  - [Step 4 — The plane mapper: crystal plane coordinates → canvas pixels](#step-4--the-plane-mapper-crystal-plane-coordinates--canvas-pixels)
  - [Step 5 — Slab In Cell: band geometry, cell outline, atom markers](#step-5--slab-in-cell-band-geometry-cell-outline-atom-markers)
  - [Step 6 — Dragging the band (cursor → slice position)](#step-6--dragging-the-band-cursor--slice-position)
  - [Step 7 — The 3D "Folded Unit Cell" scene](#step-7--the-3d-folded-unit-cell-scene)
  - [Step 8 — Element colours (shared by the slab canvas, the 3D view and the legend)](#step-8--element-colours-shared-by-the-slab-canvas-the-3d-view-and-the-legend)
  - [Step 9 — The KDE Slice canvas (projection only)](#step-9--the-kde-slice-canvas-projection-only)
  - [Step 10 — Export / screenshot rendering](#step-10--export--screenshot-rendering)
  - [Parameters and defaults](#parameters-and-defaults-1)
  - [Caveats / what this is not](#caveats--what-this-is-not-1)

---

## Structure page — KDE density slices

### What this page shows

The **Structure** page (labelled *KDE And Folded Unit Cell* in the app, `Atomic Density` in the
nav) answers one question: *where does the RMCProfile configuration put atoms inside a single
unit cell?* It takes an `.rmc6f` configuration — a supercell of $S_a \times S_b \times S_c$
crystallographic cells containing $10^4$–$10^6$ atoms — folds every atom back into one unit cell,
selects a thin slab of that cell perpendicular to a chosen direction, and renders a 2-D Gaussian
kernel-density estimate (KDE) of the atoms in that slab. Three panels are drawn from the same
state: the **KDE Slice** (the density map), **Slab In Cell** (a side view showing which atoms the
slab caught, draggable to move the slice), and **Folded Unit Cell** (a Three.js point cloud of the
folded configuration with the slab drawn as a translucent band).

The density map is a *statistical* picture of the modelled configuration, not a measurement and not
a Fourier-space quantity. It is closest in spirit to a nuclear-density / MEM map, but it is
produced by kernel smoothing of discrete atom positions from a single RMC snapshot, and it carries
the smoothing width the user chooses.

Two code paths produce it:

| Runtime | Density engine | Status |
| --- | --- | --- |
| **Server-side run source** (a Flask backend with a run *directory*) | `scipy.stats.gaussian_kde` in [`rmc_toolkits/kde.py`](../../rmc_toolkits/kde.py), served by `/api/kde/slice` | **Reference-grade.** Use this for numbers that go in a paper. |
| **Browser-loaded run** (`localRun`) | hand-written kernel sum in [`web_app/frontend/src/workers/localKdeWorker.js`](../../web_app/frontend/src/workers/localKdeWorker.js), optionally on the GPU via [`gpuKde.js`](../../web_app/frontend/src/workers/gpuKde.js) | **Visualization path.** Same estimator, float32 arithmetic on the GPU branch, a different random subsample, a different bandwidth matrix, and a cruder contour tracer. |

**The discriminator is *not* static-vs-Flask mode.** `StructurePage.jsx` branches on
`const isLocalStructure = Boolean(localRun)` and nothing else; `isStaticMode()` is used on this page
only to emit the *"Open a run folder to view the structure."* message when there is no `localRun`.
The **Demo** button in the app header (`App.jsx` → `handleToggleDemo`) is rendered unconditionally,
outside any `staticMode` branch, and it sets `localRun`. So a Flask-mode session that loads the
bundled demo run gets the **browser worker**, not `/api/kde/slice`, even though the backend is
available. `/api/kde/slice` is used only when a server-side directory is the active run source.
Static mode (GitHub Pages / no backend) is always a browser-loaded run, so it always uses the
worker.

The app says so itself: the in-app `InfoBadge` on the KDE panel and the `local-density-note` under
the canvas both read *"The Flask app uses SciPy KDE for reference-grade values"*
([`StructurePage.jsx`](../../web_app/frontend/src/components/StructurePage.jsx)).

#### Notation and units

| Symbol | Meaning | Units |
| --- | --- | --- |
| $\mathbf{f}_i = (f_{i1}, f_{i2}, f_{i3})$ | atom $i$ fractional coordinate *in the supercell box*, as written in the `.rmc6f` file | dimensionless |
| $\mathbf{N} = (N_1,N_2,N_3)$ | supercell repeat counts along $a,b,c$ | dimensionless integers |
| $\mathbf{L}$ | $3\times3$ matrix of supercell lattice vectors (rows) from the `Lattice vectors` block | Å |
| $\mathbf{A}$, rows $\mathbf{a},\mathbf{b},\mathbf{c}$ | unit-cell matrix; row $j$ is $\mathbf{L}_j$ divided by $N_j$ | Å |
| $\mathbf{x}_i = (x_i,y_i,z_i)$ | atom $i$ folded into one unit cell | dimensionless, $\in[0,1)$ |
| $\mathbf{h} = (h,k,l)$ | user-chosen slice normal, entered as three numbers in *fractional index* space | dimensionless |
| $\hat{\mathbf{h}}$ | $\mathbf{h}/\lVert\mathbf{h}\rVert_2$ | dimensionless |
| $\hat{\mathbf{u}},\hat{\mathbf{v}}$ | in-plane axes, orthonormal *in fractional space* | dimensionless |
| $d_i = \mathbf{x}_i\!\cdot\!\hat{\mathbf{h}}$ | depth of atom $i$ along the slice normal | dimensionless |
| $[d_{\min},d_{\max}]$, $\Delta_d = d_{\max}-d_{\min}$ | projection range / depth span of the unit cube along $\hat{\mathbf{h}}$ | dimensionless |
| $z_c$, $\Delta z$ | slider "Slice" and "Thickness" — **fractions of the depth span**, not Å | dimensionless |
| $\mathbf{p}_i=(u_i,v_i)$ | in-plane projection $(\mathbf{x}_i\!\cdot\!\hat{\mathbf{u}},\, \mathbf{x}_i\!\cdot\!\hat{\mathbf{v}})$ | dimensionless |
| $N_\mathrm{src}$ | unique source atoms contributing to the slab (`slabCount`) — see the warning in Step 4 | count |
| $N_\mathrm{img}$ | slab rows including periodic images | count |
| $n$ | fit points actually handed to the estimator (`fitCount`) | count |
| $m$ | periodic-image margin | fractional |
| $\mathbf{C}$ | $2\times2$ sample covariance of the fit points | fractional² |
| $f$ | bandwidth factor (`bw`) | dimensionless |
| $\mathbf{H} = f^2\mathbf{C}$ | kernel bandwidth (covariance) matrix | fractional² |
| $\rho(u,v)$ | estimated density | per unit fractional area of the slice plane |

---

### Step 1 — Load atom positions and fold the supercell into one unit cell

**Inputs.** One `.rmc6f` file found in the selected run folder, plus an optional element filter
(`all` by default).

#### Which `.rmc6f` file, when the folder holds several

Both runtimes implement the same heuristic, and it is not "the only one there":

* `web_app/backend/app.py` → `_find_rmc6f()`. If the path *is* a `.rmc6f` it is used directly.
  Otherwise the directory is globbed for `*.rmc6f` (sorted by name) and every file in the directory
  is tested against `_run_stem_from_output_name()`, a table of output-file patterns with an explicit
  priority: **0** = `<stem>-NN.log`; **1** = `<stem>-EXAFS-*_[QR]_OUTPUT.csv`, `<stem>_FT_XFQ*.csv`,
  `<stem>_[FS]Q*.csv`, `<stem>_bragg*.csv`, `<stem>_PDF*.csv`; **2** = `Frac_coord_<stem>.txt`.
  Matches are sorted by `(priority, lowercased filename)` and the first `.rmc6f` whose stem matches
  wins. If nothing matches, **the alphabetically first `.rmc6f` is used**.
* `web_app/frontend/src/browserData.js` → `chooseStructureFile()` is the parallel implementation,
  keyed on `directory + '/' + stem` so the match is per-subfolder, with the same priority table but a
  **different fallback**: it ends with `return rmc6fFiles[0]` on an unsorted
  `files.filter(f => f.name.endsWith('.rmc6f'))`, i.e. the first `.rmc6f` in directory-walk / file-picker
  order, which need not be the alphabetically first one — so a multi-model folder can resolve to
  different files on the two paths. (Python sorts: `sorted(directory.glob("*.rmc6f"))[0]`.)

**Consequence.** In a multi-configuration folder the page can silently analyse a different
`.rmc6f` than the user has in mind. The chosen path is reported (`source` in the Flask payload,
`source` in the browser structure object) but is not shown prominently in the UI.

**Math.** The `.rmc6f` header supplies $\mathbf{N}$ (`Supercell dimensions`) and the supercell
lattice matrix $\mathbf{L}$ (`Lattice vectors`, three rows). The unit-cell basis is

$$\mathbf{A}_{j} = \mathbf{L}_{j}\,/\,N_j , \qquad j = 1,2,3 ,$$

and the folding of an atom's box coordinate into one cell is the **modulo-1 of the
supercell-scaled fractional coordinate**:

$$\mathbf{x}_i = \left(\mathbf{f}_i \odot \mathbf{N}\right) \bmod 1 \in [0,1)^3 ,$$

with $\odot$ the elementwise product. This is exactly one line of code —
`unit_frac = (atom["coords"] * supercell) % 1.0` — and it is the *only* wrapping convention used
anywhere in this pipeline. NumPy's `%` returns a non-negative remainder for a positive modulus, so
the result is always in $[0,1)$. The JavaScript twin writes `((value * supercell[i]) % 1 + 1) % 1`
because JavaScript's `%` keeps the sign of the dividend; the two are numerically identical for the
same input.

Two equivalent spellings appear in the repo and both are in use:

* `rmc_toolkits/kde.py` → `load_unit_cell_positions()` folds $\mathbf{f}_i$ directly.
* `web_app/backend/app.py` → `structure()` first subtracts the per-atom cell index,
  `reduced = coords - cell_indices/supercell`, then folds. Subtracting the cell index only removes
  an integer *before* the modulo, so the two give the same $\mathbf{x}_i$. The browser parser
  ([`browserData.js`](../../web_app/frontend/src/browserData.js) → `structureFromRmc6f()`) uses the
  direct form and says so in a comment, which matters because the oldest `.rmc6f` variants carry no
  per-atom cell index at all.

`load_unit_cell_positions()` also returns Cartesian positions
$\mathbf{x}_i^{\mathrm{cart}} = x_i\mathbf{a} + y_i\mathbf{b} + z_i\mathbf{c}$ (Å) and
`cell_lengths` $=(\lVert\mathbf{a}\rVert, \lVert\mathbf{b}\rVert, \lVert\mathbf{c}\rVert)$ (Å).
**The Cartesian array is not used by the KDE endpoint** — see the units gotcha in Step 4.

**Atom-line parsing.** `rmc_toolkits/parsers.py` → `iter_rmc6f_atoms()` and its browser twin
`web_app/frontend/src/rmc6f.js` → `parseAtomLine()` both index fields *from the end* of the line so
that the modern 10-field form (`id element [type] x y z ref cellx celly cellz`) and the legacy
9-field form both parse: the reference number and three cell indices are the last four tokens, the
fractional coordinates the three before them. The browser parser additionally accepts a 5–6 field
coordinates-only form (returning `referenceNumber = null` and `cellIndices = null`) and explicitly
**rejects** 7–8-field lines, so a truncated full line is not misread as coordinates. Python has no
5–6-field branch: `iter_rmc6f_atoms()` skips any line with fewer than 9 tokens.

#### Population handed to the KDE (this differs between runtimes)

* **Server-side run source** — `/api/kde/slice` calls `load_unit_cell_positions()` itself (memoized
  by `_cached_positions`, an `lru_cache(maxsize=16)` keyed on path, mtime and element). **Every atom**
  of the selected element enters the estimate; no display sampling is applied, and the element
  filter is applied while parsing, before anything else.
* **Browser-loaded run** — the worker is posted the `points` **memo**, i.e. `structure.points`
  *filtered by `selectedElement`* (`StructurePage.jsx`, `points` `useMemo`). The element filter is
  therefore applied **client-side, after parsing**, and it is the whole mechanism by which the
  element dropdown works in this path; it is the direct counterpart of Flask's
  `load_unit_cell_positions(element=)`. The 3-D view and the Slab-In-Cell panel use the same
  filtered memo, so all three panels agree with each other.

**The browser display cap, precisely.** `structureFromRmc6f()` computes
`stride = max(1, ceil(N_atoms / maxPoints))`, keeps every `stride`-th atom **in file order**, and
then applies a hard `.slice(0, maxPoints)`; `StructurePage.jsx` passes
`STRUCTURE_MAX_POINTS = 1_000_000`. Two consequences worth naming:

1. The stride is applied to **all** atoms, before any element filter, so above the cap a minority
   element is thinned in the same global proportion (it is not given its own quota).
2. RMCProfile writes atoms grouped by reference site, so a stride can alias onto the site
   ordering — exactly the failure mode Step 5 gives as the reason the *fit* subsample is random
   rather than strided. Below one million atoms nothing is dropped and the tension does not arise.

**The Flask display sampler is a different algorithm again.** `/api/structure` does **not** stride:
`web_app/backend/app.py` → `_sample_atoms_by_site()` is a site-stratified quota sampler.
`maxPoints` is clamped to $[100,\;1{,}000{,}000]$ (`MAX_STRUCTURE_POINTS = 1_000_000`); if
`len(atoms) <= max_points` everything is kept and the reported stride is 1. Otherwise atoms are
grouped by `reference_number`, `quota = max(1, max_points // n_sites)`, each group is strided by
`max(1, len(group)//quota)` and truncated to `quota`, and the concatenation is truncated to
`max_points`. The `sampleStride` field it reports is `max(1, len(atoms)//max_points)` and **does not
describe the actual selection**. Note the split this creates in Flask mode: the **KDE uses all
atoms**, while the 3-D view, the Slab-In-Cell panel and the 50-bin auto-centring histogram (Step 4)
all use this sampled array. Above the cap, the picture and the density are drawn from different
populations.

**Outputs.** $\{\mathbf{x}_i\}$ (fractional, folded), $\mathbf{A}$ (Å),
per-element counts.

**Code.** `rmc_toolkits/kde.py` → `load_unit_cell_positions()`, `UnitCellPositions`;
`rmc_toolkits/parsers.py` → `read_cell_vectors()`, `iter_rmc6f_atoms()`;
`web_app/backend/app.py` → `_find_rmc6f()`, `_run_stem_from_output_name()`,
`_sample_atoms_by_site()`, `structure()`;
`web_app/frontend/src/browserData.js` → `chooseStructureFile()`, `structureFromRmc6f()`;
`web_app/frontend/src/rmc6f.js` → `parseAtomLine()`;
`web_app/frontend/src/workers/localStructureWorker.js`;
`StructurePage.jsx` → the `points` `useMemo`.

---

### Step 2 — Choose the slice plane: normal, in-plane axes, depth coordinate

**Inputs.** The *Normal* control (`a`, `b`, `c`, or `Custom` with three numbers, default
`[1, 1, 0]`); default direction is `c`.

**Math.** The presets are defined once on each side and agree exactly:

| Preset | $\mathbf{h}$ | $\hat{\mathbf{u}}$ | $\hat{\mathbf{v}}$ | axis labels |
| --- | --- | --- | --- | --- |
| `a` | $(1,0,0)$ | $(0,1,0)$ | $(0,0,1)$ | b, c |
| `b` | $(0,1,0)$ | $(1,0,0)$ | $(0,0,1)$ | a, c |
| `c` | $(0,0,1)$ | $(1,0,0)$ | $(0,1,0)$ | a, b |

(`SLICE_PRESETS` in `StructurePage.jsx`, `SLICE_ORIENTATIONS` in `web_app/backend/app.py`.)

Everything downstream works in **fractional index space**: the code treats the unit cell as the
unit cube $[0,1]^3$ and takes ordinary Euclidean dot products of fractional triples. The consequence
is crystallographically clean and worth stating explicitly: the set
$\{\mathbf{x} : \mathbf{x}\cdot\mathbf{h} = \mathrm{const}\}$ with $\mathbf{h}=(h,k,l)$ is exactly the
family of lattice planes with **Miller indices $(hkl)$**. So the *Custom* direction box is a Miller
index triple, not a real-space vector; preset `a` slices the $(100)$ planes, and `[1 1 0]` slices
the $(110)$ planes, for any cell metric including triclinic.

The normal is normalized in that same fractional space,
$\hat{\mathbf{h}} = \mathbf{h}/\lVert\mathbf{h}\rVert_2$, which sets only the *scale* of the depth
coordinate, not the plane family.

For a custom normal the in-plane axes are built by Gram–Schmidt against $\hat{\mathbf{h}}$, but the
two runtimes pick a **different seed vector**:

* Python `_orthogonal_axis()` seeds with the Cartesian unit vector along the *smallest-magnitude*
  component of $\hat{\mathbf{h}}$ (`np.eye(3)[argmin(|n|)]`), then $\hat{\mathbf{v}} = \hat{\mathbf{h}}\times\hat{\mathbf{u}}$.
* JavaScript `makeFreePlaneBasis()` seeds with $(1,0,0)$ when $|n_1| < 0.85$ and $(0,1,0)$
  otherwise, then $\hat{\mathbf{v}} = \hat{\mathbf{h}}\times\hat{\mathbf{u}}$.

Both give a right-handed orthonormal frame in the plane, but generally a *different* one. For
$\mathbf{h}=(1,1,0)$, Python returns $\hat{\mathbf{u}}=(0,0,1)$ while JavaScript returns
$\hat{\mathbf{u}}=(\tfrac{1}{\sqrt2},-\tfrac{1}{\sqrt2},0)$. The density field is the same up to
that in-plane rotation/reflection, but **a custom-normal slice is drawn in a different in-plane
orientation on the SciPy path than on the browser path.** The a/b/c presets are unaffected.

#### Zero and near-zero custom directions (the two runtimes disagree)

The three *Direction* inputs are `type="number"` and `updateCustomDirection()` coerces with
`Number(value)`, so **clearing a box yields 0** (`Number('') === 0`). Clearing all three gives
$\mathbf{h}=(0,0,0)$, and the two runtimes handle a zero vector completely differently:

* **JavaScript** — `normalize(vector, fallback)` returns the fallback when
  $\lVert\mathbf{v}\rVert \le 10^{-9}$. `makeSliceConfig()` calls
  `normalize(customDirection, [0, 0, 1])`, so a zero direction **silently becomes a c-slice**, with
  the label still reading `[0 0 0]`. `makeFreePlaneBasis()` carries two further fallbacks
  ($[0,1,0]$ for $\hat{\mathbf{u}}$, $[0,0,1]$ for $\hat{\mathbf{v}}$), but for a unit input normal
  the Gram–Schmidt residual is never shorter than $\sqrt{1-0.85^2}\approx0.53$, so those two are
  unreachable in practice.
* **Python** — `_normalize_vector()` **raises** `ValueError("normal must be a non-zero 3D vector")`
  for $\lVert\mathbf{v}\rVert \le 10^{-12}$, which `kde_slice_endpoint()` turns into an error
  response.

**But the app never triggers the Python raise**, because the frontend sends
`sliceConfig.normal` — already normalized, hence already the $(0,0,1)$ fallback. A cleared custom
direction therefore produces a silent c-slice in *both* modes. The `ValueError` is only reachable by
a hand-written API request with `nx=0&ny=0&nz=0`.

The **depth span** is the extent of the unit cube's corner projections,

$$d_{\min} = \min_{\mathbf{q}\in\{0,1\}^3} \mathbf{q}\cdot\hat{\mathbf{h}}, \quad
d_{\max} = \max_{\mathbf{q}\in\{0,1\}^3} \mathbf{q}\cdot\hat{\mathbf{h}}, \quad
\Delta_d = d_{\max}-d_{\min} = \frac{\lVert\mathbf{h}\rVert_1}{\lVert\mathbf{h}\rVert_2},$$

evaluated over the eight `_CUBE_CORNERS` / `CUBE_CORNERS`. For any single-axis preset $\Delta_d=1$; for
$(1,1,0)$, $\Delta_d=\sqrt2$.

The in-plane plot limits are the corresponding corner projections onto $\hat{\mathbf{u}},\hat{\mathbf{v}}$:
$x_{\mathrm{lim}} = [\min_\mathbf{q}\mathbf{q}\cdot\hat{\mathbf{u}},\ \max_\mathbf{q}\mathbf{q}\cdot\hat{\mathbf{u}}]$
and likewise for $\hat{\mathbf{v}}$. **This is the bounding rectangle of the projected cube, not the
cell cross-section** — see Step 7, where it becomes the evaluation grid.

**Outputs.** $\hat{\mathbf{h}}, \hat{\mathbf{u}}, \hat{\mathbf{v}}$, $[d_{\min},d_{\max}]$, plot
extent.

**Code.** `rmc_toolkits/kde.py` → `_plane_basis()`, `_orthogonal_axis()`, `_normalize_vector()`,
`_CUBE_CORNERS`; `StructurePage.jsx` → `makeSliceConfig()`, `makeFreePlaneBasis()`, `normalize()`,
`projectionRange()`; `localKdeWorker.js` → `makeFreePlaneBasis()`, `normalize()`.

---

### Step 3 — Restore periodicity by tiling neighbour images

**Why.** Folding into one cell (Step 1) throws away every atom's periodic neighbours. A kernel
evaluated near a cell face would then see only "half" the atoms and the density would decay
artificially toward the boundary. It would also make a slab centred at $z_c=0$ empty even when the
structure has a dense layer at $z\approx 1$.

**Math.** For every folded atom, all 26 non-trivial integer translations
$\boldsymbol{\delta}\in\{-1,0,1\}^3\setminus\{\mathbf{0}\}$ are generated and kept when the shifted
point lies in the padded cube:

$$\mathbf{x}_i + \boldsymbol{\delta} \in [-m,\ 1+m]^3, \qquad
m = \min\!\Big(0.5,\ \max\big(0.1,\ 2f,\ \Delta z\big)\Big).$$

$m$ is the **margin** in fractional units: it must cover both the kernel reach (the kernel $\sigma$
scales as $f$ times an $O(1)$ data spread) and the slab depth, so that both the in-plane density and
the depth selection wrap. With the default $f=0.03$ and $\Delta z = 0.08$, $m = 0.1$; the slider
maxima ($f\le0.15$, $\Delta z\le0.5$) push $m$ to at most $0.5$.

Every retained row carries a `source_index` / `sourceIndex` pointing back to its originating atom,
so image duplicates can be counted out later (Step 6).

`kde.py` prepends the original (unshifted) array unconditionally and then appends the surviving
images; the worker runs the margin test over all 27 offsets including $\boldsymbol{\delta}=\mathbf{0}$.
For folded inputs ($\mathbf{x}_i\in[0,1)^3$) the identity offset always passes, so the two are
equivalent — but note the **row order differs**: Python emits one block of originals followed by one
block per surviving offset, while JavaScript emits all surviving offsets of atom 0, then of atom 1,
and so on. Index-based subsampling (Step 5) therefore picks different points even before the two
random streams are taken into account.

#### The margin is a truncation, not an exact wrap

The margin test is applied **independently per axis over the full $3\times3\times3$ offset
product**, so for a population spread through the cell the expected augmented row count is

$$N_{\mathrm{aug}} \;=\; (1+2m)^3 \times (\text{folded population}).$$

Measured directly (`_augment_periodic_images` on 20 000 uniform points): **1.735× at $m=0.1$**
(default), 2.297× at $m=0.16$, and exactly **8× at $m=0.5$** (the slider-driven maximum, where every
one of the 27 offsets survives for every atom). This is the memory and compute multiplier the
augmentation costs.

More importantly, images *beyond* $m$ are simply dropped, so the periodic wrap is exact only up to
the margin. The effective kernel truncation radius is $m/\sigma$ standard deviations, where
$\sigma = f\sqrt{\lambda}$ (Step 6):

| Case | $m$ | $\sigma$ (cell-filling slab) | truncation | dropped kernel weight |
| --- | --- | --- | --- | --- |
| defaults, $f=0.03$ | 0.1 | $\approx 0.0105$ | $\approx 9.5\sigma$ | $\sim e^{-45} \approx 3\times10^{-20}$ |
| slider max, $f=0.15$ | 0.3 | $\approx 0.0525$ | $\approx 5.7\sigma$ | $\sim e^{-16} \approx 9\times10^{-8}$ |

The browser's `exponent > -60` guard (Step 7) is a *second*, independent truncation at a Mahalanobis
radius of $\sqrt{120}\approx 11\sigma$. Which of the two binds depends on $f$: with $m=0.1$ the
crossover is near $f\approx0.026$, so at the default and above **the image margin is the binding
approximation**, and below it the exponent guard is. The image truncation applies to the SciPy
reference path as well — SciPy has no exponent cutoff, but it never sees the discarded images.

**Verification.** Python: `tests/test_kde.py::test_oriented_kde_slice_wraps_density_across_the_cell_boundary`
plants a cluster hugging the $x=0$ face and asserts the $x=1$ edge of the slice sees it to within
20 % of the $x=0$ edge, **and more than $3\times$ what the non-periodic `kde_slice()` reports**;
`test_oriented_kde_slice_wraps_slab_selection_in_depth` puts all 30 atoms at $z\approx0.98$ and
asserts a slab at $z_c=0$ still finds all 30. JavaScript
(`workers/__tests__/localKdeWorker.test.js`) asserts a single atom at $x=0.02$ with $m=0.1$ yields
exactly 2 rows both mapped to source index 0, and the same 20 %-edge-agreement and 30-atom depth
wrap. Its third assertion is **not** the same as Python's: the worker has no non-periodic mode, so
instead of comparing against a truncated reference it compares the wrapped edge against a
cluster-free point in the *same* run (`rightEdge > 2 × density[36][40]`).

**Outputs.** Augmented point set — $(1+2m)^3$ times the folded population in expectation, i.e.
1.73× at the default $m=0.1$ rising to 8× at $m=0.5$, concentrated near faces/edges/corners — plus
the source-index map.

**Code.** `rmc_toolkits/kde.py` → `_augment_periodic_images()`, and the margin line inside
`oriented_kde_slice()`; `localKdeWorker.js` → `augmentPeriodicImages()` (exported for tests) and
the `margin` line inside `computeKde()`.

---

### Step 4 — Slab selection in depth, and the fraction-vs-Ångström contract

**Inputs.** Slider `zCenter` $=z_c \in [0,1]$ (step 0.001) and slider
`Thickness` $=\Delta z \in [0.01, 0.5]$ (step 0.01, default 0.08).

**Math.** The slider values are mapped onto the depth axis and an atom is in the slab iff

$$\left| \mathbf{x}_i\cdot\hat{\mathbf{h}} - \big(d_{\min} + z_c\Delta_d\big) \right| \ \le\ \frac{\Delta z\,\Delta_d}{2}.$$

Python writes this as `center_depth = depth_min + center*depth_span`,
`thickness_depth = thickness*depth_span`, then in `kde_slice()`
`half = 0.5*max(dz,1e-12)` with an inclusive two-sided comparison. JavaScript `makeSlab()`
normalizes first — `normalizedDepth = (x·n̂ − range[0]) / depthSpan` — and tests
`|normalizedDepth − zCenter| ≤ thickness/2`. **The two are algebraically identical.** Python
additionally clamps $z_c$ into $[0,1]$ and floors $\Delta z$ at $10^{-12}$; the JavaScript path relies
on the slider bounds.

**The mask is applied to the augmented set and is not clamped to the cube.** A slab at $z_c=0$ picks
up images at depth $d \approx 0^-$ from atoms whose folded depth is $\approx 1$; that is the whole point of
Step 3. Keep this in mind when reading `slabCount` and the companion panels (below, and Step 11).

#### `slabCount` counts *source atoms with an image in the slab*, not atoms in the slab

`slab_count = np.unique(source_index[mask]).size` (Python) / a `Set` of `sourceIndex` values
(JavaScript). Two readings follow, and neither is the naive one:

* An atom whose own folded coordinate is nowhere near the slab still counts, provided one of its
  periodic images landed inside. The repo's own fixture makes this concrete:
  `test_oriented_kde_slice_wraps_slab_selection_in_depth` puts all 30 atoms at $z\approx0.98$ and
  asserts `slabCount == 30` for a slab centred at $z_c = 0$.
* Conversely one atom can contribute **several rows** to the same slab when the margin admits more
  than one of its images (edge and corner cells), so $N_\mathrm{img} > N_\mathrm{src}$ in general. That inequality is exactly
  why the renormalization $\kappa = N_\mathrm{img}/N_\mathrm{src}$ of Step 6 exists.

The canvas label reads `"{slabCount} atoms in slab"`; read it as *"atoms with at least one periodic
image in the slab"*.

#### The units gotcha (read this before quoting a slab thickness)

`AGENTS.md` states: *"`z` / `dz` are cell-edge fractions at the API/slider boundary, converted to
Ångström inside `kde.py`."* **The first half is right only for the `a`/`b`/`c` presets, and the
second half is stale.** In the code as it stands, nothing is converted to Å anywhere in the KDE
slice pipeline:

* `/api/kde/slice` passes `positions.fractional_positions` — *not* `positions.positions` (the Å
  array) — into `oriented_kde_slice()`, with an inline comment saying "Keep the KDE slice in
  fractional coordinates so non-orthogonal cells can be projected through the actual cell basis in
  the frontend."
* Every subsequent quantity ($t$, the slab half-width, the KDE covariance, the evaluation grid, the
  contour coordinates) is therefore dimensionless.
* Ångströms enter **only at draw time**, when `StructurePage.jsx` maps
  $\hat{\mathbf{u}},\hat{\mathbf{v}}$ through `unitCell.unitVectors` to get the real parallelogram
  (`makeProjectedPlane(vectorFromFraction(uVector, unitCell.unitVectors), …)`).

So the honest statement of the contract is: **$z_c$ and $\Delta z$ are fractions of the projection
range of the unit cube along the chosen normal**, and they only coincide with "fraction of a cell
edge" for the `a`/`b`/`c` presets (where $\Delta_d=1$).

**Converting $\Delta z$ to Ångström.** Combining the code's definitions with the standard
interplanar spacing $d_{hkl}$ of the plane family $(hkl)=\mathbf{h}$, a depth increment
$\Delta d$ corresponds to a real-space distance $\Delta d \cdot \lVert\mathbf{h}\rVert_2 \cdot d_{hkl}$,
hence

$$\text{slab thickness (Å)} \;=\; \Delta z \; \Delta_d \; \lVert\mathbf{h}\rVert_2 \; d_{hkl}
\;=\; \Delta z \; \big(|h|+|k|+|l|\big)\, d_{hkl}.$$

Sanity checks: for the `c` preset in an orthogonal cell this is $\Delta z\cdot c$ (the naive
reading); for $\mathbf{h}=(1,1,0)$ in a cubic cell it is $\Delta z\cdot 2 \cdot a/\sqrt2 = \sqrt2\,a\,\Delta z$,
which is $\Delta z$ times the full $[110]$ diagonal of the cell — as it must be, since $z_c$ sweeps
$0\to1$ across the whole cube. This conversion is *not* performed anywhere in the app; the readout
shows only the dimensionless $z_c$ and $\Delta z$ values.

#### Auto-centring (it overwrites the slider, and more often than you would guess)

An effect in `StructurePage.jsx` histograms the atoms into **50 depth bins** and sets
$z_c = (\mathrm{argmax}+0.5)/50$, "so the view is populated on load (the geometric midpoint can fall in
a gap between atomic layers)". The details matter for reproducing the behaviour:

* **Triggers.** The dependency array is `[points, sliceConfig, pointDepth]`. `points` changes on
  structure load and on every element-filter change; `sliceConfig` is rebuilt by
  `makeSliceConfig(sliceDirection, customDirection)` and so changes on every normal change **and on
  every keystroke in a custom-direction box** — typing a single digit re-runs the search and
  overwrites $z_c$.
* **Binning.** `bin = max(0, min(49, floor(depth·50)))` — a clamped floor, so depth exactly 1 lands
  in bin 49. The argmax scan uses a strict `count > counts[best]`, so **ties resolve to the lowest
  bin**.
* **Population.** It histograms the element-filtered *display* array (strided in the browser,
  site-quota-sampled in Flask — Step 1), **never** the periodic-augmented set. The "densest layer"
  is therefore found with no wrapping at all, on a possibly sampled population, while the density it
  centres is computed on the full wrapped set.

#### The plane-section polygon (`planePolygon` / `planeVertices`)

The cross-section of the unit cube by the plane $\mathbf{x}\cdot\hat{\mathbf{h}} = d_{\mathrm{c}}$
(with $d_{\mathrm{c}}$ = `center_depth`) is computed identically in three places and is used for four
different things: the drawn cell outline, the clip region for the density blit, the point set whose
projection defines the drawing plane bounds, and the 3-D slab caps. The algorithm:

1. For each of the 12 cube edges $(\mathbf{p}_0,\mathbf{p}_1)$ compute
   $d_0 = \mathbf{p}_0\!\cdot\!\hat{\mathbf{h}} - d_{\mathrm{c}}$ and
   $d_1 = \mathbf{p}_1\!\cdot\!\hat{\mathbf{h}} - d_{\mathrm{c}}$.
2. Push $\mathbf{p}_0$ verbatim when $|d_0| \le 10^{-9}$, and $\mathbf{p}_1$ when $|d_1| \le 10^{-9}$
   (a corner lying on the plane).
3. When $d_0 d_1 < 0$ push the crossing $\mathbf{p}_0 + t(\mathbf{p}_1-\mathbf{p}_0)$ with
   $t = d_0/(d_0-d_1)$.
4. Deduplicate: a vertex is dropped if it is within Euclidean distance $10^{-8}$ of one already
   kept. If fewer than **3** unique vertices survive, return an **empty** polygon.
5. Order the survivors about the centroid $\mathbf{c}$ of the unique vertices by the angle
   $\operatorname{atan2}\big(\boldsymbol{\delta}\!\cdot\!\hat{\mathbf{v}},\,\boldsymbol{\delta}\!\cdot\!\hat{\mathbf{u}}\big)$,
   with $\boldsymbol{\delta} = \mathbf{w}-\mathbf{c}$ the offset of vertex $\mathbf{w}$ from the
   centroid.

The section of a cube by a plane is always convex, so step 5 gives a simple polygon (a triangle,
quadrilateral, pentagon or hexagon). **Python re-derives $\hat{\mathbf{u}},\hat{\mathbf{v}}$ inside
`_plane_section_vertices()` from `_plane_basis(normal)` with default axes, ignoring any
`u_axis`/`v_axis` passed in.** For the `b` preset that default $\hat{\mathbf{v}}$ is $(0,0,-1)$
against the preset's $(0,0,1)$, so the vertex cycle comes out reversed; since the polygon is convex
and closed this only flips the winding and nothing downstream depends on it.

`planePolygon` is the same vertex list projected onto the *actual*
$(\hat{\mathbf{u}},\hat{\mathbf{v}})$ of the slice. **Fallback:** if it is empty, `drawKdeSlice()` substitutes the full
extent rectangle $[[x_{\min},y_{\min}],[x_{\max},y_{\min}],[x_{\max},y_{\max}],[x_{\min},y_{\max}]]$,
so the panel degrades to drawing the whole bounding box rather than failing.

#### The Flask query-parameter contract

The endpoint parses its arguments by hand, and one rule surprises API consumers:

| Query arg | Default | Notes |
| --- | --- | --- |
| `dir` | `.` | resolved inside `DATA_ROOT`; `403` outside it |
| `element` | absent → all | `''` and `all` both map to `None` (no filter) |
| `orientation` | `c` | lowercased; if it is `a`, `b` or `c` the **preset wins and `nx/ny/nz` are ignored entirely**; anything else → `custom` |
| `nx`, `ny`, `nz` | `0.0, 0.0, 1.0` | only read when `orientation` is not a preset |
| `z` | `0.5` | fraction of the projection range $\Delta_d$ |
| `dz` | `0.08` | fraction of $\Delta_d$ |
| `bw` | `0.03` | plain `float()`. The `_bw_argument()` helper that accepts `"scott"`/`"silverman"` is **not** wired to this endpoint |
| `grid` | `120` | clamped to 16–400 inside `kde_slice()` |
| `levels` | `8` | contour count |
| `log` | `false` | true only for `'1'`, `'true'`, `'yes'` (case-insensitive); **any other spelling silently means linear** |

The frontend always sends **both** `orientation` and `nx/ny/nz`, so `orientation=c&nx=1&ny=1&nz=0`
returns a c-slice, not a $(110)$ slice.

**Outputs.** The slab point list $\{\mathbf{p}_i\}$ (2-D, in-plane fractional projections),
$N_\mathrm{img}$ = slab rows including images, $N_\mathrm{src}$ = unique source atoms with an image in the slab = `slabCount`.

**Code.** `rmc_toolkits/kde.py` → `oriented_kde_slice()` (depth mapping, margin,
`_plane_section_vertices()`) and `kde_slice()` (mask/`slabCount`); `localKdeWorker.js` →
`makeSlab()`, `planeSectionVertices()`; `StructurePage.jsx` → the auto-centring effect,
`pointDepth()`/`inActiveSlab()`, `planeSectionVertices()`; `web_app/backend/app.py` →
`_slice_orientation_from_request()`, `kde_slice_endpoint()`.

---

### Step 5 — Deterministic pseudo-random subsampling to 6000 fit points

**Why.** The estimator cost is $O(\mathrm{grid}^2 \times n)$. A 52 000-atom sample configuration can
put tens of thousands of points in a thick slab, and the sliders are meant to be interactive. The
cap keeps the slider responsive. The header comment in `kde.py` states the rationale directly:
*"The density estimate is stable well below the full population, and the eval cost scales with the
number of fit points."*

**Why *random* and not every $k$-th point.** RMCProfile writes atoms in a structured order
(by site/reference number, then by cell index). Taking a stride would sample that structure and
alias onto the crystal lattice — you would preferentially keep certain sites or certain cells and
the density map would show a spurious superlattice. A pseudo-random draw without replacement
removes the correlation between selection and position. (Note the tension with Step 1: the
browser's *display* cap above one million atoms **is** a stride, applied before the element filter.)

**Math.** If $N_\mathrm{img} > 6000$, draw a uniform sample of size $n = 6000$ without replacement from the
$N_\mathrm{img}$ slab rows; otherwise use all of them, $n=N_\mathrm{img}$.

The two runtimes use **different generators**, both seeded to a fixed constant so a given input
always yields the same picture:

* **Python** — `np.random.default_rng(rng_seed)` (PCG64) with `rng_seed = 0` (the parameter default;
  neither `oriented_kde_slice()` nor the Flask endpoint ever overrides it), then
  `rng.choice(N, 6000, replace=False)`.
* **JavaScript** — `randomUnit(seed)` with `seed = 0`, a 32-bit *mulberry32*-style integer hash
  (`value += 0x6D2B79F5; …; >>> 0) / 4294967296`), driving a **partial Fisher–Yates shuffle**
  (`sampleWithoutReplacement`): for `index` in $0..n-1$, swap `indices[index]` with a uniformly
  chosen entry from the remaining tail, then take the first $n$. This is a correct uniform
  without-replacement draw given a uniform generator.

Because the streams differ — and because the augmented arrays are in a different row order to begin
with (Step 3) — **the two runtimes fit different 6000-point subsets of the same slab**. They are
both unbiased, so the two density fields agree in expectation, but they are not bitwise comparable,
and (see Step 6) they do not even share the same bandwidth matrix.

**Statistical consequence.** The KDE is an average of $n$ kernels; subsampling raises the
pointwise standard error of the estimate by roughly $\sqrt{N_\mathrm{img}/n}$ relative to using all $N_\mathrm{img}$ points.
For a slab of 20 000 rows capped at 6000 that is a factor $\approx1.8$ on the Monte-Carlo noise of
the density value — visible as slightly grainier fine structure, not as a bias. The **peak
positions and the integrated normalization are unaffected**; only the variance is. `slabCount`
($N_\mathrm{src}$, unique source atoms with an image in the slab) and `fitCount` ($n$, points actually fitted)
are both returned and both printed on the canvas as `"{slabCount} atoms in slab (fit {fitCount})"`,
so the reader can always see when subsampling kicked in.

**Note.** The cap is applied *after* periodic augmentation, so image duplicates consume part of the
6000 budget.

**Verification.** `tests/test_kde.py::test_kde_slice_reports_subsampled_fit_count` builds
`MAX_KDE_FIT_POINTS + 10` coplanar points and asserts `slabCount == 6010` and `fitCount == 6000`.

**Outputs.** $n \le 6000$ fit points; `fitCount`.

**Code.** `rmc_toolkits/kde.py` → `MAX_KDE_FIT_POINTS = 6000` and the `rng.choice` branch in
`kde_slice()`; `localKdeWorker.js` → `randomUnit()`, `sampleWithoutReplacement()`, and the
`fitLimit = 6000` literal in `computeKde()`.

---

### Step 6 — The Gaussian kernel: bandwidth matrix and normalization

**Inputs.** The $n$ fit points $\mathbf{p}_i=(u_i,v_i)$ and the *Bandwidth* slider
$f \in [0.005, 0.15]$, step 0.005, **default 0.03**.

**The estimator.** Both paths compute the standard multivariate Gaussian KDE with a
**full-covariance, data-adaptive bandwidth matrix** — *not* an isotropic kernel and *not* a
fixed width in Å:

$$\rho(\mathbf{p}) \;=\; \frac{\kappa}{n}\sum_{i=1}^{n}
\frac{1}{2\pi\sqrt{\det \mathbf{H}}}\,
\exp\!\left[-\tfrac{1}{2}\,(\mathbf{p}-\mathbf{p}_i)^{\!\top}\mathbf{H}^{-1}(\mathbf{p}-\mathbf{p}_i)\right],$$

$$\mathbf{H} \;=\; f^{2}\,\mathbf{C}, \qquad
\mathbf{C} \;=\; \frac{1}{n-1}\sum_{i=1}^{n}(\mathbf{p}_i-\bar{\mathbf{p}})(\mathbf{p}_i-\bar{\mathbf{p}})^{\!\top},$$

with $\kappa$ the periodic-image correction defined below. This is exactly SciPy's convention:
`gaussian_kde` with a **scalar** `bw_method` sets `kde.factor = bw` and then
`self.covariance = self._data_covariance * self.factor**2` with
`_data_covariance = atleast_2d(cov(self.dataset, rowvar=1, bias=False, aweights=self.weights))`,
which with the default uniform weights $w_i = 1/n$ is exactly the $n-1$ divisor. `gaussian_kde`
then evaluates $\sum_i w_i \mathcal{N}(\mathbf{p};\mathbf{p}_i,\mathbf{H})$ with those same uniform
weights (SciPy 1.13.1, `scipy/stats/_kde.py`). The JavaScript `makeKernel()` reproduces it term for term:
it forms $\mathbf{C}$ with the same $n-1$ divisor, scales by $f^2$, inverts the $2\times2$ matrix in
closed form, and sets `normalizer = imageFactor / (2π·sqrt(det)·samples.length)`. **The browser
normalization is exact, not an approximation.**

**$\mathbf{C}$ is estimated from the *subsampled* fit points, not from all $N_\mathrm{img}$ slab rows.** Both
runtimes do the subsample first and hand the reduced array to the covariance step (`slab = slab[choice]`
immediately precedes `gaussian_kde(slab.T, bw_method=bw)`; `makeKernel(samples, …)` receives the
output of `sampleWithoutReplacement`). Two consequences: the smoothing width is itself
**seed- and subsample-dependent**, and because Python and JavaScript draw different 6000-point
subsets (Step 5) **they use slightly different $\mathbf{H}$ for the same slab**. The divergence
between the runtimes is therefore not confined to the kernel centres; it is in the kernel width too.
(For $n=6000$ the sampling error on each covariance entry is $O(n^{-1/2})\approx1.3\,\%$, so
$\sigma$ differs by a few tenths of a percent — small, but not zero.)

#### Scott / Silverman: neither is used here

SciPy's *default* would be Scott's rule, $f_{\mathrm{Scott}} = n_{\mathrm{eff}}^{-1/(d+4)} = n^{-1/6}$
for $d=2$; Silverman's is $\left(n(d+2)/4\right)^{-1/(d+4)}$. **The Structure page never uses
either.** The Flask endpoint reads the bandwidth as an unconditional `float(request.args.get("bw", 0.03))`,
so `bw` is always a user-set *constant covariance factor* that replaces the rule-of-thumb. (There is
a `_bw_argument()` helper in `app.py` that *does* accept the strings `"scott"`/`"silverman"`, but it
is wired only to `/api/pca/kde` — the [PCA Ellipsoid](pca-ellipsoid.md) page — not to
`/api/kde/slice`.)

This matters quantitatively. For $n=6000$, $f_{\mathrm{Scott}} = 6000^{-1/6} \approx 0.235$, whereas
the page's default is $f = 0.03$ — about $8\times$ narrower. **The default deliberately
under-smooths relative to the automatic rule**, because the point of the map is to resolve
individual atomic sites rather than to produce a statistically optimal density estimate. The slider
maximum, 0.15, is still below Scott's value.

#### What the bandwidth means physically

Because $\mathbf{H}=f^2\mathbf{C}$ is tied to the *spread of the slab points*, the kernel width is
not an absolute length. Along a principal axis of $\mathbf{C}$ with eigenvalue $\lambda$ the kernel
standard deviation is

$$\sigma \;=\; f\sqrt{\lambda} \quad \text{(fractional units)} .$$

A slab that fills the cell has $\mathbf{C}\approx \tfrac{1}{12}\mathbf{I}$ before augmentation
($\sqrt\lambda\approx0.289$); with the $m=0.1$ margin the padded support widens it to
$\sqrt\lambda \approx 0.35$. So the default $f=0.03$ corresponds to $\sigma \approx 0.009$–$0.010$
in fractional units — roughly **0.09–0.11 Å for a 10.4 Å cell** (the GNSe sample). Three honest
corollaries:

1. The smoothing width **changes when you change the element filter, the slab thickness, or the
   normal**, because all of those change $\mathbf{C}$. The same `bw = 0.03` is a different physical
   width on different slices.
2. Because the periodic images inflate $\mathbf{C}$ slightly, adding them widens the kernel a few
   percent relative to a non-augmented fit. This is a real (small) side effect of Step 3.
3. In a non-orthogonal cell, $\mathbf{H}$ is a covariance in *fractional* space, so a kernel that is
   near-circular in the computation is an ellipse in Å after the affine map to the real cell
   (Step 10). The smoothing is anisotropic in real space for any non-cubic cell.

#### Periodic-image renormalization $\kappa$

`gaussian_kde` divides by *every* fit point, images included, so the raw estimate integrates to 1
over the **whole padded plane** — the cell plus its $m$-wide collar — and therefore to *less* than 1
over the cell itself, with an amplitude that sags toward the faces. Both paths correct it by

$$\kappa \;=\; \frac{N_\mathrm{img}}{N_\mathrm{src}} \;\ge\; 1 \qquad (\text{slab rows} / \text{unique source atoms}),$$

applied as `density *= slab_total / slab_count` in Python (only when `slab_total > slab_count > 0`)
and folded into the JavaScript normalizer as `imageFactor`. Because the subsample is unbiased, using
the pre-subsample $N_\mathrm{img}$ with the post-subsample $n$ in the denominator is consistent: the sampling
fraction cancels in expectation. For a roughly uniform slab the image count per source atom is the
same $(1+2m)^2$ factor by which the padded in-plane support exceeds the cell, which is exactly what
has to be undone.

**The resulting units.** $\rho$ is a probability density **per unit fractional area of the slice
plane**, normalized so that $\int_{\mathrm{cell}}\rho\,\mathrm{d}u\,\mathrm{d}v \approx 1$: the
**cell** carries unit mass, and the amplitude no longer decays toward the faces. It is **not** one
unit of mass per slab atom (that would integrate to $N_\mathrm{src}$), **not** atoms Å⁻², **not** atoms Å⁻³, and
it is **not divided by the slab thickness** — thickening the slab pulls in more atoms but the field
is renormalized, so absolute values are not comparable between different $\Delta z$, different
elements, or different bandwidths. The per-atom amplitude is $1/N_\mathrm{src}$ of the field.

#### Degenerate-slab handling (the two paths differ)

| Condition | `kde.py` | `localKdeWorker.js` |
| --- | --- | --- |
| slab rows $< 5$ | zero grid, `fitCount = 0` | zero grid, `fitCount = 0` |
| $<3$ distinct $(u,v)$ pairs | zero grid (`has_enough_unique_points`) | no such check |
| points collinear (`matrix_rank < 2`) | zero grid (`has_two_dimensional_spread`) | ridge-regularized and drawn anyway |
| singular / near-singular $\mathbf{H}$ | `suppress(LinAlgError, ValueError)` → zero grid | `k00,k11 += 1e-8`; if $\det \le 10^{-12}$, inflate both diagonals by $\max(c_{00},c_{11},10^{-4})f^2+10^{-6}$ and zero the cross term |

So the browser path adds a $10^{-8}$ ridge to the diagonal of $\mathbf{H}$ **always**, not only in
the degenerate case. For the default $f=0.03$ and a cell-filling slab, $H_{00}\approx7.5\times10^{-5}$,
so the ridge is a $\sim1.3\times10^{-4}$ relative perturbation — negligible but non-zero, and it means
the browser kernel is never *exactly* $f^2\mathbf{C}$.

**Parameter coercion: a silent substitution, not a floor.** The browser writes
`const factor = Math.max(Number(bandwidth) || 0.03, 1e-4)`. The `|| 0.03` short-circuit fires for
`0`, `NaN`, `null` and any non-numeric value, so **`bw = 0` becomes 0.03, not $10^{-4}$**; the
$10^{-4}$ floor is only reachable for values in $(0, 10^{-4})$. The same idiom appears twice more
and inconsistently: the periodic margin in `computeKde()` uses the **un-floored**
`2 * (Number(bandwidth) || 0.03)`, so for a sub-$10^{-4}$ bandwidth the margin and the kernel
disagree about $f$; and the grid uses `Number(gridSize) || 120`, so `gridSize = 0` becomes 120, not
the 16 lower clamp. **Python applies none of this.** A request with `bw=0` reaches `gaussian_kde`,
raises inside the `with suppress(np.linalg.LinAlgError, ValueError)` block, and returns an
**all-zero density grid with `fitCount = 0` and HTTP 200** — verified directly: `bw=0` on 200
coplanar points gives `vmin = vmax = 0`, `slabCount = 200`, `fitCount = 0`, zero contours, and no
error reaches the client.

`tests/test_kde.py::test_kde_slice_handles_degenerate_slab_without_error` pins the Python behaviour:
five perfectly collinear points give `slabCount = 5`, `fitCount = 0`, `vmin = vmax = 0`.

**Code.** `rmc_toolkits/kde.py` → `kde_slice()` (the `gaussian_kde(slab.T, bw_method=bw)` call and
the `density *= slab_total / slab_count` rescale); `localKdeWorker.js` → `covariance()`,
`makeKernel()`, and the `margin`/`grid` coercions in `computeKde()`.

---

### Step 7 — Evaluate on the grid (CPU loop and WGSL compute shader)

**Inputs.** The kernel, the fit points, and the *Grid* selector (options **80 / 120 / 160 / 220**,
default **120**).

**The grid.** A uniform $G\times G$ lattice spanning the projected-cube extent, endpoints inclusive:

$$u_p = x_{\min} + p\,\frac{x_{\max}-x_{\min}}{G-1},\qquad
v_q = y_{\min} + q\,\frac{y_{\max}-y_{\min}}{G-1},\qquad p,q = 0..G-1 ,$$

stored row-major as `density[q][p]` (row index = $v$, column index = $u$). Python builds it with
`np.linspace` + `np.meshgrid`; JavaScript with explicit `xStep`/`yStep` — identical node positions.

#### The grid is the bounding box of the projected cube, not the cell cross-section

$[x_{\min},x_{\max}]\times[y_{\min},y_{\max}]$ comes from the projections of **all eight cube
corners** onto $(\hat{\mathbf{u}},\hat{\mathbf{v}})$ (Step 2). For the `a`/`b`/`c` presets the
cross-section *is* that rectangle (the unit square) and nothing is wasted. For any other normal the
cross-section is a smaller triangle/quadrilateral/pentagon/hexagon **inside** the rectangle, so a
substantial fraction of the $G^2$ nodes lie **outside the cell**: for $\mathbf{h}=(1,1,1)$ at
mid-depth the section is a regular hexagon of area $1.299$ inside a $1.633\times1.414 = 2.309$
rectangle, i.e. **≈44 % of the nodes are outside**.

Those nodes are clipped away at draw time (`ctx.clip()` against `planePolygon`, Step 10) and are
never shown — **but they are included in `vmin`/`vmax`**, and therefore in both the colour
normalization (Step 10) and the eight contour levels (Step 9). In log mode they pull `vmin` down to
the $10^{-12}$ floor, i.e. $-12$, which stretches the colour scale over a range the user cannot see
and compresses the visible contrast. Measured on `np.random.default_rng(0).random((3000, 3))` through
`oriented_kde_slice(pts, center=0.5, thickness=0.08, normal=h, bw=0.03, grid=120, log=True)`
(numpy 2.0.2 / scipy 1.13.1): `vmin = −11.55`, `vmax = +0.91` for $\mathbf{h}=(1,1,0)$ and
`vmin = −12.00` (the log floor exactly), `vmax = +0.68` for $\mathbf{h}=(1,1,1)$.

**Grid clamping differs**: Python `grid = int(max(16, min(grid, 400)))`, JavaScript
`Math.max(16, Math.min(Number(gridSize) || 120, 260))`. Both admit every value the UI offers, but a
hand-crafted API request for `grid=300` is honoured by Flask and would be clipped to 260 by the
browser. `tests/test_kde.py::test_kde_slice_clamps_grid_and_empty_slab` pins the lower clamp
(`grid=2` → 16).

**Reference path (Flask).** One `gaussian_kde.__call__` on the ravelled $G^2$ sample points, which
dispatches to SciPy's compiled `gaussian_kernel_estimate` using the Cholesky factor of $\mathbf{H}$.
Full float64. No distance cutoff — every kernel contributes to every node.

**CPU path (browser).** `computeDensityCpu()` is the direct $O(G^2 n)$ triple loop:

```js
const exponent = -0.5*(inv00*dx*dx + 2*inv01*dx*dy + inv11*dy*dy);
if (exponent > -60) sum += Math.exp(exponent);
density[y][x] = sum * kernel.normalizer;
```

Note the **exponent cutoff at $-60$**: terms with $e^{-60}\approx 8.8\times10^{-27}$ are dropped,
i.e. the kernel is truncated at a Mahalanobis radius $\sqrt{120}\approx11\sigma$. SciPy applies no
such cutoff. The truncation is far below float64 resolution relative to the peak, and it is usually
looser than the periodic-image truncation of Step 3, so it is not a practical difference — but it is
a real difference in the formula evaluated.

**GPU path (browser).** `gpuKde.js` runs the identical sum as a WGSL compute shader, one invocation
per grid cell:

* Workgroup size `(8, 8)`; dispatch $\lceil G/8\rceil$ workgroups in each dimension with an
  in-shader bounds guard `if (gid.x >= grid || gid.y >= grid) { return; }`.
* Bindings: a 48-byte uniform buffer packed as three `vec4` lanes — `(inv00, inv01, inv11,
  normalizer)`, `(xMin, yMin, xStep, yStep)`, `(grid, sampleCount, pad, pad)` read through a
  `Uint32Array` view of the same `ArrayBuffer`; a read-only storage buffer of tightly packed
  `vec2<f32>` samples; a read-write storage buffer of $G^2$ `f32` outputs, copied to a
  `MAP_READ` buffer and read back with `mapAsync`.
* The shader body is line-for-line the CPU expression, including the `e > -60.0` guard.

**When the GPU is used.** Only when the work is big enough to amortize device setup, buffer
uploads, and the readback round-trip:

$$\texttt{shouldUseGpu} \iff G \times G \times n \;\ge\; \mathbf{2{,}000{,}000} \quad (\texttt{GPU\_MIN\_WORK}).$$

At the default $G=120$ that threshold is crossed at $n \ge 139$ fit points, so a normal slab
(thousands of atoms) uses the GPU when WebGPU is present.

**Fallback guarantee, stated precisely.** The device/pipeline promise is created once per worker
and cached; *any* failure — no `navigator.gpu`, `requestAdapter()` returning null, `requestDevice()`
or pipeline creation rejecting, a runtime error inside `computeDensityGpu`, or sub-threshold work —
resolves to `null` and `computeKde()` falls through to `computeDensityCpu()` (`density = mapped ??
computeDensityCpu(args)`). A lost device clears the cached promise so a later message can
re-initialize. The result object reports which one ran via `backend: 'gpu' | 'cpu'`.

The repo's own wording — `AGENTS.md`: *"fall back to the CPU loop with identical output"*;
`gpuKde.js`: *"devices without WebGPU behave exactly as before"* — is true **structurally** (same
formula, same grid, same normalizer, same cutoff, reshaped to the same nested JS array) but not
**bitwise**, and the float32 narrowing is broader than the accumulator alone. Everything crosses the
boundary as `f32`:

* the sample coordinates (`new Float32Array(sampleCount * 2)`);
* all four kernel parameters and the grid geometry (`paramFloats[0..7]` = `inv00, inv01, inv11,
  normalizer, xMin, yMin, xStep, yStep`);
* the **node positions themselves**, which the shader reconstructs as
  `xMin + f32(gid.x) * xStep` rather than receiving them from the JS loop, so the grid coordinates
  differ in the last bits too;
* the accumulator (`var sum : f32`), and WGSL `exp` need not match `Math.exp` to the last bit.

For a sum of up to 6000 positive float32 terms the expected relative difference is of order
$10^{-6}$–$10^{-5}$ — invisible in an 8-bit colormap, but it should not be described as identical
arithmetic.

**Outputs.** `density[G][G]` (nested plain JS numbers / Python list of lists), plus the CPU/GPU
backend flag in static mode.

**Code.** `localKdeWorker.js` → `computeDensityCpu()`, `computeKde()`; `gpuKde.js` → `KDE_WGSL`,
`GPU_MIN_WORK`, `shouldUseGpu()`, `getGpu()`, `computeDensityGpu()`; `rmc_toolkits/kde.py` →
`kde_slice()`.

---

### Step 8 — Optional $\log_{10}$ compression

**Inputs.** The *Log scale* switch — **on by default** (`useState(true)`).

**Math.** When enabled, the whole grid is replaced in place by

$$\tilde\rho = \log_{10}\!\left(\rho + 10^{-12}\right),$$

the $10^{-12}$ floor preventing $\log(0)$ where the kernel sum underflows. Empty regions therefore
sit at exactly $-12$. Identical constant and identical formula in both runtimes
(`np.log10(density + 1e-12)` / `Math.log10(density[y][x] + 1e-12)`).

`vmin`/`vmax` are computed **after** the log transform in both paths, so the colormap and the
contour levels operate on log density when the switch is on. The canvas prints `log10 density` as
an overlay label in that case.

**Code.** `rmc_toolkits/kde.py` → `kde_slice()` (`if log:` branch); `localKdeWorker.js` →
`computeKde()` min/max loop.

---

### Step 9 — Contour extraction

**Level selection (identical in both runtimes).** Levels are **evenly spaced fractions of the
observed range** — *not* quantiles, and not fractions of the max:

$$\ell_k \;=\; v_{\min} + \frac{k}{K+1}\,(v_{\max}-v_{\min}), \qquad k = 1..K, \quad K = 8 .$$

Python writes it as `np.linspace(finite_min, finite_max, n_levels+2)[1:-1]`; JavaScript as
`vmin + (levelIndex/(levels+1))*(vmax - vmin)`. These are the same $K$ interior levels. The
frontend always requests `levels: 8`, and 8 is also the default in both functions. Because the
levels are anchored to $v_{\min}$/$v_{\max}$ of *this slice*, **contour levels are not comparable
between slices**, and in log mode they are equally spaced in $\log_{10}\rho$ (i.e. geometrically
spaced in $\rho$). As Step 7 notes, $v_{\min}$/$v_{\max}$ are taken over **all $G^2$ nodes including
those clipped away outside the cell cross-section**, so for an oblique normal the eight levels are
spread over a range wider than anything the user can see.

**Tracing — this is where the two paths genuinely differ.**

* **Flask / reference** — `contourpy.contour_generator(grid_x, grid_y, density)` (contourpy ships
  with matplotlib; used directly to avoid pyplot global state). This is contourpy's `serial`
  marching-squares generator with default settings, which returns `LineType.Separate`: a list of
  $(m,2)$ float arrays, each a **stitched polyline** with proper saddle-cell disambiguation.
  Polylines with fewer than 2 points are dropped. Emitted as
  `{"level": …, "lines": [[[x,y], …], …]}`.
* **Static / browser** — `extractContours()` is a hand-written marching-squares pass over each
  $2\times2$ cell of the grid. The four corners are visited in the cycle
  **0 → 1 → 2 → 3 → 0** = (lower-left, lower-right, upper-right, upper-left), and an edge
  $(a,b)$ is declared crossed by the **half-open** test

  ```js
  (a.value < level && b.value >= level) || (b.value < level && a.value >= level)
  ```

  i.e. a corner exactly *at* the level counts as "above". The crossing point is placed by **linear
  interpolation**, `t = (level − a.value)/(b.value − a.value)` (with `t = 0.5` when the two corner
  values differ by $\le 10^{-12}$). Because this counts sign changes around a **closed** cycle the
  count is always even — 0, 2 or 4; 1 and 3 are impossible. Two crossings → one segment. Four
  crossings (a saddle) → **arbitrarily** paired as `[e0,e1]` and `[e2,e3]` with no disambiguation by
  the cell-centre value. Segments are emitted as independent 2-point polylines — they are never
  stitched into curves.

**A level that produces nothing is dropped from the output.** Python appends only `if polylines:`,
JavaScript only `if (lines.length)`. So `contours.length < K` is normal and is not an error; an API
consumer must not assume eight entries, and must read `contour.level` rather than inferring it from
the index.

Both produce the same *set* of crossing points; the browser version can connect a saddle cell the
"wrong" way and produces many short strokes instead of continuous curves. On screen this is
indistinguishable at typical grid sizes (segments are ~1 px), but the browser contour data is not
suitable for extracting a closed iso-line.

**A behavioural discrepancy in log mode.** The Python guard is

```python
if n_levels <= 0 or not np.isfinite(density).any() or float(density.max()) <= 0:
    return []
```

That `density.max() <= 0` test is applied **after** the log transform. So whenever the peak
$\log_{10}\rho \le 0$ — i.e. peak linear density $\le 1$ per unit fractional area — the reference
path emits **zero contours**, while the browser path (whose only guard is `vmax > vmin`) still
draws all 8. Verified directly against the code with a fixed fixture — 200 points from
`np.random.default_rng(0)` at $z=0.5$, `kde_slice(..., z_center=0.5, dz=0.1, xlim=ylim=(0,1),
grid=120, log=True)`:

| bandwidth | peak $\log_{10}\rho$ | contours, `log=True` | contours, `log=False` |
| --- | --- | --- | --- |
| `bw=0.5` | $+0.027$ | 8 | 8 |
| `bw=2.0` | $-0.443$ | **0** | 8 |

(The peak is stable across `grid` 16–220; `oriented_kde_slice` on the same points gives $+0.021$ and
$-0.421$, the small offset coming from its periodic augmentation.) With the slider's maximum
bandwidth of 0.15 and a normally populated slab the peak is comfortably $>1$, so this rarely bites —
but it is a real mode-dependent difference, not a rendering choice.

**Rendering.** `drawKdeSlice()` strokes every polyline in data coordinates through the same
plane mapper used for the cell outline, at `lineWidth = 1`, in a theme-dependent colour
(`rgba(21,34,50,0.72)` in light theme, `rgba(230,236,244,0.76)` in dark). The *Contours* switch is a
pure client-side toggle — turning it off does not recompute anything.

**Code.** `rmc_toolkits/kde.py` → `_contour_segments()`; `localKdeWorker.js` → `extractContours()`;
`StructurePage.jsx` → `drawKdeSlice()` contour loop.

---

### Step 10 — Density → colour, and the affine map onto the real cell

**Normalization: per-slice min–max, linear, no shared scale.** `drawKdeSlice()` computes

$$\hat\rho_{pq} = \frac{\rho_{pq} - v_{\min}}{v_{\max}-v_{\min}}, \qquad
\text{LUT index} = \operatorname{clamp}\big(\operatorname{round}(255\,\hat\rho_{pq}),\,0,\,255\big),$$

where $v_{\min},v_{\max}$ are the min and max **of this slice only**, returned by the engine
(post-log if log scale is on, and computed over the whole bounding-box grid — Step 7). Consequences
the reader must know:

* The mapping is **linear in whatever quantity is in the grid** — linear in $\rho$ with the log
  switch off, linear in $\log_{10}\rho$ with it on (the default).
* The scale is **re-normalized on every recompute**. Moving the slider, changing the element,
  changing the bandwidth, or changing the grid all rescale the colours. **Colours are not comparable
  between two screenshots.** There is no colorbar and no numeric legend anywhere on the panel — only
  the text overlay giving `slabCount`, `fitCount`, $z_c$, $\Delta z$ and `bw`.

#### `"No atoms in this slab"` really means *no drawable density*

The draw gate is `density && grid > 0 && kde.vmax > kde.vmin`; when it fails the canvas prints
`"No atoms in this slab"` (or `"Computing KDE..."` while a request is in flight). But **every**
degenerate bail-out of Step 6 returns an all-zero grid, hence $v_{\min}=v_{\max}$, hence the same
message: fewer than 5 slab rows, fewer than 3 distinct $(u,v)$ pairs, a rank-deficient spread, or a
suppressed `LinAlgError` (including `bw = 0`). In those cases the overlay simultaneously prints
`"{slabCount} atoms in slab (fit 0)"` with a non-zero `slabCount` — the repo's own fixture
`test_kde_slice_handles_degenerate_slab_without_error` produces exactly that state
(`slabCount = 5`, `fitCount = 0`, `vmin = vmax = 0`). The two statements on screen contradict each
other and there is no indication that the estimator declined the slab; **read `fitCount = 0` with a
non-zero `slabCount` as "the estimator bailed", not "the slab is empty".**

**The colormaps are 5-anchor approximations.** [`colormaps.js`](../../web_app/frontend/src/colormaps.js)
defines five maps — `viridis`, `magma`, `seismic`, `reds`, `greys` (default **viridis**) — each as a
list of **five RGB anchors**, expanded by piecewise-**linear** interpolation into a 256-entry
`Uint8ClampedArray` LUT that is cached per name (`getLut`). The parameterization, for reproducibility:
for entry $i$, $t = i/255$, `scaled = t·4`, `lower = min(4, floor(scaled))`, `upper = min(4, lower+1)`,
`frac = scaled − lower`. The four segments are therefore $255/4 = 63.75$ entries wide and the last
anchor is hit exactly only at $i = 255$. The interpolated float is written straight into a
`Uint8ClampedArray`, which **rounds half-to-even** rather than truncating. (`sampleColormap()` is a
separate exported helper with its own clamp; this page does not use it — `StructurePage.jsx` imports
only `getLut`.)

These are *approximations of* the matplotlib maps of the same name, not the real 256-entry tables:
`viridis` is anchored at `(68,1,84) → (59,82,139) → (33,145,140) → (94,201,98) → (253,231,37)`. They
are close enough to read but are **not** perceptually uniform in the way the true matplotlib LUTs
are, and `greys` here runs dark → light (the opposite sense to matplotlib's `Greys`). Do not use a
screenshot of this panel to read off matplotlib-calibrated colour values.

**Geometry: from fractional grid to the real (possibly oblique) cell.** The density was computed in
fractional space; the drawing puts it on the true cell parallelogram:

1. $\hat{\mathbf{u}},\hat{\mathbf{v}}$ (fractional) are mapped through the unit-cell basis:
   $\mathbf{u}_{\text{Å}} = u_1\mathbf{a}+u_2\mathbf{b}+u_3\mathbf{c}$, likewise $\mathbf{v}_{\text{Å}}$
   (`vectorFromFraction(uVector, unitCell.unitVectors)`).
2. `makePlane()` embeds those two Å vectors in 2-D preserving both lengths and the angle between
   them: $\mathbf{u}\mapsto(\lVert \mathbf{u}\rVert,0)$,
   $\mathbf{v}\mapsto(\lVert\mathbf{v}\rVert\cos\theta,\ \lVert\mathbf{v}\rVert\sin\theta)$ with
   $\cos\theta = \mathbf{u}\!\cdot\!\mathbf{v}/(\lVert\mathbf{u}\rVert\lVert\mathbf{v}\rVert)$
   clamped to $[-1,1]$. **So the drawn parallelogram has the correct real-space aspect ratio and
   shear.**
3. `makePlaneMapper()` fits that parallelogram into the canvas with an 18 px padding and uniform
   isotropic scaling (`Math.min` of the two fit factors), centred. It also exposes `invert()`, used
   by the Slab-In-Cell drag handler.
4. The $G\times G$ density is written into an offscreen `ImageData` of exactly $G\times G$ pixels,
   then blitted with `ctx.transform(...)` built from the images of $(x_{\min},y_{\min})$,
   $(x_{\max},y_{\min})$ and $(x_{\min},y_{\max})$, drawn into the unit square, with
   `imageSmoothingEnabled = true` (browser bilinear interpolation) and clipped to the plane-section
   polygon (`planePolygon`, from `_plane_section_vertices()` / `planeSectionVertices()`, Step 4).

**Where the real-space basis comes from — and how it fails.** `unitCell` is a `useMemo` in
`StructurePage.jsx` that computes `unitVectors[j] = structure.latticeVectors[j] /
max(structure.supercell[j], 1e-12)` from the **structure** payload. The `unitVectors` and
`cellLengths` fields that `/api/kde/slice` returns are **never read** by the frontend. If
`structure.latticeVectors` or `structure.supercell` is missing, the memo silently falls back to the
identity basis `[[1,0,0],[0,1,0],[0,0,1]]` with lengths $(1,1,1)$ — the panels then draw a 1 Å cubic
cell with the wrong aspect ratio and shear, and **no warning is shown**.

**`--panel-aspect` is not the aspect of the drawn parallelogram.** The CSS custom property comes
from `slicePanelGeometry`, which runs `makeProjectedPlane` over the projections of **all eight cube
corners**; `drawKdeSlice` builds its *own* `makeProjectedPlane` over the **plane-section polygon**.
For an oblique normal those two bounding boxes differ, so the panel's CSS box is sized from the
cube's bounding rectangle while the drawing inside it is fitted to the cross-section. The three
panels use: KDE panel → `planeAspect`; Slab In Cell → `sideAspect` (the
$(\hat{\mathbf{u}},\hat{\mathbf{h}})$ corner projection); Folded Unit Cell → `max(planeAspect, 1)`.

**Numerical guards in the plane/mapper code.** Each can silently change the output, so they are
listed here rather than left implicit:

| Guard | Location | Effect |
| --- | --- | --- |
| $\max(\lVert\mathbf{u}\rVert\lVert\mathbf{v}\rVert,\ 10^{-12})$ | `makePlane` | avoids 0/0 in $\cos\theta$ for a degenerate basis |
| $\cos\theta$ clamped to $[-1,1]$, then $\sin\theta=\sqrt{\max(0,1-\cos^2\theta)}$ | `makePlane` | keeps $\sin\theta$ real under rounding |
| $\max(\mathrm{y-span},\ 10^{-9})$ in the aspect ratio | `makePlane`, `makeProjectedPlane` | a flat plane reports a huge, finite aspect instead of `Infinity` |
| zero span replaced by 1 | `makePlaneMapper` | a degenerate box still maps |
| $\lvert\det\rvert < 10^{-12}$ → `invert()` returns `{uFraction: 0, vFraction: 0}` | `makePlaneMapper` | **the drag handler silently pins to the plane origin instead of erroring** |
| `kde === null` → extent $[-0.5,0.5,-0.5,0.5]$, $\hat{\mathbf{u}},\hat{\mathbf{v}}$ from `sliceConfig`, `planePolygon` = the extent rectangle | `drawKdeSlice` | the empty panel still draws a plausible outline |

**A half-cell registration offset (rendering only).** The density samples sit at grid *nodes*
$p/(G-1)$, but the blit places them as image *pixels* whose centres sit at $(p+0.5)/G$. The two
sequences differ by at most $0.5/G$ of the plane extent (0.42 % at $G=120$, ~0.04 Å for a 10.4 Å
cell), sliding from $+0.5/G$ at one edge to $-0.5/G$ at the other. The **contours are drawn in data
coordinates and are exact**, so at high zoom the contour lines can sit up to half a grid cell away
from the colour feature they enclose. This affects the picture only — the returned `density` array
is unaffected.

**Canvas resolution.** The KDE canvas is measured with `getBoundingClientRect()` and sized in CSS
units to `max(320, floor(width)) × max(260, floor(height))` (the Slab In Cell canvas uses
`220 × 260`). The backing store is that size multiplied by `window.devicePixelRatio || 1`, with a
matching `ctx.setTransform(dpr, 0, 0, dpr, 0, 0)` so the draw code always works in CSS pixels. A
`ResizeObserver` on both canvases coalesces resize events through a single `requestAnimationFrame`
and bumps a `sizeTick` state, which forces a re-measure and redraw — this is what stops a
first-visit panel measured before layout settles from staying at its minimum size.

**Export.** The *Save* menu (`PANEL_SAVE_OPTIONS`) offers two entries, and only one re-renders:

* **`png` — "PNG image (1×)"** takes the live canvas verbatim via `saveCanvasAsPng(canvas, name)`,
  so it captures the current backing store at the current `devicePixelRatio`.
* **`png3x` — "High-res PNG (3×)"** re-runs the same `drawKdeSlice()` into an offscreen canvas at
  `scale = 3` with `ctx.setTransform(3,0,0,3,0,0)`, so the high-res PNG is a genuine re-render, not
  an upscaled bitmap.

`save2dPanel()` re-measures with the **same minimum floors** (320×260 / 220×260), so a panel
displayed smaller than its minimum exports at 3× the *minimum*, not 3× its displayed size. In both
cases the density grid itself is **not** re-computed at higher resolution — the $G\times G$ array is
whatever the last request returned.

**Code.** `StructurePage.jsx` → `drawKdeSlice()`, `makePlane()`, `makeProjectedPlane()`,
`makePlaneMapper()`, `save2dPanel()`, the `unitCell` and `slicePanelGeometry` memos, and the
canvas-sizing / `ResizeObserver` effects; `colormaps.js` → `ANCHORS`, `buildLut()`, `getLut()`.

---

### Step 11 — Companion panels (briefly)

The other two panels share the same slice state but do **no** density estimation.

* **Slab In Cell** — a parallel (generally **oblique/axonometric**) side view in the
  $(\hat{\mathbf{u}}, \hat{\mathbf{h}})$ plane, again mapped through the real cell vectors. It drops
  the $v$ component, i.e. projects along $\mathbf{e}_v$, which is orthographic only when the cell
  metric makes $\mathbf{e}_v \perp \mathrm{span}(\mathbf{e}_u,\mathbf{e}_h)$; see §4d of the second
  half. Every displayed atom is plotted as a 1–2 px rectangle; atoms inside the slab get their
  element colour and 2 px, atoms outside get
  `rgba(166,176,188,0.22)` and 1 px. Points are strided at
  `stride = max(1, floor(points.length / min(points.length, 1_000_000)))`, i.e. no additional
  thinning below one million atoms. The blue band is draggable: `makePlaneMapper().invert()` maps
  the cursor back to plane coordinates, and the drag sets $z_c$ live (clamped to $[0,1]$).
* **Folded Unit Cell** — a Three.js `THREE.Points` cloud, one draw call per element, positions
  $\mathbf{x}_i$ mapped through a *normalized* basis (each unit-cell vector divided by the longest
  one, so the model is unit-scaled but keeps its shape) and centred. The slab is drawn as **two
  translucent plane-section caps** — one at the slab's front depth, one at its back — each
  triangulated as an independent fan by `makeSlabGeometry()`. There are **no side-wall triangles**:
  the only thing joining the two caps is a set of `LineSegments` from `makeSectionEdgeGeometry()`,
  and those cap-to-cap lines are emitted **only when the two sections have equal vertex counts**. It
  is not a closed solid band.

#### The drawn band and the density disagree at the cell faces

Both the 2-D side view and the 3-D band clamp the slab to the cube:

$$\mathrm{depthStart} = \max\big(d_{\min},\, d_{\mathrm{c}} - \tfrac{\Delta z \Delta_d}{2}\big), \qquad
\mathrm{depthEnd} = \min\big(d_{\max},\, d_{\mathrm{c}} + \tfrac{\Delta z \Delta_d}{2}\big).$$

The KDE's depth mask (Step 4) is **not** clamped and deliberately picks up wrapped periodic images,
and the atom highlighting in the side view is decided by `inActiveSlab()`, which tests only the
*unwrapped* folded depth. So for $z_c$ near 0 or 1 the picture **understates the selection**: the
drawn band is narrower than the depth range actually sampled, and atoms that contribute to the
density (and to `slabCount`) through their images are drawn in the grey "outside" colour. Python's
`slabVertices` is computed from the same clamped faces — and is not consumed by the frontend at all,
which recomputes the band itself.

#### 3-D panel parameters (display-only, but not reproducible without them)

| Item | Value |
| --- | --- |
| Point sprite | `THREE.PointsMaterial`, `size = 0.018`, `sizeAttenuation: true` (normalized-basis units) |
| Renderer | `antialias: true`, `preserveDrawingBuffer: true`, `pixelRatio = min(devicePixelRatio, 2)` |
| Camera | `PerspectiveCamera(fov 45)`; `near = r/100`, `far = 20r` re-derived after the bounds fit |
| $r$ | `max(Box3(cell corners).getBoundingSphere().radius, 0.5)` — the sphere of the **AABB** of the corners, not of the corner set (§7.5) |
| Controls | `OrbitControls`, `enableDamping`, `dampingFactor = 0.08`, pan on, `minDistance = 0.35r`, `maxDistance = 8r` |
| Initial camera | `sphere.center + (1.7r, 1.45r, 1.55r)`, looking at `sphere.center` |
| Slab material | `#4f8cff`, `opacity 0.12`, `DoubleSide`, `depthWrite: false` |
| Cell / slab edges | `#737c86` and `#8c96a3` (opacity 0.95) `LineSegments` |
| Camera persistence | position/target/zoom saved to `cameraStateRef` on unmount and restored, so slider changes do not reset the view |

Both panels use the same per-element palette (`atomColors.js` → `buildElementColors`, which sorts
the distinct element labels before assigning colours) shown in the legend **below** the 3-D view
(the `atom-legend` block is rendered after the Three.js mount element and is styled with
`border-top`).

---

### Request lifecycle, caching and determinism

* **Debounce.** Slider changes are debounced before a recompute: **160 ms** on the Flask path (with
  an `AbortController` cancelling any in-flight request) and **80 ms** on the browser-worker path.
  Worker results are matched to a monotonically increasing request id so a stale reply is discarded.
* **Backend cache.** `_cached_positions()` is an `lru_cache(maxsize=16)` keyed on
  `(path, mtime, element)`, so re-slicing the same file does not re-parse it. The KDE itself is
  recomputed per request.
* **Determinism.** Given the same file, element, normal, $z_c$, $\Delta z$, $f$, grid and log flag,
  each runtime returns a bit-reproducible result on the CPU path (the subsample seed is fixed at 0
  in both). The GPU path is deterministic per device but float32.
* **Error handling differs, and one path deliberately keeps stale pixels on screen.**
  * *Browser worker.* `self.onmessage` wraps `computeKde` in `try/catch` and posts
    `{ id, error: error.message || 'Browser KDE computation failed' }`; a worker-level failure
    (`worker.onerror`) sets `'Browser KDE worker failed'`. Both set `kde = null`, which makes the
    canvas fall back to the $[-0.5,0.5]$ extent and print `"No atoms in this slab"`.
  * *Flask.* A failed request surfaces `err.response?.data?.error || 'KDE computation failed'` and
    clears `kde`. **Cancellations are swallowed** (`axios.isCancel(err) || err.code ===
    'ERR_CANCELED'` skips both `setKde(null)` and `setKdeError`), so an aborted request — every
    debounced slider move — intentionally leaves the *previous* density visible rather than blanking
    the panel.

#### The two JSON payloads are not the same shape

| Key | SciPy path | Browser worker | Notes |
| --- | --- | --- | --- |
| `density`, `extent`, `grid`, `bw`, `log`, `slabCount`, `fitCount`, `vmin`, `vmax`, `contours` | ✓ | ✓ | same meaning |
| `center`, `thickness` | ✓ | ✓ | the raw slider fractions in **both** — this is what the UI reads |
| `normal`, `uVector`, `vVector`, `planeVertices`, `planePolygon` | ✓ | ✓ | `uVector`/`vVector` differ for custom normals (Step 2) |
| `z`, `dz` | `center_depth`, `thickness_depth` | `zCenter`, `thickness` | **different meanings** — see below |
| `depth`, `depthThickness`, `depthRange` | ✓ | — | absolute depth-projection units |
| `slabVertices` | ✓ | — | the two clamped slab faces; **not consumed by `StructurePage.jsx`** |
| `cellLengths`, `unitVectors`, `orientation`, `source`, `element` | ✓ (added by the endpoint) | — | `cellLengths`/`unitVectors` are ignored by the frontend (Step 10) |
| `browserKde: true`, `backend: 'gpu' \| 'cpu'` | — | ✓ | worker-only |

**`z`/`dz` mean different things.** Flask returns `z = center_depth = ` $d_{\min} + z_c\Delta_d$ and
`dz = thickness_depth = ` $\Delta z\,\Delta_d$. The `dz` form is a clean scaling, but `z` carries the
$d_{\min}$ offset, which vanishes only when the normal has no negative components — so it reduces to
$z_c\Delta_d$ for any all-positive normal and to plain $z_c$ for the a/b/c presets, while for e.g.
$\mathbf{h}=(1,-1,0)$ it is shifted by $d_{\min} = -1/\sqrt2$. The browser worker returns the raw
slider fractions. The UI reads `center`/`thickness`, which *are* the slider fractions in both
runtimes, so nothing on screen is wrong — but an API consumer reading `z`/`dz` must know which
runtime produced the payload.

---

### What the test suite actually checks (and what it does not)

The Python tests in `tests/test_kde.py` are:

| Test | What it pins |
| --- | --- |
| `test_load_unit_cell_positions_filters_by_element` | GNSe sample: 52 000 atoms total, 4 000 Ga; `cell_lengths = (10.4116, 10.4116, 10.4116)` Å; folded Cartesian positions inside the cell (skipped when the gitignored sample is absent) |
| `test_load_unit_cell_positions_preserves_nonorthogonal_basis` | a triclinic 1×1×1 fixture round-trips `unit_vectors`, `fractional_positions`, `positions` |
| `test_kde_slice_returns_density_grid_and_contours` | grid/extent/contour shape of a normal call |
| `test_kde_slice_clamps_grid_and_empty_slab` | `grid=2` → 16, and an empty slab |
| `test_kde_slice_handles_nonempty_structure_with_empty_slab` | a populated structure whose slab catches nothing |
| `test_kde_slice_handles_degenerate_slab_without_error` | 5 collinear points → `slabCount 5`, `fitCount 0`, `vmin = vmax = 0` |
| `test_kde_slice_reports_subsampled_fit_count` | 6010 → `fitCount 6000` |
| `test_oriented_kde_slice_wraps_density_across_the_cell_boundary` | edge agreement to 20 %, and $>3\times$ the non-periodic reference |
| `test_oriented_kde_slice_wraps_slab_selection_in_depth` | 30 atoms at $z\approx0.98$ found by a slab at $z_c=0$ |
| `test_oriented_kde_slice_supports_axis_and_custom_normals` | preset and $(110)$ normals produce a unit normal, a non-empty `planePolygon`/`planeVertices`, and the expected `slabCount` |

`web_app/frontend/src/workers/__tests__/localKdeWorker.test.js` covers the image-margin count, the
cell-boundary wrap and the depth wrap for the worker (with the different third assertion noted in
Step 3).

**Two gaps the parity table below cannot paper over.**

1. **No test compares the Python and JavaScript density fields numerically.** The two suites run in
   different languages on different fixtures; nothing cross-checks a grid value, a `vmin`, or a
   contour coordinate between runtimes. Every "Exact" in the parity table is derived by **reading the
   code**, not by measurement.
2. **The GPU path has no automated coverage at all.** No test file imports `gpuKde.js`, and WebGPU is
   absent under vitest, so `computeDensityGpu`, the WGSL shader and the fallback chain are exercised
   only by hand in a real browser.

---

### Parameters and defaults

| Parameter | UI control | Default | Range / options | Units | Where enforced |
| --- | --- | --- | --- | --- | --- |
| Element | select | `all` | elements found in the file | — | Flask: `load_unit_cell_positions(element=)` server-side; browser: the `points` `useMemo` client-side |
| Normal $\mathbf{h}$ | menu | `c` = $(0,0,1)$ | `a`, `b`, `c`, `Custom` | Miller indices $(hkl)$, dimensionless | `SLICE_PRESETS` / `SLICE_ORIENTATIONS` |
| Custom direction | 3 number inputs | `[1, 1, 0]` | any (step 0.1) | Miller indices, dimensionless | `StructurePage.jsx` → `customDirection`; zero vector → silent $(0,0,1)$ |
| Slice centre $z_c$ | range slider | auto-set to densest of 50 depth bins; state default 0.5 | 0 – 1, step 0.001 | fraction of depth span $\Delta_d$ | slider; Python clamps to $[0,1]$ |
| Thickness $\Delta z$ | range slider | **0.08** | 0.01 – 0.5, step 0.01 | fraction of depth span $\Delta_d$ | slider; Python floors at $10^{-12}$ |
| Bandwidth $f$ | range slider | **0.03** | 0.005 – 0.15, step 0.005 | dimensionless covariance factor | slider; JS substitutes 0.03 for any falsy value then floors at $10^{-4}$; Python neither |
| Grid $G$ | select | **120** | 80 / 120 / 160 / 220 | nodes per side | clamp 16–400 (Py), 16–260 (JS, after substituting 120 for any falsy value) |
| Contour levels $K$ | none (fixed) | **8** | request param `levels` | count | frontend always sends 8; empty levels are dropped from the output |
| Colormap | select | `viridis` | viridis, magma, seismic, reds, greys | — | `colormaps.js` |
| Contours | switch | on | on/off | — | client-side only |
| Log scale | switch | **on** | on/off | — | `kde.py` / worker |
| Fit-point cap | none | **6000** | fixed | count | `MAX_KDE_FIT_POINTS`, `fitLimit` |
| Subsample seed | none | **0** | fixed | — | `rng_seed=0`; `randomUnit(0)` |
| Periodic margin $m$ | none | $\min(0.5,\max(0.1,2f,\Delta z))$ | derived | fractional | both paths; JS uses the **un-floored** $f$ here |
| Augmentation factor | none | $(1+2m)^3$ | derived | ratio | 1.73× at $m=0.1$, 8× at $m=0.5$ |
| Log floor | none | $10^{-12}$ | fixed | density | both paths |
| Exponent cutoff | none | $-60$ (browser only) $\Rightarrow 11\sigma$ | fixed | — | `localKdeWorker.js`, `gpuKde.js` |
| Kernel ridge | none | $10^{-8}$ on diag (browser only) | fixed | fractional² | `makeKernel()` |
| GPU work threshold | none | $G^2 n \ge 2{,}000{,}000$ | fixed | work units | `GPU_MIN_WORK` |
| Display atom cap | none | 1 000 000 (clamped to $\ge100$ in Flask) | fixed | atoms | `STRUCTURE_MAX_POINTS`, `MAX_STRUCTURE_POINTS` |
| Plane-section tolerances | none | $10^{-9}$ (on-plane corner), $10^{-8}$ (dedup), $\ge3$ vertices | fixed | fractional | `_plane_section_vertices()` / `planeSectionVertices()` |
| Canvas minimum size | none | 320×260 (KDE), 220×260 (slab) | fixed | CSS px | `StructurePage.jsx`; also applied to the 3× export |
| Debounce | none | 160 ms (Flask) / 80 ms (worker) | fixed | ms | `StructurePage.jsx` |

---

### Python vs JavaScript: exact parity table

Derived by code reading — see the gaps noted above; no cross-runtime numerical test exists.

| Stage | Agreement |
| --- | --- |
| Unit-cell folding | **Exact** (same modulo convention, sign-corrected in JS) |
| Periodic image tiling + margin | **Equivalent for folded inputs** (Python: unconditional originals + 26 shifted offsets; JS: all 27 offsets margin-tested). **Row order differs**, so index-based subsampling picks different points |
| Slab selection in depth | **Exact** (algebraically identical, inclusive both ends) |
| `slabCount` (unique source atoms with an image in the slab) | **Exact** |
| Subsample size (6000) | **Exact**; the **selected subset differs** (PCG64 vs mulberry32, and a different row order) |
| Bandwidth matrix $\mathbf{H}=f^2\mathbf{C}$ | Same formula, but $\mathbf{C}$ is fitted to the **subsampled** points, so the two runtimes use **different $\mathbf{H}$**; JS additionally adds a $10^{-8}$ diagonal ridge and coerces $f$ (0/NaN → 0.03, then a $10^{-4}$ floor) |
| Kernel normalization $1/(2\pi n\sqrt{\det\mathbf{H}})$ | **Exact** (the browser constant is analytically correct, not a fudge) |
| Periodic renormalization $\kappa=N_\mathrm{img}/N_\mathrm{src}$ | **Exact** |
| Evaluation grid nodes | **Exact** on the CPU paths; grid clamp maxima differ (400 vs 260); the GPU path recomputes node positions in `f32` |
| Kernel sum | SciPy: full float64, no cutoff. Browser CPU: float64 with $e<-60$ cutoff. Browser GPU: **float32 inputs, parameters, node grid and accumulator**, same cutoff |
| Degenerate slabs | Python bails to a zero grid on rank/uniqueness tests; JS regularizes and draws |
| $\log_{10}$ transform + $10^{-12}$ floor | **Exact** |
| Contour level values | **Exact** (same formula; both drop levels that yield no polylines) |
| Contour tracing | contourpy stitched polylines with saddle handling vs. per-cell 2-point segments with arbitrary saddle pairing |
| Contours in log mode | Python suppresses all contours when peak $\log_{10}\rho\le0$; JS does not |
| In-plane axes for a **custom** normal | **Differ** (different Gram–Schmidt seed → in-plane rotation/reflection) |
| In-plane axes for a/b/c presets | **Exact** |
| Zero / near-zero custom normal | **Differ**: JS falls back to $(0,0,1)$ at $\lVert\mathbf{h}\rVert\le10^{-9}$; Python raises at $\le10^{-12}$. The app never hits the raise because it sends the already-normalized fallback |
| Input population | Flask: **all** atoms of the element, filtered while parsing. Browser: the element-filtered display array, globally strided (and hard-truncated) above 1 000 000 atoms **before** filtering |
| Display sampler (companion panels) | Flask `/api/structure`: site-stratified quota sampler `_sample_atoms_by_site()`, so above the cap the panels and the KDE use different populations. Browser: the same strided array the KDE uses |
| Element label case | Python `.capitalize()`s the element token, so `SE`, `se` and `Se` **merge into one entry**; the browser parser keeps the raw token, so they stay separate. The dropdown labels, the per-element counts, the legend colours (`buildElementColors` sorts the distinct labels) and hence the KDE population for a filtered element can all differ between runtimes for the same file |

---

### Caveats / what this is not

1. **This is a smoothed picture of one RMC configuration, not a measured density.** Nothing here is
   fitted to data, refined, or error-propagated. The map inherits every artefact of the underlying
   RMCProfile run.
2. **The browser path is a visualization path.** The SciPy path served by `/api/kde/slice` is the
   reference. If a number is going into a figure caption or a paper, take it from a Flask session
   pointed at a run **directory** — note that loading the bundled **Demo** run, even in Flask mode,
   switches the page to the browser worker.
3. **The bandwidth is not a length.** $f$ multiplies the *sample covariance of the slab points*, so
   the physical smoothing width changes with the element filter, the slab thickness, the slice
   normal, and even the periodic margin. It is also computed from the ≤6000 *subsampled* points, so
   the smoothing width is seed-dependent and differs slightly between runtimes. `bw = 0.03` is
   roughly $8\times$ narrower than Scott's rule at $n=6000$ — the map is deliberately under-smoothed
   to resolve sites, which means fine structure in the map can be sampling noise rather than real
   density.
4. **The estimate is subsampled.** At most 6000 slab points (images included) are fitted. The
   estimate stays unbiased but its pointwise noise grows as $\sqrt{N_\mathrm{img}/6000}$. The canvas always
   prints both counts.
5. **The density units are fractional, not Å⁻².** With $\kappa = N_\mathrm{img}/N_\mathrm{src}$ applied, the field integrates
   to about **1 over the whole cell** — the cell carries unit mass; it is *not* 1 per slab atom. It
   is not divided by the slab thickness, so values are not comparable across different $\Delta z$,
   elements, bandwidths or slices.
6. **The colour scale is per-slice.** No colorbar, no shared normalization, no numeric legend. Two
   screenshots of this panel cannot be compared quantitatively. For an oblique normal the scale is
   additionally set by grid nodes that lie **outside** the drawn cross-section and are clipped from
   the display.
7. **The colormaps are 5-anchor approximations** of the matplotlib maps of the same name and are not
   perceptually uniform.
8. **Everything geometric happens in fractional space.** The kernel is Euclidean-isotropic in
   fractional coordinates and therefore anisotropic in Å for any non-cubic cell; the custom
   "direction" is a Miller index triple $(hkl)$, not a real-space vector; and the slab thickness in
   Å must be reconstructed by hand as $\Delta z(|h|+|k|+|l|)d_{hkl}$.
9. **The slider auto-jumps.** Changing the element filter or the normal — including typing a single
   digit into a custom-direction box — re-runs the 50-bin "densest layer" search and overwrites $z_c$.
   That search runs on the unwrapped, display-sampled population.
10. **The GPU result is float32** in its inputs, its kernel parameters, its reconstructed node grid
    and its accumulator. Structurally identical to the CPU loop, numerically equal to about
    $10^{-6}$–$10^{-5}$ relative — fine for a picture, not a bit-for-bit guarantee, and untested.
11. **In log mode the Flask path can silently drop all contours** when the peak density is below 1
    per unit fractional area, while the browser path draws them.
12. **The periodic wrap is exact only out to the margin $m$.** Images farther than $m$ from the cube
    are discarded, which truncates the kernel at $\approx9.5\sigma$ at the defaults and
    $\approx5.7\sigma$ at $f=0.15$. This applies to the SciPy path too.
13. **`slabCount` is "atoms with at least one image in the slab"**, and one atom can contribute
    several rows near an edge or corner. The drawn band and the highlighted atoms in the side view
    are clamped/unwrapped and so **understate** the selection for $z_c$ near 0 or 1.
14. **`"No atoms in this slab"` also means "the estimator declined this slab"** — fewer than 5 rows,
    fewer than 3 distinct in-plane points, a rank-deficient spread, or a suppressed `LinAlgError`
    (e.g. `bw=0` on the API). Look at `fitCount`.
15. **A missing lattice block degrades silently.** If `structure.latticeVectors` or
    `structure.supercell` is absent the frontend draws a 1 Å cubic cell with no warning.
16. **The page may not be analysing the file you think.** In a folder with several `.rmc6f`
    configurations, the one paired with a recognised output file wins; with no match Flask falls back
    to the alphabetically first one and the browser to the first in directory-walk order, so the two
    paths can pick different models (Step 1).
17. **A zero-atom browser parse is still listed as an open issue** (`AGENTS.md`, *Current known
    issues*, 2026-06-18) — but **its stated rationale is stale and should not be quoted**: it blames
    a parser that "assumes one exact atom-line format", whereas `parseAtomLine()` now handles the
    5–6-field coords-only form, the 9-field legacy form and the 10+-field modern form by indexing
    from the end (Step 1). If a zero-atom failure mode survives, it is somewhere else — the `Atoms:`
    header match, or the lattice/supercell header parse — not in the atom-line field layout.
18. **No error bars, no resolution function, no thermal deconvolution.** The width of a blob in this
    map is the convolution of the true site spread with the KDE kernel; the kernel is *not*
    deconvolved. For quantitative displacement parameters use the
    [PCA Ellipsoid](pca-ellipsoid.md) page (`rmc_toolkits/pca_kde.py`), which fits the displacement
    cloud directly.


## Structure page — Slab In Cell projection and the 3D unit-cell view

### What this page shows

The **Atomic Density** tab (nav label in [App.jsx](../../web_app/frontend/src/App.jsx); the page's own
heading reads *KDE And Folded Unit Cell*) is rendered by
[StructurePage.jsx](../../web_app/frontend/src/components/StructurePage.jsx). It puts three panels side
by side, all driven from one shared slice definition:

1. **KDE Slice** — the in-plane density map of the selected slab (the density estimator itself is
   documented in the KDE section; this section documents the *projection* that puts it on screen).
2. **Slab In Cell** — a side view down one in-plane crystal direction, showing where the slab sits
   inside the unit cell and which atoms fall in it. The highlighted band is draggable.
3. **Folded Unit Cell** — a Three.js scene: every atom of the (element-filtered) model folded into a
   single unit cell, the cell edge frame, and the slab rendered as two translucent cross-sections.

Everything on this page operates on **atoms folded into one unit cell**. The model is an RMCProfile
supercell of $N_1 \times N_2 \times N_3$ cells; all of them are superimposed. None of the three
panels shows a single physical cell of the configuration — they show the *ensemble* of all cells.

**Two data paths, and what actually selects between them.** `StructurePage` branches on
`isLocalStructure = Boolean(localRun)`, **not** on the deployment mode:

- **browser path** — a local run is loaded. The `.rmc6f` is parsed by
  [localStructureWorker.js](../../web_app/frontend/src/workers/localStructureWorker.js) and the density
  by [localKdeWorker.js](../../web_app/frontend/src/workers/localKdeWorker.js), both in Web Workers.
- **backend path** — `localRun` is `null`. The page calls `GET /api/structure` and
  `GET /api/kde/slice` on the Flask server ([app.py](../../web_app/backend/app.py)).

Only the **folder picker** is static-mode-only (`fsAccess = staticMode && supportsFileSystemAccess()`
in `App.jsx`, with a `webkitdirectory` input as the static fallback). The **Demo** button is rendered
unconditionally and calls `setLocalRun(...)`, so a *Flask deployment showing the bundled demo run
takes the browser path* — none of the backend-path statements below apply to it. Conversely, a static
deployment with no run loaded shows only the error `Open a run folder to view the structure.` This
section therefore says "browser path" / "backend path" and never "static mode" / "Flask mode" when
describing which code runs.

#### Notation and units

| Symbol | Meaning | Units |
| --- | --- | --- |
| $\mathbf{L}_i$ | supercell (box) lattice vectors, $i=1,2,3$ | Å |
| $N_i$ | supercell multiplicity along axis $i$ | — |
| $\mathbf{A}_i = \mathbf{L}_i / N_i$ | unit-cell lattice vectors | Å |
| $\ell_i = \lVert\mathbf{A}_i\rVert$, $\ell_{\max}=\max_i \ell_i$ | unit-cell edge lengths | Å |
| $\mathbf{f}=(f_1,f_2,f_3)$ | atom coordinate as stored in `.rmc6f` — fraction of the **box** | — |
| $(n_1,n_2,n_3)$ | per-atom supercell indices from `.rmc6f` | — |
| $\mathbf{x}=(x_1,x_2,x_3)$ | folded coordinate, fraction of **one unit cell**, $x_i\in[0,1)$ | — |
| $\mathbf{c}_j$, $j=1\ldots8$ | the eight unit-cube corners (`CUBE_CORNERS`) | fractional |
| $\mathbf{c}$ | body centre of the *normalized* cell (Step 2) — never a cube corner | normalized cell units |
| $\hat{\mathbf{h}}$ | slice normal, components along the fractional axes (unit length in that space) | — |
| $\hat{\mathbf{u}},\hat{\mathbf{v}}$ | in-plane basis, also in fractional-component space | — |
| $\mathbf{e}_u,\mathbf{e}_v,\mathbf{e}_h$ | Å-space images of $\hat{\mathbf{u}},\hat{\mathbf{v}},\hat{\mathbf{h}}$, i.e. $\sum_i \hat u_i\mathbf{A}_i$ etc. | Å |
| $u=\hat{\mathbf{u}}\cdot\mathbf{x}$, $v=\hat{\mathbf{v}}\cdot\mathbf{x}$ | in-plane coordinates of an atom | — |
| $\mathbf{R}=u\mathbf{e}_u+v\mathbf{e}_v+d\mathbf{e}_h$ | the atom's true Å position | Å |
| $d=\hat{\mathbf{h}}\cdot\mathbf{x}$ | out-of-plane (depth) coordinate of an atom, along the slice normal | — |
| $[d_{\min},d_{\max}]$, $\Delta_d = d_{\max}-d_{\min}$ | projection range of the unit cube along $\hat{\mathbf{h}}$ | — |
| $\tilde d = (d-d_{\min})/\Delta_d$ | normalized depth (`pointDepth`) | — |
| $z_c$, $\Delta z$ | slice centre and slab thickness, as fractions of $\Delta_d$ | — |
| $\delta = \Delta z\,\Delta_d$; $d_\mathrm{start}, d_\mathrm{end}$ | band depth thickness and its cell-clipped edges | — |
| $\mathbf{P},\mathbf{Q}$ | the two Å-space vectors handed to `makePlane` | Å |
| $\mathbf{p},\mathbf{q}$ | their flattened 2D images | Å |
| $X_{\min}\ldots Y_{\max}$, $\Delta X, \Delta Y$ | bounds and spans of the flattened plane | Å |
| $W_{px}, H_{px}$ | canvas width / height | CSS px |
| $k$ | canvas scale factor | CSS px / Å |
| $o_x,o_y$ | centring offsets of the fitted content | CSS px |
| $\rho$ | KDE density sample on the grid (see the KDE section for its normalization) | probability density per in-plane fractional unit² (integrates to ≈1 over the cell — KDE Step 6) |
| $v_{\min},v_{\max}$ | min / max of the density grid actually drawn (after the $\log_{10}$ toggle) | as $\rho$ |
| $R_{sph}$ | camera framing radius (Step 7.5) | normalized cell units |

---

### Step 1 — Parse the `.rmc6f` configuration and fold it into one cell

**Inputs:** one `.rmc6f` file, *chosen* from the run folder by the heuristic in 1a.

**Outputs (the `structure` object both paths produce):**
`{ source, totalAtoms, sampledAtoms, sampleStride, elements, elementCounts, atomIndices, supercell,
latticeVectors, points[] }`, where each point carries `element`, `referenceNumber`, the raw box
coordinates `boxX/boxY/boxZ`, and the folded unit-cell coordinates `x, y, z`. `elements` is **sorted**
in both paths (`sorted(counts.keys())` in `app.py`, `Object.keys(counts).sort()` in
[browserData.js](../../web_app/frontend/src/browserData.js)), which is what makes the element dropdown
order and the colour assignment (Step 8) deterministic. The browser path additionally returns `basis`
(one circular-mean site per reference number, with a per-site rms displacement `dispA`) and `moves`;
`basis` is what feeds the symmetry card of `ModelSummary`, so **that card only appears on the browser
path**. `ModelSummary`, rendered at the top of this page, displays `source` (basename), `totalAtoms`,
the cell lengths/angles derived from `latticeVectors`, and per-element site counts read from
`atomIndices`.

Two independent implementations exist and the page uses whichever data path it is on.

#### 1a. Which `.rmc6f` is read (a heuristic, not a given)

A run folder can hold several models, and the choice decides what all three panels show.

`app.py` → `_find_rmc6f()`: if the resolved target is itself a `.rmc6f` file, use it. Otherwise glob
`*.rmc6f` (error if none), then walk the directory's other files in case-insensitive name order and
classify each by `_run_stem_from_output_name()`, which extracts a run stem at one of three
priorities:

| priority | pattern |
| --- | --- |
| 0 | `<stem>-NN.log` (2+ digits) |
| 1 | `<stem>-EXAFS-*_Q_OUTPUT.csv` / `_R_OUTPUT.csv`, `<stem>_FT_XFQ<n>.csv`, `<stem>_FQ<n>.csv`, `<stem>_SQ<n>.csv`, `<stem>_bragg*.csv`, `<stem>_PDF*.csv` |
| 2 | `Frac_coord_<stem>.txt` |

The `(priority, filename)` list is sorted and the first `.rmc6f` whose stem matches wins; with no
match, the **alphabetically first** `.rmc6f` is used.

`browserData.js` → `chooseStructureFile()` uses the same priority ladder but matches on
`dirname + '/' + stem`, so an output file only claims an `.rmc6f` sitting in the *same* subfolder; its
fallback is the first `.rmc6f` in directory-walk order, which is not necessarily alphabetical. A
folder holding several models can therefore resolve to a different model on the two paths.

#### 1b. Metadata

Both parsers scan for two headers and are otherwise position-independent:

- a line whose first token is `Supercell` → $N_i$ from the **last three** tokens
  (`Supercell dimensions:  10 10 10`);
- a line whose first token is `Lattice` → the **next three lines** parsed as the rows
  $\mathbf{L}_1,\mathbf{L}_2,\mathbf{L}_3$ in Å.

Python: `rmc_toolkits/parsers.py` → `read_cell_vectors()`. JavaScript:
[browserData.js](../../web_app/frontend/src/browserData.js) → `readCellVectors()`. The two agree
exactly; both raise/throw if either header is missing. Note that the `Cell (Ang/deg): a b c α β γ`
line, when present, is **ignored** — the cell geometry always comes from the lattice-vector rows, so
a triclinic cell is handled by construction.

#### 1c. Atom lines

Atom records start after a line whose first token is `Atoms:`.

Python — `rmc_toolkits/parsers.py` → `iter_rmc6f_atoms()` — requires **≥ 9 whitespace tokens** and
indexes from the *end* of the line, so both the current and the older layouts parse:

```
current (10 fields):  id  element  [type]  f1 f2 f3  ref  n1 n2 n3
older    (9 fields):  id  element          f1 f2 f3  ref  n1 n2 n3
```

i.e. `coords = parts[n-7:n-4]`, `reference_number = parts[n-4]`, `cell_indices = parts[n-3:n]`. The
element is stored **capitalized** (`parts[1].capitalize()`); everything between the element and the
coordinates is joined into `type_label`.

JavaScript — [rmc6f.js](../../web_app/frontend/src/rmc6f.js) → `parseAtomLine()` — is more tolerant. It
accepts three shapes:

- `n < 5` → rejected;
- `5 ≤ n ≤ 6` → **coords-only** legacy form (`id element [label] f1 f2 f3`): the last three tokens are
  the coordinates; `referenceNumber` and `cellIndices` come back `null`;
- `7 ≤ n ≤ 8` → rejected as a truncated full line;
- `n ≥ 9` → full form, indexed from the end exactly as Python does.

Non-finite coordinates, reference number, or cell indices reject the line. The element token is
**not** capitalized. `web_app/frontend/src/__tests__/rmc6f.test.js` pins all five branches.

> **Cross-runtime differences (real).**
>
> 1. A coords-only (5–6 field) `.rmc6f` yields **zero atoms on the backend path** — `iter_rmc6f_atoms()`
>    does `if n < 9: continue`, so such a line is never yielded — but parses fine in the browser.
> 2. An element written `GA` or `se` becomes `Ga`/`Se` server-side but stays `GA`/`se` in the browser,
>    which changes the element filter list, the legend, and the colour lookup between the two paths.
> 3. **The capitalization is not even uniform inside one backend payload.** `elements` and
>    `elementCounts` come from `iter_rmc6f_atoms()` and are capitalized, but `atomIndices` comes from
>    `read_atom_indices()`, which keys on the **raw** token `parts[1]` (and skips any line with fewer
>    than 5 tokens or a non-integer `parts[-4]`). For a file whose element column is not already
>    title-case, `atomIndices['Se']` is undefined while `elementCounts['Se']` exists, so the
>    per-element **site counts** in `ModelSummary` read 0 (the *total* site count still sums all
>    `atomIndices` values and stays correct). The CLI element filter in `read_structure()` goes
>    through `read_atom_indices()` too, so it matches the raw casing.

#### 1d. Folding

Backend path ([app.py](../../web_app/backend/app.py) → `structure()`):

$$\mathrm{reduced}_i = f_i - \frac{n_i}{N_i}, \qquad x_i = \big(\mathrm{reduced}_i \cdot N_i\big) \bmod 1$$

Browser path ([browserData.js](../../web_app/frontend/src/browserData.js) → `structureFromRmc6f()`):

$$x_i = \big(((f_i N_i) \bmod 1) + 1\big) \bmod 1$$

**When both run, they are algebraically identical**: subtracting $n_i/N_i$ removes an integer after
multiplication by $N_i$, so it cannot change the value mod 1, and the extra `+1 %1` in JS only fixes
JavaScript's sign-preserving `%` for negative coordinates (Python's `%` is already non-negative). To
floating point they agree exactly.

They are **not** equally applicable, though: only the JS form is index-free. The backend form
evaluates `atom["coords"] - (atom["cell_indices"] / supercell)`, which needs the cell-index columns,
and its parser drops any atom line with fewer than 9 tokens — so a file without per-atom cell indices
produces **no atoms at all** on the backend path (consistent with the box in 1c), not a differently
folded set. `load_unit_cell_positions()` in `rmc_toolkits/kde.py` uses the index-free form
`(coords * supercell) % 1.0`, but it reads atoms through the same 9-token parser, so it inherits the
same restriction.

#### 1e. Subsampling (the two paths use *different* strategies)

Both are capped at `maxPoints = 1 000 000` (`STRUCTURE_MAX_POINTS` in `StructurePage.jsx`,
`MAX_STRUCTURE_POINTS` in `app.py`, where the query argument is clamped to $[100,\,10^6]$).

- **Backend** — `app.py` → `_sample_atoms_by_site()`. First: `if len(atoms) <= max_points: return
  atoms, 1` — an **early return** with stride 1 and no grouping at all. Above the cap, atoms are
  grouped by `reference_number` (crystallographic site); each group gets
  `quota = max(1, max_points // n_sites)` and is strided by `max(1, len(group) // quota)`, keeping
  `group[::stride][:quota]`; the concatenation is finally truncated with `sampled[:max_points]`. This
  is **site-stratified**: every site keeps ≥ 1 atom and equal representation. The reported
  `sampleStride` is only a summary, `max(1, N_atoms // max_points)`, not the stride actually applied
  to any group.
- **Browser** — `structureFromRmc6f()`: a single global stride
  $s_\mathrm{stride}=\max(1,\lceil N_\mathrm{atoms}/\mathrm{maxPoints}\rceil)$, keeping every $s_\mathrm{stride}$-th atom, then truncating
  to `maxPoints`. A rare element can be lost entirely this way; the site-stratified backend path
  cannot lose one.

With the default cap and a typical run ($\sim$52 000 atoms for the repository's GNSe sample) **no
subsampling happens on either path** — on the backend it is the `len(atoms) <= max_points` early
return that does it (the quota code never executes), and in the browser $s=1$. The difference only
bites for models above one million atoms.

`elements` / `elementCounts` are always computed over **all** atoms, not the sampled subset, on both
paths.

#### 1f. How the data reaches the page, and what re-triggers what

The load effect (deps `[directory, localRun]`) has five branches, in order:

1. `localRun.structure` present → used directly. *Never populated by the current `App.jsx`*:
   `makeRunFromEntries()` sets `structure: null`.
2. `localRun.structureFile` present → posted to `localStructureWorker` with
   `maxPoints = STRUCTURE_MAX_POINTS`.
3. `localRun` present with neither → the error `localRun.structureError` (`'No model structure
   detected'`), else `'No structure data available in this folder'`.
4. no `localRun` and `isStaticMode()` → the error `'Open a run folder to view the structure.'`
5. otherwise → `GET /api/structure?dir=<directory>&maxPoints=<STRUCTURE_MAX_POINTS>`.

Two module Workers (structure and KDE) are created per mount whenever `isLocalStructure`, and
terminated on cleanup. Both use **monotonic request ids**: each post increments
`localStructureRequestRef` / `localKdeRequestRef`, and `worker.onmessage` drops any reply whose
`event.data.id` is not the current one, so a slow earlier slice can never overwrite a newer one.

The KDE effect debounces **80 ms** on the browser path and **160 ms** on the backend path, and the
backend request carries an `AbortController` signal that is aborted on every dependency change
(cancellation errors are swallowed via `axios.isCancel` / `ERR_CANCELED`). Its dependency list is
`[structure, isLocalStructure, points, directory, selectedElement, sliceDirection, sliceConfig,
zCenter, thickness, bandwidth, gridSize, logScale]`.

Consequence to keep in mind when reading the rest of this section: **one slider tick of $z_c$
re-triggers a KDE request, both 2D draw effects, and a complete teardown/rebuild of the Three.js
scene** (Step 7), the last with camera-state restore so the viewpoint survives.

#### 1g. The element filter — the first transformation applied to every panel

```js
points = selectedElement === 'all' ? structure.points
                                   : structure.points.filter(p => p.element === selectedElement);
```

a `useMemo` on `[structure, selectedElement]` using **strict equality** on the parsed symbol. That
filtered array is what drives the Slab In Cell markers (Step 5), the 3D point clouds (Step 7), the
50-bin auto-centre histogram (Step 3d), and — on the browser path — the point set posted to the KDE
worker.

On the **backend path the KDE panel is not filtered client-side at all**: `selectedElement` is sent as
the `element=` query argument and applied server-side by `rmc_toolkits/kde.py` →
`load_unit_cell_positions()`, which skips atoms whose symbol differs from it
(`element in (None, '', 'all')` means no filter) and is memoized per `(path, mtime, element)` by
`_cached_positions`. Both sides of that comparison come from `iter_rmc6f_atoms()`, so they are
consistently capitalized; but the same control therefore acts through **two different mechanisms at
two different points of the pipeline**, and the symbol it matches is the capitalized one on the
backend path and the file's raw token on the browser path.

The element colour map and legend are *not* filtered — see Step 8.

#### 1h. Related but not used by this page: the `Frac*.txt` conversion

`rmc_toolkits/parsers.py` → `frac_lines_from_rmc6f()` writes the classic `Frac_coord_*.txt` file
(exposed as `POST /api/convert/frac` → `write_frac_from_rmc6f()`): a 5-line header followed by one
row per atom,

$$\mathrm{reduced}_i = f_i - \frac{n_i}{N_i} \quad\text{printed as}\quad \texttt{RN  x  y  z  Nx  Ny  Nz}$$

with the coordinates formatted to **5 decimal places** of a *box* fraction. `read_structure()` reads
that file back, skipping exactly the first 5 lines, and re-expands
$\mathbf{x} = (\mathrm{reduced}\cdot\mathbf{N}) \bmod 1$, optionally converting to Cartesian
$\mathbf{r}=\sum_i x_i\mathbf{A}_i$ (`mode="cartesian"`, the default; `mode="fractional"` returns
$\mathbf{x}$).

Two things to know:

- **Precision loss.** 5 decimals of a box fraction quantizes positions to $10^{-5}\lVert\mathbf{L}_i\rVert$;
  for the GNSe sample ($\lVert\mathbf{L}\rVert = 104.116$ Å) that is $\approx 1.0\times10^{-3}$ Å.
  Neither `/api/structure` nor `kde.py` uses this path — both read the `.rmc6f` directly at full
  precision — so the web page is unaffected. Only the package/CLI `read_structure()` consumers see
  the quantization.
- **`RmcStructure.atom_types` does not hold element symbols.** `read_structure()` appends `parts[0]`,
  which is the *reference number* column of the Frac file. `tests/test_parsers.py` →
  `test_read_structure_loads_full_folded_unit_cell` asserts `len(set(atom_types)) == 52`, i.e. the
  number of reference sites, confirming the field's actual content. The element filter still works,
  because `read_structure()` maps element → reference numbers through `read_atom_indices()`.

---

### Step 2 — Unit-cell basis and its normalization

**Input:** `latticeVectors` ($\mathbf{L}_i$), `supercell` ($N_i$). **Output:** the `unitCell` memo in
`StructurePage.jsx`.

$$\mathbf{A}_i = \frac{\mathbf{L}_i}{\max(N_i,\,10^{-12})}, \qquad
\ell_i = \lVert\mathbf{A}_i\rVert, \qquad
\mathbf{b}_i = \frac{\mathbf{A}_i}{\max(\ell_1,\ell_2,\ell_3,\,10^{-9})}, \qquad
\mathbf{c} = \tfrac12\sum_i \mathbf{b}_i$$

- `unitVectors` = $\mathbf{A}_i$ in **Å** — used by every 2D projection, so the 2D panels carry a real
  metric.
- `basis` = $\mathbf{b}_i$, the same cell scaled so its **longest edge is exactly 1** — used only by
  the Three.js scene, so the scene is resolution- and material-independent of the physical cell size.
  Note the $10^{-9}$ floor inside the normalization: it is what protects the 3D basis from a
  degenerate (zero-length) cell, alongside the $10^{-12}$ floor on $N_i$.
- `center` = $\mathbf{c}$, the body centre of the normalized cell, subtracted from every 3D vertex so
  the scene is centred on the origin.

If `latticeVectors` or `supercell` is missing the memo degrades to the identity cell
($\mathbf{A}_i=\hat{e}_i$, lengths $1$), which keeps the page from crashing but silently draws a cube.

Conversion from a fractional triple to Å (or to normalized scene units) is
`vectorFromFraction(fraction, basis)` $=\sum_i \mathrm{fraction}_i \cdot \mathrm{basis}_i$.

---

### Step 3 — Defining the slice: normal, in-plane basis, depth range

**Inputs:** the `Normal` control (`a`, `b`, `c`, or `Custom` with three numeric fields, default
$[1,1,0]$). **Output:** `sliceConfig = { key, label, normal, u, v, uLabel, vLabel, range }`, built by
`makeSliceConfig()`.

#### 3a. Presets

| key | $\hat{\mathbf{h}}$ | $\hat{\mathbf{u}}$ (label) | $\hat{\mathbf{v}}$ (label) |
| --- | --- | --- | --- |
| `a` | $[1,0,0]$ | $[0,1,0]$ (b) | $[0,0,1]$ (c) |
| `b` | $[0,1,0]$ | $[1,0,0]$ (a) | $[0,0,1]$ (c) |
| `c` | $[0,0,1]$ | $[1,0,0]$ (a) | $[0,1,0]$ (b) |

These are `SLICE_PRESETS` in `StructurePage.jsx` and match `SLICE_ORIENTATIONS` in
[app.py](../../web_app/backend/app.py) exactly; `_slice_orientation_from_request()` passes both $\hat{\mathbf u}$
and $\hat{\mathbf v}$ through to the estimator for a preset, so for `a`/`b`/`c` the two runtimes share
the in-plane frame. Note that the `b` preset triad is **left-handed**
($\hat{\mathbf{u}}\times\hat{\mathbf{v}} = -\hat{\mathbf{h}}$), so a `b`-normal view is mirrored
relative to a right-handed convention. Both runtimes share the convention, so they agree with each
other.

#### 3b. Custom normal

$\hat{\mathbf h} = \texttt{normalize}(\texttt{customDirection},\,[0,0,1])$, then
`makeFreePlaneBasis(normal)` builds an in-plane pair by Gram–Schmidt **in fractional-component
space**:

$$\mathbf{r} = [1,0,0]\ \text{if}\ |h_1| < 0.85,\qquad \mathbf{r} = [0,1,0]\ \text{otherwise}$$

$$\hat{\mathbf{u}} = \widehat{\mathbf{r} - (\mathbf{r}\cdot\hat{\mathbf{h}})\hat{\mathbf{h}}},\qquad
\hat{\mathbf{v}} = \widehat{\hat{\mathbf{h}}\times\hat{\mathbf{u}}}$$

`normalize()` falls back to a fixed vector when the length is $\le 10^{-9}$ (`[0,1,0]` for
$\hat{\mathbf{u}}$, `[0,0,1]` for $\hat{\mathbf{v}}$).

Because the construction is done on fractional components, $\hat{\mathbf{u}}$ and $\hat{\mathbf{v}}$
satisfy $\hat{\mathbf{h}}\cdot\hat{\mathbf{u}} = \hat{\mathbf{h}}\cdot\hat{\mathbf{v}} = 0$, which is
exactly the condition for a lattice direction to **lie in** the plane family with Miller indices
$\propto \hat{\mathbf{h}}$. So the two axes genuinely span the crystallographic plane. They are
**not** orthogonal in Å space for a non-cubic cell — Step 4 handles that correctly.

> **Degenerate input.** The same `normalize(customDirection, [0,0,1])` fallback applies to the
> **normal itself**: entering all zeros (clearing a field of the number input yields `Number('') = 0`)
> silently reverts the view to the **c-axis slice** while the panel label, built from the raw
> `customDirection`, still reads `[0 0 0]`. `updateCustomDirection` stores `Number(event.target.value)`
> with no finiteness check; the `length ≤ 1e-9` guard would not catch a NaN component either, though a
> `type="number"` input reports `''` (hence `0`) rather than `NaN` for unparseable text, so the
> all-zero case is the reachable one.

> **Cross-runtime difference (real).** `rmc_toolkits/kde.py` → `_plane_basis()` / `_orthogonal_axis()`
> picks its seed axis as the Cartesian axis with the **smallest** $|h_i|$
> (`np.eye(3)[argmin(|h|)]`), whereas `makeFreePlaneBasis()` picks $x$ unless $|h_1| \ge 0.85$. For
> $\hat{\mathbf{h}} \propto [1,1,0]$ the JS basis is $\hat{\mathbf{u}}=[0.7071,-0.7071,0]$,
> $\hat{\mathbf{v}}=[0,0,-1]$ while the Python basis is $\hat{\mathbf{u}}=[0,0,1]$,
> $\hat{\mathbf{v}}=[0.7071,-0.7071,0]$ — i.e.
> $(\hat{\mathbf u},\hat{\mathbf v})_\mathrm{Py} = (-\hat{\mathbf v},\,\hat{\mathbf u})_\mathrm{JS}$,
> a 90° rotation. The KDE canvas draws using the
> `uVector`/`vVector` the server returns, so it stays internally consistent, but on the **backend path
> with a custom normal the KDE panel and the Slab In Cell panel do not share an in-plane
> orientation**. On the browser path the worker is handed `sliceConfig.u/v`, so the two panels agree.
> Presets are unaffected (3a).

#### 3c. Depth coordinate and range

`projectionRange(normal)` evaluates $\hat{\mathbf{h}}\cdot\mathbf{c}_j$ over the eight unit-cube
corners `CUBE_CORNERS` and returns $[d_{\min}, d_{\max}]$. For the presets this is $[0,1]$; for
$\hat{\mathbf{h}}\propto[1,1,0]$ it is $[0,\sqrt2]$.

The slider quantities are **fractions of that range**, not of a cell edge:

$$\tilde d(\mathbf{x}) = \frac{\hat{\mathbf{h}}\cdot\mathbf{x} - d_{\min}}{\Delta_d}
\quad(\texttt{pointDepth}), \qquad
\text{in slab} \iff \big|\tilde d - z_c\big| \le \tfrac{\Delta z}{2}\quad(\texttt{inActiveSlab})$$

with $\Delta_d$ replaced by 1 if it evaluates to 0.

The same **depth convention** is used by the browser KDE worker
([localKdeWorker.js](../../web_app/frontend/src/workers/localKdeWorker.js) → `makeSlab`) and by the
server (`rmc_toolkits/kde.py` → `oriented_kde_slice`, which converts to absolute depth as
$d_c = d_{\min} + z_c\Delta_d$ and $\delta = \Delta z\,\Delta_d$ before selecting, after clamping
$z_c$ to $[0,1]$ and $\Delta z$ to $\ge 10^{-12}$). The formula is the same in all three.

> **But the predicate is applied to different point sets.** Both KDE implementations first tile
> periodic images — `_augment_periodic_images()` / `augmentPeriodicImages()`, keeping every image of
> every atom whose fractional coordinates fall inside $[-m,\,1+m]^3$ with
> $m = \min(0.5,\ \max(0.1,\ 2\,\mathrm{bw},\ \Delta z))$ — and run the slab test over that **augmented** cloud,
> so images with $\tilde d$ outside $[0,1]$ can enter the slab. Their reported `slabCount` is then the
> number of **unique source atoms** contributing (`sources` Set / `np.unique(source_index[mask])`).
> `inActiveSlab()` in `StructurePage.jsx` runs over the **folded points only**, with no tiling. The
> atoms the KDE integrates are therefore not the atoms the Slab In Cell panel highlights; see Step 5.5.

#### 3d. Auto-centring the slice on the densest layer

On every change of `points`, `sliceConfig` or `pointDepth` an effect histograms $\tilde d$ into
**50 equal bins**,

$$\mathrm{bin} = \max\!\big(0,\ \min(49,\ \lfloor \tilde d\cdot 50\rfloor)\big), \qquad
z_c \leftarrow \frac{\text{argmax bin} + 0.5}{50}$$

so the view opens on an atomic layer rather than in a gap. Details that matter:

- the bin index is **clamped** to $[0,49]$, so an out-of-range depth (only possible for a degenerate
  range) lands in the end bin rather than crashing;
- the scan uses a strict `>`, so **ties resolve to the lowest bin**;
- the effect returns immediately when the point set is empty (`if (!points.length) return;`),
  leaving the previous $z_c$ untouched;
- the histogram runs over the **element-filtered** `points` (Step 1g), so "the densest layer" means
  the densest layer *of the selected species*.

Consequence: changing the element filter or the normal **resets a hand-set slice position**; dragging
the band (Step 6) does not, because dragging changes neither dependency.

---

### Step 4 — The plane mapper: crystal plane coordinates → canvas pixels

This is the geometric core shared by the KDE Slice and Slab In Cell canvases.

#### 4a. Isometric flattening of an oblique plane (`makePlane`)

**Input:** two Å-space vectors $\mathbf{P},\mathbf{Q}$ (e.g. $\mathbf{e}_u,\mathbf{e}_v$).
**Output:** a 2D basis $(\mathbf{p},\mathbf{q})$, the corner bounds of the unit parallelogram, and its
aspect ratio.

$$\cos\theta = \mathrm{clamp}\!\left(\frac{\mathbf{P}\cdot\mathbf{Q}}{\max(\lVert\mathbf{P}\rVert\lVert\mathbf{Q}\rVert,\,10^{-12})},-1,1\right),\quad
\sin\theta = \sqrt{\max(0,\,1-\cos^2\theta)}$$

$$\mathbf{p} = \big(\lVert\mathbf{P}\rVert,\;0\big),\qquad
\mathbf{q} = \big(\lVert\mathbf{Q}\rVert\cos\theta,\;\lVert\mathbf{Q}\rVert\sin\theta\big)$$

This preserves the Gram matrix exactly ($\mathbf{p}\cdot\mathbf{p}=\mathbf{P}\cdot\mathbf{P}$,
$\mathbf{q}\cdot\mathbf{q}=\mathbf{Q}\cdot\mathbf{Q}$, $\mathbf{p}\cdot\mathbf{q}=\mathbf{P}\cdot\mathbf{Q}$),
so the 2D picture is a true **isometry** of the plane $\mathrm{span}(\mathbf{P},\mathbf{Q})$: lengths
and angles measured on the canvas are real Å lengths and real angles, up to the single uniform scale
$k$. $\sin\theta \ge 0$ always, so the flattening is orientation-fixing and never mirrors. Note it
preserves **only** the in-plane metric — nothing about the third direction (see 4d).

`makeProjectedPlane(P, Q, uvPoints)` runs `makePlane` and then re-derives the bounds from an explicit
list of plane coordinates $(s_j,w_j)$ mapped as $s_j\mathbf{p}+w_j\mathbf{q}$. The three call sites
pass different lists, and that is where the CSS box and the drawn content can disagree:

| call site | points passed | drives |
| --- | --- | --- |
| `drawKdeSlice` | the exact section polygon `planePolygon` (Step 9.2) | the KDE canvas fit |
| `drawSlab` | the 8 projected cube corners **+ the 4 band corners** | the slab canvas fit |
| `slicePanelGeometry` | the 8 projected cube corners only | the CSS `--panel-aspect` |

**Two aspects, not one.** `slicePanelGeometry` is a memo returning `{ planeAspect, sideAspect }`, both
computed as $\Delta X / \max(\Delta Y, 10^{-9})$ from the bounding box of the eight projected cube
corners in the **local** basis:

- $\texttt{planeAspect}$ = aspect of $\texttt{makeProjectedPlane}\big(\mathbf U,\mathbf V,\{(\hat{\mathbf u}\cdot\mathbf c_j,\ \hat{\mathbf v}\cdot\mathbf c_j)\}\big)$ → `--panel-aspect` on the **KDE** panel;
- $\texttt{sideAspect}$ = aspect of $\texttt{makeProjectedPlane}\big(\mathbf U,\mathbf H,\{(\hat{\mathbf u}\cdot\mathbf c_j,\ \hat{\mathbf h}\cdot\mathbf c_j)\}\big)$ → `--panel-aspect` on the **Slab In Cell** panel;
- the **3D** panel uses $\max(\texttt{planeAspect},\,1)$.

[StructurePage.css](../../web_app/frontend/src/components/StructurePage.css) consumes `--panel-aspect` as
`aspect-ratio` on `.kde-canvas`, `.slab-panel canvas` and `.three-mount`. Because the CSS aspect comes
from the cube corners in the local basis while the KDE canvas fits the exact section polygon in the
*server's* basis and the slab canvas fits corners **plus** band corners, the CSS box is only an
approximation of the drawn extent — the mapper letterboxes the content, never stretches it.

**A caveat on the slab fit.** The point list `drawSlab` passes covers the projected cube corners and
the band corners, but *not* the corners of the rectangle it actually strokes as the cell outline
($(u_{\min},d_{\min})\ldots(u_{\max},d_{\max})$, Step 5.4), which are generally not among the eight
projected corners. For a custom normal on an oblique cell those rectangle corners can fall outside the
fitted bounds — e.g. $\hat{\mathbf h}\propto[1,1,0]$ with a 120° angle between $\mathbf U$ and
$\mathbf H$ puts $(u_{\max},d_{\min})$ at plane-$X\approx0.707$ while the fitted $X_{\max}\approx0.35$
— so part of the cell outline is drawn outside the 18 px padded box and can be clipped by the canvas
edge.

#### 4b. Fit-and-centre (`makePlaneMapper`)

**Inputs:** a plane from 4a, the canvas $W_{px}$/$H_{px}$ in **CSS pixels**, and `padding = 18` CSS px
(the default, and the value passed explicitly at both call sites).

$$\Delta X = X_{\max}-X_{\min},\quad \Delta Y = Y_{\max}-Y_{\min}\quad(\text{each replaced by }1\text{ if }0)$$

$$k = \min\!\left(\frac{W_{px} - 2\cdot 18}{\Delta X},\ \frac{H_{px} - 2\cdot 18}{\Delta Y}\right)\ \text{[px/Å]},\qquad
o_x = \frac{W_{px} - k\,\Delta X}{2},\quad o_y = \frac{H_{px} - k\,\Delta Y}{2}$$

$$\texttt{map}(u,v):\quad
\mathbf{r} = s\,\mathbf{p} + w\,\mathbf{q},\qquad
X = o_x + (r_x - X_{\min})\,k,\qquad
Y = o_y + (Y_{\max} - r_y)\,k$$

The single $k$ for both axes means **no anisotropic stretching**; the $Y$ flip puts $+w$ upward. The
offsets recentre the fitted content, so the effective padding is symmetric and at least 18 px on the
tight axis.

#### 4c. `invert()` — screen → plane coordinates

Used by the drag handler. It undoes the affine map analytically:

$$r_x = \frac{X-o_x}{k} + X_{\min},\qquad r_y = Y_{\max} - \frac{Y-o_y}{k}$$

$$\det = p_x q_y - q_x p_y \;=\; \lVert\mathbf{P}\rVert\,\lVert\mathbf{Q}\rVert\sin\theta,\qquad
s = \frac{r_x q_y - q_x r_y}{\det},\qquad
w = \frac{p_x r_y - r_x p_y}{\det}$$

If $|\det| < 10^{-12}$ (the two Å vectors are parallel — a degenerate cell) it returns
$(u,v)=(0,0)$ rather than throwing. Because the canvas 2D transform is set to
`setTransform(dpr,0,0,dpr,0,0)`, all mapper arithmetic is in CSS pixels, which is also what
`getBoundingClientRect()` yields — so the inversion needs no device-pixel-ratio correction.

#### 4d. What the projection is, exactly

Write the folded coordinate in the triad that is orthonormal *in fractional-component space*,
$\mathbf{x} = u\,\hat{\mathbf{u}} + v\,\hat{\mathbf{v}} + d\,\hat{\mathbf{h}}$. The true Å position
is $\mathbf{R} = u\mathbf{e}_u + v\mathbf{e}_v + d\mathbf{e}_h$. Then:

- **KDE Slice canvas** draws $u\mathbf{e}_u + v\mathbf{e}_v = \mathbf{R} - d\mathbf{e}_h$ — a parallel
  projection **along $\mathbf{e}_h$**. In-plane distances and angles are true Å (4a), and inside a thin
  slab $d$ is nearly constant, so the projection is essentially a rigid translation of the layer.
  $\mathbf{e}_h$ is *not* in general the geometric normal of the drawn plane, though: orthogonality
  holds in fractional-component space
  ($\hat{\mathbf h}\cdot\hat{\mathbf u} = \hat{\mathbf h}\cdot\hat{\mathbf v} = 0$), but in Å space
  $\mathbf{e}_u\cdot\mathbf{e}_h = \sum_{ij} \hat u_i \hat h_j\,(\mathbf{A}_i\cdot\mathbf{A}_j)$, which
  vanishes only for a metric that makes it vanish (a cubic cell, or an axis normal in an orthogonal
  cell). For a general cell this projection is oblique too — the in-plane metric stays exact
  regardless.
- **Slab In Cell canvas** draws $u\mathbf{e}_u + d\mathbf{e}_h = \mathbf{R} - v\mathbf{e}_v$ — a parallel
  projection **along the in-plane direction $\mathbf{e}_v$**. Distances along $\mathbf{e}_u$ and
  $\mathbf{e}_h$ and the angle between them are true; $\mathbf{e}_v$ is generally *not* perpendicular to
  the drawing plane in a triclinic cell, so this is an **oblique (axonometric) projection**, not an
  orthographic one.

---

### Step 5 — Slab In Cell: band geometry, cell outline, atom markers

**Code:** `StructurePage.jsx` → `drawSlab(ctx, width, height)`, invoked by an effect that sizes the
canvas ($W_{px} = \max(220, \mathrm{rect.width})$, $H_{px} = \max(260, \mathrm{rect.height})$, backing
store $\times$ `devicePixelRatio`) and stores the returned geometry in `slabGeometryRef`.

**Draw order per frame:** `clearRect` → fill the whole canvas with `--canvas-bg` → cell outline → band
fill → band stroke → **atom markers** → labels. The markers are drawn *after* the band, so the in-slab
colours sit on top of the blue tint rather than being blended with it, and the labels sit on top of
everything.

**Step 5.1 — plane setup.** $\mathbf{e}_u = \texttt{vectorFromFraction}(\hat{\mathbf{u}}, \mathbf{A})$
and $\mathbf{e}_h = \texttt{vectorFromFraction}(\hat{\mathbf{h}}, \mathbf{A})$, both in Å. The eight
cube corners are projected to plane coordinates
$(u,d)_j = (\hat{\mathbf{u}}\cdot\mathbf{c}_j,\ \hat{\mathbf{h}}\cdot\mathbf{c}_j)$, giving
$s_{\min}=u_{\min}$, $s_{\max}=u_{\max}$.

**Step 5.2 — band depth, clipped to the cell.** With $\Delta_d = d_{\max}-d_{\min}$ (or 1 if zero):

$$d_c = d_{\min} + z_c\Delta_d,\qquad
\delta = \Delta z\,\Delta_d,\qquad
d_\mathrm{start} = \max\!\big(d_{\min},\, d_c - \tfrac{\delta}{2}\big),\qquad
d_\mathrm{end} = \min\!\big(d_{\max},\, d_c + \tfrac{\delta}{2}\big)$$

The drawn band is therefore clipped at the cell faces; since every folded atom already lies inside
the cell, the clipping never hides in-slab atoms, but the on-canvas band can be visually thinner than
the nominal `d = thickness` label near $z_c\to 0$ or $1$ — and, because atoms are *not* wrapped
either, the slab genuinely samples less than its nominal thickness there (Step 6).

**Step 5.3 — bounds.** `makeProjectedPlane(U, H, [...8 cube corners, 4 band corners])`. Including the
band corners guarantees the highlighted band is always inside the fitted view; the cell **rectangle**
corners are not in the list, so they are not guaranteed to be (see the caveat in Step 4a).

**Step 5.4 — outlines.** Both polygons are rectangles **in plane coordinates**, mapped through the
oblique 2D basis, so they render as **parallelograms**:

- cell outline: $(u_{\min},d_{\min}) \to (u_{\max},d_{\min}) \to (u_{\max},d_{\max}) \to (u_{\min},d_{\max})$,
  stroked in `--border-strong`, 1 px;
- band: $(u_{\min},d_\mathrm{start}) \to (u_{\max},d_\mathrm{start}) \to (u_{\max},d_\mathrm{end}) \to (u_{\min},d_\mathrm{end})$,
  filled `rgba(79, 140, 255, 0.18)` and stroked `#74a7ff`.

For the three axis presets this rectangle **is** the exact silhouette of the unit cell projected along
$\hat{\mathbf{v}}$ (the cube's shadow along a cell axis is a full cell face). For a **custom normal**
the true silhouette of a cube projected along an arbitrary direction is a hexagon; the code draws the
axis-aligned bounding rectangle in $(u,d)$ instead — so the outline is a **bounding parallelogram,
not the exact cell cross-section**. (The KDE panel, by contrast, uses the exact polygon; see Step 9.)

**Step 5.5 — atoms and depth cueing.** The loop is

```js
const sampleLimit = Math.min(points.length, SLAB_CANVAS_MAX_POINTS); // 1e6
const stride = Math.max(1, Math.floor(points.length / sampleLimit)); // == 1 in practice
```

Since `points.length` is already capped at $10^6$ upstream, `stride` is always 1 — the canvas draws
**every** point it was given; the only subsampling is the one in Step 1e. (For an *empty* point set
`sampleLimit` is 0 and `stride` evaluates to `NaN`; harmless only because the loop body is never
entered.)

Each atom is placed at $\texttt{map}(\hat{\mathbf{u}}\cdot\mathbf{x},\ \hat{\mathbf{h}}\cdot\mathbf{x})$
and drawn with `ctx.fillRect`:

| condition | colour | marker |
| --- | --- | --- |
| `inActiveSlab(point)` | `elementColors[element]` (Step 8), else `#8A8F98` | 2 × 2 CSS px |
| otherwise | `rgba(166, 176, 188, 0.22)` | 1 × 1 CSS px |

Depth cueing is therefore **binary** (in-slab vs. out-of-slab), not a continuous depth fade, and
marker size carries no element or distance information. `fillRect` anchors the marker's *top-left*
corner at the projected point, so markers sit up to 1 px right of / below the true position — a
sub-pixel bias, irrelevant for reading the figure but present. Markers are drawn in draw order
(file order), so later atoms overwrite earlier ones; there is no depth sorting.

> **Selection here is not the KDE's selection.** `inActiveSlab()` tests the **single folded copy** of
> each atom, $|\tilde d(\mathbf x) - z_c| \le \Delta z/2$, with **no periodic images**. Both KDE
> implementations tile neighbour images within $m = \min(0.5,\max(0.1,2\,\mathrm{bw},\Delta z))$ first and then
> select, and their overlay figure `slabCount` counts *distinct source atoms* of the images that
> landed in the slab (Step 3c). So the number of highlighted markers here and the
> `N atoms in slab` figure on the KDE panel are **different quantities**, and they diverge most when
> the slab is clipped at a cell face ($z_c$ near 0 or 1): the KDE density wraps and stays correct,
> while this panel and the 3D view highlight only the part of the layer that lies inside the cell.

**Step 5.6 — labels.** Four `fillText` labels in `--text`, 12 px:

- the in-plane axis label `sliceConfig.uLabel` near $(u_{\max}, d_{\min})$, at
  $\big(\min(W_{px}-24,\ X+4),\ \min(H_{px}-8,\ Y+14)\big)$;
- the normal label near $(u_{\min}, d_{\max})$, at $\big(\max(8,\ X-12),\ \max(16,\ Y-6)\big)$;
- the two band labels at a **fixed $x = 10$** (the canvas's left edge, not the band's), vertically
  anchored to $\texttt{map}(u_{\min},\ d_\mathrm{start})$ — which, because the mapper flips $Y$, is the
  band's **lower** edge on screen: `<label>=<zCenter>` at $\max(30,\ Y-6)$ and `d=<thickness>` at
  $\min(H_{px}-16,\ Y+18)$.

The panel header separately shows the clamped interval
$[\max(0, z_c - \Delta z/2),\ \min(1, z_c + \Delta z/2)]$.

**Cost.** The whole canvas is rebuilt on every $z_c$ change — one `fillRect` per point, up to $10^6$
per frame — so dragging the band re-issues that loop on every `pointermove`. This is the practical
limit on interactivity for large models.

---

### Step 6 — Dragging the band (cursor → slice position)

**Code:** the pointer-handler effect in `StructurePage.jsx`, keyed on `[structure]` (it must not be
`[]`: the canvas only exists after a structure loads — see the note in
[AGENTS.md](../../AGENTS.md)). The band geometry published by the last `drawSlab` is read from
`slabGeometryRef`; the live slice position is mirrored into `zCenterRef` so the handlers do not need
to re-subscribe on every slider tick.

1. `planeCoordsAt(event)` = `geometry.invert(clientX − rect.left, clientY − rect.top)` → $(u,d)$.
2. `overBand` hit test, in plane coordinates:
   $u_{\min} \le u \le u_{\max}$ **and** $d_\mathrm{start} \le d \le d_\mathrm{end}$. Because the test
   is done in plane coordinates rather than on screen, it is exact for an oblique/parallelogram band.
   Note that these are the **clipped** depths of Step 5.2, so the grabbable region shrinks with the
   drawn band as $z_c\to0$ or $1$; and at the minimum thickness ($\Delta z = 0.01$ of the depth range) the
   band is only about 2 CSS px tall on a 260 px canvas, i.e. nearly ungrabbable.
3. `zCenterAt(coords)` $= (d - d_{\min})/\Delta_d$ — the inverse of Step 3c.
4. On `pointerdown` inside the band: store `offset = zCenter − zCenterAt(coords)` (grab-point
   preservation, so the band does not jump under the cursor), capture the pointer, set the
   `grabbing` cursor.
5. On `pointermove` while dragging: $z_c \leftarrow \mathrm{clamp}(\texttt{zCenterAt} + \texttt{offset},\,0,\,1)$.
   Only the **centre** is clamped — nothing keeps the slab inside the cell, and nothing wraps it. At
   $z_c=0$ or $1$ only half the nominal thickness contains folded atoms, so the highlighted set and
   the 3D cross-sections are sampled asymmetrically there (the KDE panel is not, Step 3c).
   Without a drag, the cursor is `grab` over the band and `default` elsewhere.
6. `pointerup` / `pointercancel` release the capture; `pointerleave` resets the cursor.

`touch-action: none` on the canvas (CSS) keeps touch drags from scrolling the page. Thickness is not
changed by dragging.

---

### Step 7 — The 3D "Folded Unit Cell" scene

**Code:** the large Three.js effect in `StructurePage.jsx`, keyed on
`[points, unitCell, zCenter, thickness, themeVars, sliceConfig, elementColors]` — i.e. **the whole
scene is torn down and rebuilt** whenever the slice moves. Camera state is preserved across rebuilds
(Step 7.5).

> **The effect returns early when `points.length === 0`** (`if (!mount || points.length === 0) return
> undefined;`) — reachable when the file parsed to zero atoms (e.g. a coords-only `.rmc6f` on the
> backend path, Step 1c) or when subsampling left the selected element with no sampled atoms. No scene
> is built, and because the previous run's cleanup disposes the renderer **without clearing the
> mount** (`mount.replaceChildren` runs only on a successful build), the previously rendered canvas
> stays on screen as a frozen image. `modelExportRef.current` is `null` in that state, so the 3× PNG
> export resolves to `null` and silently saves nothing, while the 1× export happily saves the stale
> canvas.

**7.1 Renderer and scene.** `THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true })`,
pixel ratio $\min(\mathrm{devicePixelRatio}, 2)$, sized to the mount's client box.
`preserveDrawingBuffer` is required so the canvas can be read back for PNG export at any time.
Background = the CSS `--canvas-bg` colour. **There are no lights in the scene** — every material used
is unlit (`PointsMaterial`, `MeshBasicMaterial`, `LineBasicMaterial`), so there is no shading, no
specular highlight, and no ambient-occlusion depth cue.

**7.2 Atoms.** Points are grouped by element, and each group becomes one `THREE.Points` object with a
plain position buffer:

$$\mathbf{q} = \sum_i x_i\,\mathbf{b}_i \;-\; \mathbf{c}
\qquad (\mathbf{b}_i,\mathbf{c}\ \text{from Step 2})$$

Material: `PointsMaterial({ color: elementColors[element] ?? '#8A8F98', size: 0.018,
sizeAttenuation: true })`.

Two properties of that buffer are worth stating: positions are written into a **`Float32Array`**, so
every 3D coordinate is quantized to single precision ($\sim10^{-7}$ relative, $\approx10^{-6}$ Å for a
10 Å cell) while the 2D canvases stay in double precision; and **no subsampling happens here** — every
element-filtered point, up to the $10^6$ cap of Step 1e, is uploaded, one buffer per element.

> **These are not spheres.** Atoms are camera-facing square point sprites of world size **0.018
> normalized cell units** — 1.8 % of the longest unit-cell edge, i.e. $\approx 0.19$ Å for a 10.4 Å
> cell — scaled with distance by `sizeAttenuation`. The size is a single constant: it encodes **no**
> ionic/covalent/van-der-Waals radius and does not vary by element. Only the colour is
> element-specific.

> **There is no bonding.** The scene contains no bond cylinders, no neighbour search, and no
> distance criterion of any kind. Nothing in `StructurePage.jsx` computes interatomic distances.

**7.3 Cell frame.** `CUBE_CORNERS` (8 fractional corners) → `cellCorners(basis, center)` →
$\mathbf{b}$-space, origin-centred. `makeCellEdgeGeometry` emits the 12 edges of `CUBE_EDGES` as line
segments; material `LineBasicMaterial({ color: '#737c86' })`. For a triclinic cell this is the correct
oblique parallelepiped, because the edges are built from the actual $\mathbf{A}_i$ (normalized), not
from an axis-aligned box.

**7.4 Slab cross-sections.** The slab is drawn as its two bounding **cross-section polygons**, not as
a solid. `planeSectionVertices(normal, offset)` clips the unit cube with the plane
$\hat{\mathbf{h}}\cdot\mathbf{x} = \mathrm{offset}$:

For each of the 12 cube edges $(\mathbf{p}_0,\mathbf{p}_1)$, with
$d_j = \hat{\mathbf{h}}\cdot\mathbf{p}_j - \mathrm{offset}$:

- if $|d_j| \le 10^{-9}$, the corner itself is a vertex;
- if $d_0 d_1 < 0$, the crossing point is
  $\mathbf{p}_0 + \dfrac{d_0}{d_0-d_1}\big(\mathbf{p}_1-\mathbf{p}_0\big)$.

Vertices closer than $10^{-8}$ are merged; fewer than 3 unique vertices returns `[]` (no polygon).
The survivors are sorted about their centroid by
$\operatorname{atan2}(\boldsymbol{\delta}\cdot\hat{\mathbf{v}},\,\boldsymbol{\delta}\cdot\hat{\mathbf{u}})$,
giving a convex polygon in angular order (3–6 vertices).

Three copies of this routine exist: `StructurePage.jsx` → `planeSectionVertices()`,
`localKdeWorker.js` → `planeSectionVertices()`, and `rmc_toolkits/kde.py` →
`_plane_section_vertices()`. **The two JS copies are identical.** Python matches them on the clipping
algorithm, on the $10^{-9}$ on-plane test and on the $10^{-8}$ dedupe — the vertex *set* is the same —
but it sorts with a **different in-plane basis**: `_plane_basis(normal)` seeds from
`np.eye(3)[argmin(|h|)]` whereas both JS copies use `makeFreePlaneBasis()` (Step 3b). For a custom
normal the two bases differ (for $\hat{\mathbf h}\propto[1,1,0]$ the Python $(\hat{\mathbf u},\hat{\mathbf v})$
equals $(-\hat{\mathbf v},\hat{\mathbf u})_\mathrm{JS}$, a 90° rotation of the sort angle), so the
returned sequence is a **cyclic rotation** of the JS order: the same polygon, a different starting
vertex. For the `a`/`b`/`c` presets the bases coincide and the orders match.

The two sections at $d_\mathrm{start}$ and $d_\mathrm{end}$ (Step 5.2) are transformed into
$\mathbf{b}$-space and fed to:

- `makeSlabGeometry(sections)` — a **triangle fan per section** (`0, i, i+1`), then
  `computeVertexNormals()`. Material `MeshBasicMaterial({ color: '#4f8cff', opacity: 0.12,
  transparent: true, side: DoubleSide, depthWrite: false })`. **Only the two caps are filled — the
  side walls of the slab are never triangulated**, so the "slab" is two translucent lids, not a closed
  prism.
- `makeSectionEdgeGeometry(sections)` — each section's closed boundary loop, plus vertex-to-vertex
  "rungs" between the two sections **only when both polygons have the same vertex count**. When the
  slab straddles a cell corner the two cross-sections can have different vertex counts (e.g. 4 and 5)
  and the connecting rungs silently disappear. Material `LineBasicMaterial({ color: '#8c96a3',
  opacity: 0.95 })`.

**7.5 Camera and framing.** `PerspectiveCamera(45°, W/H, 0.01, 20)`, then re-normalized to the
geometry. The framing radius is built in two steps — a `THREE.Box3` fitted to the 8 cell corners, then
*that box's* bounding sphere:

$$R_{sph} = \max\big(\text{radius of the bounding sphere of the axis-aligned bounding box of the 8 cell corners},\ 0.5\big)$$

$$\texttt{near} = R_{sph}/100,\quad \texttt{far} = 20R_{sph},\quad
\texttt{minDistance} = 0.35R_{sph},\quad \texttt{maxDistance} = 8R_{sph}$$

For an oblique (or merely non-axis-aligned) cell the AABB step makes $R_{sph}$ strictly larger than
the minimal bounding sphere of the corners, so the framing is more conservative than a tight fit.

For a first build the camera is placed at $\mathbf{s}_\mathrm{centre} + R_{sph}\,(1.7,\,1.45,\,1.55)$ (an
off-axis three-quarter view; $\mathbf{s}_\mathrm{centre}\approx\mathbf{0}$ because the corners were
centred in Step 2) with the orbit target at the sphere centre. On any later rebuild the previously
saved `{ position, target, zoom }` from `cameraStateRef` is restored, so moving the slice slider does
not throw away the user's viewpoint. `OrbitControls` with `enableDamping: true`,
`dampingFactor: 0.08`, `enablePan: true`; a `requestAnimationFrame` loop calls `controls.update()` and
renders every frame.

**7.6 Resizing.** A `ResizeObserver` on the mount updates `camera.aspect`, calls
`updateProjectionMatrix()`, and `renderer.setSize(w, h)` — this is what recovers a scene first built
at 0 × 0 (measured before layout settles). The 2D canvases have the analogous guard: a
`ResizeObserver` bumps `sizeTick`, which re-runs the draw effects at the settled size.

**7.7 Teardown.** The cleanup disconnects the observer, snapshots the camera state, cancels the
animation frame, and disposes controls, renderer, all geometries and materials (including a
`group.traverse` over the per-element point clouds). It does **not** remove the canvas from the mount
— see the early-return note above.

---

### Step 8 — Element colours (shared by the slab canvas, the 3D view and the legend)

**Code:** [atomColors.js](../../web_app/frontend/src/atomColors.js) → `buildElementColors(elements)`,
called from a memo keyed on `[structure]` alone:

```js
elementColors = buildElementColors(structure.elements?.length ? structure.elements
                                                             : structure.points.map(p => p.element));
```

1. Unique element symbols are **sorted** first, so the assignment is deterministic and independent of
   atom order in the file.
2. Each element takes its colour from `ELEMENT_COLORS`, a CPK/Jmol-style table of **55 entries
   running H → Bi**. Everything above Bi is absent — that includes **all** actinides as well as Po,
   At and Rn — as are Kr, Sc, Tc, Ru, Rh, Pd, Xe, the lanthanides, Hf, Ta, Re, Os, Ir and Tl.
3. If the element is absent from the table **or its table colour is already taken**, the next unused
   entry of `FALLBACK_PALETTE` (16 qualitative colours) is used.
4. If the fallback palette is exhausted, an evenly spaced HSL hue
   `hsl((n·47) mod 360, 70%, 60%)` is generated, where `n` is the number of elements already assigned.
5. Any lookup miss at draw time falls back to `DEFAULT_ELEMENT_COLOR = '#8A8F98'`.

The invariant is: **no two elements in one structure ever share a colour**, and the same map is used
by the Slab In Cell markers, the 3D point clouds, and the legend rendered under the 3D panel.

Two consequences of the memo depending on `structure` only:

- the map is **invariant under the element filter** — every element keeps its colour when you filter,
  and the legend under the 3D panel always lists **every element in the file**, not just the displayed
  one;
- the `points`-based input is only a fallback for a payload with an empty `elements` list, but when it
  is taken it walks the full point array (up to $10^6$ entries) on every structure change.

Because the element symbol is the dictionary key, the Python/JS capitalization difference in Step 1c
can hand the two data paths different keys — and therefore, potentially, different colours — for the
same file.

---

### Step 9 — The KDE Slice canvas (projection only)

The density estimator belongs to the KDE section; what this page adds is the same projection
machinery, plus a set of fallbacks worth knowing when a panel looks empty. **Code:**
`StructurePage.jsx` → `drawKdeSlice(ctx, width, height)`.

1. `uVector`/`vVector` come from the KDE result when available (`kde.uVector || sliceConfig.u`), so
   backend-path custom slices are drawn in the server's basis (see the Step 3b warning).
2. `planePolygon` — the **exact** cross-section polygon of the unit cell at the slice centre depth
   (computed by `_plane_section_vertices` / `planeSectionVertices` from Step 7.4, expressed in
   $(u,v)$) — is used both as the drawn cell outline and, via `makeProjectedPlane`, to set the fitted
   bounds. With no KDE result the extent falls back to `[-0.5, 0.5, -0.5, 0.5]` and, if
   `planePolygon` is also absent, the outline degrades to the rectangle
   $[x_{\min},y_{\min}]\ldots[x_{\max},y_{\max}]$ of that extent. Unlike the Slab panel, the polygon
   outline is exact for oblique cells and custom normals, and its shape changes as the slice moves.
3. The heatmap is painted **only when `density && grid > 0 && kde.vmax > kde.vmin`**. The density grid
   is rendered into an offscreen `grid × grid` `ImageData` through the selected colormap LUT
   (normalized as $(\rho - v_{\min})/\mathrm{span}$ with $\mathrm{span} = v_{\max}-v_{\min}$ or 1 if that
   is 0, then clamped to a 0–255 LUT index), then blitted with an affine `ctx.transform` whose columns
   are $\texttt{map}(x_{\max},y_{\min}) - \texttt{map}(x_{\min},y_{\min})$ and
   $\texttt{map}(x_{\min},y_{\max}) - \texttt{map}(x_{\min},y_{\min})$ — i.e. the unit image square is
   mapped onto the parallelogram that the extent box occupies on screen — with `ctx.clip()` set to the
   cell polygon so only the in-cell part shows. `imageSmoothingEnabled = true`, so **the displayed
   field is a bilinear interpolation between grid cells**: visible structure finer than the grid pitch
   is interpolation, not data.
4. Contour polylines are mapped point-by-point through the same `mapper.map` (1 px, `themeVars.contour`).
   They are drawn inside the same branch as the heatmap, and they are *not* clipped to the cell polygon.
5. When the gate in (3) fails, the panel instead shows a single placeholder string in `--muted`,
   `500 13px Inter`, at $(14, 28)$: `Computing KDE...` while a request is in flight, otherwise
   `No atoms in this slab`.
6. Overlay text (drawn with a dark stroke `rgba(13, 18, 28, 0.62)`, `lineWidth 3`, under a white fill,
   so it stays legible over any colormap) reports `<slabCount> atoms in slab (fit <fitCount>)` at
   $(12,22)$, `<label>=<center>  d=<thickness>  bw=<bw>` at $(12,40)$, and `log10 density` at
   $(12,58)$ when the log toggle is on. Recall from Step 3c/5.5 that `slabCount` counts unique source
   atoms of the **periodic-image-augmented** set, not the markers highlighted in the Slab panel.

**Fixed draw order:** background fill → heatmap (clipped) → contours → cell outline → overlay text, so
the outline and the text always sit on top of the density.

---

### Step 10 — Export / screenshot rendering

**Code:** `save2dPanel()`, `captureModelBlob()`, `saveKdeSlice/saveSlab/saveModel` in
`StructurePage.jsx`; [figureExport.js](../../web_app/frontend/src/figureExport.js) →
`saveCanvasAsPng()`, `canvasToPngBlob()`, `downloadBlob()`, `sanitizeFilename()`.

Both 2D panels and the 3D panel offer exactly two options (`PANEL_SAVE_OPTIONS`): `png` (labelled 1×)
and `png3x` (labelled 3×).

- **2D, 1×** — `canvas.toBlob('image/png')` on the live canvas, whose backing store is
  $W_{px}\!\cdot\!\mathrm{dpr} \times H_{px}\!\cdot\!\mathrm{dpr}$ with an **uncapped**
  `devicePixelRatio`.
- **2D, 3×** — a fresh offscreen canvas of exactly $3W_{px} \times 3H_{px}$ device pixels with
  `setTransform(3,0,0,3,0,0)`, then the **same draw function** is re-run
  (`drawKdeSlice` / `drawSlab`). Because the draw code works in CSS-pixel units, every element —
  including the 2 px atom markers, the 18 px padding and the 12 px labels — scales proportionally;
  the result is genuinely higher-resolution, not an upscaled bitmap. Minimum logical sizes are
  320 × 260 (KDE) and 220 × 260 (slab) CSS px.
- **3D, 1×** — reads the live WebGL canvas out of the mount (possible only because of
  `preserveDrawingBuffer: true`), whose backing store is
  $W_{px}\!\cdot\!\min(\mathrm{dpr},2) \times H_{px}\!\cdot\!\min(\mathrm{dpr},2)$.
- **3D, 3×** — `renderer.setPixelRatio(3)`, `setSize(width, height, false)` (style untouched),
  re-render, `toBlob`, then restore the previous pixel ratio, size and frame → exactly
  $3W_{px} \times 3H_{px}$.

> **"1×" does not mean one device pixel per CSS pixel, and the four options give four different
> scales.** On a 2× display the 2D "1×" file is already 2× and the 3D "1×" file is 2× (capped); the
> "3×" options are then only a 1.5× increase in linear resolution over what was on screen. On a 3×
> phone display the 2D "1×" export is 3× while the 3D "1×" export is still 2×.

File names are `KDE_Slice_<normal label>.png`, `Slab_In_Cell_<normal label>.png`, and
`Folded_Unit_Cell.png`, passed through `sanitizeFilename()`, which (i) unwraps LaTeX-style
superscripts `^{…}` → `…`, (ii) collapses every run of characters outside `[A-Za-z0-9_.-]` to a
single `_`, (iii) strips leading/trailing `_`, and (iv) returns the literal `figure` if nothing
survives.

The normal label for a custom direction is built with
`Number(v).toLocaleString(undefined, { maximumFractionDigits: 2 })`, so both the on-canvas label and
the filename are **locale-dependent**: `KDE_Slice_[1 1 0]` is written as `KDE_Slice__1_1_0.png` (note
the doubled underscore — the `_` before `[` is itself a word character and is not absorbed into the
run), a direction component of 1.5 gives `Slab_In_Cell__1.5_1_0.png` in an en-US locale but
`Slab_In_Cell__1_5_1_0.png` where the decimal separator is a comma.

---

### Parameters and defaults

| Parameter | Code location | Default | Range / values | Units |
| --- | --- | --- | --- | --- |
| `selectedElement` | `StructurePage.jsx` state | `all` | `all` + elements found in the file | — |
| `sliceDirection` | state / `NORMAL_OPTIONS` | `c` | `a`, `b`, `c`, `custom` | — |
| `customDirection` | state | `[1, 1, 0]` | any 3 reals (number inputs, step 0.1); all-zero ⇒ `[0,0,1]` | components on the fractional axes |
| `zCenter` ($z_c$) | state + auto-centre effect | 0.5, then the densest of 50 depth bins | 0–1, slider step 0.001 | fraction of $\Delta_d$ |
| `thickness` ($\Delta z$) | state | 0.08 | 0.01–0.5, step 0.01 | fraction of $\Delta_d$ |
| `bandwidth` | state | 0.03 | 0.005–0.15, step 0.005 | SciPy `bw_method` factor (dimensionless) |
| `gridSize` | state | 120 | 80 / 120 / 160 / 220 | grid cells per axis |
| `colormap` | state | `viridis` | `COLORMAP_NAMES` | — |
| `showContours` / `logScale` | state | on / on | boolean | — |
| KDE contour levels | request param `levels` | 8 | — | — |
| Periodic-image margin (KDE only) | `kde.py`, `localKdeWorker.js` | $\min(0.5,\max(0.1,2\,\mathrm{bw},\Delta z))$ | — | fractional cell |
| KDE fit-point cap | `MAX_KDE_FIT_POINTS` / `fitLimit` | 6000 | — | atoms (incl. images) |
| `STRUCTURE_MAX_POINTS` | `StructurePage.jsx` | 1 000 000 | — | atoms |
| `MAX_STRUCTURE_POINTS` | `app.py` | 1 000 000 | request clamped to [100, 10⁶] | atoms |
| `SLAB_CANVAS_MAX_POINTS` | `StructurePage.jsx` | 1 000 000 | — | atoms (never binding in practice) |
| Auto-centre bins | `StructurePage.jsx` | 50 | bin clamped to [0, 49]; ties → lowest bin | depth bins |
| `padding` | `makePlaneMapper` | 18 | — | CSS px |
| Canvas minimum size | draw effects | 320×260 (KDE), 220×260 (slab) | — | CSS px |
| Backing resolution (2D) | draw effects | `devicePixelRatio`, uncapped | — | device px / CSS px |
| `--panel-aspect` | `slicePanelGeometry` | `planeAspect` (KDE), `sideAspect` (slab), `max(planeAspect, 1)` (3D) | — | — |
| In-slab / out-of-slab marker | `drawSlab` | 2×2 / 1×1 | — | CSS px |
| Out-of-slab colour | `drawSlab` | `rgba(166,176,188,0.22)` | — | — |
| Band fill / stroke | `drawSlab` | `rgba(79,140,255,0.18)` / `#74a7ff` | — | — |
| KDE placeholder text | `drawKdeSlice` | `Computing KDE...` / `No atoms in this slab` at (14, 28), `500 13px Inter`, `--muted` | — | CSS px |
| 3D point size | `PointsMaterial` | 0.018, `sizeAttenuation: true` | — | normalized cell units ($\ell_{\max}=1$) |
| 3D position precision | `Float32Array` | single precision | ~10⁻⁷ relative | ≈10⁻⁶ Å for a 10 Å cell |
| 3D cell edge / slab edge / slab face colour | Three.js materials | `#737c86` / `#8c96a3` (α 0.95) / `#4f8cff` (α 0.12) | — | — |
| Camera FOV | `PerspectiveCamera` | 45° | — | degrees |
| Camera near / far | after framing | $R_{sph}/100$ / $20R_{sph}$ | $R_{sph}\ge0.5$, from the AABB of the corners | normalized cell units |
| Camera offset | first build | $R_{sph}\,(1.7, 1.45, 1.55)$ | — | normalized cell units |
| Orbit min / max distance | `OrbitControls` | $0.35R_{sph}$ / $8R_{sph}$ | — | normalized cell units |
| Orbit damping | `OrbitControls` | 0.08 | — | — |
| Renderer pixel ratio | `WebGLRenderer` | $\min(\mathrm{dpr},2)$; 3 during export | — | — |
| Export pixel dimensions | `save2dPanel` / `captureModelBlob` | 2D 1×: $W\!\cdot\!\mathrm{dpr}\times H\!\cdot\!\mathrm{dpr}$; 2D 3×: $3W\times3H$; 3D 1×: $W\!\cdot\!\min(\mathrm{dpr},2)\times\ldots$; 3D 3×: $3W\times3H$ | — | device px |
| Plane-section tolerances | `planeSectionVertices` | $10^{-9}$ on-plane, $10^{-8}$ dedupe | — | fractional |
| `invert()` determinant floor | `makePlaneMapper` | $10^{-12}$ | — | Å² |
| `normalize()` zero floor | `StructurePage.jsx` | $10^{-9}$ | — | — |
| Degenerate-span guards | `pointDepth`, `makePlaneMapper`, `unitCell`, `drawKdeSlice` | depth span `\|\| 1`; $\Delta X,\Delta Y$ `\|\| 1`; $\ell_{\max}=\max(\ell_i,10^{-9})$; $N_i \ge 10^{-12}$; density span `\|\| 1` | `stride` in `drawSlab` is `NaN` for an empty point set (loop skipped) | — |
| Debounce before KDE request | effects | 80 ms (browser worker) / 160 ms (backend) | backend request aborted via `AbortController` | ms |
| Frac file coordinate precision | `frac_lines_from_rmc6f` | 5 decimals of a box fraction | — | ≈10⁻³ Å for a 104 Å box |

---

### Caveats / what this is not

- **Every panel superimposes all supercell cells.** The "Folded Unit Cell" and both 2D panels show
  $N_1N_2N_3$ overlaid copies of the cell, not one physical cell. Cluster sizes, apparent site
  splitting, and marker density are ensemble properties.
- **The highlighted-marker count and the KDE's `N atoms in slab` are different quantities.** The Slab
  In Cell markers, `inActiveSlab()` and the 3D cross-sections test the single folded copy of each
  atom; both KDE implementations tile periodic images within
  $\min(0.5,\max(0.1,2\,\mathrm{bw},\Delta z))$ first and report unique *source* atoms. The two figures diverge
  most when the slab is clipped at a cell face, where the KDE density wraps correctly and the other
  panels show only the in-cell half of the layer.
- **The 3D view has no bonds and no chemically meaningful atom radii.** Atoms are 0.018-unit square
  point sprites in a fixed size for every element; there is no neighbour search, bond-length
  criterion, or coordination analysis anywhere in `StructurePage.jsx`. Do not read coordination
  polyhedra off this view.
- **The 3D scene is unlit.** All materials are `Basic`/`Points`/`Line` — flat colour, no shading, no
  depth-based fading. Perceived depth comes only from perspective and `sizeAttenuation`.
- **The 3D "slab" is two lids, not a solid.** `makeSlabGeometry` fills only the two cross-section
  polygons; the side walls are never triangulated, and the connecting edge rungs are drawn only when
  the two cross-sections happen to have the same vertex count.
- **An empty point set freezes the 3D panel instead of clearing it.** The Three.js effect returns
  early when `points.length === 0` and the previous teardown does not clear the mount, so a stale
  frame stays visible; the 3× export silently produces nothing in that state.
- **Slab In Cell depth cueing is binary.** In-slab vs. out-of-slab, 2 px vs. 1 px. There is no
  continuous depth ramp, and no z-ordering — later atoms simply overwrite earlier ones.
- **The Slab In Cell cell outline is a bounding parallelogram for custom normals**, and it can be
  drawn outside the fitted view. It is the exact projected cell face for the `a`/`b`/`c` presets, but
  for an arbitrary normal it is the bounding box in $(u,d)$ whereas the true silhouette is
  generally a hexagon — and because its corners are not in the list `makeProjectedPlane` fits, part of
  it can fall outside the 18 px padded box and be clipped by the canvas edge on an oblique cell. The
  KDE panel's outline (the exact section polygon) does not have either limitation.
- **Both 2D projections are oblique for a general cell.** The Slab panel projects along the in-plane
  direction $\mathbf{e}_v$, and the KDE panel along $\mathbf{e}_h$, which is the Å image of the fractional
  normal and not the plane's geometric normal unless the cell metric makes
  $\mathbf U\cdot\mathbf H = \mathbf V\cdot\mathbf H = 0$. In both panels the *in-plane* lengths and
  angles are true Å (the Gram matrix is preserved exactly); a length measured across the figure in a
  general direction is not.
- **The two runtimes pick different in-plane axes for a custom normal.** See Step 3b and 7.4. The
  consequences are a relative rotation between the KDE panel and the Slab panel on the backend path, a
  cyclic rotation of the section-polygon vertex order between Python and JS, and a possible mismatch
  between the CSS panel aspect ratio (always computed from the local basis and the cube corners) and
  the drawn content (letterboxed, never distorted, because the mapper fits isotropically).
- **The `b` preset triad is left-handed** ($\hat{\mathbf{u}}\times\hat{\mathbf{v}}=-\hat{\mathbf{h}}$),
  so the `b`-normal view is mirrored relative to a right-handed convention. Both runtimes share this,
  so they agree with each other but not with a right-handed drawing.
- **Changing the element filter or the normal resets the slice position** to the densest depth bin of
  the *selected species*, discarding a hand-placed or dragged slice.
- **The band drag only clamps the centre.** The slab is never wrapped, so at $z_c = 0$ or $1$ only
  half its nominal thickness contains atoms; and because the hit test uses the clipped band, the
  grabbable region shrinks near the cell faces and is ~2 px tall at $\Delta z = 0.01$.
- **Entering an all-zero custom direction silently gives the c-axis slice**, while the panel label
  still reads `[0 0 0]`.
- **Slice position and thickness are fractions of the projection range, not of a cell edge.** For a
  custom normal the range is longer than 1 (e.g. $\sqrt2$ for $[1,1,0]$), so a thickness of 0.08 is
  $0.08\sqrt2$ in depth units, and its physical value in Å depends on the cell.
- **The KDE heatmap is bilinearly interpolated for display** (`imageSmoothingEnabled = true`), so
  detail finer than the `grid × grid` pitch is an artefact of the blit, not of the estimate.
- **Which `.rmc6f` is read is a heuristic** (Step 1a), and the two runtimes tie-break differently, so
  a folder holding several models can show different models in the two paths.
- **Legacy-file support differs between the two paths.** Coords-only (5–6 field) `.rmc6f` files parse
  in the browser and yield zero atoms through the Flask API; element capitalization differs, and it is
  not even uniform inside one backend payload — `atomIndices` keeps the file's raw casing while
  `elements`/`elementCounts` are capitalized, which can zero the per-element site counts shown above
  the panels (Step 1c).
- **Subsampling strategies differ between the two paths** (site-stratified vs. global stride; only the
  global stride can drop a rare element entirely). Both are inert at the default 10⁶-atom cap for
  typical models. The 3D view adds no subsampling of its own and uploads every filtered point.
- **`Frac*.txt` quantizes coordinates to 5 decimals of a box fraction** (~10⁻³ Å for the sample cell).
  The web page never uses that path — it reads the `.rmc6f` directly — but package/CLI users of
  `read_structure()` do.
- **There is no live theming.** `App.jsx` renders `<StructurePage … theme="light" />`, a constant, and
  `themeVars` is a `useMemo` keyed only on `theme` — so `--canvas-bg`, `--text`, `--muted` and
  `--border-strong` are read from the document **once at mount** and never refreshed. Nothing in the
  app ever sets `data-theme`, so the `:root[data-theme='dark'|'light']` blocks in `index.css` are
  inert and the base (dark) `:root` block always wins, while `themeVars.contour` is pinned to the
  light value `rgba(21, 34, 50, 0.72)`. In practice the canvases are always drawn on the dark
  `--canvas-bg` (`#10141a`) with dark contour lines — which are dark-on-dark wherever they fall
  outside the coloured heatmap.
- **Dragging the band redraws everything.** One $z_c$ change re-issues the KDE request, redraws both
  2D canvases (up to 10⁶ `fillRect` calls) and rebuilds the entire Three.js scene, on every
  `pointermove`.
