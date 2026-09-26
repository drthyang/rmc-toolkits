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
| **Browser-loaded run** (`localRun`) | hand-written kernel sum in [`web_app/frontend/src/workers/localKdeWorker.js`](../../web_app/frontend/src/workers/localKdeWorker.js), optionally on the GPU via [`gpuKde.js`](../../web_app/frontend/src/workers/gpuKde.js) | **Visualization path.** Same estimator and the same kernel $\mathbf{H}=f^2\mathbf{C}$ (parity-tested against Python goldens to $10^{-6}$ of the peak on slabs below the 6000-point fit cap), but float32 arithmetic on the GPU branch, a different random subsample above the cap, and a cruder contour tracer. |

**The discriminator is *not* static-vs-Flask mode.** `StructurePage.jsx` branches on
`const isLocalStructure = Boolean(localRun)` and nothing else; `isStaticMode()` is used on this page
only to emit the *"Open a run folder to view the structure."* message when there is no `localRun`.
The **Demo** button in the app header (`App.jsx` → `handleToggleDemo`) is rendered unconditionally,
outside any `staticMode` branch, and it sets `localRun`. So a Flask-mode session that loads the
bundled demo run gets the **browser worker**, not `/api/kde/slice`, even though the backend is
available. `/api/kde/slice` is used only when a server-side directory is the active run source.
Static mode (GitHub Pages / no backend) is always a browser-loaded run, so it always uses the
worker.

> **Read blob shapes with care.** The smoothing kernel is $f^2\times$ the covariance of the slab's
> atoms, so its shape follows how the slab's *sites* are laid out, not how the atoms move: an
> isotropic Ga site in the element-filtered GaNb₄Se₈ layer is drawn 2 : 1 elongated at the defaults.
> The map prints the kernel's $\sigma$ in Å and flags sub-grid and strongly anisotropic kernels. See
> [*The kernel's shape follows the slab's site layout*](#the-kernels-shape-follows-the-slabs-site-layout-read-this-before-reading-blob-shapes).

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
| $\mathbf{h} = (h,k,l)$ | Miller indices of the sliced plane family, entered as three numbers (*Plane (h k l)*) | dimensionless |
| $\hat{\mathbf{h}}$ | $\mathbf{h}/\lVert\mathbf{h}\rVert_2$ | dimensionless |
| $\hat{\mathbf{u}},\hat{\mathbf{v}}$ | in-plane axes, orthonormal *in fractional space* | dimensionless |
| $d_i = \mathbf{x}_i\!\cdot\!\hat{\mathbf{h}}$ | depth of atom $i$ along the slice normal | dimensionless |
| $[d_{\min},d_{\max}]$, $\Delta_d = d_{\max}-d_{\min}$ | projection range / depth span of the unit cube along $\hat{\mathbf{h}}$ | dimensionless |
| $z_c$, $\Delta z$ | slider "Slice" and "Thickness" — **fractions of the depth span**, not Å | dimensionless |
| $\mathbf{p}_i=(u_i,v_i)$ | in-plane projection $(\mathbf{x}_i\!\cdot\!\hat{\mathbf{u}},\, \mathbf{x}_i\!\cdot\!\hat{\mathbf{v}})$ | dimensionless |
| $N_\mathrm{src}$ | unique source atoms contributing to the slab (`slabCount`) — see the warning in Step 4 | count |
| $N_\mathrm{img}$ | slab rows including periodic images | count |
| $n$ | slab rows (images included, subsampled to ≤6000) the kernel is summed over (`fitCount`) | count |
| $m$ | periodic-image margin | fractional |
| $\mathbf{C}$ | $2\times2$ sample covariance of the slab's **source atoms** (one row per atom, images excluded) | fractional² |
| $f$ | bandwidth factor (`bw`) | dimensionless |
| $\mathbf{H} = f^2\mathbf{C}$ | kernel bandwidth (covariance) matrix | fractional² |
| $\rho(u,v)$ | estimated density | per unit fractional area of the slice plane |

---

### Step 1 — Load atom positions and fold the supercell into one unit cell

**Inputs.** One `.rmc6f` file found in the selected run folder, plus an optional element filter
(`all` by default).

#### Which `.rmc6f` file, when the folder holds several

Both runtimes implement the same rule, and it is not "the only one there":

* `rmc_toolkits/parsers.py` → `find_run_configuration()` is the one Python rule, used by the
  backend (`web_app/backend/app.py` → `_find_rmc6f()`, which only adds "if the path *is* a
  `.rmc6f`, use it") and by the `rmc-triplets` CLI. The directory is globbed for `*.rmc6f` (sorted
  by name); a 0-byte or marker-less candidate (`rmc6f_problem()`, e.g. left by a killed run) never
  hides a usable one. Every file in the directory is tested against `run_stem_from_output_name()`,
  a table of output-file patterns with an explicit priority: **0** = `<stem>-NN.log`; **1** =
  `<stem>-EXAFS-*_[QR]_OUTPUT.csv`, `<stem>_FT_XFQ*.csv`, `<stem>_[FS]Q*.csv`, `<stem>_bragg*.csv`,
  `<stem>_PDF*.csv`; **2** = `Frac_coord_<stem>.txt`. Matches are sorted by `(priority, lowercased
  filename)` and the first `.rmc6f` whose stem matches wins. If nothing matches, **the first usable
  `.rmc6f` by name is used**.
* `web_app/frontend/src/browserData.js` → `chooseStructureFile()` is the parallel implementation,
  keyed on `directory + '/' + stem` so the match is per-subfolder, with the same priority table and
  the same tie-breaking: candidates sorted by path and outputs by (priority, lower-cased name,
  stem), all in code-point order as Python's `sorted()` compares strings (since 1.0; before it the
  browser fell back to directory-walk order and ranked outputs with `localeCompare`, so a
  multi-model folder could resolve to different files on the two paths).

**Consequence.** In a multi-configuration folder the page analyses the configuration the run's
outputs name, which may not be the one the user has in mind. The chosen path is reported (`source`
in the Flask payload, `source` in the browser structure object) but is not shown prominently in the
UI.

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

Every reader folds $\mathbf{f}_i$ directly — `rmc_toolkits/kde.py` → `load_unit_cell_positions()`,
`web_app/backend/app.py` → `structure()` and the browser parser
([`browserData.js`](../../web_app/frontend/src/browserData.js) → `structureFromRmc6f()`) — which
matters because the oldest `.rmc6f` variants carry no per-atom cell index at all. (Subtracting the
cell index first, `coords - cell_indices/supercell`, as `structure()` did before 1.0 and the
`Frac*.txt` writer still does, only removes an integer before the modulo and gives the same
$\mathbf{x}_i$.)

`load_unit_cell_positions()` also returns Cartesian positions
$\mathbf{x}_i^{\mathrm{cart}} = x_i\mathbf{a} + y_i\mathbf{b} + z_i\mathbf{c}$ (Å) and
`cell_lengths` $=(\lVert\mathbf{a}\rVert, \lVert\mathbf{b}\rVert, \lVert\mathbf{c}\rVert)$ (Å).
**The Cartesian array is not used by the KDE endpoint** — see the units gotcha in Step 4.

**Atom-line parsing.** Both runtimes share one anchored, validated grammar (since 1.0):
`rmc_toolkits/parsers.py` → `classify_rmc6f_atom_line()` / `iter_rmc6f_atoms()` /
`parse_rmc6f_atoms()` and their browser twins in `web_app/frontend/src/rmc6f.js` →
`classifyAtomLine()` / `parseRmc6fAtoms()`. A line is `id element [label]` followed by **exactly 7**
data fields (`x y z ref cx cy cz`, the full layout) or **exactly 3** (`x y z`, the legacy
coordinates-only form, with `referenceNumber`/`cellIndices` null); the label is a bracket group or
one non-numeric token. Numbers accept Fortran `D` exponents, `ref` must be a positive integer and
each cell index an integer in $[0, N_i)$. Any spelling of the `Atoms` marker is accepted, and
bare-CR files parse. A line of valid layout with a non-finite coordinate (`NaN`, `Inf`, Fortran
`****`) is **skipped and counted**; any other line is counted as unparsed. The parse report
(`Rmc6fParseReport` / `report`) compares the accepted count with the header's `Number of atoms:`
and its `warning()` / `parseWarning` text ("… of … atom lines unparsed", "parsed M of N atoms
declared in the header", "… skipped for non-finite coordinates") reaches the Model information
card, the PCA pages and the bond-angle payload. Before 1.0 both parsers indexed fields from the
end of the line, so an extra trailing field silently shifted every column in the browser, and
Python had no coordinates-only branch. `iter_rmc6f_atoms()` yields full-layout atoms by default;
`load_unit_cell_positions()` passes `include_coords_only=True`, since the KDE needs only element
and position.

#### Population handed to the KDE (this differs between runtimes)

* **Server-side run source** — `/api/kde/slice` calls `load_unit_cell_positions()` itself (memoized
  by `_cached_positions` → `_POSITIONS_CACHE`, a `_FileCache(16)` keyed on the file signature
  `(st_mtime_ns, st_ctime_ns, st_size, st_ino)` and the element). **Every atom**
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
`rmc_toolkits/parsers.py` → `find_run_configuration()`, `run_stem_from_output_name()`,
`classify_rmc6f_atom_line()`;
`web_app/backend/app.py` → `_find_rmc6f()`, `_sample_atoms_by_site()`, `structure()`;
`web_app/frontend/src/browserData.js` → `chooseStructureFile()`, `structureFromRmc6f()`;
`web_app/frontend/src/rmc6f.js` → `classifyAtomLine()`, `parseRmc6fAtoms()`;
`web_app/frontend/src/workers/localStructureWorker.js`;
`StructurePage.jsx` → the `points` `useMemo`.

---

### Step 2 — Choose the slice plane: normal, in-plane axes, depth coordinate

**Inputs.** The *Normal* control (`a`, `b`, `c`, or *Plane (hkl)* with three numbers labelled
*h*, *k*, *l*, default $(1\,1\,0)$); default `c`.

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
family of lattice planes with **Miller indices $(hkl)$**. So the custom input is a Miller index
triple, not a real-space vector; preset `a` slices the $(100)$ planes, and $(1\,1\,0)$ slices the
$(110)$ planes, for any cell metric including triclinic. The slab's real-space normal is
$h\mathbf{a}^*+k\mathbf{b}^*+l\mathbf{c}^*$, which differs from the direction
$[hkl]=h\mathbf{a}+k\mathbf{b}+l\mathbf{c}$ in any non-orthogonal cell (by 30° for $h=(1,0,0)$ in a
hexagonal cell, 20° for $(0\,0\,1)$ in a monoclinic cell with $\beta=110°$).

**The UI says so.** Before 1.0 the input was headed *Direction*, its boxes were labelled `a`, `b`,
`c`, and the value was printed in square brackets — `[1 1 0]` on the canvas and in the exported file
names — which is the notation for the real-space direction and invited slicing perpendicular to a
different vector. It is now the *Plane (h k l)* control (boxes `h`, `k`, `l`, with a tooltip on the
reciprocal-vector normal) and is printed as `(1 1 0)`: `millerPlaneLabel()` /
`millerPlaneFileLabel()` in `workers/slabSelection.js`, pinned by
`workers/__tests__/millerPlane.test.js`. A true $[uvw]$ mode (converting a real-space direction to its
plane normal through the metric, $\mathbf{h}\propto\mathbf{G}[uvw]$) is not offered.

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

The three *Plane (h k l)* inputs are `type="number"` and `updateCustomDirection()` coerces with
`Number(value)`, so **clearing a box yields 0** (`Number('') === 0`). Clearing all three gives
$\mathbf{h}=(0,0,0)$, and the two runtimes handle a zero vector completely differently:

* **JavaScript** — `normalize(vector, fallback)` returns the fallback when
  $\lVert\mathbf{v}\rVert \le 10^{-9}$. `makeSliceConfig()` calls
  `normalize(customDirection, [0, 0, 1])`, so a zero direction **silently becomes a c-slice**, with
  the label still reading `(0 0 0)`. `makeFreePlaneBasis()` carries two further fallbacks
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

| Case | $m$ | major $\sigma$ (cell-filling slab) | truncation | $e^{-(m/\sigma)^2/2}$ |
| --- | --- | --- | --- | --- |
| defaults, $f=0.03$ | 0.1 | $0.0098$ | $10.2\sigma$ | $\approx 2\times10^{-23}$ |
| slider max, $f=0.15$ | 0.3 | $0.0488$ | $6.1\sigma$ | $\approx 6\times10^{-9}$ |

(Measured on the all-element GaNb₄Se₈ `c` slab at $z_c=0.39$, $\Delta z=0.08$, whose source-atom
covariance has $\sqrt{\lambda}=0.23/0.33$; a uniform cell-filling slab has $\sqrt\lambda=0.289$.
Since 1.0 the kernel comes from the source atoms only (Step 6), so $\sigma$ no longer grows with
$m$. Before, the images entered $\mathbf{C}$ and the slider-maximum row read $\sigma=0.069$ at
$m=0.3$, i.e. $4.3\sigma$ and $8.7\times10^{-5}$ — about $10^3\times$ the $9\times10^{-8}$ this
table then claimed, which had reused the $m=0.1$ spread.) The last column is the Gaussian factor at
the margin; the one-sided mass lost beyond it is smaller still.

The browser's `exponent > -60` guard (Step 7) is a *second*, independent truncation at a Mahalanobis
radius of $\sqrt{120}\approx 11\sigma$. Which of the two binds depends on $f$: with $m=0.1$ the
crossover is near $f\approx0.028$ for this slab, so at the default and above **the image margin is
the binding approximation**, and below it the exponent guard is. The image truncation applies to the
SciPy reference path as well — SciPy has no exponent cutoff, but it never sees the discarded images.

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

**Math.** Each atom's depth is normalized to the slider's units, $\tilde d_i = (\mathbf{x}_i\cdot\hat{\mathbf{h}} - d_{\min})/\Delta_d$,
and the atom is in the slab iff

$$\big|\tilde d_i - z_c\big| \ \le\ \frac{\Delta z}{2} + \epsilon_\mathrm{face}, \qquad \epsilon_\mathrm{face} = 10^{-9}.$$

**One expression in all three places**, evaluated in the same order: `oriented_kde_slice()` builds
`normalized_depth = (_dot3(positions, normal) - depth_min) / depth_span` and `kde_slice()` tests
`np.abs(z - z_center) <= 0.5*max(dz, 1e-12) + SLAB_FACE_TOLERANCE`; the worker's `makeSlab()` and the
page's `inActiveSlab()` both call `isInSlab()` from `workers/slabSelection.js`
(`Math.abs(normalizedDepth - zCenter) <= thickness/2 + SLAB_FACE_TOLERANCE`). `_dot3()` sums
$x h_1 + y h_2 + z h_3$ left to right, like the worker's `dot()`, so presets give bit-identical
depths. Python additionally clamps $z_c$ into $[0,1]$ and floors $\Delta z$ at $10^{-12}$; the
JavaScript path relies on the slider bounds.

**Why the tolerance.** Before 1.0 Python tested `center_depth - half <= d <= center_depth + half` in
absolute depth units and the worker and the Slab-In-Cell highlight tested
`|normalizedDepth − zCenter| ≤ thickness/2`: algebraically identical, but rounded differently exactly
at the faces. An ideal or unrelaxed configuration puts every copy of a site on the same coordinate,
and slider positions hit those faces routinely — `|0.125 − 0.165|` evaluates to
`0.04000000000000001 > 0.08/2`, so for atom layers at $z=k/8$ (as in an ideal cubic cell) with
$z_c = 0.165$, $\Delta z = 0.08$ the browser dropped the whole $z=0.125$ layer that SciPy kept (24 of
4004 slider settings on such a lattice, e.g. also $z_c = 0.21, 0.54, 0.71$ at $\Delta z = 0.08$), and
in Flask mode the highlighted atoms disagreed with the density. $10^{-9}$ is far above the round-off of
a depth in $[0,1]$ ($\sim10^{-16}$) and far below any physical offset.
`tests/test_kde_slab_faces.py` and `workers/__tests__/slabFaces.test.js` check the slab counts of that
lattice against exact integer arithmetic at every $z_c$ on a 0.005 grid for four thicknesses.

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

**$z_c$ and $\Delta z$ are fractions of the projection range of the unit cube along the chosen
normal, and nothing in the KDE pipeline is converted to Å.** They equal cell-edge fractions only for
the `a`/`b`/`c` presets. (`AGENTS.md` used to say *"`z` / `dz` are cell-edge fractions at the
API/slider boundary, converted to Ångström inside `kde.py`"*; the first half holds only for the
presets and the second half was never true of this code — checked by running `oriented_kde_slice()`
with $\mathbf{h}=(1,1,1)$: `depthThickness = 0.1386 = 0.08·√3`, extents in dimensionless
fractional projections, no Å anywhere.)

* `/api/kde/slice` passes `positions.fractional_positions` — *not* `positions.positions` (the Å
  array) — into `oriented_kde_slice()`, with an inline comment saying "Keep the KDE slice in
  fractional coordinates so non-orthogonal cells can be projected through the actual cell basis in
  the frontend."
* Every subsequent quantity ($t$, the slab half-width, the KDE covariance, the evaluation grid, the
  contour coordinates) is therefore dimensionless.
* Ångströms enter **only at draw time**, when `StructurePage.jsx` maps
  $\hat{\mathbf{u}},\hat{\mathbf{v}}$ through `unitCell.unitVectors` to get the real parallelogram
  (`makeProjectedPlane(vectorFromFraction(uVector, unitCell.unitVectors), …)`).

So the contract is: **$z_c$ and $\Delta z$ are fractions of the projection range of the unit cube
along the chosen normal**, and they only coincide with "fraction of a cell edge" for the
`a`/`b`/`c` presets (where $\Delta_d=1$). Both payloads echo them unchanged as `z`/`dz` (and
`center`/`thickness`); the Flask payload also gives `depth`/`depthThickness` in absolute
depth-projection units (see *The two JSON payloads are not the same shape* below).

**Converting $\Delta z$ to Ångström.** Combining the code's definitions with the standard
interplanar spacing $d_{hkl}$ of the plane family $(hkl)=\mathbf{h}$, a depth increment
$\Delta d$ corresponds to a real-space distance $\Delta d \cdot \lVert\mathbf{h}\rVert_2 \cdot d_{hkl}$,
hence

$$\text{slab thickness (Å)} \;=\; \Delta z \; \Delta_d \; \lVert\mathbf{h}\rVert_2 \; d_{hkl}
\;=\; \Delta z \; \big(|h|+|k|+|l|\big)\, d_{hkl}.$$

Sanity checks: for the `c` preset in an orthogonal cell this is $\Delta z\cdot c$ (the naive
reading); for $\mathbf{h}=(1,1,0)$ in a cubic cell it is $\Delta z\cdot 2 \cdot a/\sqrt2 = \sqrt2\,a\,\Delta z$,
which is $\Delta z$ times the full $[110]$ diagonal of the cell — as it must be, since $z_c$ sweeps
$0\to1$ across the whole cube. For $(111)$ in the 10.395 Å GaNb₄Se₈ cell, $\Delta z = 0.08$ is 1.44 Å,
not the 0.83 Å a "cell-edge fraction" reading suggests (1.73×). The page prints this value next to
`d` on the map, `d=0.080 (1.44 Å)`: `slabThicknessAngstrom()` in `workers/slabSelection.js` computes
$\Delta z\,\Delta_d/\lVert\mathbf{A}^{-1}\hat{\mathbf{h}}\rVert$ — the depth gradient
$\mathbf{A}^{-1}\hat{\mathbf{h}} = \hat h_1\mathbf{a}^*+\hat h_2\mathbf{b}^*+\hat h_3\mathbf{c}^*$ is the same
formula written with the reciprocal vectors — pinned by `workers/__tests__/slabThickness.test.js`.

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

**Code.** `rmc_toolkits/kde.py` → `SLAB_FACE_TOLERANCE`, `_dot3()`, `oriented_kde_slice()` (depth
normalization, margin, `_plane_section_vertices()`) and `kde_slice()` (mask/`slabCount`);
`workers/slabSelection.js` → `SLAB_FACE_TOLERANCE`, `isInSlab()`; `localKdeWorker.js` →
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
both unbiased, so the two density fields agree in expectation, but they are not bitwise comparable.
They do share the bandwidth matrix: since 1.0 $\mathbf{C}$ comes from all the slab's source atoms,
before and independently of the subsample (Step 6).

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

**Inputs.** The $n$ fit points $\mathbf{p}_i=(u_i,v_i)$, the slab's $N_\mathrm{src}$ source-atom rows
$\mathbf{q}_a$ (for $\mathbf{C}$), and the *Bandwidth* slider
$f \in [0.005, 0.15]$, step 0.005, **default 0.03**.

**The estimator.** Both paths compute the standard multivariate Gaussian KDE with a
**full-covariance, data-adaptive bandwidth matrix** — *not* an isotropic kernel and *not* a
fixed width in Å:

$$\rho(\mathbf{p}) \;=\; \frac{\kappa}{n}\sum_{i=1}^{n}
\frac{1}{2\pi\sqrt{\det \mathbf{H}}}\,
\exp\!\left[-\tfrac{1}{2}\,(\mathbf{p}-\mathbf{p}_i)^{\!\top}\mathbf{H}^{-1}(\mathbf{p}-\mathbf{p}_i)\right],$$

$$\mathbf{H} \;=\; f^{2}\,\mathbf{C}, \qquad
\mathbf{C} \;=\; \frac{1}{N_\mathrm{src}-1}\sum_{a=1}^{N_\mathrm{src}}(\mathbf{q}_a-\bar{\mathbf{q}})(\mathbf{q}_a-\bar{\mathbf{q}})^{\!\top},$$

where the sum over $i$ runs over the $n$ (subsampled) slab rows, periodic images included, and the
covariance runs over the **source atoms**: $\mathbf{q}_a$ is the in-plane position of source atom
$a$'s representative row — among its rows in the slab, the one with the smallest periodic shift
(`_image_rank()` / `imageRank()`: fewest shifted axes, ties broken by the $(o_x,o_y,o_z)$ loop
order, so an atom inside the slab is represented by itself and a depth-wrapped one by its nearest
image). $\kappa$ is the periodic-image correction defined below. For a preset normal this is the
covariance of the folded in-cell $(x,y)$ of the slab's atoms.

This is SciPy's kernel with the covariance taken from a different point set than the one summed.
`gaussian_kde` with a **scalar** `bw_method` sets `kde.factor = bw` and then
`self.covariance = self._data_covariance * self.factor**2` with
`_data_covariance = atleast_2d(cov(self.dataset, rowvar=1, bias=False, aweights=self.weights))`,
which with the default uniform weights $w_i = 1/n$ is exactly the $n-1$ divisor, and evaluates
$\sum_i w_i \mathcal{N}(\mathbf{p};\mathbf{p}_i,\mathbf{H})$ with those same uniform weights (SciPy
1.13.1, `scipy/stats/_kde.py`). `kde.py` → `_FixedCovarianceKDE` subclasses it and overrides only
`_compute_covariance()` to install $\mathbf{C}$ (`np.cov` of the source-atom rows,
`_source_atom_rows()`) instead of the dataset's covariance; the density is still SciPy's compiled
sum. **No SciPy floor is declared, and the evaluator reads different attributes in different
releases**, so `_compute_covariance()` sets every one of them from $\mathbf{C}$: `covariance` and
`log_det` (all releases); `cho_cov`, the lower Cholesky factor of $\mathbf{H}$ that SciPy ≥ 1.10
whitens with; `inv_cov` $=\mathbf{H}^{-1}$, which earlier releases (1.8.1 checked) whiten with
through its own Cholesky factor (from 1.10 on SciPy defines `inv_cov` as a property that
re-estimates the covariance from the dataset, so the subclass overrides that property too); and
`_norm_factor` $=\sqrt{\det 2\pi\mathbf{H}}$ for the oldest, pure-Python `evaluate`. Because that still leans on SciPy internals,
the constructor evaluates one point against the direct formula below and raises
`ScipyKdeUnsupported` (a `RuntimeError`) when SciPy evaluates anything else — an attribute it cannot
find, or a different value. The tolerance is $10^{-6} + 10\,\kappa(\mathbf{H})\,\varepsilon$:
evaluators that whiten through $\mathrm{chol}(\mathbf{H})$ and through
$\mathrm{chol}(\mathbf{H}^{-1})$ legitimately differ by $O(\kappa\varepsilon)$ (measured
$\le 0.33\,\kappa\varepsilon$, i.e. $1.7\times10^{-6}$ for a needle near
`COVARIANCE_CONDITION_LIMIT`), while one that ignored the supplied covariance would be off by the
change in $\mathbf{C}$ itself. `kde_slice()` turns that error into a declined slab with the
Python-only `engine` message (logged as a warning), so `/api/kde/slice` answers 200 with an
explanation rather than 500. Checked against SciPy 1.8.1 (numpy 1.22), 1.13.1, 1.17.1 and 1.18.1:
the kernel, the peak and the decline decisions agree, the peak to $\le 2\times10^{-9}$ relative
(Step 6 parity fixture).
SciPy ≥ 1.10 evaluates the sum through the lower Cholesky factor
$\mathbf{L}=f\,\mathrm{chol}(\mathbf{C})$ of $\mathbf{H}$ (`cho_cov`): whitened offsets
$\mathbf{w}=\mathbf{L}^{-1}(\mathbf{p}-\mathbf{p}_i)$, kernel $e^{-|\mathbf{w}|^2/2}$, normalization
$1/(2\pi\,L_{00}L_{11})$. The JavaScript `makeKernel()` reproduces it term for term: it forms
$\mathbf{C}$ from the same source-atom rows (`makeSlab()` returns them as `atoms`) with the same
$n-1$ divisor, takes its $2\times2$ Cholesky factor, scales it by $f$, and hands the loop the whitening
matrix $\mathbf{W}=\mathbf{L}^{-1}$ (`w00`, `w10`, `w11`) and
`normalizer = (imageFactor / samples.length) / (2π·L00·L11)`. **There is no ridge and no fallback kernel:
the browser draws exactly $f^2\mathbf{C}$ or declines the slab** (see *Degenerate-slab handling*
below). `tests/generate_kde_fixture.py` writes Python goldens on the committed demo run, the
GaNb₄Se₈ sample run and synthetic slabs, and `workers/__tests__/kdeParity.test.js` requires the
worker to reproduce them to $10^{-6}$ of the peak (measured: $\le 2\times10^{-12}$), with identical
kernels, counts and decline messages.

**$\mathbf{C}$ depends on neither the periodic images nor the subsample** (since 1.0). The images
and the 6000-point cap are evaluation devices — they decide which rows the fixed kernel is summed
over — and the kernel is fitted to the atoms. Before 1.0 both runtimes fitted $\mathbf{C}$ to the
subsampled slab rows, images included, and that made the kernel depend on things that are not the
slab's atoms:

* **the margin.** $m=\min(0.5,\max(0.1,2f,\Delta z))$ decides which images exist, so moving the
  *Thickness* slider past 0.25 or the *Bandwidth* slider past 0.125 admitted the $x/y$ images of
  sites at $\tfrac14,\tfrac34$ and rewrote $\mathbf{C}$ without adding a single atom. On the
  GaNb₄Se₈ Ga layer at $z_c=0.75$ (the same 2000 atoms for every $\Delta z$, $f=0.03$) the kernel's
  principal $\sigma$ went from $0.0019\times0.110$ Å at $\Delta z\le0.2$ to $0.156\times0.191$ Å at
  $\Delta z\ge0.3$, and the peak fell from 1260 to 259; across the bandwidth step $f=0.12\to0.13$
  ($\Delta z=0.08$) the kernel went from $0.26\times0.45$ to $0.67\times0.81$ Å for an 8 % change in
  $f$. Even a handful of image rows mattered: in the Nb layer at $z_c=0.15$, ten image rows (0.25 %
  of 3988) set the kernel's minor $\sigma$ to 0.0107 Å; the source atoms alone give 0.0022 Å. (Now:
  $0.0019\times0.110$ Å and a peak of 1260–1269 for every $\Delta z$ on the Ga layer, the small
  spread being the subsample's Monte-Carlo noise once the slab exceeds 6000 rows.)
* **the subsample seed**, and hence the runtime (PCG64 vs mulberry32): $\sigma$ differed by a few
  tenths of a percent between Python and the browser above the cap.

`tests/test_kde_bandwidth_source.py` and `workers/__tests__/kdeBandwidthSource.test.js` pin the new
behaviour: the kernel is $f^2$ times the covariance of the folded source atoms, identical across a
thickness or bandwidth margin step and across subsample seeds, and a depth-wrapped atom is
represented by its nearest image; the Python file repeats the thickness check on the real Ga layer.

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

A slab whose atoms fill the cell uniformly has $\mathbf{C}\approx \tfrac{1}{12}\mathbf{I}$
($\sqrt\lambda\approx0.289$); the real all-element GaNb₄Se₈ `c` slab at $z_c=0.39$ measures
$\sqrt\lambda = 0.23$ and $0.33$. So the default $f=0.03$ corresponds to $\sigma \approx 0.007$–$0.010$
in fractional units — roughly **0.07–0.10 Å for a 10.4 Å cell**. Three honest corollaries:

1. The smoothing width **changes when you change the element filter, the slab thickness, or the
   normal**, because all of those change $\mathbf{C}$. The same `bw = 0.03` is a different physical
   width on different slices.
2. The periodic images and the subsample do **not** change $\mathbf{C}$ (see above). Before 1.0 they
   did, by far more than the "few percent" this document used to state.
3. The kernel's **shape** — not only its width — comes from $\mathbf{C}$, i.e. from how the slab's
   sites are laid out; see the next section. (This document used to attribute the anisotropy to
   the fractional frame: "Euclidean-isotropic in fractional coordinates, hence anisotropic in Å for
   any non-cubic cell". That was wrong on both counts. A full-covariance KDE is affine-equivariant —
   mapping the points by $\mathbf{M}$ maps $\mathbf{C}$ and $\mathbf{H}$ to $\mathbf{M}\mathbf{C}\mathbf{M}^\top$ and
   $\mathbf{M}\mathbf{H}\mathbf{M}^\top$ — so computing in fractional or in Cartesian coordinates draws the
   same map, and $\mathbf{H}$ is anisotropic in fractional coordinates as well, cubic cells included.)

#### The kernel's shape follows the slab's site layout (read this before reading blob shapes)

$\mathbf{H}=f^2\mathbf{C}$ with $\mathbf{C}$ the covariance of **all** the slab's atoms. A slab holds
many sites, so $\mathbf{C}$ measures how the sites are *laid out* across the slab, not how wide any
one site is, and every site is convolved with that layout-shaped kernel: the second moments of a
drawn blob are the site's own in-plane covariance **plus $\mathbf{H}$**. This is SciPy's convention
for a scalar bandwidth factor and it is kept as the reference estimator for 1.0; it is an artefact of
the method, not of the atoms, and it is largest exactly where the map is most tempting to read —
element-filtered layers holding one or two sites.

**Measured example** (`data/5K_try1/GaNb4Se8_5K.rmc6f`, cubic $a = 10.395$ Å, element Ga, preset `c`,
$z_c = 0.25$ as auto-centred, $\Delta z = 0.08$, default $f = 0.03$). The layer holds two Ga sites on
the face diagonal, at $(\tfrac14,\tfrac34)$ and $(\tfrac34,\tfrac14)$. $\mathbf{C}$ has eigenvalues
$3.5\times10^{-5}$ and $0.125$, so the kernel is a needle: $\sigma = 0.0018 \times 0.110$ Å, 60 : 1,
along $[1\bar10]$. The Ga cloud at $(\tfrac34,\tfrac14)$ is isotropic in the plane
($\sigma = 0.061/0.062$ Å, 1000 atoms), yet its blob in the map has $\sigma = 0.062 \times 0.126$ Å:
**drawn 2.0 : 1 elongated along $[1\bar10]$** (3.7 : 1 at $f = 0.06$). The Nb layer at $z_c = 0.15$
gives a 54 : 1 kernel along $[110]$ and turns a 1.24 : 1 cloud into a 2.1 : 1 blob. The same
mechanism acts in every geometry, only more mildly when the slab's sites fill the plane:

* all atoms, same cell, $(110)$ slice: kernel aspect 1.3 at $z_c = 0.5$ but 8.1 at $z_c = 0.1$;
* a hexagonal layer (three sites of a triangular net, $a = 5$ Å, isotropic 0.08 Å clouds): kernel
  1.61 : 1, so each 3-fold site is drawn as a 1.16 : 1 ellipse; the same crystal written in the
  orthohexagonal cell gives a 1.38 : 1 kernel along a different direction — the picture depends on the
  cell setting;
* even the origin matters: shifting the Ga coordinates by $(\tfrac14,\tfrac14,0)$ — the same
  structure — splits both sites across the cell faces, the folded in-cell covariance becomes round,
  and the kernel is $0.109 \times 0.109$ Å instead of $0.0018 \times 0.110$ Å.

**When the minor axis collapses below the grid the map is aliased.** The Ga needle's
$\sigma_{\min} = 0.0018$ Å is 1/47 of the $G = 120$ node spacing (0.087 Å); the peak reads
1296 / 1352 / 1372 / 1385 and the grid-summed mass 0.967 / 0.856 / 0.924 / 0.960 at
$G = 80 / 120 / 160 / 220$, where it should be 1. A slab holding a single site (an element-filtered
perovskite $B$ layer, say) has $\mathbf{C}$ = the thermal covariance and $\sigma = f\times$ the thermal
spread, $\approx 0.002$ Å.

**What the page does about it.** It does not remove the artefact. It makes it visible: the map's
overlay prints the kernel's principal $\sigma$ in Å (`kernelSigmaAngstrom()` in
`workers/slabSelection.js` maps $\mathbf{H}$ through the in-plane metric: the nonzero eigenvalues of
$\mathbf{M}\mathbf{H}\mathbf{M}^\top$ are those of $\mathbf{H}\mathbf{G}$, $\mathbf{G}=\mathbf{M}^\top\mathbf{M}$);
both engines attach a `subgrid` warning (`KDE_WARNINGS`) when $\sigma_{\min}$ is below half the larger
grid step (`KERNEL_SUBGRID_RATIO = 0.5`: a Gaussian sampled at spacing $h$ keeps its integral to
~1 % while $\sigma\ge h/2$), and an `unresolved` warning when the kernel misses the grid nodes
altogether (the map is then neither contoured nor painted, Step 9); and the page adds a note when
the kernel is more than 3 : 1 in Å
(`KERNEL_ANISOTROPY_NOTE`) saying that elongation along its long axis is an artefact. Read blob
shapes against the printed kernel, and take displacement shapes from the
[PCA Ellipsoid](pca-ellipsoid.md) page, which fits each site's cloud directly. A kernel that is a
physical length (isotropic in the plane, width in Å through the cell metric) would remove the
artefact; it is a different estimator and not part of 1.0.

**Verification.** `tests/test_kde_kernel_diagnostics.py` and
`workers/__tests__/kernelDiagnostics.test.js` pin the `subgrid` warning (a two-site needle warns at
$G = 120$; a cell-filling slab is quiet at $G = 120$ and warns at $G = 16$) and the Å readout; the
parity golden compares the warning codes between the runtimes. The numbers above come from
`oriented_kde_slice()` on the sample run (blob moments as site covariance $+\ \mathbf{H}$).

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
has to be undone — **for the `a`/`b`/`c` presets**, where each source atom has exactly one image in
the drawn cross-section (the unit square). For an oblique normal the section is a triangle to
hexagon whose area changes with $z_c$, and an atom can have zero, one or several images inside it
(the depth window $\Delta z\,\Delta_d$ can exceed the plane period $1/\lVert\mathbf{h}\rVert$), so $\kappa$
no longer makes the section integrate to 1 (below).

**The resulting units.** $\rho$ is a probability density **per unit fractional area of the slice
plane**. For the `a`/`b`/`c` presets it is normalized so that
$\int_{\mathrm{cell}}\rho\,\mathrm{d}u\,\mathrm{d}v \approx 1$: the **cell** carries unit mass, and the
amplitude no longer decays toward the faces. **For an oblique normal it does not**: integrated over
the drawn section polygon, 40 000 uniform random points give 1.004 for $(001)$ at every $z_c$, but
for $(110)$ 1.00 / 0.76 / 0.49 and for $(111)$ 0.79 / 0.60 / 0.24 at $z_c = 0.5 / 0.3 / 0.1$
($\Delta z = 0.08$, $f = 0.03$, $G = 220$) — the level of a *uniform* structure then depends on the
slice position (mean density 0.61 → 3.08 across $z_c$ for $(111)$). On the GaNb₄Se₈ sample the
integral is 1.000 for $(110)$ at $z_c=0.5$ but 1.18 for $(111)$ with $\Delta z = 0.5$ and 0.69 for
$(123)$ with $\Delta z = 0.3$. Only absolute values are affected — the colour scale is per-slice
min–max (Step 10) — but do not compare amplitudes of oblique slices with each other or with presets.
It is **not** one
unit of mass per slab atom (that would integrate to $N_\mathrm{src}$), **not** atoms Å⁻², **not** atoms Å⁻³, and
it is **not divided by the slab thickness** — thickening the slab pulls in more atoms but the field
is renormalized, so absolute values are not comparable between different $\Delta z$, different
elements, or different bandwidths. The per-atom amplitude is $1/N_\mathrm{src}$ of the field.

#### Degenerate-slab handling (identical in both runtimes)

Both runtimes run the same tests in the same order, and a slab that fails one is **declined**:
all-zero grid, `fitCount = 0`, `kernel = null`, and a `message` naming the reason (the strings are
shared verbatim: `KDE_MESSAGES` in `kde.py` and in `localKdeWorker.js`).

| Order | Condition | `message` key |
| --- | --- | --- |
| 1 | no slab rows | `empty` — "No atoms in this slab." |
| 2 | $f$ not a finite number $>0$ | `bandwidth` |
| 3 | slab rows $< 5$ | `too_few` / `tooFew` |
| 4 | source atoms at $<3$ distinct $(u,v)$ positions | `few_unique` / `fewUnique` |
| 5 | centred source-atom positions of rank $<2$ with numpy's tolerance $\sigma_{\min}\le\sigma_{\max}\max(N,2)\,\varepsilon$ | `collinear` |
| 6 | $\mathbf{C}$ not safely positive definite: $c_{00}\le0$, $c_{11}\le0$, $1-\rho^2\le10^{-10}$ (`COVARIANCE_CONDITION_LIMIT`), or a failed Cholesky | `singular` |
| (7) | Python only: the installed SciPy's `gaussian_kde` does not evaluate the supplied kernel (`ScipyKdeUnsupported`, above) | `engine` |

Tests 4–6 run on the source-atom rows that define $\mathbf{C}$ (Step 6 above), before any
subsampling.

Test 5 is `np.linalg.matrix_rank` in Python; the worker's `hasTwoDimensionalSpread()` gets the two
singular values from a twice-orthogonalized Gram–Schmidt QR of the centred columns and the
closed-form SVD of the $2\times2$ triangle, which is as accurate as numpy's SVD, and applies the
same tolerance. Test 6 exists because a slab can pass the rank test and still be collinear to within
round-off — coordinates written with a finite number of decimals put the perpendicular spread at
$\sim10^{-10}$ — and there the sign of the Cholesky pivot $c_{11}-c_{01}^2/c_{00}$ depends on the
summation order, so the two runtimes could disagree on whether to draw. $\rho$ is the in-plane
correlation coefficient of the source atoms; $1-\rho^2\le10^{-10}$ corresponds to a kernel aspect ratio
above $\sim2\times10^5$, i.e. a needle far below any grid spacing, and the limit sits about 100× above
the worst-case summation round-off for 6000 points. The worker's `cholesky2()` applies test 6 and
then fails exactly where LAPACK's `potrf` raises; Python checks test 6 on `np.cov(atoms)` and then
lets `scipy.linalg.cholesky` (inside `_FixedCovarianceKDE`) raise `LinAlgError`.

**Before 1.0 the two paths differed here**, and not only on degenerate input: the browser added a
fixed $10^{-8}$ (fractional²) ridge to the diagonal of $\mathbf{H}$ on every slab, and whenever
$\det\mathbf{H}\le10^{-12}$ it inflated both diagonals by $\max(c_{00},c_{11},10^{-4})f^2+10^{-6}$ and
zeroed the cross term. Both constants were absolute while $\mathbf{H}$ scales as $f^2$ times the
slab spread, so they fired on ordinary, full-rank element-filtered slabs: on the GaNb₄Se₈ Ga layer
($z_c=0.25$, two sites on the cell diagonal, minor eigenvalue of $f^2\mathbf{C}$ ≈ $3\times10^{-8}$ at
$f=0.03$) the ridge alone widened the minor kernel axis 15 %, and at $f\le0.015$ the inflation
branch replaced SciPy's needle by a round kernel $60\times$ wider — maps differing by 36–55 % of the
peak from the SciPy reference, and changing shape discontinuously between two slider steps. Both
constants are gone.

**Bandwidth input.** Both runtimes use $f$ exactly as given. A value that is not a finite number
$>0$ (`0`, negative, `NaN`, `±∞`, or a non-number) declines the slab with the `bandwidth` message
and echoes `bw: null`; it contributes $0$ to the periodic margin. (Before 1.0 the worker silently
substituted `0.03` for any falsy value and floored the rest at $10^{-4}$, while Python drew a
negative $f$ as $|f|$ and returned an unexplained zero grid for `0`.) Only the grid still has a
browser-side substitution: `Number(gridSize) || 120`, so `gridSize = 0` becomes 120, not the 16
lower clamp.

**Verification.** `tests/test_kde.py::test_kde_slice_handles_degenerate_slab_without_error` (five
collinear points → `slabCount = 5`, `fitCount = 0`, `vmin = vmax = 0`) and `tests/test_kde_decline.py`
(every row of the table, plus exact agreement with `scipy.stats.gaussian_kde` on a near-collinear
slab that must be drawn) pin Python; `workers/__tests__/localKdeKernel.test.js` pins the worker's
kernel to an in-test brute force ($<10^{-9}$ of the peak, including needle and single-site kernels)
and its decline rules; the synthetic cases of the parity fixture (`single-site`, `two-site-needle`,
`collinear`, `collinear-to-round-off`, `near-collinear`, `two-positions`, `three-atoms`,
`zero-bandwidth`, `unresolved-needle`) pin the two runtimes to each other.

**Code.** `rmc_toolkits/kde.py` → `KDE_MESSAGES`, `COVARIANCE_CONDITION_LIMIT`, `_valid_bandwidth()`,
`_well_conditioned()`, `_kernel_summary()`, `_source_atom_rows()`, `_FixedCovarianceKDE` (with
`ScipyKdeUnsupported`), `kde_slice()` (the decline chain, the
`_FixedCovarianceKDE(slab.T, covariance, bw)` call and the `density *= slab_total / slab_count`
rescale);
`localKdeWorker.js` → `KDE_MESSAGES`, `COVARIANCE_CONDITION_LIMIT`, `hasDistinctPoints()`,
`hasTwoDimensionalSpread()`, `covariance()`, `cholesky2()`, `kernelSummary()`, `makeKernel()`, and the
decline chain in `computeKde()`.

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

**CPU path (browser).** `computeDensityCpu()` is the direct $O(G^2 n)$ triple loop, in SciPy's
whitened (Cholesky) form:

```js
const w0 = w00*dx;
const w1 = w10*dx + w11*dy;
const exponent = -0.5*(w0*w0 + w1*w1);
if (exponent > -60) sum += Math.exp(exponent);
density[y][x] = sum * normalizer;
```

The whitened form keeps the round-off of a needle kernel inside the squares; the quadratic form
$\mathbf{d}^\top\mathbf{H}^{-1}\mathbf{d}$ used before 1.0 summed terms up to
$\mathrm{cond}(\mathbf{H})$ times larger than the result and cancelled them, which matters on the
float32 GPU branch.

Note the **exponent cutoff at $-60$**: terms with $e^{-60}\approx 8.8\times10^{-27}$ are dropped,
i.e. the kernel is truncated at a Mahalanobis radius $\sqrt{120}\approx11\sigma$. SciPy applies no
such cutoff. The truncation is far below float64 resolution relative to the peak, and it is usually
looser than the periodic-image truncation of Step 3, so it is not a practical difference — but it is
a real difference in the formula evaluated.

**GPU path (browser).** `gpuKde.js` runs the identical sum as a WGSL compute shader, one invocation
per grid cell:

* Workgroup size `(8, 8)`; dispatch $\lceil G/8\rceil$ workgroups in each dimension with an
  in-shader bounds guard `if (gid.x >= grid || gid.y >= grid) { return; }`.
* Bindings: a 48-byte uniform buffer packed by `packKdeParams()` as three `vec4` lanes — `(w00, w10,
  w11, normalizer)` (exactly the fields the CPU loop reads), `(xMin, yMin, xStep, yStep)`, `(grid,
  sampleCount, pad, pad)` read through a `Uint32Array` view of the same `ArrayBuffer`; a read-only
  storage buffer of tightly packed `vec2<f32>` samples (`packKdeSamples()`); a read-write storage
  buffer of $G^2$ `f32` outputs, copied to a `MAP_READ` buffer and read back with `mapAsync`.
* The shader body is line-for-line the CPU expression, including the `e > -60.0` guard.

**When the GPU is used.** Only when the work is big enough to amortize device setup, buffer
uploads, and the readback round-trip:

$$\texttt{shouldUseGpu} \iff G \times G \times n \;\ge\; \mathbf{2{,}000{,}000} \quad (\texttt{GPU\_MIN\_WORK}).$$

At the default $G=120$ that threshold is crossed at $n \ge 139$ fit points, so a normal slab
(thousands of atoms) uses the GPU when WebGPU is present.

**Fallback guarantee, stated precisely.** The device/pipeline promise is created once per worker
and cached; *any* failure — no `navigator.gpu`, `requestAdapter()` returning null, `requestDevice()`
or pipeline creation rejecting, a runtime error inside `computeDensityGpu`, a read-back map with a
non-finite node (since 1.0; pinned by `workers/__tests__/gpuNonFiniteFallback.test.js`), or
sub-threshold work — resolves to `null` and `computeKde()` falls through to `computeDensityCpu()`
(`density = mapped ?? computeDensityCpu(args)`). A lost device clears the cached promise so a later message can
re-initialize. The result object reports which one ran via `backend: 'gpu' | 'cpu'`.

The repo's own wording — `AGENTS.md` and `gpuKde.js`: the CPU loop *"evaluates the same kernel in
float64"* (before 1.0 `AGENTS.md` said *"with identical output"*) — is true **structurally** (same
formula, same grid, same normalizer, same cutoff, reshaped to the same nested JS array) but not
**bitwise**, and the float32 narrowing is broader than the accumulator alone. Everything crosses the
boundary as `f32`:

* the sample coordinates (`new Float32Array(sampleCount * 2)`);
* all four kernel parameters and the grid geometry (`paramFloats[0..7]` = `w00, w10, w11,
  normalizer, xMin, yMin, xStep, yStep`);
* the **node positions themselves**, which the shader reconstructs as
  `xMin + f32(gid.x) * xStep` rather than receiving them from the JS loop, so the grid coordinates
  differ in the last bits too;
* the accumulator (`var sum : f32`), and WGSL `exp` need not match `Math.exp` to the last bit.

For a sum of up to 6000 positive float32 terms the expected relative difference is of order
$10^{-6}$–$10^{-5}$ — invisible in an 8-bit colormap, but it should not be described as identical
arithmetic. For a needle kernel the float32 rounding of the node and atom *positions* ($\sim6\times10^{-8}$)
becomes a visible fraction of the minor kernel $\sigma$ and dominates.
`workers/__tests__/gpuKdeEmulation.test.js` replays `KDE_WGSL` in float32 (`Math.fround` after every
operation) on the buffers `packKdeParams()`/`packKdeSamples()` produce and compares with the CPU loop:
$2.8\times10^{-6}$ of the peak on a cell-filling slab and $2.0\times10^{-4}$ on a two-site needle
($\sigma_{\min}\approx6\times10^{-5}$). A real GPU may fuse multiply-adds and its `exp` is not
correctly rounded, so this bounds the formula, not a particular device — verify on a WebGPU browser.

**Outputs.** `density[G][G]` (nested plain JS numbers / Python list of lists), plus the CPU/GPU
backend flag in static mode.

**Code.** `localKdeWorker.js` → `computeDensityCpu()`, `computeKde()`; `gpuKde.js` → `KDE_WGSL`,
`packKdeParams()`, `packKdeSamples()`, `GPU_MIN_WORK`, `shouldUseGpu()`, `getGpu()`,
`computeDensityGpu()`; `rmc_toolkits/kde.py` → `kde_slice()`.

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

**Whether to contour is decided on the linear density, in both runtimes.** Both engines sum the
*linear* map before the log transform — `kde_slice()` as `grid_mass = sum(density)·Δu·Δv`, the worker
as `linearSum` in the same pass that applies the log — and contour only a drawn map whose mass
reaches `UNRESOLVED_MASS_LIMIT` $=10^{-6}$ (`has_density` / `resolved`). A declined slab
therefore has no contours in either mode, a resolved map has its levels in either mode, and an
`unresolved` map (below) has none in either mode.

**An unresolved map is flagged, not drawn.** The mass of a resolved map is about 1 (the density is
normalised per source atom; 0.24–1.18 on oblique sections, Step 6), and every drawn case of the
parity fixture has at least 0.13, the most aliased two-site needle at $G=32$. A needle kernel can
instead miss every grid node: on the GaNb₄Se₈ `5KAVERAGE` Nb layer ($z_c=0.15$, $\Delta z=0.08$,
$f=0.03$, $G=120$) $\sigma_{\min}=1.4\times10^{-6}$, the linear peak is $1.2\times10^{-14}$ and the
mass $1.8\times10^{-18}$, where a real map peaks at $10^2$–$10^3$. Before 1.0 both runtimes stretched
the per-slice colour scale over that round-off and drew eight contour levels through it. Now both
engines append the `unresolved` warning (after `subgrid`, which such a kernel always has too) when the
mass is below the limit, draw no contours, and `drawKdeSlice()` skips painting a map that carries
it; the payload still holds the density, unchanged. A map that underflows to all zeros is flagged the
same way. Pinned by `tests/test_kde_unresolved.py` and `workers/__tests__/unresolvedMap.test.js`
(a synthetic two-site needle between the node rows, fully underflowed and resolved variants, and the
`5KAVERAGE` layer when `data/` is present) and by the `unresolved-needle` case of the parity fixture.

**A kernel below `KERNEL_MIN_SIGMA` $=10^{-10}$ is not evaluated.** Both engines compute the kernel
summary first and, when $\sigma_{\min}<10^{-10}$ (in-plane fractional units), skip the Gaussian sum
and return the all-zero map with that summary, `fitCount` set, `message = null` and the warnings
`subgrid` + `unresolved` (grid mass 0). SciPy whitens the *absolute* coordinates, $x/\sigma$, so a
node's residual carries a round-off of $\sim\varepsilon|x|/\sigma$ whitened units: harmless at
$10^{-10}$ (the value at an atom is off by $<10^{-7}$ for $|x|\le100$), of order 1 near
$\sigma\sim10^{-13}|x|$, where SciPy returns 0 at an atom. Before 1.0 that made
`_FixedCovarianceKDE`'s self-check fail, and the slice declined with the SciPy `engine` message for
what was a bandwidth $f\lesssim10^{-14}$; below $f\sim10^{-150}$ the normaliser
$1/(2\pi\det L)$ overflows too, and the worker's map came out NaN at every node ($\infty\cdot0$).
A kernel this narrow sits at least $10^{7}$ times below any grid step, so a node farther than
$\sim40\sigma$ from every atom is an exact float64 zero anyway. `/api/kde/slice` therefore answers
200 with that flagged zero map (e.g. `bw=1e-200`), never a NaN map and never a silent one; a result
that still came out non-finite would be a 400 (`_strict_result_response`). Pinned by
`tests/test_kde_unresolved.py`, `unresolvedMap.test.js` and `tests/test_backend_validation.py`.

Before 1.0 the Python guard was `density.max() <= 0` applied **after** the log transform, so any map
whose peak linear density was $\le 1$ per unit fractional area lost every contour on the SciPy path
while the browser drew all eight. That is not rare: the cell carries unit mass, so the mean density
over an oblique cross-section of area $A>1$ is about $1/A$, and for a disordered configuration the
smoothed peak drops below 1 once $f\gtrsim0.1$. Measured on 8000 quasi-uniform points
(`tests/test_kde_contours.py`, (111) slice, $z_c=0.5$, $\Delta z=0.08$, $f=0.1$): peak
$\log_{10}\rho = -0.05$, and 0 contours on the old Python path against 8 now (and 8 in linear mode).
The worker's twin is `workers/__tests__/logContours.test.js`.

**Rendering.** `drawKdeSlice()` strokes every polyline in data coordinates through the same
plane mapper used for the cell outline, at `lineWidth = 1`, in a theme-dependent colour
(`rgba(21,34,50,0.72)` in light theme, `rgba(230,236,244,0.76)` in dark). The *Contours* switch is a
pure client-side toggle — turning it off does not recompute anything.

**Code.** `rmc_toolkits/kde.py` → `_contour_segments()`, `UNRESOLVED_MASS_LIMIT`,
`_kernel_warnings()` and the `has_density` guard in `kde_slice()`; `localKdeWorker.js` →
`extractContours()`, `UNRESOLVED_MASS_LIMIT` and the `resolved` guard in `computeKde()`;
`StructurePage.jsx` → `drawKdeSlice()` contour loop and its `unresolved` paint gate.

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
  the text overlay giving `slabCount`, `fitCount`, $z_c$, $\Delta z$, `bw` and the kernel's $\sigma$ in Å.

#### An empty canvas says why

The draw gate is `density && grid > 0 && kde.vmax > kde.vmin` and no `unresolved` warning (Step 9).
When it fails the canvas prints
`"Computing KDE..."` while a request is in flight, `"No atoms in this slab"` when `slabCount = 0`, and
`"No density drawn for this slab"` when the slab has atoms but the estimator declined it (Step 6); in
that case the payload's `message` — the same string from either runtime — is shown under the canvas
(`kde-message-note`). The overlay still prints `"{slabCount} atoms in slab (fit 0)"`, and
`fitCount = 0` with a non-zero `slabCount` always means "declined", never "empty". An `unresolved`
map (Step 9) also prints `"No density drawn for this slab"`, with a non-zero fit count and the
`unresolved` warning under the canvas.

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
* **Backend cache.** `_cached_positions()` goes through `_POSITIONS_CACHE`, a `_FileCache(16)`
  keyed on `(path, file signature, element)`, so re-slicing the same file does not re-parse it, and
  a parse of a file that changed during the read is never cached. The KDE itself is
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
| `density`, `extent`, `grid`, `bw`, `log`, `slabCount`, `fitCount`, `vmin`, `vmax`, `contours` | ✓ | ✓ | same meaning (`bw` is `null` when the bandwidth was rejected) |
| `kernel` | ✓ | ✓ | $\mathbf{H}$ as `covariance` (in-plane fractional²) plus its principal `sigmaMinor`/`sigmaMajor`; `null` when declined |
| `message` | ✓ | ✓ | why no density was drawn (Step 6), the same string in both; `null` when drawn |
| `warnings` | ✓ | ✓ | `[{code, message}]` about a drawn map, in this order — `subgrid` when the kernel is narrower than half a grid step (Step 6), `unresolved` when the grid holds less than $10^{-6}$ of the density (Step 9); `[]` otherwise |
| `center`, `thickness` | ✓ | ✓ | the raw slider fractions in **both** — this is what the UI reads |
| `normal`, `uVector`, `vVector`, `planeVertices`, `planePolygon` | ✓ | ✓ | `uVector`/`vVector` differ for custom normals (Step 2) |
| `z`, `dz` | ✓ | ✓ | the slider fractions in both (since 1.0; see below) |
| `depth`, `depthThickness`, `depthRange` | ✓ | — | absolute depth-projection units |
| `slabVertices` | ✓ | — | the two clamped slab faces; **not consumed by `StructurePage.jsx`** |
| `cellLengths`, `unitVectors`, `orientation`, `source`, `element` | ✓ (added by the endpoint) | — | `cellLengths`/`unitVectors` are ignored by the frontend (Step 10) |
| `browserKde: true`, `backend: 'gpu' \| 'cpu'` | — | ✓ | worker-only |

**`z`/`dz` are the slider fractions in both runtimes.** Before 1.0 Flask returned `z = center_depth
= ` $d_{\min} + z_c\Delta_d$ and `dz = thickness_depth = ` $\Delta z\,\Delta_d$ (so for
$\mathbf{h}=(1,-1,0)$ and $z_c=0.5$ it returned `z = 0`), while the worker returned the fractions.
Now `kde_slice()` receives normalized depths and echoes $z_c$ and $\Delta z$ (after Python's clamp
and floor); the absolute depth-projection values remain available as `depth`, `depthThickness` and
`depthRange`.

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

Added for 1.0:

| Test | What it pins |
| --- | --- |
| `tests/test_kde_decline.py` | every decline rule of Step 6 (with its `message`), the bandwidth validation, exact agreement with `scipy.stats.gaussian_kde` on a near-collinear slab, and the `kernel` summary |
| `tests/test_kde_slab_faces.py` / `workers/__tests__/slabFaces.test.js` | face atoms of an ideal $k/8$ lattice are in the slab at every slider position (exact integer reference); the face tolerance; `z`/`dz` echo the slider fractions |
| `tests/test_kde_parity_fixture.py` | the committed browser-parity golden is still what `kde.py` computes (re-run `tests/generate_kde_fixture.py` when it fails): kernel covariance to $10^{-12}$, map shape to $10^{-8}$ of the peak, and the peak to $10^{-10}+20\,\kappa(\mathbf{H})\,\varepsilon$ relative, because SciPy releases round the Gaussian sum differently and a needle kernel amplifies that by its condition number $\kappa$ (near-collinear, $\kappa=6.6\times10^6$: $6\times10^{-11}$ on SciPy 1.17, $1.7\times10^{-9}$ on 1.8, against the 1.13 golden) |
| `workers/__tests__/kdeParity.test.js` | **cross-runtime**: the worker reproduces the Python golden — demo run, GaNb₄Se₈ run (skipped when `data/` is absent), synthetic slabs — to $10^{-6}$ of the peak, with identical `slabCount`, `fitCount`, kernel and `message` |
| `workers/__tests__/localKdeKernel.test.js` | the worker's kernel equals an in-test brute-force $f^2\mathbf{C}$ mixture; its rank test and decline rules |
| `workers/__tests__/gpuKdeEmulation.test.js` | the WGSL shader, replayed in float32 on the packed buffers, against the CPU loop |
| `tests/test_kde_unresolved.py` / `workers/__tests__/unresolvedMap.test.js` | a kernel that misses every grid node is flagged `unresolved` and draws no contours in either scale, a fully underflowed map too, the same layer at a wider bandwidth is not flagged, a declined slab carries no warning, and the `5KAVERAGE` Nb layer needle is flagged (skipped without `data/`) |
| `tests/test_kde_scipy_compat.py` | `_FixedCovarianceKDE` hands the supplied kernel to every SciPy evaluator (the ≥ 1.10 one it runs on, plus in-test replicas of the pre-1.10 `inv_cov` path and the pure-Python `_norm_factor` path); `inv_cov` is $\mathbf{H}^{-1}$ and reading it leaves $\mathbf{C}$ alone; the construction check accepts an $O(\kappa\varepsilon)$ disagreement on a needle and rejects a $10^{-3}$ one; a SciPy that cannot evaluate the kernel declines with `engine` in `kde_slice()` and over `/api/kde/slice` (HTTP 200) |
| `tests/test_kde_bandwidth_source.py` / `workers/__tests__/kdeBandwidthSource.test.js` | $\mathbf{C}$ is the covariance of the folded source atoms, unchanged across a thickness or bandwidth margin step, across subsample seeds, and for depth-wrapped atoms (plus the real Ga layer in Python); `_FixedCovarianceKDE` sums SciPy's Gaussian with the supplied covariance and raises if SciPy stops honouring it |
| `tests/test_kde_contours.py` / `workers/__tests__/logContours.test.js` | log scale keeps all eight contours on an oblique disordered slice whose log peak is negative; a declined slab has none |
| `tests/test_kde_kernel_diagnostics.py` / `workers/__tests__/kernelDiagnostics.test.js` | the `subgrid` warning; the kernel's σ in Å through the cell metric |
| `workers/__tests__/millerPlane.test.js` / `slabThickness.test.js` | the custom slice is labelled and selected as the plane family $(hkl)$; the slab thickness in Å |

**What is still not covered.** The cross-runtime golden uses slabs below the 6000-point fit cap, where
both runtimes sum the same rows; above it they draw different subsamples (Step 5) and agree only
statistically. Contour polylines are not compared (the tracers differ by design, Step 9). WebGPU itself
cannot run under vitest: the emulation pins the formula and the packing, not a device's `exp` or its
fused multiply-adds.

---

### Parameters and defaults

| Parameter | UI control | Default | Range / options | Units | Where enforced |
| --- | --- | --- | --- | --- | --- |
| Element | select | `all` | elements found in the file | — | Flask: `load_unit_cell_positions(element=)` server-side; browser: the `points` `useMemo` client-side |
| Normal $\mathbf{h}$ | menu | `c` = $(0,0,1)$ | `a`, `b`, `c`, `Custom` | Miller indices $(hkl)$, dimensionless | `SLICE_PRESETS` / `SLICE_ORIENTATIONS` |
| Plane (h k l) | 3 number inputs | $(1\,1\,0)$ | any (step 0.1) | Miller indices, dimensionless | `StructurePage.jsx` → `customDirection`; zero vector → silent $(0,0,1)$ |
| Slice centre $z_c$ | range slider | auto-set to densest of 50 depth bins; state default 0.5 | 0 – 1, step 0.001 | fraction of depth span $\Delta_d$ | slider; Python clamps to $[0,1]$ |
| Thickness $\Delta z$ | range slider | **0.08** | 0.01 – 0.5, step 0.01 | fraction of depth span $\Delta_d$ | slider; Python floors at $10^{-12}$ |
| Bandwidth $f$ | range slider | **0.03** | 0.005 – 0.15, step 0.005 | dimensionless covariance factor | slider; both runtimes use it as given and decline a non-finite or non-positive value (`bandwidth` message) |
| Grid $G$ | select | **120** | 80 / 120 / 160 / 220 | nodes per side | clamp 16–400 (Py), 16–260 (JS, after substituting 120 for any falsy value) |
| Contour levels $K$ | none (fixed) | **8** | request param `levels` | count | frontend always sends 8; empty levels are dropped from the output |
| Colormap | select | `viridis` | viridis, magma, seismic, reds, greys | — | `colormaps.js` |
| Contours | switch | on | on/off | — | client-side only |
| Log scale | switch | **on** | on/off | — | `kde.py` / worker |
| Fit-point cap | none | **6000** | fixed | count | `MAX_KDE_FIT_POINTS`, `fitLimit` |
| Subsample seed | none | **0** | fixed | — | `rng_seed=0`; `randomUnit(0)` |
| Periodic margin $m$ | none | $\min(0.5,\max(0.1,2f,\Delta z))$ | derived | fractional | both paths; a rejected $f$ counts as 0 |
| Augmentation factor | none | $(1+2m)^3$ | derived | ratio | 1.73× at $m=0.1$, 8× at $m=0.5$ |
| Log floor | none | $10^{-12}$ | fixed | density | both paths |
| Exponent cutoff | none | $-60$ (browser only) $\Rightarrow 11\sigma$ | fixed | — | `localKdeWorker.js`, `gpuKde.js` |
| Covariance conditioning limit | none | $1-\rho^2\le10^{-10}$ declines | fixed | — | `COVARIANCE_CONDITION_LIMIT` in `kde.py` and `localKdeWorker.js` |
| GPU work threshold | none | $G^2 n \ge 2{,}000{,}000$ | fixed | work units | `GPU_MIN_WORK` |
| Display atom cap | none | 1 000 000 (clamped to $\ge100$ in Flask) | fixed | atoms | `STRUCTURE_MAX_POINTS`, `MAX_STRUCTURE_POINTS` |
| Plane-section tolerances | none | $10^{-9}$ (on-plane corner), $10^{-8}$ (dedup), $\ge3$ vertices | fixed | fractional | `_plane_section_vertices()` / `planeSectionVertices()` |
| Canvas minimum size | none | 320×260 (KDE), 220×260 (slab) | fixed | CSS px | `StructurePage.jsx`; also applied to the 3× export |
| Debounce | none | 160 ms (Flask) / 80 ms (worker) | fixed | ms | `StructurePage.jsx` |

---

### Python vs JavaScript: exact parity table

Derived by code reading, and for the kernel, the counts, the decline rules and the density grid **measured**
by `kdeParity.test.js` against Python goldens (slabs below the fit cap; see the test-suite section).

| Stage | Agreement |
| --- | --- |
| Unit-cell folding | **Exact** (same modulo convention, sign-corrected in JS) |
| Periodic image tiling + margin | **Equivalent for folded inputs** (Python: unconditional originals + 26 shifted offsets; JS: all 27 offsets margin-tested). **Row order differs**, so index-based subsampling picks different points |
| Slab selection in depth | **Exact** (tested): the same normalized-depth expression with the same $10^{-9}$ face tolerance in `kde.py`, the worker and `inActiveSlab` |
| `slabCount` (unique source atoms with an image in the slab) | **Exact** |
| Subsample size (6000) | **Exact**; the **selected subset differs** (PCG64 vs mulberry32, and a different row order) |
| Bandwidth matrix $\mathbf{H}=f^2\mathbf{C}$ | **Exact** (tested, $<10^{-9}$ relative): $\mathbf{C}$ from the same source-atom rows in both, before and independently of the subsample. No ridge, no substitution for $f$ |
| Kernel normalization $1/(2\pi n\sqrt{\det\mathbf{H}})$ | **Exact** (both through the Cholesky factor: $1/(2\pi n L_{00}L_{11})$) |
| Periodic renormalization $\kappa=N_\mathrm{img}/N_\mathrm{src}$ | **Exact** |
| Evaluation grid nodes | **Exact** on the CPU paths; grid clamp maxima differ (400 vs 260); the GPU path recomputes node positions in `f32` |
| Kernel sum | SciPy: full float64, no cutoff. Browser CPU: float64 with $e<-60$ cutoff; **measured $\le2\times10^{-12}$ of the peak** against SciPy on the parity fixture. Browser GPU: **float32 inputs, parameters, node grid and accumulator**, same cutoff (emulated: $\sim3\times10^{-6}$, $2\times10^{-4}$ for a needle kernel) |
| Degenerate slabs | **Identical** (tested): the same six decline tests in the same order, the same `message` strings |
| $\log_{10}$ transform + $10^{-12}$ floor | **Exact** |
| Contour level values | **Exact** (same formula; both drop levels that yield no polylines) |
| Contour tracing | contourpy stitched polylines with saddle handling vs. per-cell 2-point segments with arbitrary saddle pairing |
| Contours in log mode | **Same rule** (tested): contour iff the *linear* density has a positive maximum |
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
   switches the page to the browser worker. Below the 6000-point fit cap the browser CPU path draws
   the same kernel and reproduces the reference to $10^{-6}$ of the peak (tested); above it the two
   draw different subsamples, and the GPU branch is float32.
3. **The bandwidth is not a length.** $f$ multiplies the *covariance of the slab's atoms*, so
   the physical smoothing width changes with the element filter, the slab position and thickness
   (through which atoms are selected), and the slice normal — but no longer with the periodic
   margin or the subsample. `bw = 0.03` is
   roughly $8\times$ narrower than Scott's rule at $n=6000$ — the map is deliberately under-smoothed
   to resolve sites, which means fine structure in the map can be sampling noise rather than real
   density.
4. **The estimate is subsampled.** At most 6000 slab points (images included) are fitted. The
   estimate stays unbiased but its pointwise noise grows as $\sqrt{N_\mathrm{img}/6000}$. The canvas always
   prints both counts.
5. **The density units are fractional, not Å⁻².** With $\kappa = N_\mathrm{img}/N_\mathrm{src}$ applied, the field integrates
   to about **1 over the whole cell for the `a`/`b`/`c` presets** — the cell carries unit mass; it is
   *not* 1 per slab atom. For an oblique normal the integral over the drawn section varies with
   $z_c$ and $\Delta z$ (0.24–1.18 measured, Step 6). It is not divided by the slab thickness, so values
   are not comparable across different $\Delta z$, elements, bandwidths or slices.
6. **The colour scale is per-slice.** No colorbar, no shared normalization, no numeric legend. Two
   screenshots of this panel cannot be compared quantitatively. For an oblique normal the scale is
   additionally set by grid nodes that lie **outside** the drawn cross-section and are clipped from
   the display.
7. **The colormaps are 5-anchor approximations** of the matplotlib maps of the same name and are not
   perceptually uniform.
8. **The kernel's shape follows the slab's site layout** (Step 6), in cubic cells too: an isotropic
   site can be drawn 2 : 1 or more along the line joining the slab's sites, the shape depends on the
   cell setting and origin, and a one- or two-site slab can collapse the kernel below the grid
   (flagged `subgrid`). The overlay prints the kernel's $\sigma$ in Å. **Everything geometric
   happens in fractional space**: the custom plane is a Miller index triple $(hkl)$ (and labelled as
   one), not a real-space vector, and the slider's $\Delta z$ is a fraction of the cube's depth
   range, not of a cell edge. The overlay prints the slab thickness in Å next to it
   (`d=0.080 (1.44 Å)`, $\Delta z(|h|+|k|+|l|)d_{hkl}$, Step 4); API consumers of `z`/`dz` must
   do that conversion themselves.
9. **The slider auto-jumps.** Changing the element filter or the normal — including typing a single
   digit into a custom-direction box — re-runs the 50-bin "densest layer" search and overwrites $z_c$.
   That search runs on the unwrapped, display-sampled population.
10. **The GPU result is float32** in its inputs, its kernel parameters, its reconstructed node grid
    and its accumulator. Structurally identical to the CPU loop, numerically equal to about
    $10^{-6}$–$10^{-5}$ relative ($2\times10^{-4}$ for a needle kernel) in a float32 emulation of the
    shader — fine for a picture, not a bit-for-bit guarantee, and no test runs a real GPU.
11. **Contours need a resolved linear density, nothing more.** Log scale changes the levels (equally
    spaced in $\log_{10}\rho$), never whether contours are drawn; before 1.0 the Flask path dropped
    all of them whenever the peak density was below 1 per unit fractional area. A map whose kernel
    misses every grid node (grid mass $<10^{-6}$, flagged `unresolved`) is neither contoured nor
    painted, in both runtimes (Step 9).
12. **The periodic wrap is exact only out to the margin $m$.** Images farther than $m$ from the cube
    are discarded, which truncates a cell-filling slab's kernel at $\approx10\sigma$ at the defaults
    and $\approx6\sigma$ at $f=0.15$ (Step 3). This applies to the SciPy path too.
13. **`slabCount` is "atoms with at least one image in the slab"**, and one atom can contribute
    several rows near an edge or corner. The drawn band and the highlighted atoms in the side view
    are clamped/unwrapped and so **understate** the selection for $z_c$ near 0 or 1.
14. **A declined slab draws nothing, in both runtimes, and says why** — fewer than 5 rows, fewer than
    3 distinct in-plane points, a collinear spread, a covariance singular to round-off, or an invalid
    bandwidth (Step 6). The canvas reads "No density drawn for this slab" and the reason is printed
    under it; `fitCount` is 0.
15. **A missing lattice block degrades silently.** If `structure.latticeVectors` or
    `structure.supercell` is absent the frontend draws a 1 Å cubic cell with no warning.
16. **The page may not be analysing the file you think.** In a folder with several `.rmc6f`
    configurations, the one paired with a recognised output file wins; with no match both runtimes
    fall back to the first usable one by name (Step 1).
17. **Skipped atom lines are reported, not hidden.** The zero-atom static-mode issue (`AGENTS.md`
    *Current known issues*, 2026-06-18) is resolved in 1.0: both parsers share the validated grammar
    of Step 1, report skipped or unparsed lines against the header's atom count, and a file with no
    parseable atom is an error naming what was found (HTTP 400 from `/api/structure`, `/api/pca/*`
    and `/api/triplets`; a thrown error in the browser) rather than an empty card.
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
| $\rho$ | KDE density sample on the grid (see the KDE section for its normalization) | probability density per in-plane fractional unit² (integrates to ≈1 over the cell for the a/b/c presets only — KDE Step 6) |
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

`app.py` → `_find_rmc6f()`: if the resolved target is itself a `.rmc6f` file, use it. Otherwise
`parsers.find_run_configuration()` (the rule the `rmc-triplets` CLI uses too) globs `*.rmc6f` (error
if none), keeps the usable candidates (not 0-byte, an Atoms marker in the first 64 KiB), then walks
the directory's other files in case-insensitive name order and classifies each by
`run_stem_from_output_name()`, which extracts a run stem at one of three priorities:

| priority | pattern |
| --- | --- |
| 0 | `<stem>-NN.log` (2+ digits) |
| 1 | `<stem>-EXAFS-*_Q_OUTPUT.csv` / `_R_OUTPUT.csv`, `<stem>_FT_XFQ<n>.csv`, `<stem>_FQ<n>.csv`, `<stem>_SQ<n>.csv`, `<stem>_bragg*.csv`, `<stem>_PDF*.csv` |
| 2 | `Frac_coord_<stem>.txt` |

The `(priority, filename)` list is sorted and the first usable `.rmc6f` whose stem matches wins; with
no match, the **alphabetically first** usable `.rmc6f` is used.

`browserData.js` → `chooseStructureFile()` uses the same priority ladder and the same tie-breaking
(candidates by path, outputs by (priority, lower-cased name, stem), code-point order) but matches on
`dirname + '/' + stem`, so an output file only claims an `.rmc6f` sitting in the *same* subfolder.
Before 1.0 its fallback was the first `.rmc6f` in directory-walk order, so a folder holding several
models could resolve to a different model on the two paths.

#### 1b. Metadata

Both parsers scan for two headers and are otherwise position-independent:

- a line whose first token is `Supercell` → $N_i$ from the **last three** tokens
  (`Supercell dimensions:  10 10 10`);
- a line whose first token is `Lattice` → the **next three lines** parsed as the rows
  $\mathbf{L}_1,\mathbf{L}_2,\mathbf{L}_3$ in Å.

Python: `rmc_toolkits/parsers.py` → `read_cell_vectors()`. JavaScript:
[rmc6f.js](../../web_app/frontend/src/rmc6f.js) → `readRmc6fCellVectors()`. The two agree
exactly; both raise/throw if either header is missing, and — since 1.0 — if the supercell is not
three positive integers or the lattice is not a finite, non-singular 3×3 matrix with a finite
volume (numbers read Fortran-aware, so `D` exponents work; checks and messages in
[run-dashboard.md](run-dashboard.md), Step 1). `Supercell dimensions: 0 0 0` used to fold every
atom onto the origin and still draw a map. Note that the `Cell (Ang/deg): a b c α β γ`
line, when present, is **ignored** — the cell geometry always comes from the lattice-vector rows, so
a triclinic cell is handled by construction.

#### 1c. Atom lines

Atom records start after the `Atoms` marker line (any spelling: `Atoms:`, `Atoms :`, `atoms:`,
`Atoms (fractional coordinates):`, case-insensitive). Both runtimes classify each line with **one
anchored grammar** — Python `rmc_toolkits/parsers.py` → `classify_rmc6f_atom_line()` (used by
`iter_rmc6f_atoms()` / `parse_rmc6f_atoms()`), JavaScript
[rmc6f.js](../../web_app/frontend/src/rmc6f.js) → `classifyAtomLine()` (used by `parseRmc6fAtoms()`):

```
id  element  [label]  x y z  ref  n1 n2 n3     full layout: exactly 7 data fields
id  element  [label]  x y z                    legacy coords-only: exactly 3 data fields
```

The label is optional: a bracket group (`[1]`, or split as `[ 1]`) or one non-numeric token.
Numbers accept Fortran `D` exponents; `ref` must be a positive integer and each cell index an
integer in $[0, N_i)$ ($N$ from the `Supercell` line). A coords-only atom comes back with
`reference_number`/`cell_indices` (`referenceNumber`/`cellIndices`) null. A line of valid layout
whose coordinates are non-finite (`NaN`, `Inf`, Fortran `****`) is **skipped and counted**
separately; every other line after the marker that fits no layout is counted as unparsed. Element
tokens are capitalized the same way in both runtimes (Python `str.capitalize()`: `SE`/`se` → `Se`).
The report (`Rmc6fParseReport` / `report`) holds those counts and the header's `Number of atoms:`,
and its warning is shown on the Model information card and returned as `parseWarning` by the
structure, PCA and bond-angle paths. `iter_rmc6f_atoms()` yields full-layout atoms by default;
position-only consumers (`load_unit_cell_positions()`, the triplets loader, `/api/structure`) pass
`include_coords_only=True`. `read_atom_indices()` is built on the same iterator, so the site table
and the atom list cannot disagree. Tests: `tests/test_parsers_rmc6f_grammar.py` and
`web_app/frontend/src/__tests__/rmc6f.test.js`, plus the coords-only golden
`rmc6f_coords_only_fixture.json` (`tests/test_coords_only_fixture.py` ↔
`coordsOnlyParity.test.js`).

> **Before 1.0** the two parsers indexed fields from the *end* of the line. Python required ≥ 9
> tokens, so a coords-only file gave **zero atoms** on the backend path while the browser parsed it;
> an extra trailing field silently shifted every column in the browser; the browser left `GA`/`se`
> uncapitalized, which changed the element filter, legend and colours between the two paths; and
> `read_atom_indices()` keyed on the raw token and read `parts[-4]` of any line, so a file whose
> element column was not title-case showed 0 sites per element on the Model information card.

#### 1d. Folding

Backend path ([app.py](../../web_app/backend/app.py) → `structure()`):

$$x_i = \big(f_i N_i\big) \bmod 1$$

Browser path ([browserData.js](../../web_app/frontend/src/browserData.js) → `structureFromRmc6f()`):

$$x_i = \big(((f_i N_i) \bmod 1) + 1\big) \bmod 1$$

**They are identical**: the extra `+1 %1` in JS only fixes JavaScript's sign-preserving `%` for
negative coordinates (Python's `%` is already non-negative). To floating point they agree exactly.
(Before 1.0 the backend first subtracted $n_i/N_i$, which removes an integer after multiplication by
$N_i$ and so cannot change the value mod 1 — but needed the cell-index columns.)

Both are index-free: `/api/structure` folds a coords-only atom (no cell indices) straight from its
coordinates, and `load_unit_cell_positions()` in `rmc_toolkits/kde.py` uses the
index-free form `(coords * supercell) % 1.0` on every atom `iter_rmc6f_atoms(...,
include_coords_only=True)` yields — so a legacy coords-only file gives the same folded positions in
both runtimes (before 1.0 it gave **no atoms at all** on the backend path).

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
(`element in (None, '', 'all')` means no filter) and is memoized per `(path, file signature, element)`
by `_cached_positions`. Both sides of that comparison come from `iter_rmc6f_atoms()`, so they are
consistently capitalized; but the same control therefore acts through **two different mechanisms at
two different points of the pipeline**, and the symbol it matches is the capitalized one on the
backend path and the file's raw token on the browser path.

The element colour map and legend are *not* filtered — see Step 8.

#### 1h. Related but not used by this page: the `Frac*.txt` conversion

`rmc_toolkits/parsers.py` → `frac_lines_from_rmc6f()` writes the classic `Frac_coord_*.txt` file
(exposed as `POST /api/convert/frac` → `write_frac_from_rmc6f()`): a 5-line header followed by one
row per atom,

$$\mathrm{reduced}_i = f_i - \frac{n_i}{N_i} \quad\text{printed as}\quad \texttt{RN  x  y  z  Nx  Ny  Nz}$$

with the coordinates formatted to **5 decimal places** of a *box* fraction. `read_structure()` pairs
`Frac_coord_<stem>.txt` with the usable `<stem>.rmc6f` of the same configuration (a folder with
exactly one of each pairs them regardless of name; any other ambiguity raises; `frac_path=` /
`rmc6f_path=` choose explicitly), cross-checks the Frac cell indices against the `.rmc6f` supercell
and its reference numbers against the `.rmc6f` sites, then reads the Frac file back, skipping exactly the first 5 lines, and re-expands
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

**Inputs:** the `Normal` control (`a`, `b`, `c`, or *Plane (hkl)* with three numeric fields
$h,k,l$, default $(1\,1\,0)$). **Output:** `sliceConfig = { key, label, normal, u, v, uLabel, vLabel,
range }` (plus `fileLabel` for a custom plane), built by `makeSliceConfig()`.

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
> `customDirection`, still reads `(0 0 0)`. `updateCustomDirection` stores `Number(event.target.value)`
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
\text{in slab} \iff \big|\tilde d - z_c\big| \le \tfrac{\Delta z}{2} + 10^{-9}\quad(\texttt{inActiveSlab} \to \texttt{isInSlab})$$

with $\Delta_d$ replaced by 1 if it evaluates to 0.

The same **depth convention and the same inclusive test** — `isInSlab()` in
[`workers/slabSelection.js`](../../web_app/frontend/src/workers/slabSelection.js), with a
$10^{-9}$ face tolerance added to $\Delta z/2$ — is used by `inActiveSlab`, by the browser KDE worker
([localKdeWorker.js](../../web_app/frontend/src/workers/localKdeWorker.js) → `makeSlab`) and by the
server (`rmc_toolkits/kde.py` → `oriented_kde_slice`/`kde_slice`, on the same normalized depth, after
clamping $z_c$ to $[0,1]$ and $\Delta z$ to $\ge 10^{-12}$), so an atom exactly on a face is in the
slab in all three (KDE Step 4).

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
   `500 13px Inter`, at $(14, 28)$: `Computing KDE...` while a request is in flight,
   `No density drawn for this slab` when the slab has atoms but the estimator declined it (the
   payload's `message` is then shown under the canvas), otherwise `No atoms in this slab`.
6. Overlay text (drawn with a dark stroke `rgba(13, 18, 28, 0.62)`, `lineWidth 3`, under a white fill,
   so it stays legible over any colormap) reports `<slabCount> atoms in slab (fit <fitCount>)` at
   $(12,22)$, `<label>=<center>  d=<thickness> (<Å> Å)  bw=<bw>` at $(12,40)$ (the Å value from
   `slabThicknessAngstrom()`, KDE Step 4), `kernel σ <minor> × <major> Å`
   at $(12,58)$ when a kernel was drawn (KDE Step 6), and `log10 density` on the next line (18 px
   lower) when the log toggle is on. The engine's `warnings` and the page's kernel-anisotropy note are
   shown under the canvas. Recall from Step 3c/5.5 that `slabCount` counts unique source
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

File names are `KDE_Slice_<normal>.png`, `Slab_In_Cell_<normal>.png`, and `Folded_Unit_Cell.png`.
For a preset the name is passed through `sanitizeFilename()`, which (i) unwraps LaTeX-style
superscripts `^{…}` → `…`, (ii) collapses every run of characters outside `[A-Za-z0-9_.-]` to a
single `_`, (iii) strips leading/trailing `_`, and (iv) returns the literal `figure` if nothing
survives: `KDE_Slice_c.png`. For a custom plane `sliceFileName()` appends
`millerPlaneFileLabel()` after sanitizing the prefix, so the Miller parentheses survive and the
indices are **locale-independent** (rounded to 0.01, `.` as the decimal point):
`KDE_Slice_(1_1_0).png`, `Slab_In_Cell_(1.5_1_0).png`. The on-canvas label `(1 1 0)` still uses
`Number(v).toLocaleString(undefined, { maximumFractionDigits: 2 })`, so a component of 1.5 reads
`1,5` where the decimal separator is a comma. (Before 1.0 the custom label was `[1 1 0]` and the
file `KDE_Slice__1_1_0.png`, locale-dependent.)

---

### Parameters and defaults

| Parameter | Code location | Default | Range / values | Units |
| --- | --- | --- | --- | --- |
| `selectedElement` | `StructurePage.jsx` state | `all` | `all` + elements found in the file | — |
| `sliceDirection` | state / `NORMAL_OPTIONS` | `c` | `a`, `b`, `c`, `custom` | — |
| `customDirection` | state | `[1, 1, 0]` | any 3 reals (number inputs, step 0.1); all-zero ⇒ `[0,0,1]` | Miller indices $(h\,k\,l)$ of the sliced planes |
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
| KDE placeholder text | `drawKdeSlice` | `Computing KDE...` / `No density drawn for this slab` / `No atoms in this slab` at (14, 28), `500 13px Inter`, `--muted` | — | CSS px |
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
- **Entering an all-zero custom plane silently gives the c-axis slice**, while the panel label
  still reads `(0 0 0)`.
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
