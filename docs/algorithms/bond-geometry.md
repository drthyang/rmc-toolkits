# Bond Geometry — algorithm reference

> Part of the [Algorithms and Math Reference](../ALGORITHMS.md). Every step is anchored to the source; if this document and the code disagree, **the code wins**.

Bond angles the RMCProfile `triplets` way: name an A–B–C triplet with **B the central atom**,
bound the two bond lengths, and histogram the angle at B over every triplet in the periodic
configuration — exactly, images included, nothing subsampled. The engine section covers the
neighbour search and the three normalizations; the page section covers what the Bond Geometry
tab adds on top of the payload (it computes no geometry of its own).

## Contents

- [Bond Geometry — the triplet engine](#bond-geometry--the-triplet-engine)
  - [What this page shows](#what-this-page-shows)
  - [Step 1 — Inputs, validation, and element matching](#step-1--inputs-validation-and-element-matching)
  - [Step 2 — Folding and the image bookkeeping](#step-2--folding-and-the-image-bookkeeping)
  - [Step 3 — The linked-cell search: cells, reach, and correctness](#step-3--the-linked-cell-search-cells-reach-and-correctness)
  - [Step 4 — Which pairs are bonds](#step-4--which-pairs-are-bonds)
  - [Step 5 — Pairing bonds into angles](#step-5--pairing-bonds-into-angles)
  - [Step 6 — The histogram and its three normalizations](#step-6--the-histogram-and-its-three-normalizations)
  - [Step 7 — The summary payload](#step-7--the-summary-payload)
  - [Step 8 — The two app boundaries and their caps](#step-8--the-two-app-boundaries-and-their-caps)
  - [The `rmc-triplets` CLI](#the-rmc-triplets-cli)
  - [Parameters and defaults](#parameters-and-defaults)
  - [Parity: Python engine vs JavaScript port](#parity-python-engine-vs-javascript-port)
  - [Caveats](#caveats)
- [Bond Geometry — the page](#bond-geometry--the-page)
  - [What the page owns](#what-the-page-owns)
  - [Step 1 — Triplet seeding](#step-1--triplet-seeding)
  - [Step 2 — The compute request and the epoch guard](#step-2--the-compute-request-and-the-epoch-guard)
  - [Step 3 — The result chips](#step-3--the-result-chips)
  - [Step 4 — The angle plot and the `fit` variant](#step-4--the-angle-plot-and-the-fit-variant)
  - [Step 5 — The partial-g(r) window helper](#step-5--the-partial-gr-window-helper)
  - [Step 6 — The folded-cell bond view](#step-6--the-folded-cell-bond-view)
  - [Parameters and defaults](#parameters-and-defaults-1)
  - [Caveats](#caveats-1)

---

## Bond Geometry — the triplet engine

### What this page shows

The **Bond Geometry** tab (component
[BondGeometryPage.jsx](../../web_app/frontend/src/components/BondGeometryPage.jsx)) answers:

> At every B atom, what angle do its A and C neighbours subtend — and how is that distributed
> over the whole configuration?

The workflow mirrors `triplets_new_bonds_sinth` from the RMCProfile tool set: two atoms are
*bonded* when their distance falls inside an inclusive window you choose by eye against the
partial $g(r)$, and every (A-bond, C-bond) pair at a central B contributes one angle

$$\theta \;=\; \arccos\frac{\mathbf r_{BA}\cdot\mathbf r_{BC}}
{\lVert\mathbf r_{BA}\rVert\,\lVert\mathbf r_{BC}\rVert}\;\in\;[0°,180°].$$

The engine of record is [rmc_toolkits/triplets.py](../../rmc_toolkits/triplets.py) (NumPy; the
module docstring is the compact math reference). It is served three ways, all built on one pass:

| Boundary | Entry point | Notes |
|---|---|---|
| Flask | `/api/triplets` in [app.py](../../web_app/backend/app.py) | `cached_bond_angle_summary`, an `lru_cache(16)` keyed on (path, mtime, every parameter) |
| Browser | `kind: 'triplets'` in [pcaKdeWorker.js](../../web_app/frontend/src/workers/pcaKdeWorker.js), engine [workers/triplets.js](../../web_app/frontend/src/workers/triplets.js) | line-for-line port; answers from the worker's cached parse of the already-loaded `.rmc6f` |
| CLI | `rmc-triplets` ([triplets_cli.py](../../rmc_toolkits/triplets_cli.py)) | commented CSV + optional PNG and raw angle list |

The routing switch is `requestPca` in
[useSiteCloud.js](../../web_app/frontend/src/useSiteCloud.js), and — as on the PCA and
Orientation pages — it is **not** a Flask-vs-static-build switch: whenever an `.rmc6f` has been
loaded as a browser file (the Demo, or a picked folder), the worker answers *in both runtimes*.
Only a typed backend directory goes through HTTP. `bond_angle_summary` defines the payload
contract all three share.

Unlike the KDE and PCA pages there is **no subsampling and no randomness anywhere in this
engine**: every atom participates, every image within reach is visited, and the result is
deterministic to the last count.

---

### Step 1 — Inputs, validation, and element matching

`_triplet_core` ([triplets.py](../../rmc_toolkits/triplets.py)) receives:

- `fractional` — $(N,3)$ coordinates as **fractions of the supercell** (the `.rmc6f` storage
  convention), mapped to Cartesian ångström by the row-vector product
  $\mathbf x = \mathbf f\,\mathsf L$ with $\mathsf L$ the $(3,3)$ supercell lattice rows — the
  same convention as `pca_kde.py`.
- `elements` — matching symbols, compared after `str.capitalize()`, the same normalization
  [parsers.iter_rmc6f_atoms](../../rmc_toolkits/parsers.py) applies — so `se`, `SE` and `Se`
  name one species.
- `triplet = (A, B, C)` with **B central**, `bond12 = (rmin, rmax)` for A–B, and an optional
  `bond23` for B–C (`None` reuses `bond12`).

Validation is strict and raises `ValueError` rather than coercing: coordinates must be finite
$(N,3)$, the lattice a finite $(3,3)$ matrix, each window needs $0 \le r_\mathrm{min} <
r_\mathrm{max}$ with finite bounds, and a triplet element with no atoms in the configuration
reports the list of symbols that *are* available. The port
([workers/triplets.js](../../web_app/frontend/src/workers/triplets.js)) additionally rejects
`null`/`''` bounds explicitly, because JavaScript's `Number(null)` would otherwise silently
coerce a missing bound to `rmin = 0` where Python's `float()` raises.

Two flags fall out of the spec before any geometry runs:

- `shared_ends` — true iff A and C are the same element **and** the two windows are equal. It
  selects between the two counting rules of Step 5 and lets the B–C search be skipped entirely
  (`bonds23 = bonds12`).
- The bin count for a requested width (Step 6) is fixed here too:
  `_bin_count = max(1, floor(180/w + 0.5))` — **round half up**, matching JavaScript's
  `Math.round` exactly. Python's banker's `round()` would disagree at exact-.5 ratios (an
  8° request: $180/8 = 22.5 \to 23$ bins in both engines, not 22 vs 23), which would be a
  different *payload shape* across runtimes, not merely different numbers.

### Step 2 — Folding and the image bookkeeping

Every coordinate is folded into $[0,1)$ first:

$$\mathbf f' = \mathbf f - \lfloor\mathbf f\rfloor$$

An RMC configuration may arrive pre-wrapped or with atoms drifted outside the box; folding plus
the explicit image bookkeeping below restores the true relative geometry either way, so both
inputs give identical results (`test_wrap_invariance` in
[tests/test_triplets.py](../../tests/test_triplets.py) pins this).

From here on, a neighbour is always identified as **(atom row, integer image)** — the candidate
atom plus the whole-box shift $\mathbf m \in \mathbb Z^3$ applied to it. That pair is what makes
two bonds "the same bond" in Step 5's exclusions, and it is exact: no minimum-image convention
is assumed anywhere.

### Step 3 — The linked-cell search: cells, reach, and correctness

`_neighbor_bonds` divides the box into $n_i$ cells per lattice direction. The length scale is
the **perpendicular width** of the box along each direction (`_perpendicular_widths`):

$$w_i \;=\; \frac{V}{\lvert \mathbf a_j \times \mathbf a_k \rvert},\qquad
V = \lvert \mathbf a_1 \cdot (\mathbf a_2 \times \mathbf a_3) \rvert$$

— the distance between the two box faces spanned by the other two vectors, which is what decides
both how many cells fit along direction $i$ and how far a periodic image can reach. Then

$$n_i = \min\bigl(\max(1, \lfloor w_i / r_\mathrm{max} \rfloor),\; 64\bigr),\qquad
k_i = \left\lceil \frac{r_\mathrm{max}\, n_i}{w_i} + 10^{-9} \right\rceil$$

with $k_i$ the number of neighbour-cell layers visited along direction $i$. The cap of 64 cells
per axis (`MAX_CELLS_PER_AXIS`) only matters for tiny $r_\mathrm{max}$ in huge boxes, where finer
cells cost memory without reducing work; the $10^{-9}$ headroom (`REACH_HEADROOM`) absorbs float
rounding at exact-integer ratios at the cost of one spurious (empty) extra layer in those rare
cases.

**Why this covers every pair.** Two atoms in cells $q$ layers apart along direction $i$ are
separated by strictly more than $(q-1)$ perpendicular cell thicknesses $w_i/n_i$. So any pair
within $r_\mathrm{max}$ sits within $k_i = \lceil r_\mathrm{max} n_i / w_i \rceil$ layers — for
**any** cell shape, triclinic included. And because the layer offsets are enumerated as absolute
shifts (every $(\delta_1,\delta_2,\delta_3)$ with $\lvert\delta_i\rvert \le k_i$) rather than
wrapped indices, the search is also correct when the box is *smaller* than $r_\mathrm{max}$: an
out-of-range cell index wraps, records the whole-box shift it wrapped by,

$$\mathbf m = \frac{\text{shifted} - \text{wrapped}}{\mathbf n}\ \ (\text{integer division}),$$

and multiple images of the *same* atom become genuine distinct neighbours
(`SmallBoxImageTests`). The candidate's image is placed beside the center in fractional space,
then mapped to Cartesian:

$$\Delta\mathbf x = \bigl(\mathbf f_\text{cand} + \mathbf m - \mathbf f_\text{center}\bigr)\,\mathsf L$$

(One implementation note, recorded in the source: the fractional→Cartesian product uses
`np.einsum` rather than `@`, because NumPy on Apple's Accelerate BLAS emits spurious
divide/overflow warnings for large $(P,3)\times(3,3)$ matmuls of finite values. Results are
identical; einsum skips that code path.)

The search is fully vectorized per offset: candidates are bucketed by flattened cell index once
(`argsort` + `bincount` + prefix sums), and each of the $\prod_i (2k_i{+}1)$ offsets gathers all
its (center, candidate) pairs in one shot via `_ragged_ranks`, a flattened 0..k−1 rank within
consecutive groups.

### Step 4 — Which pairs are bonds

A candidate pair with squared distance $d^2$ survives iff

$$r_\mathrm{min}^2 \le d^2 \le r_\mathrm{max}^2
\quad\text{and}\quad d^2 > 0
\quad\text{and not}\quad(\text{same row} \wedge \mathbf m = \mathbf 0).$$

Three deliberate choices:

- **Windows are inclusive at both ends** — read the bounds off the partial $g(r)$ and atoms at
  exactly the bound count.
- **A zero-length pair is never a bond, even under `rmin = 0`.** Bitwise-coincident atoms have
  no direction; admitting the pair would put a $0/0$ NaN into every angle it joins
  (`CoincidentAtomTests`).
- **A center is never its own neighbour in the unshifted image** ($\mathbf m = \mathbf 0$), but
  its *other* images are genuine neighbours and stay — a B atom in a small box legitimately
  bonds to its own periodic copy (`test_self_image_bonds_of_the_central_element`).

Each surviving bond is stored as (center position in the selection, candidate row, image
$\mathbf m$, Cartesian vector, length) — the `_Bonds` container.

### Step 5 — Pairing bonds into angles

`_pair_angles` sorts both bond lists by central atom (stable sort, so order is deterministic)
and forms the per-center pairing:

- **Shared ends** (`shared_ends = True`; one list paired with itself): keep the strict upper
  triangle $i < j$, so each *unordered* pair of distinct bonds counts once. An octahedrally
  coordinated B with six bonds gives exactly $\binom{6}{2} = 15$ angles — $12\times 90° +
  3\times 180°$ (`OctahedronTests`, and the same invariant asserted from the JS side).
- **Distinct ends or windows**: every (A-bond, C-bond) combination counts — *ordered* (1→2, 2→3)
  assignments — **minus** combinations where both bonds reach the same atom row in the same
  image $\mathbf m$: the degenerate zero-degree "angle" of a bond with itself, which arises when
  A and C name the same element with overlapping windows
  (`SameElementDistinctWindowTests`).

The angle is then

$$\theta = \frac{180°}{\pi}\arccos\!\Bigl(\operatorname{clip}\bigl(
\tfrac{\mathbf v_1\cdot\mathbf v_2}{\lVert\mathbf v_1\rVert\lVert\mathbf v_2\rVert},\,-1,\,1\bigr)\Bigr)$$

with the clip guarding the $\pm1$ boundary against rounding. There is no tolerance anywhere
else: the whole calculation is exact geometry on float64.

### Step 6 — The histogram and its three normalizations

Angles are binned uniformly over $[0°, 180°]$ into $K = \max(1,\lfloor 180/w_\text{req} +
0.5\rfloor)$ bins of realized width $w = 180/K$ (a requested width that does not divide 180 is
adjusted to the nearest exact tiling). Alongside the raw `counts` $N_k$ the result carries:

**`density`** — a per-degree probability density with unit integral over $[0,180]$:

$$D_k = \frac{N_k}{N\,w},\qquad \sum_k D_k\, w = 1 .$$

Note that randomly oriented bonds do **not** give a flat density: there are simply fewer ways to
form an angle near 0° or 180° than near 90°, so geometry alone bows the curve as $\sin\theta$.

**`sin_corrected`** — the count fraction divided by the *exact* isotropic reference fraction per
bin:

$$S_k = \frac{N_k/N}{\bigl(\cos\theta_k - \cos\theta_{k+1}\bigr)/2},\qquad
\int_{\theta_k}^{\theta_{k+1}} \tfrac{1}{2}\sin\theta\,d\theta
= \tfrac{\cos\theta_k - \cos\theta_{k+1}}{2}.$$

For bonds pointing in independent uniformly-random directions this is flat at $1.0$ — the
RMCProfile `sinth` view — so anything above 1 is real structure, and a peak near 180° (the
octahedral *trans* angle) is no longer suppressed by geometry. Dividing by the **bin integral**
rather than by $1/\sin\theta_c$ at the bin center is what keeps the 0° and 180° bins finite,
where $1/\sin\theta_c$ diverges.

With zero angles both curves are all-zero rather than NaN.

### Step 7 — The summary payload

`bond_angle_summary` re-runs nothing: one `_triplet_core` pass feeds every panel of the page.
The dict (camelCase keys, plain lists and scalars — JSON-safe) is the payload contract shared by
the Flask route and the worker; the port's `bondAngleSummary` must be kept in sync with any
change here. On top of the three angle curves it adds:

| Key | Content |
|---|---|
| `triplet`, `bond12`, `bond23`, `sharedEnds`, `binWidth` | the resolved spec (realized bin width, not the requested one) |
| `angleCount`, `meanAngle`, `stdAngle`, `apexCount` | angle statistics and the central-atom count; means are `None` when empty |
| `lengths12`, `lengths23` | bond-length histograms **inside each window**: fixed `LENGTH_BINS = 40` bins over the window (fixed count, not width, so any window renders at the same detail), plus `count` and `meanLength`. `lengths23` is `null` under shared ends — it would duplicate `lengths12` |
| `coordination` | `coordination[n]` = how many central atoms have exactly $n$ window-1 bonds — a double `bincount` of the per-center bond counts |

### Step 8 — The two app boundaries and their caps

The engine itself is **unrestricted** — library and CLI callers can ask for anything. The two
app boundaries apply identical request caps, because one request could otherwise blow up the
image stencil (volume grows $\sim r_\mathrm{max}^3$ and the candidate pair count with it) or the
response size — and in the browser it would freeze the shared PCA worker:

| Cap | Flask `/api/triplets` | Worker `kind: 'triplets'` |
|---|---|---|
| $r_\mathrm{max} \le 15\,$Å (each window) | 400 | thrown `Error` |
| `binWidth` $\ge 0.05°$ | 400 | thrown `Error` |
| `r23Min`/`r23Max` required together | 400 | thrown `Error` |

Request parameters are flat scalars (`end1`, `apex`, `end2`, `r12Min`, `r12Max`, `r23Min`,
`r23Max`, `binWidth`) so the identical request shape works as an HTTP query string and as a
worker message. The Flask side normalizes element case **before** the cache key
(`'se'` and `'Se'` share one entry) and resolves the `bond23` default before the call, so equal
windows hit one cache entry; the cache is `lru_cache(maxsize=16)` keyed on (path, mtime, every
parameter), mirroring `pca_kde.cached_site_displacements`. Errors map to 400 (bad parameters),
403 (path escapes the data root), 404 (no `.rmc6f`).

### The `rmc-triplets` CLI

Console entry point installed by `pip install -e .`
([triplets_cli.py](../../rmc_toolkits/triplets_cli.py); module form
`python -m rmc_toolkits.triplets_cli`). Accepts an `.rmc6f` file or a run folder (first sorted
match), writes a commented CSV (`angle_deg, counts, density_per_deg, sin_corrected` with the
spec, bond counts and mean lengths in `#` headers), optionally a PNG plot (`--plot`, Agg
backend, sin-corrected + density on twin axes) and the raw angle list (`--angles-out`). Nothing
is overwritten without `--force`.

```bash
rmc-triplets data/5K_try1 --triplet Se Nb Se --bond12 2.2 2.9 --plot se_nb_se.png
rmc-triplets config.rmc6f --triplet O Ti O --bond12 1.7 2.3 --bond23 1.7 2.3 --bin-width 0.5
```

### Parameters and defaults

| Parameter | Default | Meaning |
|---|---|---|
| `triplet` (A, B, C) | — (required) | element symbols, B central; matched after `capitalize()` |
| `bond12` | — (required) | inclusive A–B window (Å), $0 \le r_\mathrm{min} < r_\mathrm{max}$ |
| `bond23` | `None` → `bond12` | inclusive B–C window; equal ends + equal windows ⇒ unordered counting |
| `bin_width` | `1.0`° | requested width; realized width is $180/\max(1,\lfloor 180/w+0.5\rfloor)$ |
| `collect_angles` | `False` | `bond_angle_distribution` only: keep the raw angle list |
| `MAX_CELLS_PER_AXIS` | 64 | linked-cell resolution cap per lattice direction |
| `REACH_HEADROOM` | $10^{-9}$ | relative headroom on the layer reach against float rounding |
| `LENGTH_BINS` | 40 | bond-length histogram bins per window (summary payload only) |
| App-boundary caps | $r_\mathrm{max}\le15$ Å, `binWidth` ≥ 0.05° | Flask route and worker only; engine and CLI unrestricted |

### Parity: Python engine vs JavaScript port

[workers/triplets.js](../../web_app/frontend/src/workers/triplets.js) is a line-for-line port of
the Python engine, and the parity is pinned by golden fixtures rather than claimed:
[tests/generate_triplets_fixture.py](../../tests/generate_triplets_fixture.py) evaluates the
Python engine on two constructed configurations and writes
[triplets_fixture.json](../../web_app/frontend/src/__tests__/fixtures/triplets_fixture.json),
which [workers/\_\_tests\_\_/triplets.test.js](../../web_app/frontend/src/workers/__tests__/triplets.test.js)
replays against the port:

| Fixture case | What it exercises |
|---|---|
| `random-triclinic` — 48 atoms, seeded RNG, lattice $[[6,0,0],[3,5,0],[1,1,7]]$; three specs incl. shared ends, distinct windows, and B = A = C | general triclinic geometry, both counting rules, 5° and 2° bins |
| `small-box-images` — 3 atoms in a 4 Å cube with windows reaching 3.5 Å | multiple periodic images of one atom as distinct neighbours |

Measured agreement, asserted per bin and per statistic:

- **Histogram counts, coordination, length-histogram counts: exact integer equality.** No
  tolerance — the two engines must produce the same integers.
- `density`, `sinCorrected`, bin centers: $10^{-9}$ (single-pass float math; covers
  summation-order noise).
- `meanAngle`, `stdAngle`, `meanLength`: $10^{-7}$.
- Sorted raw angles (head/tail samples): $10^{-5}$ — these go through `acos` twice
  (compute, then fixture rounding).

The Python engine itself is pinned to a brute-force all-images reference over a $\pm2$ image
span on random triclinic configurations (`TriclinicBruteForceTests` in
[tests/test_triplets.py](../../tests/test_triplets.py)), plus the constructed invariants:
octahedron counting, bonds through the periodic wall, wrap invariance, zero-length exclusion,
self-image bonds, and ordered-counting degeneracies. The backend route's caps and error paths
are covered by `TripletsApiTests` in [tests/test_backend_api.py](../../tests/test_backend_api.py).

**The one documented residual divergence:** libm and V8 `acos` may differ by 1 ulp, so a
*bitwise-ideal* geometry whose cosine lands exactly on a bin edge (e.g. $\cos\theta = 0.5$ in an
undisplaced average configuration) can shift one count into the neighbouring bin between
engines. Real RMC configurations — displaced by construction — never hit this.

### Caveats

- **The angle histogram is unweighted.** Every triplet counts 1; there is no amplitude,
  distance, or multiplicity weighting of any kind.
- **`coordination` counts window-1 bonds only.** Under distinct windows it is the A-coordination
  of B; the C-side coordination is not reported.
- **The realized bin width may differ from the request** (nearest exact tiling of 180°). The
  payload reports the realized width; the CSV bin centers are authoritative.
- **Inclusive windows mean boundary atoms count.** Two runs whose $g(r)$ peak touches the bound
  can differ by exactly the boundary population — intentional, but worth knowing when comparing.
- **The engine reports geometry, not chemistry.** A "bond" is a distance window and nothing
  else; there is no bond-valence, electronegativity, or connectivity analysis.

---

## Bond Geometry — the page

### What the page owns

[BondGeometryPage.jsx](../../web_app/frontend/src/components/BondGeometryPage.jsx) renders one
controls bar (triplet selects, windows, `Distinct B–C` switch, bin width, **Compute**), the
model-information card (same `ModelSummary` as the Dashboard, symmetry card omitted), a
`Triplet result` chip strip, and three equal-width panels: the angle distribution, the
partial-$g(r)$ window helper, and the folded-cell bond view
([FoldedCellPanel.jsx](../../web_app/frontend/src/components/FoldedCellPanel.jsx)). Everything
numerical comes from the Step 7 payload; the page adds selection, presentation, and the two
helper views.

### Step 1 — Triplet seeding

The three selects are seeded **once per sites payload**, not per dataset key: on a dataset
switch the new key arrives while the previous run's sites are still in state, so keying on the
data itself waits for the real element list — and a Live Data refresh keeps the user's picks
when they still apply. The seed ranks elements by total atom count from the sites table: ends =
most abundant (in practice the anion), central = next most abundant — Se–Nb–Se for GaNb₄Se₈. An
existing valid selection is never overwritten.

### Step 2 — The compute request and the epoch guard

**Compute** issues `requestPca('triplets', {end1, apex, end2, r12Min, r12Max, r23Min?, r23Max?,
binWidth})` — the B–C window included only when the split switch is on (the engine then receives
`bond23 = null` and reuses `bond12`). A dataset switch clears any previous result immediately
and bumps a `runEpoch` ref; a compute that was in flight for the old run compares its captured
epoch on resolve and can never land a stale payload on the new dataset.

### Step 3 — The result chips

Straight reads of the payload: central-atom count (`apexCount`), bonds with mean length
(`lengths12`, and `lengths23` when not shared), the coordination summary — mean bonds per B
$\sum_n n\,c_n / \sum_n c_n$, plus the modal $n$ and its share — and the angle count with
mean ± std.

### Step 4 — The angle plot and the `fit` variant

The hero plot shows `sinCorrected` or `density` (toggle; sin-corrected is the default). Both
this plot and the partial-$g(r)$ helper render through
[InteractivePlot](../../web_app/frontend/src/components/InteractivePlot.jsx) with
`variant="fit"`: the SVG `viewBox` is taken from the rendered box (ResizeObserver, rounded to
whole pixels) instead of the fixed 8:5 aspect, one user unit = one CSS pixel, and tick density
follows the box (~1 y-tick per 70 px, ~1 x-tick per 95 px, clamped). The two cards are handed
identically sized boxes — pinned header height, two reserved legend rows
([BondGeometryPage.css](../../web_app/frontend/src/components/BondGeometryPage.css)) — so the
two figures always share one aspect ratio.

### Step 5 — The partial-g(r) window helper

The helper plots the run's measured partial pair distribution from `PDFpartials.csv` (loaded via
the same plot-file path as the Dashboard; both runtimes) so the windows can be set against the
actual first shell:

- **Curves**: the A–B partial always; a second curve whenever A–B and B–C are *different bond
  types* (pair labels looked up in either order — `Ta-Se` matches `Se-Ta`). This is independent
  of the window split: Ga–Ta–Se has two shells to bracket even with one shared window, and
  Se–Ta–Se has one shell even with two windows.
- **Guides**: dashed verticals at the current bounds, following the inputs after a 400 ms
  debounce so the plot's view state does not reset per keystroke. With the split **off**, one
  neutral-grey pair covers both bonds. With it **on**, each window gets its own pair, labelled
  `A–B rmin/rmax` and `B–C rmin/rmax` **by role** (pair names would collide for same-element
  triplets) and colored to match the curve each brackets: guides consume no palette slot, so
  shell $N$ is `PLOT_PALETTE[N]` ([plotPalette.js](../../web_app/frontend/src/plotPalette.js)).
- **Crop**: the x-range is cut at $\max(6\,\text{Å},\ 2\times$ the furthest active
  $r_\mathrm{max})$ — beyond the first-shell region nothing informs a bond window.

The panel is display-only: the page never computes a $g(r)$; a run without `PDFpartials.csv`
gets an empty panel and everything else still works.

### Step 6 — The folded-cell bond view

[FoldedCellPanel.jsx](../../web_app/frontend/src/components/FoldedCellPanel.jsx) shows the same
folded unit cell as the Atomic Density page — every atom of the supercell folded into one cell
as a point cloud, colored by element, so the spread around a site is the *measured* thermal
cloud rather than a fitted ellipsoid — with the analysis' detected bonds drawn over it:

- **Cloud**: one `THREE.Points` per element; above 120 000 atoms the cloud is strided down to
  ~120 k points (display only — the engine always sees every atom).
- **Bonds**: for each computed window, average-site pairs whose distance falls inside it,
  periodic images included ($\mathbf m \in \{-1,0,1\}^3$, same-element in-cell pairs
  deduplicated) — so a stick may reach an image just outside the box, which is the real
  coordination. Drawn as thin transparent `LineSegments` (opacity 0.45), A–B in the app accent
  blue (`0x2563eb`), B–C in amber (`0xd97706`) when the windows are distinct, so a full network
  reads as a framework without hiding the cloud.
- The a/b/c gizmo, reset view, and 1×/3× PNG export follow the other 3D panels.

Note the sticks connect **average site positions** (the folded reference sites), while the cloud
shows instantaneous atoms: a stick is the average bond, not any single configuration's bond.

### Parameters and defaults

| Control | Default | Notes |
|---|---|---|
| Triplet A, B, C | seeded per sites payload | ends = most abundant element, central = next |
| A–B window | 2.0 – 3.0 Å | inclusive; string state, validated at the boundary |
| Distinct B–C | off | off ⇒ B–C reuses the A–B window and one guide pair |
| B–C window | 2.0 – 3.0 Å | only sent when the split is on |
| Bin width | 1.0° | realized width comes back in the payload |
| Angle view | sin-corrected | toggle to per-degree density |
| Guide debounce | 400 ms | `useDebounced` on all four window inputs |
| Cloud stride cap | 120 000 points | `MAX_CLOUD_POINTS`, display only |

### Caveats

- **The page computes no geometry.** Every number on it is the engine payload; the helper and
  the folded cell are presentation over `PDFpartials.csv` and the sites table respectively.
- **The bond-length histograms are in the payload but not plotted.** An earlier layout gave them
  a panel; it duplicated the first-shell peak the partial $g(r)$ already shows, clipped to the
  window. The counts and mean lengths survive in the result chips.
- **The folded-cell sticks are average-structure bonds.** They match sites within the window at
  their *average* positions; a strongly displaced site whose average distance falls outside the
  window shows no stick even though many instantaneous bonds were counted (and vice versa).
- **The helper needs `PDFpartials.csv`.** Runs without it lose the guides' context but nothing
  else; the windows still apply exactly as typed.
