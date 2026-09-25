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
`null`, `''` and whitespace-only bounds explicitly (`isBlankValue`), because JavaScript's
`Number()` turns every one of them into `0` — a missing bound would silently become
`rmin = 0` — where Python's `float()` raises. The guard only protects callers that hand the
engine the raw value, so neither app boundary converts a bound before the check (Step 8), and
the page validates its input boxes before sending anything (page Step 2).

Two flags fall out of the spec before any geometry runs:

- `same_end` — true iff A and C are the same element. It selects between the two counting rules
  of Step 5, and it makes the engine run **one** neighbour search over both windows (each bond
  tagged with the window(s) it falls in) instead of one search per end.
- `shared_ends` — `same_end` **and** equal windows. The B–C bond list is then the A–B list
  (`bonds23 = bonds12`), and the payload reports it as `sharedEnds` (one length histogram).
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

$$\Delta\mathbf x = \bigl((\mathbf f_\text{cand} - \mathbf f_\text{center}) + \mathbf m\bigr)\,\mathsf L$$

The order of operations is deliberate: the fractional difference is taken **before** the integer
shift is added, so the same bond seen from its other end,
$((\mathbf f_\text{center} - \mathbf f_\text{cand}) - \mathbf m)\,\mathsf L$, is the *exact*
negative (IEEE rounding is symmetric under negation) and has bitwise the same length. A bond is
therefore inside or outside a window from both of its ends alike. Before 1.0 the shift was added
first, and on ideal lattices with a window bound exactly on a shell the two ends could disagree
(one bond counted from one end only — an odd directed count, e.g. 2675 on a 5×5×5 two-atom
cubic lattice with $r_\mathrm{max} = a$).

(One implementation note, recorded in the source: the fractional→Cartesian product, the
squared length and the angle's dot product are written out term by term —
$(\Delta f_1 L_{1j} + \Delta f_2 L_{2j}) + \Delta f_3 L_{3j}$ — in exactly the JavaScript port's
evaluation order, so both engines produce bitwise-identical vectors, and no BLAS matmul is
involved: NumPy on Apple's Accelerate BLAS emits spurious divide/overflow warnings for large
$(P,3)\times(3,3)$ matmuls of finite values.)

The search is vectorized per offset: candidates are bucketed by flattened cell index once
(`argsort` + `bincount` + prefix sums), and each of the $\prod_i (2k_i{+}1)$ offsets gathers its
(center, candidate) pairs in one shot via `_ragged_ranks`, a flattened 0..k−1 rank within
consecutive groups. Centers are taken in blocks of $\lfloor$`SEARCH_CHUNK`$/\max(\text{atoms per
cell})\rfloor$ ($2^{20}$ candidate pairs per offset, ~100 MB of transients), so the search's
working memory is bounded whatever $r_\mathrm{max}$ and the box; the per-center bond order
(stencil order) is the same with or without blocking. The same search runs in a
**count-only** mode that stores nothing and returns, per center, how many bonds fall in each
window-membership pattern — the exact counting pass of Step 8.

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

`_pair_angles` sorts the bond lists by central atom (stable sort, so order is deterministic)
and forms the per-center pairing:

- **Same end element** (`same_end`, A = C): the one search of Step 1 lists every bond image
  $x$ once, with flags $x \in w_{12}$ and $x \in w_{23}$. The list is paired with itself over the
  strict upper triangle $i < j$, so each *unordered* pair of distinct bond images — one physical
  triplet $\{x, B, y\}$ — is considered once, and it contributes one angle iff
  $$(x \in w_{12} \wedge y \in w_{23}) \;\vee\; (y \in w_{12} \wedge x \in w_{23}).$$
  With equal windows that is every unordered pair: an octahedrally coordinated B with six bonds
  gives exactly $\binom{6}{2} = 15$ angles — $12\times 90° + 3\times 180°$ (`OctahedronTests`,
  and the same invariant asserted from the JS side). With distinct windows the rule is
  **continuous**: moving a bound changes the count only by the triplets whose bonds cross it,
  and at $w_{23} \to w_{12}$ it reduces to the shared-window count
  (`SameElementDistinctWindowTests`). Disjoint windows (short vs long bonds) pair each short
  bond with each long one once. A bond is never paired with itself, since $i < j$.
- **Different end elements** (A ≠ C): every (A-bond, C-bond) combination counts; the A and C
  atoms are necessarily different atoms.

> Before 1.0 the same-element/distinct-window case counted *ordered* (1→2, 2→3) assignments, so
> a triplet whose two bonds both lay in the overlap of the windows counted twice, and nudging a
> B–C bound by $10^{-4}$ Å off the A–B one doubled every count (223 651 → 447 317 Se–Nb–Se
> angles on the 5 K sample). The unordered rule counts each physical triplet once.

The angle is then

$$\theta = \frac{180°}{\pi}\arccos\!\Bigl(\operatorname{clip}\bigl(
\tfrac{\mathbf v_1\cdot\mathbf v_2}{\lVert\mathbf v_1\rVert\lVert\mathbf v_2\rVert},\,-1,\,1\bigr)\Bigr)$$

with the clip guarding the $\pm1$ boundary against rounding. There is no tolerance anywhere
else: the whole calculation is exact geometry on float64.

**The exact angle count comes first.** Before any angle is formed, the count per center follows
from the bond lists alone (`_angles_per_center`):

$$N_\text{angles} = \sum_B \begin{cases}
n_{12}\,n_{23} & A \ne C,\\[2pt]
\binom{a+b+c}{2} - \binom{a}{2} - \binom{b}{2} & A = C,
\end{cases}$$

with, for a same end element, $a$ bonds only in $w_{12}$, $b$ only in $w_{23}$ and $c$ in both
(equal windows: $\binom{c}{2}$). It is what the work budget of Step 8 is checked against.

**Angles are streamed, never all held.** `_stream_angles` walks the centers in consecutive runs
whose combined bond-pair count stays within `PAIR_CHUNK` $= 2^{18}$ (a single center above it is
a run of its own), forms that run's angles, adds them to the histogram, and folds their mean and
variance into the running totals with Chan's parallel update
($\delta = \bar x_\text{chunk} - \bar x$, $M_2 \mathrel{+}= M_{2,\text{chunk}} + \delta^2 n\,n_\text{chunk}/(n+n_\text{chunk})$).
Pairing memory is therefore ~25 MB whatever the angle count; the raw list is kept only when
`collect_angles` asks for it. The JS port streams inside its pairing loop (Welford per angle)
and keeps angles only for `collectAngles`. Before 1.0 both engines materialized every angle
(~200 B/angle in NumPy index arrays, one JS array capped by V8 at $2^{27}$ elements), so a
window well inside the app's 15 Å cap needed tens to hundreds of GB in Flask and threw
`RangeError: Invalid array length` in the worker.

### Step 6 — The histogram and its three normalizations

Angles are binned uniformly over $[0°, 180°]$ into $K = \max(1,\lfloor 180/w_\text{req} +
0.5\rfloor)$ bins of realized width $w = 180/K$ (a requested width that does not divide 180 is
adjusted to the nearest exact tiling).

**Bin membership** (`_angle_bins`; `angleBin` in the port). Bins are half-open,
$[\theta_k, \theta_{k+1})$ with $\theta_k = k\,w$ (the `linspace` values), the last one closed
at 180° — numpy.histogram's convention, written out identically in both engines. Every edge is a
multiple of $w$, so for the usual widths (1, 0.5, 2, 3, 5 …) **every symmetry angle of an
undisplaced configuration sits exactly on an edge** — 60/90/120/180° of an ideal perovskite or
fcc lattice, the Nb₄ tetrahedron's 60°, a CIF-built RMCProfile start configuration — and float
noise puts the computed angle a few ulp (~$10^{-13}$°) either side of it. Left alone, that noise
split one symmetry class between two bins in an arbitrary ratio that changed under a rigid shift
of the configuration, and differently in the two engines (numpy's `arccos` and V8's
`Math.acos` differ by 1 ulp on ~17% of inputs: `arccos(0.5)` gives 59.99999999999999°, V8
60.00000000000001°). On a 3×3×3 ideal SrTiO₃, O–Sr–O at 1° bins came out
`{59: 354, 60: 294, …}` in Python and `{59: 273, 60: 375, …}` in JS. So an angle within
`EDGE_SNAP_DEG` $= 10^{-9}$° of an edge is binned **as exactly on it**:

$$k^\ast = \lfloor \theta / w + \tfrac12 \rfloor,\qquad
\lvert \theta - k^\ast w \rvert < 10^{-9}° \;\Rightarrow\; \text{bin } \min(k^\ast, K-1),$$

i.e. into the bin the edge starts. A symmetry class then lands whole in one bin, deterministic
and identical in both engines. The tolerance is four orders of magnitude above the noise and
far below any bin width or real displacement (an angle $10^{-7}$° below an edge still bins
below it). Only the binning snaps: `meanAngle`, `stdAngle` and raw angles keep the computed
values. On displaced RMC configurations nothing changes (the 5 K sample and its AVERAGE file
bin identically to `np.histogram`).

Alongside the raw `counts` $N_k$ the result carries:

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
| `lengths12`, `lengths23` | bond-length histograms **inside each window**: fixed `LENGTH_BINS = 40` bins over the window (fixed count, not width, so any window renders at the same detail), plus `count`, `uniqueBonds` and `meanLength`. `count` is the number of **B-centred bond vectors** (the histogram total); `uniqueBonds` the number of **physical bonds**, each once. They differ when the end element is the central element (A = B for `lengths12`, C = B for `lengths23`): every such bond is found from both of its ends, so `uniqueBonds = count / 2` exactly (Step 3's antisymmetry makes the halving exact; a self-image bond to $\pm\mathbf m$ is one periodic bond). Otherwise `uniqueBonds = count`. On the 5 K sample, Nb–Nb–Nb 2.6–3.4 Å has `count` 46 704 and `uniqueBonds` 23 352 — the number of distinct Nb–Nb pairs an independent `cKDTree` search finds. `lengths23` is `null` under shared ends — it would duplicate `lengths12` |
| `coordination` | `coordination[n]` = how many central atoms have exactly $n$ window-1 bonds — a double `bincount` of the per-center bond counts |

### Step 8 — The two app boundaries and their caps

The engine itself is **unrestricted** — library and CLI callers can ask for anything (and, with
Step 5's streaming, memory stays bounded; only the time grows). The two app boundaries apply
identical request caps, each bounding a different cost — and in the browser an unbounded request
would freeze the shared PCA worker:

| Cap | Bounds | Flask `/api/triplets` | Worker `kind: 'triplets'` |
|---|---|---|---|
| $r_\mathrm{max} \le 15\,$Å (each window) | the neighbour search (bond count $\sim r_\mathrm{max}^3$) | 400 | thrown `Error` |
| exact angle count $\le$ `APP_MAX_ANGLES` $= 5\times10^7$ | the pairing work (angle count $\sim r_\mathrm{max}^6$) | 400 | thrown `Error` |
| `binWidth` $\ge 0.05°$ | the response size | 400 | thrown `Error` |
| `r12Min`/`r12Max` present, `r23Min`/`r23Max` both or neither — `null`, `''` and whitespace count as missing | — | 400 | thrown `Error` |

A missing bound is an error with the same text at both boundaries ("r12Min/r12Max are required
together; missing r12Min") — never a bound of 0. Before 1.0 the worker ran `Number()` on every
bound, so a cleared minimum reached the engine as `0` and silently widened the window (Nb–Nb–Nb
3.5–4.6 Å became 0–4.6 Å on the 5 K sample: 235 883 angles, mean 107°, instead of 48 182 at
62°), while Flask rejected the same request.

The rmax cap alone does **not** bound the work: on the 52 000-atom 5 K sample a Se–Nb–Se window
2.2–15 Å forms $1.27\times10^9$ angles (Se–Se–Se 2–15 Å: $2.45\times10^9$). So both boundaries
pass the one shared budget — `APP_MAX_ANGLES` in [triplets.py](../../rmc_toolkits/triplets.py),
mirrored in [workers/triplets.js](../../web_app/frontend/src/workers/triplets.js) and pinned
equal by the parity fixture — as `max_angles` / `maxAngles`. A budgeted request first runs the
count-only search (Step 3), which stores nothing, computes the exact count of Step 5, and
refuses a spec above the budget with a message naming the count ("… would form 1,274,044,098
angles, over the limit of 50,000,000 for one request; narrow the bond windows"). Measured on the
5 K sample: refusing Se–Se–Se 2–15 Å costs ~0.2 GB and <1 s in the worker (~0.5 GB, ~8 s in
Flask, parse included); an accepted Se–Nb–Se 2.2–8 Å request ($2.9\times10^7$ angles) takes
~0.8 s / 0.4 GB in the worker and ~2 s / 0.5 GB in Flask.

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
`python -m rmc_toolkits.triplets_cli`). Accepts an `.rmc6f` file or a run folder. A run folder
resolves to **the configuration the app analyses** (`find_run_configuration`, the same rule as the
backend's `_find_rmc6f` and the browser's `chooseStructureFile`): the `.rmc6f` whose stem matches
the run's own outputs (`<stem>-NN.log` first, then `<stem>_PDFpartials.csv`, `_FQ1.csv`, …), the
first sorted file only when nothing matches — so an input supercell `GaNb4Se8.rmc6f` beside the
refined `GaNb4Se8_5K.rmc6f` no longer wins by sorting first (`'.'` < `'_'`). The CLI prints the
path it chose. It writes a commented CSV (`angle_deg, counts, density_per_deg, sin_corrected` with the
spec, physical bond counts — plus the B-centred count when the end element is the central one — and mean lengths in `#` headers), optionally a PNG plot (`--plot`, Agg
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
| `bond23` | `None` → `bond12` | inclusive B–C window; A = C ⇒ unordered counting of each triplet once (Step 5) |
| `bin_width` | `1.0`° | requested width; realized width is $180/\max(1,\lfloor 180/w+0.5\rfloor)$ |
| `collect_angles` | `False` | `bond_angle_distribution` only: keep the raw angle list (the one memory cost that grows with the angle count) |
| `max_angles` | `None` (unlimited) | refuse, before pairing, a spec whose exact angle count exceeds it; the app boundaries pass `APP_MAX_ANGLES` |
| `APP_MAX_ANGLES` | $5\times10^{7}$ | the shared app-boundary work budget (Python and JS constants must be equal) |
| `PAIR_CHUNK` | $2^{18}$ | bond pairs formed per streaming chunk (~25 MB) |
| `SEARCH_CHUNK` | $2^{20}$ | candidate pairs examined per search block and stencil offset (~100 MB) |
| `MAX_CELLS_PER_AXIS` | 64 | linked-cell resolution cap per lattice direction |
| `REACH_HEADROOM` | $10^{-9}$ | relative headroom on the layer reach against float rounding |
| `EDGE_SNAP_DEG` | $10^{-9}$° | an angle this close to a bin edge bins as exactly on it (Step 6); Python and JS constants must be equal |
| `LENGTH_BINS` | 40 | bond-length histogram bins per window (summary payload only) |
| App-boundary caps | $r_\mathrm{max}\le15$ Å, angles ≤ `APP_MAX_ANGLES`, `binWidth` ≥ 0.05° | Flask route and worker only; engine and CLI unrestricted |

### Parity: Python engine vs JavaScript port

[workers/triplets.js](../../web_app/frontend/src/workers/triplets.js) is a line-for-line port of
the Python engine, and the parity is pinned by golden fixtures rather than claimed:
[tests/generate_triplets_fixture.py](../../tests/generate_triplets_fixture.py) evaluates the
Python engine on four constructed configurations (and records the shared constants
`APP_MAX_ANGLES` and `EDGE_SNAP_DEG`, which the port must equal) and writes
[triplets_fixture.json](../../web_app/frontend/src/__tests__/fixtures/triplets_fixture.json),
which [workers/\_\_tests\_\_/triplets.test.js](../../web_app/frontend/src/workers/__tests__/triplets.test.js)
replays against the port:

| Fixture case | What it exercises |
|---|---|
| `ideal-perovskite` — undisplaced SrTiO₃ 3×3×3, fractions $(i+x)/3$; O–Ti–O, O–Sr–O, O–O–O and Sr–Ti–O (distinct windows) at 1°, 0.5° and 5° | symmetry angles exactly on bin edges (the edge snap); TiO₆ octahedra |
| `ideal-fcc` — undisplaced Cu 3×3×3 conventional cells; Cu–Cu–Cu at 1° and 3° | the case where every 60° angle used to change bin between engines |
| `random-triclinic` — 48 atoms, seeded RNG, lattice $[[6,0,0],[3,5,0],[1,1,7]]$; five specs: shared ends, different ends with distinct windows, B = A = C, and A = C with overlapping distinct windows (twice) | general triclinic geometry, both counting rules, 5°, 3° and 2° bins |
| `small-box-images` — 3 atoms in a 4 Å cube with windows reaching 3.5 Å | multiple periodic images of one atom as distinct neighbours |

Measured agreement, asserted per bin and per statistic:

- **Histogram counts, coordination, length-histogram counts: exact integer equality.** No
  tolerance — the two engines must produce the same integers.
- `density`, `sinCorrected`, bin centers: $10^{-9}$ (single-pass float math; covers
  summation-order noise).
- `meanAngle`, `stdAngle`, `meanLength`: $10^{-7}$.
- Sorted raw angles (head/tail samples): $10^{-5}$ — these go through `acos` twice
  (compute, then fixture rounding).

The work budget is pinned from both sides: `WorkBudgetTests` (exact count at the budget
accepted, one over refused, for every counting rule; streamed and blocked results identical to
unchunked ones; tracemalloc bounds on the streamed and refused paths), `TripletsBudgetApiTests`
in [tests/test_triplets_api.py](../../tests/test_triplets_api.py) (the route's 400), and the JS
`work budget` / worker-boundary suites.

The Python engine itself is pinned to a brute-force all-images reference over a $\pm2$ image
span on random triclinic configurations (`TriclinicBruteForceTests` in
[tests/test_triplets.py](../../tests/test_triplets.py)), plus the constructed invariants:
octahedron counting, bonds through the periodic wall, wrap invariance, zero-length exclusion,
self-image bonds, and the same-element distinct-window rule (continuity at touching windows,
one count per overlap triplet, disjoint shells). The backend route's caps and error paths
are covered by `TripletsApiTests` in [tests/test_backend_api.py](../../tests/test_backend_api.py).

**No residual divergence on ideal geometries.** Bond vectors, lengths and cosines are computed
in the same evaluation order in both engines (bitwise identical); only `acos` differs, by at
most 1 ulp, and the edge snap of Step 6 makes that difference invisible to the histogram. An
angle whose two engine values straddle an edge by more than float noise but less than
$10^{-9}$° away from it is measure-zero. Beyond the fixture, the Flask and worker paths were
compared on a CIF-built ideal GaNb₄Se₈ start configuration, the 5 K sample and its AVERAGE file
(21 specs): identical counts, coordination and length histograms. `IdealConfigurationTests`
pins the Python side (one bin per symmetry class at 1°, 0.5°, 5°; rigid-shift invariance; exact
numpy.histogram agreement off the edges).

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

**Compute** first turns the input boxes into a request with `tripletRequestFromInputs`
([workers/triplets.js](../../web_app/frontend/src/workers/triplets.js)): a cleared or non-numeric
box shows an error naming it ("A–B window minimum is empty — enter a number.") and nothing is
sent — it is never coerced to `0`. It then issues `requestPca('triplets', {end1, apex, end2,
r12Min, r12Max, r23Min?, r23Max?, binWidth})` — the B–C window included only when the split
switch is on (the engine then receives `bond23 = null` and reuses `bond12`). A dataset switch clears any previous result immediately
and bumps a `runEpoch` ref; a compute that was in flight for the old run compares its captured
epoch on resolve and can never land a stale payload on the new dataset.

### Step 3 — The result chips

The card's header names the triplet and the windows **the engine actually used** (the
resolved `bond12`/`bond23` of the payload) — both, labelled A–B and B–C, whenever they differ.
The chips are straight reads of the payload: central-atom count (`apexCount`), **Bonds** — the physical
bond count `uniqueBonds`, each bond once, with its mean length (`lengths12`, and `lengths23`
when not shared; a tooltip gives the B-centred count when the end element is the central one) —
the coordination summary — mean bonds per B $\sum_n n\,c_n / \sum_n c_n$, which counts a B–B
bond at both of its ends, as a coordination number should, plus the modal $n$ and its share —
and the angle count with mean ± std.

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
| A–B window | 2.0 – 3.0 Å | inclusive; string state, validated before sending (a cleared box is an error, not 0) |
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
