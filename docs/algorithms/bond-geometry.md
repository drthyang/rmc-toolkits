# Bond Geometry — algorithm reference

> Part of the [Algorithms and Math Reference](../ALGORITHMS.md). Every step is anchored to the source; if this document and the code disagree, **the code wins**.

Bond angles the RMCProfile `triplets` way: name an A–B–C triplet with **B the central atom**,
bound the two bond lengths, and histogram the angle at B over every triplet in the periodic
configuration — exactly, images included, nothing subsampled. The engine section covers the
neighbour search and the three normalizations; the page section covers what the Bond Geometry
tab adds on top of the payload (it computes no geometry of its own — only documented
presentation reductions of the payload and one closed-form reference line).

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
  - [Step 3 — The KPI rail](#step-3--the-kpi-rail)
  - [Step 4 — The angle plot and the random-bonds line](#step-4--the-angle-plot-and-the-random-bonds-line)
  - [Step 5 — The partial-g(r) window helper](#step-5--the-partial-gr-window-helper)
  - [Step 6 — The folded-cell bond view](#step-6--the-folded-cell-bond-view)
  - [Step 7 — Card states: empty, computing, stale, error](#step-7--card-states-empty-computing-stale-error)
  - [Step 8 — Layout and the colour system](#step-8--layout-and-the-colour-system)
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
| Flask | `/api/triplets` in [app.py](../../web_app/backend/app.py) | the uncached `bond_angle_summary_from_file`, memoized by `_TRIPLETS_CACHE` (a `_FileCache(16)`) keyed on (file signature, every parameter including the angle budget) |
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

From a file (`_read_configuration`, used by the Flask route and the CLI) the atoms come from the
shared `.rmc6f` grammar with `include_coords_only=True`: bond angles need only element and position,
so legacy coordinates-only lines count as well — the same atom set the browser worker takes from
`parseRmc6fAtoms()`. Lines with a non-finite coordinate are skipped and counted: the payload's
`parseWarning` names them in both runtimes (`null` when the file is clean; the `rmc-triplets` CLI
prints it on stderr, and the Model information card reports the same lines), and a file with no parseable atom is a `ValueError` (HTTP 400) *"no atoms could be
parsed — …"* naming what was found. Before 0.6.0 Flask read full-layout lines only, so a
coordinates-only file had "no atoms" there but angles in the browser.

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
therefore inside or outside a window from both of its ends alike. Before 0.6.0 the shift was added
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

$$\max(r_\mathrm{min} - \varepsilon, 0)^2 \le d^2 \le (r_\mathrm{max} + \varepsilon)^2
\quad\text{and}\quad d^2 > 0
\quad\text{and not}\quad(\text{same row} \wedge \mathbf m = \mathbf 0),
\qquad \varepsilon = \texttt{WINDOW\_TOL} = 10^{-9}\ \text{Å}.$$

Three deliberate choices:

- **Windows are inclusive at both ends** — read the bounds off the partial $g(r)$ and atoms at
  exactly the bound count. $\varepsilon$ makes that true for ideal geometries as well: a bound
  typed exactly at an ideal shell distance used to keep whichever bonds float rounding left
  inside (2400 of 3000 simple-cubic bonds at $r_\mathrm{max} = a$; a bond of exactly 1 Å dropped
  from a 0.5–1 Å window). $\varepsilon$ is four orders above the rounding noise of a length and
  far below any real distance difference (a bond $10^{-7}$ Å outside is still outside). The
  bond-length histograms clip admitted lengths into the window, so they still total `count`.
- **A zero-length pair is never a bond, even under `rmin = 0`.** Bitwise-coincident atoms have
  no direction; admitting the pair would put a $0/0$ NaN into every angle it joins
  (`CoincidentAtomTests`).
- **A center is never its own neighbour in the unshifted image** ($\mathbf m = \mathbf 0$), but
  its *other* images are genuine neighbours and stay — a B atom in a small box legitimately
  bonds to its own periodic copy (`test_self_image_bonds_of_the_central_element`).

Each surviving bond is stored as (center position in the selection, candidate row, image
$\mathbf m$, Cartesian vector, length) — the `_Bonds` container.

### Step 5 — Pairing bonds into angles

`_sort_by_center` sorts the bond lists by central atom (stable sort, so order is deterministic)
and `_Pairing` forms the per-center pairing, chunk by chunk:

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

> Before 0.6.0 the same-element/distinct-window case counted *ordered* (1→2, 2→3) assignments, so
> a triplet whose two bonds both lay in the overlap of the windows counted twice, and nudging a
> B–C bound by $10^{-4}$ Å off the A–B one doubled every count (223 651 → 447 317 Se–Nb–Se
> angles on the 5 K sample). The unordered rule counts each physical triplet once.

The angle is then

$$\theta = \frac{180°}{\pi}\arccos\!\Bigl(\operatorname{clip}\bigl(
\tfrac{\mathbf v_1\cdot\mathbf v_2}{\lVert\mathbf v_1\rVert\lVert\mathbf v_2\rVert},\,-1,\,1\bigr)\Bigr)$$

with the clip guarding the $\pm1$ boundary against rounding. Apart from the two $10^{-9}$
tolerances that make ideal configurations deterministic — `WINDOW_TOL` on the window bounds
(Step 4) and `EDGE_SNAP_DEG` on the bin edges (Step 6) — the calculation is exact geometry on
float64.

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
and keeps angles only for `collectAngles`. Before 0.6.0 both engines materialized every angle
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

For bonds pointing in independent uniformly-random directions this is flat at $1.0$, so
anything above 1 is real structure, and a peak near 180° (the octahedral *trans* angle) is no
longer suppressed by geometry.

**What the bin integral is — and is not.** By
$\cos(\theta_c - \tfrac{\Delta}{2}) - \cos(\theta_c + \tfrac{\Delta}{2}) = 2\sin\theta_c\sin\tfrac{\Delta}{2}$
(bin centre $\theta_c$, width $\Delta$ in radians), the reference is *exactly*
$\sin\theta_c\,\sin(\Delta/2)$:

$$S_k = \frac{N_k/N}{\sin\theta_c\,\sin(\Delta/2)} .$$

So `sin_corrected` is the familiar bin-centre $1/\sin\theta_c$ correction times the global
constant $1/\sin(\Delta/2)$ — the same shape, scaled so that random directions read exactly 1
(`SinCorrectionIdentityTests` pins the identity to $10^{-12}$). Neither form diverges: the bin
centres lie at $\Delta/2 \dots 180° - \Delta/2$, where $\sin\theta_c \ge \sin(\Delta/2) > 0$. Only a
per-angle weight $1/\sin\theta_i$, applied to each angle before binning, blows up at 0°/180°.
(Earlier versions of this page said the bin integral was needed to keep the end bins finite;
it is not — its benefit is the exact normalization.)

**Relation to RMCProfile's `triplets` output.** RMCProfile's TRIPLETS writes, per bin, `norm` —
the per-degree density, i.e. `density` here — and `norm/sin(theta)` $= D_k/\sin\theta_c$. With
$D_k = N_k/(N\,\Delta_\text{deg})$,

$$\texttt{norm/sin(theta)} \;=\; S_k\,\frac{\sin(\Delta/2)}{\Delta_\text{deg}}
\;\approx\; S_k\,\frac{\pi}{360}\quad(0.0087265\ \text{at } 1°),$$

a constant factor: **the same shape, not the same numbers** — `sin_corrected` is RMCProfile's
curve rescaled so that random is 1, not RMCProfile's normalization itself. The CLI writes the
exact factor for the realized bin width into its CSV header. Cross-check on the 5 K run folder,
which holds RMCProfile's own TRIPLETS output for the same configuration (`bonds_hist.pct`:
$r_\mathrm{max} = 3.5$ Å for every pair, 1000 bins of 0.18°): for Se–Nb–Se, Nb–Nb–Nb,
Se–Ga–Se, Nb–Se–Nb and Se–Nb–Nb the angle totals are identical (239 326, 47 078, 24 132,
95 731, 284 483); per-bin counts differ only by a few angles across a neighbouring edge
(RMCProfile bins in single precision; cumulative difference ≤ 4); `density` equals `norm`;
and `norm/sin(theta)` equals `sin_corrected` × $\sin(\Delta/2)/\Delta_\text{deg}$ to single
precision (`RmcProfileTripletsTests`, sample-backed, skipped without `data/`).

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

The engine itself is **unrestricted** — library and CLI callers can ask for anything. Step 5's
streaming bounds the memory the *angles* take; the bond lists themselves are stored in full and
grow as $\sim r_\mathrm{max}^3$, and the CLI has no $r_\mathrm{max}$ cap. The two app boundaries apply
identical request caps, each bounding a different cost — and in the browser an unbounded request
would freeze the shared PCA worker:

| Cap | Bounds | Flask `/api/triplets` | Worker `kind: 'triplets'` |
|---|---|---|---|
| $r_\mathrm{max} \le 15\,$Å (each window) | the neighbour search (bond count $\sim r_\mathrm{max}^3$) | 400 | thrown `Error` |
| exact angle count $\le$ `APP_MAX_ANGLES` $= 5\times10^7$ | the angles formed (count $\sim r_\mathrm{max}^6$) — and with them the pairing work, except for A = C with distinct windows (below) | 400 | thrown `Error` |
| `binWidth` $\ge 0.05°$ | the response size | 400 | thrown `Error` |
| `r12Min`/`r12Max` present, `r23Min`/`r23Max` both or neither — `null`, `''` and whitespace count as missing | — | 400 | thrown `Error` |

A missing bound is an error with the same text at both boundaries ("r12Min/r12Max are required
together; missing r12Min") — never a bound of 0. Before 0.6.0 the worker ran `Number()` on every
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

The budget counts angles, not candidate pairs. For A = C with **distinct** windows, Step 5 examines
every pair of the centre's combined bond list ($n^2$ candidates per centre in Python, $n(n-1)/2$ in
the JS loop, $n$ = bonds in either window) and keeps those with one bond in each window, so one
narrow window paired with a wide one costs more time than its angle count suggests. Every number
it produces is exact; only the run time is not bounded by `APP_MAX_ANGLES` there.

Request parameters are flat scalars (`end1`, `apex`, `end2`, `r12Min`, `r12Max`, `r23Min`,
`r23Max`, `binWidth`) so the identical request shape works as an HTTP query string and as a
worker message. The Flask side normalizes element case **before** the cache key
(`'se'` and `'Se'` share one entry) and resolves the `bond23` default before the call, so equal
windows hit one cache entry; the cache is `_TRIPLETS_CACHE` (a `_FileCache(16)` in `app.py`) keyed on
(file signature, every parameter), and a parse of a file that changed while it was read is never
cached. Errors map to 400 (bad or non-finite parameters, a file with no parseable atom), 403 (path
escapes the data root), 404 (no `.rmc6f`), 409 (the file kept changing while it was read).

### The `rmc-triplets` CLI

Console entry point installed by `pip install -e .`
([triplets_cli.py](../../rmc_toolkits/triplets_cli.py); module form
`python -m rmc_toolkits.triplets_cli`). Accepts an `.rmc6f` file or a run folder. A run folder
resolves to **the configuration the app analyses** (`find_run_configuration`, the same rule as the
backend's `_find_rmc6f` and the browser's `chooseStructureFile`): the `.rmc6f` whose stem matches
the run's own outputs (`<stem>-NN.log` first, then `<stem>_PDFpartials.csv`, `_FQ1.csv`, …), the
first sorted file only when nothing matches — so an input supercell `GaNb4Se8.rmc6f` beside the
refined `GaNb4Se8_5K.rmc6f` no longer wins by sorting first (`'.'` < `'_'`). The CLI prints the
path it chose. Outputs:

- the CSV (`--output`, default `triplets_<A-B-C>_<config>.csv`): columns `angle_deg, counts,
  density_per_deg, sin_corrected`, with the spec, physical bond counts — plus the B-centred
  count when the end element is the central one — mean lengths and the angle count in `#`
  headers;
- optionally a PNG (`--plot PATH`, Agg backend): **one** y-axis, the sin-corrected curve in
  its own units, with the per-degree density drawn dashed and **rescaled** so its peak meets the
  sin-corrected peak (legend "density, rescaled to the sin-corrected peak"). The dashed curve
  shows shape only — read density values from the CSV;
- optionally the raw angle list (`--dump-angles PATH`): one angle per line in degrees, 6
  decimals, in the engine's pairing order (unsorted). This is the one output whose memory grows
  with the angle count — the engine keeps the list only for it.

Nothing is overwritten without `--force`, and every destination is checked **before the angles
are computed** (`check_destinations`): `--output`, `--plot` and `--dump-angles` must be different
files (compared case-insensitively), none may be the configuration or a directory, each folder
must exist or be creatable, and the `--plot` extension must be a format this matplotlib can write
(no extension: PNG). A violation is one line on stderr and exit 1; `--force` relaxes only the
existing-file check. The outputs are then written through temporary files renamed into place
only after all of them succeeded, so a failed write leaves no partial set. Before 0.6.0 an
unsupported `--plot` format printed a traceback after the CSV was written, and `--dump-angles`
equal to `--output` silently replaced the histogram. `--version` prints the package version.

```bash
rmc-triplets data/5K_try1 --triplet Se Nb Se --bond12 2.2 2.9 --plot se_nb_se.png
rmc-triplets config.rmc6f --triplet O Ti O --bond12 1.7 2.3 --bond23 1.7 2.3 --bin-width 0.5
rmc-triplets data/5K_try1 --triplet Nb Nb Nb --bond12 2.6 3.4 --dump-angles nb_angles.txt
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
| `WINDOW_TOL` | $10^{-9}$ Å | a distance this close to a window bound counts as on it (inside); Python and JS constants must be equal |
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
- **Inclusive windows mean boundary atoms count** (to within `WINDOW_TOL` = $10^{-9}$ Å). Two
  runs whose $g(r)$ peak touches the bound can differ by exactly the boundary population —
  intentional, but worth knowing when comparing.
- **Ideal geometries are handled by two tolerances, both $10^{-9}$** — a distance that close to
  a window bound is on it, an angle that close to a bin edge is on it. They exist only to make
  float noise irrelevant; neither moves a real configuration's numbers.
- **App requests are refused above $5\times10^7$ angles** (`APP_MAX_ANGLES`, exact count, before
  any work on angles). The budget bounds the angles, not the candidate pairs of an A = C request
  with distinct windows (Step 8). Library and CLI calls are unrestricted: the angles stream in
  bounded memory unless the raw angle list is requested, but the stored bond lists grow as
  $\sim r_\mathrm{max}^3$.
- **`count` vs `uniqueBonds`.** `lengths.count` (and `bond12_count`) counts B-centred bond
  vectors — twice the physical bonds when the end element is the central one; `uniqueBonds`
  counts each bond once. The coordination numbers are per B and use the B-centred count.
- **The engine reports geometry, not chemistry.** A "bond" is a distance window and nothing
  else; there is no bond-valence, electronegativity, or connectivity analysis.

---

## Bond Geometry — the page

### What the page owns

[BondGeometryPage.jsx](../../web_app/frontend/src/components/BondGeometryPage.jsx) renders the
model-information card (same `ModelSummary` as the Dashboard, symmetry card omitted), one
controls bar — a `<form>`, so Enter in any field computes (triplet selects, windows,
`Distinct B–C` switch, bin width, **Compute** at the right) — and three cards: the angle
distribution as the **hero**, with the headline results in a KPI rail under its header, and
beside it the folded-cell bond view
([FoldedCellPanel.jsx](../../web_app/frontend/src/components/FoldedCellPanel.jsx)) above the
partial-$g(r)$ window helper. Everything numerical comes from the Step 7 payload; the page adds
selection, presentation (the reductions of Steps 3 and 4, each stated below), and the two helper
views. Every piece of chrome is the UI kit's (`src/ui`, see its README).

### Step 1 — Triplet seeding

The three selects are seeded **once per sites payload**, not per dataset key: on a dataset
switch the new key arrives while the previous run's sites are still in state, so keying on the
data itself waits for the real element list — and a Live Data refresh keeps the user's picks
when they still apply. The seed ranks elements by total atom count from the sites table: ends =
most abundant (in practice the anion), central = next most abundant — Se–Nb–Se for GaNb₄Se₈. An
existing valid selection is never overwritten. Each select carries its element's colour dot; B
keeps the ring; a ⇄ button swaps A and C (shown only when they differ).

### Step 2 — The compute request and the epoch guard

**Compute** (the button, or Enter in any field: the bar is a `<form noValidate>`, so the
browser's own number validation never blocks or rewrites a request) first turns the input boxes
into a request with `tripletRequestFromInputs`
([workers/triplets.js](../../web_app/frontend/src/workers/triplets.js)): a cleared or non-numeric
box is an error naming it ("A–B window minimum is empty — enter a number.") and nothing is sent —
it is never coerced to `0`. The page maps the message's leading label back to its input, marks it
`aria-invalid` (danger border) and focuses it; editing that field clears the mark and the
message. It then issues `requestPca('triplets', {end1, apex, end2, r12Min, r12Max, r23Min?,
r23Max?, binWidth})` — the B–C window included only when the split switch is on (the engine then
receives `bond23 = null` and reuses `bond12`). A dataset switch clears any previous result
immediately and bumps a `runEpoch` ref; a compute that was in flight for the old run compares its
captured epoch on resolve and can never land a stale payload on the new dataset.

The request a shown result was computed from is kept beside it. Whenever the current inputs
would send a different request — compared field by field as numbers, so `3.0` and `3.00` are the
same, and an input that does not parse counts as different — the result is **stale**: an amber
*inputs changed* chip appears in the hero header and the button reads **Update** with a dot
(both wait while a compute runs). The plot and the KPIs stay those of the shown result.

A **new configuration of the same run** — a Live Data save, picked up through the Flask
`dataEpoch` prop (App.jsx's `configEpoch`) or a browser-loaded run's changed `.rmc6f` text — bumps
the same epoch, keeps the triplet and the typed windows, reloads the element list, the Model
information card and the partials in place, and **drops** the computed distribution; the hero
header shows a *new configuration* chip and the empty-state prompt reads "New configuration —
Compute again.". The distribution is computed on demand, so it is never recomputed unasked, and a
result from the previous configuration never sits next to the new model.

### Step 3 — The KPI rail

The rail sits under the hero header and is there before Compute, its values reading "—", so
nothing below it moves when a result lands. Each tile is a straight read of the payload or one
stated reduction of it:

| Tile | Value | Sub line | Hover |
|---|---|---|---|
| **Angles** | `angleCount / apexCount`, 1 dp, "per B" | `angleCount` · realized `binWidth` · where it ran (`browser` for the worker, `server` for `/api/triplets`) | mean ± std of all angles (`meanAngle`, `stdAngle`) |
| **Coordination** | $\sum_n n\,c_n / \sum_n c_n$ from `coordination`, 2 dp, "per B" — counts a B–B bond at both of its ends, as a coordination number should | modal $n$ and its share of the central atoms · `apexCount` | — |
| **B–A bond** | `lengths12.meanLength`, 3 dp, Å | `uniqueBonds` (each physical bond once) · the resolved window `bond12` | the B-centred `count` when the end element is the central one |
| **B–C bond** | the same from `lengths23` — only when `sharedEnds` is false; two tiles on the same pair (A = C, distinct windows) are told apart as (A–B) / (B–C) | | |

Windows print at the precision they were given (2–4 decimals, so a B–C bound of 3.4001 does not
read as 3.40). The mean angle is kept off the headline on purpose: the mean of a multimodal
distribution (the Demo's Se–Ta–Se has four angle classes) is not a bond angle. It stays in the
hover text and in the CLI output.

### Step 4 — The angle plot and the random-bonds line

The hero plot shows `sinCorrected` or `density` (toggle; sin-corrected is the default) through
[InteractivePlot](../../web_app/frontend/src/components/InteractivePlot.jsx) with
`variant="fit"`: the SVG `viewBox` is taken from the rendered box (ResizeObserver, rounded to
whole pixels) instead of the fixed 8:5 aspect, one user unit = one CSS pixel. The angle axis uses
InteractivePlot's opt-in axis fields, which leave every other plot in the app untouched
(`interactivePlotAxes.test.jsx` pins Dashboard- and Auto StoG-shaped payloads to the markup the
component rendered before they existed):

- `xDomain: [0, 180]` — the whole angle range, **unpadded** (the automatic domain padded it to
  about −8…188°), which also bounds the wheel zoom;
- `xTicks` every 30° (0, 30, …, 180 — the angles a crystallographer reads: 60, 90, 120, 180),
  `xMinorStep: 10` unlabelled marks, `xGrid` vertical grid lines at the labelled ticks; after a
  zoom the nice 1-2-5 ticks take over;
- `yMin: 0` — the y axis starts at zero (no negative padding).

The data draw as a **step curve** — one flat step per bin across $\theta_c \pm w/2$, the realized
bin width from the payload — with a light same-colour area to $y=0$ (`curve: 'step'`,
`fill: true`), because the payload is a histogram, not a sampled function.

The dashed **random bonds** guide is what uniformly random bond directions give, i.e. the
isotropic reference of Step 6 drawn in the view's own units:

$$y_\text{random} = 1 \quad\text{(sin-corrected)},\qquad
y_{\text{random},k} = \frac{\cos\theta_k - \cos\theta_{k+1}}{2\,w}\ \ \text{deg}^{-1}\quad\text{(density)},$$

the second being the exact fraction $\tfrac12\int_{\theta_k}^{\theta_{k+1}}\sin\theta\,d\theta$
of random angles in bin $k$, per degree (it integrates to 1 over 0–180°, like `density`; test:
`BondGeometryLayout.test.jsx`). It is the only formula the page evaluates itself, and it is the
same bin integral the engine divides by for `sin_corrected`, so "above the line = more than
random" reads the same in both views. Axis labels: `angle at B, θ (°)` and
`sin-corrected (random = 1)` or `density (deg⁻¹)`.

### Step 5 — The partial-g(r) window helper

The helper plots the run's measured partial pair distribution from `PDFpartials.csv` (loaded via
the same plot-file path as the Dashboard; both runtimes) so the windows can be set against the
actual first shell:

- **Title**: the pair as element chips, B first — "(Ta)–(Se) partial g(r)" — or the whole
  triplet when there are two curves. The header's right side shows the **live window chip**
  (`2.00–3.00 Å`, following the inputs with the guides' debounce), or one chip per window, led by
  its bond-role dash, when the windows are split.
- **Curves**: the A–B partial always; a second curve whenever A–B and B–C are *different bond
  types* (pair labels looked up in either order — `Ta-Se` matches `Se-Ta`). This is independent
  of the window split: Ga–Ta–Se has two shells to bracket even with one shared window, and
  Se–Ta–Se has one shell even with two windows. The curves wear the bond-role colours
  (`BOND_COLORS.ab`, `.bc` — the first two plot colours).
- **Guides**: dashed verticals at the current bounds, following the inputs after a 400 ms
  debounce so the plot's view state does not reset per keystroke; a blank box draws no guide.
  With the split **off**, one neutral-grey pair covers both bonds. With it **on**, each window
  gets its own pair, labelled `A–B rmin/rmax` and `B–C rmin/rmax` **by role** (pair names would
  collide for same-element triplets) and coloured `BOND_COLORS.ab` / `.bc`, so when the bonds are
  different types each pair matches the curve it brackets (for a same-type triplet the B–C pair
  has the second colour and no curve of its own). The guides stay **out of the legend**
  (`legend: false`): the window chips name them, and the legend lists curves only.
- **Nothing is shaded**: the window is marked only by the guides. The in-app help (the A–B
  window and partial g(r) InfoBadges) states exactly these rules — the second curve follows the
  bond types, the switch only the guides — pinned by `BondGeometryPage.test.jsx`.
- **Crop**: the x-range is cut at $\max(6\,\text{Å},\ 2\times$ the furthest active
  $r_\mathrm{max})$ — beyond the first-shell region nothing informs a bond window.

The panel is display-only: the page never computes a $g(r)$. A run without `PDFpartials.csv` (or
without the pair in it) gets a slim card — the header row with a one-line note — and the folded
cell takes the freed height; everything else still works.

### Step 6 — The folded-cell bond view

[FoldedCellPanel.jsx](../../web_app/frontend/src/components/FoldedCellPanel.jsx) shows the same
folded unit cell as the Atomic Density page — every atom of the supercell folded into one cell
as a point cloud, colored by element, so the spread around a site is the *measured* thermal
cloud rather than a fitted ellipsoid — with the analysis' detected bonds drawn over it:

- **Title**: the bond as element chips — "(Ta)–(Se) bonds", or the whole triplet when B–C is a
  bond of its own — for the computed triplet (the current picks before Compute).
- **Cloud**: one `THREE.Points` per element; above 120 000 atoms the cloud is strided down to
  ~120 k points (display only — the engine always sees every atom).
- **Bonds**: for each computed window, average-site pairs whose distance falls inside it,
  periodic images included ($\mathbf m \in \{-1,0,1\}^3$, same-element in-cell pairs
  deduplicated) — so a stick may reach an image just outside the box, which is the real
  coordination. Drawn as thin transparent `LineSegments` (opacity 0.45) in the bond-role
  colours — A–B `BOND_COLORS.ab`, B–C `BOND_COLORS.bc` when it is a bond of its own
  (`!sharedEnds`) — so a full network reads as a framework without hiding the cloud.
- **Legend**: a pill inside the canvas, bottom-left: every element's colour, the triplet's
  elements in bold and the rest muted, then one bond swatch per drawn window
  (`Ta–Se 2.00–3.00 Å`); before Compute it ends with *Compute to draw bonds*.
- The a/b/c gizmo, reset view, and 1×/3× PNG export follow the other 3D panels.

Note the sticks connect **average site positions** (the folded reference sites), while the cloud
shows instantaneous atoms: a stick is the average bond, not any single configuration's bond.

### Step 7 — Card states: empty, computing, stale, error

The three cards keep their skeleton — header, KPI rail, plot toolbar, plot — before and after
Compute, so nothing moves when a result lands:

- **Empty** (a run is open, nothing computed): the hero shows the angle axis of Step 4 dimmed
  and inert (the *ghost*: same domain, ticks and the random-bonds line, at the typed bin width in
  the density view), with a prompt card centred over it: the triplet chips and window chip, one
  line ("Pick a triplet, then Compute."), a **Compute** button, and chips for the central-atom
  count, the supercell and where it will run (their hovers say the rest: every B atom in the box,
  periodic images included, exact).
- **Computing**: the first compute sweeps a shimmer over the ghost; a recompute dims the shown
  plot (inert); both show a centred *Computing A–B–C…* badge, and the button a spinner
  (`aria-busy`). Motion is off under `prefers-reduced-motion`.
- **Stale**: see Step 2.
- **Error**: the prompt's line becomes the message (`role="alert"`) and the named field is
  marked (Step 2). A run that cannot be read at all keeps the page-level banner.
- **No run**: the page-level prompt "Open a run folder with an `.rmc6f` file." and the empty
  cards.

### Step 8 — Layout and the colour system

**Layout** ([BondGeometryPage.css](../../web_app/frontend/src/components/BondGeometryPage.css)):
above 1100 px the grid is two columns, `minmax(0, 7fr) minmax(0, 5fr)`, with rows
`minmax(0, 1.25fr) minmax(0, 1fr)` and areas `"hero cell" "hero pdf"`: the hero spans both rows
(about 904 × 615 px at 1600 × 900), the folded cell and the partial $g(r)$ stack on the right.
The grid takes the height left under the model card and the controls (`flex: 1 1 0`, floor
32 rem, cap 60 rem) rather than a fixed `100vh − k` clamp, so it absorbs the model card wrapping
at 1440 px. The controls bar keeps **one row down to 1440 px with the split window on**: the
window labels are the bare bond roles (*A–B*, *B–C*; the group and the inputs keep "window" in
their accessible names) and the kit's *Enter* key hint hides below 1500 px (the button's title
names the key). Measured on the Demo run (page root `scrollHeight − clientHeight`): 0 at
1440 × 900, 1600 × 900 and 1920 × 1080, with the split off and on, for Se–Ta–Se and Ga–Ta–Se,
before and after Compute. Narrower, the bar wraps and the grid absorbs the extra row until its
floor; below the floor the page scrolls rather than squeezing the plots — at 1280 × 800 the model
card takes three rows, and the page scrolled there before too. The card
headers wrap their actions under the title on narrow cards. Without a partials file the rows
become `minmax(0, 1fr) auto`. At ≤ 1100 px everything stacks: hero
`clamp(24rem, 62vh, 36rem)`, folded cell 24 rem, partial $g(r)$ 18 rem.

**Colour system** — one meaning per colour, on all three cards:

| Colour | Marks | Where |
|---|---|---|
| element colours (`buildElementColors`) | atoms | select dots, element chips (dot, tint and the central atom's ring), 3D cloud and legend |
| `BOND_COLORS.ab` (= `PLOT_PALETTE[0]`) | the A–B bond | chip bond dashes, the A–B label bar, 3D sticks, split guides and window chip, the first partial curve, the bond KPI dash |
| `BOND_COLORS.bc` (= `PLOT_PALETTE[1]`) | the B–C bond, when it is its own | the same places, for B–C |
| `GUIDE_STROKE` (neutral grey) | references | the random-bonds line, a lone (unsplit) window's guides |

Text is never element-coloured (contrast in both themes): element colours only fill dots, tints,
rings and 3D objects.

### Parameters and defaults

| Control | Default | Notes |
|---|---|---|
| Triplet A, B, C | seeded per sites payload | ends = most abundant element, central = next |
| A–B window | 2.00 – 3.00 Å | inclusive; string state, validated before sending (a cleared box is an error, not 0) |
| Distinct B–C | off | off ⇒ B–C reuses the A–B window and one guide pair |
| B–C window | 2.00 – 3.00 Å | only sent when the split is on |
| Bin width | 1.0° | realized width comes back in the payload |
| Angle view | sin-corrected | toggle to per-degree density |
| Guide debounce | 400 ms | `useDebounced` on all four window inputs (guides and window chips) |
| Cloud stride cap | 120 000 points | `MAX_CLOUD_POINTS`, display only |

### Caveats

- **The page computes no geometry.** Every number on it is the engine payload or a reduction
  stated in Step 3 (angles per B, mean coordination, modal share); the one formula it evaluates
  itself is the random-bonds reference of Step 4, the closed-form isotropic bin fraction. The
  helper and the folded cell are presentation over `PDFpartials.csv` and the sites table
  respectively.
- **The bond-length histograms are in the payload but not plotted.** An earlier layout gave them
  a panel; it duplicated the first-shell peak the partial $g(r)$ already shows, clipped to the
  window. The counts and mean lengths survive in the bond KPI tiles.
- **The folded-cell sticks are average-structure bonds.** They match sites within the window at
  their *average* positions; a strongly displaced site whose average distance falls outside the
  window shows no stick even though many instantaneous bonds were counted (and vice versa).
- **The helper needs `PDFpartials.csv`.** Runs without it lose the guides' context but nothing
  else; the windows still apply exactly as typed.
