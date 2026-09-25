# PCA Ellipsoid — algorithm reference

> Part of the [Algorithms and Math Reference](../ALGORITHMS.md). Every step is anchored to the source; if this document and the code disagree, **the code wins**.

The second moment of each site's displacement cloud: anisotropic displacement parameters, the separable 3-D Gaussian KDE and its isosurface, and how the principal axes relate to the crystallographic axes.

## Contents

- [PCA Ellipsoid page — displacement clouds, ADP tensors, and the separable 3D KDE](#pca-ellipsoid-page--displacement-clouds-adp-tensors-and-the-separable-3d-kde)
  - [What this page shows](#what-this-page-shows)
  - [Step 1 — Parse `.rmc6f` into per-site displacement clouds](#step-1--parse-rmc6f-into-per-site-displacement-clouds)
  - [Step 2 — Select the cloud (single site, element-pooled, or all)](#step-2--select-the-cloud-single-site-element-pooled-or-all)
  - [Step 3 — Covariance and eigen-decomposition](#step-3--covariance-and-eigen-decomposition)
  - [Step 4 — ADP readouts](#step-4--adp-readouts)
  - [Step 5 — Probability-ellipsoid scaling (the $\chi^2$ quantile)](#step-5--probability-ellipsoid-scaling-the-chi2-quantile)
  - [Step 6 — Bandwidth and the sampling box](#step-6--bandwidth-and-the-sampling-box)
  - [Step 7 — The separable 3D Gaussian KDE](#step-7--the-separable-3d-gaussian-kde)
  - [Step 8 — Isosurface level selection](#step-8--isosurface-level-selection)
  - [Step 9 — Marching cubes](#step-9--marching-cubes)
  - [Step 10 — Wall projections and their contours](#step-10--wall-projections-and-their-contours)
  - [Step 11 — Density painted on the ellipsoid ("Shell")](#step-11--density-painted-on-the-ellipsoid-shell)
  - [Step 12 — The non-Gaussianity readout](#step-12--the-non-gaussianity-readout)
  - [Step 13 — PCA frame ↔ crystallographic frame (pointer)](#step-13--pca-frame--crystallographic-frame-pointer)
  - [Step 14 — The Site-ellipsoids (unit-cell) panel and site selection](#step-14--the-site-ellipsoids-unit-cell-panel-and-site-selection)
  - [Parameters and defaults](#parameters-and-defaults)
  - [Caveats / what this is not](#caveats--what-this-is-not)
- [Principal axes in the crystallographic frame](#principal-axes-in-the-crystallographic-frame)
  - [What this shows](#what-this-shows)
  - [Steps 1–2 — The eigenframe this section starts from (owned by the PCA Ellipsoid section)](#steps-12--the-eigenframe-this-section-starts-from-owned-by-the-pca-ellipsoid-section)
  - [Step 3 — Sign and handedness canonicalisation, and what a direction means](#step-3--sign-and-handedness-canonicalisation-and-what-a-direction-means)
  - [Step 4 — Per-axis excess kurtosis $\kappa_i$ (pointer)](#step-4--per-axis-excess-kurtosis-kappa_i-pointer)
  - [Step 5 — Unit-cell vectors in the shared Cartesian basis](#step-5--unit-cell-vectors-in-the-shared-cartesian-basis)
  - [Step 6 — Frame transforms between fractional and PCA coordinates](#step-6--frame-transforms-between-fractional-and-pca-coordinates)
  - [Step 7 — Orientation of each principal axis against $a$, $b$, $c$](#step-7--orientation-of-each-principal-axis-against-a-b-c)
  - [Step 8 — What the statistics panel prints: the Crystal orientation table](#step-8--what-the-statistics-panel-prints-the-crystal-orientation-table)
  - [Step 9 — What the UI shows in 3D: the PC ↔ Crystal frame switch](#step-9--what-the-ui-shows-in-3d-the-pc--crystal-frame-switch)
  - [Step 10 — What is displayed, and what is computed but not](#step-10--what-is-displayed-and-what-is-computed-but-not)
  - [Step 11 — A different displacement measure: `dispA`](#step-11--a-different-displacement-measure-dispa)
  - [Parameters and defaults](#parameters-and-defaults-1)
  - [Caveats — what this is not](#caveats--what-this-is-not-1)

---

## PCA Ellipsoid page — displacement clouds, ADP tensors, and the separable 3D KDE

### What this page shows

RMCProfile writes, for every atom in the box, the *reference site* it belongs to and the integer
*cell index* of the supercell copy it sits in. Subtracting a site's average position from each of its
copies leaves one **displacement cloud** per crystallographic site: a set of $N$ Cartesian offsets
(Å) about the average structure. The page turns that cloud into three things:

1. an **anisotropic displacement tensor** $\mathbf{U}$ (the cloud covariance) and its
   eigen-decomposition — the thermal ellipsoid drawn at a chosen enclosed probability;
2. a **3D Gaussian kernel-density estimate (KDE)** of the same cloud, shown as an isosurface, as
   density painted on the ellipsoid surface, and as three 2D "shadow" projections on the walls of a
   box (a layout after Maksim Eremenko's `PCA_KDE` utilities, credited in the code and in the
   viewport legend);
3. a **non-Gaussianity readout** (excess kurtosis) quantifying how far the cloud departs from the
   harmonic (Gaussian) model the ellipsoid assumes.

The whole analysis is *per site* (or, in the engines, per element with sites pooled). Nothing on this
page fits a model to data — it is descriptive statistics of the RMC configuration plus a smoothing
estimator.

#### Two engines, one algorithm

| | Flask (server) mode | Static (browser) mode |
|---|---|---|
| Engine | [`rmc_toolkits/pca_kde.py`](../../rmc_toolkits/pca_kde.py) (source of truth) | [`web_app/frontend/src/workers/pcaKde.js`](../../web_app/frontend/src/workers/pcaKde.js) |
| Entry points | `/api/pca/sites`, `/api/pca/kde` in [`web_app/backend/app.py`](../../web_app/backend/app.py) | `pcaKdeWorker.js` messages `{kind:'sites'}` / `{kind:'kde'}` |
| Eigensolver | `numpy.linalg.eigh` (LAPACK) | 3×3 cyclic **Jacobi** (`jacobiEigenSymmetric`) |
| $\chi^2$ quantile | `scipy.stats.chi2.ppf` | exact inverse of the closed-form CDF (`chiSquare3Quantile`, Newton) |
| Contraction | BLAS `dgemm` via `numpy` | hand-rolled triple loop with a reused buffer |
| Site reconstruction for old files | **not implemented** | fold-and-cluster fallback (`sitesByClustering`) |

A third request kind exists on both sides — `/api/pca/orientation` and the worker's
`{kind:'orientation'}`, routed to `rmc_toolkits/orientation.py` / `workers/orientation.js`. It is
**not** part of this page: it feeds the separate **Displacement Directions** page, documented in the
sibling [Displacement Directions](displacement-directions.md) reference. Nothing on the PCA Ellipsoid
page issues it. It appears here only because the same worker, the same parse cache and the same
routing function serve both pages.

#### Where the code lives

The page was split into a shared hook plus three view modules when the Displacement Directions page
became its own page; the pieces below are the current owners of each concern.

| Module | Owns |
|---|---|
| [`workers/pcaKde.js`](../../web_app/frontend/src/workers/pcaKde.js) | the browser engine: parsing, clouds, covariance, eigensolve, KDE, projections |
| [`workers/pcaKdeWorker.js`](../../web_app/frontend/src/workers/pcaKdeWorker.js) | the worker message loop: content-keyed parse cache, `sites` / `kde` / `orientation` dispatch |
| [`useSiteCloud.js`](../../web_app/frontend/src/useSiteCloud.js) | reading the local `.rmc6f` text, worker-vs-Flask routing (`requestPca`), the per-site ellipsoid table, the selected site, and `unitCell` |
| [`components/PcaKdePage.jsx`](../../web_app/frontend/src/components/PcaKdePage.jsx) | the KDE request, the main three.js viewport (isosurface, shell, walls, cameras), the statistics panel, the controls |
| [`components/SiteStructurePanel.jsx`](../../web_app/frontend/src/components/SiteStructurePanel.jsx) | the unit-cell Site-ellipsoids picker (Step 14), shared with the Displacement Directions page |
| [`components/sceneAxes.js`](../../web_app/frontend/src/components/sceneAxes.js) | the PC and a/b/c colour palettes and the `buildAxisTriad` / `buildCrystalAxes` rod builders, shared by every panel |
| [`pcaCrystalFrame.js`](../../web_app/frontend/src/pcaCrystalFrame.js) | the crystal-frame algebra (Step 13): `unitCellVectors` (via `useSiteCloud.js`) and `crystalOrientationRows` → `principalAxisOrientation` → `crystalPcaTransforms` (the *Crystal orientation* table in `PcaKdePage.jsx`) |

Which engine runs is **not** decided by the runtime mode alone. The decision lives in
[`useSiteCloud.js`](../../web_app/frontend/src/useSiteCloud.js) → `requestPca(kind, params)`, not in the
page: when the hook holds the `.rmc6f` text of a locally-loaded run (the demo or a picked folder) it
posts to the Web Worker in *either* runtime mode; otherwise it issues an HTTP GET against
`{ sites: '/api/pca/sites', orientation: '/api/pca/orientation' }[kind] ?? '/api/pca/kde'` for the
typed backend directory. `PcaKdePage.jsx` destructures `requestPca` (and `sites`, `selectedRef`,
`selectedEllipsoid`, `unitCell`, `datasetKey`, …) out of the hook and only issues the `kde` request
itself. The stats panel prints `· browser` or `· server` so you can always tell which produced the
numbers on screen (the JS result carries `browserPcaKde: true`).

#### Notation and units

| Symbol | Meaning |
|---|---|
| $\mathbf{f}_n$ | stored fractional coordinate of atom $n$ — a fraction of the **supercell**, not of one unit cell |
| $\mathbf{c}_n$ | integer cell indices of atom $n$ (which box copy it belongs to) |
| $\mathbf{N}=(N_1,N_2,N_3)$ | supercell repeat counts |
| $\mathsf{L}$ | supercell lattice matrix, **rows** = the three supercell edge vectors in Å |
| $\mathsf{A}$ | unit-cell matrix, rows $\mathbf{a},\mathbf{b},\mathbf{c}$ in Å; $\mathsf{A}_{i\cdot}=\mathsf{L}_{i\cdot}/N_i$ |
| $\delta\mathbf{f}_n$ | folded, mean-subtracted fractional offset — a fractional *displacement* (dimensionless, supercell fraction) |
| $\mathbf{u}_n$ | Cartesian displacement of atom $n$, Å |
| $\mathbf{U}$ (code: `covariance`) | 3×3 displacement covariance, Å² |
| $\lambda_1\ge\lambda_2\ge\lambda_3$ | eigenvalues of $\mathbf{U}$, Å² (code: `eigenvalues`) |
| $\mathsf{P}$ | principal axes as **rows**, orthonormal, right-handed (code: `axes`) |
| $\sigma_a=\sqrt{\lambda_a}$ | principal RMS amplitudes, Å (code: `rms`) |
| $k(p)$ | ellipsoid scale factor at enclosed probability $p$ (code: `probabilityScale`) |
| $f$ | KDE bandwidth factor (Scott/Silverman), dimensionless (code: `factor`) |
| $h_a=f\sigma_a$ | per-axis kernel bandwidth, Å (code: `bandwidth`) |
| $G$ | grid points per axis (code: `grid`) |
| $\tau(p)$ | isosurface density threshold enclosing mass fraction $p$, Å⁻³ |

---

### Step 1 — Parse `.rmc6f` into per-site displacement clouds

**Inputs.** One `.rmc6f` file: a `Supercell dimensions` line, a `Lattice vectors (Ang):` line
followed by three rows, and an `Atoms:` section.

**Header.** `rmc_toolkits/parsers.py` → `read_cell_vectors()` takes the last three tokens of the
`Supercell` line as $\mathbf{N}$ and the three lines after `Lattice` as the rows of $\mathsf{L}$
(Å). The JS twin is `pcaKde.js` → `readCellVectors()`; both scan the whole file and keep the *last*
match. Unit-cell vectors are $\mathsf{A}_{i\cdot} = \mathsf{L}_{i\cdot}/N_i$
(`SiteDisplacements.unit_vectors`; JS `pcaCrystalFrame.js` → `unitCellVectors()`).

**Atom lines.** `parsers.py` → `iter_rmc6f_atoms()` indexes from the *end* of the line, so any number
of label columns between element and coordinates is tolerated:

```
… x y z ref cellx celly cellz     ← last 4 fields are ref + cell indices,
                                     the 3 before them are the coordinates
```

Lines with fewer than 9 whitespace-separated fields are **silently skipped** by the Python parser.
The browser parser `web_app/frontend/src/rmc6f.js` → `parseAtomLine()` additionally accepts the old
5–6-field coordinates-only form, returning `referenceNumber = null` and `cellIndices = null`.
**Non-finite coordinates.** Python's `iter_rmc6f_atoms` accepts `NaN`/`inf` tokens, and one such
coordinate used to poison its site's mean and fail the *batched* eigensolve for every site
(`LinAlgError: Eigenvalues did not converge`, an HTTP 500 for the whole page).
`load_site_displacements()` now raises a `ValueError` naming the atom and its reference site
(`/api/pca/sites` returns it as a 400), and `site_ellipsoids()` refuses a `SiteDisplacements` whose
cloud is non-finite, naming the sites. The browser parser (`rmc6f.js` → `parseAtomLine()`) instead
drops such a line silently, so the static app still shows that site with one copy fewer; making the
two *parsers* agree is outside this engine.

Element symbols are normalised identically in both engines — Python's `str.capitalize()` in
`iter_rmc6f_atoms`, and `capitalizeElement()` in `siteDisplacementsFromRmc6f()` (`SE` → `Se`) — so
a site and its pooling carry the same label in either runtime.

**The displacement convention** (`pca_kde.py` → `load_site_displacements()`; JS
`sitesByReferenceNumber()` + `buildSite()`):

$$\mathbf{o}_n \;=\; \mathbf{f}_n - \frac{\mathbf{c}_n}{\mathbf{N}},\qquad
\mathbf{m}_s=\frac{1}{2\pi}\arg\!\sum_{n\in s}e^{2\pi i\,\mathbf{o}_n},\qquad
\mathbf{o}_n \leftarrow \mathbf{o}_n - \mathrm{round}\big(\mathbf{o}_n-\mathbf{m}_{s(n)}\big)$$

componentwise. The subtraction removes the origin of the atom's own box copy, leaving
$\mathbf{o}_n=(\mathbf{s}+\boldsymbol{\delta}_n)/\mathbf{N}$ — the site's position $\mathbf{s}$ inside
the unit cell plus the atom's displacement, as a supercell fraction — known only **modulo one
supercell period**, because the stored coordinate wraps at the box edge (an atom that drifted across
the outer edge is stored as $\approx 0.999$ when it belongs at $\approx -0.001$). The unwrap is taken
about each site's own **circular mean** $\mathbf{m}_s$ per axis (period 1; `_circular_site_centres`
in Python, `circularMean` in JS), which no wrap of an individual offset can move: every offset is
shifted by the whole number of periods that brings it within half a period of its site's centre. JS
rounds half to even (`roundHalfEven`) so a tie folds exactly as `np.round` does.

**Why not a fold about zero.** Up to 1.0 both engines applied $\mathbf{o}\leftarrow\mathbf{o}-\mathrm{round}(\mathbf{o})$,
a half-box minimum image about the copy's *origin*. That is only safe while
$(s_i+\delta_i)/N_i<\tfrac12$ — always true for $N_i\ge3$, but not for a site at $s_i\approx\tfrac12$ in
a one-cell-thick box ($N_i=1$, a common special position) or at $s_i\approx1$ with $N_i=2$. There the
copies on either side of the cut landed a whole box edge apart, the site mean fell between the two
lumps, $U$ came out $10^3$–$10^4$ times too large and the site marker was misplaced — silently, in both
engines, and in the Displacement Directions page that reads the same clouds. About the site centre
only a genuine displacement of half a *supercell* could fold. On the repository's GaNb₄Se₈ and GTS
runs the two rules give bit-identical results.

`tests/test_pca_kde.py::test_supercell_boundary_wrap_folds` and the JS `pcaKde.test.js` "folds an
atom that drifted across the supercell boundary" pin the edge wrap (four atoms straddling the edge
of a $1\times1\times1$, 8 Å box end up within 0.1 Å of each other); `tests/test_pca_regressions.py`
`SiteCentredFoldTests` and `pcaKdeRegressions.test.js` pin the special positions ($N_i=1$ at
$s=\tfrac12$, $N_i=2$ at $s=0.99$, and a $1\times8\times8$ box) with nonzero spread.

Each site $s$ is then centred on its own mean and mapped to Cartesian Å through the **full supercell
lattice**:

$$\bar{\mathbf{o}}_s=\frac{1}{n_s}\sum_{n\in s}\mathbf{o}_n,\qquad
\delta\mathbf{f}_n=\mathbf{o}_n-\bar{\mathbf{o}}_{s(n)},\qquad
\mathbf{u}_n=\delta\mathbf{f}_n\,\mathsf{L}\quad\mathrm{i.e.}\quad u_{nk}=\sum_i \delta f_{ni}\,\mathsf{L}_{ik}.$$

Using $\mathsf{L}$ (not $\mathsf{A}$) is correct precisely because $\delta\mathbf{f}$ is a *supercell*
fraction. Python does the per-site means with three `np.bincount` reductions; JS accumulates per
cloud in `buildSite`. Both are exact, unweighted arithmetic means, so `displacements` is centred to
round-off (asserted at `atol=1e-9` in Python, `<1e-9` in JS).

The site's position inside the average unit cell is reported as
$\;\mathrm{mod}(\bar{\mathbf{o}}_s\odot\mathbf{N},\,1)$ (`site_fractional`, serialised as
`siteFractional`), which is what `SiteStructurePanel.jsx` uses to place each marker (Step 14).

**Outputs.** `SiteDisplacements` (Python dataclass) / `{referenceNumbers, sites, latticeVectors,
supercell, reconstructed}` (JS): per-site element label, copy count, centred Cartesian cloud (Å),
fractional site position, and each row's **own** species (`atom_elements` / `atomElements`).

**Mixed-occupancy sites.** A reference number can be carried by atoms of more than one species — a
solid solution, or the result of RMCProfile swap moves. Both engines apply one rule
(`_site_compositions` / `siteComposition`): the site's `element` label is its **majority species**,
ties going to the alphabetically first (code-point order, never locale); `elementCounts` holds the
full composition (species → copies, in name order) and `mixed` is true when it has more than one
entry. Both appear in every `sites` row and in a per-site KDE payload. The site's cloud is still all
of its atoms about their common mean position, so $\mathbf U$, $\kappa$ and the KDE describe the
site as a whole. The payload's `elements` list (Python `SiteDisplacements.species`, JS
`siteSpecies()`) is every species present, minority species of mixed sites included. (Before 1.0
Python labelled such a site by the *last* atom in file order and JS by the *first*, so the same file
showed `#1 In` in one runtime and `#1 Ga` in the other, and the minority species vanished from the
element list.) The page shows a mixed site as its composition — `Ga0.75In0.25` in the picker and the
viewport heading — with a `mixed · Ga 48 · In 16` tag in the Summary column.

**Caching.** Python: `cached_site_displacements()` is `functools.lru_cache(maxsize=8)` keyed on
`(path, st_mtime)`. JS: `pcaKdeWorker.js` → `parseCached()` keys on a **cheap content signature** —
FNV-1a over every 64th character plus the total length, concatenated with the cluster threshold. That
is a heuristic, not a hash of the full text; two different files of identical length whose every-64th
characters agree would collide and reuse the wrong parse. The key is deliberately content-based
rather than path-based so two runs sharing a filename do not alias. The worker holds **one** cache
slot (`cache = { key, parsed }`), so alternating between two files re-parses on every switch.

The worker itself is now **one instance for the whole app**: `useSiteCloud.js` keeps it in a
module-level `sharedWorker` created lazily by `getWorker()` and **never terminated**, with requests
correlated by a monotonically increasing `id`. This is a real behavioural change from the previous
in-page worker, which was created in a `useEffect` on mount and `terminate()`d on unmount: the parse
cache now survives navigation and is shared between the PCA Ellipsoid and Displacement Directions
pages, so switching pages (or tabs within a page) re-uses the parsed clouds instead of re-parsing
the configuration. The cost is that the parsed model and one worker thread stay resident for the
lifetime of the tab.

#### Step 1b — Site reconstruction for coordinates-only files (browser only)

Unless **more** than half the atom lines carry both a reference number and cell indices — the test is
a strict `tagged.length > atoms.length / 2`, so an exact half also falls through — `pcaKde.js` →
`siteDisplacementsFromRmc6f()` takes the fallback path `sitesByClustering()`:

1. fold every atom into one unit cell: $u_{ni} = \mathrm{mod}(x_{ni}N_i,1)$;
2. per element, **periodic single-linkage cluster** by minimum-image Cartesian distance under the
   full unit-cell metric (`clusterPeriodic`), with a union–find over a uniform bin grid. Bins are
   sized by each axis's *perpendicular* width (cell volume ÷ opposite-face area), so a fractional
   step of $1/\mathrm{bins}$ always spans at least the threshold even for an oblique cell; the
   minimum-image distance searches all 27 neighbouring images because an oblique cell's nearest copy
   can be diagonal;
3. unwrap each cluster about its **circular mean** per axis
   ($\bar{u}=\arg(\sum e^{2\pi i u})/2\pi$) rather than about an arbitrary member, so a wide cluster
   (e.g. an orientationally disordered rotor shell) is not split;
4. sort clusters by (element, folded centroid $x$, $y$, $z$) and assign synthetic reference numbers
   $1,2,\dots$;
5. report `copiesPerCell = round(N_1N_2N_3)` alongside each site's member count, and set
   `reconstructed: true`.

Default threshold `DEFAULT_CLUSTER_THRESHOLD = 1.5` Å (UI slider 0.4–2.5 Å, step 0.1, shown only for
such files). A count equal to `copiesPerCell` is labelled *clean site*; larger is *merged /
disordered*; smaller is *fragmented — raise the distance*. **The Flask backend has no equivalent**:
`iter_rmc6f_atoms` drops sub-9-field lines and unconditionally reads a reference number plus three
cell indices from the last four fields, so a coordinates-only file yields **no atoms at all** and
`load_site_displacements()` raises *"… does not contain any atoms"* in server mode. The
`reconstructed` flag, `copiesPerCell`, the Cluster-distance slider, the `({count}/{copiesPerCell})`
suffix on the site picker and the clean/merged/fragmented tag therefore exist only in the browser
path.

---

### Step 2 — Select the cloud (single site, element-pooled, or all)

`pca_kde.py` → `displacement_cloud()` / JS `sitePcaKde()`:

- `reference_number=` → the rows of that one site;
- `element=` (not `""`/`"all"`) → **every atom of that species**, matched case-insensitively
  against each atom's own element (`atom_elements` / `atomElements`), not against the site label —
  so a mixed site contributes only its matching atoms, still centred on the site's all-species
  mean. JS: `displacementCloud()`, which `sitePcaKde()` calls;
- neither → every atom in the configuration.

Pooling across sites is meaningful *because each site was already centred on its own average
position* (Step 1): the union describes how an element moves, independent of where its sites sit in
the cell. It is not meaningful as a single "ellipsoid" if the pooled sites have different
orientations — the pooled covariance is then a symmetry-averaged tensor, not any one site's ADP.
The current UI only ever requests a single `referenceNumber`; element pooling is reachable through
the API/engine only.

**Subsampling.** `pca_kde_volume()` caps the cloud fed to the KDE at
`MAX_PCA_FIT_POINTS = 20000` points (same constant in both engines). Python draws
`np.random.default_rng(rng_seed).choice(n, limit, replace=False)` then sorts the indices; JS runs a
partial Fisher–Yates with a `mulberry32(seed)` PRNG then sorts. Both are deterministic with default
seed 0, but they draw **different subsets** — so above 20 000 points the two engines agree only
statistically, not to round-off. A single site in an $n\times n\times n$ box has $n^3$ copies — 1000
for a 10³ box, but 21 952 for 28³ and 27 000 for 30³ (a normal ~54k-atom box for a two-atom cell) —
so the cap binds for **per-site clouds in boxes of 28 or more cells per edge** as well as for
element-pooled clouds. There the KDE volume, its mass levels, its axes and eigenvalues come from
20 000 of the copies, and the Flask and static apps use different ones (a 27 000-copy Gaussian site:
eigenvalues differ by up to 0.6%; the audit measured up to 6% of $v_\mathrm{max}$ in the density).
The result reports both `count` (total) and `fitCount` (used), and the panel's metadata line prints
`fit {fitCount}/{count}`.

Note `site_ellipsoids()` (the ADP table) **never** subsamples — it always uses every atom of every
site.

---

### Step 3 — Covariance and eigen-decomposition

**Covariance.** For a site with $n$ copies and already-centred displacements,

$$U_{ab}=\frac{1}{n-1}\sum_{n} u_{na}u_{nb},\qquad a,b\in\{x,y,z\},\quad[\text{Å}^2].$$

This is the unbiased ($\mathrm{ddof}=1$) estimator, chosen to match what
`scipy.stats.gaussian_kde` uses internally. `site_ellipsoids()` computes all sites at once from six
`np.bincount` reductions (one per independent component) and divides by $\max(n-1,1)$;
`pca_kde_volume()` uses `np.cov(centered, rowvar=False, bias=False)`. JS: `covariance3()` accumulates
the six components and divides by $\max(n-1,1)$. `test_site_ellipsoids_batch_matches_direct_covariance`
asserts the batched form equals `np.cov(..., bias=False)` to `rtol=1e-9`. The $\max(\cdot,1)$ means a
one-atom "site" yields a zero tensor rather than a division by zero — see the degenerate-site
walkthrough at the end of Step 4.

**Where the mean is removed.** The two engines differ in bookkeeping, not in result. Python centres
**once**, in `load_site_displacements()` (`centered = offsets - site_mean[site_index]`), so
`site_ellipsoids()` and its kurtosis pass operate on already-centred arrays and subtract nothing
further. The JS `covariance3()` and `excessKurtosisPca()` re-derive the mean of the array they are
handed — numerically $\approx\mathbf{0}$ — and subtract it again. Every formula in this section is
written with the explicit $-\bar{\mathbf{u}}$; in Python that term is identically zero by
construction, and the results agree to round-off.

**Eigen-decomposition.** $\mathbf{U}=\sum_a\lambda_a\,\mathbf{p}_a\mathbf{p}_a^{\!\top}$ with
$\lambda_1\ge\lambda_2\ge\lambda_3\ge0$.

- Python `_eigen_decomposition()`: `np.linalg.eigh`, reverse to descending, clip
  $\lambda\leftarrow\max(\lambda,0)$ (round-off can make a near-zero eigenvalue slightly negative),
  transpose so axes are **rows**.
- JS `eigenDecomposition()`: `jacobiEigenSymmetric()` — cyclic Jacobi rotations over the pairs
  $(0,1),(0,2),(1,2)$, at most **50 sweeps**, early exit when
  $|a_{01}|+|a_{02}|+|a_{12}| \le 10^{-15}\,\lVert\mathbf{U}\rVert_F$, and an individual rotation
  skipped when $|a_{pq}|<10^{-300}$ (which would otherwise divide by zero forming
  $\theta=(a_{qq}-a_{pp})/2a_{pq}$). The convergence test is **relative** to the Frobenius norm, so a
  matrix of any scale is diagonalised to the same precision. (Before 1.0 it was an absolute
  $<10^{-18}$, which returned a round-off covariance of $\sim10^{-28}$ Å² — an average configuration —
  undiagonalised with identity axes while `eigh` returned rotated noise, so the two engines disagreed
  on axes, anisotropy and the `degenerate` flag for the same file.)

  Cyclic Jacobi was chosen over a closed-form $3\times3$ eigenvalue formula because it stays accurate
  on the near-degenerate clouds (planar or linear disorder) that trip cubic-root solvers, and a few
  sweeps of a $3\times3$ matrix cost nothing next to the KDE.

Both then apply the same **sign and handedness canonicalisation** (`_canonical_axes` in Python, the
sign + determinant block of `eigenDecomposition` in JS): each axis is flipped so its
largest-magnitude component is positive, then the third axis is negated if $\det\mathsf{P}<0$. It
exists because an eigenvector is defined only up to a sign and the sign LAPACK returns is
build-dependent. **The rule, its order-of-operations consequence for PC3, the one (unreachable)
divergence between the engines, and what a signed axis does and does not mean physically are set out
in "Principal axes in the crystallographic frame", Step 3** — that section owns direction semantics.
Here it is enough that the frame is orthonormal and right-handed:
`test_axes_are_orthonormal_and_right_handed` asserts $\mathsf{P}\mathsf{P}^{\!\top}=\mathsf{I}$ to
`atol=1e-9` and $\det\mathsf{P}=1$ to 9 places. For a near-isotropic cloud the axes remain arbitrary
within the degenerate subspace — no convention can fix that, and the code says so.

**Agreement between the two eigensolvers.** Both compute the exact symmetric decomposition to machine
precision, so eigenvalues agree at the $10^{-15}$ relative level; eigenvectors agree only up to the
arbitrariness of a (near-)degenerate subspace. This is **not** asserted by any cross-engine test —
the JS test (`eigenDecomposition` "recovers a known symmetric matrix") checks a hand-built
rotated `diag(4,2,1)` to 6 decimals, and the Python tests check orthonormality and ordering. There is
no golden-file parity test between `pca_kde.py` and `pcaKde.js` (unlike the
[Auto StoG](auto-stog.md) engines, which do have one — `autoscale_fixture.json`).

---

### Step 4 — ADP readouts

Everything in the *Displacement statistics* panel comes from Step 3 (`site_ellipsoids()` /
`siteEllipsoids()`):

| readout | formula | units | code key |
|---|---|---|---|
| Covariance $\mathbf{U}$ | $U_{ab}$ above, in Cartesian $x,y,z$ | Å² | `covariance` |
| $\lambda_a$ | eigenvalues, descending | Å² | `eigenvalues` |
| Principal axes | $\mathbf{p}_a$, unit Cartesian, one per row | — | `axes` |
| RMS amplitude | $\sigma_a=\sqrt{\lambda_a}$ | Å | `rms` |
| Semi-axes | $k(p)\sigma_a$ | Å | `semiAxes` |
| $U_\mathrm{iso}$ | $U_\mathrm{eq}=\tfrac13\sum_a\lambda_a=\tfrac13\,\mathrm{tr}\,\mathbf{U}$ | Å² | `uIso` |
| $B_\mathrm{iso}$ | $8\pi^2U_\mathrm{eq}$ | Å² | `bIso` |
| $\langle u^2\rangle^{1/2}$ | $\sqrt{U_\mathrm{eq}}$ | Å | `rmsIso` (not displayed) |
| Anisotropy | $\sqrt{\max(\lambda_1,10^{-30})/\max(\lambda_3,10^{-30})}$; `null` for a zero-spread site | — | `anisotropy` |
| Zero-spread flag | $\lambda_1 <$ `ZERO_SPREAD_VARIANCE` $=10^{-8}$ Å² (RMS $<10^{-4}$ Å on every axis) | — | `zeroSpread` |
| Degenerate flag | `zeroSpread` **or** $\lambda_3/\max(\lambda_1,10^{-30}) < 10^{-6}$ (`DEGENERATE_RATIO`) | — | `degenerate` |

$\mathrm{tr}\,\mathbf{U}$ is rotation-invariant, so $U_\mathrm{eq}=\mathrm{tr}\,\mathbf{U}/3$ is the
standard equivalent isotropic displacement parameter and $B_\mathrm{iso}=8\pi^2U_\mathrm{eq}$ the
standard Debye–Waller conversion. Both engines compute exactly these two lines and nothing else.

Note that **both** the numerator and the denominator of `anisotropy` carry the $10^{-30}$ floor, and
that the same floored "largest" also serves as the denominator of the degeneracy ratio. The
unfloored expression $\sigma_1/\sigma_3$ agrees with it for every real site but not in the limit —
see the walkthrough below.

**Zero spread.** In an RMCProfile `*AVERAGE.rmc6f`, or an ideal starting configuration, every copy of
a site sits at the same offset and the covariance is floating-point round-off
($\sim10^{-26}$–$10^{-29}$ Å²; e.g. `data/5K_try1/GaNb4Se8_5KAVERAGE.rmc6f`). A ratio test cannot see
that — ratios of noise are noise — so both engines apply one **absolute** floor,
`ZERO_SPREAD_VARIANCE = 1e-8` Å² (RMS $10^{-4}$ Å, five orders below any physical ADP and far above
round-off or the $\sim10^{-5}$ Å quantisation of a 7-decimal coordinate in a 100 Å box). Such a site
reports `zeroSpread: true`, `degenerate: true`, `axes: null`, `anisotropy: null`,
`excessKurtosis: [null, null, null]` and `nonGaussianity: null` (JSON `null`, never a `NaN` token);
$\mathbf{U}$, $\lambda_a$, `rms`, $U_\mathrm{iso}$ and $B_\mathrm{iso}$ keep their (round-off) values.
The panel prints *no displacement* in the anisotropy cell and explains the missing axes, and the KDE
request raises *"displacement cloud has zero spread"* (Step 6).

#### Step 4b — What the statistics panel actually prints

`PcaKdePage.jsx` renders these numbers in the *Displacement statistics* panel; everything in it is
read from `selectedEllipsoid` — the row of the `sites` payload matching the selected reference
number — i.e. from the **full-cloud** statistics, never from the KDE volume's subsampled fit.

*Principal axes table*, one row per axis, tinted with `PC_CSS_COLORS`
(`['#d64545','#3fa34d','#3f7fd6']` in [`sceneAxes.js`](../../web_app/frontend/src/components/sceneAxes.js),
the CSS twins of the 3D triad colours):

| Column | Value | Units | Decimals |
|---|---|---|---|
| row label | `PC1`/`PC2`/`PC3` with a coloured dot | — | — |
| `x`, `y`, `z` | components of $\mathbf{p}_a$ (`axes[a]`) | — (unit vector) | 3 |
| `λ (Å²)` | $\lambda_a$ | Å² | 4 |
| `RMS (Å)` | $\sigma_a$ | Å | 3 |
| `κ` | $\kappa_a$ (Step 12b); `—` when the axis is not resolved (`axisResolved[a]` false) or undefined | — | 2 |

Beside it a *Covariance U (Å²)* table prints $\mathbf{U}$ as a $3\times3$ Cartesian matrix to 4
decimals. The *Summary* block prints $U_\mathrm{iso}$ (4), $B_\mathrm{iso}$ (3), anisotropy (2) and
non-Gaussianity (2); the anisotropy cell appends `· degen.` when the `degenerate` flag is set, and
non-Gaussianity falls back to the KDE result's `nonGaussianity` when the site record lacks one.
Non-finite values render as an em dash (`numberFormat` returns `'—'`).

A fourth column, *Crystal orientation*, prints where each principal axis points in the crystal:
∠a/∠b/∠c, the angles to the unit-cell edges (1 decimal, degrees, the closest edge shaded), and the
direct-lattice direction $[u\,v\,w]$ it runs along (2 decimals, largest component $\pm1$, components
with $\lvert\cdot\rvert<5\times10^{-3}$ printed as `0.00`, never `-0.00`, by `uvwFormat`). Each row is
the sign-canonical representative from `crystalOrientationRows()` — the sense that makes the closest
edge acute. The column is omitted when the payload carries no lattice metadata or the cell is
singular, and a row is dimmed when its axis is not resolved (Step 12b). The algebra, the sign fold
and what $[u\,v\,w]$ does and does not mean are in "Principal axes in the crystallographic frame",
Steps 7–8.

The panel heading carries a metadata line, rendered exactly as

```
Volume {grid}³ · fit {fitCount}/{count} · captured mass {mass·100}%[ · iso @ {isoPercent}% mass]{ · browser| · server}
```

with `fitCount`/`count` through `toLocaleString()`, captured mass to **1 decimal**, and the
`· iso @ …` fragment conditional on `Number.isFinite(kde.massLevels[isoPercent]?.level)`. The whole
line renders only when a `kde` payload exists.

**A degenerate site, end to end.** For a site with a single box copy: $\mathbf{U}=\mathbf 0$, so
$\lambda_a=0$, `rms` $=$ `semiAxes` $=0$, $U_\mathrm{iso}=B_\mathrm{iso}=0$; `zeroSpread` and
`degenerate` are true; `axes`, `anisotropy`, every $\kappa_a$ and the non-Gaussianity are `null`, so
the Summary prints *no displacement* and *Non-Gaussianity* `—`, and the Principal-axes table shows a
one-line explanation instead of axes. The KDE request for that site **throws** — `pca_kde_volume()` /
`pcaKdeVolume()` raise *"a 3D KDE needs at least four points"* for $n<4$ and *"displacement cloud has
zero spread"* for $\lambda_1<10^{-8}$ Å² (Step 6) — and the page shows the red `pca-badge is-error`
overlay in the viewport. Because the sites table and the KDE volume are two independent requests,
the tables render while the 3D view shows only the error badge, and the metadata line above (which
needs `kde`) does not render at all.

**What the code does *not* compute.** There is no conversion of the Cartesian $\mathbf{U}$ to the
crystallographic $U^{ij}$ / $U_\mathrm{cif}$ basis (components on the reciprocal-cell axes) and no
$\beta_{ij}$ form. The tensor shown in the UI is labelled "Covariance U (Å²) … in Cartesian (x, y, z)
axes", which is accurate. Any comparison against CIF `U_11 … U_23` values must be done by the reader
— **the app never does this** — with

$$\mathbf{U}_\mathrm{cif}=\mathsf{D}^{-1}\,\mathsf{A}^{-\top}\,\mathbf{U}_\mathrm{cart}\,\mathsf{A}^{-1}\,\mathsf{D}^{-1},
\qquad \mathsf{D}=\mathrm{diag}\big(\lvert\mathbf a^\ast\rvert,\lvert\mathbf b^\ast\rvert,\lvert\mathbf c^\ast\rvert\big),$$

where $\mathsf{A}$ has rows $\mathbf a,\mathbf b,\mathbf c$ (Å, the notation table), the rows of
$\mathsf{A}^{-\top}$ are the reciprocal vectors $\mathbf a^\ast,\mathbf b^\ast,\mathbf c^\ast$ (no
$2\pi$), and $\mathsf{D}$ their lengths. Derivation: a Cartesian displacement is
$\mathbf u=\mathsf{A}^{\!\top}\Delta\mathbf f$, so $\langle\Delta\mathbf f\,\Delta\mathbf f^{\top}\rangle
=\mathsf{A}^{-\top}\mathbf{U}_\mathrm{cart}\mathsf{A}^{-1}$, and the CIF convention defines
$\langle\Delta f_i\Delta f_j\rangle=a^\ast_i a^\ast_j\,U^{ij}$. For an orthogonal cell aligned with the
Cartesian axes ($\mathsf A=\mathrm{diag}(a,b,c)$, $\mathsf D=\mathsf A^{-1}$) $\mathbf U_\mathrm{cif}=
\mathbf U_\mathrm{cart}$; for anything else the components differ. (Before 1.0 this page gave
$\mathsf{A}^{-1}\mathbf{U}_\mathrm{cart}\mathsf{A}^{-\top}$ — transposes in the wrong order and no
reciprocal-length normalisation, so it returned dimensionless $\langle\Delta f\Delta f\rangle$-like
numbers about $U/a^2$, e.g. $2.5\times10^{-4}$ for a 0.012 Å² $U_{11}$ in a 7 Å cell.)

Note also that `pca_kde_volume()` returns its own `uIso`/`bIso`/`anisotropy` computed from the
**floored** eigenvalues (Step 6) and from the possibly-subsampled cloud; the UI table reads the
`site_ellipsoids()` values, which are unfloored and use every atom.

---

### Step 5 — Probability-ellipsoid scaling (the $\chi^2$ quantile)

For a trivariate Gaussian, the squared Mahalanobis radius
$r^2=\sum_a q_a^2/\lambda_a$ is $\chi^2$-distributed with 3 degrees of freedom. The surface enclosing
probability $p$ is therefore the ellipsoid with semi-axes

$$k(p)\,\sigma_a,\qquad k(p)=\sqrt{F^{-1}_{\chi^2_3}(p)} .$$

**Python** (`probability_scale()`) evaluates $\sqrt{\texttt{scipy.stats.chi2.ppf}(p,3)}$ exactly, and
rejects $p\notin(0,1)$. The crystallographic 50% convention gives
$k(0.5)=\sqrt{2.3659739}=1.5381722$, asserted to 5 places by
`ProbabilityScaleTests::test_fifty_percent_is_crystallographic_constant`.

**JavaScript** (`probabilityScale()` → `chiSquare3Quantile()`) inverts the closed-form
$\chi^2_3$ CDF,

$$F(x)=\mathrm{erf}\!\big(\sqrt{x/2}\big)-\sqrt{2x/\pi}\,e^{-x/2}=P\big(\tfrac32,\tfrac x2\big),\qquad
f(x)=\sqrt{x/2\pi}\,e^{-x/2},$$

by safeguarded Newton iteration (bracket $[0,2^m]$, bisection whenever a Newton step leaves the
bracket, stop at a relative step of $4\times10^{-16}$). $F$ and its complement are evaluated by
`chiSquare3Tails()` without cancellation: the positive series of the regularised incomplete gamma
function $P(a,t)$ below $t=a+1$ and the modified-Lentz continued fraction of $Q(a,t)=1-P$ above it
(Numerical Recipes `gser`/`gcf`, $a=\tfrac32$, $t=x/2$); the residual is taken from the lower tail for
$p\le\tfrac12$ and from the upper tail (with the exact $1-p$) above, so both ends keep full relative
precision. The result equals `scipy.stats.chi2.ppf(p, 3)` to $\sim10^{-15}$;
`pcaKdeProbabilityScale.test.js` pins it at $10^{-10}$ against scipy on $p=0.01,0.02,\dots,0.99$ plus
0.6827, 0.995, 0.9973 and 0.999, and checks monotonicity out to $p=10^{-12}$ and $1-10^{-12}$.

*History.* Before 1.0 the browser interpolated a 12-node table linearly in $q=k^2$, clamping outside
$[0.10, 0.9973]$. $q(p)$ is convex above $p\approx0.2$, so the chords overestimated — up to
$+3.5\%$ in every semi-axis at $p=0.97$ ($+1.2\%$ at 0.85) — and the $p=0.6827$ node held
$q=3.5058779$ instead of $3.5268222$ ($-0.30\%$ in $k$). The drawn ellipsoid and the KDE shell sampled
on it therefore differed between the static and the Flask app at every non-node slider value; they
now agree.

The scale factor multiplies only the **drawn** ellipsoid (`semiAxes`); no other quantity depends on
it.

---

### Step 6 — Bandwidth and the sampling box

**Bandwidth factor** (`_bandwidth_factor()` / `bandwidthFactor()`), for $d=3$ and $n$ = `fitCount`:

$$f_\mathrm{Scott}=n^{-1/(d+4)}=n^{-1/7},\qquad
f_\mathrm{Silverman}=\left[\frac{n(d+2)}{4}\right]^{-1/7}=\left(\tfrac{5n}{4}\right)^{-1/7},$$

then multiplied by `bw_scale` (default 1.0; exposed by the API as `bwScale`, **not** exposed in the
UI). A positive number may be passed instead of a rule name and is used as $f$ directly. These are
exactly SciPy's `gaussian_kde` covariance factors for unweighted points. For a 10³ box,
$n=1000\Rightarrow f_\mathrm{Scott}=0.3728$; for $n=216$, $f=0.4640$; for the 20 000-point cap,
$f=0.2430$.

**Bandwidth matrix.** SciPy's estimator uses $\mathsf{H}=f^2\mathbf{C}$ with $\mathbf{C}$ the data
covariance. In the PCA frame $\mathbf{C}=\mathrm{diag}(\lambda_a)$, so

$$\mathsf{H}=\mathrm{diag}(h_a^2),\qquad h_a=f\sqrt{\lambda_a}=f\sigma_a\;[\text{Å}].$$

**Eigenvalue floor.** A cloud that is flat along one direction drives $\lambda_3\to0$ and makes
$\mathsf{H}$ singular, so `pca_kde_volume()` floors
$\lambda_a\leftarrow\max(\lambda_a,\lambda_1\cdot\texttt{EIGENVALUE\_FLOOR\_RATIO})$ with
`EIGENVALUE_FLOOR_RATIO = 1e-8` — bounding $\mathrm{cond}(\mathsf{H})$ at $10^8$ (a bandwidth ratio
of $10^4$) without visibly moving a well-conditioned site. The pre-floor ratio is what sets the
`degenerate` flag. If $\lambda_1<$ `ZERO_SPREAD_VARIANCE` ($10^{-8}$ Å², Step 4) the call raises
*"displacement cloud has zero spread (RMS below 1e-4 A on every axis)"* — before 1.0 the test was
$\lambda_1\le0$, which a round-off covariance passes, so an average configuration drew a
$10^{-14}$ Å "cloud" with $v_\mathrm{max}\approx10^{42}$ Å⁻³; fewer than 4 points raises *"a 3D KDE
needs at least four points"*.

**Input validation.** `pca_kde_volume()` / `pcaKdeVolume()` reject non-finite points
(*"displacement cloud contains non-finite coordinates"*) and require `extent`, `bw_scale` and a
numeric `bw` to be positive **and finite**, and `grid` finite. Before 1.0 the checks were `<= 0`,
which NaN passes: a NaN `extent`, `bwScale` or `bw` returned an all-NaN volume that Flask serialised
as bare `NaN` tokens (HTTP 200, invalid JSON for `JSON.parse`), and JS accepted `Infinity`. Both
routes now answer a bad parameter with HTTP 400.

**Box half-widths.** The KDE convolves the cloud with the kernel, so the *estimate*'s variance along
axis $a$ is $\lambda_a+h_a^2=\lambda_a(1+f^2)$. The sampling box is sized on that broadened width:

$$w_a=\texttt{extent}\cdot\sigma_a\sqrt{1+f^2}\quad(\texttt{halfWidths}),\qquad
\texttt{boxHalfWidths}=\begin{cases}(\max_b w_b)\,(1,1,1) & \texttt{cubicBox}\\ (w_1,w_2,w_3) & \text{otherwise.}\end{cases}$$

The volume, `mass`, `massLevels` and the PC wall projections are **always** sampled on the per-axis
box $\pm w_a$, on whose $G$ nodes every axis's kernel $h_a=f\sigma_a$ is resolved
($\Delta_a/h_a=2\,\texttt{extent}\sqrt{1+f^2}/(f(G-1))\approx0.6$ at the UI defaults, whatever the
anisotropy). `cubicBox` only sizes the **display** box, `boxHalfWidths`: the cube the page draws the
shadow box and its walls on. The page always requests `cubicBox: true`, with `extent` from a slider
(2–5 σ, step 0.5, default **4**); the engine defaults are `extent = 3.0`, `cubic_box = False`. For a
Gaussian the fraction of mass inside a $\pm e\sigma$ box per axis is $\mathrm{erf}(e/\sqrt2)^3$:
0.870 at $e=2$, **0.992** at $e=3$, 0.99981 at $e=4$.

*Why the cube is display-only (changed in 1.0).* Until 1.0 `cubicBox` also sampled every axis over
the cube, with the same $G$ nodes, so the grid spacing along a short axis grew to
$\Delta=2\,\texttt{extent}\,\sigma_1\sqrt{1+f^2}/(G-1)$ while its kernel stayed $f\sigma_3$. Once
$\sigma_1/\sigma_3\gtrsim3$ ($G=40$) the narrow kernels fell between nodes and the rectangle-rule mass,
the HDR levels $\tau(p)$ (the isosurface) and the shell all aliased: 36% captured mass at
$\sigma_1/\sigma_3=20$, and for a planar or collinear cloud **0%** on the UI's even grids (no node on
the plane, no isosurface) or $\approx2\times10^5$% on an odd one (one node on the plane weighted by the
wide $\Delta$). The volume is now independent of `cubicBox` —
`CubicDisplayBoxTests` / `pcaKdeRegressions.test.js` assert identical densities and unit mass for
those clouds.

**Grid.** $G$ points per axis on `np.linspace(-w_a, +w_a, G)`, clamped to $8\le G\le128$ in both
engines. UI options 24/32/40/48/56/64, default **40**; engine default 48. Grid spacing
$\Delta_a=2w_a/(G-1)$; `cellVolume` $=\prod_a\Delta_a$ (Å³).

---

### Step 7 — The separable 3D Gaussian KDE

**The estimator.** For kernel centres $\mathbf{y}_m$ (the cloud, mean-subtracted and projected onto
the principal axes, $\mathbf{q}_m=(\mathbf{y}_m-\bar{\mathbf{y}})\mathsf{P}^{\!\top}$) SciPy's
`gaussian_kde` evaluates at a point $\mathbf{q}$

$$\hat{\rho}(\mathbf{q})=\frac{1}{n(2\pi)^{3/2}|\mathsf{H}|^{1/2}}
\sum_{m=1}^{n}\exp\!\left[-\tfrac12(\mathbf{q}-\mathbf{q}_m)^{\!\top}\mathsf{H}^{-1}(\mathbf{q}-\mathbf{q}_m)\right].$$

**Why the PCA frame makes it separable.** $\mathsf{H}=f^2\mathbf{C}$ is diagonal in the eigenbasis of
$\mathbf{C}$ — which is exactly the frame this analysis wants anyway. The quadratic form collapses to
$\sum_a (q_a-q_{ma})^2/h_a^2$ and the Gaussian **factorises**:

$$\hat{\rho}(\mathbf{q})=\frac{1}{n(2\pi)^{3/2}h_1h_2h_3}
\sum_m \prod_{a=1}^{3}\exp\!\left[-\frac{(q_a-q_{ma})^2}{2h_a^2}\right].$$

Sampling on a grid **aligned with those axes** (which is what `axisCoords` are) means each factor
depends on only one grid index. Define three kernel matrices, $K^{(a)}\in\mathbb{R}^{G\times n}$,

$$K^{(a)}_{i m}=\exp\!\left[-\tfrac12\!\left(\frac{g^{(a)}_i-q_{ma}}{h_a}\right)^{\!2}\right]
\qquad(\texttt{\_kernel\_matrix} / \texttt{kernelMatrix}),$$

and the volume is a tensor contraction over the *sample* index:

$$\rho_{ijk}=\underbrace{\frac{1}{n(2\pi)^{3/2}h_1h_2h_3}}_{\texttt{norm}}
\sum_{m=1}^{n}K^{(1)}_{im}K^{(2)}_{jm}K^{(3)}_{km}\qquad[\text{Å}^{-3}].$$

**This is not an approximation.** It is an exact algebraic rearrangement of the same estimator; only
the *evaluation points* are chosen to be axis-aligned.

**Cost.** Building the three kernel matrices needs $3nG$ exponentials instead of $nG^3$ for direct
evaluation ($n=1000$, $G=48$: $1.44\times10^5$ vs $1.11\times10^8$; the `pca_kde.py` module docstring
rounds the same comparison to `2e5` rather than `1.1e8`). The remaining contraction is still $O(nG^3)$ multiply–adds, but they are dense
linear-algebra operations rather than transcendental calls plus a 3×3 quadratic form. Both engines
assemble one PC3 slice at a time so no $n\times G^3$ temporary is ever materialised:

- Python: `density[:, :, k] = (kernels[0] * kernels[2][k]) @ kernels[1].T` — a $(G\times n)\cdot(n\times G)$
  BLAS `dgemm` per $k$;
- JS: the same, with an explicit `xzRow` buffer of $K^{(1)}_{im}K^{(3)}_{km}$ reused across all $j$.

**Memory layout.** `density` is flat, **C order over (PC1, PC2, PC3)**:
`index = (i*grid + j)*grid + k`. Python returns it as a Python list (via JSON); JS returns a
`Float64Array`.

**Verification — what the tests actually assert.**

- `tests/test_pca_kde.py::test_volume_matches_scipy_gaussian_kde`: for `bw ∈ {"scott", "silverman",
  0.35}` the returned volume is rebuilt into Cartesian points
  ($\mathbf{x}=\bar{\mathbf{x}}+\mathbf{q}\,\mathsf{P}$) and compared against
  `scipy.stats.gaussian_kde(cloud.T, bw_method=bw)` evaluated at exactly those points, with
  `rtol=1e-9, atol=1e-12`. That is the round-off-equality claim, against the reference implementation
  itself.
- `pcaKde.test.js` "matches the brute-force anisotropic 3D KDE to round-off": the JS engine is
  compared against a **full multivariate** KDE written inside the test (explicit 3×3 inverse of
  $\mathsf{H}=f^2\mathbf{C}$, one exponential per (grid point, sample) pair) with max absolute and
  max relative error both $<10^{-9}$. Note this reference is the test's own brute force, *not* SciPy —
  the browser engine is never compared to SciPy directly.

`axes` and `mean` are returned precisely so the caller can undo the rotation:
$\mathbf{x}_\mathrm{cart}=\texttt{mean}+\mathbf{q}\cdot\texttt{axes}$.

---

### Step 8 — Isosurface level selection

`_iso_levels()` / `isoLevels()` produce two families of thresholds from the sampled volume.

**Mass levels (the ones the page uses).** Sort all $G^3$ density values descending into
$\rho_{(1)}\ge\rho_{(2)}\ge\dots$, form the running quadrature sum
$S_r=\Delta V\sum_{t\le r}\rho_{(t)}$ with $\Delta V=$ `cellVolume`, and let $S_\infty$ be the total
(`mass`). For each requested $p$,

$$\tau(p)=\rho_{(r^\ast)},\qquad r^\ast=\min\{r: S_r\ge p\,S_\infty\}$$

(`np.searchsorted(..., side="left")` in Python; an equivalent lower-bound binary search in JS). The
surface $\{\rho=\tau(p)\}$ is the **highest-density region** — the smallest volume containing
fraction $p$ of the density in the box — which is exactly what makes it comparable with the $p\%$
ellipsoid. Requested probabilities default to `np.linspace(0, 1, 101)` (i.e. $p=0.00,0.01,\dots,1.00$),
so `massLevels[i]` is the level for $p=i/100$ and the UI indexes it directly with the integer
percentage.

Two honest caveats: $p$ is a fraction of the **captured** mass, not of the estimator's unit total; and
$S_r$ is a plain rectangle-rule sum over grid nodes, so $\tau$ inherits an $O(\Delta^2)$ quadrature
error.

**`mass`.** $S_\infty=\Delta V\sum_{ijk}\rho_{ijk}$ — the fraction of the KDE's unit mass the box
holds, i.e. the truncation error of the box. Displayed as "captured mass NN.N%".
`test_mass_recovers_probability_normalization` asserts $>0.99$ at `extent=4, grid=64` (JS: same at
`extent=4, grid=56`).

**Density levels.** `densityLevels[i] = vmin + p*(vmax − vmin)` — a linear ramp on the density range.
Computed and returned by both engines but **not used by the UI**.

**UI wiring.** The slider *Level* (1–99%, default **50** — the same level as the ellipsoid's default
probability, because the two surfaces are only comparable at equal $p$; the pre-1.0 default of 25%
against a 50% ellipsoid put a Gaussian site's surface ~16% *inside* the ellipsoid, the very cue the
tooltip calls anharmonic) picks `kde.massLevels[isoPercent].level`,
with a fallback of `vmax * (1 − isoPercent/100)` when mass levels are unavailable
(`PcaKdePage.jsx`, the `massLevel` line in the scene-rebuild effect). `test_mass_levels_bracket_the_cloud`
pins the monotonicity: higher enclosed probability ⇒ lower density threshold.

---

### Step 9 — Marching cubes

[`web_app/frontend/src/workers/marchingCubes.js`](../../web_app/frontend/src/workers/marchingCubes.js)
→ `marchingCubes(field, nx, ny, nz, isoLevel)` contours the scalar field. It is the classic
Lorensen–Cline algorithm with the public-domain Paul Bourke lookup tables: a 256-entry `EDGE_TABLE`
(bitmask of intersected cube edges) and a 256-row `TRI_TABLE` (triangles as edge-index triples).
Three.js's own `MarchingCubes` helper is a metaball generator that rebuilds a field from point
sources and cannot contour an existing volume, which is why this module exists.

Per cell $(i,j,k)$ over $[0,n_x-2]\times[0,n_y-2]\times[0,n_z-2]$:

1. classify the 8 corners: bit $c$ set when `value < isoLevel` (inside = *below* the level; the
   tables are consistent with that convention);
2. for each edge flagged by `EDGE_TABLE[cubeIndex]`, place the vertex by **linear interpolation**
   $t=(\tau-v_a)/(v_b-v_a)$, falling back to $t=0.5$ when $|v_b-v_a|<10^{-12}$;
3. interpolate a per-vertex normal from central-difference field gradients at the two corners
   (clamped at the borders), normalised and **negated** — the field increases inward for a density
   blob, so $-\nabla\rho$ points outward;
4. emit triangles as a non-indexed soup of positions in **grid-index space**.

`PcaKdePage.jsx` maps each vertex to Cartesian with `pcaVertexToCartesian()`, which assumes the
uniform spacing `axisCoords[a][1] − axisCoords[a][0]`:
$\mathbf{x}=\texttt{mean}+\sum_a\big(g^{(a)}_0+f_a\Delta_a\big)\mathbf{p}_a$.

**Discrepancy worth naming:** the page then calls `geometry.computeVertexNormals()` and **discards**
the gradient normals `marchingCubes` computed and returned. The rendered shading therefore uses
face-averaged normals, not the analytic field gradient. The module's own doc comment advertises
gradient normals as its output; they are exercised by `marchingCubes.test.js` ("normals point outward
from a density blob", >95% of vertices with $\mathbf{n}\cdot\hat{\mathbf{r}}>0$) but not by the app.

Other assertions in that test file: a sphere field contours to vertices within 0.6 grid cells of the
true radius, and a field that never crosses the level returns `count === 0`.

Rendering: `MeshPhongMaterial`, colour = the active colormap sampled at $t=0.72$, `opacity 0.55`,
`DoubleSide`.

---

### Step 10 — Wall projections and their contours

#### 10a. PC-frame walls — the analytic marginal

Because the 3D estimator is separable, integrating out one axis analytically leaves the 2D KDE with
the *same* bandwidth sub-block. `_projection()` / `projection()` therefore computes, for the plane
spanned by axes $(u,v)$,

$$\rho^{(uv)}_{ij}=\frac{1}{n\,2\pi h_u h_v}\sum_m K^{(u)}_{im}K^{(v)}_{jm}
\qquad[\text{Å}^{-2}],$$

a single $(G\times n)\cdot(n\times G)$ matrix product per plane. Three planes are returned: `pc12`,
`pc13`, `pc23`, each with `density` (indexed `[first][second]`), `extent` (the two axes' coordinate
ranges, Å), `axes`, `bandwidth`, `vmax`.

This is the **honest shadow of the displayed volume**, not an independently bandwidth-optimised 2D
estimate — the marginal keeps the 3D rule's $h_u,h_v$. Both suites verify it against a fine Riemann
sum of the volume over the dropped axis:
`test_projection_is_true_marginal_of_the_volume` (Python, `grid=96, extent=5`, `rtol=2e-3, atol=1e-6`)
and the JS twin (`grid=80, extent=5`, max relative error $<3\times10^{-3}$). The tolerance is loose
because the check itself is a discrete sum, not because the marginal is approximate.

**Placement** (`makeProjectionWall`): on the far face of the **display** box, at $-1.06\,W_t$ along
the remaining axis $t=3-u-v$ ($W$ = `boxHalfWidths`, the cube), pushed 6% outward so it sits just past
the cloud, with basis $(\mathbf{p}_u,\mathbf{p}_v,\mathbf{p}_u\times\mathbf{p}_v)$. Each face gets two
coplanar quads (`MeshBasicMaterial`, `DoubleSide`, `opacity 0.96`, `depthWrite: false`): a
$2W_u\times2W_v$ backing in the colormap's zero colour (`renderOrder 0`), and the density texture sized
to the range the projection was **evaluated** on — its `extent`, the per-axis box $2w_u\times2w_v$ —
centred on it (`renderOrder 1`). A short axis therefore shows its true, narrow marginal inside a
full-size face instead of an under-sampled one stretched over it.
Texture (`projectionTexture`): `density/(vmax || 1)` clipped to $[0,1]$ through the colormap LUT, so
**each wall is normalised to its own maximum** (walls are not comparable to each other in absolute
density). The texture row index is written flipped, `(nSecond − 1 − j)`, to cancel `CanvasTexture`'s
vertical flip on upload, so the second axis increases along the plane's local $+Y$. A wireframe box
at $\pm W_a$ (`makeBoundingBox`) completes the shadow-box. Contour lines map their grid indices through
the same per-projection range (`projectionRange`).

**Contours** (`makeProjectionContours`): marching *squares* (`contourSegments`) at the fixed
fractions

$$\{0.10,\;0.25,\;0.40,\;0.55,\;0.70,\;0.85\}\times v_\mathrm{max}^{(uv)},$$

i.e. **relative levels of that plane's own peak**, *not* enclosed-mass levels — do not read them as
probability contours. Each cell contributes 0, 1, or 2 segments (the 4-crossing ambiguous case is
emitted as two segments without a saddle-point disambiguation, so a saddle may be connected the
"wrong" way). Lines are drawn at $-1.05\,w_t$ (just inside the wall) with `renderOrder = 2` so the
translucent wall cannot paint over them.

#### 10b. Crystal-frame walls — re-binned, not analytic

When the *Frame* switch is set to **Crystal**, `orthonormalCrystalFrame()` (in `PcaKdePage.jsx`)
builds a Gram–Schmidt orthonormal frame from the unit cell —
$\mathbf{e}_0\parallel\mathbf{a}$, $\mathbf{e}_1$ in the $\mathbf{a}$–$\mathbf{b}$ plane,
$\mathbf{e}_2=\mathbf{e}_0\times\mathbf{e}_1$, with $\mathbf{c}$ never referenced. It is used **only**
to draw the box, its walls and the crystal-mode default camera; the a/b/c rods and the
look-down-a/b/c cameras use the true (possibly oblique) cell edges. What its three faces do and do
not mean for an oblique cell, and its missing singularity guard, are covered in "Principal axes in
the crystallographic frame", Step 9b.

`projectDensityOntoFrame(kde, frame, half, nBins)` then re-bins the *same* sampled volume. For every
**grid node** $(i,j,k)$ — `axisCoords` are `linspace(-w_a, +w_a, G)`, so the samples are nodes, with
the first and last sitting exactly on the box faces, **not** voxel centres — it maps the node to a
Cartesian offset $\mathbf{c}$ from the cloud mean through `axes`, takes the three in-frame
coordinates $\mathbf{c}\cdot\mathbf{e}_t$, converts each to a continuous bin coordinate

$$\texttt{toBin}(c)=\frac{c+\texttt{half}}{\texttt{span}}\,(n_\mathrm{bins}-1),\qquad
\texttt{half}=\max_a w_a,\ \ \texttt{span}=2\,\texttt{half},\ \ n_\mathrm{bins}=G,$$

and **bilinearly splats** the weight $w=\rho_{ijk}\Delta V$ into three 2D accumulators (the
$\mathbf{e}_0$–$\mathbf{e}_1$, $\mathbf{e}_0$–$\mathbf{e}_2$ and $\mathbf{e}_1$–$\mathbf{e}_2$
planes, integrating out the remaining axis). Nodes with $\rho\le0$ are skipped. Each accumulator is
finally divided by $\texttt{binArea}=(\texttt{span}/(n_\mathrm{bins}-1))^2$ to give a surface density
in Å⁻². Three defensive fallbacks guard the arithmetic: `cellVolume || 1`, `span || 1`,
`binArea || 1`. So:

- these walls are a **numerical** projection of the discretised volume, whereas the PC walls are the
  analytic marginal;
- density that falls outside $\pm\max_a w_a$ in the crystal frame is silently dropped. The PCA
  sampling box's corners reach up to $\lVert\mathbf w\rVert\le\sqrt3\max_a w_a$ along a crystal
  direction, so corner cells can be lost. At the default `extent = 4` those corners carry a
  negligible fraction of the mass, but the crystal-frame walls do not integrate to the same total
  as the PC-frame ones;
- the bilinear splat adds a small smoothing of order one bin, so the crystal walls read blockier and
  slightly more diffuse than the PC ones;
- $\Delta V$ is the **node-spacing product** $(2w_a/(G-1))^3$, so $\sum\rho\,\Delta V$ is a
  node-sampled Riemann sum in which the edge nodes are over-weighted by roughly $2\times$ per axis —
  the same convention `mass` uses in both engines (Step 8), but not a cell-centred quadrature.

The crystal-frame projections are shaped exactly like the PC ones (`{pc12, pc13, pc23}` with
`density`, `axes`, `vmax`), so the same wall/contour builders draw them unchanged — including the
same per-plane `vmax` normalisation and the same relative contour fractions. Two naming traps: the
keys stay `pc12`/`pc13`/`pc23` in crystal mode, where they actually mean the
$\mathbf{e}_0$–$\mathbf{e}_1$, $\mathbf{e}_0$–$\mathbf{e}_2$ and $\mathbf{e}_1$–$\mathbf{e}_2$
planes; and `kde.cellVolume` is the KDE grid's node-spacing volume (Å³), not the crystallographic
unit-cell volume.

**Cost.** `crystalDisplay` is an $O(G^3)$ **main-thread** `useMemo` keyed on `[kde, unitCell]` only —
not on `axisFrame`. It therefore re-runs for every new volume even while the viewport is in PC mode,
where its result is unused.

---

### Step 11 — Density painted on the ellipsoid ("Shell")

`makeEllipsoidKdeSurface()` is the default density view. It takes a unit `SphereGeometry(1, 96, 64)`,
maps each vertex $\hat{\mathbf{s}}$ to the ellipsoid surface point
$\mathbf{q}=(k\sigma_1\hat s_x,\,k\sigma_2\hat s_y,\,k\sigma_3\hat s_z)$ in the PCA frame (using
`selectedEllipsoid.semiAxes`), bakes the world position
$\mathbf{x}=\texttt{mean}+\sum_a q_a\mathbf{p}_a$, and samples the KDE volume there by **trilinear
interpolation** on the flat grid (`sampleFieldTrilinear` in `marchingCubes.js`).

**Outside the sampled box there is no data.** The shell reaches $k(p)\sigma_a$ along each axis, the
box only $\texttt{extent}\cdot\sigma_a\sqrt{1+f^2}$, so a small Box with a high Level pokes the PC1
tips out ($k(p)>\texttt{extent}\sqrt{1+f^2}$: $p\gtrsim0.79$ at Box 2, $0.93$ at 2.5, $0.99$ at 3). The
sampler returns `NaN` there and the page paints those vertices neutral grey (`NO_DATA_RGB`), leaves
them out of the colour stretch, and adds a legend note naming the Box that would cover the shell
(`shellBoxNeeded`, $\lceil 2\,\texttt{extent}\max_a k\sigma_a/w_a\rceil/2$). Before 1.0 the sampler
clamped to the border node and painted the box-face density — nearer the centre, so 1.4–2.3× too
dense — as bright false caps along PC1.

Colouring is stretched to the **shell's own** min/max density (not the global $0..v_\mathrm{max}$), then
a user contrast gain is applied symmetrically about the mid-tone:

$$t=\mathrm{clip}\!\left[\,0.5+\left(\frac{\rho_v-\rho_\mathrm{min}^\mathrm{shell}}
{\rho_\mathrm{max}^\mathrm{shell}-\rho_\mathrm{min}^\mathrm{shell}}-0.5\right)\cdot\texttt{contrast},\;0,\;1\right],
\qquad \texttt{contrast}\in[0.5,3],\ \text{default }1.$$

Interpretation: for a Gaussian *population* an iso-probability ellipsoid is a level set of the
density, so an infinite sample would paint one colour; a systematic pattern (caps along an axis, a
band) marks where the real density departs from the harmonic reference. A **finite** cloud does not
paint one colour: at $n=1000$ the KDE on the 50% ellipsoid of a truly Gaussian cloud varies by about
$2\times$ ($(\max-\min)/\text{mean}\approx0.7$) from sampling noise alone, and **because the range is
auto-stretched, that noise is rendered with the full colour range**. Read the magnitude from the
numeric non-Gaussianity, not from the shell's colours; the Shell tooltip says so.

The wireframe ellipsoid itself (`showEllipsoid`) is a unit sphere transformed by
$\mathsf{P}^{\!\top}\mathrm{diag}(k\sigma_a)$ with translation `mean`, drawn as a translucent
wireframe (opacity 0.4) in a user-chosen colour.

---

### Step 12 — The non-Gaussianity readout

Code: `_shape_statistics()` in Python (called by `site_ellipsoids()`, batched — one `einsum` projects
every atom onto its own site's axes, then `bincount` moments — and by `pca_kde_volume()`), and its
line-for-line port `shapeStatistics()` in JS (called by `siteEllipsoids()` and `pcaKdeVolume()`). All
moments are population ($1/n$) moments of the cloud projected onto its principal axes,

$$q_{ma}=(\mathbf{u}_m-\bar{\mathbf{u}})\cdot\mathbf{p}_a,\qquad
m_2^{(a)}=\frac1n\sum_m q_{ma}^2,\qquad m_4^{(a)}=\frac1n\sum_m q_{ma}^4 .$$

**12a. The headline number — Mardia's multivariate kurtosis** (`nonGaussianity`, the
*Non-Gaussianity* summary row):

$$b_2=\frac1n\sum_m\Big(\sum_{a\in D}\frac{q_{ma}^2}{m_2^{(a)}}\Big)^{2}
=\frac1n\sum_m\big[(\mathbf{u}_m-\bar{\mathbf{u}})^{\!\top}\mathbf{S}^{-1}(\mathbf{u}_m-\bar{\mathbf{u}})\big]^2,
\qquad
\texttt{nonGaussianity}=\frac{3\,\big(b_2-d(d+2)\big)}{d(d+2)}\;\overset{d=3}{=}\;\frac{b_2-15}{5},$$

with $\mathbf{S}$ the ML ($1/n$) covariance and $D$ the $d$ defined axes (12c; $d=3$ for every
ordinary site). $b_2$ is **affine invariant** — it cannot depend on how the eigenvectors of a
(near-)isotropic site happen to be oriented — and for any elliptical distribution (Gaussian,
Student-$t$, any scale mixture of normals) $b_2=d(d+2)(1+\kappa/3)$, so the normalised value equals
the excess kurtosis $\kappa$ along **every** direction. For a non-elliptical cloud it is a
rotation-invariant summary: with independent components it is $\tfrac15\sum_a\kappa_a$ (three
fifths of the axial mean), so a feature confined to one axis (a split along $x$) is diluted by $1/5$ —
read that axis's own $\kappa$ (12b) when it is resolved. Small-sample bias: for a Gaussian
$\mathbb E[b_2]=d(d+2)(n-1)/(n+1)$, i.e. $-6/n$ in the normalised value ($-0.006$ at $n=1000$).

*Why not the mean of the per-axis kurtoses (the pre-1.0 readout).* When two or three eigenvalues
are equal up to sampling error — every cubic site (isotropic $\mathbf U$), the in-plane pair of
every uniaxial one — the eigenvectors are simply the directions of largest and smallest *sample*
variance, and those correlate with the fourth moments: PC1 leans toward the biggest outliers. On the
repository's GTS_250K Ga site 8 ($\bar{4}3m$, eigenvalue gaps 2.1% and 6.9%) the old mean was 2.50,
while the same cloud gives 1.38 on the cube axes and anything from 1.36 to 3.58 over random
orthonormal frames. `test_non_gaussianity_is_affine_invariant` (Python and JS) pins the new value to
$10^{-9}$ under an arbitrary linear map and under $10^{-7}$ perturbations that re-orient a degenerate
frame; the old readout moved by 0.33 on the same test.

**12b. Per-axis excess kurtosis** (`excessKurtosis`, the `κ` column):
$\kappa_a=m_4^{(a)}/\big(m_2^{(a)}\big)^2-3$, **shown only when the axis is resolved**
(`axisResolved`). An axis is resolved when its eigenvalue is separated from each neighbour's by more
than `AXIS_RESOLUTION_SIGMAS = 3` standard errors of the gap,

$$\big(m_2^{(a)}-m_2^{(a+1)}\big)\;>\;3\sqrt{\mathrm{SE}_a^2+\mathrm{SE}_{a+1}^2},\qquad
\mathrm{SE}_a=\sqrt{\big(m_4^{(a)}-(m_2^{(a)})^2\big)/n},$$

PC1 needing the 1–2 gap, PC3 the 2–3 gap and PC2 both. $\mathrm{SE}_a$ is the asymptotic standard
error of a sample variance along a fixed direction. The factor 3 is calibrated by simulation: for a
truly degenerate pair the statistic has median $\approx1$ (eigenvalue repulsion), a 95th percentile
of $\approx2$ and exceeds 3 in under 1% of samples at $n=216$–$8000$ for Gaussian, $t_8$ and
$\langle111\rangle$-split clouds — a "gap larger than its sampling error" (factor 1) rule would call
half of all degenerate pairs resolved. The price is power: at $n=1000$ a 17% eigenvalue gap is
resolved about a third of the time, and heavy tails (large $m_4$) widen the error further — on the
5 K GaNb₄Se₈ run, whose clouds have $\kappa\approx5$–$35$, PC2 is never resolved. The panel prints
`—` for an unresolved $\kappa_a$, dims that axis's crystal-orientation row, and adds a one-line note
naming the degenerate pair. The raw values stay in the payload with their flags.

**12c. Definedness — one rule, no floor constant, in both engines.** $\kappa_a$ (and axis $a$'s
share of $b_2$) exists only when the site has spread ($\lambda_1\ge$ `ZERO_SPREAD_VARIANCE`) and axis
$a$ has not collapsed ($\lambda_a\ge10^{-6}\lambda_1$, the `DEGENERATE_RATIO` of the degenerate flag);
otherwise $m_2^{(a)}$ is round-off, the ratio is $0/0$, and $\kappa_a$ is reported as `null`. A planar
cloud therefore has no $\kappa_3$ and its non-Gaussianity is the $d=2$ Mardia value; a zero-spread
site has neither. (Before 1.0 the engines guarded the denominator with different floors —
$\max(m_2,10^{-30})^2$ in Python, $\max(m_2^2,10^{-30})$ in JS — and printed unrelated noise, e.g.
$-0.74$ against exactly $-3.00$, for the same frozen site.)

**How to read the sign.** For the multivariate value and for a resolved axis alike:

- $\approx0$ — harmonic/Gaussian motion; the ellipsoid is a faithful description;
- $>0$ — a peaked, heavy-tailed well, or a minority off-centre component (a strongly unequal split).
  The covariance ellipsoid then *overstates* the concentrated core, which is what makes a KDE
  isosurface look tighter than its ellipsoid;
- $<0$ — flat-topped or **bimodal**. An equal-occupancy double well $N(\pm d,s^2)$ has excess
  kurtosis $-2d^4/(s^2+d^2)^2$ along the split — always negative, tending to $-2$ as $d/s$ grows,
  whether or not the two wells are resolved. A symmetric split site and a heavy-tailed single well
  therefore give **opposite** signs (`test_symmetric_split_site_is_platykurtic` pins $-0.50$ along the
  split for $d=0.15$, $s=0.08$ Å). Before 1.0 the UI tooltip said "positive = split sites", which was
  the wrong way round.

Kurtosis is computed on the raw cloud, not on the KDE, is blind to skew, and — as a fourth moment —
is dominated by the tails: a handful of outlying copies can set it.

The panel reads the **sites-table** values (full cloud); the KDE payload carries its own, computed
on the possibly-subsampled fit with the same function.

Test evidence: a Gaussian cloud gives $|\texttt{nonGaussianity}|<0.3$ (Python) / $<0.4$ (JS); a
Student-$t_3$ cloud gives $>1.0$; a synthetic Gaussian-displaced site gives $|\cdot|<0.5$; an elliptical
scale mixture recovers its analytic 1.6875 within 0.12; `axisResolved` is `[F,F,F]` for an isotropic,
`[T,F,F]` for a uniaxial and `[T,T,T]` for a triaxial Gaussian cloud (both engines).

The per-site table (including `nonGaussianity`) is published upward to the AI-assistant context
(`web_app/frontend/src/llm/context/runContext.js` → `pcaContext()`), which ranks sites by
non-Gaussianity and ships a `note` string defining the quantity so the model does not misread it.
That module is outside this engine; its note still describes the pre-1.0 mean-of-axes readout and
"split site" sign (tracked as a hand-off).

---

### Step 13 — PCA frame ↔ crystallographic frame (pointer)

[`web_app/frontend/src/pcaCrystalFrame.js`](../../web_app/frontend/src/pcaCrystalFrame.js) is a pure
module (no three.js, no DOM) holding the frame algebra: `unitCellVectors()`
($\mathsf{A}_{i\cdot}=\mathsf{L}_{i\cdot}/N_i$), `crystalPcaTransforms()`
($\texttt{fracToPca}=\mathsf{P}\mathsf{A}^{\!\top}$ and its inverse) and
`principalAxisOrientation()` (direction cosines and angles to $\mathbf{a},\mathbf{b},\mathbf{c}$, the
crystallographic direction $[uvw]$, and the dominant crystal axis). **All of that is documented in
the next section, "Principal axes in the crystallographic frame", Steps 5–7**, which owns the
crystal-frame direction material; it is not repeated here.

Two facts this page depends on:

1. **No re-referencing is ever needed.** Everything in the scene lives in one shared Cartesian basis
   — the clouds were mapped through the supercell lattice $\mathsf{L}$ (Step 1) and the unit-cell
   vectors are that same lattice divided by the supercell counts — so a dot product between a
   principal axis and a cell edge is immediately meaningful and the relations are exact.
2. **Two entry points are wired into the running app.** `unitCellVectors` is imported by
   [`useSiteCloud.js`](../../web_app/frontend/src/useSiteCloud.js), which hands the resulting
   `unitCell` to the page; and `PcaKdePage.jsx` imports `crystalOrientationRows`, which calls
   `principalAxisOrientation` → `crystalPcaTransforms`, to render the *Crystal orientation* table
   (Step 4b; next section, Step 8). The transformation matrices themselves (`fracToPca`,
   `pcaToFrac`) and the direction cosines are computed on that path but not printed. What the
   *Crystal* frame toggle changes is separate: the drawn axis triad, the shadow-box orientation and
   its re-binned wall projections (Step 10b), and the camera-snap directions — all built from
   `unitCell` and `orthonormalCrystalFrame()` inside the page. The next section, Step 10, lists
   what is displayed and what is not.

---

### Step 14 — The Site-ellipsoids (unit-cell) panel and site selection

The picker is no longer part of `PcaKdePage.jsx`. It lives in
[`components/SiteStructurePanel.jsx`](../../web_app/frontend/src/components/SiteStructurePanel.jsx), a
self-contained component that owns its own three.js scene, its own `Reset view` and PNG save menu,
and its own `a b c` gizmo toggle. `PcaKdePage.jsx` renders it as
`<SiteStructurePanel sites selectedRef onSelectSite selectedEllipsoid elementColors loadingSites />`
— it receives the sites table and a selection callback and computes nothing statistical itself. The
same component is mounted by the Displacement Directions page, which is why it was extracted. The
extraction was behaviour-preserving: the geometry below is unchanged from when it was inline.

**Which site is selected** is decided in [`useSiteCloud.js`](../../web_app/frontend/src/useSiteCloud.js),
not in either page: after a `sites` response the hook keeps the current reference number if it still
exists, otherwise prefers the first site whose `count === copiesPerCell` (a *clean* reconstructed
site, so a coordinates-only file defaults to a genuine thermal site rather than a disordered shell),
otherwise the first site in the list. `copiesPerCell` is null for files with real site columns, so
the preference is inert for them.

Geometry of the scene, all from `sites` (Step 1 + Step 3):

- **Unit-cell vectors** are recomputed inline as `lattice[i][k] / supercell[i]` — the component does
  *not* call `unitCellVectors()` and has no guard against a zero repeat count. Identical output for
  every real file; see the next section, Step 5, for the five copies of this division.
- **cell wireframe**: the 12 edges of the parallelepiped from $\mathsf{A}$-mapped corners
  $\{0,1\}^3$ (the pairs at Hamming distance 1), grey, opacity 0.4;
- **site position**: `siteFractional` mapped through $\mathsf{A}$;
- **each marker** is a unit `SphereGeometry(1, 20, 16)` transformed by
  $\mathsf{P}^{\!\top}\mathrm{diag}(\mathrm{semi})$ with
  $\mathrm{semi}_a=\max(\sigma_a\,s,\;0.15\,r_0)$, where $r_0=0.05\times$ (shortest cell edge) and the
  **global magnification** is $s=r_0/\overline{\sigma}$. Precisely, $\overline{\sigma}$ is the
  arithmetic mean of **all three RMS values of all sites pooled** (`positions.flatMap(p => p.rms)`,
  i.e. $3S$ numbers), not the selected site's own mean — a per-configuration constant, so changing
  the site set (e.g. moving the Cluster-distance slider on a reconstructed file) rescales *every*
  marker. **This is a display magnification, not a probability level** — real RMS amplitudes (~0.1 Å)
  are invisible at cell scale, so the panel shows relative size and anisotropy between sites, not
  absolute size. (Contrast with the main viewport's ellipsoid, which is a true $k(p)\sigma$ surface
  in Å.) A site lacking `rms`/`axes` falls back to a uniform sphere of radius $r_0$;
- **colour and selection state**: unselected markers take their element colour at opacity 0.45; the
  selected marker is opaque, takes the contrast highlight `SELECTED_ATOM_COLOR = 0xff7a1a`
  (with a $0.35\times$ emissive of the same colour), and is wrapped in two concentric translucent
  glow shells at $1.9\times$ and $1.35\times$ its own extent $\max(\mathrm{semi},r_0)$;
- **selected-site PC triad**: `buildAxisTriad` from
  [`sceneAxes.js`](../../web_app/frontend/src/components/sceneAxes.js), fed the site's own `axes`, rod
  length $2.3\times$ and radius $0.1\times$ that marker's extent, in the PC1/PC2/PC3 colours
  `TRIAD_COLORS` (`0xd64545`, `0x3fa34d`, `0x3f7fd6`). This panel is the one place where the PC axes and
  the a/b/c rods are visible together — the main viewport's *Frame* switch is exclusive;
- **bonds**: straight lines between in-cell sites only (no periodic images), with cutoff
  $\min(3.4,\max(2.0,1.3\,d_\mathrm{nn}))$ Å where $d_\mathrm{nn}$ is the shortest site–site distance
  above 0.4 Å — a display heuristic, not a bonding analysis;
- **optional a/b/c gizmo** at the cell origin (`buildCrystalAxes` from `sceneAxes.js`, length
  $0.5\times$ shortest cell edge, radius $0.035\times$ that length), on by default, identified by the
  gold/teal/orchid `CELL_AXIS_COLORS` and the legend under the canvas. The rods are **normalised**, so
  they carry the cell's orientation only, never the relative lengths of $a$, $b$, $c$; a cell vector
  with $\lVert\mathbf{a}_j\rVert^2<10^{-12}$ is skipped, silently drawing two rods while the legend
  still shows three.

**Picking.** A `Raycaster` over `sitesGroup` on `pointerup`, gated by a 5-pixel pointer-travel budget
so an orbit drag is not read as a click; the hit mesh carries `userData.referenceNumber`, which is
passed to `onSelectSite`. Hovering swaps the cursor between `pointer` and `grab`. The camera frames
the cell **once** (radius $1.7\times$ the longest cell edge) and then keeps the user's orientation
across rebuilds; `Reset view` restores that framing.

**The dropdown picker** in the controls bar (still in `PcaKdePage.jsx`) is the other way to select a
site. Each option reads ``#{referenceNumber} {label} — U={uIso} Å²`` (4 decimals; `label` is the
element, or the composition such as `Ga0.75In0.25` for a mixed site — `siteLabel()`), with
`` ({count}/{copiesPerCell})`` appended only when `copiesPerCell` is set — i.e. only for
browser-reconstructed, coordinates-only files. The matching tag renders **not** next to the picker
but at the top of the *Displacement statistics* summary column (`pca-site-tag`, in the
`pca-stats-col--summary` block), reading `{count}/{copiesPerCell} copies · <status>` where the status
is `clean site` ($\texttt{count}=\texttt{copiesPerCell}$), `merged / disordered` (larger) or
`fragmented — raise the distance` (smaller). The Flask `/api/pca/sites` payload carries no
`copiesPerCell`, so neither the suffix nor the tag ever renders in server mode (Step 1b). No
direction information appears in the picker.

---

### Parameters and defaults

| Parameter | UI control | UI default | Engine default | Range / values | Units |
|---|---|---|---|---|---|
| `grid` $G$ | Grid select | **40** | 48 | 24/32/40/48/56/64 in UI; clamped to [8, 128] in both engines | points per axis |
| `bw` | Bandwidth select | `scott` | `scott` | `scott`, `silverman`, or a positive float (API/engine only) | — |
| `bwScale` | *not exposed* | 1.0 | 1.0 | $>0$ | — |
| `extent` | Box slider | **4.0** | 3.0 | 2–5, step 0.5 | broadened σ |
| `cubicBox` | *forced* | **true** | `False` | display box only (`boxHalfWidths`); the volume is always per-axis | — |
| `probability` $p$ | Level slider (ellipsoid) | **0.5** | 0.5 | 0.10–0.99, step 0.01 | — |
| `isoPercent` | Level slider (isosurface) | **50** (= ellipsoid default) | — | 1–99, step 1 | % of captured mass |
| `projections` | Projections toggle | on | `True` | — | — |
| `clusterThreshold` | Cluster slider (reconstructed files only) | **1.5** | 1.5 (`DEFAULT_CLUSTER_THRESHOLD`) | 0.4–2.5, step 0.1 | Å |
| `shellContrast` | Contrast slider | 1.0 | — | 0.5–3.0, step 0.1 | gain |
| colormap / shell colormap | selects | `viridis` | — | `COLORMAP_NAMES` | — |
| `MAX_PCA_FIT_POINTS` | — | — | **20 000** | — | points |
| `EIGENVALUE_FLOOR_RATIO` | — | — | **1e-8** | — | fraction of $\lambda_1$ |
| `DEGENERATE_RATIO` | — | — | **1e-6** | — | $\lambda_3/\max(\lambda_1,10^{-30})$ |
| `ZERO_SPREAD_VARIANCE` | — | — | **1e-8** (both engines) | — | Å² on $\lambda_1$; below it `zeroSpread` |
| kurtosis definedness | — | — | $\lambda_1\ge10^{-8}$ Å² and $\lambda_a\ge10^{-6}\lambda_1$, else `null` (both engines) | — | — |
| `AXIS_RESOLUTION_SIGMAS` | — | — | **3** (both engines) | — | standard errors of an eigenvalue gap for `axisResolved` |
| Jacobi sweeps / tolerance | — | — | 50 sweeps; off-diagonal sum $\le10^{-15}\lVert\mathbf U\rVert_F$ (**relative**); rotation skipped at $\lvert a_{pq}\rvert<10^{-300}$ | — | — |
| `rng_seed` | — | — | 0 | — | — |
| $k(0.5)$ | — | — | 1.5381722 | — | — |
| contour fractions | — | — | 0.10, 0.25, 0.40, 0.55, 0.70, 0.85 | — | × plane's $v_\mathrm{max}$ |
| wall / contour offsets | — | — | $-1.06\,w_t$ / $-1.05\,w_t$; wall opacity 0.96, `renderOrder = 2` on contours | — | along the third axis from `mean` |
| Site-panel marker scale | — | — | $r_0/\overline{\sigma}$, $r_0=0.05\times$ shortest cell edge, semi-axis floor $0.15\,r_0$ | — | display magnification |
| Site-panel gizmo / triad | `a b c` toggle (on) | — | gizmo $0.5\times$ shortest edge (radius $0.035\times$); PC triad $2.3\times$ marker extent (radius $0.1\times$) | — | Å |
| PC / cell axis colours | — | — | `#d64545 #3fa34d #3f7fd6` / `#e0a419 #18a3a0 #b15ad8` | `sceneAxes.js` | — |
| table precision | — | — | axes 3, λ 4, RMS 3, κ 2, $U_\mathrm{iso}$ 4, $B_\mathrm{iso}$ 3, anisotropy 2, non-Gaussianity 2, captured mass 1 | — | decimals; non-finite → `—` |

Derived quantities and their units: `covariance`, `eigenvalues`, `uIso`, `bIso` in Å²; `rms`,
`semiAxes`, `bandwidth`, `halfWidths`, `axisCoords`, `mean` in Å; `cellVolume` in Å³; `density`,
`vmin`, `vmax`, `massLevels[].level` in Å⁻³; projection `density` in Å⁻²; `mass`, `anisotropy`,
`excessKurtosis`, `nonGaussianity`, `factor` dimensionless.

---

### Caveats / what this is not

1. **The KDE is kernel-broadened; the ellipsoid is not.** The estimator's variance along each axis is
   $\lambda_a(1+f^2)$, while the drawn ellipsoid uses $\lambda_a$. With Scott's rule that is a
   $\sqrt{1+f^2}$ inflation of $+10.2\%$ at $n=216$, $+6.7\%$ at $n=1000$, $+3.8\%$ at $n=8000$. For a
   *perfectly Gaussian* cloud the $p\%$ mass isosurface therefore sits **outside** the $p\%$ ellipsoid
   by exactly that factor (`KernelBroadeningGuidanceTests` pins it). A surface that hugs or falls
   inside the ellipsoid is the anharmonic signal; a surface slightly outside it may be nothing but the
   bandwidth. Neither engine corrects for this; the isosurface tooltip states it with the current
   factor ($\sqrt{1+f^2}$ from the loaded volume), and the two levels default to the same 50%.
2. **A KDE is a smoother, not a model.** It has no physics in it: no harmonic approximation, no
   temperature, no separation of static from dynamic disorder. Everything shown is the RMC
   configuration's *static snapshot* of positions.
3. **Round-off equality is to the estimator, not to truth.** The separable volume equals
   `scipy.stats.gaussian_kde` (Python) and a brute-force multivariate Gaussian KDE (JS) to $10^{-9}$.
   Whether that estimator is the right description of the cloud is a separate question, and Scott's
   rule is a rule of thumb calibrated for Gaussian data — it over-smooths multimodal (split-site)
   clouds.
4. **The two engines are not proven equal to each other.** Each is pinned against its own reference.
   There is no golden-file parity test between `pca_kde.py` and `pcaKde.js`. The known systematic
   difference is the different subsample draw above 20 000 points (Step 2). Eigenvector signs follow the same canonicalisation rule
   in both engines (with one unreachable difference in how the handedness flip is written — next
   section, Step 3), so they do not differ.
5. **Server mode cannot read coordinates-only `.rmc6f` files.** The fold-and-cluster site
   reconstruction exists only in the browser engine. Reconstructed sites are a *heuristic* grouping
   controlled by a distance knob; a `162/27` count means the "ellipsoid" is really a disordered shell
   and should be read from the KDE, not as an ADP.
6. **`U` is Cartesian.** No $U_\mathrm{cif}$/$U^{ij}$/$\beta_{ij}$ conversion is performed anywhere;
   the tabulated tensor is not directly comparable to CIF ADP components unless the cell is
   orthogonal and axis-aligned. Convert with
   $\mathbf{U}_\mathrm{cif}=\mathsf{D}^{-1}\mathsf{A}^{-\top}\mathbf{U}_\mathrm{cart}\mathsf{A}^{-1}\mathsf{D}^{-1}$
   (Step 4b).
7. **Non-Gaussianity is Mardia's multivariate kurtosis**, normalised to the marginal excess
   kurtosis of an elliptical distribution (Step 12). It is frame-independent but carries no
   information about skewness or the direction of a feature, and a one-axis feature is diluted by
   $1/5$. Its sign separates a heavy-tailed well ($>0$) from a symmetric split site ($<0$); per-axis
   $\kappa$ is shown only for resolved axes.
8. **Element-pooled clouds mix orientations.** The engine supports `element=`, and pooling is valid
   because each site is pre-centred, but the resulting single ellipsoid is only meaningful when the
   pooled sites are symmetry-equivalent *and* similarly oriented.
9. **Mass fractions are box-relative and quadrature-limited.** `mass` reports how much of the unit
   total the box captured (≈0.992 at `extent=3` for a Gaussian); isosurface
   levels enclose $p$ of *that* captured mass, from a rectangle-rule sum on the grid.
10. **Wall projections are per-plane normalised** (each texture is `density/vmax` of its own plane)
    and their contour levels are fixed fractions of that plane's peak — they are **not** enclosed-mass
    contours and are not comparable between walls. The crystal-frame walls are additionally a
    bilinear re-binning of the sampled volume that discards density outside the cube half-width, so
    they do not carry the same total as the PC-frame analytic marginals.
11. **The drawn ellipsoid mixes two computations.** `semiAxes` come from the sites table (full cloud,
    unfloored eigenvalues) while the orientation and centre come from the KDE result (possibly
    subsampled, floored eigenvalues). They coincide whenever the cloud has at most 20 000 copies, and
    diverge slightly above the cap — for a pooled cloud, or for a single site in a box of $\ge28$
    cells per edge (the UI's only request is per-site, so the second case is reachable from the
    page).
12. **Marching-cubes normals are recomputed by the renderer**, discarding the field-gradient normals
    the module produces; the 4-crossing (saddle) case in the 2D contour code is emitted without
    disambiguation. Both affect appearance only, not the reported numbers.
13. **The browser parse cache key is a sampled signature** (FNV-1a over every 64th character plus the
    length), not a full hash — a theoretical collision would reuse the wrong parsed configuration.
    The cache holds one slot and lives in a worker that is never torn down, so it is shared with the
    Displacement Directions page for the lifetime of the tab.
14. **The Site-ellipsoids markers are magnified by a per-configuration constant**, not drawn at a
    probability level. Shapes and orientations are exact; absolute sizes are not, and changing the
    site set rescales every marker at once (Step 14).
15. **Directions appear only in the statistics panel.** The site picker carries no direction
    information; the panel prints the Cartesian axis components and, in the *Crystal orientation*
    column, the angles to $a$/$b$/$c$ and $[u\,v\,w]$ of the sign-canonical representative
    (Step 4b). An unresolved axis's orientation is sampling noise and is dimmed.


## Principal axes in the crystallographic frame

### What this shows

This section continues the PCA Ellipsoid page. The previous section built, for each
crystallographic site, a **displacement cloud** — one Cartesian offset per supercell copy of the
site — reduced it to a covariance tensor $C\equiv\mathbf{U}$, diagonalised it, and used the result
to draw thermal ellipsoids and a separable 3D KDE. This section takes the *same* eigenframe and
asks the geometric question the previous one deliberately left alone: **what does a principal axis
mean as a direction, and how does it sit relative to the crystallographic $a$, $b$, $c$ axes?**

That covers: the sign and handedness canonicalisation and why it does *not* make an axis an arrow;
the unit-cell vectors in the shared Cartesian basis; the fractional ↔ PCA transform matrices; the
direction cosines and angles to $a$, $b$, $c$; the crystallographic direction $[u\,v\,w]$; the
dominant-axis rule; the crystal-frame display frame used by the viewport's shadow box; what the app
does and does not surface of all this; and one unrelated site-displacement scalar, `dispA`, that is
easy to confuse with the PCA amplitudes.

Everything upstream of the eigenframe — parsing the `.rmc6f`, building and folding the clouds, the
covariance, the eigensolvers, the ADP readouts, the ellipsoid scaling, the bandwidth, the KDE, the
isosurface, the wall projections, the shell, the non-Gaussianity and the site picker — belongs to
the previous section, "PCA Ellipsoid page — displacement clouds, ADP tensors, and the separable 3D
KDE", and is cross-referenced rather than repeated. Its "Where the code lives" table names the
current owner of each module, including the shared hook
[`useSiteCloud.js`](../../web_app/frontend/src/useSiteCloud.js) that routes every request and the
`sites`/`kde` payload shapes used below.

> **Not the Displacement Directions page.** The app has a nav tab of that name (renamed from
> "Orientation" in commit `bf039be`); it is documented in the sibling
> [Displacement Directions](displacement-directions.md) reference. That page asks a *different*
> question of the same clouds — the
> solid-angle distribution of the individual displacement vectors, hex-binned on a Goldberg sphere,
> served by `/api/pca/orientation` and the worker's `{kind:'orientation'}`. Nothing in this section
> describes that page; the two share only the site-cloud plumbing and the site picker.

The crystal-frame geometry itself — cell vectors, frame transforms, direction cosines — lives in a
pure-JavaScript module,
[`web_app/frontend/src/pcaCrystalFrame.js`](../../web_app/frontend/src/pcaCrystalFrame.js), with no
Python counterpart. The page's *Crystal orientation* table uses it through `crystalOrientationRows`
(Step 8); what is displayed and what is only computed is listed in
[Step 10](#step-10--what-is-displayed-and-what-is-computed-but-not).

#### Notation used in this section

| Symbol | Meaning | Units |
| --- | --- | --- |
| $N$ | number of box copies of the site (atoms in the cloud) | — |
| $\mathbf{u}_m$ | Cartesian displacement vector of copy $m$ (see the note below) | Å |
| $\bar{\mathbf u}$ | cloud mean (`mean`) | Å |
| $C$ | $3\times3$ cloud covariance (`covariance`), the ADP tensor $U$ | Å² |
| $\hat{\mathbf e}_i$ | $i$-th **principal** axis, a unit Cartesian vector; row $i$ of `axes`; $i=1,2,3$ | — |
| $P$ | $3\times3$ matrix whose **rows** are $\hat{\mathbf e}_1,\hat{\mathbf e}_2,\hat{\mathbf e}_3$ | — |
| $\lambda_i$ | $i$-th eigenvalue of $C$ (`eigenvalues`), the variance along $\hat{\mathbf e}_i$ | Å² |
| $\sigma_i=\sqrt{\lambda_i}$ | RMS amplitude along $\hat{\mathbf e}_i$ (`rms`) | Å |
| $\kappa_i$ | excess kurtosis along $\hat{\mathbf e}_i$ (`excessKurtosis`) | — |
| $\mathbf L_i$ | supercell lattice vector $i$ (`latticeVectors` row $i$) | Å |
| $n_i$ | supercell repeat count along axis $i$ (`supercell`) | — |
| $\mathbf a_1,\mathbf a_2,\mathbf a_3$ | unit-cell vectors $a,b,c$; rows of $M$ | Å |
| $\ell_i=\lVert\mathbf a_i\rVert$ | unit-cell edge length (`cellLengths`, `cellEdgeA`) | Å |
| $M$ | unit-cell matrix, rows $a,b,c$, in the shared Cartesian basis | Å |
| $\mathbf r$ | an arbitrary Cartesian vector (a displacement, not a position) | Å |
| $\mathbf f$ | a displacement in **unit-cell fractional** coordinates | — |
| $\mathbf q$ | a displacement in the PCA frame | Å |
| $[u\,v\,w]$ | crystallographic direction indices (`crystalDirection`) | — |
| $\theta_{ij}$ | angle between $\hat{\mathbf e}_i$ and the cell edge $\mathbf a_j$ | degrees |
| $\mathbf e_0,\mathbf e_1,\mathbf e_2$ | orthonormal **display** frame built from the cell (Step 9b), **0-indexed** — *not* the principal axes $\hat{\mathbf e}_i$ | — |
| $L$ | drawn rod length of an axis triad (Step 9a) | Å |

> **Symbol collision, stated once.** $\hat{\mathbf e}_1,\hat{\mathbf e}_2,\hat{\mathbf e}_3$
> (hatted, 1-indexed) are the **principal axes of the displacement cloud**.
> $\mathbf e_0,\mathbf e_1,\mathbf e_2$ (unhatted, 0-indexed) are the **orthonormal display
> frame derived from the unit cell** in Step 9b. "The $\mathbf e_0$–$\mathbf e_1$ face" is a
> crystal-frame face, never a PC plane.

> **What $\mathbf u_m$ is, exactly.** $\mathbf u_m$ is the Cartesian offset of copy $m$ from its
> own box copy's origin, already centred on the site mean $\bar{\mathbf u}$. Its construction —
> the supercell fold, the per-site mean subtraction, the mapping through $\mathsf L$, and the
> fold-and-cluster fallback for coordinates-only files — is the previous section's Step 1, and
> which engine subtracts the mean where is its Step 3. Formulas here keep the explicit
> $-\bar{\mathbf u}$; in Python that term is identically zero by construction.

Everything below lives in **one shared orthonormal Cartesian basis in Ångström**: the basis in
which the `.rmc6f` lattice vectors are written. The displacement clouds are mapped into it
through the supercell lattice, and the unit-cell vectors are derived from the same lattice. That
shared basis is the reason no re-referencing is ever needed between the PCA axes and the
crystal axes.

---

### Steps 1–2 — The eigenframe this section starts from (owned by the PCA Ellipsoid section)

The cloud covariance and its eigen-decomposition are derived in the previous section, Step 3, and
the readouts built from them in its Step 4. In brief, so the algebra below has its symbols:

$$C=\frac{1}{\max(N-1,1)}\sum_{m}\bigl(\mathbf u_m-\bar{\mathbf u}\bigr)\bigl(\mathbf u_m-\bar{\mathbf u}\bigr)^{\mathsf T}
\ \ [\text{Å}^2],\qquad
C\,\hat{\mathbf e}_i=\lambda_i\hat{\mathbf e}_i,\qquad \lambda_1\ge\lambda_2\ge\lambda_3\ge0 .$$

$C$ is symmetric, positive semi-definite, and **is** the anisotropic displacement tensor $U$ for
that site. `axes` returns the eigenvectors as **rows** ($P$); `eigenvalues` are clamped at 0 from
below; `rms` is $\sigma_i=\sqrt{\lambda_i}$ in Å. Solvers, flooring constants, the degeneracy flag,
the subsample cap and the two engines' test coverage: previous section, Steps 2–4 and 6.

One convention matters for everything that follows: **the ordering is by variance, descending**.
PC1 is the softest (largest-amplitude) direction, PC3 the stiffest. Nothing else orders them — PC1
is *not* "the $a$-like axis", and the index carries no crystallographic meaning at all.

The other thing worth carrying forward is a negative result about the tests. Each engine's KDE is
pinned against its own reference (SciPy for Python, an inline brute-force multivariate KDE for the
worker), and **neither is compared to the other**: there is no cross-engine assertion on
eigenvectors, eigenvalues or densities anywhere in the suite. On a near-degenerate cloud two
different — and equally valid — eigenframes would both pass. The **"PCA-KDE is separable, not
approximate"** bullet in [AGENTS.md](../../AGENTS.md) — "both assert exact agreement against the full
estimator" — reads as a cross-engine claim and is not one. So a *direction* read off one runtime is
not guaranteed to match the other's beyond what the shared canonicalisation of Step 3 enforces.

### Step 3 — Sign and handedness canonicalisation, and what a direction means

An eigenvector is defined only up to sign, and LAPACK's choice depends on the build. Both
engines apply the same two-step canonicalisation (`_canonical_axes` in Python,
`eigenDecomposition` in JS) so a rerun on another machine prints the same table:

1. For each axis, find the component of largest magnitude and flip the whole axis if that
   component is negative. (Ties/zeros default to $+1$.)
2. Compute $\det P$; if it is negative, negate $\hat{\mathbf e}_3$, leaving a right-handed frame
   $\hat{\mathbf e}_1\times\hat{\mathbf e}_2=\hat{\mathbf e}_3$.

**The order matters and step 2 wins.** Negating $\hat{\mathbf e}_3$ for handedness makes its
largest-magnitude component negative again. So the sign rule of step 1 is guaranteed only for
$\hat{\mathbf e}_1$ and $\hat{\mathbf e}_2$; the printed PC3 row violates it whenever the
pre-flip frame was left-handed — roughly half the time.

> **One unreachable divergence between the engines.** Python implements the flip as
> $\hat{\mathbf e}_3 \mathrel{*}= \operatorname{sign}(\det P)$, so $\det P = 0$ would *zero* the
> third axis; JS tests `if (det < 0)` and leaves it untouched. Orthonormal eigenvectors give
> $\lvert\det P\rvert = 1$, so this case cannot arise — but the two are not literally the same
> operation.

**Read the consequence carefully.** The canonicalisation makes the *printed numbers*
reproducible; it does not make the *physical direction* signed. A principal axis is an
undirected line: the cloud's covariance is invariant under $\mathbf u\mapsto-\mathbf u$, so
$\hat{\mathbf e}_i$ and $-\hat{\mathbf e}_i$ describe the same displacement direction. A table
row of `(0.577, 0.577, 0.577)` means the atoms move along $\pm$ that Cartesian line, not toward
$+x+y+z$ specifically. Any physical statement drawn from a principal axis must be invariant
under that sign flip.

A second limitation is intrinsic: when two eigenvalues are (nearly) equal, the axes spanning
that degenerate subspace are arbitrary — any rotation within the subspace is an equally valid
eigenbasis. No sign convention can repair that; for a near-isotropic site, PC1 and PC2
individually mean nothing and only the plane they span does.

**Test evidence.** `test_axes_are_orthonormal_and_right_handed` asserts
$PP^{\mathsf T}=\mathsf I$ to `atol=1e-9` and $\det P = 1$ to 9 places, i.e. it pins the *frame*,
not the individual axes: an arbitrary rotation inside a degenerate subspace passes it unchanged.

### Step 4 — Per-axis excess kurtosis $\kappa_i$ (pointer)

Each axis also carries a shape number, the excess kurtosis $\kappa_i$ of the cloud's projection
onto $\hat{\mathbf e}_i$ — reported as `excessKurtosis[i]` with a resolution flag
`axisResolved[i]`. It is **derived in the previous section, Step 12**, including the population
($N$) moments, the definedness rule, the resolution test and how to read the sign; the headline
`nonGaussianity` is Mardia's rotation-invariant kurtosis, not an average of the $\kappa_i$.

The things to carry into a *direction* statement: $\kappa_i$ means something only when
`axisResolved[i]` is true — inside a degenerate eigenvalue pair the axis itself is sampling noise.
A resolved $\kappa_i$ is a **marginal** number attached to one axis: $\kappa_1>0$ says the
distribution along PC1 is peaked and heavy-tailed, $\kappa_1<0$ that it is flat-topped or bimodal
(a symmetric split along that line); it says nothing about the sign of the displacement (Step 3).
Combine it with the KDE isosurface before attributing a mechanism to a direction.

### Step 5 — Unit-cell vectors in the shared Cartesian basis

**Inputs:** `latticeVectors` — the three **supercell** vectors $\mathbf L_i$ (Å, rows) from the
`.rmc6f` `Lattice vectors` block — and `supercell`, the repeat counts $n_i$ from the
`Supercell dimensions` line. Both are returned by `/api/pca/sites` and by the worker's `sites`
response, so the frontend has them in either runtime.

**Operation:** divide out the repeats,

$$\mathbf a_i=\frac{\mathbf L_i}{n_i},\qquad i=1,2,3\ \ (\mathbf a_1=a,\ \mathbf a_2=b,\ \mathbf a_3=c),$$

producing $M$, the $3\times3$ unit-cell matrix with the cell edges as **rows**, in Å.

**Guards:** returns `null` unless both arguments are length-3 arrays; a zero repeat count falls
back to $n_i=1$ (the test `if (Math.abs(supercell[i]) > 0)` guards only against zero — a negative
count would be used as-is and would flip the vector).

**Why this makes the PCA axes and cell vectors directly comparable.** The displacement clouds are
built as fractional offsets within the supercell and mapped to Cartesian Å through the *same*
$\mathbf L_i$ (previous section, Step 1). Dividing that lattice by the repeat counts
therefore yields cell edges expressed in the identical Cartesian basis as the eigenvectors of
$C$. No change of basis, no origin shift, no re-referencing: a dot product between
$\hat{\mathbf e}_i$ and $\mathbf a_j$ is immediately meaningful.

**Code:** `unitCellVectors` in [`pcaCrystalFrame.js`](../../web_app/frontend/src/pcaCrystalFrame.js);
test `unitCellVectors` in
[`pcaCrystalFrame.test.js`](../../web_app/frontend/src/__tests__/pcaCrystalFrame.test.js).

> **Five implementations of the same division exist in the frontend.**
>
> 1. `unitCellVectors` in [`pcaCrystalFrame.js`](../../web_app/frontend/src/pcaCrystalFrame.js)
>    line 95 — guard `Math.abs(n) > 0`.
> 2. `conventionalCell` in [`symmetryModel.js`](../../web_app/frontend/src/symmetryModel.js)
>    line 15 — guard `Math.max(n, 1)`.
> 3. `conventionalCell` in
>    [`llm/context/pairCorrelations.js`](../../web_app/frontend/src/llm/context/pairCorrelations.js)
>    line 19 — same guard, and a *deliberate* re-implementation: its comment records that "the llm
>    module must not import host modules" (the import boundary in `src/llm/README.md`), so this
>    copy is a design decision rather than drift.
> 4. `unitVec` inside `sitesByClustering` in
>    [`workers/pcaKde.js`](../../web_app/frontend/src/workers/pcaKde.js) line 617 — the PCA engine
>    itself recomputes
>    `latticeVectors.map((row, i) => row.map((value) => value / Math.max(supercell[i], 1)))` for the
>    clustering metric (guard `Math.max(n, 1)`).
> 5. An inline `lattice.map((row, i) => row.map((value) => value / supercell[i]))` at
>    [`SiteStructurePanel.jsx`](../../web_app/frontend/src/components/SiteStructurePanel.jsx)
>    line 170 — the component that renders the Site-ellipsoids scene — with **no guard at all**.
>
> They agree for every real file; the guards differ only on malformed metadata. A sixth, *scalar*
> variant is not counted here because it computes a different quantity: `cellEdgeA` in
> [`browserData.js`](../../web_app/frontend/src/browserData.js) line 516 divides each row's
> **norm** by $n_i$ (guard `Math.max(n, 1)`), giving cell-edge lengths $|\mathbf a_i|$ rather than
> the vectors. Note the `SiteStructurePanel` copy moved with the panel when it was extracted from
> `PcaKdePage.jsx`; it was not consolidated onto `unitCellVectors`, even though the hook now hands
> the same component's host page a computed `unitCell`.

### Step 6 — Frame transforms between fractional and PCA coordinates

**Inputs:** $P$ (PCA axes as rows, dimensionless) and $M$ (cell vectors as rows, Å).

**Conventions, stated unambiguously.** Both matrices are arrays of rows, `m[row][col]`. Because
the *rows* of $M$ are the cell vectors, a fractional column vector $\mathbf f$ maps to Cartesian
as

$$\mathbf r=M^{\mathsf T}\mathbf f\qquad\Longleftrightarrow\qquad r_k=\sum_i f_i\,M_{ik},$$

i.e. $\mathbf r = f_1 a + f_2 b + f_3 c$. Its inverse is $\mathbf f = M^{-\mathsf T}\mathbf r$.
Similarly, because the rows of $P$ are the (orthonormal) principal axes, $P^{-1}=P^{\mathsf T}$
and a Cartesian vector maps to PCA coordinates as $\mathbf q = P\,\mathbf r$, with
$\mathbf r = P^{\mathsf T}\mathbf q$.

**Operation:** compose the two,

$$\texttt{fracToPca}=P\,M^{\mathsf T},\qquad
\texttt{pcaToFrac}=\bigl(P\,M^{\mathsf T}\bigr)^{-1}=M^{-\mathsf T}P^{\mathsf T}.$$

**Units.** `fracToPca` carries **Å per unit-cell fraction**: feed it a fractional *displacement*
(dimensionless, referred to the **unit** cell, not the supercell) and it returns that
displacement's components along PC1/PC2/PC3 in Å. `pcaToFrac` is its inverse, in unit-cell
fractions per Å.

**Return shape.** `crystalPcaTransforms` returns **three** keys —
`{ fracToPca, pcaToFrac, unitCell }` — echoing the cell matrix back alongside the two transforms.

**Null return, and the hole in the guard.** The function's only validation is

```js
if (!Array.isArray(pcaAxes) || pcaAxes.length !== 3 || !unitCell) return null;
```

plus a `null` from `invert3`. So it returns `null` when `pcaAxes` is not a length-3 array, when
`unitCell` is falsy, or when $P M^{\mathsf T}$ is singular — `invert3` declares singularity at
$\lvert\det\rvert<10^{-12}$ (or a non-finite determinant) and returns `null` rather than emitting
`NaN`s, so a collapsed cell degrades gracefully. It does **not** check that `unitCell` is a
3-element array of 3-element rows, nor that the rows of `pcaAxes` have length 3: a truthy but
malformed `unitCell` (wrong length, short rows, a plain object) **throws** inside
`transpose3`/`multiply3` instead of returning `null`. Since $P$ is orthonormal
($\lvert\det P\rvert=1$), a legitimate singularity here means a degenerate **cell**, not a
degenerate cloud.

**Code:** `crystalPcaTransforms` in
[`pcaCrystalFrame.js`](../../web_app/frontend/src/pcaCrystalFrame.js), built on the local
`transpose3` / `multiply3` / `matVec3` / `determinant3` / `invert3` helpers (adjugate-over-
determinant inverse). Tests cover the cubic case (where $P M^{\mathsf T}=8P$ for an 8 Å cell),
an oblique monoclinic cell ($\beta=105°$) round-tripping to the identity, and the collapsed-cell
`null`.

### Step 7 — Orientation of each principal axis against $a$, $b$, $c$

`principalAxisOrientation(pcaAxes, unitCell)` builds the transforms of Step 6 **first** and
returns `null` if they fail:

```js
const transforms = crystalPcaTransforms(pcaAxes, unitCell);
if (!transforms) return null;
```

That gating is worth stating, because the direction cosines and angles of 7a need no matrix
inversion at all and stay perfectly well defined for a collapsed or malformed cell — yet they
are discarded along with $[u\,v\,w]$. A caller handing in a degenerate cell gets nothing rather
than the angles. The function then inverts $M$ a **second** time, independently of `pcaToFrac`
(see 7b).

It returns `{ perAxis, cellLengths, fracToPca, pcaToFrac }`, with `perAxis[i]` reporting three
things for each principal axis.

#### 7a. Direction cosines and angles

$$\cos\theta_{ij}=\hat{\mathbf e}_i\cdot\hat{\mathbf a}_j,\qquad
\hat{\mathbf a}_j=\frac{\mathbf a_j}{\lVert\mathbf a_j\rVert},\qquad
\theta_{ij}=\frac{180}{\pi}\arccos\bigl(\mathrm{clamp}(\cos\theta_{ij},-1,1)\bigr).$$

`cellLengths` — the un-normalised $\ell_j=\lVert\mathbf a_j\rVert$ in Å — is returned alongside.

Two preconditions the code assumes rather than enforces:

- **The cell edges are normalised, with a fallback.** `cellUnit` divides by
  `cellLengths[i] || 1`, so a **zero-length** cell edge is left as the zero vector and silently
  reports $\cos\theta_{ij}=0$, $\theta_{ij}=90°$ for that edge, rather than failing.
- **The PCA axes are taken on trust.** `principalAxisOrientation` never normalises `pcaAxes`; it
  dots the caller's rows as-is. Feed it non-unit rows and `dot3(axis, u)` is not a cosine, and
  the $[-1,1]$ clamp then *hides* the error by pinning the result at $0°$ or $180°$. The clamp is
  a round-off guard only **given** unit input. Both engines do supply unit rows.

These angles are always well defined, including for a triclinic cell, because they measure the
angle between two ordinary Cartesian directions. What is *not* generally true for an oblique
cell is $\sum_j\cos^2\theta_{ij}=1$; that identity needs the $\hat{\mathbf a}_j$ to be mutually
orthogonal.

#### 7b. Crystallographic direction $[u\,v\,w]$

$$\begin{bmatrix}u\\v\\w\end{bmatrix}_i \;\propto\; M^{-\mathsf T}\,\hat{\mathbf e}_i,
\qquad\text{normalised so } \max_k\bigl\lvert[u\,v\,w]_{i,k}\bigr\rvert = 1 .$$

**Why $M^{-\mathsf T}$ and not $M^{-1}$.** Step 6 fixed $\mathbf r=M^{\mathsf T}\mathbf f$, so the
fractional components of a Cartesian vector are $\mathbf f = M^{-\mathsf T}\mathbf r$ — full
stop. Concretely, the **rows** of $M^{-\mathsf T}$ are the reciprocal-lattice vectors
$a^\ast,b^\ast,c^\ast$ defined by $\mathbf a_i\cdot\mathbf a^\ast_j=\delta_{ij}$ (no $2\pi$),
so $M^{-\mathsf T}\hat{\mathbf e} = (a^\ast\!\cdot\hat{\mathbf e},\,b^\ast\!\cdot\hat{\mathbf e},\,c^\ast\!\cdot\hat{\mathbf e})$
— exactly the coefficients of $\hat{\mathbf e}$ in the $a,b,c$ basis. Contracting with $M^{-1}$
instead would dot $\hat{\mathbf e}$ against the *columns* of $M^{-\mathsf T}$, which are not
lattice vectors of anything; contracting with $M$ itself would give the covariant components
$(a\!\cdot\hat{\mathbf e},\,b\!\cdot\hat{\mathbf e},\,c\!\cdot\hat{\mathbf e})$, which are
proportional to the Miller indices $(hkl)$ of the plane *perpendicular* to $\hat{\mathbf e}$ —
a different object.

The three contractions are **not** simply "the same for a cubic cell": $[u\,v\,w]$ and $(hkl)$
point the same way only when the cell metric is isotropic ($a=b=c$, all angles $90°$), while
$M^{-1}$ and $M^{-\mathsf T}$ additionally coincide whenever $M$ happens to be **symmetric** — as
it is for the axis-aligned orthorhombic cells `.rmc6f` files normally carry,
$M=\mathrm{diag}(a,b,c)$. Conversely, for a cubic cell written in a rotated orientation
($M=\lambda R$ with $R$ a non-symmetric rotation) $M\hat{\mathbf e}$ and
$M^{-\mathsf T}\hat{\mathbf e}$ are both $\propto R\hat{\mathbf e}$ while
$M^{-1}\hat{\mathbf e}\propto R^{\mathsf T}\hat{\mathbf e}$ is not.

**What the code does, as opposed to what is equivalent.** Analytically the normalised
$[u\,v\,w]_i$ is column $i$ of `pcaToFrac`, rescaled. The implementation does **not** read
`pcaToFrac`. It performs a second, independent inversion and transposes it,

```js
const invT = transpose3(invert3(unitCell));   // M⁻ᵀ
const frac = invT ? matVec3(invT, axis) : [0, 0, 0];
```

so the numbers come from a different determinant/adjugate evaluation and agree with `pcaToFrac`
only to round-off. The apparent `invT ? … : [0,0,0]` guard is **dead code**: a `null` from
`invert3` would already have thrown a `TypeError` inside `transpose3` one line earlier.

**What $[u\,v\,w]$ means physically:** the axis is parallel to the lattice vector
$u\,a+v\,b+w\,c$. Nothing more. In particular:

- It is **not** a unit vector. Its Cartesian length is $\lVert M^{\mathsf T}[u\,v\,w]^{\mathsf T}\rVert$ Å.
- It is **not** rounded or reduced to small integers. A PC axis $30°$ from $a$ in the $a$–$b$
  plane of a cubic cell prints as $[1,\,0.577,\,0]$, not $[2\,1\,0]$ or $[\sqrt3\,1\,0]$.
- Because $\hat{\mathbf e}_i$ is sign-ambiguous (Step 3), so is $[u\,v\,w]_i$: it stands for
  $\pm[u\,v\,w]$.
- The max-magnitude normalisation makes the largest component exactly $\pm1$; the divisor falls
  back to 1 if all three components vanish.

#### 7c. Dominant (closest) crystal axis

$$\text{dominant}(i)=\arg\max_j\bigl\lvert\cos\theta_{ij}\bigr\rvert,$$

reported as `{ index, label: 'a'|'b'|'c', angleDeg }` where `angleDeg` is $\theta_{ij}$ for that
$j$. Because the **absolute value** of the cosine is compared, a direction and its negative are
treated as the same axis: an axis at $170°$ to $b$ is "closest to $b$", and
`principalAxisOrientation` then reports `angleDeg` $=170°$. The page never shows that raw value — it
goes through `crystalOrientationRows` (Step 8), which prints the $10°$ representative. Ties (two
equal $\lvert\cos\rvert$) resolve to the lowest index by the strict `>` in the scan.

**Code:** `principalAxisOrientation` in
[`pcaCrystalFrame.js`](../../web_app/frontend/src/pcaCrystalFrame.js). Tests assert $0°/90°$ and
$[1\,0\,0]$/$[0\,0\,1]$ for axis-aligned PCs in a cubic cell, $30°/60°$ for a rotation about
$z$, and — in a $\beta=110°$ monoclinic cell — that a PC axis placed exactly along the $c$ edge
recovers $[0\,0\,1]$.

#### The oblique-cell subtlety, stated plainly

For a non-orthogonal cell the two outputs answer **different** questions and must not be
conflated:

- $\theta_{ij}$ is a genuine angle in real space between two Cartesian directions. Always
  meaningful, always in $[0°,180°]$.
- $[u\,v\,w]$ is a **fractional (crystallographic) direction**. It is not a unit vector, its
  components are not angles, and the lattice direction $[u\,v\,w]$ is in general **not**
  perpendicular to the lattice plane $(u\,v\,w)$. That perpendicularity holds only in a cubic
  cell; even in an orthorhombic cell $[110]\not\perp(110)$ unless $a=b$.

---

### Step 8 — What the statistics panel prints: the Crystal orientation table

The *Displacement statistics* panel and the site picker are specified in the previous section
(Step 4b and Step 14 respectively) — columns, formulas, decimal places, the metadata line. As
**direction information** the panel carries two things:

- The *Principal axes* table prints $\hat{\mathbf e}_i$ as $x,y,z$ to 3 decimals, i.e. the axis as
  a unit vector **in the `.rmc6f` Cartesian frame**.
- The *Crystal orientation* table (`PcaKdePage.jsx`, `crystalOrientation` memo) prints, per
  principal axis, the angles ∠a/∠b/∠c to the unit-cell edges and the direction $[u\,v\,w]$ — the
  output of `crystalOrientationRows(selectedEllipsoid.axes, unitCell)`.

`crystalOrientationRows` is `principalAxisOrientation` (Step 7) folded to **one representative per
axis**. An eigenvector's sign is arbitrary (Step 3), and negating an axis sends every angle
$\theta\to180°-\theta$, every cosine to its negative and $[u\,v\,w]\to[-u\,-v\,-w]$. Of the $\pm$ pair
it returns the sense whose **closest crystal edge is acute** (`cosines[dominant] >= 0`), flipping
the cosines, angles and $[u\,v\,w]$ together so a row stays internally consistent (the `negate`
helper keeps a zero at $+0$), and recomputes `dominant.angleDeg` from the flipped angles. So the
shaded cell is always the smallest angle ($\le90°$) and $[u\,v\,w]$ points along it — how a
direction is quoted. It returns `null` when `principalAxisOrientation` does (non-array axes, a falsy
or singular cell), and the page then omits the column; it is also omitted when `unitCell` is `null`
(no lattice metadata) or `axes` is `null` (a zero-spread site).

Formatting in the table: angles `numberFormat(deg, 1)` + `°`, the closest edge's cell shaded
(`is-diagonal`); $[u\,v\,w]$ through `uvwFormat`, 2 decimals, largest component exactly $\pm1$,
$\lvert v\rvert<5\times10^{-3}$ printed as `0.00`. A row whose axis is not resolved from a neighbour
(`axisResolved`, previous section Step 12b) is dimmed and italic, with a tooltip saying its
direction is sampling noise.

- The site picker labels carry element (or composition), $U_\mathrm{iso}$ and, for reconstructed
  files, the copies-per-cell ratio — no direction content.
- The *Covariance $U$* table prints $C$ in the same Cartesian frame, with no conversion to the
  crystallographic $U^{ij}$ basis anywhere in the app (the conversion is given in the previous
  section, Step 4b).

`pcaCrystalFrame.test.js` → `crystalOrientationRows` pins the fold: an already-acute axis is left
untouched; an axis $170°$ from $a$ is reported as the $10°$ representative with every angle
$\theta\to180°-\theta$ and $[u\,v\,w]$ exactly negated (a zero component staying $+0$); in a
$\beta=110°$ monoclinic cell the dominant angle is never obtuse and equals that edge's printed angle;
and a collapsed cell returns `null`.

### Step 9 — What the UI shows in 3D: the PC ↔ Crystal frame switch

A two-button **Frame** toggle (`PC` | `Crystal`) in the viewport header drives the axis triad,
the shadow-box wall projections, the camera-snap buttons, and **Reset view** together
(`selectFrame`). The `Crystal` button is disabled when `unitCell` is `null` (no lattice/supercell
metadata). Default is `pc`. Introduced in **v0.4.0** (2026-07-16) — see
[docs/CHANGELOG.md](../CHANGELOG.md).

#### 9a. Axis rods — and which frame's vectors each uses

The main viewport draws **exactly one** triad at a time:

| Mode | Function | Vectors used | Colours | Length / radius |
| --- | --- | --- | --- | --- |
| `pc` | `buildAxisTriad(mean, kde.axes, L, 0.012·L)` | the principal axes $\hat{\mathbf e}_i$ | PC1 `#d64545` red, PC2 `#3fa34d` green, PC3 `#3f7fd6` blue | $L=\max(\texttt{halfWidths})$ |
| `crystal` | `buildCrystalAxes(mean, unitCell, L, 0.01·L)` | the **true** cell edges $\mathbf a_j$ | $a$ `#e0a419` gold, $b$ `#18a3a0` teal, $c$ `#b15ad8` orchid | $L=\max(\texttt{halfWidths})$ |

Both builders **normalise** each direction before drawing, so all three rods are the same length.
The crystal rods therefore convey the cell's *orientation* only — never the relative lengths of
$a$, $b$, $c$. Axis identity is by rod colour, matched to a legend under the canvas; there are no
3D letter labels. Rods start at the cloud mean and extend outward one half-box.

`buildCrystalAxes` ([`sceneAxes.js`](../../web_app/frontend/src/components/sceneAxes.js) line 53)
additionally **skips** any cell vector with $\lVert\mathbf a_j\rVert^2 < 10^{-12}$ (`continue`),
so a degenerate cell silently draws **two** rods instead of three while the $a/b/c$ legend below
the canvas is unchanged. `buildAxisTriad` has no such skip.

Because the switch is exclusive, **the main viewport cannot show the PC triad and the $a/b/c$
rods at the same time**, and it prints no angle between them: comparing the two frames means
flipping the switch and judging by eye. The one place both frames are visible together is the
**Site ellipsoids** panel
([`SiteStructurePanel.jsx`](../../web_app/frontend/src/components/SiteStructurePanel.jsx); geometry in
the previous section, Step 14), where the opt-in $a/b/c$ gizmo at the cell origin coexists with the
selected site's own PC triad. Both panels call the same two builders from
[`sceneAxes.js`](../../web_app/frontend/src/components/sceneAxes.js), which is why PC1/PC2/PC3 and
$a$/$b$/$c$ read as the same objects wherever they appear.

#### 9b. `orthonormalCrystalFrame` — the display frame for the shadow box

$$\mathbf e_0=\hat{a},\qquad
\mathbf e_1=\frac{b-(b\cdot\mathbf e_0)\,\mathbf e_0}{\lVert\cdot\rVert},\qquad
\mathbf e_2=\frac{\mathbf e_0\times\mathbf e_1}{\lVert\cdot\rVert}.$$

A Gram–Schmidt orthonormalisation of $(a,\,b)$, completed by a cross product: $\mathbf e_0$ along
$a$, $\mathbf e_1$ in the $a$–$b$ plane, $\mathbf e_2=\mathbf e_0\times\mathbf e_1$. **$c$ is
never referenced** — the function reads only `unitCell[0]` and `unitCell[1]`. Rows are the frame
axes, matching the `axes` row convention. For a **right-handed** orthogonal cell this is exactly
$(\hat a,\hat b,\hat c)$; for a left-handed one $\mathbf e_2=-\hat c$.

**No singularity guard.** The local helper is
`norm = (v) => { const n = Math.hypot(...) || 1; … }`, so it returns the **zero vector** rather
than failing. If $a\parallel b$ (or $a$ is zero) the Gram–Schmidt residual is zero, giving
$\mathbf e_1=\mathbf 0$ and hence $\mathbf e_2=\mathbf e_0\times\mathbf 0=\mathbf 0$: the crystal
shadow box, its three wall projections and the crystal default camera all collapse onto a line,
silently. The `Crystal` toggle is gated only on `unitCell !== null`, which such a cell still
satisfies.

**It is used only for the crystal-mode shadow box, its wall projections, and the crystal-mode
default camera framing.** The $a/b/c$ rods and the look-down-$a$/$b$/$c$ cameras use the true
cell edges. Consequence for an oblique cell, which the code comments do not spell out: of the
three box faces, only the $\mathbf e_0$–$\mathbf e_1$ plane is a genuine crystallographic plane
(it spans $a$ and $b$). The $\mathbf e_0$–$\mathbf e_2$ face contains $a$ but is perpendicular to
the $a$–$b$ plane, so it is not the $a$–$c$ plane unless $\beta=90°$; the
$\mathbf e_1$–$\mathbf e_2$ face is the plane perpendicular to $a$, not the $b$–$c$ plane.

#### 9c. What gets drawn on the frame's three faces (pointer)

The density painted on the crystal-mode shadow box is the *same* KDE volume as in PC mode,
re-binned onto the three faces of the frame above by `projectDensityOntoFrame` — a bilinear splat
of each grid node's mass into 2D accumulators, divided by the bin area. The algorithm, its
node-sampled (not cell-centred) quadrature, the density it silently drops outside the cube
half-width, the per-plane `vmax` normalisation, the wall/contour constants and the main-thread cost
are all in the previous section, Step 10b.

Two things belong here because they are statements about the **frame**, not about the density:

- the returned keys are still `pc12`/`pc13`/`pc23` in crystal mode, where they actually mean the
  $\mathbf e_0$–$\mathbf e_1$, $\mathbf e_0$–$\mathbf e_2$ and $\mathbf e_1$–$\mathbf e_2$ planes —
  and by 9b only the first of those is a genuine crystallographic plane for an oblique cell;
- the box is sized from the KDE half-widths, not from the cell, so **the crystal-mode shadow box is
  not a unit cell** and its faces carry no cell dimension.

#### 9d. Cameras

All framings go through `placeMainCamera`, which flushes OrbitControls' damped orbit velocity
first (so a reframe issued during a glide lands exactly on the requested pose), then places the
camera at radius $4.3\times\max(\texttt{halfWidths})$ (falling back to 1 when that is zero) from
the cloud mean, sets near/far to radius/100 and radius$\times$100, and matches the orthographic
frustum to the perspective `CAMERA_FOV = 45°` so the projection toggle preserves framing.

- **Default / Reset view** (`frameMainCamera`) — looks down the body diagonal of the active
  frame, $\mathbf d\propto\mathbf e_0+\mathbf e_1+\mathbf e_2$, with row 2 of that frame as the
  up vector. The active frame is `kde.axes` in PC mode and `crystalDisplay.frame` (the
  orthonormal frame of 9b) in crystal mode. Applied automatically on a new site or new run, and
  on every frame switch.
- **Look down PC$n$** (`frameAlongAxis`) — view direction is $\hat{\mathbf e}_n$; the screen-up
  axis is `VIEW_UP_AXIS = [2, 2, 1]`, i.e. PC3 is up when looking down PC1 or PC2, PC2 is up
  when looking down PC3. No orthogonalisation is needed: the PC axes are orthonormal already.
- **Look down $a$/$b$/$c$** (`frameAlongCellAxis`) — view direction is the **true** normalised
  cell edge $\hat{\mathbf a}_n$. The up vector starts as a neighbouring edge
  (`CELL_VIEW_UP_AXIS = [2, 2, 1]`: $c$ up for $a$ and $b$, $b$ up for $c$) and is
  **Gram–Schmidt-orthogonalised against the view direction**,
  $\mathbf{up}\leftarrow\mathbf{up}-(\mathbf{up}\cdot\mathbf d)\,\mathbf d$, then normalised.
  This is what keeps an oblique cell from producing a degenerate or tilted framing. If
  `up.lengthSq() < 1e-10` — i.e. the residual's **squared** length falls below $10^{-10}$ Å²,
  equivalently a length below $10^{-5}$ Å — the code falls back to $(0,1,0)$, or to $(1,0,0)$
  when that is within $\lvert\cos\rvert>0.99$ of the view direction, and re-orthogonalises.

The button row's label switches between `PC` (buttons `1`, `2`, `3`) and `Cell` (buttons `a`,
`b`, `c`) with the frame.

#### 9e. Degenerate clouds — what the viewport shows (pointer)

A one-copy site, or any site of an average/ideal configuration, has no spread
($\lambda_1<10^{-8}$ Å²): the panel prints *no displacement*, `—` for the non-Gaussianity, and no
axes, while the KDE request throws and the viewport shows only an error badge. The full walkthrough
is in the previous section, Step 4 and 4b.

The consequence *for directions* is the part that belongs here: the eigenvectors of a round-off
covariance are whatever the solver returns, so both engines return `axes: null` for such a site and
no crystal-orientation row is drawn. A planar or linear cloud keeps its axes (the spanned plane or
line is well defined) and is flagged by `degenerate` at $\lambda_3/\lambda_1<10^{-6}$. A near-tie
between two non-zero eigenvalues, which makes the axes inside that pair individually meaningless,
is flagged per axis by `axisResolved` (previous section, Step 12b).

### Step 10 — What is displayed, and what is computed but not

The module's application importers are
[`useSiteCloud.js`](../../web_app/frontend/src/useSiteCloud.js), which calls `unitCellVectors` inside
the `unitCell` `useMemo` as `unitCellVectors(sites.latticeVectors, sites.supercell)` (guarded by
`sites?.latticeVectors && sites?.supercell`), and `PcaKdePage.jsx`, which calls
`crystalOrientationRows` → `principalAxisOrientation` → `crystalPcaTransforms`. Because
`useSiteCloud` is the shared hook for **both** the PCA Ellipsoid page and the Displacement Directions
page, both receive the same `unitCell`.

**Displayed** (PCA Ellipsoid page, *Crystal orientation* table, Step 8): the angle $\theta_{ij}$
between each principal axis and $a$, $b$, $c$; the closest edge (shaded); the $[u\,v\,w]$ direction —
all for the acute-sense representative.

**Computed but not displayed:**

- the direction cosines themselves (only their angles are printed);
- the `fracToPca` / `pcaToFrac` matrices;
- `cellLengths` (the panel prints no $a$, $b$, $c$ magnitudes);
- the raw, un-folded `principalAxisOrientation` output (an axis at $170°$ to $b$ — the table shows
  $10°$).

The AI-assistant run context carries none of them. `pcaContext` in
[`runContext.js`](../../web_app/frontend/src/llm/context/runContext.js) emits per site only
`ref`, `element`, `U_iso_A2`, `rms_axes_A`, `anisotropy` (documented there as `rms1/rms3`),
`non_gaussianity`, and `degenerate` — **no axis directions at all**. The model therefore never
sees which way a site moves, only how much and how anisotropically.

**And it sees only the top of the list.** `pcaContext` rounds `U_iso_A2`, each entry of
`rms_axes_A`, `anisotropy` and `non_gaussianity` to **3 significant figures** (`roundSig`),
sorts by `non_gaussianity` then `U_iso_A2` descending, and emits at most `MAX_SITES = 12` rows,
appending `sites_omitted` with the remainder. On a 52-site run the assistant is shown 12 sites,
not 52. (`symmetryContext` applies the same cap; see Step 11.)

What a user can observe about direction today is therefore (a) the Cartesian components in the
Principal axes table, (b) the Crystal orientation table, and (c) the visual comparison enabled by
the PC ↔ Crystal switch.

### Step 11 — A different displacement measure: `dispA`

The app carries a **second, unrelated** site displacement number, and the two must not be
confused.

**Where:** `structureFromRmc6f` in
[`browserData.js`](../../web_app/frontend/src/browserData.js) — the static / local-file structure
parser, **not** the PCA engine. It runs in the browser only; the Flask structure path returns no
`basis`, so `dispA` simply does not exist in Flask mode.

**It also requires the reference-number column.** The accumulator loop starts with
`if (referenceNumber === null) return;`, so for a **coordinates-only** `.rmc6f` no per-site
accumulators are built at all, `structure.basis` comes back empty, `describeSymmetry` returns
`null` (`if (!structure?.basis?.length || !structure?.latticeVectors) return null;` in
[`symmetryModel.js`](../../web_app/frontend/src/symmetryModel.js)), and no `mean_disp_A` /
`max_disp_A` ever reaches the assistant — even though the PCA page still reconstructs sites for
that same file by folding and clustering.

**Operation.** For each reference site, and separately for each of the three cell axes $i$,
accumulate the circular statistics of the folded within-unit-cell fraction. The fold is a
**positive** modulo, not a plain one, because real `.rmc6f` files carry negative coordinates:

$$w_{i,m}=\bigl((x_{i,m}\,n_i \bmod 1)+1\bigr)\bmod 1\ \in[0,1),$$

$$C_i=\sum_m\cos(2\pi w_{i,m}),\quad S_i=\sum_m\sin(2\pi w_{i,m}),\quad
R_i=\frac{\sqrt{C_i^2+S_i^2}}{N}\in[0,1].$$

$R_i$ is the resultant length; the site's mean fractional position is
$\operatorname{atan2}(S_i,C_i)/2\pi$ folded into $[0,1)$. The circular standard deviation follows
as

$$\sigma^{\text{frac}}_i=\frac{\sqrt{-2\ln \max(R_i,\,10^{-6})}}{2\pi}\quad[\text{unit-cell fractions}],$$

and the axes are combined in quadrature after scaling by the **unit-cell edge lengths**
$\ell_i=\lVert\mathbf L_i\rVert/\max(n_i,1)$ (Å):

$$\texttt{dispA}=\sqrt{\sum_{i=1}^{3}\bigl(\sigma^{\text{frac}}_i\,\ell_i\bigr)^2}\quad[\text{Å}].$$

An axis with $R_i\ge1$ (zero spread, or a single copy) contributes nothing; the $10^{-6}$ floor
caps $\sigma^{\text{frac}}_i$ at $\approx0.837$ cell fractions.

**How it differs from the PCA amplitudes.**

| | `dispA` | PCA `rms` / `eigenvalues` |
| --- | --- | --- |
| Estimator | circular std $\sqrt{-2\ln R}$, population (divide by $N$) | sample variance, unbiased (divide by $N-1$) |
| Frame | the three **fractional cell axes**, then scaled by edge length | Cartesian Å, then rotated to the cloud's own eigenframe |
| Metric | quadrature sum — **ignores the cell metric cross terms** $\mathbf a_i\cdot\mathbf a_j$, so it is not a true Cartesian magnitude for an oblique cell | exact; $C$ is a Cartesian tensor |
| Output | one scalar per site | a $3\times3$ tensor: three amplitudes **and** three directions |
| Runtime | browser structure parser only (no Flask equivalent) | both runtimes |
| Needs reference-number column | **yes** — empty `basis` without it | no (browser worker reconstructs by clustering) |
| Where used | per-Wyckoff-orbit `mean_disp_A` / `max_disp_A` in the AI-assistant context | Principal axes / Summary tables, drawn ellipsoids, assistant `pca` block |

For an **orthogonal** cell and a spread small compared with the cell, the circular std tends to
the linear std and

$$\texttt{dispA}\;\approx\;\sqrt{\lambda_1+\lambda_2+\lambda_3}\;=\;\sqrt{3\,U_\text{iso}},$$

i.e. `dispA` is roughly the total RMS displacement magnitude, the trace information of $C$ with
all directional content discarded. The unit test in
[`browserData.test.js`](../../web_app/frontend/src/__tests__/browserData.test.js) pins the estimator:
two copies at $\pm0.02$ cell fractions on a 10 Å edge give
$\sqrt{-2\ln\cos(2\pi\cdot0.02)}/2\pi\times10 \approx 0.2003$ Å, against the 0.2 Å a linear std
would give.

**Consumers.** `dispA` never appears in a UI panel. It flows
`browserData.structureFromRmc6f` → `structure.basis[]` →
[`symmetryModel.js`](../../web_app/frontend/src/symmetryModel.js) `describeSymmetry` (which attaches
each Wyckoff orbit's member indices into `basis`) → `symmetryContext` in
[`runContext.js`](../../web_app/frontend/src/llm/context/runContext.js), which averages it over each
orbit into `mean_disp_A` and takes the maximum into `max_disp_A`, both rounded to **2 significant
figures**, sorts the orbits by `mean_disp_A` descending, and — like `pcaContext` — emits at most
`MAX_SITES = 12` orbits with a `sites_omitted` count for the rest. That whole chain is gated on
the user having opened the AI Assistant page
([`Dashboard.jsx`](../../web_app/frontend/src/components/Dashboard.jsx), `wantAssistantData`).

---

### Parameters and defaults

Engine-side constants — the covariance denominator, eigenvalue ordering and clamp,
`EIGENVALUE_FLOOR_RATIO`, `DEGENERATE_RATIO`, `MAX_PCA_FIT_POINTS`, `rngSeed`, the grid clamp, the
Jacobi sweep budget and tolerances, `ZERO_SPREAD_VARIANCE` and the kurtosis definedness rule, the anisotropy floors — and every
UI default of the PCA Ellipsoid page are tabulated in the previous section, "Parameters and
defaults". The table below lists only what this section owns.

| Name | Value | Where | Meaning |
| --- | --- | --- | --- |
| eigenvalue ordering | descending $\lambda_1\ge\lambda_2\ge\lambda_3$ | `_eigen_decomposition`, `eigenDecomposition` | defines PC1/PC2/PC3; carries no crystallographic meaning |
| sign rule | largest-\|component\| positive, then $\det P>0$ | `_canonical_axes` (Py), `eigenDecomposition` (JS) | reproducibility only; PC3 may violate step 1 |
| handedness flip | `*= sign(det)` (Py) vs `if (det < 0)` (JS) | same | differs only for the unreachable $\det P = 0$ |
| `invert3` singularity | $\lvert\det\rvert<10^{-12}$ or non-finite → `null` | [`pcaCrystalFrame.js`](../../web_app/frontend/src/pcaCrystalFrame.js) | collapsed-cell guard |
| cosine clamp | $[-1,1]$ before `Math.acos` | `principalAxisOrientation` | round-off guard **only for unit `pcaAxes`** |
| cell-edge normalisation fallback | divisor `cellLengths[i] \|\| 1` | `principalAxisOrientation` | a zero-length edge reports $\cos=0$, $\theta=90°$ |
| angle unit | degrees, `RAD_TO_DEG = 180/π` | `principalAxisOrientation` | — |
| $[u\,v\,w]$ normalisation | $\max_k\lvert\cdot\rvert=1$ (divisor falls back to 1) | `principalAxisOrientation` | no integer reduction |
| supercell division guard | $\lvert n_i\rvert>0$ else 1 | `unitCellVectors` | zero-repeat guard only |
| `dispA` resultant floor | `1e-6` | `structureFromRmc6f` | caps $\sigma^\text{frac}$ at $\approx0.837$ |
| assistant site cap | `MAX_SITES = 12`, plus `sites_omitted` | `pcaContext`, `symmetryContext` in [`runContext.js`](../../web_app/frontend/src/llm/context/runContext.js) | `pcaContext` rounds to 3 s.f., `symmetryContext` to 2 s.f. |
| default reference frame | `pc` | `axisFrame` state in [`PcaKdePage.jsx`](../../web_app/frontend/src/components/PcaKdePage.jsx) | `crystal` disabled when `unitCell` is `null` |
| camera FOV / distance | `45°` / $4.3\times\max(\texttt{halfWidths})$ (`\|\| 1`) | `CAMERA_FOV`, `placeMainCamera` | — |
| `VIEW_UP_AXIS` | `[2, 2, 1]` | `frameAlongAxis` | PC3 up for PC1/PC2, PC2 up for PC3 |
| `CELL_VIEW_UP_AXIS` | `[2, 2, 1]` | `frameAlongCellAxis` | $c$ up for $a$/$b$, $b$ up for $c$ |
| cell-axis up-vector fallback | `up.lengthSq() < 1e-10` (length $<10^{-5}$ Å) → $(0,1,0)$, or $(1,0,0)$ if $\lvert\cos\rvert>0.99$ | `frameAlongCellAxis` | degenerate-framing guard |
| crystal-rod zero skip | $\lVert\mathbf a_j\rVert^2<10^{-12}$ → rod omitted | `buildCrystalAxes` in [`sceneAxes.js`](../../web_app/frontend/src/components/sceneAxes.js) | fewer than three rods, legend unchanged |
| PC triad colours | `#d64545`, `#3fa34d`, `#3f7fd6` | `TRIAD_COLORS` / `PC_CSS_COLORS` in [`sceneAxes.js`](../../web_app/frontend/src/components/sceneAxes.js) | PC1/PC2/PC3, shared by every 3D panel and the table |
| cell axis colours | `#e0a419`, `#18a3a0`, `#b15ad8` | `CELL_AXIS_COLORS` in the same module | $a$/$b$/$c$, deliberately unlike the PC palette |

### Caveats — what this is not

- **A principal axis is a line, not an arrow.** The sign canonicalisation (largest component
  positive, then a handedness flip on PC3) exists for reproducibility across LAPACK builds and
  browsers, not because the direction is physically signed. $\hat{\mathbf e}_i$ and
  $-\hat{\mathbf e}_i$ are the same result. And because the handedness flip runs second, the
  "largest component positive" rule is guaranteed only for PC1 and PC2.
- **Near-degenerate eigenvalues make individual axes meaningless.** For an isotropic or
  near-isotropic site, the axes within the degenerate subspace are arbitrary; only the subspace
  is determined. `axisResolved` (previous section, Step 12b) flags such axes — the panel then prints
  no $\kappa$ and dims their crystal-orientation row. The `degenerate` flag catches only the extreme case
  ($\lambda_3/\max(\lambda_1,10^{-30})<10^{-6}$), not the far more common near-tie between
  $\lambda_1$ and $\lambda_2$. Check the printed $\lambda$ column before reading a direction.
- **The crystal-frame direction the app shows is a folded representative.** The *Crystal
  orientation* table prints angles to $a$/$b$/$c$ and $[u\,v\,w]$ for the sense of each axis that
  makes its closest edge acute (`crystalOrientationRows`, Step 8); the opposite sense is the same
  axis. The fractional↔PCA matrices are computed but not shown (Step 10), and an unresolved axis's
  row is sampling noise (dimmed).
- **$[u\,v\,w]$ is not $(hkl)$.** Even when it is displayed elsewhere or computed by hand from
  `pcaToFrac`: unless the cell metric is isotropic, the lattice direction $[u\,v\,w]$ is not
  perpendicular to the lattice plane $(u\,v\,w)$, and $[u\,v\,w]$ is not a unit vector.
- **The crystal-mode shadow box is not a cell.** It is an orthonormal Gram–Schmidt frame anchored
  to $a$ and $b$ (never $c$), sized from the KDE half-widths. Only its
  $\mathbf e_0$–$\mathbf e_1$ face is a true crystallographic ($a$–$b$) plane for an oblique
  cell, and the rods that carry the real $a$/$b$/$c$ directions are normalised, so relative
  cell-edge lengths are never shown.
- **Covariance is a second moment, not a mechanism.** A large $\lambda_1$ along some direction is
  equally consistent with a soft harmonic mode, a double-well split site, and static strain
  disorder. $\kappa_i$ and the KDE isosurface are what separate those cases; the ellipsoid alone
  cannot.
- **`dispA` is not a PCA amplitude.** It is a circular-statistics scalar on the fractional cell
  axes, with the cell metric's cross terms ignored and all directional content discarded
  (Step 11). It never appears in a UI panel, and it exists only in the browser structure parser.
- **Where the numbers themselves come from is the previous section's problem.** The estimator
  choices that set the *magnitudes* — the subsample cap and the two engines' different draws, the
  browser's Jacobi solver, the KDE bandwidth broadening,
  the display magnification of the Site-ellipsoids markers — are all documented in "PCA Ellipsoid
  page", Caveats. None of them changes a *direction*, except through the eigenframe of a
  subsampled cloud.
