# Displacement Directions — algorithm reference

> Part of the [Algorithms and Math Reference](../ALGORITHMS.md). Every step is anchored to the source; if this document and the code disagree, **the code wins**.

The direction-space counterpart to the ellipsoid: keep which way each atom moved, discard how far, and estimate the probability density over solid angle. The engine section covers the sphere tiling and histogram; the view section covers the rendering.

## Contents

- [Displacement Directions — the orientation-distribution engine](#displacement-directions--the-orientation-distribution-engine)
  - [What this page shows](#what-this-page-shows)
  - [Step 1 — From displacements to directions](#step-1--from-displacements-to-directions)
  - [Step 2 — The frame (`cartesian` | `pca`)](#step-2--the-frame-cartesian--pca)
  - [Step 3 — Choosing the resolution](#step-3--choosing-the-resolution)
  - [Step 4 — The tiling: icosahedron → geodesic subdivision → Goldberg dual](#step-4--the-tiling-icosahedron--geodesic-subdivision--goldberg-dual)
  - [Step 5 — Cell assignment: spherical Voronoi in O(1)](#step-5--cell-assignment-spherical-voronoi-in-o1)
  - [Step 6 — The histogram](#step-6--the-histogram)
  - [Step 7 — Antipodal (inversion) asymmetry](#step-7--antipodal-inversion-asymmetry)
  - [Step 8 — The orientation tensor](#step-8--the-orientation-tensor)
  - [Outputs](#outputs)
  - [Parameters and defaults](#parameters-and-defaults)
  - [Parity: Python engine vs JavaScript port](#parity-python-engine-vs-javascript-port)
  - [Caveats](#caveats)
- [Displacement Directions — the sphere view, axis views, and site picker](#displacement-directions--the-sphere-view-axis-views-and-site-picker)
  - [What this page shows](#what-this-page-shows-1)
  - [Step 1 — Page composition: what the page owns vs. what the view owns](#step-1--page-composition-what-the-page-owns-vs-what-the-view-owns)
  - [Step 2 — Loading a site: the shared `useSiteCloud` hook](#step-2--loading-a-site-the-shared-usesitecloud-hook)
  - [Step 3 — The request: which options reach the engine, and which never do](#step-3--the-request-which-options-reach-the-engine-and-which-never-do)
  - [Step 4 — Cell mesh: the typed-array layout](#step-4--cell-mesh-the-typed-array-layout)
  - [Step 5 — Amplitude relief: the exact radius map](#step-5--amplitude-relief-the-exact-radius-map)
  - [Step 6 — Colour: which scalar, and the contrast transfer](#step-6--colour-which-scalar-and-the-contrast-transfer)
  - [Step 7 — The colorbar and the summary strip](#step-7--the-colorbar-and-the-summary-strip)
  - [Step 8 — Cell outlines](#step-8--cell-outlines)
  - [Step 9 — Hover picking and the tooltip](#step-9--hover-picking-and-the-tooltip)
  - [Step 10 — Frames and axis rods](#step-10--frames-and-axis-rods)
  - [Step 11 — The "Axis views" mini panel](#step-11--the-axis-views-mini-panel)
  - [Step 12 — The site picker: `SiteStructurePanel`](#step-12--the-site-picker-sitestructurepanel)
  - [Step 13 — Reset and Save (PNG export)](#step-13--reset-and-save-png-export)
  - [Step 14 — Computed but not displayed](#step-14--computed-but-not-displayed)
  - [Parameters and defaults](#parameters-and-defaults-1)
  - [Caveats](#caveats-1)

---

## Displacement Directions — the orientation-distribution engine

### What this page shows

The **Displacement Directions** tab (internal page id `orientation` in
[App.jsx](../../web_app/frontend/src/App.jsx), component
[OrientationPage.jsx](../../web_app/frontend/src/components/OrientationPage.jsx) +
[OrientationView.jsx](../../web_app/frontend/src/components/OrientationView.jsx)) answers one
question about a crystallographic site:

> Of all the RMC box copies of this site, **which way** did the atoms move — and is any
> direction preferred beyond what counting noise alone would produce?

It works on exactly the same input as the **PCA Ellipsoid** page: the per-site displacement
cloud $\{\Delta\mathbf{r}_i\}_{i=1}^{N}$, in Cartesian ångström, produced by
`load_site_displacements` in [pca_kde.py](../../rmc_toolkits/pca_kde.py) — an atom's offset from the
*average* position of its own reference site, folded over the supercell boundary, centred on the
site mean and mapped to Cartesian Å through the supercell lattice $\mathsf{L}$. **That
construction is owned by *PCA Ellipsoid page* → Step 1** (and its Step 1b for the browser-only
fold-and-cluster fallback on coordinates-only files); nothing about it is re-derived here. The two
pages then take that cloud apart along orthogonal lines:

| | PCA Ellipsoid page | Displacement Directions page |
|---|---|---|
| Statistic | second moment $\mathbf{U} = \langle \Delta\mathbf{r}\,\Delta\mathbf{r}^{\mathsf T}\rangle$ | direction density $p(\mathbf{u})$ over solid angle |
| Keeps | how far and along which axes, jointly | only which way; $\lVert\Delta\mathbf{r}\rVert$ discarded (or used as a weight) |
| Shape | an ellipsoid — always convex, always centrosymmetric | an arbitrary function on $S^2$ |

The complementarity is the whole point, and it is structural, not stylistic. A covariance
matrix has six free numbers and reproduces exactly one kind of surface. Three physically
common signatures survive in $p(\mathbf{u})$ and are provably invisible in $\mathbf{U}$:

1. **Discrete off-centre spots.** A site hopping between, say, eight $\langle 111\rangle$
   off-centre positions produces eight sharp lobes on the sphere. Its covariance is that of a
   single fat, smooth, convex ellipsoid — the eight lobes and one blob have the same second
   moment.
2. **Odd-order anharmonicity / antipodal asymmetry.** $\mathbf{U}$ is built from
   $\Delta\mathbf{r}\,\Delta\mathbf{r}^{\mathsf T}$, which is invariant under
   $\Delta\mathbf{r} \to -\Delta\mathbf{r}$. Any $+\mathbf{u}$ vs $-\mathbf{u}$ imbalance —
   static off-centring, a one-sided potential — is annihilated by that product. This page
   therefore **never folds the map antipodally**, and reports `antipodalAsymmetry` explicitly.
3. **Any non-convexity at all** — a girdle with a hole in it, a cone, a ring.

Conversely, this page throws away the radial information the ellipsoid is made of. It is a
complement, not a replacement; the `amplitude`/`amplitude2` weights and the *relief* display
bring a controlled amount of that radial information back.

The engine of record is [rmc_toolkits/orientation.py](../../rmc_toolkits/orientation.py) (NumPy,
served through the Flask route `/api/pca/orientation` in
[app.py](../../web_app/backend/app.py)). The JavaScript port
[workers/orientation.js](../../web_app/frontend/src/workers/orientation.js) runs inside the shared
PCA worker ([pcaKdeWorker.js](../../web_app/frontend/src/workers/pcaKdeWorker.js), request
`kind: 'orientation'`). The routing switch is in `requestPca`
([useSiteCloud.js](../../web_app/frontend/src/useSiteCloud.js)) and it is **not** a
Flask-vs-static-build switch: whenever a `.rmc6f` has been loaded as a browser file (the Demo,
or a picked folder), the worker answers — *in both runtimes*, Flask backend running or not.
Only a typed backend directory goes through the HTTP route. The two engines produce the same
JSON shape up to the small divergences listed under "Parity" below.

---

### Step 1 — From displacements to directions

**Input.** $\Delta\mathbf{r}_i \in \mathbb{R}^3$, Cartesian ångström, $i = 1 \dots N_{\text{tot}}$
(`totalPoints`). For a single site, $N_{\text{tot}}$ is the number of supercell copies of that
site; for an element-pooled cloud it is the sum over that element's sites.

**Which cloud.** `site_orientation_histogram` / `siteOrientationHistogram` select it through
`displacement_cloud` ([pca_kde.py](../../rmc_toolkits/pca_kde.py)) and the matching `else`-chain in
`siteOrientationHistogram`:

| Selector | Cloud |
|---|---|
| `reference_number` given | that one site's copies |
| `element` given (not `""`/`"all"`) | every site of that element, pooled |
| **neither** (or `element` in `{"", "all"}`) | **every site of every element, pooled** |

`displacement_cloud` is the *same* selector the ellipsoid engine uses, and the reason pooling is
meaningful at all — each site is separately mean-centred before the union is taken — is set out in
*PCA Ellipsoid page* → Step 2. Two differences matter here: that page's KDE caps its fit at
`MAX_PCA_FIT_POINTS = 20000` points while this engine subsamples nothing (Step 2, detail 4 below,
and caveat 2 of the page section), and this engine's *default* selection is the pooled-everything
row above.

The UI always sends a `referenceNumber`, so the last row is only reachable from a script, a
hand-made HTTP request, or a direct library call — but it is the *default* of both engines, and
it mixes species with different masses and different $\mathbf{U}$ into one cloud. That is almost
never the intended read (caveat 11).

**Amplitude.** $a_i = \lVert \Delta\mathbf{r}_i \rVert$ (Å).

**Direction.**

$$\mathbf{u}_i \;=\; \frac{\Delta\mathbf{r}_i}{\lVert \Delta\mathbf{r}_i \rVert} \;\in\; S^2 ,
\qquad \lVert\mathbf{u}_i\rVert = 1 \ \text{(dimensionless)} .$$

**The amplitude floor.** $\mathbf{u}_i$ is undefined at $\Delta\mathbf{r}_i = \mathbf{0}$ and
numerically meaningless just above it: an atom sitting 10⁻¹² Å from its site mean has a
"direction" that is pure round-off from the mean subtraction, not physics. The engine
therefore always drops points below

```python
NEGLIGIBLE_AMPLITUDE = 1e-9   # Angstrom   orientation.py:91
```

**User cutoffs.** Two optional, cumulative cuts sharpen a weak pattern:

- `min_amplitude` / `minAmplitude` — an absolute floor in Å. **Not validated:** a negative value
  is accepted, has no effect (the $10^{-9}$ floor dominates), and is then hidden by the
  `max(threshold, 0)` in the reported `amplitudeCutoff`. In JS a non-numeric string becomes
  `NaN` via `Number(minAmplitude)`, which makes every `>` comparison false and surfaces as the
  survival error below, not as a type error.
- `min_amplitude_quantile` / `minAmplitudeQuantile` $\in [0, 1)$ — the $q$-quantile of the
  $\{a_i\}$ themselves (NumPy's default linear-interpolation quantile; the JS port
  reimplements that exact rule in `quantile`). The quantile is taken over **all $N_{\text{tot}}$
  amplitudes, before any cut** — including the ones below the $10^{-9}$ Å floor — not over the
  survivors.

The effective threshold and keep-mask are

$$t \;=\;\begin{cases}
\texttt{min\_amplitude} & q = 0 \quad\text{(the default)}\\[2pt]
\max\big(\texttt{min\_amplitude},\; Q_q(\{a_i\}_{i=1}^{N_{\text{tot}}})\big) & q > 0
\end{cases}
\qquad
\text{keep } i \iff a_i > \max(t,\; 10^{-9}\,\text{Å}) .$$

The $q = 0$ branch is a genuine short circuit, not a simplification:
[orientation.py:588](../../rmc_toolkits/orientation.py) and
[orientation.js:509](../../web_app/frontend/src/workers/orientation.js) both guard with
`if quantile > 0`, so $Q_0 = \min_i a_i$ is never evaluated — which matters, because with
`min_amplitude = 0` a literal $\max(0, Q_0)$ would drop the single smallest point.

Note the **strict** inequality: an atom exactly on its site mean is always dropped. The
surviving count is `usedPoints` ($N$ from here on), and `rejectedPoints` $= N_{\text{tot}} - N$.

*Rationale for the quantile knob.* Near-zero displacements are isotropically distributed by
construction (their direction is noise), so including them dilutes any real anisotropy toward
uniform. A modest quantile cut therefore raises contrast without inventing structure — but it
also discards data, so it moves the Poisson noise floor up. Both are reported.

**Errors.** `orientation_histogram` / `orientationHistogram` raise six errors in all
(`ValueError` in Python, `Error` in JS), four of them input validation before any work:

| Message | Raised when |
|---|---|
| `vectors must be a numeric array with shape (N, 3)` — JS: `vectors must be an array of [x, y, z] rows` | wrong input shape |
| `weight must be one of ('count', 'amplitude', 'amplitude2')` | unknown weight |
| `frame must be one of ('cartesian', 'pca')` | unknown frame |
| `min_amplitude_quantile must lie in [0, 1)` | quantile outside $[0,1)$ |
| `no displacement vectors survive the amplitude cutoff` | $N = 0$ after the cut |
| `displacement weights sum to zero` | $\sum_m M_m \le 0$ |

The last is unreachable through the ordinary path — every surviving $a_i > 10^{-9}$ Å, so every
weight is strictly positive — but it guards a caller who reaches the function directly.
`goldberg_tiling` / `goldbergTiling` adds `frequency must lie in [1, 64]`. `site_orientation_histogram`
adds `Unknown reference number …` / `Unknown element …` from the site selection. The Flask route
maps all of them: `ValueError` → **400**, `PermissionError` → **403**, `FileNotFoundError` →
**404**, anything else → **500** ([app.py](../../web_app/backend/app.py),
`pca_orientation_endpoint`).

Code: `orientation_histogram` in [orientation.py](../../rmc_toolkits/orientation.py);
`orientationHistogram` in [orientation.js](../../web_app/frontend/src/workers/orientation.js).

---

### Step 2 — The frame (`cartesian` | `pca`)

`FRAMES = ("cartesian", "pca")`.

- **`cartesian`** (default) — the directions stay in the Cartesian ångström frame in which the
  `.rmc6f` supercell lattice vectors are expressed (`displacements = centered @ lattice_vectors`,
  [pca_kde.py](../../rmc_toolkits/pca_kde.py) `load_site_displacements`). It is the *same* basis the
  unit-cell vectors $\mathbf{a},\mathbf{b},\mathbf{c}$ are derived in —
  $\mathbf{a}_i = \mathbf{L}_i/n_i$, see *Principal axes in the crystallographic frame* → Step 5 —
  which is why no re-referencing is ever needed between a direction on this sphere and a cell edge.
  This is the crystal frame,
  so the sphere's axes are directly comparable to $\mathbf{a}$, $\mathbf{b}$, $\mathbf{c}$ — but
  see the obliquity note below for what "comparable" means in a non-orthogonal cell. The view
  draws the a/b/c rods in this mode when the file carries cell metadata, and falls back to the
  site's PC axes when it does not (Step "Outputs", `pcaAxes`).
- **`pca`** — every direction is rotated into the cloud's own principal axes before binning:
  $\mathbf{u}' = \mathbf{A}\,\mathbf{u}$, where the rows of $\mathbf{A}$ are the three PCA axes
  in descending eigenvalue order. PC1 lands on local $+x$, PC2 on $+y$, PC3 on $+z$, which
  makes different sites (and different elements) directly superimposable. The rods drawn in this
  mode are the literal identity triad, not `pcaAxes` — in the rotated frame they *are* PC1/2/3.

Five details matter for interpretation:

1. **The PCA rotation is fitted on *every* vector, before any amplitude cut**, so the frame does
   not move when a quantile cut is active. It is computed only if $N_{\text{tot}} \ge 4$; below
   that $\mathbf{A} = \mathbf{I}$, which means **`frame="pca"` silently degrades to exactly the
   `cartesian` map for $N_{\text{tot}} < 4$, with nothing in the payload flagging it** (`frame`
   still echoes `"pca"`).
2. **$\mathbf{A}$ comes from the shared `_eigen_decomposition`** (Step 8), including the sign and
   handedness canonicalisation, so PC1 here *is* PC1 there. What a signed axis does and does not
   mean is owned by *Principal axes in the crystallographic frame* → Step 3.
3. **The frame agrees with the ellipsoid page to round-off, not bit-for-bit.** The covariance
   here is `np.cov(centered, rowvar=False, bias=False)`
   ([orientation.py:584-585](../../rmc_toolkits/orientation.py); `covariance3(vectors)` at
   [orientation.js:505](../../web_app/frontend/src/workers/orientation.js), which re-derives the mean
   itself) — the unbiased ($n-1$) estimator of the *re-centred* vectors — while `site_ellipsoids`
   builds the same tensor from a raw `bincount` outer-product sum divided by $\max(n-1, 1)$ with
   **no** mean subtraction. Both estimators, and the fact that the Python cloud is already centred
   once in `load_site_displacements` while the JS helpers subtract a numerically-zero mean again,
   are *PCA Ellipsoid page* → Step 3 ("Where the mean is removed"). The consequence here is only
   that the summation orders differ in the last bits: measured
   $\max|\mathbf{A}_{\text{orient}} - \mathbf{A}_{\text{ellipsoid}}| = 5\times10^{-16}$ on a
   synthetic $8^3$ site (and `np.array_equal` is `False`).
4. **For large element-pooled clouds the two frames genuinely differ.** The orientation engine
   never subsamples; the KDE path caps its fit at `MAX_PCA_FIT_POINTS = 20000` points with
   `rng_seed = 0`, and the two engines draw *different* subsets above that (*PCA Ellipsoid page*
   → Step 2, "Subsampling"). Above 20 000 pooled points the PCA-KDE page's frame is a subsample
   fit and this page's is the full fit — a real difference, not round-off. (Per-site
   `site_ellipsoids`, which is what the site picker and the panel header show, is not
   subsampled either.)
5. **Switching the frame re-bins the data; it does not relabel it.** The tiling is a pure
   function of $\nu$ (`@lru_cache` keyed on `frequency` alone, [orientation.py:378](../../rmc_toolkits/orientation.py))
   and the rotation is applied to the *directions* before `assign_cells`
   ([orientation.py:597-611](../../rmc_toolkits/orientation.py)). So the same data binned in the two
   frames lands in different cells, and `counts`, `expected`, `zScore`, `peakCell`,
   `emptyFraction`, `significance` and `antipodalAsymmetry` all change. Measured on a 3000-point
   cloud at $\nu = 6$: `counts` differ, `peakCell` 123 → 328, `emptyFraction` 0.0608 → 0.0525,
   `significance` 3.343 → 3.354, $\mathcal{A}$ 0.1753 → 0.1747. Only the point-level quantities
   (`totalPoints`, `usedPoints`, `meanAmplitude`, `rmsAmplitude`) are frame-invariant.

**Which frame each output is in.** `centers`, `polygons`, `peakDirection` and `orientationAxes`
are in the **request's** frame — PCA-rotated when `frame="pca"`, because `tensor` at
[orientation.py:661](../../rmc_toolkits/orientation.py) is built from the already-rotated
`directions`. `pcaAxes` is always in the **original Cartesian** frame (it comes from line 585,
pre-rotation). A consumer that wants a PCA-frame vector back in crystal coordinates must
rotate by $\mathbf{A}^{\mathsf T}$:
$\mathbf{v}_{\text{cart}} = \mathbf{A}^{\mathsf T}\mathbf{v}_{\text{pca}}$, with
$\mathbf{A} = $ `pcaAxes` (rows).

**Obliquity.** The Cartesian frame is whatever the `.rmc6f` lattice vectors express, so for a
non-orthogonal cell $\mathbf{a}$, $\mathbf{b}$, $\mathbf{c}$ are oblique to each other and to the
sphere's own $x/y/z$. The a/b/c rods the view draws are then not an orthonormal triad, and
"the lobe points along $\mathbf{a}$" is a statement about a direction, not about a coordinate —
the general form of that trap (a Cartesian component is not a fractional index, and
$[u\,v\,w]$ is the contraction that converts one to the other) is *Principal axes in the
crystallographic frame* → "The oblique-cell subtlety, stated plainly".
Separately, the tiling's own orientation is fixed by the icosahedron construction (§4.1 puts
2-fold axes on Cartesian $x$, $y$, $z$) and bears no relation to the crystal's point group, so
symmetry-equivalent lobes land in cells of different area and shape and get different
$\Omega_m$, `expected` and $z$. Both points are in the Caveats.

---

### Step 3 — Choosing the resolution

The tiling resolution is the geodesic **frequency** $\nu$, and the cell count is
$C = 10\nu^2 + 2$ (derived in Step 4). Passing `frequency=None` (the UI's *Auto*) invokes the
over-binning guard:

```python
def recommended_frequency(n_points, *, target_per_cell=12, max_frequency=24):
    cells = max(12.0, n_points / target_per_cell)
    frequency = round(sqrt(max(cells - 2.0, 10.0) / 10.0))
    return clip(frequency, MIN_FREQUENCY, min(max_frequency, MAX_FREQUENCY))
```

i.e. invert $C = 10\nu^2 + 2$ at the cell count that gives `DEFAULT_TARGET_PER_CELL = 12`
points per cell:

$$\nu_{\text{rec}} \;=\; \operatorname{round}\!\left(\sqrt{\frac{\max\!\big(N/12,\,12\big) - 2}{10}}\,\right),
\qquad \text{clipped to } [1, 24] .$$

**Where the floors bite.** `cells` is pinned at 12 for every $N \le 144$, and the inner
$\max(\cdot, 10)$ pins the radicand at 1, so $\nu_{\text{rec}} = 1$ there. Rounding then holds
$\nu = 1$ well past that: measured transition points are

| $N$ range | Auto $\nu$ | cells |
|---|---|---|
| $\le 293$ | **1** | 12 (the dodecahedron) |
| 294 – 774 | 2 | 42 |
| 775 – 1493 | 3 | 92 |
| 1494 – 2454 | 4 | 162 |
| 2455 – 3653 | 5 | 252 |
| $\ge 10\,854$ | $\ge 10$ | $\ge 1002$ |

Two consequences. **Auto can pick a resolution the manual control cannot**: the dropdown starts
at $\nu = 2$, so $\nu = 1$ is reachable only through Auto. And a typical site with a few hundred
copies gets $\nu = 1$ or 2 from Auto — twelve or forty-two cells — against the UI's fixed default
of $\nu = 10$ (1002 cells). That gap is the whole content of caveat 7.

**Why 12.** At $n$ points per cell the fractional Poisson noise on an isotropic map is
$1/\sqrt{n}$; $n = 12$ gives $\approx 29\%$ cell-to-cell scatter. That is small enough that a
genuine $1.5\times$ lobe stands out and large enough that noise does not manufacture one.

**The failure mode it guards.** Push $\nu$ past the data and every cell holds 0 or 1 points.
The map becomes Poisson confetti — and confetti on a sphere looks exactly like structure to the
eye. This is the single easiest way to over-read this kind of plot, which is also why `zScore`
and `significance` are computed and displayed.

**Rounding is load-bearing for cross-runtime agreement.** Python's `round()` is round-half-to-
even. At exactly 774 surviving points the square root is exactly 2.5, so Python gives $\nu = 2$
while a naive `Math.round` would give 3 — two different tilings for the same data in the two
runtimes. The JS port therefore implements `roundHalfEven` explicitly, and *both* test suites
pin the same three values (`recommendedFrequency(774) == 2`, `(300) == 2`, `(12000) == 10`).

Bounds: `MIN_FREQUENCY = 1` ([orientation.py:85](../../rmc_toolkits/orientation.py)),
`MAX_FREQUENCY = 64` (line 86, $C = 40\,962$ cells). `recommended_frequency` never exceeds 24 on
its own; the UI dropdown offers `auto, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24` and **defaults to
$\nu = 10$** (1002 cells), not Auto.

---

### Step 4 — The tiling: icosahedron → geodesic subdivision → Goldberg dual

`goldberg_tiling(frequency)` / `goldbergTiling(frequency)` builds the bin geometry. It depends
on nothing but $\nu$, so it is cached — `@lru_cache(maxsize=8)` in Python (per process),
an unbounded `Map` in JS (per worker). Measured build cost on the dev machine: 0.01 s at
$\nu = 10$, 0.04 s at $\nu = 24$, 0.25 s at $\nu = 64$ (Python); 0.01 / 0.03 / 0.11 s (JS).
*(The module docstring's "~1 s to tessellate" at $\nu = 64$ is a conservative overestimate.)*

#### 4.1 The icosahedron

The 12 vertices are the cyclic permutations of $(0, \pm 1, \pm \varphi)$ with the golden ratio

$$\varphi = \frac{1 + \sqrt 5}{2} \approx 1.618\,034 ,$$

generated in `_icosahedron()` as the three families $(0, s_1, s_2\varphi)$,
$(s_1, s_2\varphi, 0)$, $(s_2\varphi, 0, s_1)$ over $s_1, s_2 \in \{+1,-1\}$, then normalised to
the unit sphere.

The 20 faces are **derived, not tabulated**: the edge length is the minimum pairwise vertex
separation *excluding self-pairs* (`distances[distances > 1e-9].min()`,
[orientation.py:165](../../rmc_toolkits/orientation.py)); two vertices are adjacent when

$$\big|\,d_{ab} - \ell_{\text{edge}}\,\big| \;<\; 10^{-9}\,\max(\ell_{\text{edge}},\,1)$$

— a *relative* tolerance, not an absolute one ([orientation.py:166](../../rmc_toolkits/orientation.py);
[orientation.js:80,85](../../web_app/frontend/src/workers/orientation.js) is identical) — and a face
is any mutually-adjacent triple. The construction asserts it found exactly 20. Each triple is
then wound counter-clockwise as seen from outside by **swapping $b \leftrightarrow c$ when**

$$\big((\mathbf{v}_b-\mathbf{v}_a)\times(\mathbf{v}_c-\mathbf{v}_a)\big)\cdot(\mathbf{v}_a+\mathbf{v}_b+\mathbf{v}_c) \;<\; 0 ,$$

keeping $(a,b,c)$ otherwise. The strict `< 0` matters: an exact 0 keeps the original order in
both engines. (It cannot occur for an icosahedron, but the *parity* of the two ports depends on
the comparison being written identically.) Deriving rather than hard-coding makes winding and
completeness properties of the code instead of a table that could be transcribed wrong — and it
is why the Python and JS engines produce *identical cell indices* (they enumerate in the same
order).

#### 4.2 Frequency-$\nu$ geodesic subdivision

For each face $f$ with corners $(\mathbf{A}, \mathbf{B}, \mathbf{C})$ and integers
$i, j \ge 0$, $k = \nu - i - j \ge 0$:

$$\mathbf{P}(f,i,j) \;=\; \frac{k\,\mathbf{A} + i\,\mathbf{B} + j\,\mathbf{C}}{\nu}\;
\Big/\; \left\lVert \frac{k\,\mathbf{A} + i\,\mathbf{B} + j\,\mathbf{C}}{\nu} \right\rVert .$$

This is the standard **Class I** geodesic sphere: subdivide the flat triangle uniformly, then
project radially. Points shared between faces are merged by an exact **combinatorial** key —
`("v", corner_id)` for a corner, `("e", min, max, position)` for an edge point (with the
position flipped when the edge is traversed the other way, in `_edge_key`), `("f", face, i, j)`
for a face interior point — never by rounding coordinates, so the merge cannot fail on a
floating-point tie.

**Vertex count.** Each face carries $(\nu+1)(\nu+2)/2$ lattice points, 20 faces give
$10(\nu+1)(\nu+2)$; the 12 corners are each counted 5 times (subtract $12 \times 4 = 48$) and
the 30 edges each carry $\nu-1$ interior points counted twice (subtract $30(\nu-1)$):

$$V \;=\; 10(\nu+1)(\nu+2) - 48 - 30(\nu-1) \;=\; \boxed{10\nu^2 + 2} .$$

The construction asserts this count and raises if it is not met.

**Triangles.** Per face, $\nu(\nu+1)/2$ "upward" triangles $\{(i,j), (i{+}1,j), (i,j{+}1)\}$ and
$\nu(\nu-1)/2$ "downward" triangles $\{(i{+}1,j),(i{+}1,j{+}1),(i,j{+}1)\}$, i.e. $\nu^2$ per
face and $F_\triangle = 20\nu^2$ overall, all inheriting the parent face's counter-clockwise
winding (the $(i,j)$ parameter plane maps $\mathbf{A}\to(0,0)$, $\mathbf{B}\to(\nu,0)$,
$\mathbf{C}\to(0,\nu)$, which preserves orientation). With $E_\triangle = 30\nu^2$ this
triangulation has $\chi = (10\nu^2+2) - 30\nu^2 + 20\nu^2 = 2$, as a sphere must.

#### 4.3 The Goldberg dual

**Cells are the geodesic *vertices*, not its triangles.** Each of the $C = 10\nu^2+2$ geodesic
vertices becomes one cell centre $\mathbf{c}_m$ (unit vector), and the polygon drawn around it
has one vertex per incident triangle, placed at that triangle's **spherical circumcentre**.

For a counter-clockwise triangle of unit vectors $(\mathbf{p}_0,\mathbf{p}_1,\mathbf{p}_2)$,

$$\mathbf{q} \;=\; \frac{(\mathbf{p}_1-\mathbf{p}_0)\times(\mathbf{p}_2-\mathbf{p}_0)}
{\lVert (\mathbf{p}_1-\mathbf{p}_0)\times(\mathbf{p}_2-\mathbf{p}_0)\rVert}$$

is the unit normal of the plane through the three points; since they all lie in that plane,
$\mathbf{q}\cdot\mathbf{p}_0 = \mathbf{q}\cdot\mathbf{p}_1 = \mathbf{q}\cdot\mathbf{p}_2$, so
$\mathbf{q}$ is **equiangular** from all three — the spherical circumcentre — and the CCW
winding makes it point outward. That is exactly the definition of a spherical Voronoi vertex
for those three centres, which is what ties the drawn polygon to the assignment rule (Step 5).

The result is a **Goldberg polyhedron** $\mathrm{GP}(\nu, 0)$: hexagons plus exactly 12
pentagons — the family of geodesic domes, fullerene cages and H3-style geospatial grids.

#### 4.4 Why exactly 12 pentagons, always

Not an implementation shortcut, and not removable. Let the tiling have $p$ pentagons and $h$
hexagons, with every vertex trivalent (three cells meet at each dual vertex — true here because
each geodesic triangle contributes one polygon vertex to each of its three corner cells):

$$F = p + h, \qquad E = \tfrac{5p + 6h}{2}, \qquad V = \tfrac{2E}{3} = \tfrac{5p+6h}{3} .$$

Then

$$\chi \;=\; V - E + F \;=\; \frac{5p+6h}{3} - \frac{5p+6h}{2} + (p+h)
\;=\; -\frac{5p+6h}{6} + p + h \;=\; \frac{p}{6} .$$

For a sphere $\chi = 2$, hence $p = 12$ — **independent of $h$, and therefore of $\nu$**. Set
$p = 0$ (all hexagons) and you get $\chi = 0$: a torus, not a sphere. The 12 pentagons sit at
the vertices of the parent icosahedron.

Both engines assert this at construction time and raise
`dual tiling is not a Goldberg polyhedron (bad cell degrees)` otherwise. Both test suites check
$C = 10\nu^2+2$ and exactly 12 degree-5 cells at $\nu \in \{1,2,3,5,8\}$. At $\nu = 1$ the dual
is the regular **dodecahedron**: 12 cells, all pentagons, no hexagons (verified).

#### 4.5 Tangent basis and polygon vertex ordering

Polygon vertices arrive in triangle order, which is not the order a renderer or an area formula
needs. `_angular_order` / `angularOrder` sorts each cell's incidences counter-clockwise about
its own centre.

At a unit normal $\mathbf{n} = \mathbf{c}_m$, `_tangent_basis` builds

$$\mathbf{e}_1 = \frac{\mathbf{r} - \mathbf{n}(\mathbf{r}\cdot\mathbf{n})}
{\lVert \mathbf{r} - \mathbf{n}(\mathbf{r}\cdot\mathbf{n})\rVert}, \qquad
\mathbf{e}_2 = \mathbf{n} \times \mathbf{e}_1 ,$$

with $\mathbf{r}$ the Cartesian axis **least** aligned with $\mathbf{n}$ (`argmin |n|`), which
keeps the projection well conditioned. The frame is right-handed:
$\mathbf{e}_1 \times \mathbf{e}_2 = \mathbf{e}_1\times(\mathbf{n}\times\mathbf{e}_1) = \mathbf{n}$.
Each incident point $\mathbf{p}$ is projected into the tangent plane,
$\mathbf{p}_\perp = \mathbf{p} - \mathbf{n}(\mathbf{p}\cdot\mathbf{n})$, and ordered by

$$\theta = \operatorname{atan2}\big(\mathbf{p}_\perp\cdot\mathbf{e}_2,\; \mathbf{p}_\perp\cdot\mathbf{e}_1\big) \in (-\pi, \pi] ,$$

ascending — counter-clockwise seen from outside. Python does this with one stable
`np.lexsort((angle, owners))`; JS sorts each owner's list with a stable comparator. Ties (which
do not occur for distinct circumcentres) fall back to incidence order in both.

`polygons` is stored padded to six vertices; a pentagon repeats its last vertex, and `sizes`
carries the true count (5 or 6). The JSON returned to the browser is trimmed back to `sizes`.

#### 4.6 Exact solid angle of each cell

Cell areas are **not** approximated. `_spherical_polygon_areas` / `sphericalPolygonArea` sums a
signed triangle fan from vertex $\mathbf{v}_0$, using the **van Oosterom–Strackee** form of
**Girard's theorem** (spherical excess). For unit vectors $\mathbf{a},\mathbf{b},\mathbf{c}$:

$$\Omega(\mathbf{a},\mathbf{b},\mathbf{c}) \;=\; 2\,\operatorname{atan2}\!\Big(
\mathbf{a}\cdot(\mathbf{b}\times\mathbf{c}),\;\;
1 + \mathbf{a}\cdot\mathbf{b} + \mathbf{b}\cdot\mathbf{c} + \mathbf{c}\cdot\mathbf{a}\Big)
\quad\text{[steradian]},$$

$$\Omega_m \;=\; \sum_{i=1}^{n_m - 2} \Omega(\mathbf{v}_0, \mathbf{v}_i, \mathbf{v}_{i+1}) .$$

This is algebraically Girard's $\Omega = A + B + C - \pi$, rewritten as a tangent-half-angle so
that it stays accurate for the small, near-equilateral cells here, where differencing three
interior angles against $\pi$ loses most of its significant digits. **L'Huilier's formula is not
used** (it needs the three arc lengths and is aimed at a different conditioning problem).

The padded pentagon vertex is harmless by construction: with $\mathbf{b} = \mathbf{c}$ the
numerator $\mathbf{a}\cdot(\mathbf{b}\times\mathbf{b})$ is *exactly* 0 and the denominator
$2(1+\mathbf{a}\cdot\mathbf{b}) > 0$, so `atan2(0, +) = 0` — the degenerate fan step contributes
exactly nothing.

Both engines assert $\sum_m \Omega_m = 4\pi$ and raise otherwise — but **with different
tolerances**, one of the few places the two ports already diverge:

| Engine | Test | Effective bound on $\lvert\delta\rvert = \lvert\sum\Omega_m - 4\pi\rvert$ |
|---|---|---|
| Python, [orientation.py:426](../../rmc_toolkits/orientation.py) | `np.isclose(sum, 4*pi, rtol=1e-9)` | $10^{-8} + 10^{-9}\cdot 4\pi \approx 2.26\times10^{-8}$ (NumPy's default `atol = 1e-8` is included, and contributes about 44 % of the bound) |
| JS, [orientation.js:351](../../web_app/frontend/src/workers/orientation.js) | `abs(total - 4*pi) > 1e-9 * 4*pi` | $10^{-9}\cdot 4\pi \approx 1.26\times10^{-8}$ (purely relative) |

The Python check is therefore about $1.8\times$ looser. Nothing has ever tripped either, but a
tiling that tripped the JS assertion and passed the Python one would be a silent cross-runtime
divergence.

#### 4.7 Adjacency and the antipode map

Every ordered pair $(v,w)$ occurs exactly once across all counter-clockwise triangles, so the
directed edge list *is* the adjacency list with no de-duplication needed. It is sorted with the
same `_angular_order` routine, and the engine asserts each cell's neighbour count equals its
polygon degree. `neighbors` is a $(C, 6)$ table padded with $-1$ for the 12 pentagons.

The tiling is centrosymmetric (the icosahedron is, and radial subdivision preserves it), so
$-\mathbf{c}_m$ is exactly another cell centre. `antipode[m]` is resolved by running $-\mathbf{c}_m$
through the ordinary assignment path (one code path, and it doubles as an assertion): the
construction raises if $\max_m \lVert \mathbf{c}_{\text{antipode}(m)} + \mathbf{c}_m \rVert_\infty > 10^{-9}$.
Both test suites additionally check that the map is an involution.

#### 4.8 How equal-area is it, really

Not perfectly, and this is the honest number. Cells are largest near the icosahedral face
centres and smallest at the 12 pentagons, because the Class I construction spaces points evenly
on the *flat* triangle and the radial (gnomonic) projection compresses angular spacing away from
the face centre. Measured $\Omega_m$ relative to the mean $4\pi/C$:

| $\nu$ | cells | min / mean | max / mean | max / min | pentagon / mean |
|---|---|---|---|---|---|
| 4 | 162 | 0.752 | 1.083 | 1.44 | 0.752 |
| 8 | 642 | 0.662 | 1.172 | 1.77 | 0.662 |
| **10** (UI default) | **1002** | **0.645** | **1.184** | **1.84** | **0.645** |
| 16 | 2562 | 0.619 | 1.198 | 1.93 | 0.619 |
| 24 | 5762 | 0.605 | 1.204 | 1.99 | 0.605 |
| 64 | 40962 | 0.588 | 1.206 | 2.05 | 0.588 |

The 12 pentagons are always the smallest cells, and the spread widens monotonically with $\nu$.
Area standard deviation is $\approx 12\%$ of the mean at both $\nu = 6$ and $\nu = 10$.

> **Test-coverage note.** `test_area_spread_is_modest` (Python) and its JS twin assert
> $\max\Omega/\min\Omega < 2$ — but only at $\nu = 8$ (measured 1.77). The invariant is **false**
> at $\nu \ge 40$ (2.03 at 40, 2.05 at 64), which `MAX_FREQUENCY = 64` permits. The assertion is
> a sanity check at a working resolution, not a property of the construction.

Two consequences:

- **This is still far more uniform than an equal-angle $(\theta,\varphi)$ grid**, whose cells
  shrink as $\sin\theta$ and collapse to zero-area slivers at the poles. That is the comparison
  the design is making.
- **Raw counts are emphatically not a density.** Plot $n_m$ directly and the 12 pentagons show
  up as $\approx 35\%$ dark spots at $\nu = 10$ — i.e. the map prints the parent icosahedron
  onto your data. Dividing by each cell's own $\Omega_m$ (Step 6) removes this exactly, which is
  why the exact solid-angle computation of §4.6 is not a nicety.

> **Documentation discrepancy.** The `orientation.py` module docstring says cells "range over
> roughly +/-10% at usual frequencies". The measured range at the UI's default $\nu = 10$ is
> $-36\%$ to $+18\%$ (ratio 1.84). The $\pm 10\%$ figure is about right for the *hexagons'*
> standard deviation, but understates the pentagons and the full spread. Nothing downstream
> depends on the claim — every count is divided by its own exact $\Omega_m$ — but the docstring
> number should be corrected.

---

### Step 5 — Cell assignment: spherical Voronoi in O(1)

**The rule.** A direction $\mathbf{u}$ belongs to the cell whose centre it is nearest to in
angle:

$$m(\mathbf{u}) \;=\; \arg\max_m\; \mathbf{u}\cdot\mathbf{c}_m
\qquad(\text{maximising } \cos\theta \equiv \text{minimising } \theta).$$

Cell boundaries are therefore the **spherical Voronoi diagram** of the centres — and, by §4.3,
the polygon drawn for cell $m$ is *exactly* that Voronoi region, because its vertices are the
circumcentres (Voronoi vertices) of the incident triangles. The drawn shape, the binned region
and the divided-by area are the same object. Both suites verify this non-circularly by Monte
Carlo (see "Tests" below).

A literal $\arg\max$ costs $O(C)$ per point — 1002 dot products per atom at $\nu = 10$. `_assign`
/ `assignOne` does it in $O(1)$ in three stages.

**(a) Locate the icosahedral face — exact ray/cone test.** For face $f$ with corner matrix
$\mathbf{M}_f = [\mathbf{A}\,|\,\mathbf{B}\,|\,\mathbf{C}]$ (corners as *columns*), solve
$\mathbf{u} = \lambda_0\mathbf{A} + \lambda_1\mathbf{B} + \lambda_2\mathbf{C}$ by the
precomputed inverse `face_inverse[f]` $= \mathbf{M}_f^{-1}$:

$$\boldsymbol{\lambda}^{(f)} = \mathbf{M}_f^{-1}\mathbf{u} .$$

$\mathbf{u}$ lies in face $f$'s cone iff all $\lambda^{(f)}_a \ge 0$; for every other face at
least one is negative. So

$$f^\star = \arg\max_f\ \min_a \lambda^{(f)}_a$$

picks the containing face with **no epsilon**, and breaks exact edge ties deterministically
(first maximum wins, in both engines).

**(b) Invert the gnomonic map to a lattice index.** Normalising
$\boldsymbol{\beta} = \boldsymbol{\lambda}/\sum_a\lambda_a$ gives the barycentric coordinates of
the point where the ray $\mathbf{u}$ pierces the *plane* of the face triangle — a **gnomonic
(central) projection**, which is linear in barycentric coordinates and therefore free. The
lattice points of §4.2 sit at barycentric $\big(k/\nu, i/\nu, j/\nu\big)$, so

$$(k, i, j) \;\approx\; \nu\,\boldsymbol{\beta}$$

rounded to integers that still sum to $\nu$, by the **largest-remainder rule**: floor all three,
then hand the leftover $\nu - \sum \lfloor\cdot\rfloor$ units to the components with the largest
fractional parts (Python: stable `argsort` of the negated remainders; JS: stable `Array.sort` —
same tie behaviour). `lattice[f*, i, j]` is the seed cell.

For any valid unit direction $\sum_a \beta_a = 1$, so the leftover is in $\{0,1,2\}$ and the sum
is restored exactly. Both implementations, however, only ever hand out **at most 3** units
(Python: `slot < deficit` with `slot ∈ {0,1,2}`, [orientation.py:353-354](../../rmc_toolkits/orientation.py);
JS: `for (slot = 0; slot < 3 && deficit > 0; ...)`,
[orientation.js:259-262](../../web_app/frontend/src/workers/orientation.js)), so a degenerate input
whose barycentrics do not sum to 1 is silently under-filled rather than rejected.

Four clamps and guards keep the path total, all mirrored in both engines:

| Guard | Code | Effect |
|---|---|---|
| barycentric denominator floored at $10^{-300}$ | `chosen.sum(...)` → `np.maximum(..., 1e-300)`, [orientation.py:343](../../rmc_toolkits/orientation.py) / [orientation.js:255](../../web_app/frontend/src/workers/orientation.js) | a zero $\boldsymbol{\lambda}$ cannot divide by zero |
| $\boldsymbol{\beta}$ clipped to $[0,1]$ before scaling by $\nu$ | `nu * np.clip(barycentric, 0.0, 1.0)`, :347 / :256 | keeps `floor` inside the lattice |
| `_normalize` divides by $\max(\lVert\mathbf{v}\rVert, 10^{-300})$ | :145 / :45 | a zero vector normalises to zero instead of raising |
| $-1$ lattice slot degrades to cell 0 | :361 / :265 | pathological input gives a valid cell, not an index error |

Together these mean a **zero vector is assigned, not rejected**: all $\lambda = 0$, the first
face wins the `argmax`, $\boldsymbol{\beta} = \mathbf{0}$, and the leftover-3 rule lands it
deterministically on `lattice[0][1][1]` (for $\nu \ge 3$). It can only get there through a direct
`assign_cells` call — the histogram path drops zero vectors at the $10^{-9}$ Å floor.

The seed is only approximate because the gnomonic map is not an isometry: equal planar spacing
is *not* equal angular spacing, and the distortion grows toward the triangle corners.

**(c) Greedy walk on the cell-adjacency graph.** From the seed, repeatedly move to whichever
neighbour has a larger $\mathbf{u}\cdot\mathbf{c}$, stopping when none does.

*Why it provably reaches the true nearest cell.* The cell-adjacency graph here is the geodesic
triangulation, which is the **Delaunay** graph dual to the spherical Voronoi diagram of the
centres. Delaunay graphs have the "no false local optimum" property: for any query point
$\mathbf{u}$ and any site $s$ that is *not* the nearest site to $\mathbf{u}$, $s$ has a Delaunay
neighbour strictly closer to $\mathbf{u}$. Hill-climbing therefore cannot stall anywhere except
at the global nearest site.

*Why it terminates immediately.* The gnomonic seed and the truth differ only near cell
boundaries, never by more than a cell or two, independent of $\nu$. Measured on 20 000 uniform
random directions:

| $\nu$ | seed already correct | improving rounds used | mismatches vs brute force |
|---|---|---|---|
| 6 | 97.5 % | 1 | 0 / 20000 |
| 10 | 98.1 % | 1 | 0 / 20000 |
| 24 | 98.8 % | 1 | 0 / 20000 |
| 64 | 98.9 % | 1 | 0 / 20000 |

Both engines cap the walk at **8 rounds**. In every case measured, one improving round plus one
confirming round suffices, at every resolution.

**Two entry points, one of which normalises.** `assign_cells` / `assignCells` (the public
wrapper) calls `_normalize` first, so callers may pass unnormalised vectors
([orientation.py:452-462](../../rmc_toolkits/orientation.py)). The internal `_assign` / `assignOne`
does **not** — it assumes unit input. Two consequences worth knowing: the antipode map of §4.7 is
built by calling `_assign` directly on `-centers` (already unit, line 434), and the histogram
path calls `assign_cells` on rows that `orientation_histogram` has *already* divided by their
amplitude (line 611), so the re-normalisation there is redundant work, not a correctness
requirement.

---

### Step 6 — The histogram

#### 6.1 Weights: what each cell accumulates

`WEIGHTS = ("count", "amplitude", "amplitude2")`. With $w_i$ the weight of surviving point $i$
and $m$ its cell:

$$M_m \;=\; \sum_{i \,:\, m(\mathbf{u}_i) = m} w_i ,
\qquad
w_i = \begin{cases}
1 & \texttt{"count"} \quad\text{(default)}\\[2pt]
a_i \ \text{[Å]} & \texttt{"amplitude"}\\[2pt]
a_i^2 \ \text{[Å}^2\text{]} & \texttt{"amplitude2"}
\end{cases}$$

- **`count`** is the orientation distribution proper — every atom votes once, however far it
  moved. This is the only weight for which "density" means "probability of a direction".
- **`amplitude`** biases toward the directions the big excursions took.
- **`amplitude2`** is the **angular decomposition of the mean-square displacement**. This is
  exact and worth stating precisely: since $\Delta\mathbf{r} = a\,\mathbf{u}$,
  $\sum_i a_i^2 \mathbf{u}_i\mathbf{u}_i^{\mathsf T} = \sum_i \Delta\mathbf{r}_i\Delta\mathbf{r}_i^{\mathsf T}$,
  so integrating this map against $\mathbf{u}\mathbf{u}^{\mathsf T}$ reproduces
  $\mathbf{U}/\operatorname{tr}\mathbf{U}$ — the very tensor the ellipsoid page draws. Both sides
  are *uncentred* second moments, so the identity is exact for the already-centred site clouds
  this page consumes (measured $3\times10^{-16}$; see Step 8 for the one thing that breaks it).

Alongside the weighted mass the engine always keeps the **raw counts**
$n_m = \#\{i : m(\mathbf{u}_i) = m\}$ and the per-cell amplitude sum $S_m = \sum_i a_i$.

#### 6.2 Density, and why the division by $\Omega_m$

$$\rho_m \;=\; \frac{M_m}{\big(\sum_{m'} M_{m'}\big)\,\Omega_m}
\qquad [\text{sr}^{-1}] ,
\qquad\text{so}\qquad \sum_m \rho_m \Omega_m = 1 .$$

$\rho$ is a genuine density over solid angle: it integrates to 1 over the sphere, and it is
resolution-independent in the limit (doubling $\nu$ does not halve it). Dividing by each cell's
*own* $\Omega_m$ — not by the mean $4\pi/C$ — is what removes the $-36\%$ pentagon artefact of
§4.8 exactly.

The display quantity is the dimensionless **enhancement**:

$$\boxed{\;\mathcal{E}_m \;=\; 4\pi\,\rho_m\;}
\qquad
\begin{cases}
\mathcal{E} = 1 & \text{everywhere, for an isotropic site}\\
\mathcal{E}_m = 1.8 & \text{"this direction is 1.8}\times\text{ more likely than chance"}\\
\mathcal{E}_m < 1 & \text{depletion}
\end{cases}$$

`vmin`/`vmax` are $\min_m \mathcal{E}_m$ and $\max_m \mathcal{E}_m$; the colourbar runs
$0 \to \texttt{vmax}$ with a tick drawn at the isotropic $1\times$.

#### 6.3 Poisson significance against the isotropic null

**Null model.** The $N$ surviving directions are i.i.d. uniform on the sphere. Then
$n_m \sim \mathrm{Binomial}(N, \Omega_m/4\pi) \approx \mathrm{Poisson}(e_m)$ with

$$e_m \;=\; \frac{N\,\Omega_m}{4\pi} \qquad \text{(reported as } \texttt{expected} \text{),}$$

and, using the Poisson identity $\operatorname{Var} = \operatorname{E}$,

$$z_m \;=\; \frac{n_m - e_m}{\sqrt{\max(e_m,\,10^{-12})}} .$$

Three things to note:

1. **$z$ always comes from the raw counts $n_m$, never from $M_m$ and never after smoothing.**
   Smoothing correlates neighbouring cells, so a post-smoothing $z$ would overstate the evidence
   by design. This means that under `amplitude`/`amplitude2` weighting, **the colour and the
   significance measure different things** — the colour shows weighted mass per steradian, the
   hover $z$ shows how surprising the *head count* is.
2. The exact multinomial standard deviation is $\sqrt{e_m(1 - \Omega_m/4\pi)}$, smaller than
   $\sqrt{e_m}$ by a factor $\sqrt{1-1/C}$ — 0.05 % at $\nu = 10$. Using the Poisson form is
   therefore very slightly conservative.
3. A whole-map summary is reported as

$$\texttt{significance} \;=\; \sqrt{\tfrac{1}{C}\textstyle\sum_m z_m^2}\ \ [\sigma],$$

   the RMS $z$. $\approx 1$ means the entire pattern is consistent with pure counting noise;
   well above 1 means real directional structure. This is the number that stops an over-binned
   map from being read as physics.

#### 6.4 Neighbour smoothing

`smoothing` is a count of diffusion passes on the cell-adjacency graph
(`_smooth` / `smooth`). Each pass, with $\alpha = \texttt{SMOOTHING\_ALPHA} = 0.5$
([orientation.py:95](../../rmc_toolkits/orientation.py)) and $d_m \in \{5,6\}$ the cell's degree:

$$M_m' \;=\; (1-\alpha)\,M_m \;+\; \alpha \sum_{m'\,\in\,\mathcal{N}(m)} \frac{M_{m'}}{\max(d_{m'},\,1)} .$$

The $\max(d, 1)$ is in both engines ([orientation.py:500](../../rmc_toolkits/orientation.py),
[orientation.js:431](../../web_app/frontend/src/workers/orientation.js)); it never binds for a real
tiling, where every degree is 5 or 6.

The pass count is **truncated toward zero**, not rounded: `int(passes)` / `Math.trunc(smoothing)`.
A fractional `smoothing = 2.9` runs 2 passes. A negative value runs **none** (the whole block is
guarded by `if smoothing > 0`) but is still echoed verbatim into the result as
`"smoothing": int(smoothing)` — so a payload can report `smoothing: -3` on a completely unsmoothed
map.

Each cell keeps $1-\alpha$ of its mass and splits the rest evenly over its own neighbours, so
**the total is invariant by construction** — the smoothed map still integrates to the same
number of points (or the same $\sum a_i$, or $\sum a_i^2$) as the raw one. Both suites assert
$\sum_m M_m$ is unchanged to 9 decimal places, and that smoothing lowers the peak without moving
the lobe.

Applied to: `mass`, the per-cell amplitude sum $S_m$, and a float copy of the counts (the
"count field"). **Not** applied to `counts` itself, hence not to `zScore`.

> **Caveat — smoothing is mass-conserving but not area-weighted.** The diffusion moves mass
> *per cell*, taking no account of the $\pm$36 %/+18 % area spread of §4.8. A perfectly
> isotropic input (mass $\propto \Omega_m$) does not stay perfectly isotropic. Measured at the UI
> default $\nu = 10$, an exactly-isotropic map reads:
>
> | passes | resulting $\mathcal{E}$ range |
> |---|---|
> | 1 | 0.978 – 1.035 |
> | **2 (UI default)** | **0.976 – 1.060** |
> | 4 | 0.955 – 1.098 |
> | 12 (slider max) | 0.903 – 1.202 |
>
> The excess sits in the ring of hexagons *around* each of the 12 pentagons (mean 1.060 there
> versus 0.999 elsewhere at 2 passes), so heavy smoothing paints a faint icosahedral pattern
> onto an isotropic site. At the default 2 passes it is a $\pm 6\%$ effect — well below the
> Poisson noise of any real map — but at 8–12 passes it approaches the size of a weak real lobe.
> Read `zScore`, which is untouched by smoothing, when in doubt.

#### 6.5 Per-cell mean amplitude (the relief)

$$\langle a \rangle_m \;=\; \begin{cases} S_m / n^{\text{f}}_m & n^{\text{f}}_m > 10^{-12}\\ 0 & \text{otherwise}\end{cases}
\qquad [\text{Å}]$$

where $S_m$ and the count field $n^{\text{f}}_m$ are both taken *after* smoothing. Smoothing the
numerator and denominator separately — rather than smoothing the ratio — keeps a smoothed relief
the honest mean of the same smoothed population the colour shows, and keeps every value inside
the raw data's amplitude range (asserted by `test_cell_mean_amplitude_smoothing_stays_bounded`).
Empty cells report exactly 0, never NaN.

This is the **radial-relief** quantity. `reliefFactors` in
[orientationSphere.js](../../web_app/frontend/src/orientationSphere.js) turns it into a per-cell
scale factor

$$\phi_m \;=\; \begin{cases}
1 & \langle a\rangle_m \le 0 \ \text{ or } \ \langle a\rangle \le 0 \quad\text{(empty cell: exactly neutral)}\\[2pt]
\mathrm{clamp}_{[0.3,\,2.2]}\big(1 + \texttt{relief}\cdot(\langle a\rangle_m/\langle a\rangle - 1)\big) & \text{otherwise,}
\end{cases}$$

with $\langle a\rangle$ = `meanAmplitude`. Note the empty-cell rule: a cell with
`cellMeanAmplitude == 0` is pinned to **exactly 1**, not to the lower clamp 0.3 — an empty
direction sits on the unit sphere rather than denting in.

The surface is **not** then bulged cell-by-cell, which would tear it open at every cell border.
`vertexRadii` keys each polygon vertex by its rounded coordinates (circumcentres are shared by up
to three cells, an identity lost across the JSON boundary) and displaces it to the **mean of the
relief factors of the cells sharing it**, so the relief surface stays crack-free while each cell
keeps one flat colour.

Because the relief is computed independently of the colour weighting, **shape and colour carry
separate information**: colour = how often this direction, shape = how far when it happens.

---

### Step 7 — Antipodal (inversion) asymmetry

This is the readout that exists because the ellipsoid cannot have one. Using the exact
`antipode` map of §4.7:

$$\mathcal{A} \;=\; \frac{1}{2N}\sum_{m=1}^{C} \big| n_m - n_{\text{antipode}(m)} \big|
\;=\; \frac{1}{N}\sum_{\text{unordered pairs}} \big|n(\mathbf{u}) - n(-\mathbf{u})\big| .$$

The factor $\tfrac12$ is because each unordered antipodal pair appears twice in the sum over
cells. Range and meaning:

- $\mathcal{A} = 0$ — perfectly inversion-symmetric cloud (what a harmonic site gives).
- $\mathcal{A} = 1$ — fully one-sided: every occupied cell's antipode is empty.

**$\mathcal{A}$ is bounded above by 1 by construction.** Since
$|n_m - n_{\text{ant}(m)}| \le n_m + n_{\text{ant}(m)}$ and $\sum_m n_m = N$,

$$\mathcal{A} \;\le\; \frac{1}{2N}\sum_m \big(n_m + n_{\text{ant}(m)}\big) \;=\; \frac{2N}{2N} \;=\; 1 .$$

This matters for the flag below, because the null floor has no such bound.

**The null floor.** Pure counting noise produces a non-zero $\mathcal{A}$, so a bare value is
unreadable. For two i.i.d. $\mathrm{Poisson}(\mu)$ cells with $\mu = N/C$,
$X - Y \approx \mathcal{N}(0, 2\mu)$ and $\operatorname{E}|X-Y| = \sqrt{4\mu/\pi} = 2\sqrt{\mu/\pi}$.
Summing over $C$ cells and applying the $1/(2N)$ prefactor:

$$\mathcal{A}_{\text{null}} \;=\; \frac{1}{2N}\,C\cdot 2\sqrt{\frac{N/C}{\pi}}
\;=\; \boxed{\sqrt{\frac{C}{\pi N}}} \qquad (\texttt{antipodalAsymmetryNull}).$$

The Gaussian step assumes $\mu = N/C \gg 1$. **At the UI's default operating point it is not:**
$\nu = 10$, $C = 1002$, and a 216-copy site gives $\mu \approx 0.22$. The formula is still what
both engines report; read it as an order-of-magnitude floor there, not as a calibrated
expectation.

**The $3\times$ flag is structurally unreachable at high resolution.** The UI flags the pair red
when $\mathcal{A} > 3\,\mathcal{A}_{\text{null}}$
([OrientationView.jsx:447](../../web_app/frontend/src/components/OrientationView.jsx)). Because
$\mathcal{A} \le 1$ while $\mathcal{A}_{\text{null}} = \sqrt{C/\pi N}$ can exceed $1/3$ — or
exceed 1 outright — the condition can only ever fire when

$$3\sqrt{\frac{C}{\pi N}} < 1 \iff N > \frac{9C}{\pi} .$$

| $\nu$ | $C$ | $N$ needed for the flag to be *possible* |
|---|---|---|
| 1 (Auto floor) | 12 | 35 |
| 2 | 42 | 121 |
| 6 | 362 | 1 038 |
| 8 | 642 | 1 840 |
| **10 (UI default)** | **1002** | **2 871** |
| 24 | 5762 | 16 507 |

Measured, to make the failure concrete: a **perfectly one-sided** 216-point cloud at
$\nu = 10$ reports `usedPoints = 216`, $\mathcal{A} = 1.000$ (its maximum),
$\mathcal{A}_{\text{null}} = 1.215$, flag **False**, and `significance` $= 1.35\sigma$ — the page
reads a maximally asymmetric cloud as consistent with noise. The same cloud on *Auto*
($\nu = 1$) gives $\mathcal{A} = 0.991$, $\mathcal{A}_{\text{null}} = 0.133$, flag **True**.
**So the antipodal readout, exactly like the map itself, must be read on Auto resolution for
typical site copy counts** (caveat 7).

Both suites assert that a deliberately one-sided cloud clears $3\times$ the floor while a
centrosymmetric Gaussian lands between $0.2\times$ and $2\times$ it — at test resolutions where
$N > 9C/\pi$ holds.

Note $\mathcal{A}$ is built from raw counts, so it is unaffected by smoothing and by the weight
choice.

---

### Step 8 — The orientation tensor

$$\mathbf{T} \;=\; \frac{\sum_i w_i\,\mathbf{u}_i\mathbf{u}_i^{\mathsf T}}{\sum_i w_i}
\qquad\text{(dimensionless, } \operatorname{tr}\mathbf{T} = 1 \text{ since } \lVert\mathbf{u}_i\rVert = 1).$$

Computed **from the points, not from the binned map**, so it is resolution-independent — it does
not change when you move the $\nu$ slider. It uses the same weights $w_i$ and the same frame as
the histogram, and only the surviving points.

Eigen-decomposition gives $\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge 0$ with
$\lambda_1+\lambda_2+\lambda_3 = 1$ (`orientationEigenvalues`) and the corresponding axes as
rows (`orientationAxes`). Reading them (the standard Woodcock classification):

| Eigenvalue pattern | Shape on the sphere | Physical reading |
|---|---|---|
| $\lambda_1 \approx \lambda_2 \approx \lambda_3 \approx \tfrac13$ | uniform | isotropic site |
| $\lambda_1 \gg \lambda_2 \approx \lambda_3$ | two lobes on one axis (bipolar) | one soft direction — a rod |
| $\lambda_1 \approx \lambda_2 \gg \lambda_3$ | a great-circle **girdle** | one *stiff* direction, free in the plane — a disc |

The scalar summary reported is

$$\texttt{orientationAnisotropy} = 3\lambda_1 - 1 \in [0, 2],$$

0 for a uniform sphere, 2 for a perfect single axis.

**Exact identity worth knowing.** With `weight="amplitude2"` and no amplitude cut,

$$\mathbf{T} \;=\; \frac{\sum_i a_i^2\,\mathbf{u}_i\mathbf{u}_i^{\mathsf T}}{\sum_i a_i^2}
\;=\; \frac{\sum_i \Delta\mathbf{r}_i\Delta\mathbf{r}_i^{\mathsf T}}{\operatorname{tr}\sum_i \Delta\mathbf{r}_i\Delta\mathbf{r}_i^{\mathsf T}}
\;=\; \frac{\mathbf{U}}{\operatorname{tr}\mathbf{U}} ,$$

where $\mathbf{U}$ is the ADP tensor `site_ellipsoids` reports. Two things about that identity
are easy to get wrong:

- **Any $n$ vs $n-1$ convention cancels identically.** It is a scalar factor on both numerator
  and denominator of $\mathbf{U}/\operatorname{tr}\mathbf{U}$, so it cannot contribute a residual
  of any size. (Measured: for a synthetic cloud, $|\mathbf{T} - \mathbf{M}/\operatorname{tr}\mathbf{M}|$
  is $1.1\times10^{-15}$ using $\mathbf{M} = \mathbf{c}^{\mathsf T}\mathbf{c}$ **and**
  $1.1\times10^{-15}$ using $\mathbf{M}/(n-1)$ — bit-for-bit the same.)
- **$\mathbf{T}$ is built from the *uncentred* second moment** $\sum_i \Delta\mathbf{r}\Delta\mathbf{r}^{\mathsf T}$
  ([orientation.py:661](../../rmc_toolkits/orientation.py)), and so is `site_ellipsoids`'
  $\mathbf{U}$ (a raw `bincount` outer-product sum). The only thing that can break the identity
  is therefore a **non-zero cloud mean**, which `np.cov` would subtract and neither of these
  does. Against `np.cov` the same synthetic cloud gives $4.2\times10^{-5}$ — the entire
  discrepancy is the mean term.

For a real site the cloud is already centred by `load_site_displacements`, so the identity holds
exactly: measured $3\times10^{-16}$ (and $9.7\times10^{-17}$ against `site_ellipsoids`'
`covariance` on a synthetic $8^3$ site). In the `pca` frame $\mathbf{T}$'s eigenvalues are then
just the ellipsoid's normalised eigenvalues. So `amplitude2` weighting is the setting in which
this page and the ellipsoid page are looking at *literally the same tensor*, just decomposed
differently.

#### The shared `_eigen_decomposition`, and why sharing it is mandatory

Both the PCA frame (Step 2) and $\mathbf{T}$ go through the *same* routine the ellipsoid engine
uses — `_eigen_decomposition` / `_canonical_axes` in
[pca_kde.py](../../rmc_toolkits/pca_kde.py), imported by `orientation.py`
([orientation.py:80](../../rmc_toolkits/orientation.py)), and `eigenDecomposition` in
[pcaKde.js](../../web_app/frontend/src/workers/pcaKde.js), imported by `orientation.js`. It sorts
eigenvalues **descending** (LAPACK `eigh` returns ascending), clamps them at $\ge 0$, returns the
axes as **rows** (`eigh` gives columns), and canonicalises the axes in two steps: flip each axis so
its largest-magnitude component is positive, then flip the third axis if the frame came out
left-handed.

**Those two steps, their order, and their limits are owned by *Principal axes in the
crystallographic frame* → Step 3** and are not restated here. Three of its findings change how the
outputs of *this* page must be read, so they are recorded rather than assumed:

- **Handedness wins over the sign rule.** Negating $\hat{\mathbf e}_3$ for handedness makes its
  largest-magnitude component negative again, so the "largest component positive" property holds
  for PC1 and PC2 only; a PC3 row that violates it is expected, not a bug. `pcaAxes` and
  `orientationAxes` inherit this exactly.
- **The two engines write the handedness flip differently** — Python
  `vectors[..., 2] *= sign(det)` ([pca_kde.py](../../rmc_toolkits/pca_kde.py) `_canonical_axes`) vs
  JS `if (det < 0) axes[2] = -axes[2]` ([pcaKde.js:121](../../web_app/frontend/src/workers/pcaKde.js)).
  For orthonormal eigenvectors $|\det| = 1$, so the $\det = 0$ case that would *zero* Python's
  third axis is unreachable; the operations are still not literally the same.
- **A canonicalised axis is still an undirected line.** $\mathbf{T}$ is built from
  $\mathbf{u}\mathbf{u}^{\mathsf T}$ and is blind to $\mathbf{u} \to -\mathbf{u}$ (caveat 8), so a
  sign printed in `orientationAxes` or `pcaAxes` is a reproducibility convention, never a claim
  that the atoms move *toward* that end. Only the map itself and `antipodalAsymmetry` (§7) carry
  that information.

The solvers behind the routine — LAPACK `eigh` on one side, the 50-sweep cyclic Jacobi with its
$10^{-18}$ absolute off-diagonal tolerance and $10^{-300}$ rotation skip on the other — are
described in *PCA Ellipsoid page* → Step 3; their consequence for this page is in "Parity" below.

**Why sharing is required, not merely tidy.** The entire purpose of `frame="pca"` and of the
`orientationAxes` output is to be *compared against the ellipsoid's axes* — "the direction
distribution peaks along PC1", "the girdle normal is PC3". If the two engines canonicalised
signs differently, PC1 here could be $-$PC1 there and every such statement would silently invert.
There is no way to detect that from the numbers alone, so the two views must share one
convention by sharing one function.

Honest limitation, inherited: for a **near-isotropic** cloud the axes within a degenerate
subspace are arbitrary, and no convention can fix that. In the `pca` frame such a site's sphere
is arbitrarily rotated about the degenerate axis — and, because the two runtimes use different
eigensolvers (below), they can pick *different* arbitrary rotations for the same data.

---

### Outputs

Every field of the returned object, with units. "Shown" = surfaced anywhere in the UI (verified
by grepping [OrientationView.jsx](../../web_app/frontend/src/components/OrientationView.jsx) and
[OrientationPage.jsx](../../web_app/frontend/src/components/OrientationPage.jsx), excluding tests and
the workers themselves). "Frame" is `req` for the request's frame (PCA-rotated when
`frame="pca"`), `cart` for the original Cartesian frame, `—` for scalars and per-cell scalars.

| Field | Type / units | Frame | Meaning | Shown |
|---|---|---|---|---|
| `frequency`, `cellCount`, `pentagonCount` | int | — | $\nu$, $C = 10\nu^2+2$, always 12 | **no** — the Resolution dropdown prints `ν=n (10n²+2)` from its own constant, client-side; **in Auto mode the resolved $\nu$ and $C$ are never surfaced**. `cellCount`'s only frontend use is a hover bounds guard (OrientationView.jsx:452) |
| `totalPoints`, `usedPoints`, `rejectedPoints` | int | — | before / after the amplitude cut | no |
| `amplitudeCutoff` | Å | — | $\max(t, 0)$ — **excludes** the $10^{-9}$ floor | no |
| `weight`, `smoothing`, `frame` | str / int / str | — | echo of the request | `weight`, `smoothing` (footer line) |
| `pcaAxes` | 3×3, rows | **cart** | PCA axes from the full cloud, *pre*-rotation | only as the fallback reference triad when `frame="cartesian"` **and** the file carries no unit-cell metadata. The PCA frame draws the identity basis; the crystal frame draws a/b/c (OrientationView.jsx:389-399) |
| `centers` | $(C,3)$ unit vectors | req | cell centre directions | yes (hover) |
| `polygons` | ragged $(C, n_m, 3)$ | req | Voronoi cell outlines, CCW from outside (omitted if `geometry=false`) | yes (mesh) |
| `neighbors` | $(C,6)$, $-1$ padded | — | adjacency (omitted if `geometry=false`) | no |
| `areas` | $(C,)$ sr | — | exact cell solid angles | no |
| `sizes` | $(C,)$ 5 or 6 | — | polygon degree | no |
| `antipode` | $(C,)$ | — | index of the cell at $-\mathbf{c}_m$ | no |
| `counts` | $(C,)$ int | — | **raw** head count per cell (never smoothed) | yes (hover) |
| `mass` | $(C,)$ | — | weighted, **post**-smoothing | no |
| `density` | $(C,)$ sr⁻¹ | — | from smoothed `mass`; integrates to 1 | no |
| `enhancement` | $(C,)$ dimensionless | — | $4\pi\rho$, **smoothed**; 1 = isotropic | yes (colour) |
| `expected`, `zScore` | $(C,)$ | — | isotropic null and its $z$, both from **raw** counts | `zScore` (hover) |
| `vmin`, `vmax` | dimensionless | — | enhancement range | `vmax` (colourbar) |
| `meanCount`, `emptyFraction` | float | — | $N/C$; fraction of cells with **raw** count 0 | no |
| `meanAmplitude`, `rmsAmplitude` | Å | — | $\langle a\rangle$, $\sqrt{\langle a^2\rangle}$ over survivors | `meanAmplitude` (relief scale) |
| `cellMeanAmplitude` | $(C,)$ Å | — | per-cell $\langle a\rangle_m$, smoothed numerator **and** denominator | yes (relief + hover) |
| `antipodalAsymmetry`, `…Null` | dimensionless | — | §7, from raw counts | yes |
| `orientationTensor`, `orientationEigenvalues` | 3×3 / 3 | req / — | §8 | no |
| `orientationAxes` | 3×3, rows | **req** | §8 — PCA-frame axes when `frame="pca"` | no |
| `orientationAnisotropy` | dimensionless | — | $3\lambda_1 - 1$ | yes |
| `peakCell`, `peakEnhancement` | int / float | — | argmax over the **smoothed** enhancement | `peakEnhancement` only |
| `peakDirection` | unit vector | **req** | `centers[peakCell]` | yes |
| `peakZScore` | float | — | the **raw** $z$ at that cell | yes |
| `significance` | $\sigma$ | — | RMS $z$ over cells, raw counts | yes |
| `recommendedFrequency` | int | — | what Auto would have chosen | no |
| `referenceNumber`, `element`, `siteFractional` | int / str / 3 fractional | — | site tag added by `site_orientation_histogram`; `siteFractional` is $\mathrm{mod}(\text{site mean}\times\text{supercell},\,1)$, i.e. the site's mean position reduced into one unit cell ([pca_kde.py:172](../../rmc_toolkits/pca_kde.py)) | **no** — the panel header's element and site number come from `selectedEllipsoid` (the `/api/pca/sites` payload), not from this result; `siteFractional` here is read by nothing (the site panel reads `site.siteFractional` from the sites payload instead) |
| `source` | str | — | Flask only: resolved `.rmc6f` path | no |
| `browserOrientation` | `true` | — | JS only: marks the browser engine | yes ("· browser") |

**Smoothed vs raw, in one line.** Smoothed: `mass`, `density`, `enhancement`, `vmin`, `vmax`,
`cellMeanAmplitude`, `peakCell`, `peakEnhancement`. Raw: `counts`, `expected`, `zScore`,
`peakZScore`, `significance`, `emptyFraction`, `antipodalAsymmetry`. Two consequences worth
naming: with smoothing on, `peakCell` (smoothed argmax) and `peakZScore` (raw $z$ *at* that cell)
describe different populations; and `emptyFraction` does **not** shrink when smoothing is applied,
even though after one pass essentially every cell holds mass.

**Computed but never displayed:** `frequency`, `cellCount`, `pentagonCount`, `mass`, `density`,
`expected`, `areas`, `sizes`, `antipode`, `neighbors`, the whole orientation tensor and its axes,
`peakCell`, `recommendedFrequency`, `meanCount`, `emptyFraction`, `rmsAmplitude`,
`amplitudeCutoff`, `referenceNumber`, `element`, `siteFractional`, and the point-budget triple
(`totalPoints`, `usedPoints`, `rejectedPoints`). They are all in the API/worker payload, so a
script or notebook can use them.

---

### Parameters and defaults

**Engine parameters.** The "Transport" column records what each parameter is actually reachable
through: **HTTP** = read by `pca_orientation_endpoint` in [app.py](../../web_app/backend/app.py);
**worker** = forwarded by the `kind: 'orientation'` branch of
[pcaKdeWorker.js](../../web_app/frontend/src/workers/pcaKdeWorker.js); **UI** = sent by
`OrientationView`'s fetch. Two parameters are library-only, and one is browser-only.

| Parameter | Key (JS) | Transport | Range | Engine default | UI default | Effect |
|---|---|---|---|---|---|---|
| `frequency` | `frequency` | HTTP + worker + UI | 1 – 64, or `None`/absent = auto | `None` (auto) | **10** (1002 cells) | geodesic frequency; $C = 10\nu^2+2$ |
| `weight` | `weight` | HTTP + worker + UI | `count` \| `amplitude` \| `amplitude2` | `count` | `count` | what each cell accumulates |
| `min_amplitude` | `minAmplitude` | HTTP + worker, **never sent by the UI** | **not validated** — negatives are accepted and reported as `amplitudeCutoff = 0` | `0.0` | — | absolute short-displacement cut |
| `min_amplitude_quantile` | `minAmplitudeQuantile` | HTTP + worker + UI | $[0, 1)$, validated | `0.0` | `0` (slider caps at 0.5) | quantile short-displacement cut; skipped entirely at $q = 0$ |
| `smoothing` | `smoothing` | HTTP + worker + UI | int, truncated toward 0; $\le 0$ = no passes | `0` | **2** (slider max 12) | neighbour-diffusion passes |
| `frame` | `frame` | HTTP + worker + UI | `cartesian` \| `pca` | `cartesian` | `cartesian` | crystal frame vs PCA frame |
| `geometry` | `geometry` | HTTP + worker + UI | bool | `True` | `true` | include `polygons` + `neighbors` |
| `target_per_cell` | `targetPerCell` | **library call only** — neither the Flask route nor the worker forwards it | int $\ge 1$ | `12` | not exposed | auto-resolution target |
| `reference_number` / `element` | `referenceNumber` / `element` | HTTP + worker; UI sends `referenceNumber` only | — | both `None` → **all sites of all elements pooled** | from the site picker | site selection (`element` pools all sites of one element; `""`/`"all"` normalise to the pooled default in both transports) |
| — | `clusterThreshold` | **browser worker only** — no Python equivalent | 0.4 – 2.5 Å, slider step 0.1 | `DEFAULT_CLUSTER_THRESHOLD = 1.5` Å ([pcaKde.js:472](../../web_app/frontend/src/workers/pcaKde.js)) | 1.5 Å, control shown only when `sites.reconstructed` | fold-and-cluster distance used to *reconstruct* sites when the file has no reference-site/cell columns; part of the worker's parse cache key, so changing it re-derives every displacement and therefore every number on the page |

**The `clusterThreshold` path has no server counterpart.** The fold-and-cluster site
reconstruction (`sitesByClustering` in [pcaKde.js](../../web_app/frontend/src/workers/pcaKde.js)) is
browser-only, and the Python parser cannot read a coordinates-only `.rmc6f` at all, so such a file
simply cannot be served by `/api/pca/orientation`. The algorithm, the clustering guards and the
clean/merged/fragmented labelling are *PCA Ellipsoid page* → Step 1b; caveat 12 below records what
those reconstructed sites do to *this* page's numbers.

**Module constants** ([orientation.py](../../rmc_toolkits/orientation.py) lines 85–103, mirrored at
[orientation.js](../../web_app/frontend/src/workers/orientation.js) lines 24–31):

| Constant | Value | Unit | Role |
|---|---|---|---|
| `MIN_FREQUENCY` | 1 | — | $\nu=1$ is the dodecahedron (12 cells, all pentagons) |
| `MAX_FREQUENCY` | 64 | — | 40 962 cells; 0.25 s (Py) / 0.11 s (JS) to tessellate |
| `NEGLIGIBLE_AMPLITUDE` | $10^{-9}$ | Å | always-on floor; below it a direction is round-off |
| `SMOOTHING_ALPHA` | 0.5 | — | fraction of mass a cell exports per pass |
| `DEFAULT_TARGET_PER_CELL` | 12 | points/cell | auto-resolution target ($\approx 29\%$ Poisson scatter) |
| `recommended_frequency(max_frequency=)` | 24 | — | auto-resolution never exceeds this |
| greedy-walk cap | 8 | rounds | measured need: 1 |

**Display-only options** (owned by `OrientationPage.jsx`, never sent to the engine):
`colormap` = `viridis`, `contrast` = 1.0 (0.5 – 3.0), `relief` = 0.5, `showOutline` = true,
`showAxes` = true. The contrast transfer (`colorCoordinate` in
[orientationSphere.js](../../web_app/frontend/src/orientationSphere.js)) is

$$t \;=\; \mathrm{clamp}_{[0,1]}\!\Big(\pi_0 + c\,\big(\mathrm{clamp}_{[0,1]}(\mathcal{E}/\mathcal{E}_{\max}) - \pi_0\big)\Big),
\qquad \pi_0 = \mathrm{clamp}_{[0,1]}\!\big(1/\mathcal{E}_{\max}\big),
\qquad t = 0 \ \text{ if } \mathcal{E}_{\max} \le 0 .$$

All **three** `clamp01` calls are in the code, and the pivot clamp is not decorative: on a wholly
depleted map ($\mathcal{E}_{\max} < 1$) the isotropic level sits off the top of the bar and
$1/\mathcal{E}_{\max}$ would exceed 1. The pivot is the isotropic level, so this is a symmetric
gain about $1\times$; $c = 1$ is exactly the plain linear mapping. The colourbar gradient is
painted through the identical transfer, so bar and sphere always agree. It changes colours only,
never the numbers.

---

### Parity: Python engine vs JavaScript port

**Which runtime runs what.** The split is *not* Flask-build vs static-build; it is whether a
`.rmc6f` **text** has been loaded into the browser.

| Data source | Path | Engine |
|---|---|---|
| A **typed backend directory** (no browser file loaded) | `GET /api/pca/orientation` → `site_orientation_histogram` | [orientation.py](../../rmc_toolkits/orientation.py) — **source of truth** |
| A **locally-loaded `.rmc6f`** (the Demo, or a picked folder) — *in either runtime, Flask running or not* | worker message `{kind: 'orientation'}` → `siteOrientationHistogram` | [orientation.js](../../web_app/frontend/src/workers/orientation.js) |

The switch is in `requestPca(kind, params)`
([useSiteCloud.js](../../web_app/frontend/src/useSiteCloud.js)): if a local `.rmc6f` text is loaded,
the shared worker answers; otherwise axios hits
`{ sites: '/api/pca/sites', orientation: '/api/pca/orientation' }[kind] ?? '/api/pca/kde'`. So a
developer running the full Flask stack and opening a folder through the file picker is reading
the **JavaScript** numbers, not the Python ones. The Flask
route caches only the *parse* (`cached_site_displacements`, keyed on path + mtime, `maxsize=8`)
and the *tiling* (`goldberg_tiling`, `maxsize=8`) — the histogram itself is recomputed every
request. The browser side caches the *tiling* per worker (an unbounded `Map`) and the *parse* in
the single **app-lifetime** worker `useSiteCloud.js` creates and never terminates — one cache slot,
keyed on a cheap content signature plus the cluster threshold, shared with the PCA Ellipsoid page;
see *PCA Ellipsoid page* → Step 1, "Caching". The practical consequence for this page is that
navigating away and back does not re-parse the configuration, and that alternating between two
runs re-parses on every switch.

**What the two test suites assert.** They are *parallel* suites — the same properties,
independently computed in each language — not a shared-golden parity suite.

| Assertion | [tests/test_orientation.py](../../tests/test_orientation.py) | [orientation.test.js](../../web_app/frontend/src/workers/__tests__/orientation.test.js) |
|---|---|---|
| $C = 10\nu^2+2$, exactly 12 pentagons, all degrees ∈ {5,6} | $\nu \in \{1,2,3,5,8\}$ | same |
| Centres are unit vectors | $\nu = 6$ | $\nu = 6$ |
| $\sum \Omega_m = 4\pi$, all $\Omega_m > 0$ | $\nu \in \{1,4,9\}$ | $\nu = 6$ |
| Area ratio $\max/\min < 2$ | $\nu = 8$ | $\nu = 8$ |
| Neighbour lists symmetric | $\nu = 4$ | $\nu = 4$ |
| Antipode exact + involutive | $\nu = 5$ | $\nu = 5$ |
| Polygon vertices within the cell ($\mathbf{v}\cdot\mathbf{c} > 0.9$) | $\nu = 6$ | — |
| Out-of-range $\nu$ rejected | yes | yes |
| Assignment vs brute-force $\arg\max$ | 4000 dirs, $\nu\in\{2,6,11\}$, > 99.99 % | 2000 dirs, same $\nu$, > 99.99 % |
| Every centre self-assigns | $\nu = 7$ | $\nu = 7$ |
| Unnormalised input accepted | yes | — |
| **Monte-Carlo area cross-check** (assignment Voronoi fractions × $4\pi$ = analytic $\Omega_m$) | 200 000 dirs, $\nu = 3$, rtol 12 % | 60 000 dirs, $\nu = 2$, tol 15 % |
| `recommendedFrequency` at 774 / 300 / 12000 = 2 / 2 / 10 | **pinned** | **pinned (shared verbatim)** |
| Isotropic cloud: $\langle\mathcal{E}\rangle \approx 1$, significance < 1.5, \|anisotropy\| < 0.05 | 20 000 pts, $\nu=8$ | same |
| Lobed cloud peaks on the long axis, $\mathcal{E}_{\text{peak}} > 2$, significance > 3 | yes | yes |
| One-sided cloud: peak > 4× its antipode, $\mathcal{A} > 0.5$ and $> 3\mathcal{A}_{\text{null}}$ | yes | yes |
| Symmetric cloud: $0.2\,\mathcal{A}_{\text{null}} < \mathcal{A} < 2\,\mathcal{A}_{\text{null}}$ | yes | yes |
| Smoothing conserves $\sum M_m$ (9 dp), lowers the peak, keeps the lobe | yes | yes |
| `cellMeanAmplitude` tracks per-direction travel; empty cells = 0, all finite | yes | yes |
| Smoothed `cellMeanAmplitude` stays within the raw amplitude range | yes | — |
| `amplitude2` moves the peak from the many-small axis to the few-large axis | yes | yes |
| Quantile cut drops the right fraction | 0.25 → 750 ± 2 of 1000 | 0.25 → 745–755 |
| `pca` frame puts an oblique lobe on $\pm x$ | yes | yes |
| Auto-$\nu$ grows with $N$ | yes | yes |
| Orientation tensor equals the direct $\frac1N\sum\mathbf{u}\mathbf{u}^{\mathsf T}$; eigenvalues sum to 1 | atol $10^{-9}$ | 9 dp |
| `geometry` toggle; polygon lengths equal `sizes` | yes | yes |
| Bad input rejected (2-column, all-zero, bad weight, bad frame, quantile = 1) | yes | yes |
| Site tagging (`referenceNumber`, `element`, `totalPoints` = 512, lobe on $x$) | synthetic 8³ `.rmc6f` | same synthetic file, JS parser |
| Worker routing `{kind:'orientation'}` returns the histogram shape | — | $\nu=3$ → 92 cells, `browserOrientation === true` |
| HTTP route `/api/pca/orientation` | [test_backend_api.py](../../tests/test_backend_api.py): 200, `cellCount` 162 at $\nu=4$, `totalPoints` 216, no `polygons` when `geometry=false`, lobe on $x$; bad `weight` → **400** | — |

**What is *not* covered, and where the two could diverge.**

- **No cross-engine golden test exists.** Nothing in CI compares Python output to JS output for
  the same input. `recommendedFrequency(774|300|12000)` are the only values pinned in both
  suites. In particular, if the two ever enumerated the icosahedron, its faces, the geodesic
  lattice, or the angular ordering in a different order, `centers[i]` would mean different cells
  in the two runtimes and **both suites would still pass** — every assertion above is
  index-agnostic. (Contrast [autoScale.js](../../web_app/frontend/src/workers/autoScale.js), which
  *is* parity-tested against committed Python goldens.)
- **Measured for this document, they do agree.** Running both engines on identical inputs on the
  dev machine: at $\nu = 4$, `centers`, `areas`, `sizes`, `antipode` and even the *per-cell
  polygon vertex order* match to $\le 9\times 10^{-16}$; on a shared 3000-point cloud at
  $\nu = 6$, `smoothing=2`, in both frames and both weightings, `counts` and `peakCell` are
  bit-identical and every float field agrees to $\le 3\times 10^{-13}$
  (`enhancement`, `zScore`, `density`, `cellMeanAmplitude`, `pcaAxes`, `orientationTensor`,
  `significance`, `antipodalAsymmetry`). This is a verification, not a regression guard.
- **Different eigensolvers.** Python uses LAPACK `numpy.linalg.eigh`; the JS port uses a 3×3
  cyclic **Jacobi** rotation (`jacobiEigenSymmetric`); the sweep budget, the *absolute*
  $10^{-18}$ convergence test and why Jacobi was chosen over a closed-form cubic are in
  *PCA Ellipsoid page* → Step 3.
  They agree to round-off on well-conditioned clouds (measured $\le 10^{-15}$ above), but for a
  **near-degenerate** covariance the axes within the degenerate subspace are arbitrary and the
  two can return different ones. This affects `frame="pca"` (the sphere is then rotated
  differently in the two runtimes) and `orientationAxes`. It does not affect the `cartesian`
  frame at all.
- **Cosmetic-only differences:** Python reports `pentagonCount` from the actual count, JS
  hard-codes 12 (both after the same construction assertion); JS adds `browserOrientation:
  true`, Flask adds `source`; Python's tiling cache is bounded (`maxsize=8`), the JS `Map` is
  unbounded (a user who visits every frequency at $\nu \le 64$ retains all of them).
- **Empty input errors differently.** For a zero-row cloud, Python's
  `np.asarray([])` has `ndim == 1`, so it fails the shape check and raises
  `vectors must be a numeric array with shape (N, 3)`; JS's `vectors.some(...)` passes trivially
  on an empty array, so the run reaches `no displacement vectors survive the amplitude cutoff`
  instead. Same class of failure, different message.
- **One shape divergence.** With `element = "all"` the JS engine stamps `result.element = "all"`
  onto the payload (`else if (element)` — the string is truthy,
  [orientation.js:716-718](../../web_app/frontend/src/workers/orientation.js)) while Python omits the
  key entirely (`elif element not in (None, "", "all")`,
  [orientation.py:734](../../rmc_toolkits/orientation.py)). Both transports normalise `"all"` to
  `null`/`None` before the call, so it is reachable only from a direct library call — but it does
  mean "both produce the same JSON shape" is true of the wire, not of the functions.
- **Two places the tolerances already differ.** The $\sum\Omega_m = 4\pi$ assertion is a mixed
  absolute+relative bound in Python (`np.isclose`'s default `atol = 1e-8` is included) and purely
  relative in JS — see §4.6. And the eigensolvers differ, as above.
- **Deliberately identical, easy to break:** the half-to-even rounding in `recommendedFrequency`
  (§3), the stable tie-break in the largest-remainder lattice rounding, the strict `< 0` face
  winding test (§4.1), the relative adjacency tolerance $10^{-9}\max(\ell, 1)$ (§4.1), the
  first-maximum tie-break in the face `argmax`, the `argmin |n|` choice of tangent reference axis,
  and the 8-round walk cap. Each of these is a place where an "obvious simplification" in one
  engine silently desynchronises the two.

---

### Caveats

1. **Cells are not equal-area** — $-36\%$ to $+18\%$ at $\nu = 10$ (§4.8). This is fully
   compensated in `density`/`enhancement`/`expected` by dividing by each cell's exact
   $\Omega_m$, so nothing downstream is biased. But *never plot `counts` directly*: the raw
   histogram carries a $\approx 35\%$ dark spot at each of the 12 pentagons.
   The module docstring's "±10 %" understates this.
2. **Smoothing introduces a small icosahedral artefact** (§6.4): it conserves total mass exactly
   but is not area-weighted, so an exactly isotropic map reads 0.976–1.060 at the UI default of
   2 passes and 0.903–1.202 at 12. `zScore` and `antipodalAsymmetry` are computed from raw
   counts and are immune.
3. **`zScore` ignores the weight, and mixes smoothed with raw.** Under
   `amplitude`/`amplitude2` the colour shows weighted mass per steradian while the significance
   readout still describes the head count — deliberate (there is no clean Poisson null for a
   weighted sum), but the two must be read as separate statements. Within one readout, note also
   that `peakCell`/`peakEnhancement` are argmax/value over the **smoothed** map while
   `peakZScore` is the **raw** $z$ at that cell, and `emptyFraction` is a raw-count quantity that
   does not shrink under smoothing.
4. **`amplitudeCutoff` under-reports, and `min_amplitude` is unvalidated.** It is
   `max(threshold, 0)`, so with no cutoff requested it reads 0 even though points below
   $10^{-9}$ Å were dropped — `usedPoints` can be less than `totalPoints` with
   `amplitudeCutoff == 0`. A **negative** `min_amplitude` is likewise accepted, has no effect, and
   is reported as 0. And the quantile is taken over all $N_{\text{tot}}$ amplitudes before any
   cut, not over the survivors — so a repeated quantile cut is not idempotent.
5. **The Delaunay guarantee is empirical here.** The greedy walk's correctness and the
   "polygon = Voronoi cell" identity both rest on the geodesic triangulation being the Delaunay
   triangulation of its own vertices. The code does not prove this; the Monte-Carlo area test
   ($\nu = 3$ / $\nu = 2$) validates it at low resolution, and brute-force comparison confirms it
   at $\nu = 6, 10, 24, 64$ (0 mismatches in 20 000 directions each, measured for this document).
   The changelog records only the unqualified phrase "brute-force-verified", without naming
   frequencies. The walk is also hard-capped at 8 rounds, so a hypothetical pathological tiling
   would silently return a near-nearest cell rather than fail.
6. **Boundary ties are resolved arbitrarily.** A direction exactly equidistant from two centres
   lands in whichever the walk reaches first. Measure-zero for real data; it is why the
   brute-force parity tests assert $> 99.99\%$ rather than 100 %.
7. **Over-binning is still possible on purpose — and it silently disables the asymmetry flag.**
   The UI defaults to a fixed $\nu = 10$, not to Auto, so a site with only a few hundred copies is
   over-binned out of the box (1002 cells, $< 1$ point per cell); Auto would have picked $\nu = 1$
   or 2 for the same data (§3). The 2× default smoothing hides this visually. Worse, because
   $\mathcal{A} \le 1$ while $\mathcal{A}_{\text{null}} = \sqrt{C/\pi N}$ is not bounded, the
   `A > 3·null` red flag **cannot fire at all** unless $N > 9C/\pi$ — 2 871 copies at $\nu = 10$
   (§7). A perfectly one-sided 216-point site at the default resolution reports
   $\mathcal{A} = 1.000$ against a floor of 1.215 and is displayed as unremarkable. Always check
   `significance` (should be $\gg 1$) and the hover `z` before believing a lobe, and switch the
   Resolution control to *Auto* — for both the map **and** the antipodal readout.
8. **The orientation tensor is antipodally blind, like the ellipsoid.**
   $\mathbf{u}\mathbf{u}^{\mathsf T} = (-\mathbf{u})(-\mathbf{u})^{\mathsf T}$, so a "cluster"
   in $\mathbf{T}$ means a preferred *axis*, not a preferred *direction*. Only the map itself and
   `antipodalAsymmetry` can tell $+\mathbf{u}$ from $-\mathbf{u}$.
9. **Near-isotropic sites have arbitrary PCA frames** (§8), and the two runtimes can pick
   different ones. Use `frame="cartesian"` when the site is close to isotropic.
10. **This is a visualisation of an RMC ensemble, not a refined model.** The directions come from
    one RMCProfile configuration; a lobe reflects that configuration's atom placements, which are
    constrained by, not uniquely determined by, the data. Two independent runs on the same data
    should be compared before a discrete-hop claim is made.
11. **Element pooling mixes sites — and the *default* mixes elements.** `element=` pools every
    site of that element. Because each site's cloud is separately mean-centred, the pooled cloud
    is meaningful — but if the sites are symmetry-related with *different orientations*, their
    lobes superimpose into a symmetry-averaged pattern. Per-site is the safer read. Worse, with
    **neither** `reference_number` nor `element` supplied (or with `element` in `{"", "all"}`)
    both engines pool **every site of every element** into one cloud, so every readout then mixes
    species with different masses and different $\mathbf{U}$. The UI never does this — it always
    sends a `referenceNumber` — but a script or a hand-made HTTP request gets it by default.
12. **Reconstructed sites are clusters, not crystallographic sites.** On a file with no
    reference-site/cell columns the browser engine rebuilds "sites" by folding every atom into one
    unit cell and single-linkage-clustering same-element atoms within `clusterThreshold`
    (default 1.5 Å, slider 0.4 – 2.5) — the procedure is *PCA Ellipsoid page* → Step 1b. What it
    costs *here*: the displacement cloud is then relative to a **cluster centroid**, a single
    cluster can merge two distinct crystallographic sites (or split a wide disordered shell), and
    moving the slider re-derives every number on the page — the tiling, the counts, the tensor,
    all of it. A merged cluster is the worst case for this page in particular, because a
    superposition of two clouds about a spurious mean produces a bimodal direction map and an
    `antipodalAsymmetry` measured from a site centre that does not exist. The page flags this only
    indirectly, by showing the Cluster control at all. This path exists **only in the browser
    engine** (see "Parameters and defaults").
13. **The tiling is not aligned to the crystal, and the crystal frame may be oblique.** The
    icosahedron of §4.1 is built with fixed 2-fold axes on Cartesian $x$, $y$, $z$, and
    `goldberg_tiling` is cached on $\nu$ alone — so the tiling never rotates with the data or with
    the space group. Symmetry-equivalent lobes therefore land in cells of different area and
    shape and receive different $\Omega_m$, `expected` and $z$; two directions that *must* be
    equivalent by symmetry will not report identical statistics. Separately, for a non-orthogonal
    cell the a/b/c rods the view draws are oblique to each other and to the sphere's own axes, so
    "comparable to $\mathbf{a}$, $\mathbf{b}$, $\mathbf{c}$" means comparable as directions, not
    as an orthonormal coordinate system.
14. **`frame` re-bins rather than relabels.** Switching between `cartesian` and `pca` moves the
    data across a *fixed* tiling, so `counts`, `expected`, `zScore`, `peakCell`, `emptyFraction`,
    `significance` and `antipodalAsymmetry` all change for identical data (Step 2, detail 5).
    Do not compare per-cell numbers across frames. And below four points `frame="pca"` is
    silently the identity, i.e. exactly the cartesian map, with nothing in the payload saying so.


## Displacement Directions — the sphere view, axis views, and site picker

### What this page shows

The **Displacement Directions** tab (nav label in [App.jsx](../../web_app/frontend/src/App.jsx);
internal page id `orientation`) is the direction-counterpart of the **PCA Ellipsoid** tab. Both read
the *same* per-site displacement cloud — for every atom belonging to one reference site,
$\Delta\mathbf r = \mathbf r_{\text{atom}} - \langle \mathbf r\rangle_{\text{site}}$ in Cartesian
Ångström. The PCA page asks *how far and in what shape*; this page throws the length away and keeps
only the unit direction

$$\mathbf u = \frac{\Delta\mathbf r}{|\Delta\mathbf r|},$$

then displays the probability density of $\mathbf u$ over solid angle, hex-binned on a Goldberg
sphere.

This section documents **everything between the engine's JSON and the pixels**: page composition,
the mesh/relief/colour construction, the frames and axis rods, the fixed-angle mini views, and the
PNG export. It does *not* re-derive anything owned elsewhere:

| Owned by | What |
| --- | --- |
| *Displacement Directions — the orientation-distribution engine* | the Goldberg tiling, direction→cell assignment, `enhancement`, `zScore`, smoothing, `antipodalAsymmetry`, the orientation tensor, and the full parameter/parity/caveat tables for the engine |
| *PCA Ellipsoid page* → Steps 1, 1b, 2, 3 | the displacement cloud itself, the coordinates-only site reconstruction, cloud selection, the covariance and its eigen-decomposition (`axes`, `rms`, `uIso`) |
| *PCA Ellipsoid page* → Step 14 | `SiteStructurePanel` — the unit-cell site picker this page mounts as its third panel, shared verbatim with that page |
| *Principal axes in the crystallographic frame* → Steps 3, 5, 9a, 9b | what a canonicalised axis direction means, `unitCellVectors`, the shared `sceneAxes.js` rod builders, and `orthonormalCrystalFrame` |

Quantities are used here as named engine outputs, with their definition
recalled in one line where the display depends on it.

#### Which runtime computes what

| Path | Trigger | Engine | Source of truth |
| --- | --- | --- | --- |
| Browser worker | a **locally loaded run** — the Demo, or a folder picked in the browser — in *both* Flask and static mode | [`workers/orientation.js`](../../web_app/frontend/src/workers/orientation.js) via [`pcaKdeWorker.js`](../../web_app/frontend/src/workers/pcaKdeWorker.js) | JS port of `orientation.py`, kept in sync |
| Flask API | a **typed backend directory** (`?dir=…`), no local file loaded | `GET /api/pca/orientation` in [app.py](../../web_app/backend/app.py) → `site_orientation_histogram` in [orientation.py](../../rmc_toolkits/orientation.py) | Python, source of truth |

The two engines are written to agree cell-for-cell (same construction order, so cell indices match;
`recommendedFrequency` was fixed to round half-to-even like Python's `round()`). Everything in *this*
section — mesh, relief, colours, cameras, picker — is browser-side JavaScript/Three.js in both
modes; the backend never renders anything for this page.

The response carries a boolean `browserOrientation`, set only by the JS engine
([`workers/orientation.js:680`](../../web_app/frontend/src/workers/orientation.js)). The legend under
the sphere prints `· browser` or `· server` from it, so which engine produced the picture is always
visible on screen.

---

### Step 1 — Page composition: what the page owns vs. what the view owns

[`OrientationPage.jsx`](../../web_app/frontend/src/components/OrientationPage.jsx) is the host. It owns
**all user options as React state** and renders them in one top controls bar (the same
`.pca-controls` markup and CSS the PCA Ellipsoid page uses, in
[`PcaKdePage.css`](../../web_app/frontend/src/components/PcaKdePage.css)), grouped into three clusters:

1. **Site and resolution** — site `<select>`, `Resolution` (ν), and a `Cluster` slider that appears
   only when the loaded file had to have its sites reconstructed.
2. **Weighting and cutoff** — `Weight`, `Min |Δr|`, `Smoothing`.
3. **Appearance** — `Amplitude height`, `Colormap`, `Contrast`, `Cell borders`, `Axes`.

Most controls carry an `InfoBadge` popover with a plain-language explanation. Three do not —
**Colormap** ([`OrientationPage.jsx:243`](../../web_app/frontend/src/components/OrientationPage.jsx)),
**Cell borders** (`:270`) and **Axes** (`:275`) are bare `<span className="control-name">` labels,
on the grounds that they are self-explanatory.

[`OrientationView.jsx`](../../web_app/frontend/src/components/OrientationView.jsx) owns **no options**;
they all arrive as props. It owns the fetch, the two WebGL renderers, and the readouts, and it
returns exactly **two** panels: the mini "Axis views" column and the sphere. Its per-panel actions
live in that sphere panel's header, not the controls bar: the **Crystal | PCA** frame toggle,
**Reset view** and **Save**.

The third panel is **not** OrientationView's.
[`SiteStructurePanel.jsx`](../../web_app/frontend/src/components/SiteStructurePanel.jsx) is a sibling
component rendered directly by `OrientationPage`
([`OrientationPage.jsx:308`](../../web_app/frontend/src/components/OrientationPage.jsx)), and it owns
its own header actions: **a b c**, **Reset view**, **Save**, plus a `Loading…` chip while
`loadingSites`.

#### The 3 : 6.5 : 6.5 grid

`OrientationView` returns a **React fragment** — no DOM element of its own — so its two panels
become direct children of the page's grid container, sitting alongside `SiteStructurePanel`:

```
.orient-layout   (CSS grid, gap 0.85rem)
├─ .orient-mini-panel        "Axis views"        minmax(0, 3fr)
├─ .orient-main-panel        the sphere          minmax(0, 6.5fr)
└─ .pca-unitcell-panel       "Site ellipsoids"   minmax(0, 6.5fr)
```

Two CSS details make this work
([`PcaKdePage.css`](../../web_app/frontend/src/components/PcaKdePage.css), the
`/* ── Orientation page ── */` block):

* `.orient-mini-panel, .orient-main-panel, .orient-layout .pca-unitcell-panel { grid-area: auto; }`
  neutralises the *named* grid areas (`viewport`, `unitcell`) those same class names carry on the
  PCA Ellipsoid page, so on this page they flow into columns 1–3 instead.
* The row is **height-clamped**, not flex-grown:

  ```css
  .orient-layout { flex: none; height: clamp(20rem, calc(100vh - 17.75rem), 60rem); }
  ```

  The `17.75rem` subtrahend is the app header + controls bar + footer + page padding. The clamp is
  authoritative (`flex: none` prevents a flex-grow from stretching the row past what fits), so on a
  16:9 monitor the whole page reaches the footer with no page scroll — verified at 1600×900. Below
  `1200px` viewport width the grid collapses to one column with `height: auto` and explicit panel
  orders: sphere (`order: 1`, `clamp(20rem, 55vh, 34rem)`), picker (`order: 2`,
  `clamp(18rem, 45vh, 30rem)`), axis views last (`order: 3`, fixed `12rem`).

> **Comment/code discrepancy.** Three comments (`OrientationView.jsx` lines 12 and 454, and
> `PcaKdePage.css` line 937) say the view's *root renders with `display: contents`*. There is no such
> rule for it — `display: contents` exists only on `.pca-side`, the PCA page's wrapper. `OrientationView`
> returns a fragment, which is strictly better (no element at all). The behaviour described is
> correct; the mechanism named is not.

The page is mounted once per session and kept alive: `App.jsx` renders it under
`visitedPages.orientation` and hides it with `.workspace-page.is-hidden { display: none; }` when
another tab is active. So switching tabs does **not** unmount the WebGL contexts, and returning to
the page shows the last camera and the last result with no refetch. (The unmount path still disposes
every geometry/material and calls `forceContextLoss()` on both renderers — a leak fixed in the
2026-07-24 review pass, when the view was a tab inside the PCA page and really did unmount.)

While hidden, the page's ancestor is `display: none`, so `mount.clientWidth/clientHeight` are 0. The
three renderers handle that in **two different ways**, and only one of them is a real early return:

* `renderMinis` reads `mount.clientWidth`/`clientHeight` **raw** and genuinely bails:
  `if (!width || !height) return;`
  ([`OrientationView.jsx:170-172`](../../web_app/frontend/src/components/OrientationView.jsx)).
* The two `ResizeObserver` handlers (sphere and picker) instead **fall back to the size captured at
  scene construction**: `const w = mount.clientWidth || width` with `width` itself
  `mount.clientWidth || 640` / `clientHeight || 520` for the sphere (`OrientationView.jsx:226-227,
  278-285`) and `|| 260` / `|| 260` for the picker (`SiteStructurePanel.jsx:47-48, 116-123`). Because
  the page is first mounted while *visible*, those construction fallbacks never fire in practice and
  `width`/`height` are real non-zero pixels — which makes the handlers' own
  `if (w === 0 || h === 0) return;` guard **dead code**. While hidden, the handler proceeds and
  re-applies the *initial* dimensions.

The end result is the one that matters — no canvas is ever resized to nothing — but the mechanism is
a fallback, not a bail-out: a hidden page's canvases sit at their construction-time size until the
next resize that happens while they are visible.

---

### Step 2 — Loading a site: the shared `useSiteCloud` hook

[`useSiteCloud.js`](../../web_app/frontend/src/useSiteCloud.js) is the only data plumbing both this page
and the PCA Ellipsoid page use. It is called here as
`useSiteCloud({ directory, localRun, clusterThreshold })` — note `probability` is **not** passed, so
it takes its default `0.5`.

**One worker for the whole app.** `getWorker()` lazily creates a single module-level
`sharedWorker` (`workers/pcaKdeWorker.js`) and *never* tears it down, so this page and the PCA
Ellipsoid page share one parsed model and one single-slot, content-keyed parse cache; the cache
key, its collision risk and the behavioural change from the older per-mount worker are in *PCA
Ellipsoid page* → Step 1, "Caching". What follows from it here: switching tabs re-uses the cached
clouds instead of re-parsing a multi-megabyte configuration, and the parse survives navigating away
from the page entirely.

**Request routing** (`requestPca(kind, params)`):

* If `rmc6fText` is non-null — i.e. `localRun.structureFile.sourceFile` exists and its `.text()` has
  resolved for *that exact* `File` object — the request is posted to the worker with a monotonically
  increasing integer `id`, and a one-shot `message` listener resolves on the matching id. This is the
  path for the Demo run and any browser-picked folder, **in both runtime modes**.
* Otherwise it is an `axios.get` against the Flask API with `{ dir: directory || '.', ...params }`.
  The endpoint map is `{ sites: '/api/pca/sites', orientation: '/api/pca/orientation' }`, defaulting
  to `/api/pca/kde` for any other kind.

The text is tagged with the `File` it came from (`loadedText = { file, text }`) and
`rmc6fText` is only considered current when `loadedText.file === localFile`; on a dataset switch it
is `null` until the new text loads, so a request can never run against the previous model's text.
`ready` is `!localFile || Boolean(rmc6fText)` and gates the fetch.

**Read failure.** If `localFile.text()` rejects, the hook resets `loadedText` to
`{ file: null, text: null }` *and* sets `sitesError` to the literal string
`'Could not read the structure file.'`
([`useSiteCloud.js:58-60`](../../web_app/frontend/src/useSiteCloud.js)). `OrientationPage` renders that
as a `.pca-error-banner` **above the panel grid**
([`OrientationPage.jsx:285`](../../web_app/frontend/src/components/OrientationPage.jsx)). It is the only
error string this page owns; every other message is passed through verbatim from an engine or from
axios.

**Site table.** The sites effect requests `kind: 'sites'` with `{ probability, clusterThreshold }`
and stores the result. Its dependency array is
`[requestPca, localFile, rmc6fText, probability, clusterThreshold]`, and `requestPca` is itself a
`useCallback` on `[rmc6fText, directory]`, so the table is refetched on **five** triggers:

1. mount;
2. a new local `File` being picked;
3. that file's text resolving (`rmc6fText` flipping from `null` to the contents);
4. a `directory` change (the Flask path — via `requestPca`'s identity);
5. a `probability` or `clusterThreshold` change.

Before requesting, the effect early-returns with `setSites(null)` when `localFile && !rmc6fText`
([`useSiteCloud.js:97`](../../web_app/frontend/src/useSiteCloud.js)). A dataset switch therefore blanks
the site table — and with it the picker, the site `<select>` and the sphere panel's header title —
for at least one render while the new text loads.

Selection is sticky, and the "prefer a clean reconstructed site" default that decides the *first*
selection lives in the hook rather than in either page — see *PCA Ellipsoid page* → Step 14,
"Which site is selected".

**Unit cell — computed twice, differently.** `unitCell` (used by the sphere, its rods and the mini
views) is `unitCellVectors(sites.latticeVectors, sites.supercell)` from
[`pcaCrystalFrame.js`](../../web_app/frontend/src/pcaCrystalFrame.js), imported by the hook at
[`useSiteCloud.js:17`](../../web_app/frontend/src/useSiteCloud.js) and called in the `unitCell`
`useMemo` at lines 130–135: the supercell lattice rows divided element-wise by the repeat counts,
giving $\mathbf a,\mathbf b,\mathbf c$ **as rows, in the same Cartesian Å basis** as every
displacement and every PCA axis in the app. That shared basis is why the a/b/c rods can be drawn in
the sphere's frame without any re-referencing. The operation, its `null`-unless-both-length-3 shape
guard, its `Math.abs(supercell[i]) > 0` divisor fallback (which repairs a **zero** repeat count but
passes a **negative** one through signed, flipping that edge vector), and the fact that three
separate copies of this division exist in the frontend are all *Principal axes in the
crystallographic frame* → Step 5.

The consequence specific to this page is that its two 3D panels sit on opposite sides of that
guard. `SiteStructurePanel` does **not** call `unitCellVectors`. It re-derives the cell inline and unguarded,
`const unit = lattice.map((row, i) => row.map((value) => value / supercell[i]))`
([`SiteStructurePanel.jsx:170`](../../web_app/frontend/src/components/SiteStructurePanel.jsx)), after a
guard that only checks `if (!handle || !sites) return;`. So on the same degenerate metadata that the
sphere silently repairs (divisor 1), the picker divides by zero and propagates `Infinity`/`NaN` — a
NaN `edgeLength`, NaN marker scales and NaN camera framing, i.e. an empty picker panel. The two
panels do not agree here; the "falls back to 1" behaviour is `unitCellVectors`' alone.

> **Runtime asymmetry — the `Cluster` control is browser-only.** The fold-and-cluster site
> reconstruction (for old `.rmc6f` files with no reference-site/cell columns) exists only in
> [`workers/pcaKde.js`](../../web_app/frontend/src/workers/pcaKde.js) — the algorithm is *PCA
> Ellipsoid page* → Step 1b; `rmc_toolkits/pca_kde.py` has no
> equivalent and its `/api/pca/sites` payload carries neither `reconstructed` nor `copiesPerCell`.
> Consequences on the Flask path: the `Cluster` slider never renders, the "prefer a clean site"
> default selection degenerates to "first site", the `(count/copiesPerCell)` suffix is missing from
> the site labels, and the `clusterThreshold` query parameter this page sends to
> `/api/pca/orientation` is silently ignored by the endpoint (it is not read). None of this is an
> error — a backend directory's file is expected to have proper site columns — but the two runtimes
> are not feature-identical here.

---

### Step 3 — The request: which options reach the engine, and which never do

`OrientationView`'s fetch effect (`useEffect` with deps
`[requestPca, ready, selectedRef, frequency, weight, frame, smoothing, minQuantile, clusterThreshold]`)
issues:

```js
requestPca('orientation', {
  referenceNumber: selectedRef,
  frequency: frequency === 'auto' ? null : frequency,
  weight, frame, smoothing,
  minAmplitudeQuantile: minQuantile,
  clusterThreshold,
  geometry: true
})
```

So the *engine* options are exactly this effect's dependency list, and the *display* options are
exactly the mesh-rebuild effect's. Printing both side by side:

```js
// fetch effect  (OrientationView.jsx:160)
[requestPca, ready, selectedRef, frequency, weight, frame, smoothing, minQuantile, clusterThreshold]
// mesh-rebuild effect  (OrientationView.jsx:403)
[result, colormap, contrast, relief, showOutline, showAxes, frame, unitCell, renderMinis]
```

* **Refetch** — site, resolution ν, weight, frame, smoothing, `Min |Δr|` quantile, cluster distance.
* **Client-side only, no refetch** — amplitude height (relief), colormap, contrast, cell borders,
  axes. Changing those re-runs only the mesh-rebuild effect.

Note `frame` appears in **both** lists: flipping Crystal | PCA refetches *and* swaps the rod triad.
`unitCell` is in the rebuild list only — a late-arriving cell (the sites table resolving after the
histogram) redraws the rods without re-requesting the map.

The rebuild is **total, not incremental**. It disposes all three groups and re-runs `buildCellMesh`
unconditionally (`OrientationView.jsx:339-403`), so a pure colormap or contrast change re-triangulates
all 3 996 triangles and reallocates two 36 k-float arrays even though only the `color` attribute
changed; toggling `Cell borders` rebuilds the mesh too. That dominates the `computeVertexNormals()`
waste flagged in Step 4 by an order of magnitude.

Three engine parameters are **never exposed by this page**: `minAmplitude` (the absolute Å cutoff —
only the quantile form is wired up), `element` (element-pooled clouds, which both engines support),
and `targetPerCell` (the 12-per-cell over-binning target). Each takes its engine default (0, `null`,
12).

> **There is always an amplitude floor, even at `Min |Δr| = 0`.** The engine takes
> `cutoff = max(threshold, NEGLIGIBLE_AMPLITUDE)` with
> `NEGLIGIBLE_AMPLITUDE = 1e-9` Å ([`workers/orientation.js:26, 512`](../../web_app/frontend/src/workers/orientation.js))
> and the keep test is **strict**, `amplitude[i] > cutoff`. So an atom sitting exactly on its site
> mean has no direction and is dropped, at every setting. `usedPoints` is therefore not in general
> equal to the site's atom count even at the default.

`geometry: true` is always requested, so the cell polygons cross the worker/HTTP boundary on **every**
request even though they depend only on ν. At the default ν = 10 that is 1002 ragged polygons of 5–6
unit vectors ≈ 6000 triples of `Number`s per response.

**Serialisation across the two paths.** On the worker path the params object is posted as-is. On the
Flask path it goes through `axios.get(url, { params })`, which **omits `null`/`undefined` keys
entirely** — so `frequency: null` ("Auto") arrives as an *absent* `frequency` query key, and
`pca_orientation_endpoint` resolves that to `None` → the auto path
(`frequency=int(frequency_raw) if frequency_raw not in (None, "") else None`,
[`app.py:648-652`](../../web_app/backend/app.py)). That silent omission is the only reason the two
runtimes agree on Auto. Likewise `geometry: true` is serialised as the string `"true"` and parsed by
`request.args.get("geometry", "true").lower() in ("1", "true", "yes")`.

Three guards live in this effect:

* The result is stored with `setHover(null)` first, because a hover index from a *previous*, larger
  tiling would be out of range for a smaller one.
* An error whose message matches `/Unknown reference number/i` is swallowed (no error banner, current
  sphere retained). A cluster-distance change transiently renumbers the sites, so the in-flight
  request can legitimately name a site that no longer exists; the sites effect re-selects a valid one
  a moment later. **Every other error reaches the badge and clears the sphere** (`setResult(null)`),
  including the engine's own `'no displacement vectors survive the amplitude cutoff'` (thrown when
  `used < 1`) and its validation errors on `minAmplitudeQuantile ∈ [0,1)`, `weight ∈ WEIGHTS`,
  `frame ∈ FRAMES` and `frequency ∈ [1, 64]`.
* The **third guard is less benign**: the effect opens with
  `if (!ready || selectedRef == null) { setResult(null); return; }`
  ([`OrientationView.jsx:126`](../../web_app/frontend/src/components/OrientationView.jsx)), which clears
  the result but **not** `error` or `loading`. A red badge (or a stuck `Computing…`) from the previous
  request therefore survives over an emptied panel until the next successful load.

Both badges are **overlaid on the canvas**, not in the panel header: `{loading && <div
className="pca-badge">Computing…</div>}` and `{error && <div className="pca-badge
is-error">{error}</div>}` sit inside `<div className="pca-canvas orient-canvas" ref={mountRef}>`
(`OrientationView.jsx:535-537`), absolutely positioned at the viewport's top-left corner
(`.pca-badge`, `PcaKdePage.css:438`). The `<h3>` header holds only the title, the atom count, the
frame toggle, Reset view and Save. (It is the *picker*, `SiteStructurePanel`, that puts a `Loading…`
chip in its header.)

---

### Step 4 — Cell mesh: the typed-array layout

All the pure geometry lives in
[`orientationSphere.js`](../../web_app/frontend/src/orientationSphere.js) — deliberately free of
Three.js and DOM so it is unit-testable
([`orientationSphere.test.js`](../../web_app/frontend/src/__tests__/orientationSphere.test.js)).

> **Read Steps 4–6 together.** The helpers are documented one concern at a time, which is *not* the
> order they run in. The mesh-rebuild effect
> ([`OrientationView.jsx:359-402`](../../web_app/frontend/src/components/OrientationView.jsx)) executes:
> 1. `reliefFactors` → per-cell radial factors (Step 5);
> 2. `vertexRadii` → the `radii` Map (Step 5);
> 3. `buildCellMesh(polygons, enhancement, colormap, vmax, radii, contrast)` (this step for the
>    layout, Step 6 for the `contrast` transfer that runs *inside* it);
> 4. `buildCellOutline` (Step 8);
> 5. `buildAxisRods` (Step 10);
> 6. `renderMinis()` (Step 11).
>
> So `radii` and `contrast` are fully determined *before* the mesh call, even though this step meets
> them only as the bare defaults `radii = null, contrast = 1`.

**Input.** `result.polygons` is a ragged list, one entry per cell, of $n_c \in \{5, 6\}$ unit vectors
wound counter-clockwise as seen from outside. (The engines pad internally to 6 by repeating the last
vertex; both slice back to the true size before emitting, so no degenerate vertex reaches the
renderer.) `result.enhancement` is one dimensionless number per cell.

**`buildCellMesh(polygons, values, colormap, vmax, radii = null, contrast = 1)`** emits a
**non-indexed** triangle soup — there is no index buffer:

| Array | Type | Length | Layout |
| --- | --- | --- | --- |
| `positions` | `Float32Array` | $9T$ | per triangle: $x_0y_0z_0\,x_1y_1z_1\,x_2y_2z_2$ |
| `colors` | `Float32Array` | $9T$ | same layout, RGB in $[0,1]$ |
| `triangleCell` | `Int32Array` | $T$ | triangle index → cell index |
| `triangles` | `Number` (scalar) | — | $T$ itself — the loop bound for iterating the mesh, and the value the unit test asserts against $\sum_c(n_c-2)$ |

with the triangle count

$$T \;=\; \sum_c (n_c - 2) \;=\; \Big(\sum_c n_c\Big) - 2C \;=\; (6C - 12) - 2C \;=\; 4C - 12,$$

using $C = 10\nu^2 + 2$ cells of which exactly 12 are pentagons. At the default $\nu = 10$:
$C = 1002$ (990 hexagons + 12 pentagons), $T = 3996$ triangles, 11 988 emitted vertices,
$|positions| = |colors| = 35\,964$ floats.

**Fan triangulation.** For a polygon $P = [v_0, \dots, v_{n-1}]$ the fan is

$$(v_0, v_i, v_{i+1}),\qquad i = 1,\dots,n-2,$$

i.e. always anchored on $v_0$. Because the polygons are convex spherical cells and CCW from outside,
every fan triangle inherits that winding, and the mesh material uses `THREE.FrontSide`.

**Vertices are duplicated per cell on purpose.** Each of a cell's $3(n_c-2)$ emitted vertices is
written the *same* RGB, so the cell is a single flat colour with no interpolation across it — the
crisp hex-bin look — without needing flat-shading qualifiers or per-face materials. The test
`flat-colors each cell` asserts exactly this, and that the colour set covers all `cellCount` cells.

**Consumption** (`OrientationView.jsx`, mesh-rebuild effect):

```js
geometry.setAttribute('position', new THREE.BufferAttribute(mesh.positions, 3));
geometry.setAttribute('color',    new THREE.BufferAttribute(mesh.colors, 3));
geometry.computeVertexNormals();
material = new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.FrontSide });
sphere.userData.triangleCell = mesh.triangleCell;
```

The scene contains **no lights** and the material is `MeshBasicMaterial` (unlit), so the sphere is
pure colour with no shading gradient — every pixel of a cell is exactly its colormap sample. Note
`computeVertexNormals()` is called but its output is never used by an unlit material; it is dead work
(one extra pass over 36 k floats per rebuild).

---

### Step 5 — Amplitude relief: the exact radius map

This is the page's second, independent visual channel: **colour says how *often* atoms move that
way; radius says how *far* they move when they do.**

**Input.** `result.cellMeanAmplitude[c]` $= A_c$ (Å) — the mean $|\Delta\mathbf r|$ of the atoms that
landed in cell $c$, 0 for an empty cell — and `result.meanAmplitude` $= \bar A$ (Å), the mean
$|\Delta\mathbf r|$ over **all surviving** vectors. Both come from the engine. When smoothing is on,
$A_c$ is smoothed as a *ratio of smoothed sums* (numerator and denominator diffused separately) so it
stays the mean of the same smoothed population the colour shows; $\bar A$ is never smoothed.

**Per-cell radial factor** — `reliefFactors(cellMeanAmplitude, globalMean, relief, clampLow = 0.3, clampHigh = 2.2)`:

$$f_c \;=\;
\begin{cases}
1, & \neg(A_c > 0) \ \text{ or } \ \neg(\bar A > 0),\\[4pt]
\operatorname{clip}\!\left(1 + \rho\left(\dfrac{A_c}{\bar A} - 1\right),\; 0.3,\; 2.2\right), & \text{otherwise,}
\end{cases}$$

The guard is written in the code as `if (!(mean > 0) || !(globalMean > 0)) return 1;`
([`orientationSphere.js:46`](../../web_app/frontend/src/orientationSphere.js)) — a **negated positivity
test**, not a sign test. That distinction is load-bearing: it catches `NaN`, `undefined` and `null`
as well as $\le 0$, so a *missing* `cellMeanAmplitude` entry (a truncated or malformed response)
degrades to the neutral radius 1 exactly like an empty cell, instead of producing a NaN vertex.

Here $\rho$ is the **Amplitude height** slider, dimensionless, range $[0, 1]$ in steps of 0.05,
**default 0.5**. The map has three properties worth stating explicitly:

* $A_c = \bar A \Rightarrow f_c = 1$ for any $\rho$ — a cell whose movers travel the site-average
  distance sits exactly on the unit sphere. The unit sphere is the *neutral* surface, not the floor.
* $\rho = 0 \Rightarrow f_c \equiv 1$ — a perfect sphere, verified in the unit test.
* $\rho = 1$ makes the radius literally $A_c/\bar A$, clamped to $[0.3, 2.2]$. An **empty cell is
  neutral (1), not 0** — a hole in the direction distribution is shown by colour, never by a dent.

The clamp bounds 0.3 and 2.2 are defaults of `reliefFactors`; `OrientationView` calls it with only
three arguments, so they always apply.

**Per-vertex radius** — `vertexRadii(polygons, factors)`. A per-*cell* radius would tear the surface
at every shared edge, so each polygon vertex takes the **arithmetic mean of the factors of the cells
incident on it**:

$$r(v) \;=\; \frac{1}{m_v}\sum_{c \,\ni\, v} f_c ,$$

where $m_v$ is the number of cells that listed $v$ (3 for a Goldberg dual). Vertex identity is lost
across the worker/JSON boundary, so vertices are keyed by rounded coordinates:

```js
const vertexKey = (v) => `${v[0].toFixed(9)},${v[1].toFixed(9)},${v[2].toFixed(9)}`;
```

Nine decimals is safe on both counts: bit-identical sources agree exactly, and adjacent circumcentres
are separated by roughly $0.7/\nu$ (radians on the unit sphere, hence the same in the Cartesian
components) — ~$3\times10^{-2}$ at the finest tiling this page can actually request (ν = 24) and
~$10^{-2}$ at the engine's own bound `MAX_FREQUENCY = 64`, which the UI never offers. Either way the
separation is seven orders of magnitude above the rounding, so it cannot merge two different
vertices. The test `produces a crack-free relief surface` asserts that every emitted
copy of a given direction lands at an identical radius to 9 decimal places, while the radii still
vary by more than 0.2 across the sphere.

`buildCellMesh` then scales each fan vertex by `radii.get(vertexKey(v)) ?? 1` — the `?? 1` means a
vertex missing from the map degrades to the unit sphere rather than collapsing to the origin.

The call site is guarded: relief is computed only when `relief > 0 && result.cellMeanAmplitude`,
otherwise `radii` is `null` and every scale is exactly 1.

> **Approximation, stated honestly.** Because a vertex radius is the *mean* of three cell factors,
> the surface never actually attains $f_c$ anywhere inside cell $c$ unless all three of a vertex's
> cells agree. A single isolated long-amplitude cell therefore bulges visibly **less** than its
> factor suggests, and a sharp amplitude step is rendered as a ramp one cell wide. The relief is a
> continuous piecewise-linear interpolant of the cell factors, not a stepped extrusion of them. The
> colour channel has no such smoothing — cells stay flat and exact.

---

### Step 6 — Colour: which scalar, and the contrast transfer

**The mapped scalar is `enhancement`**, not density and not the z-score:

$$\text{enhancement}_c \;=\; 4\pi \, \rho_c , \qquad \rho_c = \frac{m_c}{\big(\sum_k m_k\big)\,\Omega_c},$$

with $m_c$ the cell's (possibly smoothed) weighted mass and $\Omega_c$ its exact solid angle in
steradian. It is dimensionless, **exactly 1 everywhere for an isotropic site**, and reads directly as
"this direction is 1.8× more likely than chance". Definition and derivation belong to the engine
section; what matters here is that the display normalisation is built around the value 1.

**Normalisation range: $[0, v_{\max}]$, one-sided, NOT symmetric about 1.** `vmax` is
`result.vmax` = $\max_c$ enhancement over the tiling. **At contrast $\gamma \ge 1$**, zero maps to the
bottom of the colormap and the peak cell to the top (asserted by the test
`maps the peak cell to the top of the colormap`, which calls
`colorCoordinate(result.vmax, result.vmax)` at the *default* contrast and so does not cover
$\gamma < 1$). The isotropic level 1 sits at the *fraction* $1/v_{\max}$ along the bar — near the
bottom whenever the map is strongly peaked. There is no diverging/symmetric mode; a `seismic`
colormap is offered but **at $\gamma = 1$** its white midpoint lands at enhancement $v_{\max}/2$,
not at 1.

Both statements fail below $\gamma = 1$, which the Contrast slider reaches (min 0.5): the transfer
compresses toward the isotropic pivot, so neither end of the colormap is used. With $v_{\max} = 4$
and $\gamma = 0.5$ the peak maps to $t = 0.625$ and enhancement 0 maps to $t = 0.125$ — the bar's
extremes are simply never painted, and seismic's white sits at $v = 3.0$.

**The contrast transfer** — `colorCoordinate(value, vmax, contrast = 1)`:

$$p \;=\; \operatorname{clip}_{[0,1]}\!\left(\frac{1}{v_{\max}}\right),
\qquad
t(v) \;=\; \operatorname{clip}_{[0,1]}\!\Big(p + \gamma\big(\operatorname{clip}_{[0,1]}(v/v_{\max}) - p\big)\Big),$$

The early return is `if (!(vmax > 0)) return 0;` — again a negated positivity test, so a **missing or
`NaN` `vmax`** paints the whole sphere the colormap floor rather than producing NaN colours. (There
is no equivalent guard on `value`: see the LUT note below.) Here $t \in [0,1]$ is the colormap
coordinate, $\gamma$ is the
**Contrast** slider (dimensionless, range 0.5–3.0 step 0.1, **default 1**), and $p$ is the colour
coordinate of the isotropic level — the physically meaningful pivot.

Properties (all unit-tested):

* $\gamma = 1$ reduces to the plain linear $v/v_{\max}$, so the control is backward-compatible and
  callers that omit it are unchanged.
* $v = 1$ is the **fixed point** for every $\gamma$: the isotropic tone never moves.
* $\gamma > 1$ pushes super-isotropic cells up and sub-isotropic cells down (faint lobes *and*
  depletions gain contrast); $\gamma < 1$ flattens both toward the isotropic tone.
* The outer `clip` means high contrast saturates rather than wrapping.

Contrast is **purely a colour gain**: it does not touch the relief radii, the z-scores, or any
readout in the summary strip.

**The LUT** is [`colormaps.js`](../../web_app/frontend/src/colormaps.js). Available maps: `viridis`
(default), `magma`, `seismic`, `reds`, `greys` — 5 anchor RGB triples each. `buildLut(name)` walks
them piecewise-linearly into a 256-entry `Uint8ClampedArray` (cached per name in `lutCache`), with
$S = \text{anchors.length} - 1 = 4$ segments; for LUT entry $i$, $t = i/255$,

$$u = tS,\quad \ell = \min(S, \lfloor u \rfloor),\quad h = \min(S, \ell+1),\quad
\text{lut}[i] = A_\ell + (A_h - A_\ell)(u - \ell).$$

Three details the formula hides:

* **Ties-to-even rounding, not `Math.round`.** Assigning a float into a `Uint8ClampedArray` clamps to $[0,255]$ and
  rounds half-to-even — so LUT entries can differ by 1 from a naive `round()` of the interpolant.
* **Unknown-name fallback.** `ANCHORS[name] || ANCHORS.viridis` — an unrecognised colormap silently
  renders as viridis rather than throwing or rendering black.
* **Index clamp.** `sampleColormap(name, value)` is
  `Math.max(0, Math.min(255, Math.round(value * 255))) * 3`, so an out-of-range coordinate saturates
  at an end of the map rather than reading out of bounds.

The clamp does **not** rescue `NaN`: `Math.round(NaN)`, `Math.min`/`Math.max` of `NaN` are all `NaN`,
so `lut[NaN]` is `undefined` and the returned triple is `[undefined, undefined, undefined]` → NaN
colour channels in the attribute buffer. The only way to reach that here is an `undefined`
`values[cell]` (a short or ragged `enhancement` array), since `colorCoordinate` already absorbs a bad
`vmax`. Nothing guards it. `buildCellMesh` divides the sampled bytes by 255 into the float colour
attribute.

---

### Step 7 — The colorbar and the summary strip

#### 7a. The colorbar

`colorbarGradient(colormap, stops = 24, contrast = 1, vmax = 1)` returns a CSS
`linear-gradient(to right, …)` string used as the `background` of `.orient-colorbar`.
`OrientationView` calls it as `colorbarGradient(colormap, 24, contrast, result?.vmax ?? 1)` inside a
`useMemo`.

For $s = 0 \dots 24$ (**25 colour stops**):

$$\text{base} = \frac{s}{24}, \qquad
\text{colour} = \mathrm{LUT}\big[\operatorname{round}\!\big(t(\text{base}\cdot v_{\max})\cdot 255\big)\big],\qquad
\text{position} = (100\cdot \text{base})\%\ \text{(1 decimal)} .$$

The crucial property is that the bar's geometric axis is **linear in enhancement from 0 to
$v_{\max}$**, while each stop's *colour* is painted through the **identical** `colorCoordinate`
transfer the cells use. So the bar is an exact legend at any contrast — the test
`tracks the cell color transfer under contrast` re-derives every stop's RGB from
`sampleColormap(colormap, colorCoordinate(percent/100 · vmax, vmax, contrast))` and asserts byte
equality.

The row is labelled `0` on the left and `<vmax> × isotropic` on the right (1 decimal). When
$v_{\max} > 1$ a 2 px vertical tick (`.orient-colorbar-marker`, `title="1× = isotropic"`) is placed
at

$$\text{left} = \left(\frac{100}{v_{\max}}\right)\!.\mathtt{toFixed(2)}\ \%,$$

i.e. the CSS `left` is written to 2 decimal places, marking the isotropic level. When
$v_{\max} \le 1$ — a map with no cell above chance anywhere — the tick is omitted rather than pinned
to the right edge.

**One formatter behind every printed number.** The right-hand colorbar label, all four summary
statistics and three of the four tooltip lines go through the module-local

```js
const numberFormat = (value, digits = 2) => (Number.isFinite(value) ? value.toFixed(digits) : '—');
```

([`OrientationView.jsx:50-51`](../../web_app/frontend/src/components/OrientationView.jsx); the host page
has its own copy with `digits = 4` for the site labels). So **any non-finite value — `NaN`,
`Infinity`, a missing field — prints as an em-dash `—`, never as `NaN`**. The tables below give only
the digit counts; the `—` fallback applies to every one of them. The two exceptions are noted where
they occur: the tooltip's atom count is printed raw, and the tick position uses `.toFixed(2)`
directly.

#### 7b. The summary strip

Directly under the colorbar, `.orient-summary` prints four statistics. **All four are engine
outputs — none is recomputed in the browser** — and each is recalled here with its definition and the
formatting the UI applies. $N$ is the number of surviving displacement vectors (`usedPoints`), $C$
the cell count, $n_c$ the **raw** (unsmoothed) count in cell $c$, and $e_c = N\,\Omega_c/4\pi$ its
isotropic expectation.

| Readout | Engine field | Definition | Format |
| --- | --- | --- | --- |
| `peak N.NN× at [x, y, z] (z = N.N)` | `peakEnhancement`, `peakDirection`, `peakZScore` | the cell with the largest `enhancement`; its centre direction and its raw-count Poisson $z$ | 2 dp, direction 2 dp, $z$ 1 dp |
| `anisotropy N.NN` | `orientationAnisotropy` | $3\lambda_1 - 1$ of the orientation tensor $T = \langle \mathbf u\mathbf u^{\mathsf T}\rangle$ (weighted by the selected weight), $\lambda_1$ its largest eigenvalue. $T = I/3$ for a uniform sphere, so the value is **0 for isotropic, 2 for a perfect single axis**. Computed from the vectors, not the bins, so it is resolution-independent. | 2 dp |
| `± asymmetry N.NN (noise floor N.NN)` | `antipodalAsymmetry`, `antipodalAsymmetryNull` | $\dfrac{1}{2N}\sum_c \lvert n_c - n_{\bar c}\rvert$ over cells, $\bar c$ the exact antipodal cell — equivalently $\sum_{\text{pairs}}\lvert n(\mathbf u) - n(-\mathbf u)\rvert / N$: **0 for an inversion-symmetric cloud, 1 for a fully one-sided one**. The floor is the level pure Poisson noise alone produces, $\sqrt{C/(\pi N)}$, from $\mathbb E\lvert X-Y\rvert \approx 2\sqrt{m/\pi}$ for two iid Poisson($m$) cells. | both 2 dp |
| `map significance N.Nσ` | `significance` | $\sqrt{\big\langle z_c^2 \big\rangle_c}$ with $z_c = (n_c - e_c)/\sqrt{e_c}$ — the RMS per-cell Poisson $z$ against the isotropic null. **≈1 means the whole map is consistent with counting noise**; well above 1 means real directional structure. | 1 dp |

The only *client-side* logic in the strip is the red flag:

```js
const asymmetrySignificant = result
  ? result.antipodalAsymmetry > 3 * result.antipodalAsymmetryNull
  : false;
```

which adds `.is-flagged` (colour `var(--danger)`) to the asymmetry stat. The factor **3** is a UI
constant ([`OrientationView.jsx:447`](../../web_app/frontend/src/components/OrientationView.jsx)), not
an engine threshold. A genuine +u/−u imbalance is the headline finding of this page — it is exactly
the physics (static off-centring, odd-order anharmonicity) that a covariance ellipsoid is
structurally blind to — which is why it gets the only coloured state in the strip. See Caveat 1 for
when that flag cannot fire at all.

The strip does **not** show the resolution actually used, the number of vectors that survived the
amplitude cut, or the tensor's eigenvectors; see Step 14.

---

### Step 8 — Cell outlines

`buildCellOutline(polygons, inflate = 1.002, radii = null)` builds a flat `Float32Array` of
`LineSegments` positions: for every cell, every edge $(v_i, v_{(i+1) \bmod n})$ contributes 6 floats.

$$E \;=\; \sum_c n_c \;=\; 6C - 12 \quad(\text{6000 segments at } \nu = 10,\ 36\,000 \text{ floats}).$$

Each endpoint is scaled by `(radii.get(key) ?? 1) × inflate`, so the borders **follow the relief
surface** and sit 0.2 % above it — enough to defeat z-fighting with the mesh, small enough to read as
coincident. Interior edges are shared by two cells and are therefore emitted **twice**; this is
deliberate (harmless for a line overlay, far simpler than deduplication) and doubles the segment
count over the true edge set.

`OrientationView` draws them with
`LineBasicMaterial({ color: 0x10151c, transparent: true, opacity: 0.35 })` when the **Cell borders**
switch is on (**default on**). The unit tests pin both the unscaled case (every outline vertex at
radius `inflate`) and the relief case (every vertex at `1.5 × 1.01` for uniform factors 1.5).

---

### Step 9 — Hover picking and the tooltip

Picking is a **raycast against the mesh**, not a nearest-centre lookup on the cell table. In the
scene-setup effect:

```js
pointer.x =  ((e.clientX - rect.left) / rect.width)  * 2 - 1;
pointer.y = -((e.clientY - rect.top)  / rect.height) * 2 + 1;
raycaster.setFromCamera(pointer, camera);
const hit  = raycaster.intersectObjects(meshGroup.children, false)[0];
const cell = hit ? hit.object.userData.triangleCell[hit.faceIndex] : null;
```

Only `meshGroup` is tested (outlines and axis rods are not pickable), non-recursively, and only the
nearest hit is used. Three.js reports the intersected `faceIndex`; because the geometry is
non-indexed and built in fan order, `faceIndex` **is** the triangle ordinal, so the `triangleCell`
lookup is O(1) and exact — including on the relief surface, since the raycast runs against the
displaced positions.

An **effect-scoped `lastCell`** suppresses React state churn when the pointer stays inside the same
cell. It is a plain `let lastCell = null;` declared inside the scene-setup `useEffect` body
([`OrientationView.jsx:251`](../../web_app/frontend/src/components/OrientationView.jsx)) and captured by
the `onPointerMove`/`onPointerLeave` closures — *not* module state. That matters: module scope would
be shared by every mounted instance and would survive a remount, whereas this deduplication state is
per-scene and resets whenever the scene is rebuilt. Note the `else if (cell != null)` branch still
calls `setHover` on every pointer move over the sphere so the tooltip can follow the cursor. Leaving
the canvas clears both.

The tooltip (`.orient-tooltip`, positioned at `hover.x + 14`, `hover.y + 12` in canvas-local px)
shows four lines:

| Line | Value | Provenance |
| --- | --- | --- |
| `[x, y, z]` | `result.centers[cell]`, 2 dp, via `formatDirection` | the **cell centre**, not the ray direction |
| `N.NN× isotropic` | `result.enhancement[cell]`, 2 dp | smoothed if smoothing > 0 |
| `n atoms · z = ±N.N` | `result.counts[cell]` **printed raw — no `numberFormat`, so no `—` fallback and no rounding**; `result.zScore[cell]` 1 dp | **raw, unsmoothed counts** |
| `⟨\|Δr\|⟩ = N.NNN Å` | `result.cellMeanAmplitude[cell]`, 3 dp; shown only when > 0 | smoothed if smoothing > 0 |

Two honest notes on this readout:

* **Mixed provenance in one tooltip.** The colour/enhancement and ⟨|Δr|⟩ lines reflect the smoothed
  field; the count and z lines are the raw counts (by design — smoothing correlates neighbours and
  would overstate the evidence). With the default 2× smoothing they refer to different populations of
  the same cell.
* **The direction is a Cartesian cell centre, not a crystallographic index.** In the Crystal frame it
  is $[x,y,z]$ in the shared Cartesian Å basis; in the PCA frame it is
  $[\text{PC1}, \text{PC2}, \text{PC3}]$ components. Neither is converted to $[uvw]$ — and neither
  is anything on the PCA
  Ellipsoid page: `principalAxisOrientation`, which computes `crystalDirection`, direction cosines
  and the dominant crystal axis, is implemented and unit-tested but imported by nothing except its
  own test (*Principal axes in the crystallographic frame* → Step 10). No page in the app prints a
  $[uvw]$ today. Its angular precision is the cell size: at ν = 10 a cell spans
  $\Omega = 4\pi/1002 = 0.01254$ sr, an equivalent cap half-angle of
  $\arccos(1 - \Omega/2\pi) \approx 0.063\ \text{rad} = 3.6°$; at ν = 24 it is 1.5°.

A second guard sits between the raycast and the render:

```js
const hoverCell = hover && result && hover.cell < result.cellCount ? hover.cell : null;
```

so a hover captured against a larger tiling can never index past a newly-shrunk result.

---

### Step 10 — Frames and axis rods

#### The Crystal | PCA toggle

The toggle in the sphere panel header sets the `frame` prop, which is an **engine** parameter, so
flipping it refetches:

* **Crystal** (`frame: 'cartesian'`, the default) — directions stay in the shared Cartesian Å frame
  the whole app uses, the same basis in which $\mathbf a,\mathbf b,\mathbf c$ are expressed.
* **PCA** (`frame: 'pca'`) — the engine rotates each unit direction into the site's own principal
  frame, $\mathbf u' = P\,\mathbf u$ with $P$ the principal-axis matrix (rows = axes). After the
  rotation PC1 is $+x$, PC2 is $+y$, PC3 is $+z$, which makes different sites directly comparable.
  The fit uses **all** vectors, before any amplitude cut, and goes through the *same*
  `_eigen_decomposition` the PCA Ellipsoid page's axes come from, so PC1 here is PC1 there — see the
  **engine** section, Step 2 (which frame each output is in, and the measured $5\times10^{-16}$
  round-off between the two covariance estimators) and its Step 8 (why sharing the routine is
  mandatory). What a canonicalised sign does and does not mean — including that the
  largest-component-positive rule holds for PC1/PC2 but not PC3, and that a principal axis is an
  undirected line either way — is *Principal axes in the crystallographic frame* → Step 3.
  Below four vectors the engine never fits at all and `frame: 'pca'` is silently the identity; see
  the `pcaAxes` fallback note below.

#### Rod geometry

`buildAxisRods(axes, colors, radius, length)` (local to `OrientationView`) draws, for each axis row,
**two** cylinder meshes — one for $+\hat{\mathbf d}$ and one for $-\hat{\mathbf d}$ — because a
direction map has no preferred end of an axis:

```js
dir.normalize();                                          // rows are used for orientation only
rod = Mesh(CylinderGeometry(radius, radius, length, 10),
           MeshBasicMaterial({ color, transparent: true, opacity: 0.9 }));
rod.position = dir * (sign * length / 2);
rod.quaternion.setFromUnitVectors(TRIAD_UP, dir * sign);   // TRIAD_UP = (0,1,0)
```

Each rod therefore spans from the sphere centre out to $\pm\,length$, giving a full through-axis of
$2\,length$. Rows with $|\mathbf d|^2 < 10^{-12}$ are skipped.

**Three rod builders, three behaviours.** This page draws rods from three different functions, and
they are not interchangeable:

| Builder | Where | Signs drawn | Span | Degenerate-row guard |
| --- | --- | --- | --- | --- |
| `buildAxisRods` | local to `OrientationView`, the sphere's triad | **both** $+\hat{\mathbf d}$ and $-\hat{\mathbf d}$ (2 meshes/axis) | $[-L,\,+L]$ through the centre | `lengthSq() < 1e-12` → skip that axis |
| `buildCrystalAxes` | [`sceneAxes.js:49-64`](../../web_app/frontend/src/components/sceneAxes.js), the picker's a/b/c gizmo | **one-sided**, $+$ only | $[\text{origin},\ \text{origin} + L\hat{\mathbf d}]$ | `lengthSq() < 1e-12` → skip |
| `buildAxisTriad` | [`sceneAxes.js:30-43`](../../web_app/frontend/src/components/sceneAxes.js), the picker's selection triad | **one-sided**, $+$ only | same | **none** — a zero row normalises to `NaN`, giving a NaN quaternion and a NaN position |

All three `normalize()` their input rows, so only the **orientation** of an axis is carried; lengths
are set by the caller. All three position at $\text{origin} + \hat{\mathbf d}\,L/2$ and orient with
`quaternion.setFromUnitVectors(TRIAD_UP, dir)`. `buildCrystalAxes` and `buildAxisTriad` are the
shared builders — the PCA Ellipsoid page's viewport calls the same two, and their signatures and
the "normalised rods carry orientation, never $|a|,|b|,|c|$" consequence are *Principal axes in
the crystallographic frame* → Step 9a. `buildAxisRods` is local to `OrientationView` and exists
only because a direction map has no preferred end of an axis.

`TRIAD_UP` and both palettes are shared via
[`sceneAxes.js`](../../web_app/frontend/src/components/sceneAxes.js) so PC1/2/3 and a/b/c read as the
same objects on every panel in the app:

| Palette | Values | Used for |
| --- | --- | --- |
| `TRIAD_COLORS` / `PC_CSS_COLORS` | `0xd64545` / `0x3fa34d` / `0x3f7fd6` (red, green, blue) | PC1, PC2, PC3 |
| `CELL_AXIS_COLORS` / `CELL_AXIS_CSS` | `0xe0a419` / `0x18a3a0` / `0xb15ad8` (gold, teal, orchid) | a, b, c |
| `SELECTED_ATOM_COLOR` | `0xff7a1a` | selected site in the picker |

The crystal palette is deliberately unrelated to the PC tricolor so the two frames never read as the
same thing.

#### Only the selected frame's rods are drawn

When the **Axes** switch is on (**default on**), exactly one triad appears:

| `frame` | `unitCell` | Rods drawn | Colors | radius | length |
| --- | --- | --- | --- | --- | --- |
| `'pca'` | any | identity $x,y,z$ (= PC1/2/3 after the rotation) | `TRIAD_COLORS` | 0.006 | 2.7 |
| `'cartesian'` | present | `unitCell` rows $\mathbf a,\mathbf b,\mathbf c$, normalised | `CELL_AXIS_COLORS` | 0.0045 | 2.35 |
| `'cartesian'` | absent | `result.pcaAxes \|\| []` | `TRIAD_COLORS` | 0.006 | 2.7 |

The third row has two fallbacks stacked on it:

* The call is `buildAxisRods(result.pcaAxes || [], …)`, so a response without `pcaAxes` draws **no
  triad at all** while the legend still prints PC1/PC2/PC3 swatches.
* `pcaAxes` is not always the site's principal axes. The engine initialises it to the **identity** and
  only replaces it when the site has at least 4 vectors
  (`let pcaAxes = [[1,0,0],[0,1,0],[0,0,1]]; if (totalPoints >= 4) pcaAxes = eigenDecomposition(…)`,
  [`workers/orientation.js:503-506`](../../web_app/frontend/src/workers/orientation.js)). For a
  1–3-atom site the rods are therefore the **world** axes wearing the PC tricolor.

The legend under the sphere follows the same rule for its swatches: a/b/c only when
`frame === 'cartesian' && unitCell`, otherwise PC1/PC2/PC3. Note that this choice is driven by the
**`frame` prop and `unitCell`, not by the response** (`result.frame` is returned and never read), so
during an in-flight refetch after a frame flip the swatches can already describe a frame the
displayed sphere was not computed in.

The legend's **trailing note**, by contrast, *is* response-derived, so it always describes what was
actually computed: `direction distribution` for `result.weight === 'count'`, else
`weighted by |Δr|` / `weighted by |Δr|²`; then ` · smoothed N×` when `result.smoothing > 0`; then
` · browser` or ` · server` from `result.browserOrientation`.

> **The oblique-cell approximation.** The a/b/c rods carry only the cell **orientation**: each row is
> normalised, so all three rods have equal on-screen length regardless of $|a|, |b|, |c|$, and for a
> non-orthogonal cell they are correctly non-perpendicular on screen. Nothing on this page
> orthonormalises the crystal frame — unlike the PCA Ellipsoid page's viewport, whose shadow box is
> built on the Gram–Schmidt frame `orthonormalCrystalFrame` (*Principal axes in the crystallographic
> frame* → Step 9b, including its unguarded $a\parallel b$ collapse). The consequence here is that
> the three "axis views" of an oblique cell are **not** an orthogonal set of projections (see
> Step 11).

> **Inconsistency worth flagging.** In the third row above (Crystal frame with no cell metadata) the
> sphere and the legend show **PC axes in the PC tricolor**, but the Axis-views mini panel labels its
> three panes `x`, `y`, `z` in grey and snaps down the identity axes — because `miniAxes` falls back to
> the identity while the rod builder falls back to `result.pcaAxes`. The mini views and the rods
> disagree about what the three axes are in that one case. The common case (cell metadata present) is
> consistent.
>
> The sub-case **inverts** the blame: when the site has fewer than 4 vectors, `pcaAxes` *is* the
> identity, so the mini panes' `x, y, z` labels are then correct and it is the rods' PC tricolor that
> is misleading.

---

### Step 11 — The "Axis views" mini panel

The left panel is **one extra `WebGLRenderer`** (a second GL context) rendering the **same scene**
object through three scissored viewports. Cross-renderer scene sharing is safe because Three.js
uploads GPU resources per renderer.

**Which three axes.** `miniAxes` (a `useMemo` on `[frame, unitCell]`) picks:

| Condition | `axes` | `labels` | `colors` |
| --- | --- | --- | --- |
| `frame === 'pca'` | identity | `PC1, PC2, PC3` | `PC_CSS_COLORS` |
| `frame === 'cartesian'` and `unitCell` | `unitCell` rows | `a, b, c` | `CELL_AXIS_CSS` |
| `frame === 'cartesian'`, no `unitCell` | identity | `x, y, z` | grey `#8a97a8` ×3 |

**Camera placement.** `renderMinis` reuses **one** static camera for all three panes. It is
constructed as `new THREE.PerspectiveCamera(45, 1, 0.01, 100)` — FOV 45°, **aspect 1**, near 0.01,
far 100 ([`OrientationView.jsx:301`](../../web_app/frontend/src/components/OrientationView.jsx)) — and
only receives its real aspect `width / paneHeight` inside the render loop, so the constructor's
aspect never reaches the screen. For pane $i \in \{0,1,2\}$:

$$\hat{\mathbf d}_i = \frac{\mathbf{axes}_i}{|\mathbf{axes}_i|},\qquad
\mathbf{up}_i^{\text{raw}} = \mathbf{axes}_{\,\mathrm{upIndex}[i]},\qquad
\mathrm{upIndex} = [2, 2, 1],$$

$$\mathbf{up}_i = \frac{\mathbf{up}_i^{\text{raw}} - (\mathbf{up}_i^{\text{raw}}\cdot\hat{\mathbf d}_i)\,\hat{\mathbf d}_i}
{\big|\;\cdot\;\big|},\qquad
\mathbf{p}_i = 5.6\,\hat{\mathbf d}_i,\qquad \text{lookAt}(0,0,0).$$

In words: **looking down axis 1 or axis 2, axis 3 points up the screen; looking down axis 3, axis 2
points up** — the same screen-up convention as the **PCA Ellipsoid** page's snap views
(`VIEW_UP_AXIS` and `CELL_VIEW_UP_AXIS`, both `[2, 2, 1]`, in
[`PcaKdePage.jsx:493, 507`](../../web_app/frontend/src/components/PcaKdePage.jsx); those cameras are
described in *Principal axes in the crystallographic frame* → Step 9d, which also documents the
same Gram–Schmidt orthogonalisation of the up vector for oblique cells). *(Not the Atomic
Density page: that nav label maps to page id `structure` → `StructurePage.jsx`, which has no snap
views and no camera-up convention of this kind. The source comment in `renderMinis` says only "the
density page's snap views", which is ambiguous and points at the wrong tab.)* The single
Gram–Schmidt step is what makes this well defined for an oblique cell: "c up" while looking down a
really means *the component of c perpendicular to a*. If that component degenerates
($|\mathbf{up}|^2 < 10^{-10}$, i.e. two cell edges parallel) it falls back to world $(0,1,0)$.

**Distance 5.6 is not arbitrary.** With a 45° vertical FOV the visible half-height at that distance
is $5.6\tan 22.5° = 2.32$, just clearing the relief clamp of 2.2, so a fully-inflated relief surface
can never clip a pane edge.

**Scissoring.** Every invocation first calls `renderer.setSize(width, height)` — `renderMinis` is not
only a camera/viewport pass, it re-establishes the drawing buffer, which is why the mini
`ResizeObserver` can simply call it with no other work. Then `paneHeight = floor(height/3)`; WebGL's
viewport origin is bottom-left, so pane $i$ (counting from the top of the panel) is placed at

$$y_i = \text{height} - (i+1)\cdot \text{paneHeight},$$

with `setViewport` and `setScissor` identical and `setScissorTest(true)` around the loop.
`camera.aspect` is set to `width / paneHeight` per pane.

> **The `floor` leaves a remainder.** The three panes occupy
> $y \in [\text{height} - 3\,\text{paneHeight},\ \text{height}]$, so up to **2 px at the bottom of the
> canvas** is covered by no viewport and no scissor rect: it keeps whatever was last cleared there.
> And because the CSS click zones are sized as a fraction of the *full* mount
> (`.orient-multiview-zone { top: i·100/3 %; height: 33.334%; }`,
> [`PcaKdePage.css:1063-1070`](../../web_app/frontend/src/components/PcaKdePage.css)) while pane $i$ is
> rendered at $[i\cdot\text{paneHeight},\ (i{+}1)\cdot\text{paneHeight}]$ from the top, the buttons and
> the rendered panes drift out of registration by up to 2 px at the third pane. Both are cosmetic at
> the panel's size, but they are real.

**Cost.** The cameras never move, so `renderMinis` is *not* on the animation loop. It runs only from
(a) the end of the mesh-rebuild effect and (b) a `ResizeObserver` on the mini mount. It reads only
refs (`sceneRef`, `miniRef`, `miniMountRef`, `miniAxesRef`), never state, so it is a single stable
`useCallback` with an empty dependency list.

> **Stale-pane edge case.** The rebuild effect disposes the three groups and then early-returns when
> `!result?.polygons`, *before* reaching `renderMinis()`. So on an error or a cleared result the main
> viewport goes empty (its rAF loop keeps rendering the now-empty scene) while the three mini panes
> keep showing the last successfully rendered image until the next successful rebuild.

**Click to snap.** Three transparent `<button class="orient-multiview-zone">` elements are absolutely
positioned over the canvas (`top: i·100/3 %`, `height: 33.334%`, `z-index: 1` above the canvas's
`z-index: 0`), each carrying the axis label in that axis's colour. Clicking pane $i$ calls
`snapMain(i)`, which applies the **same** $\hat{\mathbf d}$/up rule to the main camera:

```js
controls.target.set(0, 0, 0);
camera.up.copy(up.normalize());
camera.position.copy(dir).multiplyScalar(3.05);
controls.update();
```

so the main viewport reproduces exactly the mini's viewing direction and screen-up, at distance 3.05
instead of 5.6. Hover paints a ~38 % dark scrim over the pane
(`color-mix(in srgb, #0b0f16 38%, transparent)`); the selector is written as
`.orient-multiview .orient-multiview-zone:hover` specifically to out-specify the global
`button:hover` rule, and the scrim is that strong because a weaker dark overlay on the light panel
just greys toward white.

---

### Step 12 — The site picker: `SiteStructurePanel`

[`SiteStructurePanel.jsx`](../../web_app/frontend/src/components/SiteStructurePanel.jsx) is **shared
verbatim** with the PCA Ellipsoid page — it was extracted from `PcaKdePage` when this page was
created, behaviour-preserving — and **its scene is specified in *PCA Ellipsoid page* → Step 14**:
the inline unguarded unit-cell division (Step 2 above for what that costs on this page), the cell
wireframe from the eight $\{0,1\}^3$ corners at Hamming distance 1, the `siteFractional` marker
placement, the marker transform $\mathsf{P}^{\!\top}\operatorname{diag}(\mathrm{semi})$ on a shared
`SphereGeometry(1, 20, 16)`, the element-colour/selected-marker states with their two glow shells,
the PC triad, the bond heuristic and the a/b/c gizmo. Every constant it uses is in this section's
Constants table. What follows is only what is *not* in that specification, plus the wiring specific
to this page: the panel is a **sibling of `OrientationView`**, rendered directly by
`OrientationPage` (Step 1), with its own header actions and its own `Loading…` chip.

The symbols used below: $\sigma_i =$ `site.rms[i]` $= \sqrt{\lambda_i}$ (Å) is the RMS displacement
along principal axis $i$; $r_0 = 0.05\times$ the shortest unit-cell edge; the global magnification
is $k = r_0/\overline{\sigma}$ with $\overline{\sigma}$ the arithmetic mean of *all* $\sigma$ over
*all* sites and all three axes; and each semi-axis is $s_i = \max(\sigma_i k,\ 0.15\,r_0)$.

**$k$ is not always a magnification.** The code has two fallbacks
([`SiteStructurePanel.jsx:222-226`](../../web_app/frontend/src/components/SiteStructurePanel.jsx)):

```js
const rmsValues = positions.flatMap((p) => p.rms || []);
const meanRms = rmsValues.length ? rmsValues.reduce(…) / rmsValues.length : baseRadius;
const ellipsoidScale = baseRadius / (meanRms || baseRadius);
```

If **no** site in the table carries `rms`, $\overline{\sigma}$ falls back to $r_0$; and a
$\overline{\sigma}$ of exactly 0 is likewise replaced by $r_0$. In either case $k = r_0/r_0 = 1$ and
the semi-axes become the **raw RMS in Ångström** — roughly 0.1 Å against a several-Å cell, i.e.
effectively invisible markers rather than a magnified footprint. Worse, the first branch is a
*whole-table* property but the marker choice is *per site*: sites that do carry `rms` are still drawn
as ellipsoids at $k = 1$ while sites without fall back to a uniform $r_0$ sphere, so the two marker
families are then on incomparable scales.

Two consequences the shared specification does not draw out, and the panel is honest about neither
on screen:

* `useSiteCloud` requests the sites table with `probability = 0.5` and both engines return
  `semiAxes` $= k_{50\%}\sigma$ with $k_{50\%} = \sqrt{\chi^2_{3,\,0.5}} \approx 1.538$ (*PCA
  Ellipsoid page* → Step 5) — **and the panel never reads `semiAxes`.** The 50 % scaling is
  computed on both the Python and JS sides and discarded here; what is drawn is 1σ times the global
  magnification $k$, so only **relative** sizes and **shapes** are meaningful.
* The `max(..., 0.15 r_0)` floor stops a near-degenerate third axis collapsing into an invisible
  sliver, at the cost of over-stating the thickness of a genuinely planar site.

Marker materials, the two-light rig, the glow-shell scales and the triad proportions are unchanged
from the shared specification; their exact values are in the Constants table below.

**Selection highlight — two independent sources of the same axes.** The glow shells and the
selected-site PC triad are sized off `extent = max(s_1, s_2, s_3, r_0)`, so both scale with the
marker however soft or anisotropic the site is.

The triad is gated on `if (selectedEllipsoid?.axes)`
([`SiteStructurePanel.jsx:312`](../../web_app/frontend/src/components/SiteStructurePanel.jsx)) and is
oriented from the `selectedEllipsoid` **prop**, while the ellipsoid *body* under it is oriented from
`positions[].axes`, read out of the `sites` table. Both normally resolve to the same table entry
(`selectedEllipsoid` is a `useMemo` on `[sites, selectedRef]` in `useSiteCloud`), so in practice they
agree — but they are two independent sources, and the visible consequence is that a selected site
whose `selectedEllipsoid` lacks `axes` gets the shells and no triad even though its marker is a
correctly-oriented ellipsoid. Note also `buildAxisTriad` has **no** degenerate-row guard (Step 10), so
a zero axis row here produces a NaN rod rather than a skipped one.

**Bond criterion.** The adaptive cutoff
$d_{\text{cut}} = \min(3.4,\ \max(2.0,\ 1.3\,d_{\text{near}}))$ Å over in-cell site representatives
only, with $d_{\text{near}}$ the shortest site–site distance above 0.4 Å, is *PCA Ellipsoid page*
→ Step 14. Stated once more because it is easy to over-read on this page: it is a **display
heuristic with no chemistry in it** — no covalent radii, no element-pair table, no periodic images
— so it will miss long bonds beyond 3.4 Å, miss every bond that crosses a cell boundary, and can
connect non-bonded second neighbours in a loosely packed cell.

> **The degenerate branch.** `nearest` is initialised to `Infinity` and the minimisation is over the
> pairs with $d_{ij} > 0.4$ Å only. When that set is **empty** — a cell with a single site, or sites
> that all overlap within 0.4 Å — $d_{\text{near}}$ stays `Infinity` and
> $d_{\text{cut}} = \min(3.4, \max(2.0, \infty)) = 3.4$ Å. The "adaptive" cutoff silently pins to its
> ceiling rather than being undefined. No bonds are drawn in that case either way, but the cutoff is
> no longer material-adaptive.

**The a/b/c gizmo.** Toggled by the `a b c` header button (**default on**), drawn by
`buildCrystalAxes` at the cell origin with length $0.5\,$`edgeLength` and radius $0.035\times$ that.
Two properties matter for reading this page next to the sphere: the edge vectors are individually
normalised, so the gizmo shows the cell **orientation** only (never the relative lengths of
$a$, $b$, $c$), and unlike the sphere's rods these are **one-sided** — they run from the origin
corner outward, not through it (Step 10). This panel is also the one place where a/b/c and the PC
triad are visible *together*; the sphere's frame toggle is exclusive.

**Selection propagation.** The pick itself — a `Raycaster` over `sitesGroup` on `pointerup`, gated
by a 5 px pointer-travel budget so an orbit drag is not read as a click, with the cursor swapped by
the same pick on `pointermove` — is *PCA Ellipsoid page* → Step 14. Three details specific to the
wiring here: the test is `if (moved > 5) return;`, so exactly 5 px still counts as a click; the
callback is reached as `selectRef.current?.(hit.userData.referenceNumber)`, held in a ref and
refreshed by its own effect so the once-mounted scene never has to rebuild to see a new handler; and
that callback is `setSelectedRef` from `useSiteCloud`, so a click here re-selects the site for the
sphere, the panel header title, the axis views and the export filename in one step.

**Framing.** The once-only framing (radius $1.7\times$ the longest cell edge, orientation preserved
across rebuilds, `Reset view` restoring it) is likewise specified there. The numbers behind it:
`defaultRadius = span * 1.7 || 1` with `span = max_i |unit_i|` — the `|| 1` catching a zero or NaN
span — camera at `center + (R, 0.7R, R)`,
`near = R/100`, `far = 100R`. The camera is *constructed* as
`PerspectiveCamera(40, w/h, 0.01, 1000)`, so FOV 40° throughout but near/far only become $R/100$ and
$100R$ at that first framing: a scene that is never framed (no `sites`) keeps the constructor's
0.01/1000. `resetStructureView` re-applies the same framing with its own `radius = defaultRadius || 1`
guard. Later rebuilds keep the user's orientation and only re-set `controls.target` to the cell
centre. A second effect re-fits the canvas
to its mount whenever `sites`/`selectedRef`/`selectedEllipsoid` change and the renderer size differs
from the mount by ≥ 1 px — because the three panels share a row whose content lands in stages.

**Site labels** are built by the host page, not the panel, and use the *same* format as the PCA
Ellipsoid page's dropdown — `#<ref> <element> — U=<uIso, 4 dp> Å²` plus
` (<count>/<copiesPerCell>)` on the browser path — specified in *PCA Ellipsoid page* → Step 4b.
One difference: this page's `<select>` carries **no** clean/merged/fragmented tag under it, so the
`(count/copiesPerCell)` suffix is the only reconstruction diagnostic on screen here (Caveat 5).

---

### Step 13 — Reset and Save (PNG export)

**Reset view** (sphere panel) restores the default camera exactly:
`controls.target = (0,0,0)`, `camera.up = (0,1,0)`, `camera.position = (2.5, 1.85, 2.5)`
(|r| ≈ 3.99). The picker's **Reset view** restores its stored `center`/`defaultRadius` framing.

**Save** is a `SaveMenu` with two options, identical on both panels. Both buttons are `disabled` until
there is something to export — `!result` on the sphere, `!sites` on the picker — as is each panel's
**Reset view**.

| id | Label | Hint | Path |
| --- | --- | --- | --- |
| `png` | Standard PNG | 1× | `renderer.render(...)` then `saveCanvasAsPng(renderer.domElement, name)` |
| `png3x` | High quality PNG | 3× | `setPixelRatio(3)` → `setSize(w, h, false)` → render → `domElement.toBlob` → restore ratio → `setSize` → render again |

Both renderers are created with `preserveDrawingBuffer: true` precisely so `toBlob` can read the
canvas at any time; both paths force a synchronous `render` immediately before the capture rather
than trusting whatever the rAF loop last left there. In the 3× path `renderer.getSize()` returns the
**logical (CSS) size**, so `setSize(size.x, size.y, false)` with `pixelRatio = 3` yields a drawing
buffer of exactly $3w \times 3h$; the `false` third argument leaves the CSS size alone so the page
does not visibly jump. The restore-and-re-render at the end puts the on-screen canvas back.

Filenames go through `sanitizeFilename`
([`figureExport.js:77-83`](../../web_app/frontend/src/figureExport.js)), which is a four-step pipeline,
not a single substitution:

1. an empty/absent name becomes the literal `'figure'`;
2. `^{…}` superscript markup is unwrapped (`.replace(/\^\{([^}]+)\}/g, '$1')`) — inherited from the
   chart exporters, inert for this page's names;
3. every run of non-`[\w.-]` characters collapses to a single `_`;
4. leading and trailing `_` are stripped, and if nothing survives the result is again `'figure'`.

Applied to:

* sphere: `Orientation_<element>_site<referenceNumber>.png`, or `Orientation_sphere.png` when no site
  is selected;
* picker: `PCA_Site_ellipsoids.png` — note this name is **not** parameterised by run or site.

> **The two paths disagree on a failed encode.** `toBlob` can hand back `null`. The 3× path guards it
> (`if (blob) downloadBlob(…)`) and therefore **silently does nothing** — no file, no message. The 1×
> path goes through `canvasToPngBlob`, which **rejects** with `new Error('Could not encode the
> figure')`; neither `saveView` nor `saveStructureView` catches, so that rejection escapes to the
> `SaveMenu` caller. Same failure, two different user-visible outcomes.

**How the export differs from the screen** — all of these are consequences of exporting the WebGL
canvas alone:

* **Everything outside the canvas is missing**: the colorbar and its 1× tick, the whole summary strip
  (peak, anisotropy, ± asymmetry, significance), the legend line (including the `· browser` / `·
  server` provenance tag), the panel title, and any hover tooltip. A saved sphere carries **no scale
  and no numbers**; the reader cannot recover $v_{\max}$ from the image.
* **The background is transparent**, not the panel colour — both renderers use `alpha: true` and the
  scene has no background. On a white slide the dark cell outlines and rods survive; on a dark slide
  the figure needs its own backdrop.
* **"Standard PNG" resolution depends on the display.** It captures the live buffer, which was sized
  with `setPixelRatio(window.devicePixelRatio || 1)` — so it is 2× on a Retina screen and 1× on a
  standard one. Only the 3× option is resolution-deterministic.
* **The Axis views panel has no Save button at all**; the mini renderer is never exported.
* The 3× render happens at the panel's current aspect ratio, so the exported framing is whatever the
  panel's on-screen shape was — a narrow browser window bakes a narrow figure.

---

### Step 14 — Computed but not displayed

The engines return considerably more than this page shows. From the orientation response, the UI
reads only `polygons`, `enhancement`, `vmax`, `cellMeanAmplitude`, `meanAmplitude`, `centers`,
`counts`, `zScore`, `cellCount`, `peakEnhancement`, `peakDirection`, `peakZScore`,
`orientationAnisotropy`, `antipodalAsymmetry`, `antipodalAsymmetryNull`, `significance`, `weight`,
`smoothing`, `browserOrientation`, and `pcaAxes` (fallback rods only). **Never rendered anywhere:**

| Field | What it is | Why its absence matters |
| --- | --- | --- |
| `recommendedFrequency` | the ν the ~12-per-cell guard would choose | With **Auto** selected the resolution control just reads "Auto" and no panel shows the ν or cell count that was actually used. The 2026-07-24 changelog removed the "N cells (ν=…)" summary line as "already in the Resolution control" — true for a manual ν, **false for Auto**. |
| `usedPoints`, `rejectedPoints`, `amplitudeCutoff` | vectors surviving the `Min \|Δr\|` cut, and the Å threshold that quantile resolved to | The header prints `selectedEllipsoid.count` — **all** atoms at the site — so raising `Min \|Δr\|` changes the map without changing the displayed atom count, and the resolved Å cutoff is never shown. |
| `orientationTensor`, `orientationEigenvalues`, `orientationAxes` | $\langle \mathbf u\mathbf u^{\mathsf T}\rangle$, its eigenvalues, its canonicalised axes | Only the scalar $3\lambda_1 - 1$ is shown. The tensor's *axes* go through the same `_eigen_decomposition` as the site's PCA axes and are therefore directly comparable to them (engine section, Step 8) — and are never drawn or tabulated. The Woodcock girdle-vs-rod reading the three eigenvalues support is unavailable in the UI for the same reason. |
| `density`, `mass`, `areas`, `expected`, `sizes`, `antipode`, `neighbors` | per-cell engine internals | `neighbors` is requested (via `geometry: true`) and transferred on every request, then discarded. |
| `vmin`, `meanCount`, `emptyFraction`, `rmsAmplitude`, `pentagonCount`, `totalPoints`, `frequency` | map-quality summaries | `emptyFraction` in particular would immediately expose over-binning (see Caveats). |
| `peakCell` | the **cell index** behind `peakDirection` / `peakEnhancement` / `peakZScore` | The summary strip names the peak by direction but never highlights that cell on the sphere, though the index to do so is right there (and the unit test uses it). |
| `frame` | the engine's echo of the frame it actually computed in | Never read. The legend's a/b/c-vs-PC swatches are chosen from the **`frame` prop** instead (Step 10), which is why they can lead the displayed sphere during a refetch. |
| `source` (Flask path only) | the absolute server-side path of the `.rmc6f` the map came from (`app.py:660`) | Transferred on every HTTP response, never shown — the page cannot tell the user *which file on the server* it plotted. |
| `referenceNumbers`, `totalAtoms` (sites table) | the raw reference-number list and the total atom count | Both returned by `summarizeSites`; the page rebuilds what it needs from `sites.sites`. |
| `semiAxes`, `probability` (sites table) | the 50 %-probability ellipsoid scaling, $k(0.5) = 1.5381722$ (*PCA Ellipsoid page* → Step 5) | Computed on both runtimes; the picker draws 1σ × a global magnification instead (Step 12). |
| `excessKurtosis`, `nonGaussianity` (sites table) | per-axis excess kurtosis and its mean | Both are printed by the PCA Ellipsoid page's statistics panel (*PCA Ellipsoid page* → Step 4b); neither is surfaced here, so a strongly anharmonic site looks the same as a harmonic one in this page's chrome. |

Two client-side computations are also unused: `geometry.computeVertexNormals()` on an unlit material
(Step 4), and the `datasetKey` / `rmc6fText` values returned by `useSiteCloud` that this page does not
destructure.

---

### Parameters and defaults

**Engine options (owned by `OrientationPage`, sent with the request)**

| Control | State | Default | Range / options | Unit |
| --- | --- | --- | --- | --- |
| Site | `selectedRef` | first "clean" site, else first reference number | reference numbers in the file | — |
| Resolution | `frequency` | **10** (1002 cells) | `auto`, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24 → $10\nu^2+2$ cells | — |
| Weight | `weight` | `count` | `count` \| `amplitude` (\|Δr\|) \| `amplitude2` (\|Δr\|²) | — |
| Min \|Δr\| | `minQuantile` | **0** | 0 – 0.5, step 0.05 (displayed as %) | quantile of \|Δr\| |
| Smoothing | `smoothing` | **2** | 0 – 12, step 1 | diffusion passes |
| Frame | `frame` | `cartesian` (**Crystal**) | `cartesian` \| `pca` | — |
| Cluster | `clusterThreshold` | 1.5 | 0.4 – 2.5, step 0.1; visible only when `sites.reconstructed` (browser path) | Å |
| — | `geometry` | always `true` | — | — |
| — | `probability` (sites) | 0.5 (hook default) | — | — |

**Display options (client-side only, no refetch)**

| Control | State | Default | Range | Notes |
| --- | --- | --- | --- | --- |
| Amplitude height | `relief` | **0.5** | 0 – 1, step 0.05 | $\rho$ in Step 5 |
| Colormap | `colormap` | `viridis` | viridis, magma, seismic, reds, greys | 256-entry LUT |
| Contrast | `contrast` | **1.0** | 0.5 – 3.0, step 0.1 | $\gamma$ in Step 6 |
| Cell borders | `showOutline` | on | — | inflate 1.002 |
| Axes | `showAxes` | on | — | selected frame only |

**Constants**

| Constant | Value | Where |
| --- | --- | --- |
| Relief clamps | 0.3 (low) / 2.2 (high) | `reliefFactors` defaults, `orientationSphere.js:44` |
| Vertex key rounding | 9 decimals | `vertexKey`, `orientationSphere.js:36` |
| Negligible amplitude | $10^{-9}$ Å, strict `>` | `NEGLIGIBLE_AMPLITUDE`, `workers/orientation.js:26` |
| Engine frequency bounds | `MIN_FREQUENCY = 1`, `MAX_FREQUENCY = 64` (UI reaches only 2–24) | `workers/orientation.js:24-25` |
| Auto-frequency cap | `maxFrequency = 24`, `targetPerCell = 12` | `recommendedFrequency`, `workers/orientation.js:413` |
| LUT | 5 anchors → 256-entry `Uint8ClampedArray`, index clamped to `[0, 255]`, unknown name → viridis | `colormaps.js:27-57` |
| Outline inflation | 1.002 | `OrientationView` call site |
| Outline material | `0x10151c`, opacity 0.35 | `OrientationView` |
| Colorbar stops | 24 (25 colours) | `colorbarGradient` call |
| Asymmetry flag | `antipodalAsymmetry > 3 × antipodalAsymmetryNull` | `OrientationView.jsx:447` |
| Sphere mount fallback | 640 × 520 px (dead in practice — see Step 1) | `OrientationView.jsx:226-227` |
| Sphere camera | FOV 45°, near 0.01, far 100, default `(2.5, 1.85, 2.5)` | scene effect |
| OrbitControls | damping 0.12, `enablePan: false` | both panels |
| Snap distance | 3.05 | `snapMain` |
| Mini distance / FOV | 5.6 / 45° | `renderMinis` |
| Mini camera construction | `PerspectiveCamera(45, 1, 0.01, 100)` — aspect 1 until the loop | `OrientationView.jsx:301` |
| Mini pane height | `floor(height/3)` rendered vs `33.334 %` clicked | `renderMinis`; `PcaKdePage.css:1068` |
| Mini hover scrim | `color-mix(in srgb, #0b0f16 38%, transparent)` | `PcaKdePage.css:1091` |
| Rod cylinder segments | 10 radial | all three rod builders |
| Rod skip threshold | $\lvert\mathbf d\rvert^2 < 10^{-12}$ (`buildAxisRods`, `buildCrystalAxes`; **none** in `buildAxisTriad`) | `OrientationView.jsx:59`; `sceneAxes.js:53` |
| Gram–Schmidt degeneracy | $\lvert\mathbf{up}\rvert^2 < 10^{-10}$ → world $(0,1,0)$ | `renderMinis`, `snapMain` |
| PC rods | radius 0.006, length 2.7 (each sign) | `OrientationView` |
| a/b/c rods | radius 0.0045, length 2.35 (each sign) | `OrientationView` |
| Tooltip offset | +14 px x, +12 px y from the pointer, canvas-local | `OrientationView.jsx:539` |
| Picker mount fallback | 260 × 260 px | `SiteStructurePanel.jsx:47-48` |
| Picker marker geometry | `SphereGeometry(1, 20, 16)`, shared by every marker and shell | `SiteStructurePanel.jsx:256` |
| Picker marker radius $r_0$ | 0.05 × shortest cell edge | `SiteStructurePanel` |
| Picker semi-axis floor | 0.15 × $r_0$ | `SiteStructurePanel` |
| Unselected marker material | `MeshPhongMaterial` in the element colour (`DEFAULT_ELEMENT_COLOR` if unlisted), opacity 0.45, shininess 20, emissive `0x000000` | `SiteStructurePanel.jsx:262-269` |
| Selected marker material | opaque, shininess 70, emissive 0.35 × `0xff7a1a` | `SiteStructurePanel.jsx:262-269` |
| Picker lighting | `AmbientLight(0xffffff, 0.7)` + `DirectionalLight(0xffffff, 0.7)` at $(1,1,1)$ — the only lit scene on the page | `SiteStructurePanel` |
| Picker length scale / origin | `edgeLength` $=\min_i\lVert\mathbf{unit}_i\rVert$ (Å); `center = toCartesian([0.5, 0.5, 0.5])` | `SiteStructurePanel` |
| Cell wireframe | 12 edges from the eight $\{0,1\}^3$ corners at Hamming distance 1; `LineSegments`, `0x8a97a8`, opacity 0.4 | `SiteStructurePanel.jsx:196` |
| Bond line | `0x9aa3ad`, opacity 0.7 | `SiteStructurePanel.jsx:249` |
| Bond cutoff | `min(3.4, max(2.0, 1.3 × nearest))` Å, pairs with $d > 0.4$ Å; pins to 3.4 Å when no pair exceeds 0.4 Å | `SiteStructurePanel` |
| Gizmo | length 0.5 × shortest edge, radius 0.035 × length | `SiteStructurePanel` |
| Highlight shells | two unlit shells, scale 1.9 @ opacity 0.18 and 1.35 @ 0.28, `depthWrite: false`, both off `extent = max(s_1, s_2, s_3, r_0)` | `SiteStructurePanel` |
| Picker triad | length 2.3 × extent, radius 0.1 × extent | `SiteStructurePanel` |
| Picker framing | radius `1.7 × longest cell edge \|\| 1`, FOV 40°; near/far `R/100` / `100R` after first framing, else the constructed 0.01 / 1000 | `SiteStructurePanel` |
| Click travel budget | ≤ 5 px (`if (moved > 5) return`) | `SiteStructurePanel.jsx:95` |
| Grid | 3 : 6.5 : 6.5, gap 0.85rem, `clamp(20rem, 100vh − 17.75rem, 60rem)`, breakpoint 1200 px | `PcaKdePage.css` |
| Export | 1× = `devicePixelRatio`; 3× = `setPixelRatio(3)` | both panels |

---

### Caveats

**1. The shipped default resolution over-bins by the engine's own criterion — and that silently
disables one readout.** The `Resolution` default is ν = 10 (1002 cells), chosen for a legible picture
together with the 2× smoothing default. The engine's own guard, `recommendedFrequency`, targets ~12
points per cell; for a typical site with ~1000 copies (one per unit cell of a 10×10×10 supercell) it
returns **ν = 3** (92 cells). At ν = 10 the expected count per cell is ≈ 1. Three consequences:

* The colour field is only legible because of the smoothing; the raw map at ~1 count/cell is Poisson
  confetti. Turning smoothing to 0 at the default ν shows that directly.
* Per-cell `zScore` and the `map significance` (RMS of z) are computed from raw counts against
  $\text{expected} \approx 1$, so `significance` sits near 1 (its pure-noise value) unless the
  anisotropy is strong. The readout is **resolution-dependent** and is not comparable between two
  maps at different ν.
* **The ± asymmetry flag can be unreachable.** The null floor is
  $\text{null} = \sqrt{C/(\pi N)}$ with $C$ cells and $N$ surviving vectors, and the UI flags red only
  when $\text{asymmetry} > 3\,\text{null}$. But `antipodalAsymmetry` is bounded above by 1 by
  construction. At ν = 10 and $N = 1000$, $\text{null} = \sqrt{1002/(\pi\cdot1000)} = 0.565$ and the
  threshold is $1.695 > 1$ — **no possible data can trigger the flag**. Reachability requires
  $3\sqrt{C/(\pi N)} < 1$, i.e. $C < \pi N/9 \approx 0.35\,N$; at $N = 1000$ that means $C < 349$, and
  the largest tiling the Resolution dropdown offers under that bound is **ν = 5** ($C = 252$), since
  ν = 6 already gives 362. Users looking for off-centring must lower the resolution (or select a site
  with many more copies) before the flag means anything. The numeric value and its floor are always printed, so
  the information is on screen — only the red flag is dead at the default.

**2. Every surviving vector is binned — there is no point budget, and no RNG.** Unlike the PCA-KDE
path, which explicitly caps its fit with `subsample(points, maxFitPoints, rngSeed)`
([`workers/pcaKde.js:196, 326`](../../web_app/frontend/src/workers/pcaKde.js)), the orientation engine
subsamples **nothing**: every vector above the amplitude cutoff enters the histogram, however many
there are. Two consequences. First, the map is fully **deterministic** given the options — re-running
it cannot shift a lobe, so a feature that moves between runs is a data change, not sampling noise.
Second, a very large site costs time and memory linearly, with no ceiling. Readers arriving from
*PCA Ellipsoid page* → Step 2 ("Subsampling") should not carry that page's cap — or the fact that
its two engines draw *different* 20 000-point subsets — over to this one.

Relatedly, the **reachable resolution range is ν = 2 … 24**: `recommendedFrequency` caps at
`maxFrequency = 24` and `FREQUENCY_OPTIONS` stops at 24, so the engine's `MAX_FREQUENCY = 64` is a
validation bound this page can never reach.

**3. Colour is one-sided, so depletion is compressed.** The scale runs $0 \to v_{\max}$ with the
isotropic level at $1/v_{\max}$. On a strongly peaked map ($v_{\max} = 8$, say) everything from
"completely empty" to "isotropic" is squeezed into the bottom 12.5 % of the colormap. The Contrast
knob helps (it pivots about the isotropic level and pushes depletions down as well as lobes up), but
there is no symmetric/diverging normalisation, and `vmax` is a single extreme cell — one noisy peak
compresses the entire rest of the map.

**4. The relief is a smoothed interpolant, and its two channels are not independent in practice.**
Vertex radii average three cell factors (Step 5), so single-cell features under-bulge. And while
colour (frequency) and radius (mean amplitude) are formally independent quantities, both are computed
from the same smoothed population when smoothing > 0, and with `weight = amplitude` or `amplitude2`
the colour itself becomes amplitude-weighted — at which point the two channels are strongly
correlated by construction. The "independent information" claim in the UI copy holds for
`weight = count` and is weakened for the other two weights.

**5. On a reconstructed file, "site" is a clustering result, not a crystallographic fact.** When the
`.rmc6f` carries no reference-site/cell columns, the browser engine rebuilds the sites by folding
every atom into one unit cell and grouping same-element atoms within the user-chosen `Cluster`
distance (0.4–2.5 Å, default 1.5 Å; Step 2, and *PCA Ellipsoid page* → Step 1b for the clustering
itself). At a threshold that is too generous, two genuinely
distinct crystallographic positions merge into a single "site" — and then the whole page is
describing a **superposition of two clouds about a spurious mean**. The direction map becomes bimodal
for a reason that has nothing to do with anharmonicity, and the ± asymmetry readout in particular is
uninterpretable, because the "site centre" it is measuring displacements from does not exist.

The only on-screen signal is the `(count/copiesPerCell)` suffix in the site label — the same cue the
hook's default-selection heuristic uses — so before reading a map on a reconstructed file, check that
$\text{count}/\text{copiesPerCell} \approx 1$, and re-check the map at two or three cluster distances
to confirm it is stable. Note this diagnostic **does not exist on the Flask path**, which returns
neither `reconstructed` nor `copiesPerCell`.

**6. Nothing on this page is a publication number.** The map is a visualisation of a binned, usually
smoothed histogram. The Python `orientation.py` is the source of truth; the browser port is asserted
to agree (same construction order, half-to-even rounding of the auto frequency, brute-force-verified
assignment), but the browser path is a *visualisation* path in the same sense as the static-mode KDE.
If a number is going into a paper, re-derive it from the package.

**7. Frames.** In the Crystal frame the sphere lives in the shared Cartesian Å basis — **not** in
fractional coordinates — and the a/b/c rods carry orientation only, normalised to equal on-screen
length. A hover direction is therefore a Cartesian triple, never a $[u\,v\,w]$; converting one to
the other requires the cell matrix and is not done anywhere in the app (*Principal axes in the
crystallographic frame* → Steps 7b and 10). For an oblique cell the three Axis views are
consequently not mutually orthogonal projections,
and the "c up" of the "down a" view is only the component of c perpendicular to a. In the PCA frame
the sphere's own $x,y,z$ *are* PC1/2/3 by construction (the engine rotated the directions), so the
crystal rods are correctly not drawn — but that also means a map in the PCA frame carries no
crystallographic reference at all, and for a near-isotropic site the PCA frame itself is arbitrary
within its degenerate subspace (engine section, caveat 9).

**8. The mini panes can go stale**, and in the no-cell-metadata case their labels (`x, y, z`) disagree
with the rods actually drawn (the site's PC axes — unless the site has fewer than 4 vectors, in which
case the labels are right and the rod *colours* are wrong). Both are described in Steps 10–11.

**9. Exports lose the scale.** A saved PNG contains the canvas only: no colorbar, no $v_{\max}$, no
peak/anisotropy/asymmetry/significance numbers, and a transparent background. Record the settings
alongside the image, or screenshot the panel instead.
