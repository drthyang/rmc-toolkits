# Agent Guide

Onboarding for AI agents and new contributors. Human users want [README.md](README.md) /
[QuickStart.md](QuickStart.md). This file is the "pick up where we left off" record: architecture,
key files, conventions, and current state. Full chronological history lives in
[docs/CHANGELOG.md](docs/CHANGELOG.md); forward plans in [docs/ROADMAP.md](docs/ROADMAP.md).
Per-page math and data flow — every formula, default, threshold, and approximation, anchored to the
function that implements it — lives in [docs/ALGORITHMS.md](docs/ALGORITHMS.md); read the relevant
page section before changing an engine, and update it when the math changes.

## What this project is

Post-processing for **RMCProfile** modeling outputs, in three layers. RMCProfile performs atomistic
configuration optimization under experimental constraints; avoid calling those runs
"refinements" (Rietveld-style parameter refinement is a different workflow). Since 2026-07-17 the
app also covers the *pre*-processing step: **Auto StoG** (`rmc_toolkits/scaling.py` + the
`rmc-autoscale` CLI + the Auto StoG tab) automatically scales measured total-scattering S(Q) and
writes the classic stog/RMCProfile-ready file family — see
[docs/STOG_SCALING_PLAN.md](docs/STOG_SCALING_PLAN.md) for the verified math and validation record.

1. **`rmc_toolkits/`** — pure-Python package (parsing, plots, KDE). The source of truth; new app
   code should call into this, not the legacy scripts.
2. **`web_app/`** — Flask API (`backend/app.py`) + React/Vite SPA (`frontend/`).
3. **`src/`** — original standalone research scripts, kept for CLI workflows.

The same React app ships in two runtime modes:
- **Flask mode** — backend serves the built SPA and provides server-side file browsing, SciPy KDE,
  conversion, and Live Data.
- **Static mode** (`VITE_STATIC_MODE=true`) — GitHub Pages build. No backend; the browser parses
  files locally and computes KDE in a Web Worker (WebGPU + CPU fallback).

## Architecture map

```
rmc_toolkits/
  parsers.py     RMC CSV/log, STOG parsing (incl. stog.inp / STOG xy / .dat headers), .rmc6f metadata + atom iteration, Frac*.txt conversion, structure loading (read_structure pairs Frac/rmc6f by stem), shared .rmc6f atom-line grammar (classify_rmc6f_atom_line / parse_rmc6f_atoms / Rmc6fParseReport), robust .log reader (read_chi_log), Rwp column roles (rwp_columns / fit_rwp), the one run-folder rule (find_run_configuration)
  plots.py       plot-kind detection, matplotlib figures, Rwp/chi metrics, PNG serialization
  kde.py         unit-cell position loading + server-side KDE slice (_FixedCovarianceKDE: gaussian_kde with the slab-atom kernel H = bw²·C), decline messages / warnings (+ contours) (source of truth)
  pca_kde.py     per-site RMC displacement clouds → PCA/thermal-ellipsoid stats + separable 3D KDE volume (source of truth)
  orientation.py displacement-direction distribution: Goldberg (hex + 12 pentagon) sphere tiling, solid-angle histogram, enhancement/z-score, orientation tensor; calibrated peak (Poisson + Šidák), map (exact-moment X²), antipodal (conditional binomial) and Bingham tests (source of truth)
  transforms.py  Keen-2001 conversions, sine-FT pair, Lorch, omitted-low-Q correction (cancellation-free low-Q basis), Fourier filter (r = 0-safe, gpdf_slope_at_zero), stog low-r enforcement (pure numpy, no I/O)
  scaling.py     Auto StoG engine: level sweep, closed-form (a,b) fit (Huber IRLS, rows scaled by √w), self-consistent filter loop, FZ amplitude mode + its conditioning (fz_limit_fit: a_fz standard error, a_fz_reliable), first-shell detection (smallest-r shell, either sign) + data-located low-r window (_place_low_r_window, first_shell_candidates), automatic enforcement cutoff at the first-shell foot (first_shell_foot / auto_enforcement_cutoff), aliasing limit (alias_limit = π/max ΔQ), estimate_rho0 (density from amplitude-criteria concordance), diagnostics (no I/O)
  scattering.py  Faber-Ziman coefficients from a chemical formula (Sears table, 89 elements)
  scaling_cli.py rmc-autoscale CLI: stog.inp/--data → engine → classic stog output family + provenance JSON (all scaling file I/O); resolve_coefficients (⟨b⟩²/⟨b²⟩ from one source), usable_sigma, refuse_failed_fit, stog_inp_closest_approach = the stog.inp line-22 r0 rule — all shared with the API
  triplets.py    bond-angle (triplet) distribution engine, RMCProfile triplets-style: A-B-C with B central, per-bond length windows, linked-cell periodic search with explicit image shifts (triclinic-safe), exact per-bin sin(θ) correction; bond_angle_summary is the JSON payload contract shared by /api/triplets and workers/triplets.js, bond_angle_summary_from_file the uncached file entry point (source of truth); angles streamed (PAIR_CHUNK) with an exact-count work budget (APP_MAX_ANGLES) checked by a count-only search at both app boundaries; EDGE_SNAP_DEG / WINDOW_TOL make ideal geometries deterministic
  triplets_cli.py rmc-triplets CLI: .rmc6f or run folder (same .rmc6f as the app) → angle histogram CSV (+ optional PNG plot / raw-angle dump)

web_app/backend/app.py    Flask API; data-root guard; `_number()` numeric-parameter validation (bad input → 400); StrictJSONProvider (NaN → null in data series) + `_strict_result_response()` / `_require_finite_scaling()` (a non-finite computed result → 400); file-signature LRU caches (`_FileCache` keyed on `_file_signature()` = (st_mtime_ns, st_ctime_ns, st_size, st_ino); a torn read is never cached → 409) for KDE, PCA, triplets and scaling; /api/scaling/preview|run share the CLI writer

web_app/frontend/src/
  App.jsx                        shell, run-folder selection, page nav, Live Data (`configEpoch` → the analysis pages' `dataEpoch` prop)
  browserData.js                 static-mode local file parsing + run assembly (chooseStructureFile mirrors parsers.find_run_configuration, code-point tie-breaks)
  rmc6f.js                       shared .rmc6f atom-line grammar + parse report (classifyAtomLine, parseRmc6fAtoms, readRmc6fCellVectors) — mirror of parsers.py; keep in sync
  plotDomain.js                  one-pass axis domains, null-tolerant hover search (nearestFiniteIndex), plot payload checks
  symmetry.js                    browser-only space-group finder: strain-tested lattice rotations, least-squares-refined operations, closure walk, tolerance ladder, orbits
  spaceGroupSymbol.js / spaceGroupTable.js   H-M symbol, standard-setting search, centring from the full translation lattice, 230-group table
  wyckoff.js / wyckoffTable.js   Wyckoff letters (1731 ITA positions) read in the naming cell
  symmetryModel.js               structure → finder glue, 2000-site cap and 384-operation budget, orbitLabel()
  siteLabel.js                   mixed-occupancy site label (composition) shared by the PCA Ellipsoid and Displacement Directions pages
  colormaps.js                   colormap LUTs for the KDE canvas
  orientationSphere.js           pure helpers for the orientation sphere (cell mesh/outline typed arrays, relief radii, colorbar gradient) — unit-tested, no Three.js
  pcaCrystalFrame.js             pure crystallographic-frame math (3×3 algebra, unit-cell vectors from the supercell lattice, fractional ⟷ PCA transforms, per-PC angles to a/b/c + [u v w], `crystalOrientationRows` (the Crystal orientation table), `projectVolumeOntoFrame` (crystal-frame wall marginals by line integrals through the volume)) — unit-tested, no Three.js
  useSiteCloud.js                shared hook: .rmc6f text loading, worker/API request routing, per-site ellipsoid table, selected site (one app-lifetime worker → shared parse cache across the PCA Ellipsoid and Orientation pages); `dataEpoch` in its request dependencies reloads in place on a Flask Live Data save
  api.js                         frontend API base URL config (VITE_API_BASE_URL)
  llm/                           experimental AI assistant — local LLM (Ollama/LM Studio) or cloud (OpenAI/Gemini); see its README
    context/                     dashboard state → compact LLM context JSON (symmetry + per-site displacements, pair correlations, char budget)
    provider/client.js           OpenAI-compatible client (models, SSE streaming incl. reasoning, connection hints)
    prompts/                     shared system prompt + chat/watchdog message builders
    watchdog/                    convergence heuristics (source of truth) + LLM-narrated badge hook
    useAssistant.js              shared hook: settings, connection probe/auto-connect, run context
    components/                  AssistantPage (chat-only) + connection bar, settings drawer, ChatView (Thinking panel), WatchdogBadge
  components/
    AutoStogPage.jsx             Auto StoG tab — pre-processing, fully client-side in BOTH runtimes and independent of the run folder: page-local S(Q) upload (± stog.inp) → grouped params (fieldsets w/ descriptions) → worker auto-scale (+ rho0 self-consistency estimate when rho0 is empty) → readout + S(Q)/GK/D(r) plots → zip export. Does NOT call /api/scaling/* (those remain for API/CLI use)
    Dashboard.jsx                all-plots run dashboard
    ModelSummary.jsx             Model information + Detected SG cards (parse warning, move counters, tolerance ladder)
    InteractivePlot.jsx          browser-native SVG plot renderer (hover, legend, drag-zoom)
    PlotViewer.jsx               PNG plot rendering + metadata
    StructurePage.jsx            KDE slice, Slab In Cell, Three.js 3D view  ← most complex component
    PcaKdePage.jsx               PCA Ellipsoid tab: site picker, ADP table, Three.js isosurface + ellipsoid + wall projections; unit-cell picker via SiteStructurePanel
    OrientationPage.jsx          Orientation tab (own workspace page — not a PCA product): owns the options (top controls bar, PCA-page style) + the 3:6.5:6.5 equal-height grid (Axis views : sphere : SiteStructurePanel), height viewport-clamped for 16:9
    BondGeometryPage.jsx        Bond Geometry tab (beside Atomic Density): triplet + window controls, Model information card (ModelSummary with showSymmetry={false} — no Detected SG card; same structure source as Dashboard/StructurePage), result chips, angle-distribution panel (sin-corrected|density toggle), bond-length step histogram, partial-g(r) window helper; 16:9 viewport-clamped grid with flush card edges; compute-on-demand via useSiteCloud requestPca('triplets') → worker or /api/triplets
    OrientationView.jsx          renders display:contents → its two panels drop into the page grid: the Axis-views mini panel (three fixed-angle a/b/c | PC1/2/3 views, click to snap) and the sphere panel (flat-shaded Goldberg cells, amplitude relief, only the selected frame's axis rods, header Crystal|PCA toggle + Reset + Save, per-cell hover, colorbar + asymmetry/significance strip)
    SiteStructurePanel.jsx       clickable unit-cell site picker (thermal-ellipsoid markers, bonds, a/b/c gizmo) shared by the PCA Ellipsoid and Orientation pages
    sceneAxes.js                 shared axis palettes (PC tricolor, a/b/c) + triad/rod builders for every Three.js panel
    FileExplorer.jsx             file navigation
  workers/
    localKdeWorker.js            static-mode KDE worker (same kernel, decline rules, warnings and slab test as kde.py; GPU-or-CPU density map, contours); parity-tested against Python goldens (kdeParity.test.js ← tests/generate_kde_fixture.py)
    gpuKde.js                    WGSL compute-shader density map + shouldUseGpu heuristic + cached device init
    slabSelection.js             shared slab membership (isInSlab, SLAB_FACE_TOLERANCE), Miller-plane labels, kernel σ and slab thickness in Å — pure, used by the KDE worker and StructurePage
    pcaKde.js                    static-mode PCA-KDE engine (JS port of pca_kde.py): 3×3 Jacobi eigensolver, per-site clouds, separable volume + projections
    autoScale.js                 static-mode Auto StoG engine (JS port of scaling.py + transforms.py + stog parsers + Faber-Ziman); parity-tested against Python goldens (autoScale.test.js)
    autoScaleWorker.js           off-thread runner for autoScale.js (transferable buffers)
    orientation.js               static-mode displacement-orientation engine (JS port of orientation.py: Goldberg hex+pentagon sphere tiling, solid-angle histogram); parity-tested against Python goldens in orientationFixes.test.js (incl. JS regularizedGamma/normalQuantile pinned to scipy)
    triplets.js                  static-mode bond-angle engine (port of triplets.py; parity-tested against Python goldens in triplets_fixture.json — regenerate with tests/generate_triplets_fixture.py)
    pcaKdeWorker.js              static-mode PCA-KDE worker (parses clouds once, answers 'sites'/'kde'/'orientation'/'triplets' requests off-thread)
    marchingCubes.js             isosurface extraction over a scalar field (Lorensen-Cline tables) for the Three.js KDE surface; `sampleFieldTrilinear` (NaN outside the grid — never clamp)
```

## Key conventions & gotchas

- **`z` / `dz` are fractions of the unit cube's projection range along the slice normal** at the
  API/slider boundary. They equal cell-edge fractions only for the `a`/`b`/`c` presets; for (111)
  the range is √3. The KDE works entirely in fractional coordinates: nothing in `kde.py` or
  `/api/kde/slice` converts to Ångström, and Å enter only at draw time (`StructurePage.jsx`,
  through `unitCell.unitVectors`). The real slab thickness is `dz·(|h|+|k|+|l|)·d_hkl`, printed on
  the map (`slabThicknessAngstrom()` in `workers/slabSelection.js`). Both payloads echo `z`/`dz` as
  given; the Flask payload's `depth`/`depthThickness` are absolute depth-projection units. Keep
  that contract when touching KDE code (docs/algorithms/notation.md §3b).
- **Structure KDE: one kernel, one slab test, two runtimes.** `kde.py` (source of truth) and
  `localKdeWorker.js`/`gpuKde.js` draw the same kernel `H = bw²·C`, where `C` is the covariance of
  the slab's *source atoms* (one row per atom, periodic images excluded; `_source_atom_rows` /
  `makeSlab`). They decline the same slabs with the same `message` (`KDE_MESSAGES`,
  `COVARIANCE_CONDITION_LIMIT = 1e-10`) and attach the same `warnings`: `subgrid` (minor σ below
  half a grid step), and `unresolved` when the grid-summed linear density is below
  `UNRESOLVED_MASS_LIMIT = 1e-6`. An `unresolved` map is neither contoured nor painted. A kernel with
  σ below `KERNEL_MIN_SIGMA = 1e-10` (in-plane fractional units, e.g. `bw = 1e-200`) is not
  evaluated at all: both engines return the finite zero map with its kernel, flagged `subgrid` +
  `unresolved` (HTTP 200). Slab membership is `|d − z_c| ≤ dz/2 + SLAB_FACE_TOLERANCE (1e-9)` in
  `kde.py` and `workers/slabSelection.js`; the page's Slab-In-Cell highlight uses the same helper.
  `kde.py` evaluates through `_FixedCovarianceKDE`, a `scipy.stats.gaussian_kde` subclass that sets
  the covariance attributes of every SciPy release (`cho_cov` for ≥ 1.10, `inv_cov`/`_norm_factor`
  for older ones) and checks itself at construction. On a SciPy it cannot drive, it declines with the
  Python-only `engine` message rather than failing the request. `kdeParity.test.js` checks the
  worker against Python goldens: re-run `PYTHONPATH=. python tests/generate_kde_fixture.py` whenever
  `kde.py` changes. The kernel's shape follows the slab's site layout; that is a known, documented
  artefact (docs/algorithms/structure.md), not a bug to regularise away.
- **Rwp is `None`/`null` when it is undefined, never `0`.** `rwp()` exists twice — `parsers.py`
  (source of truth) and `browserData.js` (static-mode port) — and the two must agree exactly; keep
  them in sync. It sums only the points where *both* the observed and fitted columns are finite
  (a masked region reaches the readers as NaN). The observed column is the EXPERIMENT: RMCProfile
  writes (x, calculated, experimental), and `fit_rwp()`/`fitRwp()` pick the roles with
  `rwp_columns()`/`rwpColumns()` (a header naming exp/obs vs calc/rmc/fit wins, positional order
  otherwise). It returns the unavailable sentinel for the two degenerate cases: no such point, and a
  zero denominator. Neither is a fit quality, and `0` in that slot reads as a *perfect* one — the
  dashboard chip renders the sentinel as "—".
- **KDE fit is subsampled to 6000 slab points** (deterministic pseudo-random, to avoid RMC
  atom-order aliasing); the kernel is fitted to all the slab's source atoms, not the subsample.
  Static-mode KDE is a *visualization* path — the server-side SciPy path is the reference for
  publication values (below the cap the browser CPU branch matches it to 1e-6 of the peak).
- **GPU KDE must always degrade gracefully.** Missing `navigator.gpu`, no adapter, device/shader
  error, lost device, a non-finite read-back, or sub-threshold work all fall back to the CPU loop,
  which evaluates the same kernel in float64. GPU is used only when `grid*grid*samples >= 2_000_000`.
- **PCA-KDE is separable, not approximate.** `pca_kde.py` samples the 3D KDE on a grid aligned with
  the cloud's principal axes; with SciPy's `H = factor²·C` bandwidth, `C` and `H` are both diagonal
  in that frame, so the Gaussian factorizes into three 1D kernels and the volume is their tensor
  product (`N·3·grid` exponentials, contracted via BLAS, instead of `N·grid³`). The result equals
  `scipy.stats.gaussian_kde` to round-off — `tests/test_pca_kde.py` and the JS
  `pcaKde.test.js` both assert exact agreement against the full estimator. `pcaKde.js` is a straight
  port (a 3×3 Jacobi eigensolver with a relative stopping test stands in for `eigh`); keep the two in
  sync when touching either.
- **PCA displacement convention**: an atom's offset is `coords − cellIndices/supercell`, unwrapped
  over the *supercell* period about its own site's circular mean per axis (`o −= round(o − centre)`,
  `_circular_site_centres` / `circularMean`). Never fold about zero: that tears a site at x = ½ in a
  one-cell-thick box, or near x = 1 in a two-cell box, into halves a box edge apart. The offset is
  then mean-subtracted per reference site and mapped to Cartesian Å through the full supercell
  `latticeVectors`. Clouds pooled by element are meaningful because each site is already centred on
  its own average position. Pooling selects atoms by their OWN element (`atom_elements` /
  `atomElements`, also in the browser orientation histogram via `displacementCloud`), and a
  mixed-occupancy site is labelled by its majority species (ties to the alphabetically first) with
  `elementCounts` + `mixed`; the pages name it by its composition (`siteLabel.js`).
- **PCA statistics conventions (1.0)**: `nonGaussianity` is Mardia's multivariate excess kurtosis,
  (b₂ − 15)/5. It is rotation- and affine-invariant, and it equals the marginal κ of any elliptical
  distribution. A symmetric split site is NEGATIVE. Per-axis κ means something only where
  `axisResolved` holds (eigenvalue gap > 3 SE, `AXIS_RESOLUTION_SIGMAS`). A site with λ₁ <
  `ZERO_SPREAD_VARIANCE` (1e-8 Å², e.g. an `*AVERAGE.rmc6f`) is `zeroSpread`: its axes, anisotropy
  and κ are null and the KDE refuses it. `cubicBox` only sizes the display box (`boxHalfWidths`); the
  volume, mass, iso levels and PC walls are always sampled on the per-axis box. `_shape_statistics`
  (pca_kde.py) and `shapeStatistics` (pcaKde.js) are line-for-line ports; keep them in sync. The
  browser χ²₃ quantile is exact (`chiSquare3Quantile`).
- **Orientation histogram bins are hexes + exactly 12 pentagons, area-normalized.** A sphere cannot
  be tiled by hexagons alone (Euler), so `orientation.py` uses a Goldberg polyhedron (10ν²+2 cells)
  and divides every count by the cell's exact solid angle — never plot raw counts, or the 12
  pentagons print the icosahedron onto the map. `enhancement = 4π·density` is 1 for an isotropic
  site. `zScore` (and `peakZScore`) come from raw counts only (never after smoothing) and are
  *local, uncorrected* values, never significances. Every significance readout is a calibrated
  one-sided normal deviate, with tails floored at 1e-300: `peakSignificance` (exact Poisson tail of
  the peak cell's raw count, Šidák-corrected over all C cells), `mapSignificance` (Pearson's X²
  against a gamma matched to its exact multinomial mean, variance and skewness — not χ²_{C−1}, which
  is badly anti-conservative below ~0.1 atom per cell; `null` below `MAP_TEST_MIN_PAIRS` = 0.1
  expected coincident pairs), `antipodalAsymmetryZ` and `orientationAnisotropySignificance`
  (Bingham χ²₅); the legacy `significance` is an RMS z, not σ. The map is deliberately **never
  antipodally folded**: a +u/−u imbalance about the site mean (skewness from odd-order anharmonicity
  or unequally occupied opposite wells) is the signal the ellipsoid cannot show; a coherent
  off-centring shared by every copy moves the site mean and is invisible here. `antipodalAsymmetry`
  is tested against its exact inversion-symmetric null (per antipodal pair X ~ Bin(T, ½): `…Null` ±
  `…NullSd`, flag `antipodalAsymmetrySignificant` = z > 3). Exact Voronoi ties go through a
  centrosymmetric hemisphere fold and then to the lowest cell index (`ASSIGN_TIE_TOL`); tied peak
  cells go to the lowest index within `PEAK_TIE_RTOL = 1e-9`; neighbour/polygon cycles start at their
  smallest index. `recommended_frequency` floors to ≥ 12 points per cell. `workers/orientation.js` is
  a straight port: keep the two in sync, with the same construction order so cell indices agree, and
  regenerate the shared golden values in `tests/test_orientation_fixes.py` ↔
  `orientationFixes.test.js` together.
- **Principal axes in the crystal frame: angles are Cartesian, `[u v w]` is fractional.** The PCA
  axes and the unit-cell vectors live in the same Cartesian basis, so ∠a/∠b/∠c (`pcaCrystalFrame.js`)
  are honest angles between directions for any cell, oblique included. `[u v w] = M⁻ᵀ·axis` is a
  *direct-lattice* direction — in an oblique cell it is NOT normal to the like-indexed (h k l) plane,
  so never present it as one. An eigenvector's sign is arbitrary, so a direction and its negative are
  the same axis: `crystalOrientationRows` picks the sense that makes the closest crystal axis acute,
  and the PCA Ellipsoid table shows that representative.
- **`StructurePage.jsx` canvases render conditionally** (`{structure && (...)}`). Effects that attach
  listeners to those canvases must depend on `structure` (or the canvas ref), not `[]` — otherwise
  they run at mount before the canvas exists and never attach. (This was the slab-drag bug, fixed
  2026-06-21.)
- **Slab In Cell drag**: cursor→slice mapping inverts the 2D plane projection in `makePlaneMapper`
  (`invert()`); the band geometry is published each render into `slabGeometryRef` for the pointer
  handlers. Drag updates `zCenter` live.
- **Three.js atom palette** is a Nature-style scheme; each element gets a distinct color shared by
  the slab and 3D views via a legend above them. Every Three.js view calls
  `renderer.forceContextLoss()` after `dispose()` on teardown, so rebuilds and remounts do not pile
  up WebGL contexts.
- **rho0 self-consistency is criteria concordance, not a new fit.** The density-limit amplitude
  depends on rho0 (C2 rows scale with `-4·pi·rho0·r`) — degenerate with the scale from low-r alone;
  the Q→0 Faber-Ziman amplitude does not (needs only ⟨b²⟩ + the measured level). `estimate_rho0`
  (scaling.py, JS port in `autoScale.js` — keep in sync) root-finds
  `a_fz/a_density(rho0) = 1` by fixed-point iteration. It REQUIRES a composition, confines rho0 to
  `RHO0_PHYSICAL_RANGE` = [0.005, 0.25] Å⁻³ and accepts a concordant root only where
  `density_limit_satisfied` holds (else `converged=False` with a `reason` — and `stopped` when a trial
  density could not be fitted — which the CLI and page quote when they refuse it). Both engines run
  the same deterministic iteration (since 1.0 the JS auto loop filters with the same S(0) target as
  Python), so the iterated rho0 agrees to round-off: `autoScale.test.js` asserts 1e-10 relative with
  `converged`/`iterations` equal (measured ≤ 1e-13 on the parity fixture).
- **Auto StoG first shell and low-r enforcement.** `detect_first_peak_onset` returns the FIRST
  coordination shell (smallest r, either sign: inverted negative-b shells such as Ti–O or Mn–Sn
  count), not the strongest |g| feature; a flank that reaches the search start is not a shell.
  Without `r0`/`r_fit_max`, `autoscale` locates the density-limit window from two trial fits
  (`START_WINDOW_WIDTHS`) and refits on [lo, onset − 0.25]; a candidate onset is accepted only when
  its own refit re-detects it (±0.15 Å, at most `MAX_WINDOW_REFITS` = 4 refits), and `r0_detected`
  is the onset the window was built from, so the window top is always `r0_detected − 0.25`. It
  raises instead of fitting across a first shell that starts within ~0.55 Å of `r_cutoff` (bonds
  < ~1.75 Å at the default 1.0 need a lower `r_cutoff`), when no shell can be confirmed, and when a
  refit gives a ≤ 0 — a fit with a non-physical scale or a window across the first shell is never
  returned. The automatic enforcement cutoff is `auto_enforcement_cutoff` = min(foot, anchor − 0.25)
  with anchor = min(detected onset, given r0) — never the onset itself (that zeroed 6–9 % of the
  first shell) and never above a given r0. `scaling._place_low_r_window` ↔ `autoScale.js
  placeLowRWindow` and `first_shell_candidates` ↔ `firstShellCandidates` are straight ports: keep
  `autoScale.js` in exact parity (fixture `expected.detector|window|enforcement`) and regenerate with
  `PYTHONPATH=. python tests/generate_autoscale_fixture.py` from a worktree (without it the editable
  install imports the main checkout). The real-data sweep regression is opt-in:
  `RMC_TOOLKITS_FULL_SWEEP=1 python -m unittest tests.test_stog_a_placement` (~2–5 min).
- **Huber IRLS scales rows by √w** in `_solve_affine` and `fz_limit_fit` (and their JS twins), so the
  fixed point is Huber's M-estimator (c = 1.345); `a_fz_rel_se` is its sandwich standard error.
  Before 1.0 rows were scaled by w (a redescending estimator). `a_fz_reliable = True` is necessary,
  not sufficient.
- **Classic-named Auto StoG outputs follow the Fortran conventions.** `scale.gr` / `<stem>.gr` hold
  the unfiltered g(r) (≈1 at large r), and `scale_ft.gr` / `<stem>_ft.gr` hold g_filtered(r) plus a
  third column r·[g−1]. Only the `_rmc` files are Keen G_K / D / F_K. Before 1.0 they held g−1 and
  4πρ0 r(g−1).
- **⟨b⟩² and ⟨b²⟩ come from one source.** `scaling_cli.resolve_coefficients` (CLI + API) and JS
  `resolveCoefficients` (page) pair a composition's ⟨b²⟩ only with an agreeing ⟨b⟩² (within 2 %);
  the CLI prints its warnings and `/api/scaling/*` return them as `warnings`. `ScalingConfig` /
  `makeConfig` reject ⟨b²⟩ < ⟨b⟩² (S(0) > 0).
- **Q order and failed fits.** `crop_sq` / `cropSq` drop non-finite rows, sort to ascending Q and
  raise on duplicate or overlapping Q (concatenated banks). The public transforms raise on
  non-increasing grids. An auto fit with a ≤ 0 carries `fit_failure`, is never converged, and is
  never written (CLI `refuse_failed_fit`, /api/scaling/run, page worker). Despike runs once, and
  `n_despiked` is the true count.
- **Triplets constants are mirrored.** `APP_MAX_ANGLES` (5×10⁷), `EDGE_SNAP_DEG` and `WINDOW_TOL`
  (both 1e-9) exist in both `rmc_toolkits/triplets.py` and `workers/triplets.js` and are pinned
  equal through `triplets_fixture.json`. Regenerate the fixture with
  `PYTHONPATH=. python tests/generate_triplets_fixture.py` — run as a script it otherwise imports
  whichever `rmc_toolkits` is installed (e.g. an editable install of another checkout). When A and C
  are the same element each physical triplet counts once (one bond in each window, either
  assignment). `lengths.count` / `bond12_count` are B-centred (2× the bonds when the end element is
  the central one); `uniqueBonds` / `unique_bonds12` count each bond once. The budget bounds the
  angles formed, not the n² candidate pairs of an A = C request with distinct windows.
- **One `.rmc6f` atom-line grammar in both runtimes** (parsers.py ⟷ rmc6f.js): `id element [label]`
  then exactly 7 (`x y z ref cx cy cz`) or 3 (coords-only) fields, validated (Fortran `D` exponents,
  cell indices inside the supercell, any Atoms-marker spelling, bare-CR files); non-finite lines are
  skipped and counted; the count is compared with the header's `Number of atoms:`
  (`Rmc6fParseReport` / `report`, surfaced as `parseWarning` by /api/structure, /api/pca/*,
  /api/triplets and the worker). Element tokens are capitalized alike (`SE` → `Se`).
  `iter_rmc6f_atoms` yields full-layout atoms by default; position-only consumers (KDE loader,
  triplets loader, /api/structure) pass `include_coords_only=True`. Zero parseable atoms is an error
  naming what was found (HTTP 400). The committed demo run is exercised by
  `tests/test_parsers_demo_run.py` and `__tests__/demoRun.test.js`; coords-only parity by
  `rmc6f_coords_only_fixture.json` (`tests/generate_coords_only_fixture.py`); Flask ⟷ browser plot
  parity by `__tests__/fixtures/plot_parity_fixture.json` (regenerate with
  `python tests/generate_plot_parity_fixture.py`).
- **One run-folder rule.** Which `.rmc6f` a folder means is `parsers.find_run_configuration()` —
  used by the backend's `_find_rmc6f()` and the `rmc-triplets` CLI — and `chooseStructureFile()` in
  the browser: the configuration a run output names (by priority, then lower-cased name), else the
  first by name; 0-byte / marker-less candidates never hide a usable one; ties break in code-point
  order in both runtimes.
- **Backend parameters and caches.** Every numeric query/JSON parameter goes through `_number()` in
  `app.py`: it must be finite, integral where required, and inside its per-parameter range, and a
  violation is a ValueError → HTTP 400 in every route. Every parsed-file cache is a `_FileCache`
  keyed on `_file_signature()`, never `st_mtime` alone. docs/REFERENCE.md lists every route with its
  ranges; update it when you add a route (tests/test_backend_contract.py enforces this).
- **Non-finite numbers: data series → `null`, computed results → 400.** `app.py` installs
  `StrictJSONProvider`: a NaN in a parsed data series (a masked CSV region, a NaN log row) reaches the
  browser as JSON `null` (`allow_nan=False` guards every response), and the frontend treats null as a
  gap (`plotDomain.js` `nearestFiniteIndex` / `plotPayloadError`). A KDE-slice, PCA-KDE or
  orientation *result* that comes out NaN/Infinity for finite but extreme parameters is refused with
  HTTP 400 by `_strict_result_response()`, which runs the stdlib encoder with `allow_nan=False`
  itself, so it does not depend on the installed provider; the scaling routes do the same through
  `_require_finite_scaling()`. Never send a computed result through the nulling path.
- **Flask-mode Live Data reloads the analysis pages in place.** `App.jsx` watches the `.rmc6f`
  signature in `/api/files` on every Load and every Live Data poll, and bumps `configEpoch` when it
  changes. The Atomic Density, Bond Geometry, PCA Ellipsoid and Displacement Directions pages take
  it as a `dataEpoch` prop in the dependencies of their backend fetches (`useSiteCloud`'s
  `requestPca` and site table, the structure and partials requests), so a page never mixes two
  configurations and keeps the picks that still apply; Bond Geometry drops its computed
  distribution (never recomputed unasked). The pages are not remounted (docs/algorithms/notation.md
  §3c).
- **The Detected SG card never shows a wrong number.** Every reported operation set is a closed
  group, and a symbol or ITA number is accepted only when it belongs to the detected class and
  centring in a standard setting. Otherwise the card shows the crystal class, a "≥" lower bound,
  "undetermined" or "not analysed" (number null, no Wyckoff letters). Wyckoff labels are
  `wyckoffMultiplicity` + letter (naming-cell multiplicity, `orbitLabel()`), not the given-cell orbit
  size — the assistant context uses the same rule. The finder is browser-only (main thread) and has
  no Python source of truth.
- **Backend data-root guard**: relative paths resolve under `RMC_TOOLKITS_DATA_ROOT` (default repo
  root); absolute paths are rejected unless inside the root or a natively-picked folder.
- **Package root exports**: the public engine API (parsers, KDE, PCA, orientation, triplets,
  scaling, transforms) is re-exported from `rmc_toolkits`; `tests/test_package_api.py` checks that
  every `__all__` name resolves. Add new public functions there.
- **`src/llm/` import boundary**: the AI assistant module receives run data **only as props**
  (`runName`, `plotFiles`, `rValueFile`, `structure`, `symmetry`, `liveData`) and must not import
  from the rest of the app except `figureExport.js` (`downloadBlob`/`sanitizeFilename`). Cell math
  is duplicated from `ModelSummary.jsx` on purpose. This keeps the module extractable — don't
  "clean up" the duplication by adding host imports. The R-value series it receives is **ln(χ²) of
  the last .log column** (named in `plotData.chiColumn`, e.g. `X_ray_(R)1` — one fit term, not the
  total; browserData applies `Math.log`); the context builder labels it so the model reads it
  correctly, and reads `non_gaussianity` as Mardia's kurtosis (sites ranked by its magnitude).

## Run & test

```bash
# Backend (venv with numpy scipy flask flask-cors matplotlib)
source .venv/bin/activate
RMC_TOOLKITS_PORT=5050 python web_app/backend/app.py

# Frontend dev
cd web_app/frontend
VITE_API_BASE_URL=http://localhost:5050 npm run dev

# Static-mode dev (no backend)
VITE_STATIC_MODE=true npm run dev

# Tests (tests that need the gitignored data/ runs skip when they are absent)
MPLCONFIGDIR=/tmp/rmc_toolkits_matplotlib python -m unittest discover -s tests

# Lint frontend
cd web_app/frontend && npx eslint src

# Frontend unit tests (vitest — engine ports + Python goldens, components, src/llm)
cd web_app/frontend && npm test
```

CI (`.github/workflows/tests.yml`) runs the Python suite on a matrix — Python 3.9 with the
dependency floors pinned exactly (numpy 1.22.4, SciPy 1.8.1, matplotlib 3.6.3, contourpy 1.0.7),
3.11 and 3.13 with the latest releases — plus the frontend lint, vitest and build, on every push/PR
to `main`. `.github/workflows/pages.yml` deploys the static dashboard. The `rmc_toolkits` package is
pip-installable (`pip install -e .`, see `pyproject.toml`; floors numpy ≥ 1.22, scipy ≥ 1.8,
matplotlib ≥ 3.6, contourpy ≥ 1.0.7) and exposes `__version__` (1.0.0).

The repo's sample data lives in `data/` (GaNb₄Se₈ runs, gitignored). Point the run folder at a
subdirectory containing a `.rmc6f` (e.g. `data/5K_try1`) to exercise the KDE/3D page. The committed
demo run (`web_app/frontend/public/demo`, GaTa₄Se₈ 250 K) backs the parser, KDE-golden and
plot-parity tests, so those run in CI.

> The machine's Anaconda Python has a broken numpy and no Flask — always use the dedicated `.venv`.
> Golden generators (`tests/generate_*.py`) import whichever `rmc_toolkits` is installed: run them as
> `PYTHONPATH=. python tests/generate_….py` from the checkout under test.

## Current known issues

- **iPhone Safari static mode** (2026-06-17): unreliable after selecting a local run folder; desktop
  static mode works. Likely in the mobile folder-selection / file-enumeration path (the atom-line
  parser, once suspected, was replaced by the shared grammar in 1.0 and reports what it skips). Keep
  the local Flask workflow as the supported path for mobile until a unified run-source abstraction
  lands.
- `src/RMC_3D.py` imports Mayavi and runs work at import time. `src/STOG_plot.py` also has
  top-level plotting, but STOG should stay hidden until a dedicated preprocessing workflow returns.
- The GaNb₄Se₈ runs and the `data/stog_tests` family (Mn₃Sn, FeCoSn) are gitignored. Their tests —
  Auto StoG real-data detection/enforcement/ρ₀, the PCA/orientation real-data checks, the 5 K
  triplets cross-checks — skip in CI; only the committed demo run and synthetic fixtures run there.
- With Live Data **off**, a Flask analysis page picks up a newly saved configuration only on its
  next request (press Load, or switch Live Data on); notation.md §3c.
- Maintainer decisions deferred beyond 1.0 are listed in [docs/CHANGELOG.md](docs/CHANGELOG.md)
  (v1.0.0, "Deferred beyond 1.0").

## Next best steps

1. Commit trimmed fixtures of the GaNb₄Se₈ and `stog_tests` runs so the real-data tests run in CI
   (the committed demo run already covers the parsers, the KDE golden and plot parity).
2. Refactor `src/RMC_plot.py` and `src/RMC_3D.py` into thin wrappers with no import-time work.
   Defer `src/STOG_plot.py` unless preprocessing becomes a visible workflow again.
3. ~~Make browser `.rmc6f` parsing tolerant + diagnostic~~ **DONE in 1.0** — one validated grammar
   in both runtimes with a parse report against the header count.
4. Symmetry finder: a Python port or an spglib cross-check in CI (the finder is browser-only), and
   moving it into a Web Worker so the 2000-site / 384-operation caps can be raised.
5. Structure KDE: the physical-kernel option (isotropic in the real plane, width in Å) that removes
   the slab-layout kernel artefact; oblique-slice normalisation; an `[uvw]` input mode.
6. χ² history: plot every χ² column and the total, and have the watchdog classify on the total (today
   the chart and the badge follow the last `.log` column only).
7. Displacement Directions: decide the default resolution (Auto instead of ν = 10), and remove or
   rename the legacy `significance` / `peakZScore` payload fields before an API freeze.
8. PCA: identical subsamples in both engines above 20 000 copies (port the JS draw to Python, or drop
   the per-site cap).
9. Bond Geometry: pair A = C distinct-window requests as (w12 bond, w23 bond) so the work, not only
   the angle count, is bounded by `APP_MAX_ANGLES`.
10. Add the z-distribution histogram and global x-z projection panels to match `src/RMC_KDE.py`;
    `/api/project/scan` summaries.
11. ~~Thermal-ellipsoid ("PCA_KDE") view~~ **DONE** — `rmc_toolkits/pca_kde.py` (source of truth) +
    `workers/pcaKde.js` (static mode), `/api/pca/sites` + `/api/pca/kde`, and the **PCA Ellipsoid**
    tab (`PcaKdePage.jsx` + `marchingCubes.js`). Possible follow-ups: element-pooled clouds in the
    picker (the engine already supports `element=`), a per-site ellipsoid overlay in the main
    structure view, PNG/CSV export of the volume, U_cif/β_ij beside the Cartesian tensor, and
    comparing two sites side by side. The displacement-*orientation* view (`orientation.py` +
    `workers/orientation.js` + `/api/pca/orientation`, UI in `OrientationPage.jsx` /
    `OrientationView.jsx`) landed 2026-07-24. Reference: Maksim Eremenko's PCA_KDE utilities at
    <https://github.com/MaximEremenko/Utilities/tree/main/RMCProfileUtilities/PCA_KDE>.

See [docs/ROADMAP.md](docs/ROADMAP.md) for the full phased plan.
