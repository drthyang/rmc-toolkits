# Development Log

Chronological record of notable changes, newest first. For current architecture and conventions see
[AGENTS.md](../AGENTS.md); for forward plans see [ROADMAP.md](ROADMAP.md).

## Unreleased

Auto StoG page redesign + **Phase 5: the static-mode engine** — the hosted app now runs
Auto StoG entirely in the browser.

- **Static-mode Auto StoG (plan Phase 5)**: `src/workers/autoScale.js` is a straight JS
  port of the Python engine — transforms (trapezoid sine FT, Lorch, omitted-low-Q
  correction, Fourier filter, first-peak zeroing), the level sweep (prefix-sum OLS with
  numpy `.astype(int)` edge parity), the Huber-IRLS closed-form fit, sweep/joint/FZ
  amplitude modes, despiking, diagnostics, the stog parsers (`stog.inp`, xy files, `::`
  headers, Fortran-style writer), and the Faber-Ziman calculator (Sears table + formula
  parser). `autoScaleWorker.js` runs it off-thread with transferable buffers. **Parity is
  tested, not assumed**: `tests/generate_autoscale_fixture.py` freezes Python golden
  numbers into `src/__tests__/fixtures/`, and `autoScale.test.js` (11 vitest tests)
  asserts the JS fit matches to 1e−6 relative (level sweep to 1e−9, manual pipeline
  samples to 1e−9). Verified live: the in-browser fit on the FeCoSn 199 K run reproduces
  the backend exactly (a = 1.1839, b = −0.19804, 3 iterations, level 1.0120 ± 0.017).
  Static export packs the classic 9-file family into a zip (client-side `writeStogXy`,
  provenance JSON included); the Auto StoG tab now shows in both runtimes and
  `browserData.isSupportedFile` admits `.inp`.
- **Page redesigned on the app's design language, laid out for 16:9** (was hardcoded-hex
  cards in a narrow sidebar): a PcaKdePage-style horizontal controls bar (SOURCE picker +
  micro-labeled parameter fields + Auto-scale + Advanced pill), an expandable advanced
  bar (windows/grid, sweep-vs-joint, density-vs-FZ, toggle pills, fixed-(a, b) expert
  run), a stat-card readout strip (correction vs hand values, convergence with the
  per-iteration a-trajectory, high-Q level ± uncertainty, fit quality, density-limit
  verdict, concordance), then a full-width S(Q) card over side-by-side G_K(r)/D(r)
  cards — all on `index.css` tokens (borders, shadows, accent, pills, tabular numerals).
- **InteractivePlot gains guide-line support** (additive; Dashboard untouched): series
  with `role: 'guide'` render dashed/muted outside the palette rotation and are skipped
  by hover snapping (dashed legend swatches); `defaultHidden` series start muted;
  `initialYDomain` sets the un-zoomed default view (used to keep the G_K low-r level
  readable instead of the first peak); swapping `plotData` now resets zoom/hover/hidden
  state (previously stale zoom survived a re-fit).
- **Plot content**: measured-unscaled S(Q) ships default-hidden (one legend click away),
  the S → 1 asymptote, the measured level (drawn only over its admissible window), and
  the S(0) Faber-Ziman target marker join S(Q); theory lines −⟨b⟩² and −4πρ₀⟨b⟩²r anchor
  the G_K/D(r) cards.
- **Form & workflow**: session persistence (source + all settings survive a reload via
  sessionStorage), the Formula field is labeled *(neutron)* with an x-ray tooltip, an
  inline guard disables Auto-scale when FZ mode lacks ⟨b²⟩, and micro-labels no longer
  pass through `text-transform: uppercase` (which had turned ρ₀ into a capital-rho
  P-lookalike).

Auto StoG Phases 3–4 — the Flask scaling API and the **Auto StoG** page.

- **`/api/scaling/preview` + `/api/scaling/run`** (`web_app/backend/app.py`): POST endpoints
  driving the Phase-1 engine server-side. `preview` resolves a classic `stog.inp` or a bare
  data file (+ overrides; `.dat`-header ρ0/r0 and `formula` defaults mirror the CLI), runs the
  auto-fit (or a fixed manual scaling) behind a per-(path, mtime, config) LRU cache, and
  returns the full series (raw/scaled/filtered S(Q), G_K, D(r), enforced variants), theory
  guide values, diagnostics, and provenance; `{"inspect": true}` is the cheap no-compute form
  used to pre-fill the page. `run` writes the classic output family through the *same* writer
  as the CLI (`ft.dat` included) with the no-clobber guard mapped to HTTP 409, outputs
  restricted to the configured data roots. 5 new backend tests (synthetic run under
  `results/`; 20-test backend module green).
- **The Auto StoG tab** (`AutoStogPage.jsx` + CSS): first in the tab row (Dashboard remains
  the default page), Flask mode only — static mode hides the tab and the page shows a
  pointer to the local app. Automation-first per the plan: pick a source file (stog.inp
  pre-fills everything; data files pre-fill from the `.dat` header), one **Auto-scale**
  button, then a diagnostics readout (a/b beside the stog.inp hand values, convergence +
  level ± uncertainty, low-r rms, one-sided density-limit verdict, amplitude concordance
  with a "check ρ₀" hint on discord) over three InteractivePlot charts with theory
  guide-lines — S(Q) (asymptote + measured level), G_K(r) (−⟨b⟩²), D(r) (−4πρ₀⟨b⟩²r) —
  the latter two zoomed to the low-r region, enforced curves overlaid when enforcement is
  on. Advanced panel: r-windows, Lorch/despike/robust/low-Q-correction/σ toggles, sweep vs
  joint architecture, density vs FZ amplitude criterion, and fixed-(a, b) expert runs.
  Export card writes the RMCProfile-ready family (default `autoscale/` beside the input,
  Force required to overwrite, 409 surfaced as a hint). Verified end-to-end against the
  FeCoSn 199 K run in the live app: page auto-fit reproduces the CLI exactly
  (a = 1.1839, b = −0.1980, level 1.0120 ± 0.017, concordance 1.06).
- File browser (`/api/files`) now also lists `*.inp`, `*.sq`, and `*.dat` so scaling
  sources are pickable.
- Docs-on-completion pass: ROADMAP Phase 7 "Deferred preprocessing" → **active Auto StoG**
  feature; AGENTS.md architecture map gains the four scaling modules, the endpoints, and
  the page. Remaining stretch: plan Phase 5 (static-mode Web-Worker port).

Auto StoG Phases 1–2 — automatic total-scattering data scaling engine (`rmc_toolkits.scaling` +
`rmc_toolkits.transforms`) and the `rmc-autoscale` CLI, plus the app rename to
**RMCProfile Workbench**.

- **Faber-Ziman amplitude mode** (idea: Tsung-Han Yang): `ScalingConfig(
  amplitude_criterion="fz")` / CLI `--amplitude fz` implements the "subtract the measured
  high-Q level, scale Q→0 onto S(0) = 1 − ⟨b²⟩/⟨b⟩², shift the level back to 1" procedure —
  closed form on top of the level sweep, no self-consistent loop (the criterion is
  filter-independent), requires ⟨b²⟩. Because it never touches ρ0, it is the natural
  cross-check for the density-limit amplitude: on the FeCoSn validation data a ±10% ρ0
  error moves the density amplitude ~1:1 while the fz amplitude is bit-identical, so the
  concordance diagnostic turns a wrong `NUMBER_DENSITY ::` into a measurable discord. In
  fz mode `diagnostics_summary` suppresses the (vacuous) self-concordance and the density
  residuals act as the independent check.
- **FeCoSn 199 K validation + robustness campaign** (`data/stog_tests/199K`, local; script
  `data/stog_tests/robustness_199K.py` + results JSON): a third complete classic-Fortran
  run now validates the stack — manual parity `scale.fq` 6.8e−14 max|Δ|, `ft.dat`/
  `scale_ft.sq` 2.1e−5 rms, enforced `scale_ft_rmc.gr` 9.6e−5 rms. 63-case perturbation
  study: affine pre-corruption invariance exact (≤2.6e−11%); recovered scale stable ±3%
  over Qmax ∈ [18, 28] and Qmin ∈ [0.5, 1.6]; noise graceful (±0.3% at σ = 0.005, ±2.6%
  at σ = 0.05); spikes −19% unflagged→flagged, `--despike` restores to +4.3%; every
  catastrophic case (rolloff Qmax, starved Qmin, spikes) raised
  `density_limit_satisfied=False`. All three auto criteria (sweep+density, joint, fz)
  land 8–13% above the colleague's hand scaling while improving the honest low-r
  residual (auto 0.110 vs hand 0.165, −33%), with density/fz concordance 4.6%. The
  x-ray fixture tests now select the first available FeCoSn run (`100K` or `199K` — same
  stog parameterization), and a CLI-level x-ray parity test joins the suite. 139 tests.
- **Auto StoG Phase 2 — the `rmc-autoscale` CLI** (`rmc_toolkits/scaling_cli.py`; console
  entry via `[project.scripts]`, module form `python -m rmc_toolkits.scaling_cli`): drop-in
  replacement for an interactive classic-stog session. Reads a classic `stog.inp` — or
  `--data FILE --qmin --qmax` with `--formula`-computed coefficients (scattering.py) and
  ρ0/r0 pre-filled from the `.dat` `NUMBER_DENSITY ::`/`MINIMUM_DISTANCES ::` header — and
  auto-fits (a, b) by default; `--manual` reruns the stog.inp hand scaling and
  `--scale`/`--offset` fix them explicitly. Writes the seven classic output files (scaled
  S(Q), unfiltered g−1, filtered S(Q), filtered g−1 + D(r) companion column, and the RMC
  FK/GK/D(r)) plus a provenance JSON carrying the full configuration, fit history, and
  `diagnostics_summary`. Safety per plan §5: outputs default into an `autoscale/` directory
  beside the input and nothing is overwritten without `--force`, so the tool can never
  silently clobber the real STOG outputs a `stog.inp` sits beside. The RMC files get the
  exact Fortran `first_peak_zero` enforcement by default in stog.inp mode (`--no-enforce`
  opts out); the honest pre-enforcement low-r rms is always printed. `write_stog_xy` gains
  an optional third column for the `scale_ft.gr` layout. The classic fixed-name `ft.dat`
  Fourier-filter correction is written too (data-grid), so the output family matches a
  stog/pystog session file-for-file. Tests: `tests/test_scaling_cli.py`
  (9: synthetic auto/manual/data-mode end-to-end, no-clobber guard, error surfaces,
  module-entry smoke, skip-if-absent Fortran parity). Full suite: 135 tests green.

- **New `rmc_toolkits/transforms.py`**: Keen-2001-convention conversions (S(Q) ⟷ F(Q) ⟷ F_K(Q);
  g(r) ⟷ G_PDF(r) ⟷ G_K(r) ⟷ D(r)), the trapezoid sine-FT pair, Lorch window, the analytic
  omitted-low-Q correction (affine-basis form), the classic stog/pystog Fourier filter, and the
  stog low-r enforcement stage. Discretization validated against pystog 0.6.7 and a complete
  classic Fortran stog run (`data/stog_tests/stog_59438`, local-only): filter correction and
  filtered S(Q) agree to ~6e-4 rms; enforced RMC outputs match the Fortran files to 1e-9.
- **New `rmc_toolkits/scaling.py`**: the auto-scaler. Affine correction `S_corr = a·S + b`
  (multiply convention; classic stog's `yoffset/yscale` map via `a = 1/yscale`), fitted by a
  closed-form linear least-squares against the high-Q asymptote (Keen Eq. 21) and the low-r
  density limit (Eqs. 15/29 in g-space), inside a self-consistent loop with the Fourier filter.
  Recovers known (a, b) on synthetic data to ~0.3% (the omitted-low-Q correction, on by
  default, is what makes that possible — 8% bias without it). `diagnostics_summary` reports the
  honest pre-enforcement residuals and flags datasets whose missing low-Q information makes the
  absolute scale unrecoverable from self-consistency (the 59438 example is such a case: its
  filtered outputs violate the Krogh-Moe sum rule ~26× even at the expert's hand scaling).
- **Parsers**: `StogInput`/`read_stog_inp` (classic 23-line stog.inp, with explicit
  `NotImplementedError` on unexercised variants), `read_stog_xy` (count-header/NaN-tolerant),
  `read_dat_header` (`TITLE ::` / `NUMBER_DENSITY ::` / `MINIMUM_DISTANCES ::`), and
  `write_stog_xy`.
- **Tests**: `tests/test_transforms.py` + `tests/test_scaling.py` (30 tests): synthetic
  round-trips and known-scale recovery always run; Fortran-run parity tests skip cleanly when
  the local example is absent. Full suite: 95 tests green.
- **New `rmc_toolkits/scattering.py` — Faber-Ziman coefficient calculator**: bound coherent
  neutron scattering lengths for 89 natural elements (NIST NCNR / Sears 1992; real part for
  the complex-b absorbers B, Cd, Dy, Eu, Gd, In, Sm), a chemical-formula parser (decimals,
  parentheses: `"Sr0.5Ba0.5TiO3"`, `"Al2(SO4)3"`), and `faber_ziman()` returning ⟨b⟩² (the
  stog "Faber-Ziman coefficient") and ⟨b²⟩ in both barns and fm² — the ecosystem mixes units
  (pystog's argon example is fm²; classic stog inputs are barns). Cross-validated against
  pystog's argon config (3.644 fm², exact). Per-element overrides support isotopic samples;
  null-matrix compositions (⟨b⟩ ≈ 0) are rejected with a clear error.
- **App rename**: "RMCProfile Run Monitor" → **"RMCProfile Workbench"** (header, tab title,
  READMEs), reflecting the multi-tool scope ahead of the Auto StoG page.
- **Adversarial review hardening** (15-agent verified review; 11 confirmed findings fixed):
  Q ≤ 0 grid rows are cropped and `fourier_filter` rejects non-positive grids (NaN poisoning);
  the omitted-low-Q correction returns exactly zero when data start at Q = 0 (pystog parity —
  no double counting); the Lorch-branch removable singularity at r = π/Qmax gets its analytic
  limit; `np.trapezoid` shim for NumPy < 2.0; `nr`/`rmax`/`yscale` validation;
  `read_stog_xy` picks the modal column count (numeric headers can't eat the data);
  `fk_to_sq` exported. Most important: the diagnostics flag is now the one-sided
  `density_limit_satisfied` — verification *demonstrated numerically* that a smooth missing-
  low-Q deficiency is silently absorbed into a ~9–21% biased scale with all residuals clean,
  so **False proves the absolute scale is unrecoverable, but True does not certify it**.
- **Level sweep — a criterion-driven answer to "what Q is high enough?"** (idea:
  Tsung-Han Yang). `level_sweep()` searches every candidate high-Q window (both edges swept,
  O(1) per-window fits via prefix sums); a window is *admissible* iff its slope is
  statistically zero given its own fit noise (no hand-set tolerance), the minimum-variance
  admissible window wins, and the level spread across all admissible windows is the honest
  level uncertainty. End artifacts exclude themselves — the criterion independently
  rediscovered both experts' hand cuts (FeCoSn: 24.5 vs hand 26 with the rolloff onset
  caught earlier; PG3: 28.8 vs hand 28). `autoscale` now defaults to the **sweep-anchored
  architecture** (`c1_mode="sweep"`): offset tied to the measured level (`b = 1 − a·level`),
  leaving the density limit a single amplitude dof — the "shift by the level, then scale"
  decomposition, which removes the 2-dof level/amplitude trade-off pathologies and converges
  ~4× faster. No flat window → `asymptote_found=False` and automatic fallback to the joint
  fit. Reported in provenance and `diagnostics_summary`.
- **Dual amplitude criteria with a concordance diagnostic**: alongside the density-limit
  amplitude, `amplitude_from_fz_limit()` independently estimates the scale from the Q→0
  Faber-Ziman limit (Keen Eq. 21: `S(0) = 1 − ⟨b²⟩/⟨b⟩²`, robust low-Q extrapolation;
  requires `b_sq_avg`). `diagnostics_summary` reports both, their ratio, and an
  `amplitudes_concordant` flag — agreement is evidence the absolute scale is trustworthy;
  disagreement *quantifies* how much the data cannot decide it (FeCoSn: 12% discord).
- **Robust high-Q level fitting**: Huber IRLS re-weighting of the joint fit (default on;
  per-block MAD scaling so C1 and C2 are each protected against isolated outliers), optional
  per-point `sigma` 1/σ-weighting of the high-Q rows, an experimental `c1_slope_nuisance`
  term absorbing linear tail drift in the level estimate, and opt-in rolling-median
  `despike` for detector-glitch contamination — measured to restore clean 0.3% recovery
  under tail spikes that otherwise ring through the transform into the low-r window (a
  channel row re-weighting cannot reject: ~80% scale error without despiking). Despike stays
  OFF by default because it also flags real Bragg maxima on crystalline data (12% of points
  on the 59438 benchmark); the dropped count is reported in provenance (`n_despiked`).
- **Second validation dataset — FeCoSn 100 K x-ray run** (`data/stog_tests/100K`, local):
  exercises the normalized-S(Q) conventions (⟨b⟩² = 1, flat −1 level, `1.0 0 0`
  enforcement). Fortran parity to 7.9e−6 rms through the filter stage; the auto-scaler lands
  within 8% of the expert's hand values and improves the low-r residual by 26% — the
  density limit is satisfiable here (Qmin = 0.5 Å⁻¹), unlike the neutron 59438 case, and the
  one-sided diagnostic correctly reports True.
- **Exact Fortran final-step semantics** (from `stog_new3.f90`, located during verification):
  `first_peak_zero()` implements the real ripple-removal rule — zero g(r) where r ≤ cutoff
  *and* outside the first-peak window [rmin, rmax] — which degenerates to the flat −⟨b⟩²
  replacement for the validation example's parameters. The mysterious second stog.inp
  "yoffset" is the Fortran's global "Add values" knob (`y·(1+vadd)+vadd`); still rejected as
  unsupported when nonzero.
- Plan: [STOG_SCALING_PLAN.md](STOG_SCALING_PLAN.md) (build phases, verified math spec,
  validation results).

## v0.4.0 — 2026-07-16

PCA Ellipsoid: PC ⟷ crystal reference-frame switch, and Site-ellipsoids crystal axes + reset/save.

- **The main viewport can be viewed in the principal-axis (PC) frame or the crystallographic (a, b, c) frame.**
  A **Frame** switch (PC | Crystal) in the plot header flips the axis triad, the shadow-box wall projections,
  the look-down camera buttons, and **Reset view** together — PC1/PC2/PC3 in PC mode, a/b/c in crystal mode. In
  crystal mode the box + walls switch to an orthonormal frame built from the unit cell and the SAME 3D KDE
  density is re-binned onto the a/b/c planes (`projectDensityOntoFrame`), so the wall shadows are the honest
  crystal-plane projections of the displayed density; the a/b/c rods and the look-down-a/b/c camera use the true
  cell edges. Reset view frames whichever box is active corner-on, and switching frames snaps to that frame's
  default view. a/b/c are keyed by rod color in the viewport legend (no 3D letter labels). The crystal-frame directions come from `src/pcaCrystalFrame.js`
  (`unitCellVectors`), which derives the unit-cell vectors in the shared Cartesian basis the PCA axes and
  density already live in (`src/__tests__/pcaCrystalFrame.test.js`).
- **The Site ellipsoids panel gains crystallographic axes + reset/save.** An opt-in a/b/c gizmo at the unit-cell
  origin (its own **a b c** toggle, keyed by color in the panel legend) and **Reset view** + **Save** (PNG,
  1× / 3×) controls mirroring the main panel.

Old coordinates-only `.rmc6f` files reconstruct thermal-ellipsoid sites; PCA controls regrouped with a shell contrast knob.

- **The oldest `.rmc6f` files (coordinates only) now drive PCA-KDE and Atomic Density.** That format drops the
  reference-site and per-atom cell columns entirely (`id element [label] x y z`, 6 fields), so there was no
  per-site grouping — the file failed to load at all. `parseAtomLine` now also accepts those short lines, and
  when a file carries no reference numbers the PCA parser reconstructs sites by folding every atom into one unit
  cell and clustering per element (periodic, full cell-metric minimum image, circular-mean unwrap of each
  cluster). Each reconstructed site reports its member count against the one-per-cell expectation from the
  supercell: `27/27` is a clean crystallographic site, while `162/27` flags atoms that do not resolve into
  separate sites at the chosen distance — close sites, or an orientationally-disordered group such as a rotor
  shell, whose "ellipsoid" is a shell best read from the KDE. A **Cluster** distance knob (shown only for such
  files) tunes the grouping; the site list and a badge surface the count and a clean / merged-or-disordered
  label, and the default selection prefers a clean site. Atomic Density needed only the coordinate fold, which
  no longer depends on the cell columns. New synthetic, oblique-cell, and real SF6 190 K (rotor-phase) cases in
  `src/__tests__/rmc6f.test.js`. (The Flask `/api/pca` backend does not reconstruct; the browser worker path does.)
- **PCA Ellipsoid controls regrouped for a clearer layout.** The wireframe and the density painted on its
  surface now share one **Ellipsoid** group (Wireframe · Level · Color · Shell · Colormap · Contrast), with
  every toggle labeled by what it toggles. A new **Contrast** knob stretches the KDE-shell colormap around its
  mid-tone — a single control over the effective vmin/vmax — so faint departures from the harmonic ellipsoid
  stand out. Turning on the **Isosurface** now also clears the wireframe and shell for a clean volume view.
- **The Displacement-statistics panel no longer clips.** It now sizes to its content (so the full covariance
  matrix, all three principal axes, and the summary line always show), and the Site-ellipsoids 3D view below it
  re-fits to the remaining space even though its PCA arrives asynchronously from a worker.

Dashboard plots labeled by the fit function declared in the run-control `.dat`.

- The RMCProfile run-control `<stem>.dat` records the correction / fit-function form per dataset
  (`> DATA_TYPE :: G(r)` with `> FIT_TYPE :: D(r)` means the data is fit as D(r)). The dashboard now uses
  that fit type as the plot's heading and y-axis label, so a `.gr` file fit as D(r) reads **D(r)**, not
  G(r) (F(Q), S(Q), … likewise). `parseRunSettings` already extracted it; a new `fitTypeByFilename` maps it
  by file name and the run assembly pairs it onto the matching plot file.
- Finding the right `.dat` among the many in a run folder (chi2.dat, optimization.dat, …) uses the existing
  stem match (`<rmc6f-stem>.dat`) first, then a capped content scan of the other `.dat` files as a fallback,
  reading only each file's head so a large data `.dat` is never read in full.
- Any `.gr` / `.sq` / `.fq` STOG data file now loads as a dashboard plot, not just the default `scale_ft.*`
  names — runs commonly use descriptive data-file names (e.g. `PMN_300k_rmc_..._v2.gr`).

Robust `.rmc6f` parsing, dashboard box-zoom, Atomic Density render fix, and PCA polish.

- **Older `.rmc6f` files now load.** The only structural difference from the current format is the per-atom
  type label (`id element [type] x y z ref cx cy cz`, 10 fields); older 2018-era files omit it (9 fields).
  Atom-line parsing now indexes the reference number, cell indices, and coordinates from the END of the line,
  tolerating any number of label columns. A new shared `web_app/frontend/src/rmc6f.js` (`parseAtomLine`) backs
  both the structure and PCA browser parsers; Python `iter_rmc6f_atoms` (used by `pca_kde` and `kde`) gets the
  same treatment. Covered by `src/__tests__/rmc6f.test.js` and `OldFormatRmc6fTests` in `tests/test_parsers.py`.
- **Dashboard charts support box zoom.** Drag a rectangle to zoom into that region on both axes (a thin
  horizontal or vertical drag still zooms just that axis); series are clipped to the plot area, and Reset
  zoom / double-click restore the full view.
- **Atomic Density first-visit render fix.** The KDE-slice / slab canvases and the 3D model measured their size
  once on mount with no observer, so a first visit before layout settled could leave the 2D panels below
  display resolution and the 3D model at 0×0. They now re-measure via a `ResizeObserver` (2D redraw; 3D updates
  the camera + renderer in place), matching the pattern the PCA viewport already uses.
- **PCA Ellipsoid polish:** a **Black** wireframe option; the wireframe is drawn opaque while the KDE shell is
  on (crisp cage over the colored surface); the **KDE shell** toggle's `?` moved out of the switch `<label>`
  so the toggle click is no longer swallowed by the help badge, and its popover opens leftward (`align="end"`);
  Site-selector markers are drawn as each site's calculated thermal ellipsoid, with the selection glow + triad
  scaled to the marker so a soft site can't outgrow them; and the U_iso / B_iso subscripts render correctly
  (the stat label was an inline-flex box that dropped `vertical-align` on `<sub>`).

PCA Ellipsoid: selectable / KDE-projected ellipsoid, Site selector; taller dashboard charts.

- **KDE shell** (new, optional): a "KDE shell" toggle paints the KDE density onto the p% ellipsoid surface
  (per-vertex trilinear sample of the volume, current colormap), so anharmonic departures from the harmonic
  ellipsoid read as hot/cold patches on the shell. The coloring is stretched to the shell's own density range
  (not the global 0..vmax), so the small variation across a near-iso-probability shell reads at full contrast
  (a Gaussian cloud stays near-flat, as it should). It shows the same density as the isosurface from the
  outside and would occlude it, so the two are mutually exclusive — enabling one switches the other off, and
  both-off is allowed (wireframe + wall projections only).
- **Ellipsoid wireframe color** is now chosen from the controls bar (swatch + dropdown: Amber, White, Silver,
  Cyan, Violet), kept transparent so it never crowds the isosurface. The default is Amber (`#ff7a1a`) — the
  same warm tone as the selected-site highlight in the Site selector, so the cage reads as "this site."
- The **Unit cell** side panel is renamed **Site selector**, since its job is picking the reference site to
  analyze (click an atom to load its PCA-KDE).
- Dashboard charts fill the panel width on wide (≥ 1500 px, 16:9) screens instead of being capped short and
  letterboxed with empty side margins: the 8:5 plot's height cap goes 360 px → 540 px, so it grows to fill the
  card at the same aspect (e.g. 1080p: 576×360 with side gaps → 618×386 filling). Fills fully through 1440p.

PCA Ellipsoid statistics detail and a 16:9 layout; Demo button beside the page tabs.

- The **Displacement statistics** panel now shows the full covariance tensor U (Å², Cartesian x/y/z with the
  variance diagonal highlighted) and a **Principal axes** table — each PCA component (eigenvector) as a row
  in Cartesian x/y/z with its eigenvalue λ (Å²), RMS amplitude (Å), and per-axis excess kurtosis κ. PC1/PC2/
  PC3 rows are color-keyed to the 3D triad and the viewport legend. The old single-line "RMS axes" and
  "Eigenvalues" rows fold into that table (now per-PC columns); U_iso / B_iso / Anisotropy / Non-Gaussianity
  stay as the scalar summary. Both the Flask and browser-worker paths already return `covariance` / `axes` /
  `eigenvalues` / `rms` / `excessKurtosis`, so the detail appears in every runtime mode.
- The page now fills the screen height instead of sprawling wide. The 3D viewport takes three-quarters of the
  width on the left; the statistics and unit-cell panels stack as two **equal-height** rows in the remaining
  quarter, with the viewport spanning both so all four outer edges line up. The grid flexes to fill the page's
  leftover height, so the view roughly fills a 16:9 screen. The four scalar figures sit in a 2×2 grid and the
  statistics body scrolls within its row if a small window can't show every table; below 980 px it stacks into
  one column at natural heights and the page scrolls.
- The header **Demo** toggle moves into the left cluster just after the AI Assistant tab (beside the page
  tabs), leaving Live Data + the run-folder field as the right-hand data-source group.

PCA Ellipsoid: reset the main view on a new dataset, and a deterministic camera reset.

- Loading a different run now returns the main 3D panel to the default body-diagonal view instead of
  inheriting the previous model's orbit. The reframe is keyed on the run's identity (App's `runId`),
  so it fires even when the new run's first site shares a reference number with the old one, while a
  Live Data refresh of the same run and any slider/layer tweak leave the view untouched.
- Root-caused and fixed why the reset sometimes landed rotated off the default: `frameMainCamera` set
  the pose and then called `OrbitControls.update()` once, which applied the control's leftover damped
  orbit velocity to the new pose. It now flushes that velocity first (a damping-off `update()` that
  zeroes the pending delta), so the reset is deterministic regardless of prior momentum — the same
  reason it also fixes the "Reset view" button landing slightly off right after a drag. This let the
  reframe live in one place (the scene rebuild) and removed the visibility-tracking scaffolding.

PCA displacements in the AI Assistant context.

- The run context sent to the local model now carries a `pca_displacements` section: one entry per
  reference site from the PCA Ellipsoid analysis (isotropic ADP `U_iso_A2`, the three principal RMS
  amplitudes `rms_axes_A`, `anisotropy`, and mean-excess-kurtosis `non_gaussianity`), ranked most
  non-Gaussian first so the anharmonic / split sites lead and survive budget trimming. This is
  information the model previously couldn't see — symmetry gives only mean displacement per Wyckoff
  orbit, not the anisotropy or the anharmonicity. The system prompt gains a matching bullet so the
  model interprets the kurtosis correctly, and the character-budget ladder trims the least-anharmonic
  PCA sites (then their explanatory note) alongside the other evidence.
- Wiring: the PCA Ellipsoid page publishes its computed per-site table upward (`onSitesChange`), App
  holds it, and the AI Assistant page threads it into `buildRunContext`. The section is present once
  the PCA Ellipsoid page has been opened for the run (its analysis runs there), and follows the active
  dataset. Covered by four new cases in `src/llm/__tests__/runContext.test.js`.

PCA Ellipsoid page refinements.

- The main 3D view has a **Reset view** button in the panel header that returns the camera to the
  default body-diagonal framing after you orbit or zoom away — the same framing applied automatically
  when the site changes, now available on demand.
- Fixed the PCA Ellipsoid not recomputing for a newly loaded dataset. The static-mode worker cached
  parsed displacement clouds by file *path*, so a different run reused the previous model's clouds. The
  cache is now content-addressed (keyed on a signature of the `.rmc6f` text), and the page ties the
  loaded text to the file it came from, so requests never run against the previous dataset. Covered by
  `src/workers/__tests__/pcaKdeWorker.test.js`.

PCA / thermal-ellipsoid KDE computation engine.

- New engine turns per-site RMC displacement clouds into anisotropic displacement tensors (the
  thermal ellipsoids) and a smooth 3D probability density. An atom's offset from the average
  structure is `coords − cellIndices/supercell` folded over the supercell boundary, mean-subtracted
  per reference site, then mapped to Cartesian Å; the covariance of each site's cloud is its ADP,
  and a Gaussian KDE of the cloud is the density the ellipsoid only approximates. PCA and KDE are
  standard tools; the specific per-site analysis and the shadow-box visualization below follow
  **Maksim Eremenko's PCA_KDE utilities**
  (<https://github.com/MaximEremenko/Utilities/tree/main/RMCProfileUtilities/PCA_KDE>) — we followed
  that approach and reimplemented it independently into the toolkit's dual-mode architecture, not a
  port of his code (his `KDE.js` evaluates a full multivariate Gaussian KDE; ours factorizes it).
- **The KDE is separable, and exact.** Sampling on a grid aligned with the cloud's principal axes
  makes SciPy's `H = factor²·C` bandwidth diagonal, so the 3D Gaussian factorizes into three 1D
  kernels and the volume becomes their tensor product: `N·3·grid` exponentials contracted through
  BLAS instead of `N·grid³` kernel evaluations. Measured **36–51× faster than a naïve
  `scipy.stats.gaussian_kde` sweep of the same volume, with max abs error ~1e-13** (machine
  precision) — it is the same estimator, not an approximation. A 52-site, 52 000-atom configuration
  builds all clouds in ~80 ms and all 52 ellipsoids in ~1 ms; a single-site 48³ volume solves in
  ~7 ms server-side.
- `rmc_toolkits/pca_kde.py` is the source of truth: `load_site_displacements`, `site_ellipsoids`
  (batched, one pass), `pca_kde_volume` / `site_pca_kde` (volume + PC-plane projections + iso
  thresholds by enclosed probability mass and by raw density). Axes are sign- and
  handedness-canonicalized for reproducibility; a floored eigenvalue keeps flat/linear clouds
  non-singular and flags them `degenerate`.
- `web_app/frontend/src/workers/pcaKde.js` is a straight JS port for static mode (a 3×3 Jacobi
  eigensolver stands in for `eigh`), driven by `pcaKdeWorker.js` off the main thread — ~96 ms for a
  1000-point 48³ volume in-browser, no GPU needed. `tests/test_pca_kde.py` and
  `src/workers/__tests__/pcaKde.test.js` both assert the separable volume equals the full
  brute-force estimator to round-off (15 + 14 tests), plus site extraction, supercell-boundary
  folding, mass normalization, and known-anisotropy recovery.
- Flask exposes `/api/pca/sites` (per-site ellipsoid table) and `/api/pca/kde` (one site's or one
  element's volume), both behind a per-(path, mtime) LRU cache.
- **PCA Ellipsoid page** (new top-level tab; the "KDE / 3D" tab is renamed "Atomic Density"). A site picker over all reference sites, an
  ellipsoid summary table, a Three.js scene showing the KDE isosurface (extracted by a self-written
  marching-cubes module — `workers/marchingCubes.js`, since Three's `MarchingCubes` only builds
  metaballs) nested with the p% thermal-ellipsoid wireframe and PC-axis triad, and the three
  PC-plane KDE projections as heatmaps. An isosurface-mass slider sweeps the enclosed-probability
  threshold. Works in both runtimes: Flask endpoints in backend mode, the `pcaKdeWorker` +
  `pcaKde.js` engine off the main thread in static mode (verified on the bundled Demo run).
- **Non-Gaussianity (excess kurtosis)** is reported per site. It quantifies why a KDE isosurface can
  sit inside its harmonic ellipsoid: the covariance is inflated by fat tails while the isosurface
  tracks the peaked core. Verified on `data/5K_try1` — Nb sites read ~5–7 (strongly anharmonic,
  isosurface well inside the ellipsoid), near-isotropic Ga sites read ~1 (isosurface ≈ ellipsoid),
  and a synthetic Gaussian cloud reads ~0. Marching cubes has its own sphere/normals tests; the
  isosurface-vs-ellipsoid scaling was validated to agree with theory (KDE surface ≈ 1.06× the
  ellipsoid for a true Gaussian).

Periodic boundary conditions for the KDE slice.

- The KDE under-counted near cell boundaries: folding atoms into one unit cell drops their
  periodic neighbors, so the density decayed toward faces/edges/corners (an edge blob showed
  roughly half the interior amplitude, a corner roughly a quarter), and a slab centered on the
  z=0 face missed the atoms just below z=1 entirely.
- Both KDE paths now tile periodic images from the 26 neighbor cells within a margin that covers
  the kernel reach and the slab depth (`min(0.5, max(0.1, 2*bw, thickness))`). In
  `rmc_toolkits/kde.py`, `_augment_periodic_images` feeds `oriented_kde_slice`, and `kde_slice`
  reports `slabCount` as unique source atoms and rescales the density (SciPy's `gaussian_kde`
  divides by every fit point, images included). In `localKdeWorker.js` the same augmentation runs
  before `makeSlab`, and the image factor rides in the kernel normalizer, so the WebGPU path
  inherits the fix unchanged. `kde_slice` called directly (no `source_index`) keeps the old
  truncated behavior for generic point clouds.
- Verified on `snao.rmc6f` at a z=0 slab: slab population 3110 → 6034 atoms, corner density
  20.8 → 55.4, edge-column peaks now match interior peaks (~73 vs ~75), opposite edges agree to
  ~5% (fit-subsampling noise). New tests: in-plane wrap symmetry and depth-wrap slab selection in
  `tests/test_kde.py`, plus a mirrored vitest suite for the worker (whose `onmessage` registration
  is now guarded so tests can import the module outside a worker).

Source data file names on figures (user-experience feedback).

- Every dashboard plot card shows its source file name in monospace under the title — the two
  otherwise-identical `EXAFS Q-space` cards are now distinguishable — with the full path on hover.
- The collapsed R-value strip shows its source log name(s); a combined multi-log strip lists every
  parsed log.
- The Model information title cell now displays the source `.rmc6f` file name (previously hidden
  in a tooltip), which also covers the KDE / 3D panels since they all derive from that one file.
- File names sit outside the card `h3`, so "Save all figures" export names are unchanged.

Math rendering in AI Assistant replies.

- Assistant Markdown now typesets LaTeX: integrated `remark-math` + `rehype-katex` so
  dollar-delimited math (e.g. `$R_{wp}$`, `$\chi^2$`, `$\text{\AA}$`) in replies renders as real
  math instead of literal text — including inline math inside GFM table cells, which is how models
  tend to format result summaries. KaTeX runs *after* `rehype-sanitize` (raw → sanitize → katex) so
  model HTML stays fully sanitized while KaTeX emits its own trusted markup, and KaTeX fonts are
  bundled locally so the static GitHub Pages build needs no CDN.

Demo run and refreshed screenshots.

- Added a header **Demo** toggle that loads a bundled GaTa4Se8 250 K example run (under
  `web_app/frontend/public/demo/`) so first-time visitors see a populated dashboard; a second click
  clears it.
- Recaptured the Dashboard and KDE/3D screenshots against the demo run and added an AI Assistant demo
  GIF (a local Ollama model summarizing the run as a table with inline LaTeX math). Renamed the screenshots to
  `assets/rmc-toolkits-dashboard-demo.png` / `assets/rmc-toolkits-kde-demo.png` (new paths bust the
  stale README image cache) and removed the unused legacy screenshots from `assets/`.

Relicensed to AGPLv3.

- Switched the project license from MIT to the **GNU Affero General Public License v3.0** to keep
  the code open while requiring that modified versions — including those run as a hosted network
  service (AGPL §13) — release their source. Updated `LICENSE`, `pyproject.toml` (license field +
  trove classifier), the README badge/notice, and the in-app footers.
- Added an **About & documentation** link in the Dashboard and KDE/3D footers pointing to the repo,
  which also surfaces the source for AGPL §13 network-use compliance.

RMCProfile-first positioning.

- README, QuickStart, roadmap, and agent notes now frame the current app as an RMCProfile
  modeling-output dashboard. STOG is treated as deferred preprocessing/legacy parser support, not
  a current RMC modeling feature.
- The dashboard excludes STOG plot kinds from the visible plot list and assistant context while
  keeping the underlying parser support in place.

RMCProfile wording cleanup.

- Replaced RMCProfile "refinement" language in docs and assistant prompts with RMC modeling /
  atomistic configuration optimization wording, to avoid implying Rietveld-style parameter
  refinement.
- Renamed the assistant's move-counter context from `refinement` to
  `configuration_optimization`.

Gentler AI Assistant startup.

- The automatic on-load connection probe no longer shows a red "Connection failed" before the user
  has set anything up. The status reads a neutral **"Connect a model to start"** (grey dot), and the
  settings drawer stays clean; the full error + setup hint appear only after the user actively
  presses **Test**.

Assistant "Beta" badge and a cleaner chat speaker label.

- The header badge now reads **Beta** (was "Experimental") in **amber** — the complement of the
  app's cobalt accent, and the same palette as the cloud-provider notices — so the beta reminder
  stands out instead of blending into the blue theme.
- Chat messages are labeled **Assistant** rather than the raw model id (which stays visible in the
  header's model switcher), for a cleaner transcript.

Run history and run-control settings in the AI Assistant's context.

- **`configuration_optimization` block** from the `.rmc6f` header: moves generated/tried/accepted
  plus the derived **acceptance ratio** and **accepted moves per atom** — the standard gauge of
  whether the configuration has been sampled long enough — and accumulated running time.
- **`run_settings` block** from the RMCProfile run-control file. The correct `.dat` is the one whose
  basename matches the chosen structure stem (so `chi2.dat`, `optimization.dat`, … are never picked
  up). Extracts title/material/phase/temperature, **minimum distances labeled per element pair**
  (hard closest-approach constraints — the model is told a g(r) peak pinned there may be
  constraint-limited), max move sizes per element, time/save limits, flags, and the fitted-data
  blocks. Static mode only, like the rest of the browser-parsed context.
- System prompt teaches both blocks; 9 new tests (header counters, `.dat` parser, stem-matched
  selection, derived stats, pair labeling).

Chat rendering fixes and a window-filling chat box.

- Markdown now parses **and sanitizes** inline HTML (`rehype-raw` + `rehype-sanitize`), so `<br>`
  line breaks inside table cells render instead of collapsing onto one line; scripts, event
  handlers, and images are still stripped (no injection surface, no external loads from model
  output) — covered by tests.
- The AI Assistant chat box now **fills the window height** (message log scrolls internally, the
  composer stays pinned at the bottom) and the column is a bit wider (880 → 1000 px).

Rendered Markdown in AI Assistant chat replies.

- Assistant answers now render **Markdown** — GitHub-flavored **tables**, lists, code blocks, inline
  code, headings — via `react-markdown` + `remark-gfm`, so a "summarize this as a table" request
  shows a real table instead of raw `| --- |` text. Rendering builds a React tree with raw HTML
  disallowed, so the no-injection property is preserved (covered by a test); user messages stay
  plain. (eslint: enable `ignoreRestSiblings` for the `{ node, ...props }` component overrides.)

App identity — a minimalist pair-distribution "wave" mark.

- Replaced the placeholder "R" brand with a single shared g(r)-wave SVG (flat brand blue, `#2563eb`):
  a base-aware favicon (`public/favicon.svg`), the header brand mark, and the AI Assistant
  avatar / empty-state icon. Removed the default `vite.svg`.

ChatGPT-style "Thinking" indicator for reasoning models.

- Reasoning models (e.g. qwen3 via Ollama) stream their chain-of-thought in a separate field; the
  streaming client now yields structured content/reasoning chunks and the chat shows a collapsible,
  shimmering **"Thinking"** panel (live reasoning) that collapses to **"Thought for Ns"**
  (re-expandable), plus a bouncing-dots indicator for non-reasoning models slow to their first token.
  Reasoning is preserved on each turn.

Local-distortion evidence in the AI Assistant's run context.

- The context JSON gains a **`symmetry` block** — space group, the symmetry-vs-tolerance **ladder**
  (distortion magnitude + character), and the Wyckoff-orbit **sites ranked by rms displacement** —
  plus **`pair_correlations`**: each partial-PDF pair's first g(r) peaks next to the
  nearest-neighbour distance the average structure predicts, so the model can point at the sites
  and pairs participating in short-range correlations.
- Per-site rms displacements (`dispA`) fall out of the circular-mean accumulators
  `structureFromRmc6f` already had — derived per axis from the resultant length, scaled to Å; the
  symmetry orbits now carry `members` (basis indices) so the context can aggregate per orbit.
- New `llm/context/pairCorrelations.js` (average-structure NN distances over periodic images +
  a smoothed local-maxima peak finder); context budget raised to ~4.5k chars with an explicit trim
  order (extra peaks → ladder middle rungs → history → low-displacement sites → datasets, each
  recorded with an `*_omitted` count); system prompt teaches the model how to read the new blocks.
- Tests: 51 passing (peak finder, NN distances, displacement math on synthetic .rmc6f fixtures,
  symmetry-block assembly, trim order).

AI Assistant rework (shipped as PRs #3/#4, catching the log up): moved from a dashboard card to a
dedicated **AI Assistant page** (chat-only — Summary/Report tabs removed, everything happens in the
chat), modern chat UI with a persistent connection bar (status dot, model switcher, settings gear,
auto-connect), redesigned Connection Settings, and optional **cloud providers** (OpenAI, Gemini)
with Bearer API-key auth and explicit data-leaves-your-device warnings; local Ollama/LM Studio
remains the default. Ollama CORS + Safari setup documented in-app and in QuickStart.

Experimental AI Assistant — local LLM in the data-monitor pipeline.

- Added `web_app/frontend/src/llm/`, a self-contained experimental module connecting the dashboard
  to a **local LLM** (Ollama `/v1` or LM Studio, both OpenAI-compatible) **directly from the
  browser** — no server, no API keys; run-derived summaries go only to the model server the user
  runs. Four features: one-click **run summary/assessment** (streamed), **chat Q&A** with the run
  context injected, one-click **Markdown run report** (deterministic metrics tables + clearly
  labeled AI narrative), and a **live convergence watchdog** badge on the R-value card
  (slope heuristics are the source of truth; the LLM only writes the note; piggybacks on the
  existing Live Data poll).
- Pipeline seams are explicit as a learning artifact (see the module README): context builder with
  history downsampling + a ~3k-char budget → prompt templates → hand-rolled SSE streaming client →
  streaming UI with a "context sent to the model" inspector. Connection test translates failures
  into actionable CORS/setup hints.
- Dashboard mounts the collapsible **AI Assistant** card (zero network activity until used) and the
  watchdog badge; strict import boundary (props in, one `figureExport` helper out) keeps the module
  extractable to its own repository.
- Added **vitest** (first frontend unit tests, 43 passing) covering the context builder,
  convergence heuristics, SSE parsing, prompt templates, and report assembly; CI now runs
  `npm test` between lint and build.

## v0.3.0 — 2026-07-02

Symmetry analysis in the browser.

- Added a **Detected SG** card beside *Model information* on both the Dashboard and KDE/3D pages: a
  client-side, table-free **FINDSYM-like space-group finder** (ported from the RMC-phonon-dynamics
  project — no spglib/WASM). From the folded `.rmc6f` basis it detects the space group (H–M symbol +
  number), point group, and operation count, and draws an interactive **tolerance ladder** — each
  brick is the space group that holds over a range of atomic-position tolerance; click a rung to
  select it. New frontend modules `symmetry.js` (the finder) and `symmetryModel.js` (structure →
  conventional cell + basis glue); `browserData.structureFromRmc6f` now also returns a
  circular-mean per-reference-number basis.
- The space-group **number table is complete for every producible symbol** (all point-group ×
  allowed-centering combinations), and the headline always shows the point group even when a number
  is unavailable.
- The **Detected SG tolerance selection persists** when switching between the Dashboard and KDE/3D
  pages (shared via context).
- KDE/3D: **contours and log-scale density are on by default.**

## v0.2.0 — 2026-06-26

Figure export, axis labels, and browser-first positioning.

- Added per-figure **Save** controls. Dashboard charts (inline SVG) export as **PNG** (raster) or
  true-vector **SVG**; the KDE/3D panels (canvas + WebGL) export **PNG** at native or **3×**
  resolution. The earlier raster-PDF attempt was dropped — a PDF of a chart isn't true vector, and
  SVG is the proper vector format (convertible to PDF offline if needed).
- Added **Save all figures** in the dashboard's *Loaded N plot files* header. It bundles every
  visible chart into a single `.zip` (PNG or SVG), avoiding the browser's multi-download blocking.
- Added supporting modules: `figureExport.js` (SVG/canvas rasterization + standalone-SVG
  serialization with inlined computed styles), a shared `SaveMenu` component, and a dependency-free
  store-method ZIP writer (`zipArchive.js`, with CRC-32).
- High-resolution panel export re-renders the same drawing onto a 3× offscreen canvas (KDE slice,
  slab) or re-renders the Three.js scene at a higher pixel ratio (`preserveDrawingBuffer`), so the
  output is genuinely higher-resolution rather than upscaled.
- Gave proper vertical-axis labels to plots that previously showed a generic `data`: PDF/partials →
  `G(r)`, x-ray/neutron S(Q) → `S(Q)`, Bragg → `Intensity`. Applied in both the Flask
  `/api/plot/data` endpoint and the browser static-mode parser.
- Reframed the README and roadmap around the hosted browser app (GitHub Pages) as the primary,
  no-install way to use the dashboard; the Flask backend and Python package are now positioned as
  optional local/advanced paths. Added a Phase 7 plan for interactive STOG reduction and guided
  RMCProfile run setup.
- Capitalized **RMCProfile** consistently in the app title, header, and messages.

## v0.1.0 — 2026-06-21 (first tagged release)

First public, releasable version: the reusable Python package, Flask API, interactive dashboard,
server-side + browser (WebGPU/CPU) KDE, and Three.js structure viewer with a draggable slice band.

RMCProfile EXAFS dataset output parsing and large-structure display (added late in the v0.1.0 cycle):

- Added plot-kind detection for RMCProfile EXAFS dataset output files (`*-EXAFS-*_Q_OUTPUT.csv` and
  `*-EXAFS-*_R_OUTPUT.csv`) in both the Python package and browser static mode.
- Added `read_exafs_csv`, which handles Q-output files with a descriptive title row before the column
  header and R-output files with real/imaginary/modulus transform columns.
- Added backend and frontend chart labels for these dataset outputs:
  - Q-space: `k (Å^{-1})` vs `χ(k) k²`.
  - R-space: `r (Å)` vs `FT[χ(k) k²]`.
- Added dashboard file badges/order and tests for parser, plot metadata, file listing, and
  `/api/plot/data`.
- Raised the structure endpoint and Structure page point limits to 1,000,000 and raised the Slab In
  Cell canvas draw cap to match, so moderate structures such as `data/RMC/snao.rmc6f` render all
  returned atoms instead of every other point.

Release engineering added in this version:

- **License:** added an MIT `LICENSE`.
- **Packaging:** added `pyproject.toml` so `rmc_toolkits` is pip-installable (`pip install -e .`) and
  exposes `rmc_toolkits.__version__`.
- **CI:** added `.github/workflows/tests.yml` — runs the Python test suite plus frontend lint/build
  on every push and PR to `main`.
- **Tests:** sample-data-backed tests now skip cleanly when the gitignored GNSe example dataset is
  absent, so the suite is green on a fresh clone and in CI (17 run, 16 skip without the sample).

## 2026-06-21 — Draggable slice band + structure-view color work

- **Slab In Cell drag-to-move slice.** Users can now grab the highlighted band in the *Slab In Cell*
  panel and drag it along the slice axis to set the slice position (`zCenter`) live; thickness is
  unchanged. Implemented by adding an `invert()` to `makePlaneMapper` (recovers plane coordinates
  from a screen point), publishing the band geometry each render into `slabGeometryRef`, and wiring
  pointer handlers on the slab canvas (`grab`/`grabbing` cursors, pointer capture,
  `touch-action: none`).
- **Bug fix:** the drag listeners initially never attached — the effect used `[]` deps but the slab
  canvas only renders after a structure loads (`{structure && ...}`), so `slabCanvasRef` was null at
  mount. Fixed by depending on `[structure]`. Verified end-to-end against `data/5K_try1`: drag moved
  SLICE 0.39 → 0.709 with the band and in-slab atoms updating.
- Earlier in the day: gave each element a distinct color across the slab and 3D views, added a shared
  atom color legend above the structure views, and relocated the atom badges.

## 2026-06-19 — WebGPU-accelerated browser KDE

- Moved the static-mode density-map hot loop (`O(grid² · samples)`, up to ~400M `exp()` per slice)
  onto the GPU. Each grid cell is independent → one WebGPU compute-shader invocation per cell.
- Added `workers/gpuKde.js`: inline WGSL shader, lazily-cached adapter/device init, `computeDensityGpu`
  (writes buffers, dispatches, reads back via `mapAsync`, reshapes to the same `density[y][x]` grid as
  the CPU path), and a `shouldUseGpu` work-size heuristic.
- Refactored `localKdeWorker.js` into `computeDensityCpu` + a GPU-or-CPU branch; made `computeKde`
  and `onmessage` async. Same `{ id, result }` message shape, so `StructurePage.jsx` needed no
  changes (its request-`id` guard already drops stale replies). Result carries `backend: 'gpu'|'cpu'`.
- GPU used only when `grid·grid·samples >= 2_000_000`; init attempted once per worker and cached.
- **Robustness is the design point:** missing `navigator.gpu`, no adapter, device/shader error, lost
  device, or sub-threshold work all fall back to the JS loop with identical output; no GPU failure
  surfaces through `worker.onerror`.
- Verified (Chromium/Metal): GPU(f32) vs CPU(f64) parity ~1.7e-6 relative across grids 16/120/260;
  density step ~59× faster at grid 120 and ~107× faster at grid 260. Still subsamples to 6000 slab
  points and uses f32 — a visualization path, not a substitute for server-side SciPy KDE.

## 2026-06-18 — Live Data and loaded-file controls

- Added optional Live Data monitoring for the local Flask dashboard: polls the selected folder for
  supported-file changes and refreshes metadata/charts without a manual reload.
- Added file mtime + size to `/api/files` so the frontend can detect changes cheaply.
- Added a collapsed `Loaded N plot files` panel; expanding shows detected files as badges that can
  hide/show their chart.
- Aligned local folder-selector copy with the hosted dashboard (`Run folder` / `Select Folder`) and
  hardened the native picker startup so a missing default folder falls back to the nearest existing
  directory.

## 2026-06-16 — Hosted static dashboard

- Added the GitHub Pages workflow (`.github/workflows/pages.yml`); Pages source set to **GitHub
  Actions** so the built React/Vite app is served instead of a Jekyll README page.
- Added static-mode local file loading (`browserData.js`): open the hosted dashboard, select a local
  run folder, and parse RMCProfile CSV/log/STOG outputs with no Python and no upload.
- Extended static mode to parse uploaded `.rmc6f`, populate the model summary, and render the slab
  projection + Three.js 3D view.
- Added `workers/localKdeWorker.js`: off-thread browser KDE capped at 6000 slab points, deterministic
  pseudo-random sampling, contour segments for the existing overlay.
- Kept Flask `/api/kde/slice` as the reference server-side SciPy path.

## 2026-06-15 — Palette and docs

- Three.js atom palette → Nature-style (Ga deep blue, Nb vermillion, Se teal).
- Refreshed `assets/rmc-toolkits-KDE.png`; updated README + frontend docs for the local preview flow.

## 2026-05-27 — UI refresh and interactive dashboard

- Bright/dark theme system via CSS variables with a persisted header toggle.
- Removed the sidebar-first workflow; header now hosts data-path input, Dashboard/KDE nav, and theme
  toggle. App defaults to the `data` path.
- Added `GET /api/plot/data` (parsed series + normalized scientific labels: `χ`, `Å`, `Q (Å⁻¹)`).
- Replaced PNG cards with `InteractivePlot.jsx` (SVG, hover readouts, legend toggles, integer x
  ticks, drag-to-zoom + reset). Simplified the Dashboard to a three-card grid.
- Reworked KDE/3D into three side-by-side panels (XY slice, slab x-z projection, Three.js model);
  2D panels preserve lattice aspect ratios. Added gray slab-edge outlines in the 3D model.

## 2026-05-27 — Real server-side KDE

- Added `rmc_toolkits/kde.py`: loads unit-cell-folded cartesian (Å) positions from a `.rmc6f`
  (optional element filter) and computes an XY `scipy.stats.gaussian_kde` density for a z-slab
  (ported from `src/RMC_KDE.py`). Returns plain arrays (grid, extent, contour polylines via
  `contourpy`, slab count, vmin/vmax); fit subsampled to 6000 points for slider responsiveness.
- Added `GET /api/kde/slice`. `z`/`dz` are cell-edge fractions converted to Å internally; loaded
  positions cached per (path, mtime, element) with `lru_cache`.
- Made the backend port configurable via `RMC_TOOLKITS_PORT` (default 5000).
- Replaced the browser box-blur "density" with the fetched real KDE grid (colormap LUT + contour
  overlay); fetches debounced with `AbortController`. Added bandwidth/colormap/grid/contour/log-scale
  controls and `colormaps.js`. Default z-slice auto-snaps to the densest band on load.
- Added `.venv/`, `__pycache__/`, `*.pyc` to `.gitignore`.

## 2026-05-27 — Package + tests foundation

- Added `rmc_toolkits/parsers.py` (RMC CSV/log/STOG, `.rmc6f`, `Frac*.txt`) and
  `write_frac_from_rmc6f` conversion.
- Added `rmc_toolkits/plots.py` (plot detection, figures, metrics, PNG serialization); backend now
  calls into the package instead of duplicating logic.
- Added backend endpoints: `/api/health`, `/api/files`, `/api/plot`, `/api/plot/metadata`,
  `/api/convert/frac`, `/api/structure`, `/api/kde/slice`; plus a data-root guard
  (`RMC_TOOLKITS_DATA_ROOT`).
- Frontend uses `VITE_API_BASE_URL` (default `http://localhost:5000`); file-path field no longer
  fires a request per keystroke.
- Added the Dashboard all-plots view and the KDE/3D page (element filter, z-slice controls, density
  KDE canvas, Three.js folded unit cell with OrbitControls + translucent slab overlay, and the
  slab-in-cell x-z projection).
- Fixed structure sampling so the sample renders all 52,000 atoms and preserves all 52 reference
  sites (sampling grouped by reference number, not raw stride).
- Added a `unittest` suite under `tests/` covering CSV/log parsing, Rwp, `.rmc6f` metadata +
  conversion, structure loading, plot detection/metadata/PNG, and KDE loading/slice.
