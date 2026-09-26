# Development Log

Chronological record of notable changes, newest first. For current architecture and conventions see
[AGENTS.md](../AGENTS.md); for forward plans see [ROADMAP.md](ROADMAP.md).

## v1.0.0 — 2026-09-25

The first stable release. Every engine's math and physics was audited end to end, and the
audit's roughly 100 defects were fixed in both runtimes. They included wrong numbers, Python ↔
browser disagreements, inputs that crashed or silently produced garbage, and docs that described
code that no longer existed. Nine groups fixed them test-first (Auto StoG ×2, structure KDE, PCA
ellipsoid, displacement directions, bond angles, parsers and dashboard, symmetry finder, Flask
API). Each group was reviewed independently, then merged and integrated. Several numbers the app
reports change *because they were wrong*: read **Upgrading from 0.5.0** before comparing 1.0
output with 0.5.0.

### Upgrading from 0.5.0

Each 0.5.0 → 1.0 example below comes from running both releases on the same input. The inputs are
the committed demo run (`web_app/frontend/public/demo`, GaTa₄Se₈ 250 K, 52 000 atoms), the
maintainer's Mn₃Sn and FeCoSn 199 K total-scattering runs (`data/stog_tests`, not in the
repository), or a synthetic cell where stated. The reasons are in
[Correctness fixes](#correctness-fixes), engine by engine.

#### Numbers that change

| Quantity | Where | 0.5.0 → 1.0 | Cause |
| --- | --- | --- | --- |
| Rwp | Dashboard chip; `metrics.rwp` of `/api/plot/*`; `make_plot` | Demo F(Q) 4.569 % → 4.565 %, x-ray G(r) 1.189 % → 1.186 %. A fit with calc = 0.7 × expt: 42.9 % → 30.0 %. The change grows with the amplitude mismatch between the curves. | Divided by the experimental curve, not the calculated one ([parsers](#parsers-and-run-dashboard)) |
| Structure KDE peak height and contour levels | Atomic Density map, both runtimes; `vmax` and `contours` of `/api/kde/slice`; `oriented_kde_slice` | Demo at the default dz 0.08 and bw 0.03 (c, a and Se-only slices): 0 to +23 %, most between +3 and +16 % (c-slice z = 0.55: 115.8 → 134.2). Thick or sparse slabs rise up to +64 % (z = 0.5, dz = 0.2: 49.6 → 80.8; Se at z = 0.3, dz = 0.1, 11 atoms: 220.6 → 362.4). The integrated density is unchanged, and a slab whose atoms all sit on one layer does not move (z = 0.25: 596.4). | The kernel is fitted to the slab's own atoms, without periodic images or the subsample, so it is narrower ([Structure KDE](#structure-kde)) |
| Auto StoG a, b: first shell and fit window | `rmc-autoscale`, `/api/scaling/*`, Auto StoG page; density-limit fits without r₀ | Mn₃Sn (composition given, ρ₀ 0.063049 Å⁻³) at Qmin 1.0: 59438 Q 1.0–27 (a, b) = (0.335, 0.666) → (1.105, −0.102); 300 K 0.331 → 0.752; 500 K 0.272 → 0.936. At Qmin 0.82, 59438 Q 0.82–28: 0.362 → 1.289. | 0.5.0 took the second shell (onset 3.43–3.53 Å) for the first, so its window crossed the first shell ([Auto StoG](#auto-stog)) |
| Auto StoG a, b: Huber weighting | same | Where 0.5.0 already found the first shell, the scale drops: 300 K Q 0.82–28 1.309 → 1.208 (−7.7 %), 500 K 1.567 → 1.397 (−10.8 %), 59438 Q 1.0–29 with `--r0 2.67` 1.199 → 1.011 (−15.6 %). | Rows scaled by √w (the Huber estimator), not by w |
| Auto StoG FZ amplitude and ρ₀ estimate | `--amplitude fz`, `--estimate-rho0` | FeCoSn 199 K (its stog.inp, ⟨b²⟩ = 1.10426): a_fz 1.25116 → 1.1784 (−5.8 %), ρ₀ 0.060692 → 0.057045 Å⁻³ (expert value 0.057329). The density-limit scale barely moves: 1.18293 → 1.18536. | Huber √w |
| Auto StoG with a formula whose ⟨b⟩² disagrees with the stog.inp | CLI `--formula`, API `formula` | FeCoSn stog.inp (⟨b⟩² = 1) with `--formula FeCoSn`: a 1.18874 → 1.18536. The formula's ⟨b²⟩ is no longer used, so `a_fz` and the concordance line (0.5.0: a_fz/a = 0.427, DISCORDANT) are gone unless you pass `--b-sq-avg`. | ⟨b⟩² and ⟨b²⟩ come from one source |
| Automatic low-r enforcement cutoff | the `_rmc` RMCProfile files; the CLI `enforcement:` line | 59438 Q 1.0–27: 3.48 → 2.43 Å; 300 K Q 1.0–27: 3.43 → 2.45 Å; 59438 Q 1.0–29 with `--r0 2.67`: 3.56 Å (above the given r₀) → 2.42 Å. | Cut at the foot of the first shell, never at the onset and never above r₀ |
| PCA `nonGaussianity` | PCA Ellipsoid table; `/api/pca/sites`, `/api/pca/kde` | Demo site 1 (Ta) 4.465 → 4.688, site 2 (Se) 0.456 → 0.398. The per-axis excess kurtosis is unchanged. A symmetric split site now reads negative. | Mean per-axis κ replaced by Mardia's (b₂ − 15)/5 ([PCA ellipsoid](#pca-ellipsoid)) |
| PCA spread of a site at x = ½ of a one-cell-thick box | PCA Ellipsoid, Displacement Directions | Synthetic 8 × 8 × 1 box (c = 10 Å), site at z = ½, σ = 0.08 Å: PC1 rms 4.98 → 0.080 Å. The demo (10 × 10 × 10 cells) is unaffected. | Offsets unwrapped about the site's circular mean |
| Displacement Directions Auto resolution | Resolution "Auto"; `/api/pca/orientation` without `frequency` | Demo sites (1000 copies each): ν = 3 (92 cells) → 2 (42 cells). 300 points: ν 2 → 1; 12 000 points: ν 10 → 9. | `recommended_frequency` floors ([Displacement Directions](#displacement-directions)) |
| Displacement Directions significance | significance strip; payload | Demo site 1 at ν = 10: the strip's "map significance 1.0σ" (the legacy RMS z `significance`) and the peak's "z = 5.4" give way to `peakSignificance` 0.82 σ and `mapSignificance` 0.14 σ. | Calibrated tests replace local z-scores |
| Antipodal-asymmetry null | same | Demo site 1 at ν = 10: `antipodalAsymmetryNull` 0.565 → 0.527 ± 0.018. The asymmetry itself (0.562) is unchanged; z = 1.96, not flagged. | The exact inversion-symmetric null |
| Bond-angle counts, A = C with distinct overlapping windows | Bond Geometry; `/api/triplets`; `rmc-triplets` | Demo Se–Ta–Se, windows 2.3–2.6 / 2.5–2.9 Å: 173 957 → 159 216 angles, mean 109.02° → 109.37°. Equal windows are unchanged (Se–Ta–Se 2.3–2.9 Å: 221 391). | Each physical triplet counts once ([bond angles](#bond-angles)) |
| Bond counts | same | New `uniqueBonds`: demo Ta–Ta–Ta 2.7–3.3 Å keeps `count` 43 624 (B-centred bond vectors, as before) and reports `uniqueBonds` 21 812. | Each physical bond counted once |
| Detected SG ladder | Dashboard Detected SG card | Demo: P1 / Pn / P2 (No. 3, given to a 4-operation set) / P3m / F-43m → P1 / Cm (No. 8) / Cmm2 (No. 35) / P-42₁m (No. 113) / F-43m (No. 216). The full group appears from 0.031 Å (was 0.038 Å). | Closed groups only, standard settings, refined translations ([symmetry finder](#symmetry-finder)) |
| `dispA` | AI Assistant context (`mean_disp_A`, `max_disp_A`) | Isotropic cloud, 0.1 Å per axis: hexagonal cell 0.190 → 0.172 Å, fcc in its 60° rhombohedral cell 0.211 → 0.172 Å. Cubic cells, the demo included, are unchanged. | Mapped through the full cell metric |

#### Files that change

Written files:

- **Classic Auto StoG outputs follow Fortran stog.** `scale.gr` / `<stem>.gr` hold the
  unfiltered g(r), which tends to 1 at large r. `scale_ft.gr` / `<stem>_ft.gr` hold the filtered
  g(r) plus a third column r·[g(r) − 1]. In 0.5.0 they held g − 1 and 4πρ₀r(g − 1). The `_rmc`
  files (Keen G_K, D, F_K) keep their conventions; their values move with a and b (above).
- **The provenance JSON** (`<stem>_provenance.json`) gains `diagnostics.a_fz_rel_se`,
  `a_fz_reliable`, `r_alias_limit` and `rmax_beyond_alias_limit`, `enforcement.source`,
  `provenance.fz_limit` (the Q→0 head fit), `provenance.fit_failure` and
  `provenance.r_alias_limit`, and `rho0_estimate.reason`, `stopped` and `q_first`. When the
  formula's ⟨b⟩² disagrees with the stog.inp value, `diagnostics.a_fz`, `amplitude_concordance`,
  `amplitudes_concordant`, `fk_qmin` and `fk_q0_theory` are absent unless ⟨b²⟩ is given.
- **Output families are written whole or not at all.** `rmc-autoscale`, `/api/scaling/run`,
  `rmc-triplets` and `/api/convert/frac` check every destination before computing and write
  through temporary files renamed into place. A failed run leaves no partial family and never
  replaces an existing file.

Files read differently:

- **stog.inp line 22** follows the classic r₀ rule: `peak_rmin` when the first-peak window starts
  inside the cutoff, else `peak_cutoff`. `2.48 2.65 3.1` now gives r₀ = 2.48 Å (0.5.0: 2.65 Å). A
  value that would leave a fit window narrower than 0.1 Å is ignored, and r₀ is detected.
- **`.rmc6f` atom lines** must be `id element [label]` followed by exactly 7 fields
  (`x y z ref cx cy cz`) or 3 (coordinates only). 0.5.0 misread other layouts silently in both
  runtimes. With one extra trailing field on every demo line, the Flask c-slice at z = 0.25 was
  empty, `/api/pca/sites` found 10 sites instead of 52, and the browser model had 10 basis sites.
  Such lines are now skipped and counted, and a file with none parseable is an error quoting the
  first line. Element names are capitalised in the browser as in Python (SE → Se).
- **A CSV with a stray non-numeric cell** fails with a line-numbered error in the browser. It used
  to plot a silent NaN.

Inputs that now stop before any output is written, with a message naming the problem (HTTP 400, a
thrown worker error or a CLI exit code). 0.5.0 returned a result for most of them, or failed
partway:

- **Auto StoG.** A first shell within ~0.55 Å of `r_cutoff` (Si–O, P–O, B–O and C–O at the
  default 1.0 Å). An unpinned density-limit fit that cannot confirm a first shell: 9 of 56 Mn₃Sn
  (run, Qmin, Qmax) configurations, among them 59438 Q 1.0–28 (0.5.0: a = 0.336) and 55537
  Q 0.82–28 (0.870). A fit with a ≤ 0. Overlapping Q banks, ⟨b²⟩ < ⟨b⟩², a negative or
  non-finite `r_cutoff`, a pinned density-limit window narrower than 0.1 Å. A negative or NaN
  Qmin, a non-finite Qmax, `--scale` NaN or 0, a non-finite `--offset`. An explicit enforcement
  cutoff that is non-finite, negative or at/beyond rmax, and a reversed first-peak window.
- **Auto StoG outputs.** An output that names the input data file or the stog.inp, even with
  `--force` (in `--data` mode, `--out-dir` = the data folder now needs `--out-stem`). Two outputs
  on one file, a directory at an output path, an output folder that cannot be created.
- **`.rmc6f` headers.** A zero, negative or fractional supercell, and a lattice with a
  non-numeric, NaN or missing row, a zero volume or an overflowing volume. `Supercell dimensions:
  0 0 0` used to fold every atom onto the origin.
- **Analysis requests.** A KDE-slice element the file lacks, a zero-spread PCA site's KDE, a
  PCA-KDE volume that captures less than 10⁻⁶ of the density, a negative, NaN or fractional
  orientation smoothing or frequency, a Bond Geometry request above 5×10⁷ angles (app boundaries only; library and CLI
  calls are unrestricted), and a cleared bond-window box.
- **Frac conversion.** A file with no full-layout atom line, and an output that is its own source
  or a directory.
- **`rmc-triplets`.** Destinations that clash, the configuration or a directory as a destination,
  a folder that cannot be created, and a plot format matplotlib cannot write.

#### API changes

| Endpoint | Added | Changed | New statuses |
| --- | --- | --- | --- |
| every `/api/*` route | — | Invalid numbers (NaN, ±∞, text, out of range) are refused by `_number()`; a NaN in a data series is `null`; a non-finite computed result is refused. | 400 where 0.5.0 answered 200 (e.g. `bw=nan` returned `"bw": NaN`, invalid JSON) or 500. 409 when the source file keeps changing during the read. An unknown `/api/*` path is a JSON 404 and a wrong method a JSON 405 with `Allow` (0.5.0: HTML pages). |
| `/api/plot/data`, `/api/plot/metadata` | `chiColumn` (`.log`) | `*_FQn.csv`: title "S(Q) (x-ray)" → "F(Q)", `yLabel` S(Q) → F(Q). `*_PDFpartials.csv`: "PDFpartials", G(r) → "Partial g(r)", g(r). `.log`: "R-value", log(χ), series "R" → "χ² history: X_ray_(R)1", ln(χ²), series "X_ray_(R)1". Classic `scale.gr` / `scale_ft.gr` are labelled g(r). `metrics.rwp` values (above). | An unsupported or unreadable file is a 400 on every plot route (0.5.0 answered some with 500). |
| `/api/structure` | `parseReport`, `parseWarning` | — | 400 for an invalid header or zero parseable atoms. |
| `/api/kde/slice` | `kernel`, `message`, `warnings`; query `ux`…`vz` (the page's in-plane frame for a custom plane) | `element` is case-insensitive (0.5.0: `element=se` drew an empty map). `bw` is `null` for a declined bandwidth. For a custom normal, `z`/`dz` echo the slider fractions: (1 1 1) at z = 0.5, dz = 0.05 echoed 0.866 / 0.0866 in 0.5.0 and echoes 0.5 / 0.05 now; `depth`/`depthThickness` keep depth units. `z` is clamped to [0, 1] and echoed clamped. A custom (1 1 0)-type map in Flask mode is drawn in the page's frame, rotated 90° from 0.5.0. | 400 for an unknown element (0.5.0: 200 with an all-zero map), an overflowing bandwidth and a non-orthogonal frame. |
| `/api/pca/sites` | per site `axisResolved`, `zeroSpread`, `elementCounts`, `mixed`; top-level `parseWarning` | `nonGaussianity` is Mardia's (above). A zero-spread site has `null` axes, anisotropy, κ and `nonGaussianity`. | 400 for an invalid `.rmc6f`. |
| `/api/pca/kde` | `axisResolved`, `boxHalfWidths`, `elementCounts`, `mixed` | `cubicBox` only sizes `boxHalfWidths`; `halfWidths` is always the per-axis box. | 400 for a zero-spread site or a volume that captures less than 10⁻⁶ of the density (0.5.0: a 200 all-zero volume). |
| `/api/pca/orientation` | `peakSignificance`, `peakPValue`, `peakLocalPValue`, `peakCount`, `peakExpected`, `peakTieCount`, `mapSignificance`, `mapPValue`, `mapChiSquare`, `mapDegreesOfFreedom`, `mapNullSd`, `mapNullSkewness`, `mapExpectedPairs`, `antipodalAsymmetryNullSd`, `antipodalAsymmetryZ`, `antipodalAsymmetrySignificant`, `orientationEffectivePoints`, `orientationAnisotropyNull`, `orientationBinghamStatistic`, `orientationBinghamPValue`, `orientationAnisotropySignificance`, `parseWarning` | `mapPValue` and `mapSignificance` are `null` on maps too sparse to test. `recommendedFrequency` and `antipodalAsymmetryNull` values change (above). The legacy `significance` (an RMS z) and `peakZScore` (a local z) stay but are no longer displayed. | 400 for negative, NaN or non-integer `smoothing` and a non-integer `frequency` (0.5.0: `smoothing=-1` was a silent no-op). |
| `/api/triplets` | `lengths12/23.uniqueBonds`, `parseWarning` | A = C counts with distinct windows (above). | 400 above 5×10⁷ angles, counted before any angle is formed. |
| `/api/scaling/preview`, `/api/scaling/run` | `warnings` (the coefficient warnings the CLI prints); in `diagnostics` `a_fz_rel_se`, `a_fz_reliable`, `r_alias_limit`, `rmax_beyond_alias_limit`, `first_shell_below_r0`; in `provenance` `fit_failure`, `fz_limit`, `r_alias_limit` | `enforce` is one tri-state. The concordance fields can be absent (see Files). JSON booleans are strict: `true`/`false`, `0`/`1` or the words 1/true/yes/on and 0/false/no/off. | 400 for every refusal listed above, a malformed `.inp` (with its own parse error), a non-boolean flag and a body nested too deeply. 404 when a stog.inp names a folder as its data file. 409 when an output exists and `force` is false (unchanged). |
| `/api/convert/frac` | `parseWarning` | `overwrite: "false"` is false (it counted as true). | 400 for the refusals listed above. |

#### Defaults and safety

- **`python web_app/backend/app.py` listens on `127.0.0.1` with debug off.** 0.5.0 listened on
  every interface (`0.0.0.0`) with Flask debug mode on, which exposed the Werkzeug interactive
  debugger (arbitrary code execution) and the data API to the network. Set
  `RMC_TOOLKITS_HOST=0.0.0.0` to listen on the network, or use Gunicorn or Docker. Set
  `RMC_TOOLKITS_DEBUG=1` for local debugging only; the startup line warns when debug meets a
  network address. `PORT` still wins over `RMC_TOOLKITS_PORT`, and a malformed value stops the
  server with an error naming the variable (`server_settings()`).
- **The Docker image no longer bundles `data/`.** The `Dockerfile` used to `COPY data ./data`,
  which failed on a clean clone (the folder is gitignored) and baked private runs into the image
  where it existed. `.dockerignore` now excludes `data/`, `node_modules`, virtualenvs and caches,
  and the image has an empty `/app/data` to mount runs on
  (`docker run -v /path/to/runs:/app/data …`).
- **PCA isosurface level: 25 % → 50 %**, the ellipsoid's probability.
- **Automatic low-r enforcement** cuts at the first-shell foot (numbers above). It was already
  automatic in 0.5.0 whenever no cutoff was given (CLI `--data`, the API without a stog.inp, a
  blank Cutoff on the page), although 0.5.0's `--help` said it was off in `--data` mode. A
  stog.inp or `--enforce-cutoff` still wins, and `--no-enforce` turns it off.
- **Displacement Directions** keeps ν = 10 as the default; only Auto moved (above).
- **`--estimate-rho0`** seeds ρ₀ = 0.05 Å⁻³ when the input gives no density.
- **Live Data in Flask mode** reloads the analysis pages in place when the `.rmc6f` changes. With
  Live Data off, press Load to pick up a newly saved configuration.
- **The static-mode workers refuse what the Flask routes refuse** (see
  [Robustness and parity](#robustness-and-parity)), so a malformed request errors in both runtimes.

### Highlights

- **Physics fixes with visible consequences.** The dashboard Rwp was divided by the calculated
  curve, reading 43 % high for calc = 0.7 × expt. PCA sites at x = ½ in a one-cell-thick box were
  torn in two, giving an rms of half the box edge (5.0 Å in a 10 Å box) instead of 0.08 Å. Auto
  StoG fitted its density-limit window across the first shell of most oxides (synthetic SrTiO₃:
  a = 5.28 for a true 10; ReO₃: a = −3.2), and its automatic low-r cleanup cut 6–9 % of the first
  shell away. The Detected SG ladder named operation sets that are not groups ("P2 No. 3" for a
  4-operation set on the bundled demo).
- **One answer in both runtimes.** The browser ports are held to committed Python goldens: Auto
  StoG to round-off (≤ 5×10⁻¹³ on every real run), the Structure KDE to 10⁻⁶ of the peak
  (measured ≤ 2×10⁻¹²) on slabs within its 6000-row fit cap, bond angles to exact integer
  histograms, and displacement directions to 10⁻⁹ relative. The cap counts the periodic-image
  rows too, so a few thousand atoms reach it at larger bandwidths. Above it the two runtimes sum
  different unbiased subsamples, and their maps differ by several percent of the peak (5.6–8.5 %
  on over-cap slabs of the bundled 52 000-atom demo run); the Flask/SciPy map is the reference.
  The shared `.rmc6f` parser and the plot readers have goldens too. PCA is the exception: each
  engine is pinned to its own reference, with no shared golden ([ALGORITHMS.md](ALGORITHMS.md)).
- **Auto StoG files match classic stog.** The classic-named files hold the Fortran stog
  functions, g(r) and r·[g(r) − 1], and the stog.inp r₀ rule is the classic one. Automatic
  low-r enforcement cuts at the foot of the first shell, and a fit with a non-physical scale is
  never written.
- **Statistics you can quote.** Displacement Directions reports calibrated significances instead
  of local z-scores. PCA non-Gaussianity is rotation-invariant (Mardia). Bond Geometry reports
  physical bond counts and counts each same-element triplet once.
- **Inputs are validated, not guessed.** Both runtimes share one validated `.rmc6f` atom-line
  grammar, with a parse report against the header's atom count; this fixes the zero-atom
  static-mode issue. Every Flask numeric parameter is range-checked (bad input is a 400, never a
  200 with invalid JSON), and parsed-file caches never serve a torn read.
- **Stable-release tooling.** Dependency floors are declared, and CI runs a Python
  3.9/3.11/3.13 matrix. Tests on the committed demo run now exercise real RMCProfile files in CI,
  where the real-data tests skip.

### Correctness fixes

#### Auto StoG

*`scaling.py`, `scaling_cli.py`, `autoScale.js`, `rmc-autoscale`*

- **First-shell detection finds the first shell.** `detect_first_peak_onset` took the flank of
  the strongest |g| feature, which skipped weak or inverted first shells: SrTiO₃'s Ti–O read at
  2.6 Å instead of 1.85 Å, and Mn₃Sn 59438 (Q 1.0–28) at 3.49 Å instead of ~2.74 Å. It now returns
  the smallest-r shell of either sign that stands out of the ripple field below it. A flank that
  reaches the search start is not a shell.
- **The density-limit window is placed from the data.** Without r₀ or `r_fit_max`, the blind
  [1.2, 2.2] Å window used to sit across the first shell. Two trial fits now propose shell onsets.
  The smallest onset is refitted on [lo, onset − 0.25] Å and accepted only when that refit
  re-detects it within 0.15 Å, so a ripple lobe is dropped (at most 4 refits). A refit with
  a ≤ 0 raises. Every fit returned on the four Mn₃Sn runs (Qmin 0.82/1.0 × Qmax 24–30) has a > 0
  and a window top of 2.40–2.50 Å.
- **Automatic low-r enforcement at the first-shell foot.** Enforcement used to zero g(r) up to
  the first-shell onset, deleting 6–9 % of the first-shell pair density. It now cuts at
  min(foot, anchor − 0.25 Å), where the anchor is min(detected onset, given r₀) and the foot is
  `first_shell_foot`. The cut never lands above a given r₀, and the first-shell coordination
  number is preserved to ≤ 0.3 %. On Mn₃Sn 59438 at Q 1.0–27 the cut is at 2.43 Å (expert
  `rmccut` 2.48 Å). A given r₀ that sits above the detected shell is flagged
  `first_shell_below_r0`.
- **Huber IRLS is the Huber estimator.** Both engines multiplied rows by the weight w, which
  minimises Σw²e² (a redescending estimator). Rows are now scaled by √w. On FeCoSn 199 K the
  density and Faber-Ziman amplitudes agree to 0.6 % (5.8 % before), and the ρ₀ self-consistency
  estimate lands 0.5 % from the expert density (5.9 % before). The √w-only table and the
  0.5.0 → 1.0 table are in [auto-stog.md](algorithms/auto-stog.md) (Step 6).
- **The page fits what the CLI fits.** The JS loop filtered with S(0) = 0 instead of the
  composition target, so page and CLI (a, b) differed by up to 2.0 % on Mn₃Sn. They now agree to
  round-off, and so does the iterated ρ₀.
- **ρ₀ self-consistency.** `estimate_rho0` confines ρ₀ to 0.005–0.25 Å⁻³ and accepts a root only
  where the density limit holds. With r₀ pinned at 2.6 Å, Mn₃Sn 300 K had "converged" at
  0.428 Å⁻³, about 50 g/cm³. Every
  non-converged exit carries a `reason`, and the CLI and page name the trial density that could
  not be fitted. `extrapolated` is judged on the first measured Q.
- **⟨b⟩² and ⟨b²⟩ come from one source.** A composition's Sears ⟨b²⟩ was paired with an
  unrelated ⟨b⟩² (for example 1 for normalised x-ray data). That invented an S(0) target (+0.55
  for FeCoSn) and a 60 %-low "converged" density. S(0) > 0 is now rejected.
- **Q order and bank overlap.** Descending-Q files are sorted. Before, they gave a negated G(r),
  and a reversed synthetic fitted a = −10.1 for a true +9.97. Duplicate or overlapping Q
  (concatenated banks) raise.
- **Despiking and σ.** `--despike` runs once, so the fit and the written files share one point
  set (2136 vs 2308 points on 59438 before) and `n_despiked` is the true count. The page applies
  the CLI's σ-column guard, and an unconverged JS run reports `max_iter` iterations.
- **Faber-Ziman conditioning.** New fields `a_fz_rel_se` (Huber sandwich error) and
  `a_fz_reliable`. A reliable a_fz is described as necessary, not sufficient, and 59438's a_fz is
  flagged unreliable.
- **Aliasing limit.** `r_alias_limit` is π/max ΔQ of the transformed grid, and an r_max beyond it
  is flagged. Despike gaps lower it to 24 Å on 59438.
- **No sliver fit windows.** A stog.inp line-22 cutoff pins r₀ only when the default window it
  leaves is at least 0.1 Å wide (`MIN_AUTO_WINDOW`, the automatic placement's floor); otherwise
  r₀ is detected. A cutoff just above r_cutoff + 0.45 Å used to pin a 1–6 point window: on
  FeCoSn 199 K, `1.46 0 0` gave a = 1.001 and `1.0 0 0` with r_cutoff 0.5 gave a = 0.672 (1.185
  with the shipped file), both "converged" with the density limit "satisfied". A pinned
  density-limit window (r₀ / r_fit_max) narrower than 0.1 Å now raises. `r_cutoff` must be finite
  and ≥ 0 in both engines (a negative value returned a "converged" a = 0.026 and was written).
- **The input is never an output.** In `--data` mode the default stem is the data file's, so
  `--out-dir` pointing at the data folder made `<stem>.sq` the measured file itself, and
  `--force` (which the refusal suggested) replaced it with the scaled, cropped S(Q); a rerun then
  silently fitted already-scaled data (a = 1).
- **Outputs are written whole or not at all.** A stog.inp declaring the FK(Q) name as `ft.dat`
  used to exit 0 with the RMCProfile input replaced by the filter correction, and a declared
  `sub/rmc.gr` or a directory at `ft.dat` failed after five or seven files were written. Every
  destination is now checked before any computation (names compared case-insensitively), the
  family goes through temporary files, and declared subfolders are created. An API `outDir` that
  is a file is a 400, not a 500 after the fit.
- **Input checks in every entry point.** `--scale nan` exited 0 with nine all-NaN RMCProfile
  files. A NaN Qmax passed Python's `qmax <= qmin`. `--enforce-cutoff nan` was reported as
  applied, and a cutoff of 1000 with rmax 50 replaced the whole G(r). `ScalingConfig` /
  `makeConfig` and `validate_enforcement` / `validateEnforcement` (CLI, API, page worker) now
  refuse these. A given r₀ with no detected shell no longer crashes with a TypeError.
- **stog readers.** Both engines read CR-only line endings, a BOM, bad bytes, Fortran D exponents
  and ASCII-only digits the same way.

#### Transforms

*`transforms.py` and its JS port*

- The Lorch low-Q correction basis is cancellation-free. Grid points 10⁻⁹–10⁻⁶ Å from π/Qmax
  got coefficients off by up to ~80× or with flipped sign. The unwindowed moments use a Taylor
  series at small Q₀r, where they were off by 100 %.
- The Fourier filter accepts an r grid that starts at 0. Before, every output Q was NaN;
  `g_filtered(0)` is now the continuous extension (`gpdf_slope_at_zero`).
- The public transforms raise on a grid that is not strictly increasing.

#### Structure KDE

*`kde.py`, `localKdeWorker.js`, `gpuKde.js`*

- **The browser draws the SciPy kernel.** The worker's 10⁻⁸ ridge and its isotropic fallback are
  gone. Both runtimes evaluate H = bw²·C exactly and decline the same slabs with the same
  `message`. Before, on the GaNb₄Se₈ Ga layer at bw ≤ 0.015, the browser map differed from SciPy's
  by 36–55 % of the peak.
- **Kernel from the source atoms.** C is the covariance of one row per slab atom. Periodic images
  and the 6000-row subsample enter only the density sum, so the thickness and bandwidth sliders
  no longer reshape the kernel. The kernel is narrower, so peaks rise in both runtimes (sizes in
  [Upgrading](#numbers-that-change)).
- **One slab test**, `|d − z_c| ≤ dz/2 + 10⁻⁹`, in `kde.py`, the worker and the Slab-In-Cell
  highlight.
- **Contours.** Flask's log mode no longer drops contours when the peak is below 1. A map whose
  kernel misses every grid node is flagged `unresolved` and is neither contoured nor painted.
- **SciPy.** `_FixedCovarianceKDE` drives SciPy 1.8 through 1.18; `/api/kde/slice` used to answer
  500 on SciPy < 1.10.
- **Honest labels.** The custom slice is labelled as the Miller plane (h k l) it is. The map
  prints the slab thickness and the kernel σ in Å and warns on sub-grid kernels.
- **Custom (hkl) slices are drawn in one frame.** For a custom plane the page now sends its
  in-plane frame to `/api/kde/slice` (`ux`…`vz`), and the route draws the map in it
  (`_custom_slice_frame`, validated orthogonal to the normal). Flask used to pick its own axes, so
  the default (1 1 0) map — and (1 0 1) — came back rotated 90° against the browser map and against
  the page's own Slab In Cell panel, letterboxed in a panel sized for the other orientation. The
  densities were always the same; the parity golden now includes (1 1 0) and (1 0 1) demo slices.

#### PCA ellipsoid

*`pca_kde.py`, `pcaKde.js`, `pcaCrystalFrame.js`*

- Offsets are unwrapped about each site's circular mean, not folded about zero. A site at x = ½
  in a one-cell-thick box, or near x = 1 in a two-cell box, was torn in two, making U 10³–10⁴×
  too large. Displacement Directions inherited the fix.
- The browser ellipsoid scale uses the exact χ²₃ quantile. The old table was up to +3.5 % off and
  had a wrong 0.6827 node.
- The volume, mass, iso levels and PC walls are always sampled on the per-axis box. On planar
  clouds the captured mass had read 0 % or ~2×10⁵ %. The shell no longer paints clamped box-face
  density, and the crystal-frame walls are line-integral marginals (no moiré; error 6–22 % → 2–3 %
  at grid 40).
- Mixed-occupancy sites are labelled by their majority species. Element pooling selects atoms by
  their own element.
- Both engines refuse a volume that captures less than 10⁻⁶ of the density (a bandwidth far below
  the node spacing, or an extent of 10⁶); it used to be an all-zero volume with no warning.
- The Cartesian → CIF U_ij formula in the docs is corrected, as is the claim that the
  20 000-point cap binds only pooled clouds (it also binds single sites in boxes of ≥ 28 cells per
  edge).

#### Displacement Directions

*`orientation.py`, `orientation.js`*

- **Calibrated significances.** The peak readout is `peakSignificance`: the exact Poisson tail
  of the peak cell, Šidák-corrected over all cells. The map test compares Pearson's X² with a
  gamma matched to its exact isotropic moments, in `mapSignificance`. It replaces an RMS z
  printed as "σ": a hemisphere-only cloud used to read "1.4 σ" and now reads > 10 σ. The map test
  is withheld below 0.1 expected coincident pairs.
- **Antipodal asymmetry.** The null is now the exact inversion-symmetric one, conditional on
  each pair's total, and carries an engine flag at 3 SD. The old √(C/πN) floor could exceed 1,
  so the red flag could never fire at ν = 10.
- **Anisotropy.** The isotropic expectation is 9/√(10πN_eff), and a Bingham test is added.
- **Ties.** Exact Voronoi ties are broken centrosymmetrically (an exactly centrosymmetric cloud
  read an antipodal asymmetry of 1.000 at odd ν),
  tied peaks resolve by one 10⁻⁹ rule, and neighbour cycles start at their smallest index.
  Before, 185 of 1002 neighbour rows differed between engines at ν = 10.
- `recommended_frequency` returns the largest ν with ≥ 12 points per cell, as its docstring
  promised.
- Non-finite rows and options are rejected by name. Before, Python failed with a LAPACK error
  and JS returned NaN axes.

#### Bond angles

*`triplets.py`, `triplets.js`, `rmc-triplets`*

- **Bounded work.** Angles stream into the histogram, and both app boundaries refuse a request
  whose exact angle count exceeds 5×10⁷ before forming any. The 15 Å cap had allowed ~10⁹-angle
  requests: tens of GB in Flask, `RangeError` in the worker.
- **One count per physical triplet.** For A = C with distinct windows the old ordered rule
  counted overlap triplets twice.
- **Ideal geometries are deterministic.** An angle within 10⁻⁹° of a bin edge bins on it, and a
  bond within 10⁻⁹ Å of a window bound is inside. A shell typed at its exact distance used to
  lose up to 60 % of its bonds.
- **Physical bond counts.** `uniqueBonds` counts each bond once; the directed count is 2× when
  A = B.
- A cleared window box is an error. It used to become rmin = 0 silently: Nb–Nb–Nb 3.5–4.6 Å
  became 0–4.6 Å.
- `rmc-triplets <run folder>` picks the same `.rmc6f` as the app.
- The docs correct the sin θ rationale and give the factor to RMCProfile's `norm/sin(theta)`.
  Angle totals match RMCProfile's own TRIPLETS output for five triplets on the 5 K run.
- **`rmc-triplets` checks its destinations before computing.** `--dump-angles` equal to
  `--output` replaced the histogram with exit 0, and an unsupported `--plot` format printed a
  traceback after the CSV was written. `--version` added.

#### Parsers and run dashboard

*`parsers.py`, `plots.py`, `browserData.js`, `rmc6f.js`*

- **Rwp** is normalised by the experimental column. The column roles come from the header when
  it names them, else from RMCProfile's (x, calc, expt) order.
- **`.rmc6f` atom lines.** Both runtimes read the fields from the end of the line, so one extra
  trailing field shifted every coordinate by one column, silently (examples under
  [Files that change](#files-that-change)). Both now share one anchored, validated grammar. It
  reads D exponents, any spelling of the `Atoms` marker, bare-CR files and coords-only lines. It
  skips and counts unparseable and non-finite lines and reports "k of n atom lines unparsed".
- **`.rmc6f` header.** `read_cell_vectors` and `readRmc6fCellVectors` read the header numbers
  like the atom lines (D exponents) and validate the supercell and the lattice, with the same
  error text in both runtimes. A NaN, collinear or 10³⁰⁰ lattice used to give 200s with NaN
  cells, an anisotropy of 10¹⁴ or all-zero angle counts, and a D-exponent lattice was a raw float
  error.
- **Structure files.** 0-byte or marker-less candidates are skipped; an empty `.rmc6f` used to
  hide six valid configurations in `data/250K_try1/supercell`. `read_structure` pairs Frac/rmc6f
  files by stem. Before, it could fold a 5×10×10 Frac file with a 10×10×10 supercell.
- **χ² history.** `.log` reads check the header's column count, drop a half-written last line and
  keep NaN rows as gaps. Static mode charts one run's logs, not every run spliced together. The
  curve is named by its column (`χ² history: X_ray_(R)1`): it is one fit term, not a total.
- **Labels.** `*_FQn.csv` is labelled F(Q), `*_PDFpartials.csv` g(r); PNG and Flask use the same
  strings.
- **dispA** goes through the full cell metric. It was +10 % in hexagonal and +23 % in
  rhombohedral cells.

#### Symmetry finder

*`symmetry.js`, `spaceGroupSymbol.js`, `wyckoff.js`*

- **Closed groups only.** A space-group number is given only to an operation set that is closed
  under composition. A set that fails is labelled "not a group" with no number (1 of 2065 rungs in
  random sweeps). Each operation's translation is least-squares refined, so the answer no longer
  depends on atom order; the demo's full-group residual fell from 0.038 to 0.031 Å.
- **Standard settings.** A symbol is shown only when it is a tabulated symbol of the detected
  class and centring, in a standard setting found from the group's own elements (any axis order,
  centred, primitive or rhombohedral cells, supercells reduced). Otherwise the card shows the
  crystal class or a "≥" lower bound, with no number. A centred group given on a primitive cell
  is named in its centred cell (C2 used to read P2 No. 3, I4 P4 No. 75, R3 P3 No. 143).
- **Lattice rotations** are tested by Cartesian strain on the τ scale. A c/a = 1.004 perovskite
  reads `P4/mmm` below 0.016 Å and `Pm-3m` above it.
- **I2₁2₁2₁ and I2₁3 are told apart** from I222 and I23. The 0.5.0 entry listed that as a known
  limitation.
- **Wyckoff letters** are read in the naming cell and paired with that cell's multiplicity. The
  assistant's labels follow.
- **Correction to the 0.5.0 entry.** The table holds 1731 positions, not 1724. The seven
  P222₁ (#17) and I2₁2₁2₁ (#24) special positions the 0.5.0 check "rejected" are correct. The
  rejection came from test fixtures that described those two groups about a shifted origin. At
  the ITA origin, P222₁'s 2a site `x,0,0` lies on the 2-fold axis along a and expands to 2
  points. All seven are restored.
- **Bounded cost.** Bases over 2000 sites, or boxes over 384 candidate operations, read
  "not analysed" instead of freezing the page. An unanalysable lattice reads "undetermined"
  instead of P1 No. 1.

#### Flask API

*`web_app/backend/app.py`*

- **One validator for every numeric parameter** (`_number()`). NaN, ±∞, text, lists and
  out-of-range values are a 400 that names the parameter. Before, they could produce a 200 with
  invalid JSON or an all-zero map. A computed result that comes out NaN/∞ is also a 400,
  including a KDE bandwidth whose kernel overflows and a scaling fit that overflows.
- **Parsed-file caches** are keyed on the full file signature (`st_mtime_ns`, `st_ctime_ns`,
  `st_size`, `st_ino`), never on a whole-second mtime. A file rewritten over sshfs or `scp -p`
  used to be served stale until the server restarted. A parse of a file that changed while it
  was being read is never cached; a file that keeps changing is a 409.
- **Live Data in Flask mode** re-checks the `.rmc6f` signature on every poll and on Load. The
  analysis pages then reload in place, keeping the picks that still apply, so a page never mixes
  two configurations. Only Bond Geometry's computed distribution is dropped. Every Three.js view
  releases its WebGL context on teardown.
- **Request edges are 4xx, never 500 or a silent 200.** The static-asset rule matched `/api/*`
  before the SPA route, so unknown paths and wrong methods got Flask's HTML pages. A malformed
  `.inp` with the default `kind: 'auto'` was re-read as S(Q) data (preview: "data mode requires
  qmin and qmax"). JSON booleans were read loosely (`"maybe"` read as false, an object as
  true, `inspect: "false"` entered inspect mode). A stog.inp whose data file is a folder and a
  JSON body nested thousands deep were 500s. `/api/kde/slice` matched `element` case-sensitively
  and drew an all-zero map for an element the file lacks. `/api/convert/frac` wrote a
  header-only Frac file with a 200 for a file with no full-layout atom line, and answered 500
  for an output that is its own source or a directory. The new statuses are in
  [API changes](#api-changes).
- `/api/scaling/*` parse `enforce` once, as a tri-state. They return the coefficient warnings
  the CLI prints and never mutate a cached result.
- [REFERENCE.md](REFERENCE.md) documents all 15 routes with their ranges and caps, and a
  contract test keeps it complete.

### Robustness and parity

- Python ↔ browser goldens: Auto StoG (`autoscale_fixture.json`, round-off), Structure KDE
  (`kde_parity_fixture.json`), bond angles (`triplets_fixture.json`), displacement directions
  (values shared between `test_orientation_fixes.py` and `orientationFixes.test.js`), the
  `.rmc6f` parser (`rmc6f_coords_only_fixture.json`) and the plot readers
  (`plot_parity_fixture.json`, including neutron, Bragg and EXAFS layouts). The WebGPU shader is
  checked by a float32 emulation.
- Non-finite input is handled one way everywhere: `.rmc6f` lines with NaN/∞ coordinates are
  skipped, counted and reported in both runtimes (no batched `LinAlgError`). NaN in a data series
  is `null`. A NaN result is a 400.
- Both runtimes resolve a flat run folder to the same configuration (`find_run_configuration`,
  `chooseStructureFile`), with code-point tie-breaks. A picked folder with subfolders can differ:
  the browser also searches the subfolders (and falls back to the first usable `.rmc6f` by full
  path), while the server looks only at the folder itself.
- The package root exports the 1.0 engine API (`first_shell_foot`, `auto_enforcement_cutoff`,
  `fz_limit_fit`, `alias_limit`, `bond_angle_summary_from_file`, …).
- The static-mode workers refuse what the Flask routes refuse (`workers/requestGuards.js`): the PCA
  worker's orientation `frequency` must be an integer and `smoothing` an integer in [0, 64] (10⁹
  passes used to pin the worker), a KDE or orientation result holding NaN/∞ is an error with
  `_strict_result_response`'s message instead of a posted NaN volume, and an unknown `kind` is an
  error instead of a silent KDE. The Auto StoG worker requires a finite, non-zero `a` and a finite
  `b` in manual mode (a NaN `a` used to post `ok: true` with all-NaN curves) and refuses a
  non-finite fit. The structure worker reads `maxPoints` as `/api/structure` does. Every worker
  answers a null message with an error, so the caller's promise always settles.
- Coordinates-only site reconstruction on the PCA page is 3.5–6.5× faster. The `.rmc6f`
  classifier's fast path parses 52 000 atoms in ~0.2 s.

### Tooling

- Dependency floors `numpy>=1.22`, `scipy>=1.8`, `matplotlib>=3.6` and `contourpy>=1.0.7`
  (`kde.py` imports contourpy directly). The suite passes on the floors (Python 3.9) and on the
  newest releases (Python 3.13).
- CI runs the Python suite on 3.9 with the floors pinned exactly, and on 3.11 and 3.13 with the
  latest releases. The frontend job runs lint, `npm test` and the build.
- The committed demo run (`web_app/frontend/public/demo/GTS_250K.*`) backs tests that run in CI:
  `tests/test_parsers_demo_run.py`, `__tests__/demoRun.test.js`, the KDE golden and the plot
  parity golden. The GaNb₄Se₈ and `stog_tests` real-data tests still skip without `data/`, and
  the full Mn₃Sn sweep is opt-in (`RMC_TOOLKITS_FULL_SWEEP=1`).
- Version 1.0.0 (Production/Stable classifier). The Auto StoG tab stays hidden in the shipped
  build (`SHOW_AUTO_STOG = false`); the engine, CLI and API are supported.

### Deferred beyond 1.0

The maintainer decisions the audit considered and did not take for 1.0 (a physical KDE kernel,
plotting every χ² column, a Python symmetry finder, Auto as the default orientation resolution,
AXIS_RESOLUTION_SIGMAS = 2, automatic FZ fallback for degenerate density limits, bounded A = C
pairing work, and more) are listed per engine in [ROADMAP.md](ROADMAP.md#1x-candidates).

## v0.5.0 — 2026-08-14

Two new analysis pages, correct space-group naming, and the math reference.

- **Bond Geometry** — bond-angle (triplet) distributions in the app, in both runtimes, on the new
  `rmc_toolkits.triplets` engine and its `rmc-triplets` CLI.
- **Displacement Directions** — the hex-binned orientation sphere.
- **Auto StoG** — composition-first absolute scaling, ρ₀ self-consistency and the `rmc-autoscale`
  CLI. The engine and CLI ship; the tab stays hidden until the page behaves.
- **Symmetry** — non-symmorphic groups named as themselves against all 230, with Wyckoff letters
  per orbit, plus the first tests for the symmetry code.
- **[docs/ALGORITHMS.md](ALGORITHMS.md)** — a code-anchored account of every operation each page
  performs, with seven per-page derivations.
- **PCA Ellipsoid** — principal axes reported in the crystallographic frame.
- Fixes: Rwp reports "unavailable" rather than a fake perfect fit; the ρ₀ fixture tolerance is
  per dataset.

**"Local Geometry" renamed to "Bond Geometry"** (2026-08-14) — the tab shipped in this release
under the name *Local Geometry*, which is generic: in the PDF/total-scattering sense every
analysis page here is local structure, so the name didn't say what distinguishes this one. The
page is bond angles, bond lengths, coordination and detected bonds — bond geometry — and the tab,
component (`BondGeometryPage.jsx`/`.css`), algorithm reference
([algorithms/bond-geometry.md](algorithms/bond-geometry.md)) and every cross-reference now carry
that name. The internal page key stays `geometry`, and the engine keeps its RMCProfile-rooted
`triplets` naming (`rmc_toolkits.triplets`, `rmc-triplets`, `/api/triplets`,
`workers/triplets.js`) — the rename is presentation-layer only. Entries below predating the
rename keep their original titles; the tag was re-cut so v0.5.0 ships the tab under its final
name.

**Local Geometry: folded-cell canvas overflowed the page on retina displays** (2026-08-14) —
`FoldedCellPanel` resized its renderer with `setSize(width, height, false)`, which skips styling
the canvas, and `.pca-structure` gives the canvas no CSS size of its own — so the canvas laid out
at its attribute size, CSS size × devicePixelRatio, twice its panel on a 2× display, and the
oversized canvas wrecked the whole tab's layout. Plain `setSize(width, height)` now lets three.js
style the canvas, as `SiteStructurePanel` always has. Also fixed in passing: the panel's *Reset
view* called `controls.reset()` with no saved state, snapping the camera to its construction
position at the origin (inside the cell); the rebuild now ends with `controls.saveState()` so
reset restores the framed view of the current data. The v0.5.0 tag was re-cut onto this fix's
merge commit; the original 2026-08-13 tag, which pointed at the release-notes commit and so
predated both this fix and the algorithm reference, was deleted before anything referenced it.

**Local Geometry algorithm reference** (2026-08-13) — new
[algorithms/bond-geometry.md](algorithms/bond-geometry.md) (created as `local-geometry.md`,
renamed with the tab), completing the per-page set: every
analysis page now has its own code-anchored derivation document. Engine section: inputs and the
half-up bin-count rule, folding + image bookkeeping, the linked-cell search with the
strictly-more-than-$(q-1)$-thicknesses covering argument, bond admission (inclusive windows,
zero-length exclusion, self-image rule), the two pairing rules, the three normalizations (with
why the bin-integral `sinth` stays finite at 0°/180°), the summary payload, the app-boundary
caps, and the CLI. Page section: seeding, the epoch guard, the `fit` plot variant, the
partial-g(r) helper rules (curves by bond type, guides by window role, palette-matched colors),
and the folded-cell bond view — including the caveat that its sticks are average-structure
bonds, not instantaneous ones. Parity is quoted from the golden-fixture tests (exact integer
counts; 1e-9 / 1e-7 / 1e-5 float tiers) rather than asserted. `ALGORITHMS.md` index row and the
README's bond-angle derivation link now point at it.

**Auto StoG: `test_estimate_rho0_near_hand_value` failed on the 100 K dataset** (2026-08-13) —
the x-ray fixture tests select the first available FeCoSn run (`100K` or `199K`), on the
stated grounds that the two share a stog parameterization. True for every assertion in
that class except one: how closely `estimate_rho0` reproduces the expert's hand density is
a property of the measured S(Q), not of the parameterization. The 4.7 % figure recorded
here and in the plan was measured on **199 K**; on **100 K** the estimate lands at 0.0640 Å⁻³,
11.7 % from the hand 0.057329 — over the hardcoded 10 % bound. Anyone holding only the
100 K run saw a red suite, and CI never caught it because both datasets are gitignored and
the test skips.

- Not an algorithm fault: on 100 K the estimate converges to 0.064033 ± 0.00002 from seeds
  spanning 0.02–0.20, concordance → 1 to 10⁻⁵ in 2–3 iterations, `extrapolated=False`. It
  is the correct root of `a_fz/a_density(ρ₀) = 1` for that data. `scaling.py` is unchanged.
- The tolerance is now per run (199 K 10 %, 100 K 13 %) rather than one bound covering
  both — a shared ~13 % bound would stop the test noticing a regression on 199 K.
- New `test_estimate_rho0_is_seed_independent` asserts the part that does *not* depend on
  which temperature is on disk: seeds spanning a 10× range must agree to <1 %. A
  seed-dependent answer would mean the fixed-point iteration is terminating on its seed
  rather than on the data, which no tolerance on the hand value would catch.
- Both figures are now documented (`docs/algorithms/auto-stog.md` Step 10, plan §ρ₀); the
  100 K entry previously recorded scale accuracy but never ρ₀, since `estimate_rho0`
  postdates it.

**Non-symmorphic space groups are named correctly; Wyckoff letters for all 230** (2026-08-06) —
the *Detected SG* panel reported the **symmorphic parent** of every group with a screw axis or
glide plane. The operations were always right; only the naming step threw the translation parts
away, so diamond came back as Fm-3m rather than Fd-3m, a Pnma structure as Pmmm, an a⁰a⁰c⁻ tilted
perovskite as I4/mmm rather than I4/mcm, and hcp as P6/mmm. Since octahedral-tilt phases are almost
all non-symmorphic (I4/mcm, Pnma, R-3c, Imma), this hit the common RMC case.

- New `spaceGroupSymbol.js`: splits each operation into its **intrinsic translation**,
  `(1/n)·Σ Rᵏ·t`, the projection of `t` onto the invariant subspace of `R`. It is
  origin-independent, which is what makes it a reliable element label — non-zero along the axis
  means a screw `n_m`, non-zero in the mirror plane means a glide `a/b/c/n/d/e`. The symbol is
  assembled positionally over each crystal system's symmetry directions. Handedness comes from
  `det[axis, v, R·v]`, so the enantiomorphic pairs (4₁/4₃, 3₁/3₂, 6₁/6₅, 6₂/6₄) stay distinct.
- New `spaceGroupTable.js`: all 230 groups, replacing a 64-number symmorphic-only table. Pre-2002
  spellings (Cmca → Cmce) resolve to the same number.
- The builder **proposes ranked candidates and the table disposes**, rather than committing to one
  spelling. Necessary because H–M plane priority (m > e > a > b > c > n > d) is not always the
  answer: I-centring puts both b- and c-glides in the same plane of I4/mcm, where ITA writes `c`.
  Non-standard axis settings are searched the same way — a structure handed over in the Pbnm
  setting is reported as Pnma.
- **R centering is detected** (obverse *and* reverse). `CENTERING_SETS` had only F/I/A/B/C, so a
  rhombohedral structure reported `P`; R-3m came out as P-3m1 — wrong centering, and order 12
  claimed for a group with 36 operations. The `R3`/`R-3m` rows of the old number table and the
  `trig: 'PR'` entry in `ALLOWED_CENTERING` had been unreachable.
- A symbol is only built from elements the positions can place. A subgroup part-way up the
  tolerance ladder keeps its cubic parent's axes, putting its 3-fold along ⟨111⟩ where the
  hexagonal-axes positions expect [001]; that used to emit the meaningless `P2m`, and now falls
  back to the crystal class (`P3m`). A symbol merely in a non-standard setting (`Pn` for #7) is
  still shown as-is. Naming also runs only for operation sets that close, so the ladder no longer
  pays for it on partial mid-transition sets.
- **Wyckoff letters for all 230 groups**, up from four (216/225/229/221). New `wyckoff.js` +
  `wyckoffTable.js` (1724 positions, packed one string per group). Multiplicity and site symmetry
  alone are often ambiguous — Pm-3m has both 3c and 3d at multiplicity 3 with site symmetry 4/mmm —
  so ties are broken on the position's coordinate form, tested against every member of the orbit
  because the detected representative is whichever atom came first, not the one International
  Tables prints. That fixes cubic SrTiO₃, whose oxygen orbit previously matched both and got no
  letter; it is now 3c. A letter is still withheld when nothing fits uniquely, since the finder
  never shifts the origin to the standard setting.
- Knock-on fix: diamond used to be identified as #225 and was therefore handed **Fm-3m Wyckoff
  letters for an Fd-3m structure**. It now resolves to #227 and gets that group's own letters.
- **First tests for the symmetry code** — there were none, and the one test mentioning a space
  group hard-codes `F-43m` as an LLM-context fixture without ever running the finder. 736 new
  tests across `__tests__/symmetry.test.js` and `__tests__/wyckoff.test.js`, built on **all 230
  space groups** given by ITA generators in `__tests__/fixtures/spaceGroups.js`. Each group is
  asserted three ways: its generator closure reproduces the general-position multiplicity, its
  symbol is built correctly from exact operations, and the finder recovers it end-to-end from a
  structure made by expanding a generic point. Every tabulated Wyckoff position is re-derived from
  the group's own operations at test time, so the committed data cannot drift from the symmetry it
  claims to describe. Plus an integration test over the bundled 52 000-atom demo model
  (GaTa₄Se₈ → F-43m, Ga 4c / Ta 16e / Se 16e ×2), and unit tests for intrinsic translations, glide
  letters, screw handedness, centering, the ladder, and site orbits.
- All reference data was generated, then **kept only where it verified computationally**: the 230
  symbols come from two independent passes that agreed on every entry; each generator set was
  accepted only when its closure size matched (crystal-class order) × (centring multiplicity); each
  Wyckoff position only when its coordinate expanded to exactly the stated multiplicity. That check
  rejected 7 wrong positions — for P222₁ the proposed 2-fold site `x,0,0` lies on no axis and
  expands to 4 points, not 2.
- Known limitation, pinned by test rather than hidden: **I222/I2₁2₁2₁ and I23/I2₁3 cannot be told
  apart.** I centring gives both members of each pair the same element types along the same
  directions (every pure 2-fold has a 2₁ screw half a cell away), and they differ only in where
  those axes sit relative to one another — which needs the origin-aware matching this finder does
  not do. It reports I222 and I23.
**Local Geometry tab — bond angles in the app, both runtimes** (2026-08-06) — the triplets engine
(merged in #36) is now a workspace page: pick A–B–C with **B central**, bound the bond lengths,
Compute. Payload contract is `rmc_toolkits.triplets.bond_angle_summary` (angle histogram with
counts / per-degree density / exact sin-corrected, bond-length histograms per window, coordination
distribution), served identically by the new Flask `/api/triplets` route (lru-cached on
path+mtime+params) and by the shared PCA worker's new `triplets` request — the same
worker-or-API routing as the PCA pages via `useSiteCloud`.

- `workers/triplets.js` is a line-for-line port of the Python engine (linked-cell periodic search
  with explicit image shifts, shared-ends unordered counting, zero-length pairs never bond),
  parity-tested against Python goldens: exact integer histograms, 1e-9 float agreement
  (`triplets_fixture.json`, regenerated by `tests/generate_triplets_fixture.py`).
  `pcaKde.js` now keeps the per-atom element + supercell-fraction list on its parse result so the
  worker answers triplets requests from the same cached parse.
- `LocalGeometryPage.jsx` follows the PCA/Orientation page pattern (top controls bar, `pca-panel`
  cards, InfoBadges): triplet selects seeded per dataset (ends = most abundant element, central =
  next), inclusive windows with an optional distinct B–C window, result chips (bonds, mean length,
  coordination mode share, mean angle), a sin-corrected|density toggle on the hero plot, and a
  window helper that plots the run's partial g(r) (`PDFpartials.csv`, both modes) with dashed
  guides at the current bounds, cropped to the first-shell region. Three equal-width full-height
  columns: angle distribution, partial g(r), folded unit cell. An earlier fourth column held a
  histogram of the bond lengths inside the window — dropped as redundant, since it was the same
  first-shell peak the partial g(r) already plots, only clipped to the window and without the
  surrounding context that makes the guides worth reading. The engine still returns those
  histograms; the result chips report their counts and mean length.
- **The partial-PDF helper carries both bonds of the triplet.** A second curve is drawn whenever
  A–B and B–C are different bond *types* — Ga–Ta–Se brackets a Ga-Ta and a Ta-Se shell, Se–Ta–Se
  only ever has one, and the pair is looked up in either label order. That is independent of the
  **distinct B–C** switch, which governs the *windows* rather than the curves: with it off one pair
  of guides covers both bonds, with it on each window gets its own pair, labelled `A–B`/`B–C` by
  role — the pair names would be identical for a same-element triplet, whose two windows still sit
  on one shared curve. The first-shell crop widens to clear whichever window reaches further.
  Two windows also take the two curve colors (`plotPalette.js`, new module: a component file that
  also exports constants breaks Fast Refresh). Guides consume no palette slot, so shell N is
  `PLOT_PALETTE[N]` — each window is drawn in the color of the shell it brackets, and a lone window
  keeps the neutral guide stroke since there is nothing to tell it apart from.
- **The 3D panel is the folded unit cell, not site ellipsoids** (`FoldedCellPanel.jsx`): every atom
  folded back into one cell as a point cloud, the same view the Atomic Density page shows, so the
  spread around a site is the measured thermal cloud rather than a fitted ellipsoid. The detected
  bonds are drawn over it as thin transparent lines (`LineSegments`, opacity 0.45) instead of
  opaque cylinders, so a full coordination network reads as a framework without hiding the cloud.
  Bond matching is unchanged — average site pairs inside the window, periodic images included.
  `SiteStructurePanel` keeps its ellipsoid view for the PCA and Orientation pages and gains a
  `title` prop rather than being renamed for one caller.
- **`InteractivePlot` gains an opt-in `fit` variant** that takes its `viewBox` from the rendered box
  (ResizeObserver, rounded to whole pixels) instead of a fixed 8:5 or 1440×320 aspect. The fixed
  aspect letterboxed inside a card of any other shape — the two plots here filled 52–57% of their
  cards, the rest dead space. Under `fit` one user unit is one CSS pixel, so margins keep their
  physical size and tick density follows the box (~1 y tick / 70px, 1 x tick / 95px, clamped);
  measured fill is now 100%. Opt-in because it needs a caller that gives the plot a definite
  height: the fixed variants are untouched and every other page renders exactly as before.
- Both plots are handed **identically sized boxes** so the same data does not read at two aspect
  ratios side by side. Equal columns are not enough on their own — it is the chrome above each
  stage that differs, so the card header is pinned to one height (one title wraps, the other does
  not) and two legend rows are reserved in both (only the partial PDF's legend wraps, once the
  split adds a second pair of guides). Measured: 492×363 in both, viewBox matching.
- The bond-angle normalization note is rewritten. It previously led with the correction rather than
  the thing being corrected; it now defines **density** first, states plainly that random bonds do
  *not* give a flat line there (fewer ways to form an angle near 0°/180° than near 90°), and only
  then explains that **sin-corrected** divides that factor out so flat 1 means random.
- Verified live in both runtimes: Flask + `data/5K_try1` Se–Ga–Se gives 4.00 Se/Ga (4-fold 100.0%)
  at 109.4° ± 3.5°; static Demo (GaTa₄Se₈ 250 K) Se–Ta–Se gives 5.95 Se/Ta (6-fold 95.8%) with the
  split-octahedron distribution. Backend: `TripletsApiTests` (octahedron: exact 12×90° + 3×180°,
  400/403 paths); frontend suite 198 tests incl. the new parity file.
- A 41-agent adversarial review confirmed 17 findings, all fixed: bin-count rounding unified to
  half-up (`_bin_count` — Python's banker's `round()` gave a different **payload shape** than
  `Math.round` at e.g. 8° bins), JS bin assignment now replicates numpy's edge-corrected algorithm,
  dataset switches clear results and epoch-guard in-flight computes, triplet seeding keys on the
  sites payload (not the dataset key racing ahead of it), both app boundaries cap `rmax ≤ 15 Å` and
  `binWidth ≥ 0.05°` and reject half-specified B–C windows (the engine itself stays unrestricted for
  CLI/library use), element case is normalized before the route cache, and assorted UI polish
  (input styling, loading states, debounced helper guides). One documented residual: libm vs V8
  `acos` can differ by 1 ulp, so a *bitwise-ideal* geometry whose cosine sits exactly on a bin edge
  (cos = 0.5 in an undisplaced average configuration) may shift one count between adjacent bins
  across engines — unreachable with real RMC data.

**Bond-angle (triplet) distribution engine + `rmc-triplets` CLI** (2026-08-05) —
`rmc_toolkits/triplets.py` (source of truth, module docstring is the math reference) +
`rmc_toolkits/triplets_cli.py`, following RMCProfile's `triplets_new_bonds_sinth` workflow: pick an
A–B–C triplet with **B the central atom**, bound the A–B and B–C bond lengths (inclusive windows),
and histogram the angle at B. Engine-only by design — offline validation comes before any Flask
route, browser port, or page.

- Neighbours come from a linked-cell search that carries **explicit periodic-image shifts** instead
  of assuming minimum-image: exact for any cell shape (triclinic included) and for boxes smaller
  than the cutoff, where several images of one atom are genuine distinct neighbours. Equivalent
  ends (same element + window) count each unordered bond pair once — an octahedron gives C(6,2)=15
  angles; distinct windows count ordered 1→2/2→3 assignments minus the bond-with-itself
  combinations. All pair formation is vectorized (ragged group cartesian products, no per-atom
  Python loop): the 52 000-atom GaNb₄Se₈ box runs in well under a second per triplet.
- Three curves per histogram: raw `counts`; `density` (per-degree, unit integral); and
  `sin_corrected` — the count fraction over the *exact* isotropic bin fraction
  `(cosθ₁−cosθ₂)/2`, i.e. the "sinth" normalization computed from the bin integral rather than
  `1/sin(θ_center)`, so the 0°/180° bins stay finite and an isotropic configuration reads flat 1.0.
- `rmc-triplets` (console script; `python -m rmc_toolkits.triplets_cli`) takes an `.rmc6f` or a run
  folder, writes a commented CSV (angle, counts, density, sin-corrected) plus an optional PNG and
  raw-angle dump, and prints bond counts, mean lengths, per-central coordination and angle stats.
  Existing outputs are never overwritten without `--force` (the `rmc-autoscale` convention).
- A 14-agent adversarial review (brute-force refutation over ~175 randomized/hand-built
  configurations) confirmed the periodic geometry and counting; its three surviving minor findings
  are fixed: zero-length pairs are never bonds (`rmin=0` + bitwise-coincident atoms used to send a
  0/0 NaN angle past the histogram), the CLI reports unwritable outputs and truncated `.rmc6f`
  headers as clean errors instead of tracebacks, and overwrites now require `--force`.
- `tests/test_triplets.py` (27 tests): exact geometry fixtures (octahedron 12×90°+3×180°), bonds
  through the periodic wall, multi-image bonds in a box smaller than the window, self-image bonds
  of the central element, wrap invariance, **exact agreement with a brute-force all-images
  reference in a skewed triclinic cell**, isotropic sin-correction flatness, normalization and
  validation errors, loader parity, CLI end-to-end. Numpy-on-Accelerate note: large `(P,3)@(3,3)`
  matmuls emit spurious divide/overflow warnings on Apple silicon, so the engine uses `einsum` for
  the displacement→Cartesian map (bit-identical result, no BLAS warning path).
- Physics check on `data/5K_try1` (GaNb₄Se₈): Se–Ga–Se with a 2.1–2.7 Å window peaks at 109.4°
  with 3.99 Se per Ga (GaSe₄ tetrahedra); Se–Nb–Se with 2.2–2.9 Å shows the distorted-octahedron
  splitting (~76/88/104/160°) at 5.78 Se per Nb. Next steps when validation settles: Flask route,
  `workers/` port (parity-tested), and a page.

**Algorithms and math reference** (2026-07-26) — `docs/ALGORITHMS.md` plus seven per-page documents
under `docs/algorithms/`: a code-anchored account of every mathematical operation each page performs
on the user's data, so a reader can audit how a plot, density map, symmetry label, scaled dataset or
direction map was produced. ~16 000 lines.

- The hub carries a page index, a pipeline diagram, a table of *where* each computation runs (Python
  package / Flask route / browser worker), and a consolidated list of every limitation and
  approximation, each linking to the section that explains it. `algorithms/notation.md` reconciles
  the symbol table across sections and lists the symbol *collisions* rather than merging them.
- Two conventions the reference commits to: **the code wins** wherever document and source disagree,
  and *reference-grade* (the float64 Python path) is distinguished from *visualization-grade* (the
  browser ports), with measured cross-engine tolerances quoted — and the gaps named where no such
  test exists (`pcaKde.js` and `orientation.js` are pinned to their own in-language references, not
  to Python; `localKdeWorker.js` has no numerical reference test at all).
- README gains a "Math Under the Hood" section with the signature equation per analysis page — fit
  residual, 2-D Gaussian KDE, ADP tensor with its $\chi^2$ probability ellipsoid, solid-angle
  enhancement — each with the caveat that makes the number honest. QuickStart, REFERENCE and AGENTS
  point at the reference from their own entry points.

**Rwp reports "unavailable" instead of a fake perfect fit** (2026-07-25) — `rwp()` returned `0`
whenever its denominator was falsy. For an observed column that is entirely NaN (a fully masked
dataset) the denominator is NaN, `!NaN` is `true`, and the dashboard chip showed `Rwp 0.000` — read
as a perfect fit over data that isn't there. A genuinely zero denominator collapsed to the same `0`.

- Both implementations now sum only the points where the observed *and* fitted values are finite,
  and return `None` / `null` for the two degenerate cases (no finite pair, zero denominator).
  `rmc_toolkits/parsers.py` stays the source of truth; `browserData.js` is the static-mode port —
  keep them in sync (see AGENTS.md). `PlotResult.metrics` is now `dict[str, float | None]`, so
  `/api/plot/metadata` serializes the sentinel as JSON `null`, matching static mode.
- `Dashboard.jsx` renders the sentinel as `Rwp —` via a single shared `renderRwpChip` helper (the
  two chip sites had duplicated the formatting). The assistant's run context already skipped
  non-finite Rwp values, so it drops the dataset's `rwp` key as before.
- Tests: all-NaN observed, partially-NaN columns on either side, and a zero-denominator column, in
  both `tests/test_parsers.py` and `src/__tests__/browserData.test.js`.

**PCA Ellipsoid: principal axes in the crystallographic frame** (2026-07-25) — the Displacement
statistics panel gains a fourth column, **Crystal orientation**: for every principal axis, its angle
to the unit-cell edges a, b, c (the closest one shaded, same accent wash the covariance diagonal
uses) and the crystallographic direction `[u v w]` it runs along. `src/pcaCrystalFrame.js` had
computed all of this since v0.4.0 but nothing imported it — only `unitCellVectors` was wired up, so
the app knew which lattice direction each displacement axis pointed along and never said so.

- Angles are between Cartesian directions (both frames already share one Cartesian basis), so they
  are well defined for any cell, oblique or not. `[u v w] = M⁻ᵀ·axis`, normalised so the largest
  component is 1, is a **direct-lattice** direction — in an oblique cell it is not normal to the
  like-indexed (h k l) plane. The panel's InfoBadge says both, since the distinction is easy to
  misread off a table of numbers.
- An eigenvector's sign is arbitrary (v and −v are one axis), so the new
  `crystalOrientationRows` reports the member of the ± pair whose closest crystal axis is the
  *acute* one — angles and `[u v w]` flip together, keeping each row internally consistent, and
  "dominant" then literally means smallest angle. Unit-tested in `pcaCrystalFrame.test.js`
  (already-acute rows untouched, a 170°-from-a axis reported as 10° with the direction negated, no
  obtuse dominant angle for an oblique cell, singular cell → null).
- The column appears only when the file carries lattice metadata; the statistics grid then switches
  to four tracks (`.has-crystal`, sized to avoid inner scrollbars down to a 1280-wide screen) and
  still stacks to one column below 980 px.
- Also: the header wordmark centers "WORKBENCH" under "RMCProfile" instead of leaving it flush-left,
  so the two lines read as one lockup.

**Displacement-orientation engine (hex-binned sphere)** (2026-07-24) — new engine for the
distribution of displacement *directions* per site: `u = Δr/|Δr|` binned in solid angle on
a Goldberg polyhedron (the dual of a frequency-ν geodesic icosahedron — hexagons plus the
12 pentagons any hexagonal tiling of a sphere must carry; 10ν²+2 cells), each cell's count
divided by its exactly-computed solid angle. Complements the ellipsoid view: discrete
hop-sites and antipodal asymmetry (odd anharmonicity) are visible here and invisible in U.

- `rmc_toolkits/orientation.py` (source of truth) + `workers/orientation.js` (static-mode
  port, keep in sync): tiling construction with derived (not hard-coded) icosahedron faces
  and combinatorial-key merging, O(1) direction→cell lookup (gnomonic seed + greedy walk,
  brute-force-verified), exact centrosymmetric antipode map.
- Histogram options: auto resolution via `recommended_frequency` (over-binning guard,
  ~12 pts/cell), weights `count | amplitude | amplitude2` (the latter's angular map is the
  decomposition of ⟨u²⟩), amplitude cutoffs (absolute + quantile — near-zero |Δr| has
  noise for a direction), mass-conserving neighbour smoothing, `cartesian | pca` frame
  (PCA axes shared with the ellipsoid engine via the same canonicalization). The map is
  deliberately never antipodally folded: the +u/−u imbalance is the signal.
- Outputs: `enhancement = 4π·density` (1 = isotropic, reads as "N× more likely than
  chance"), per-cell Poisson `zScore` vs the isotropic null (smoothing never feeds the
  z), `antipodalAsymmetry` = Σ_pairs|n(u)−n(−u)|/N with its Poisson noise floor
  `antipodalAsymmetryNull` (inversion asymmetry U is blind to), orientation tensor
  ⟨u uᵀ⟩ + eigenvalues + `orientationAnisotropy` (3λ₁−1), cell polygons for rendering.
- Wiring: `/api/pca/orientation` (Flask), `{kind: 'orientation'}` in `pcaKdeWorker.js`,
  package exports. Tests: `tests/test_orientation.py` and
  `workers/__tests__/orientation.test.js` — tiling invariants, brute-force assignment
  parity, isotropic flatness, lobe recovery, one-sided asymmetry detection, weight/frame
  behavior, API + worker routing.
- **Renamed the nav tab "Orientation" → "Displacement Directions"** (2026-07-24, same
  branch): clearer, and the direct direction-counterpart to "PCA Ellipsoid" (amplitude /
  shape). The internal route key stays `orientation`. Also: the Axis-views hover now
  darkens with a real ~38% scrim (a weak dark overlay over the light panel just greyed
  toward white, and the global `button:hover` had out-specified it), and the a/b/c
  labels are larger.
- **Orientation controls: contrast knob + defaults** (2026-07-24, same branch): new
  **Contrast** slider (0.5–3×) applies a symmetric color gain about the isotropic level
  (`colorCoordinate` in `orientationSphere.js`: `t = clamp(pivot + contrast·(v/vmax −
  pivot))`, `pivot = 1/vmax`), so faint lobes and depletions stand out; contrast = 1 is
  exactly the old linear mapping (backward compatible), and `colorbarGradient` paints
  each stop through the identical transfer so the bar and the sphere always agree
  (unit-tested). Defaults changed to **ν = 10** (was Auto) and **2× smoothing** (was 0)
  for a legible out-of-the-box map, and the smoothing slider now reaches **12×** (was 4).
  The redundant "N cells (ν=…) / N vectors" summary line is removed (both are already in
  the Resolution control and the panel header). Axis-views hover is now a translucent
  accent tint + ring instead of a near-white wash.
- **Orientation page — three-panel layout + PCA-page design parity** (2026-07-24, same
  branch): the page now matches the PCA Ellipsoid page's conventions — all options live
  in the top controls bar, and the header actions (Crystal|PCA frame toggle, Reset view,
  Save) sit in each panel's own header. Three equal-height panels in one grid at
  **3 : 6.5 : 6.5** — Axis views (the fixed-angle mini column, now its own panel) :
  sphere : site picker; `OrientationView` renders `display:contents` so its two panels
  drop straight into the grid. Panel height is viewport-clamped
  (`calc(100vh − 17.75rem)`) so the three fill a 16:9 screen down to the footer with zero
  scroll (verified 1600×900). The mini views zoom out (camera pushed to 5.6×) so a
  fully-inflated relief surface never clips the pane edges, and the main sphere's default
  view is pulled back too. The axis rods now show **only the selected coordinate frame**
  (a/b/c in Crystal, PC1/2/3 in PCA) — the mini views, main sphere, and legend all show
  one system at a time instead of overlaying both.
- **Orientation page 16:9 layout + shared site picker** (2026-07-24, same branch): the
  clickable Site-ellipsoids unit-cell picker is extracted from PcaKdePage into
  `SiteStructurePanel.jsx` (axis palettes/builders into `sceneAxes.js`) and now sits on
  the Orientation page too — sphere panel wide left, picker right
  (`.orient-page-layout` grid; the panels' named grid-areas from the PCA layout are
  neutralized there). The fixed-angle mini views moved from a bottom row to a **column
  left of the sphere** (scissor panes stacked vertically; WebGL y counts from the
  bottom). The sphere row height is viewport-clamped (`calc(100vh − 30.5rem)`, self-sized
  `flex: none` so a zero-height flex ancestor can't collapse it) — on a 16:9 monitor the
  whole page (site bar, controls, sphere + minis + picker, readouts, footer) fits with
  zero scroll, verified at 1600×900.
- **Orientation promoted to its own workspace page** (2026-07-24, same branch): the
  histogram is not a PCA product, so the Density | Orientation tab coupling is gone —
  `OrientationPage.jsx` is a top-level nav page with its own site picker, and the PCA
  Ellipsoid page is back to its pre-tab layout. The shared plumbing (text loading,
  worker/API routing, sites table, selected site) moved to the `useSiteCloud` hook with
  one app-lifetime worker, so both pages share the parse cache. The sphere view gains
  **three fixed-angle mini views** (down a/b/c, or PC1/2/3 in the PCA frame): one extra
  scissored renderer re-rendered only on scene changes (static cameras, no per-frame
  cost); clicking a mini snaps the main camera to that axis view.
- **Amplitude relief** (2026-07-24, same branch): the engines report
  `cellMeanAmplitude` — the mean |Δr| of each cell's movers (smoothing applies to the
  numerator/denominator sums, not the ratio; empty cells report 0) — and the sphere gains
  a **Relief** slider that bulges each cell radially by its mean |Δr| relative to the
  site average. Color (how often) and shape (how far) carry independent information.
  Shared polygon vertices average their cells' factors so the relief surface is
  crack-free (unit-tested); cell borders follow the surface; the hover readout adds
  ⟨|Δr|⟩.
- **Adversarial review pass** (multi-agent, 2026-07-24) — confirmed findings fixed:
  auto-resolution rounding divergence (Python `round()` is half-to-even, JS used
  `Math.round` half-up — at exactly 774 surviving points the two engines picked
  different tilings; JS now rounds half-to-even and both suites pin
  `recommendedFrequency(774) == 2`); stale hover index could crash the tooltip when a
  resolution change shrank the tiling (hover now resets with each result + bounds
  check); OrientationView unmount leaked its geometries and a WebGL context per tab
  switch (groups disposed + `forceContextLoss` in the scene cleanup); the
  density-integrates-to-one tests were circular (areas cancel by construction) and are
  replaced by a Monte-Carlo cross-check that assignment Voronoi fractions × 4π
  reproduce the analytic polygon areas. Greedy-walk assignment verified exact against
  brute force at ν = 24/40/64 (0/20000 mismatches, all centers self-assign).
- **UI: Orientation tab** in the PCA Ellipsoid main panel (`OrientationView.jsx` +
  `orientationSphere.js` pure helpers, unit-tested), sharing the page's site picker.
  Density | Orientation tabs in the panel header (the density canvas hides via CSS, never
  unmounts — its once-mounted scene effect survives tab switches). Flat-shaded Goldberg
  cells colored by enhancement, cell borders, PC1/2/3 + a/b/c axis rods, Crystal ↔ PCA
  frame toggle, per-cell hover readout (direction, ×isotropic, count, z), colorbar with a
  1× tick, and a summary strip: cells/ν, used vectors, peak (+z), orientation anisotropy,
  ± asymmetry vs its Poisson floor (flagged red when > 3× floor), map significance.
  Verified live in both runtimes (Flask `data/RMC`, browser Demo run).

**Auto StoG tab hidden** (2026-07-23) — the page is not behaving correctly yet, so it is
gated out of the UI while the work continues offline. `App.jsx` gains a `SHOW_AUTO_STOG`
constant (currently `false`) guarding the nav button and the page mount; with it off the
page is unreachable, since `activePage` defaults to `dashboard` and that button was the
only route to it. Nothing was removed — `AutoStogPage.jsx`, `workers/autoScale*.js`,
`rmc_toolkits.scaling`, the `rmc-autoscale` CLI, and `/api/scaling/*` are untouched, and
their tests still run. Flip the constant to restore the tab.

**Auto StoG independence + ρ₀ self-consistency** (2026-07-18) — the page is now true
pre-processing, decoupled from the run folder, and the number density can come from the
data itself.

- **ρ₀ self-consistency** (`estimate_rho0` in `scaling.py`, JS port in `autoScale.js`,
  CLI `--estimate-rho0`): the density-limit amplitude depends on ρ₀ (C2 rows scale with
  the density line −4πρ₀r) — degenerate with the scale from low-r alone — while the Q→0
  Faber-Ziman amplitude needs only ⟨b²⟩ + the measured level. ρ₀ is recovered as the root
  of `a_fz/a_density(ρ₀) = 1` by fixed-point iteration (`ρ ← ρ·concordance`, 2–4
  `autoscale` passes). Requires a composition; `extrapolated` flags Qmin beyond the FZ fit
  width. Validated: synthetic truth 0.05 recovered within 2% from seeds 0.02–0.2; FeCoSn
  x-ray 199 K gives 0.0600 vs the hand 0.057329 (<5%) from any seed. ρ₀ resolution order
  is now: value → `NUMBER_DENSITY ::` header → mass density + composition →
  self-consistent estimate (auto-run when ρ₀ is left empty; explicit "Estimate ρ₀" button
  too) — so composition + Q window are the only required inputs.
- **Auto StoG page decoupled from the run folder** (`AutoStogPage.jsx`): pre-processing
  now has page-local file handling — an upload/dropzone (S(Q) ± stog.inp, multi-file;
  stog output files filtered out) replaces the shared run-folder listing, and the page is
  fully client-side in BOTH runtimes (worker engine + zip export; the Flask
  `/api/scaling/*` endpoints stay for API/CLI use but the page no longer calls them).
  Q window prefills from the data's finite extent; `App.jsx` passes no props.
- **Parameters grouped with descriptions**: the Advanced panel is five fieldsets, each
  with a one-line purpose and per-field tooltips — *Amplitude & offset* (High-Q mode,
  amplitude criterion, Robust/σ/Despike), *Coefficients* (⟨b⟩², ⟨b²⟩ x-ray overrides),
  *Transform* (r-cut, grid, Lorch, low-Q corr.), *Low-r region* (r₀, fit window,
  enforcement), *Fixed scaling* (a, b).
- Tests: `estimate_rho0` synthetic + x-ray + CLI coverage (Python), fixture-backed JS
  parity (tolerance bounded by the rtol stopping rule, not single-pass precision);
  browser flow verified end-to-end in static mode (upload → estimate 0.060664 →
  auto-scale a = 1.2511, concordance 1.00 → 9-file zip).
- **Correctness review of the page** (multi-agent sweep of every description + flow +
  export against the engine; 12 findings, 8 confirmed after adversarial verification):
  manual runs now recover r₀ by detection so a checked "Enforce low-r" with an empty
  cutoff is honored (was a silent no-op and a CLI-parity break); a fixed-(a, b) run no
  longer fails the unused FZ-amplitude validation; "Run fixed" gains the missing
  `estimating` concurrency guard; three tooltips corrected to the code's real resolution
  orders (ρ₀ chain includes stog.inp above the data header; r₀ is header → stog.inp peak
  window → detection; x-ray needs ⟨b⟩² = 1 *and* ⟨b²⟩); export title lines use %g-style
  numbers (the `(browser)` marker vs the CLI's version token is intentional); provenance
  JSON now records `stogInpReference` (the loaded stog.inp's hand scaling), the fit
  `history`, and `c1ModeEffective`, matching the CLI payload.
- **Discordant-data guardrails** (found live on the Mn₃Sn neutron runs, whose low-Q
  hole breaks the density-limit amplitude): `estimate_rho0` stops on non-physical
  amplitudes instead of iterating to a bound, and an unconverged estimate is *refused*
  everywhere (CLI error, worker error, page banner) with guidance to set ρ₀ explicitly
  and use the FZ amplitude — previously the page silently fit with a garbage density
  (a = −1.02 on PG3_54139). The neutron-vs-x-ray state is now always visible: a chip
  shows the coefficients in effect and their source, warns when Advanced overrides
  shadow a typed composition, and a "Reset params" button clears a previous sample's
  settings. Folder drops report cleanly instead of failing silently. Engine + CLI
  guard tests (Mn₃Sn-fixture-gated).

**Composition-first Auto StoG** — procedure clarified end-to-end
([SCALING_PROCEDURE.md](SCALING_PROCEDURE.md)); the user now provides only the chemical
composition + [Qmin, Qmax] (+ a density when the data header has none); validated on three
new complete Fortran neutron runs (POWGEN Mn₃Sn).

- **The procedure document** (`docs/SCALING_PROCEDURE.md`): the definitive recipe for
  absolute-scale S(Q)/G(r) — inputs and defaults tables, the Keen-convention function map,
  the 7-step pipeline, how to read the one-sided density flag and the concordance trust
  metric, and material-class guidance. Research base: pystog 0.6.7 source (operation order
  confirmed read → merge → manual scale → transform → filter → Lorch → Keen conversions;
  no auto-scaler, low-r minimizer an explicit TODO) and ADDIE (density conversions,
  periodictable molar masses, `<b_coh>^2` hand-off).
- **Composition-aware omitted-low-Q correction**: the analytic [0, Qmin] correction now
  extrapolates S(Q) to the composition-derived Keen Eq. 21 target
  ``S(0) = 1 − ⟨b²⟩/⟨b⟩²`` instead of pystog's solid-state 0 — algebraically a
  one-line change (``const' = (1 − s0)·const + s0·coef``, still affine in (a, b)). For
  negative-b compositions this is O(1): Mn₃Sn has ⟨b²⟩/⟨b⟩² = 13.06 → S(0) = −12.06, and
  the composition-aware target cuts the low-r residual ~40% on the PG3 runs.
  ``ScalingConfig.s0_target`` (auto from ``b_sq_avg``; pin 0 for classic parity).
- **Data-derived first-shell r₀** (`detect_first_peak_onset` + a second refinement pass in
  `autoscale`): the dominant |g| feature's left flank (35% of peak) — |g| because
  negative-b totals can have an *inverted* first shell; peak-relative because sub-r₀
  ripples scale with the amplitude. Detected 2.73–2.77 Å on the Mn₃Sn runs and 2.53 Å on
  FeCoSn (hand-chosen classic cutoffs: 2.40–2.68). Sets the low-r fit window without a
  structural prior and becomes the default low-r **enforcement cutoff** in data mode
  (CLI + API + page; explicit `--no-enforce`/`enforce:false` still opts out);
  `diagnostics_summary` reports `r0_detected`/`window_refined` and the *effective* fit
  window from provenance.
- **ADDIE-style density toolkit** (`scattering.py` + JS port): `ATOMIC_MASS_U`
  (periodictable/CIAAW values for the full element set), `molar_mass`,
  `number_density_from_mass_density` / inverse (ρ₀ = ρ_m·N_A/10²⁴·n/M). CLI
  `--mass-density`, API `massDensity`, page "or ρ g/cm³" field. Mn₃Sn check:
  0.063049 atoms/Å³ ↔ 7.421 g/cm³.
- **Mn₃Sn neutron validation** (`data/stog_tests/{stog,stog_300K,stog_500K}`, local):
  three complete Fortran runs; parameters recovered from the outputs to ~1e−12
  (hand scalings ×2.5/×2.05/×10 — all satisfying b = 1 − a, i.e. the level-anchored
  decomposition the sweep formalizes; ⟨b⟩² 0.015407 = `faber_ziman("Mn3Sn")` exactly).
  New `Mn3SnNeutronTests` (+ detection/s0 unit tests): manual filter-stage parity ≤2e−3
  rms per run, composition/density round trips, first-shell detection, and the headline
  physics finding — the density limit is degenerate on this material (one-sided flag
  False, hand values mutually inconsistent by 5×) while the **FZ amplitude injects the
  composition's S(0) and lands at O(10)** consistently; the procedure doc encodes that
  decision rule.
- **Page**: primary bar is now Composition · Qmin · Qmax · ρ₀/mass-density with a live
  ⟨b⟩²/⟨b²⟩/S(0) chip; ⟨b⟩²/⟨b²⟩ overrides moved to Advanced (x-ray note kept); a
  collapsible "How Auto StoG works" explainer; a First-shell-r₀ readout card; and
  **full-range G_K(r) and D(r)** (no 8 Å slice — theory guide lines stay confined to the
  low-r region; the G_K default y-zoom keeps the −⟨b⟩² level readable, box-zoom for
  detail). JS engine/worker fully synced (s0-aware correction, detector, two-pass,
  converters) with new Python-golden parity tests (detection case + Mn₃Sn converters).
  Verified live on the 300 K run: composition-only inputs → a = 1.309, 7 iterations,
  r₀ 2.750 Å detected + enforced, density flag red, concordance 7.95 flagged.

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
