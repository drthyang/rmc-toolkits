# Automatic Total-Scattering Data Scaling — Build Plan

Plan for the "dedicated preprocessing module" that ROADMAP Phase 7 defers to: automatic
scale/offset determination for measured S(Q), STOG-compatible Fourier transform and Fourier
filter, and RMCProfile-ready outputs. Replaces STOG's manual `try again` loop with a minimizer
wrapped in a self-consistent loop.

**Golden reference:** D. A. Keen, *J. Appl. Cryst.* **34**, 172–177 (2001),
[doi:10.1107/S0021889800019993](https://doi.org/10.1107/S0021889800019993). All function names,
normalizations, and limits below cite its equation numbers.

**Validation example:** `data/stog_tests/stog_59438/` — a complete real STOG run (POWGEN run
59438): raw + rebinned S(Q), the actual `stog.inp`, every intermediate (`scale.fq`, `scale.gr`,
`ft.dat`, `scale_ft.*`), and the final RMC outputs (`scale_ft_rmc.{fq,gr,dr}`).

---

## 1. Verified math spec

Everything in this section is verified either against the Keen paper directly (eq. numbers) or
by exact numerical checks against the example run (marked ✓data, agreement to machine
precision unless noted).

### 1.1 Functions and conventions (Keen 2001)

| Function | Definition | Keen Eq. |
| --- | --- | --- |
| F(Q) | `Σ_ij c_i c_j b_i b_j [A_ij(Q) − 1]` (barns) | (9) |
| G(r) | `Σ_ij c_i c_j b_i b_j [g_ij(r) − 1]` | (10) |
| FT pair | `F(Q) = ρ0 ∫ 4πr² G(r) sin(Qr)/(Qr) dr`; inverse with `1/((2π)³ρ0)` prefactor | (11), (12) |
| S(Q) | `S(Q) − 1 = F(Q) / ⟨b⟩²` where `⟨b⟩² = (Σ_i c_i b_i)²` | (19) |
| D(r) | `D(r) = 4π r ρ0 G(r)` | (29) |

### 1.2 Limits — the scaling constraints

| Limit | Value | Keen Eq. | Status |
| --- | --- | --- | --- |
| S(Q→∞) | 1 | (21) | paper |
| F(Q→0) | `−⟨b²⟩ + η`, `⟨b²⟩ = Σ_i c_i b_i²`, η = compressibility term (ignorable for dense solids) | (14) | paper |
| S(Q→0) | `1 − ⟨b²⟩/⟨b⟩²` (0 only for monatomic; slightly negative otherwise) | (21) | paper |
| G(r < r0) | `−⟨b⟩²` (flat level; r0 = closest interatomic approach) | (15) | paper + ✓data |
| D(r < r0) | `−4π ρ0 ⟨b⟩² r` (straight line through 0) | (29)+(15) | paper + ✓data |

Caution: ⟨b²⟩ (Q→0 limit of F) and ⟨b⟩² (normalization and low-r level) are different numbers
whenever the sample has more than one element. Do not conflate them.

### 1.3 STOG file semantics (decoded from the example, ✓data)

Annotated `stog.inp` from the example:

```
1                       # number of data files
PG3_59438_SQ_rebin.dat  # rebinned S(Q): uniform ΔQ grid (0.01), NaN below data Qmin
1.0 28.0                # Qmin Qmax used for scaling + FT
-9 0.1                  # yoffset yscale  → S_scaled = S_raw/yscale + yoffset   [✓data, 16 digits]
0                       # Qoffset
scale.fq                # scaled S(Q)  (S-convention despite .fq name; oscillates about 1)
scale.gr                # unfiltered transform of scale.fq
50                      # rmax
5000                    # n r-points (Δr = 0.01)
N                       # Lorch window off
0.063049                # ρ0 (Å⁻³)                                              [✓data via D(r) slope]
0                       # second yoffset (0 in example; semantics TBC — see §6)
N                       # "try again" manual loop → replaced by our minimizer
Y                       # Fourier filter on
1.0                     # r_cutoff for the filter
scale_ft.sq             # filtered S(Q)
scale_ft.gr             # filtered G(r) (2 value cols: G-like, D(r) = 4πρ0 r · col1  [✓data])
0.015407                # "Faber-Ziman coefficient" = ⟨b⟩² = (Σ c_i b_i)² barns   [✓data: −⟨b⟩² is
                        #   exactly the flat low-r level of scale_ft_rmc.gr]
scale_ft_rmc.fq         # F(Q) = ⟨b⟩²[S(Q) − 1]  (→ 0 at high Q  [✓data])
scale_ft_rmc.gr         # Keen G(r)   (flat −⟨b⟩² below r0  [✓data])
scale_ft_rmc.dr         # D(r) = 4πρ0 r G(r)  (slope −4πρ0⟨b⟩² = −1.2207e−2 Å⁻¹  [✓data])
2.48 2.65 3.1           # cutoff, rmin, rmax of 1st peak (final ripple clean-up)
```

File conventions, pinned empirically (pystog cross-run, 2026-07-17): `scale.gr` and
`scale_ft.gr` column 1 hold **g(r) − 1 ≡ G_K(r)/⟨b⟩²** (dimensionless, oscillates about 0);
`scale_ft.gr` column 2 = D(r). The `rmc` outputs are Keen G_K(r) (barns), F_K(Q), and D(r).

Verified relations (all exact against the example):

- **Manual scaling formula**: `S_scaled(Q) = S_raw(Q)/yscale + yoffset`. With the example's
  `-9 0.1` this is `10·S_raw − 9`, i.e. the asymptote-preserving form `a·(S−1) + 1` with a = 10:
  the colleague amplified the oscillation about 1 tenfold without moving the asymptote.
- **Fourier filter**: `scale_ft.sq = scale.fq − (ft.dat − 1)`. `ft.dat` is the Q-space
  back-transform of the unphysical r < r_cutoff content, expressed as an S(Q)-like function
  about 1. Our implementation must reproduce `ft.dat` itself (it's in the fixture — a direct
  oracle for the filter, independent of any documentation).
- **RMC conversions**: `scale_ft_rmc.fq = ⟨b⟩²·(scale_ft.sq − 1)`;
  `scale_ft_rmc.gr` = Keen G(r); `scale_ft_rmc.dr = 4πρ0 r · scale_ft_rmc.gr`.
- These output names are already recognized by `detect_plot_kind()` → `"stog"` in
  `rmc_toolkits/plots.py` and by the frontend `readStog()` — the existing dashboard renders our
  outputs with zero changes.

### 1.4 The auto-scaling problem

Model the correction as `S_corr(Q) = a·S_meas(Q) + b` (two unknowns). **Convention: multiply.**
Classic `stog_new` *divides* (`S/yscale + yoffset`) while pystog *multiplies* (`y·Scale + Offset`)
— a known footgun. Our API uses the multiply form everywhere; the `stog.inp` reader converts
(`a = 1/yscale`, `b = yoffset` — the example's `-9 0.1` ⇒ a = 10, b = −9) and the docstrings
state it. Physics gives two
independent target regions:

- **C1 (high-Q asymptote, Eq. 21):** over a tail window `[Q_hi, Qmax]`, `S_corr` should
  oscillate about 1.
- **C2 (low-r level/line, Eq. 15/29):** after transforming, `G_corr(r)` should sit on `−⟨b⟩²`
  (equivalently `D_corr(r)` on `−4πρ0⟨b⟩²r`) over a window `[r_cut + δ, r0 − δ]` below the
  first peak.
- **C3 (Q→0, Eq. 14):** `F_corr(Q→0) → −⟨b²⟩` — used as a *diagnostic* (reported, not fitted),
  since measured data rarely reach low enough Q (example starts at Q ≈ 0.96 Å⁻¹).

Because the sine transform is linear in the data, `G[a·S + b]` is affine in (a, b):
`G_corr(r) = a·G_data(r) + (a + b − 1)·K(r)` where `G_data = FT[S_meas − 1]`-style and the
offset kernel `K(r) = FT[1]` is computed **numerically on the same Q-window/grid/window-function
as the data transform** (not analytically), so truncation effects cancel self-consistently.
Both basis transforms are precomputed once per window choice — the joint C1+C2 objective is
**linear least squares** in (a, b): one closed-form solve, no iteration needed for fixed
windows. Discretization: direct trapezoidal sine transform (`np.trapezoid`), not FFT — grids
are short and the r-grid is arbitrary. The self-consistent loop
(below) is then only about the interaction between scaling and the Fourier filter.

### 1.5 The self-consistent loop

```
parse + rebin →
  k = 0: fit (a, b) on unfiltered data      (closed-form LSQ over C1 ∪ C2)
  repeat:
    apply (a, b) → transform (+ omitted-low-Q analytic correction) →
    Fourier filter (r < r_cut) → filtered S(Q)
    re-fit (a, b) on the *filtered* data
  until |Δa|/|a| < 1e−6 and |Δb| < 1e−6   (max 20 iterations; record trajectory)
→ final transform → optional low-r enforcement + first-peak cleanup (stog parity)
→ RMC conversions → outputs + provenance JSON
```

Two stages added after the cross-validation run (§3.1):

- **Omitted-low-Q correction**: data start at Qmin (0.96 Å⁻¹ in the example) and the missing
  [0, Qmin] range biases the low-r transform. pystog carries the analytic correction
  (`_low_x_correction`, transformer.py); classic stog applies one too. Include it; without it
  the density-limit window is contaminated (observed: g−1 mean ≈ +1.1 instead of −1 in
  [1.2, 2.2] even for the hand-scaled example).
- **Low-r enforcement (final stog step, config-gated)**: classic stog *hard-replaces*
  r < cutoff (the first number of the `2.48 2.65 3.1` line) with the exact theoretical values —
  verified: `scale_ft_rmc.gr` = −0.015407 to 16 digits for every r ≤ 2.48, first deviation at
  2.49. So the published RMC files satisfy the limits *by construction*, not as evidence of a
  good fit. Implement as an explicit, clearly-labeled final stage (`enforce_low_r=True` for
  stog parity), with the *pre-enforcement* low-r residual always reported — that residual is
  the honest fit-quality metric.

Implementation: closed-form linear solve as the default inner step, wrapped in an optional
`scipy.optimize.least_squares` path (per-point weights from the data's σ column; optional
`soft_l1` robust loss for Bragg-contaminated tails; extensible to extra parameters — e.g. a
Placzek-like linear-in-Q term — without changing the architecture). Expected convergence: 2–4
iterations (the filter perturbs the low-r region only weakly once the scale is near-correct).

---

## 2. Module layout

Following the repo's flat-module package style and the modular-design convention (small,
focused, separately testable pieces):

```
rmc_toolkits/
  transforms.py    NEW — Keen-convention conversions (S↔F↔F_RMC; G↔D↔T; Eq. 19/29/30) +
                   discrete sine-FT pair on uniform grids (Eq. 11/12, 45/46 forms) +
                   Lorch window + Fourier filter. Pure numpy. No file I/O.
  scaling.py       NEW — the auto-scale engine: windows, closed-form (a,b) solve, scipy
                   minimizer wrapper, self-consistent loop, diagnostics report dataclass,
                   provenance dict. Calls transforms.py. No file I/O.
  parsers.py       EXTEND — readers: stog.inp, STOG .dat (count-header + 2–3 col + NaN),
                   .dat metadata header (TITLE::, NUMBER_DENSITY::, MINIMUM_DISTANCES:: —
                   same format as the shipped demo GTS_250K.dat); writers: STOG-format
                   outputs (count + title line + Fortran-style columns).
  plots.py         EXTEND (small) — a scaling-diagnostics figure builder (S(Q) with the
                   1-asymptote and −⟨b²⟩ marker; G(r)/D(r) with theory lines; per-iteration
                   (a, b) trajectory).
tests/
  test_transforms.py   NEW
  test_scaling.py      NEW
```

Public API sketch:

```python
@dataclass
class ScalingConfig:
    qmin: float; qmax: float
    rho0: float                  # Å⁻³   (auto-read from .dat header when present)
    b_avg_sq: float              # ⟨b⟩² barns ("Faber-Ziman coefficient" in stog.inp)
    b_sq_avg: float | None = None  # ⟨b²⟩ barns — enables the C3 diagnostic
    r_cutoff: float = 1.0        # Fourier-filter cutoff
    r0: float | None = None      # closest approach; default from MINIMUM_DISTANCES
    fit_offset: bool = True      # b in the linear correction
    q_tail_frac: float = 0.15    # C1 window = top fraction of [qmin, qmax]
    rmax: float = 50.0; nr: int = 5000
    lorch: bool = False

def autoscale(q, sq, config, sigma=None) -> ScalingResult
    # ScalingResult: a, b, iterations, per-iteration residuals, C1/C2/C3 diagnostics,
    #                final arrays for all six STOG outputs, provenance dict

def run_stog_compatible(inp_path) -> ScalingResult   # drives autoscale from a stog.inp
```

Every SPDX header: `AGPL-3.0-or-later`, © 2026 Tsung-Han Yang.

---

## 3. Validation plan (uses `data/stog_tests/stog_59438/`)

### 3.1 pystog cross-validation — already executed (2026-07-17)

pystog 0.6.7 (cloned from [github.com/neutrons/pystog](https://github.com/neutrons/pystog))
was run end-to-end against the example before any of our code exists
(script preserved at `data/stog_tests/crosscheck_pystog.py`; needs a venv with
`numpy h5py pystog`):

| Stage | Comparison | Result |
| --- | --- | --- |
| Scaling | `a·S + b` vs `scale.fq` | rms 1.4e−13 (exact) |
| Forward FT | pystog `g(r)−1` vs `scale.gr` | rel. rms 2.5e−3 of full scale → convention pinned |
| Fourier filter | pystog filter section vs `ft.dat` | rms 6.1e−4 |
| Filtered S(Q) | pystog vs `scale_ft.sq` | rms 6.1e−4 (rel. 3e−6) |
| RMC outputs | pystog conversions vs `scale_ft_rmc.*` | 2–7% — explained: classic stog's final low-r enforcement + first-peak cleanup, absent in pystog (§1.5) |

Conclusions: **pystog is a valid oracle through the filter stage** (agreement at
discretization level); the trapezoid sine-kernel discretization is compatible with the
Fortran's; the residual rmc-stage differences are the deliberately separate enforcement stage,
not transform error. The FT-discretization risk in §5 is retired.

Additional finding that reshapes V3: the hand-tuned (a = 10, b = −9) satisfies C1 well
(filtered S(Q) tail mean 1.0011) but leaves a substantial density-limit violation below the
first peak (g−1 mean ≈ +1.1 instead of −1 in [1.2, 2.2] — partly the omitted-low-Q effect,
§1.5). Classic stog then papers over it via enforcement. Therefore V3 must **not** assert
exact recovery of the hand values; it asserts the auto-fit achieves a **low-r residual ≤ the
hand-tuned one** (and reports the recovered (a, b) for human comparison).

### 3.1b Criterion study on the example (Phase-1 implementation, 2026-07-17)

Probing scale-determination criteria against the Fortran run's own filtered outputs at the
expert's scaling settled the objective choice:

| Criterion at (a=10, b=−9), filtered | Value | Target | Verdict |
| --- | --- | --- | --- |
| C1 high-Q tail mean | 0.9996 | 1 | satisfied |
| Low-r window mean (g−1, [1.2, 2.23]) | +0.77 | 0 | violated |
| Krogh-Moe sum rule ∫Q²[S−1]dQ (+S≈0 hole model) | −32.7 | −1.24 | violated ~26× |

The dataset (Qmin = 1.0 Å⁻¹, crystalline) is missing O(1) structure below Qmin, so **no
affine (a, b) can satisfy the absolute normalization criteria** — the expert's 10× encoded
external knowledge, not a self-consistency condition. Consequences implemented:

- The auto-fit objective is C1 + **pointwise** low-r g-residuals (stable on ripple-heavy
  data; binned "level matching" is available via ``c2_bins`` but proved pathological here —
  the mean level is degenerate with the low-Q-hole artifact and lets the solve run away).
- ``diagnostics_summary`` reports ``g_window_mean`` and an ``absolute_scale_reliable`` flag:
  a window mean far from 0 with a good tail means the absolute scale needs external input
  (then use ``scale_pipeline`` with a fixed a).
- The **omitted-low-Q correction defaults on** — on synthetic data with Qmin = 0.6 it cuts
  the recovered-scale bias from 8% to 0.3%. (Derivation note: pystog's linear-F correction is
  algebraically identical to assuming the solid-state limit S(0) = 0 with a linear join —
  the physically right default for dense materials.)
- Also discovered en route: the classic workflow's Fourier filter performs part of the
  normalization — the (10, −9)-scaled *unfiltered* tail sits at ~0.94 and the filter's
  correction lifts it to 1.000 — so C1 must always be evaluated on *filtered* data (the
  loop's ``delta`` term does this).

### 3.2 Test matrix

| # | Test | Pass criterion |
| --- | --- | --- |
| V1 | **Manual-path reproduction**: run the engine with the colleague's exact `stog.inp` parameters (a = 10, b = −9 fixed; no auto-fit) | reproduce `scale.fq` to ~machine precision (formula already verified exact); reproduce `ft.dat`, `scale_ft.sq/.gr`, `scale_ft_rmc.{fq,gr,dr}` within a stated rms tolerance (start loose, tighten as grid/endpoint conventions are matched — `ft.dat` is the oracle that pins the filter discretization) |
| V2 | **Limit assertions** on our outputs | rmc.gr low-r level = −0.015407 ± 1e−6; rmc.dr low-r slope = −1.2207e−2 ± 1e−6 Å⁻¹; filtered S(Q) tail mean = 1 ± tol; rmc.fq tail → 0 |
| V3 | **Blind auto-scale** (headline test): feed the raw rebinned S(Q); minimizer finds (a, b) with no hints | auto-fit low-r residual ≤ the hand-tuned (a=10, b=−9) residual; recovered (a, b) reported and expected in the same neighborhood (see §3.1 — exact recovery is *not* the criterion) |
| V4 | **Synthetic round-trips**: synthetic g(r) (Gaussian first shell) → forward FT with known ρ0 → corrupt with known (a_true, b_true) + noise → recover | (a, b) recovered to tolerance vs noise level; FT forward∘inverse = identity to round-off; Lorch on/off; monatomic sanity: ⟨b²⟩ = ⟨b⟩² |
| V5 | **pystog cross-validation** (optional, `[scaling-dev]` extra): same transforms through pystog | agreement within tolerance; not a hard dependency |

pystog specifics for V5 (verified against its source, v0.6.7): `Converter`/`Transformer`
formulas map 1:1 onto ours (`GK(r)` = Keen G(r), `FK(Q)` = Keen F(Q), `G(r)` = PDFFIT);
CLI config keys `Files[].{Filename,ReciprocalFunction,Qmin,Qmax,Y:{Scale,Offset}}`,
`RealSpaceFunction`, `NumberDensity`, `"<b_coh>^2"`, `"<b_tot^2>"`, `FourierFilter:{Cutoff}`,
`LorchFlag`. It is GPL-3.0-or-later (AGPLv3-compatible, but copyleft — keep out of base deps)
and pins `numpy>=2.1,<3` + Python ≥3.10 while this repo supports ≥3.9 — another reason it must
stay a test-only extra behind `skipUnless(find_spec("pystog"))`. Do not vendor it.

Fixture policy: **all of `data/` is gitignored**, so the local example can never be a tracked
fixture in place. Tests skip cleanly when `data/stog_tests/` is absent (existing repo pattern);
V4 synthetic tests always run in CI. If the example (or a trimmed subset) is cleared for
publication, it goes under `tests/fixtures/stog_59438/` — open decision (§6).

---

## 4. Delivery phases

- **Phase 0 — fixtures + docs**: wire the example into tests (skip-if-absent); this plan doc;
  CHANGELOG entry. **DONE 2026-07-17.**
- **Phase 1 — engine (the core)**: `transforms.py` + `scaling.py` + parser extensions + full
  test suite V1–V4. Pure numpy/scipy; **no new required dependencies**; pystog only as an
  optional dev extra for V5. **DONE 2026-07-17** — 30 new tests, 95-test suite green; see
  §3.1b for the criterion study that shaped the objective.
- **Phase 2 — CLI**: `python -m rmc_toolkits.scaling <stog.inp|args>` + a `[project.scripts]`
  console entry (`rmc-autoscale`). Drop-in replacement for the interactive stog session:
  reads the same `stog.inp`, writes the same six output files + a provenance JSON.
- **Phase 3 — Flask API**: `/api/scaling/preview` (arrays for the live view) and
  `/api/scaling/run` (write outputs), following the existing `_resolve_inside_root` guard and
  per-(path, mtime, params) LRU-cache patterns.
- **Phase 4 — frontend "Auto StoG" page** (decided 2026-07-17): a new tab named **Auto StoG**,
  placed **first** in the tab row (preprocessing precedes monitoring in the real workflow),
  while the app continues to **default to Dashboard** on load (`activePage` stays
  `'dashboard'`). UX is **automation-first — this is the app's selling point**: the primary
  action is one **Auto-scale** button (the minimizer + self-consistent loop); users never
  have to hand-tune yscale/yoffset like classic stog's `try again` loop. Manual (a, b)
  adjustment survives only as an "Advanced" affordance for experts. Layout: file pick →
  parameter form (ρ0, ⟨b⟩², windows; pre-filled from the `.dat` header when present) →
  Auto-scale → results with theory guide-lines (1, −⟨b²⟩, −⟨b⟩², −4πρ0⟨b⟩²r) and the
  fit-quality readout (pre-enforcement low-r rms, (a, b), iteration trajectory) → export.
  Backend-mode first; outputs then flow into the existing dashboard via the already-present
  `"stog"` plot kind.
- **Phase 5 (stretch) — static-mode JS port**: the engine is small (sine transforms + a 2×2
  linear solve + the filter, ~300 lines); a Web-Worker port would preserve the browser-first
  privacy guarantee. Decide after Phase 4 usage; until then the page is backend-only (the
  first such feature — precedent decision flagged in §6).
- **Docs on completion**: revise ROADMAP Phase 7 "Deferred preprocessing" → active
  preprocessing module (its stated precondition — "a dedicated preprocessing module with a
  clear user path" — is exactly what this builds); update the AGENTS.md architecture map;
  CHANGELOG.

Phase 1 is the deliverable that meets the current goal ("all math, testable, validated on the
example"); everything after is surfacing.

---

## 5. Risks

- ~~FT discretization~~ **retired** (§3.1): pystog's trapezoid sine kernel reproduces the
  Fortran intermediates to 6e−4 rms; conventions pinned empirically.
- ~~Second `yoffset` semantics~~ **resolved** (Fortran source `stog_new3.f90`): it is the
  global "Add values" knob, `y·(1+vadd)+vadd`, targeting the high-Q → 1 asymptote. Still
  parsed-but-rejected when nonzero (our fitted `b` supersedes it). `Qoffset` remains
  rejected-if-nonzero.
- ~~Final ripple-removal semantics~~ **resolved** (`stog_new3.f90` L398): zero g(r) where
  `r <= rmccut` and outside `[rmcmin, rmcmax]` — implemented exactly as
  `transforms.first_peak_zero()`; for the example's parameters (2.48, 2.65, 3.1) it
  degenerates to the flat −⟨b⟩² replacement (`enforce_low_r`).
- **Windows choice sensitivity**: defaults (tail fraction, low-r fit window) must be reported
  in provenance so results are reproducible; diagnostics flag when C1 and C2 disagree on the
  scale (symptom of background-shape problems a constant b can't fix).
- **ρ₀ provenance**: the low-r target line is −4πρ₀⟨b⟩²r — a wrong `NUMBER_DENSITY::` yields a
  confidently wrong scale. Surface ρ₀ prominently and editably; validate finite/positive; the
  supercell density from the RMC configuration should agree (future cross-check).
- **Free offset can absorb genuine normalization error**: always report the low-r rms alongside
  (a, b) so a large |b| with a good fit is visible as a warning sign; `fit_offset=False` is one
  flag away.
- **Output clobbering**: the CLI writes canonical stog names — default to a safe `--out-stem`
  and a no-clobber guard so it never silently overwrites a user's real STOG outputs.
- **Scope: neutron / constant-b only in v1.** X-ray data has Q-dependent f(Q) and Compton
  contributions that a constant (a, b) cannot model; label any x-ray use as out of scope until
  a dedicated mode exists.

## 5b. Level-sweep architecture (added 2026-07-17, maintainer's design)

The community question "what Q is high enough for the level?" is answered by a
criterion-driven two-sided window search (`level_sweep`): admissibility = slope
statistically zero given the window's own fit noise; optimality = minimum level variance;
uncertainty = level spread across admissible windows; end artifacts self-exclude
(independently reproduced both experts' hand cuts on the two validation datasets).
`autoscale` defaults to the sweep-anchored 1-dof architecture (`b = 1 − a·level`), with the
Q→0 Faber-Ziman limit as a second, independent amplitude criterion and a concordance
diagnostic. Next refinement: run the sweep on the *uncropped* data so it can inform the
Qmax choice itself (currently it operates within the configured [Qmin, Qmax]).

## 6. Open decisions (for Tsung-Han)

1. **Commit the 59438 example (or a trimmed subset) as a CI fixture?** ~1.7 MB total; it would
   live under `tests/fixtures/` (all of `data/` is gitignored). Provenance/permission to
   publish in an AGPL repo is your call.
2. **Offset default**: fit both (a, b) by default (recommended, reported with a warning when
   |b| is large), or scale-only with offset opt-in?
3. ~~⟨b⟩² / ⟨b²⟩ source~~ **RESOLVED (2026-07-17)**: `rmc_toolkits/scattering.py` computes
   both from a chemical formula (NIST/Sears table, 89 elements, formula parser, isotopic
   overrides, barn + fm² units). `ScalingConfig(b_avg_sq=faber_ziman("SrTiO3").b_avg_sq_barn,
   b_sq_avg=...b_sq_avg_barn)`.
4. **Static-mode stance** for the future page: backend-only (simplest) vs the Phase-5 JS port
   (preserves the hosted-app privacy story).

---

## 7. External verification addendum (multi-agent research, 2026-07-17)

An 11-agent research/design workflow (pystog source + docs, STOG lineage, licensing, three
independently-designed architectures, adversarial fact-checks, judge synthesis) confirmed the
plan's load-bearing assumptions:

- **pystog has no auto-scaler** (verified in `src/pystog/stog.py`): every `Y:{Scale,Offset}`
  is user-supplied (defaults 1.0/0.0). Its only low-r cost-function helpers
  (`_lowR_mean_square`) are dead code marked *"Currently not used in PyStoG workflow since was
  done manually"* with a TODO to automate — i.e. ORNL themselves flagged the gap this module
  fills. The auto-fit is new logic regardless of engine choice.
- **Judge recommendation** matches this plan's phasing: library + CLI first (scored 8.5/10),
  then the backend page (7.5/10), static-mode Web-Worker parity last (7/10) — one pure-numpy
  source of truth, pystog test-only.
- **Licensing**: pystog is GPL-3.0-or-later — compatible with this repo's AGPL-3.0-or-later,
  but kept out of runtime deps (also incompatible version pins; see §3).
- **Naming — resolved by maintainer (2026-07-17)**: the user-facing page is **"Auto StoG"**
  (deliberate nod to the classic StoG program, with "Auto" signaling the differentiator);
  the Python module stays neutral (`rmc_toolkits/scaling.py`). The app itself was renamed
  **RMCProfile Workbench** (was "Run Monitor") to reflect the multi-tool scope.

A second workflow (math verification against Keen 2001 / RMCProfile manual / pystog
transformer internals) is folded into §1 as its findings land; §1's ✓data checks are already
independently exact.
