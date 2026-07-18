# Getting S(Q) and G(r) on absolute scale — the Auto StoG procedure

The complete, robust recipe for putting measured total-scattering S(Q) — and everything
derived from it — on absolute (barns) scale, as implemented by Auto StoG
(`rmc_toolkits.scaling`, the `rmc-autoscale` CLI, `/api/scaling/*`, and the Auto StoG tab).
Verified against Keen's convention paper, the classic Fortran `stog_new3`, pystog 0.6.7, and
five complete instrument runs (POWGEN Mn3Sn ×4 incl. run 59438, FeCoSn x-ray ×2). Companion
docs: [STOG_SCALING_PLAN.md](STOG_SCALING_PLAN.md) (verified math + validation record).

**References.** D. A. Keen, *J. Appl. Cryst.* **34**, 172 (2001) — all function conventions
and limits (equation numbers below). [pystog](https://github.com/neutrons/pystog) — the
maintained reimplementation of classic StoG (its operation order is what Auto StoG
automates; it has **no** auto-scaler — every `Y.Scale/Offset` is user-supplied and its
low-r minimizer is a `TODO`). [ADDIE](https://github.com/neutrons/addie) — the ORNL
front-end whose density and Faber-Ziman conventions the inputs follow.

---

## 1. What the user must provide — and what everything else defaults to

| Input | Required? | Default / source |
| --- | --- | --- |
| **Chemical composition** | **yes** | drives ⟨b⟩², ⟨b²⟩ (Sears/NIST bound coherent lengths), the S(0) target, and mass↔number density conversion |
| **Qmin, Qmax** | **yes** | your fit/transform window — the one genuinely experimental choice (see §5) |
| ρ₀ (number density, atoms/Å³) | one of three | `NUMBER_DENSITY ::` data header → **mass density** (g/cm³, converted ADDIE-style: ρ₀ = ρ_m · N_A/10²⁴ · n/M) → explicit value |
| Everything else | no | see table below |

Defaults (all overridable, none normally touched):

| Parameter | Default | Why |
| --- | --- | --- |
| Fourier-filter cutoff | 1.0 Å | classic stog/pystog convention; removes sub-atomic artifacts |
| r grid | 5000 pts to 50 Å | classic stog defaults |
| r₀ (closest approach) | **detected from the data** (first-shell |g| flank, §3 step 5) | places the low-r fit window without a structural prior |
| Low-r fit window | [cutoff + 0.2, r₀ − 0.25] | below the first shell, above the filter |
| High-Q architecture | level sweep (`b = 1 − a·level`) | measures "what Q is flat enough" statistically |
| Amplitude criterion | density limit, FZ as cross-check | switchable to `fz` when the density limit is degenerate (§4) |
| Low-Q correction | ON, extrapolating to the **composition-derived S(0)** | §3 step 2 — this is load-bearing (53% scale bias without it on real data) |
| Robust re-weighting | Huber IRLS ON | isolated Bragg/ripple outliers cannot drag the fit |
| Lorch window | OFF | resolution first; turn on for ripple-heavy display |
| Low-r enforcement | ON at the detected r₀ (stog.inp cutoffs when present) | classic-product parity; the *pre*-enforcement residual is always reported |

## 2. The functions (Keen 2001 conventions)

- `S(Q)` — normalized structure factor, → 1 at high Q (Eq. 19/21).
- `F(Q) = Q[S(Q) − 1]`; Keen's barns-scale `F_K(Q) = ⟨b⟩²[S(Q) − 1]` (Eq. 9/19).
- `g(r)` — pair distribution function, → 1 at large r, ≡ 0 below the closest approach r₀.
- `G_K(r) = ⟨b⟩²[g(r) − 1]` (Eq. 10/16): flat **−⟨b⟩²** below r₀ (Eq. 15).
- `D(r) = 4πρ₀ r G_K(r)` (Eq. 29): straight line of slope **−4πρ₀⟨b⟩²** below r₀.
- Limits: `S(∞) = 1` (Eq. 21); `S(0) = 1 − ⟨b²⟩/⟨b⟩²` (Eq. 21, compressibility term
  ignorable for dense solids); `F_K(0) = −⟨b²⟩` (Eq. 14).
- ⟨b⟩² = (Σ cᵢbᵢ)² is the stog "Faber-Ziman coefficient" (barns); ⟨b²⟩ = Σ cᵢbᵢ² is a
  *different* number. pystog quotes fm² in places (1 barn = 100 fm²); Auto StoG computes
  both in both unit systems from the composition.

## 3. The pipeline (what one Auto-scale press runs)

Classic order (read → merge → scale → transform → filter → Lorch → Keen conversions), with
the manual "try again" scale loop replaced by physics:

1. **Read + crop** the S(Q) file (count headers, NaN padding, σ column tolerated) to
   (0, Qmin…Qmax]; optional despiking for detector glitches.
2. **Composition constants**: ⟨b⟩², ⟨b²⟩, and the Q→0 target `S(0) = 1 − ⟨b²⟩/⟨b⟩²`. The
   analytic correction for the unmeasured [0, Qmin] range extrapolates S(Q) linearly to
   *that* S(0) — pystog's correction is the special case S(0) = 0, which is badly wrong for
   negative-b compositions (Mn₃Sn: S(0) = **−12.06**; using the composition-aware target
   cut the low-r residual ~40% on the PG3 runs).
3. **Level sweep**: every candidate high-Q window is line-fitted in O(1); windows whose
   slope is statistically zero are admissible; the minimum-variance one defines the
   measured level L ± its honest spread. The offset is anchored: **b = 1 − a·L** (your
   colleagues' hand scalings all satisfied exactly this with L ≈ 1).
4. **Amplitude a** from the low-r density limit — g(r) → 0 below r₀ — solved in closed
   form (the sine transform is affine in (a, b)), inside a self-consistent loop with the
   **Fourier filter** (r < cutoff content removed and re-transformed; ft.dat is that
   correction). Converges in ~3–7 iterations.
5. **r₀ detection**: the dominant first-shell |g| feature is located and its left flank
   (35% of peak height) taken as the data's closest approach — |g| because negative-b
   compositions can have an *inverted* first shell — then the fit window is refined and
   the fit re-run once. Detected values: 2.73–2.77 Å across the Mn₃Sn runs, 2.53 Å for
   FeCoSn (matching the hand-chosen classic cutoffs 2.40–2.68).
6. **Independent cross-check**: `a_fz` from the Q→0 Faber-Ziman limit (level-subtracted
   head extrapolated to S(0)). Concordance `a_fz/a ≈ 1` is the absolute-scale trust
   metric; discord quantifies what the data cannot decide (and flags a wrong ρ₀ ~1:1).
7. **Outputs**: scaled S(Q), unfiltered g−1, filtered S(Q)/g−1(+D), and the
   RMCProfile-ready `F_K(Q)`, `G_K(r)`, `D(r)` — with classic low-r enforcement applied at
   r₀ (flags and pre-enforcement residuals reported) — plus a provenance JSON.

## 4. Reading the verdicts — when is the scale actually absolute?

- **`density_limit_satisfied` is one-sided.** False *proves* self-consistency cannot fix
  the absolute scale on this data (missing low-Q information); True is necessary, not
  sufficient — a smooth low-Q deficiency is silently absorbed into a biased scale.
- **Concordance is the trust metric.** The density-limit and FZ amplitudes share nothing
  but the data; agreement (FeCoSn: 4–6%) is strong evidence, disagreement is an alarm
  (wrong ρ₀ moves only the density amplitude; missing low-Q moves them apart).
- **Negative-b / near-null-matrix compositions (the Mn₃Sn case).** ⟨b²⟩/⟨b⟩² ≫ 1 means
  S(Q) carries an O(⟨b²⟩/⟨b⟩²) dive to S(0) that data starting at Qmin ≈ 0.8 never see:
  the density limit is degenerate (flag False on every PG3 run), and the historical hand
  scalings are mutually inconsistent (×2.5, ×2.05, ×10 for the same material). Here the
  **composition is the scale information**: use `--amplitude fz` (the level-subtract →
  pin S(0) → restore-level construction), and treat the result as the defensible one.
- The RMC-ready files satisfy the Keen limits *by construction* (enforcement); judge fit
  quality only on the reported pre-enforcement numbers.

## 5. Practical guidance

- **Qmax**: end it before detector rolloff — the level sweep will show a shrinking
  admissible window and the fit degrades loudly if rolloff enters (robustness study:
  stable ±3% over Qmax 18–28 on FeCoSn, collapse + flag at 30).
- **Qmin**: as low as the reduction allows; the correction handles the rest. Cutting real
  low-Q information (Qmin ≳ 1.5–2) starves the density limit (flagged).
- **X-ray data**: the Sears table is neutron — set ⟨b⟩² (usually 1 for normalized S(Q))
  and ⟨b²⟩ = ⟨Z²⟩/⟨Z⟩² explicitly (f(0) = Z).
- **Isotopic samples**: per-element b overrides are supported in the library
  (`faber_ziman(..., b_overrides_fm=...)`).
- ρ₀ sanity: the implied mass density is shown; Mn₃Sn's 0.063049 atoms/Å³ ↔ 7.42 g/cm³.
