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
| Low-r enforcement | ON at the foot of the detected first shell — min(foot, onset − 0.25 Å), below its rising flank — (stog.inp / explicit cutoffs when present) | classic-product parity without removing first-shell signal; the *pre*-enforcement residual is always reported |

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
5. **r₀ detection**: the *first* coordination shell — the smallest-r |g| feature that
   stands out of the ripple field below it (≥ 4× that field, or ≥ 2× while ≥ 50 % of the
   range maximum), not the tallest one — is located and its left flank (35% of its own
   height) taken as the data's first-shell onset. |g| because negative-b pairs give
   *inverted* shells (Ti–O in titanates, Mn–Sn in Mn₃Sn), often weaker than the second
   shell. Without a given r₀ the window is located, not assumed: two trial fits
   ([r_cut+0.2, +0.3] and [.., +1.0] Å) propose candidate onsets whatever the sign of
   their scale; the smallest is refitted on [r_cut+0.2, onset − 0.25] and must be
   re-detected (within 0.15 Å) on that refit's own g(r) — a candidate the refit no longer
   shows was a ripple and is dropped, a lower shell the refit uncovers is tried next, and a
   refit with a ≤ 0 stops the run (a fit with a ≤ 0 is never returned). If no shell is
   confirmed, or a confirmed one leaves < 0.1 Å of window (bonds shorter than ~1.75 Å at
   the default r_cut = 1.0: Si–O, P–O, B–O, C–O), the run stops with the r_cut to use
   instead of fitting across the shell. Confirmed onsets, composition-only, over Qmin 0.82
   and 1.0 × Qmax 24–30 (56 configurations): 2.65–2.75 Å on the four Mn₃Sn runs (window top
   2.40–2.50 Å, a > 0 in every returned fit). Nine configurations stop instead — PG3_55537 at
   (Qmin, Qmax) = (0.82, 24/25/28) and (1.0, 24–27) and 59438 at (1.0, 28) with "could not
   locate the first coordination shell" (the degenerate density limit leaves the inverted
   Mn–Sn shell just under the detector's margin), 300 K at (1.0, 25) with a non-physical
   refit scale — so quote an onset only with its Q range, and for such data use
   `--amplitude fz` or pin r₀. 2.52–2.53 Å for FeCoSn 199 K (Qmin 0.5/1.0 × Qmax 22–26).
   These are flank points of the first peak, i.e. *above* the hand-chosen classic cutoffs
   2.40–2.68 Å, which sit below it.
6. **Independent cross-check**: `a_fz` from the Q→0 Faber-Ziman limit (level-subtracted
   head extrapolated to S(0)). Concordance `a_fz/a ≈ 1` is the absolute-scale trust
   metric; discord quantifies what the data cannot decide (and flags a wrong ρ₀ ~1:1).
7. **Outputs**: scaled S(Q), unfiltered g(r) (`scale.gr`), filtered S(Q) and g(r) with
   r·[g(r)−1] (`scale_ft.sq`, `scale_ft.gr` — the Fortran stog conventions), and the
   RMCProfile-ready `F_K(Q)`, `G_K(r)`, `D(r)` — with classic low-r enforcement applied
   below the first shell (automatic cutoff = min(foot, onset − 0.25 Å): 2.43 Å on Mn₃Sn
   59438 at Qmin 1.0 / Qmax 27 vs the expert's 2.48 (the expert's Qmax 28 stops, above),
   2.40–2.50 Å on the Mn₃Sn configurations that fit over Qmin 0.82/1.0 × Qmax 24–30,
   2.27–2.28 Å on FeCoSn 199 K; the first-shell coordination number is preserved to
   ≤ 0.3 %; flags and pre-enforcement residuals reported) — plus a provenance JSON.

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
  **composition is the scale information**: `--amplitude fz` (the level-subtract →
  pin S(0) → restore-level construction) gives a = 9–16 on the 55537/55526/54139 runs —
  but only when its Q→0 extrapolation is well conditioned. On run 59438 the Bragg-dominated
  head extrapolates to within noise of the level: a_fz = 54 (Qmin 0.82), 76 (1.0), 309 (1.05),
  flagged `a_fz_reliable = False` (relative error 29–145 %). **`a_fz_reliable = True` is
  necessary, not sufficient**: the flag only says S_meas(0) − level is resolved from its
  statistical error, and a systematically biased low-Q head passes it. On two of the three
  runs above the reliable-flagged a_fz still drifts with Qmin — 55537: 9.0 → 4.1–4.7, 54139
  (500 K): 15.7 → 23.3 over Qmin 0.82–1.05 in 0.01 steps (Qmax 28), every value flagged
  reliable (relative error 11–20 % and 10–18 %); 55526 (300 K) stays within 8.9–10.7.
  Before trusting an FZ scale, re-run at a few Qmin values and check a_fz is stable, check the
  concordance with the density-limit amplitude, and compare it with other runs of the
  material / an external density. The CLI and the page print this caveat next to every
  reliable a_fz.
- The RMC-ready files satisfy the Keen limits *by construction* (enforcement); judge fit
  quality only on the reported pre-enforcement numbers.

## 5. Practical guidance

- **Qmax**: end it before detector rolloff — the level sweep will show a shrinking
  admissible window and the fit degrades loudly if rolloff enters (robustness study:
  stable ±3% over Qmax 18–28 on FeCoSn, collapse + flag at 30).
- **Qmin**: as low as the reduction allows; the correction handles the rest. Cutting real
  low-Q information (Qmin ≳ 1.5–2) starves the density limit (flagged).
- **X-ray data**: the Sears table is neutron — set ⟨b⟩² (usually 1 for normalized S(Q))
  and ⟨b²⟩ = ⟨Z²⟩/⟨Z⟩² explicitly (f(0) = Z). A composition given alongside (e.g. for the
  mass-density conversion) never supplies ⟨b²⟩ to a ⟨b⟩² from another source: the pair must
  come from one source, and ⟨b²⟩ < ⟨b⟩² (S(0) > 0) is refused.
- **Isotopic samples**: per-element b overrides are supported in the library
  (`faber_ziman(..., b_overrides_fm=...)`).
- ρ₀ sanity: the implied mass density is shown; Mn₃Sn's 0.063049 atoms/Å³ ↔ 7.42 g/cm³.
