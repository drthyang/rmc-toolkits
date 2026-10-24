---
name: displacement-anisotropy
description: >-
  Judge a site's displacements honestly: thermal ellipsoids from the PCA of its displacement
  cloud (U_iso, principal amplitudes, anisotropy), Mardia kurtosis, preferred displacement
  directions and their calibrated significance (σ), +u/−u asymmetry, and when an anisotropy or a
  direction is real rather than sampling noise. Use for "most anharmonic / most disordered
  site?", "preferred direction?", "is this site anisotropic?", "trigonal or tetragonal
  distortion?". NOT for bond windows (coordination-shells) or space groups (local-symmetry).
---

# Displacement anisotropy and directions

**Status: Draft: awaiting the user's interview.** Written from the 2026-10-10 spot check, the
review that followed, the PCA Ellipsoid and Displacement Directions algorithm notes, and the
tools' checks; the split-site step, the 15° family rule and the peak-test rule were added in
validation batch 1 (the 2026-10-14 tool review, approved by the user that day). A step marked "(to
confirm with the user)" is not yet the user's own practice: say in the answer which choice you
made.

Each atom's displacement is its offset from its own reference site's mean position, so a site's
cloud is centred by construction: a coherent off-centring shared by every copy moves the mean
and is invisible here.

## Steps

1. **Rank the sites:** `list_sites` (optionally for one element). Per site: U_iso (Å²), the
   principal RMS amplitudes, anisotropy = rms1/rms3, and Mardia's multivariate excess kurtosis
   (0 for a Gaussian; > 0 heavy-tailed or peaked; < 0 flat-topped). A site is split (two wells)
   only by the app's rule, not by a negative Mardia kurtosis alone: a resolved PC1 whose excess
   kurtosis is ≤ −0.5 and at least 3 sampling errors below zero (`list_sites` then gives its static
   displacement `delta_A`). Split sites come first (largest `delta_A`), then by |kurtosis|; the
   `largest_*_refs` name the extremes.
   **A split site is not yet disorder.** Ordered antiphase displacements (a superstructure) give a
   site exactly the same cloud as random two-site disorder: only the correlation between cells
   differs. Before calling a split site (`delta_A`) or a flat one (negative kurtosis) random
   off-centring, run `analyze_modulations`: an ordered zone-boundary point (or order at a general
   wavevector) of that element makes it order with that q, not disorder.
2. **One site at a time:** `analyze_displacement_directions` with `site_ref`. Pooling an element
   over sites mixes their local frames (each vertex of an M₄ tetrahedron sits on its own 3m
   axis), so a preference every site shares in its own frame cancels in the pooled map: an
   isotropic pooled map is not evidence against a per-site or cluster distortion.
3. **Read the significances, all calibrated one-sided σ** (≥ 3 is significant):
   `peak.significance_sigma` (the peak cell, corrected over all cells), the whole map's departure
   from isotropy, `orientation_tensor.anisotropy_significance_sigma` (Bingham), and the
   +u/−u asymmetry z (a skewness the ellipsoid cannot show). A local z-score is never a
   significance. **A preferred direction needs the peak test:** `peak.significance_sigma` ≥ 3σ,
   unrounded (2.96σ is never rounded to 3σ). The map and the tensor tests show anisotropy, which
   every anisotropic thermal ellipsoid has (an O flattened across its bond, thermal motion only:
   map 7.2σ, peak 2.3σ): say "anisotropic", not "a preferred direction".
4. **Directions:** `peak.direction_uvw` and the tensor's principal axis are direct-lattice
   `[u v w]` directions (not plane normals); the map is not antipodally folded, so +u and −u
   differ. The peak names its direction only to within its cell's radius (`peak.radius_deg`).
   Name the nearest low-index family (<100>, <110>, <111>) only when the angle to it is within
   15° (the tools' cone, 15° included: under half the 35.3° between <110> and <111>, the closest
   families) and the cell is no wider (`peak.radius_deg` ≤ 15°; the tool names no family above
   it); otherwise give the `[u v w]` as measured.
5. **Size, not only significance.** With thousands of atoms a tiny effect reaches 3σ: quote the
   enhancement (× isotropic) or the anisotropy with its σ. For a cluster-level question
   (trigonal versus tetragonal distortion of an M₄ cluster), use `analyze_clusters`: its <100>,
   <110>, <111> gyration-axis families against the shuffle null, with effect size and z. The
   shape is read at each cluster's own size, and a collective size fluctuation (breathing) is
   reported on its own: it is not a shape distortion. Where the member sites are split, the null
   already holds the distortion's amplitude, and only the pattern (axis families, trigonality)
   can show it.

## Judgment calls

- The sphere resolution (the tool, and the page, which opens at Auto, take the engine's choice
  from the point count) and whether a result must hold at two resolutions (to confirm with the
  user).
- Whether a 3–5σ direction with a small enhancement (under 1.2× isotropic) is reported as a
  finding or as a hint (to confirm with the user).
- When to report per-axis kurtosis (only on axes the PCA resolves) (to confirm with the user).
- The 15° cone for naming a low-index family (batch 1; to confirm with the user).

## Acceptance checklist

- [ ] Sites are named by `ref` and element, with U_iso, anisotropy and kurtosis from the tool.
- [ ] Directions are analysed per site, or pooling is stated with its caveat.
- [ ] "Preferred direction" only with a peak significance ≥ 3σ, quoted with the enhancement; a
      significant map or tensor test alone is called anisotropy.
- [ ] `[u v w]` given as a direction, the low-index family named only within 15°.
- [ ] A split site is called disorder only after `analyze_modulations` found no order of its
      element.
- [ ] Cluster-level claims rest on `analyze_clusters` (effect size and z).
- [ ] No invented uncertainty: every number from a tool.

## Known failure modes

- **The pooled map read as "no distortion"** (review: 16 Ta sites pooled → isotropic). The
  tool now notes the cancellation; per site and with `analyze_clusters` the Ta₄ clusters prefer
  <111> compression (19 % vs 14 %, z = 11), with no <001> signal.
- **Anisotropic thermal motion read as a preferred direction** (tool review, 2026-10-14, L3): the
  direction check counted the map test, which any thermal ellipsoid passes. *Firm rule in the
  app:* only the peak test makes a direction significant; the map test is reported as
  anisotropy.
- **Ordered antiphase displacements read as random off-centring** (validation arm A, the
  zone-boundary dev case): the PCA cloud is split exactly as for two-site disorder;
  `analyze_modulations` finds the order at (0, 0, ½) at > 8σ.
- **A direction called preferred without significance.** *Firm rule in the app:* an answer that
  claims a preferred displacement direction when every directions analysis of the question found
  none at ≥ 3σ gets a warning under it.
- **Effect size ignored:** thousands of clusters make a 1 % effect "significant"; the cluster
  tool requires both its z and its effect size to clear their thresholds.
- **Numbers differ from the page** only when the user picks another resolution there: the page
  opens at Auto, the engine's resolution and the tool's, and reads its numbers from the raw
  histogram (smoothing only redraws the map); say which you quote.

## Report

The sites in question with U_iso, anisotropy and kurtosis; per site, the strongest direction with
its enhancement and σ (or "none at ≥ 3σ"); the asymmetry where significant; for clusters, the
analyze_clusters result with effect size and z; and what the average ellipsoid cannot show.

## References

- [significance-reading](significance-reading.md) — what each σ in the direction and cluster results tests, with the Demo's values.

Further reading for Claude Code (not loaded in the app): `docs/algorithms/pca-ellipsoid.md`,
`docs/algorithms/displacement-directions.md`.
