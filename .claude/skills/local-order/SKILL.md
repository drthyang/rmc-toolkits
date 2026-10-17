---
name: local-order
description: >-
  Decide whether a configuration has local order beyond the average structure, and say "no order"
  only when every analysis that could show one came back empty: displacement directions per site,
  split sites against the superstructure test, clusters, order between cells, a lowering of the
  average structure, split coordination shells. Use for "is there local order?", "anything beyond
  thermal motion?". NOT for one site's directions (displacement-anisotropy).
---

# Local order: when to say "no order"

**Status: Draft: awaiting the user's interview.** Written in validation batch 1 (2026-10-17) from
the 2026-10-14 tool review's proposal, approved by the user that day, and the validation's dev
cases (synthetic configurations whose order is known). A step marked "(to confirm with the user)"
is not yet the user's own practice: say in the answer which choice you made.

Each kind of local order hides from every analysis but one or two. The average structure folds the
supercell into one cell; each site's displacement cloud is centred on its own mean; an RMC
configuration that the data did not shape looks like independent atoms. So one empty analysis
says nothing about the others: "no order" is the conclusion that every test the Workbench has came
back empty. What was not run is "not tested", never "absent".

## Steps

1. **Sites:** `list_sites`. Note the sites with negative kurtosis (flat or split: two positions),
   large positive kurtosis (heavy tails, an anharmonic well) or large anisotropy.
2. **Directions, per site:** `analyze_displacement_directions` with `site_ref`, one site per
   symmetry-inequivalent orbit, never pooled over sites whose local frames differ. Empty means no
   peak at ≥ 3σ. A significant map or tensor test alone is anisotropy, which thermal ellipsoids
   have: it is not order.
3. **Order between cells:** `analyze_modulations` (every element). Empty means no zone-boundary
   point and no general wavevector called ordered. A long wave (every |q_i| ≤ 0.25) is a box-scale
   correlation, not order. A split site whose element has an ordered point is order with that q (a
   superstructure), not random disorder.
4. **Clusters:** for every element whose own first shell is short and closed (`find_gr_peaks` on
   the A–A pair: dimers, trimers, tetrahedra), `analyze_clusters` in its `first_shell.window_A`.
   Empty means no shape distortion or pattern beyond independent atoms at the clusters' own size,
   no breathing, no rigidity, and a regular average cluster. An extended network (the tool says
   so) has no clusters to test.
5. **Polyhedra:** `find_gr_peaks` for each bonded pair: the whole first shell's CN and peaks against
   the average structure. A split shell the average structure does not explain is a distorted
   polyhedron (coordination-shells); a split explained by a split centre site is that site's
   disorder (step 1–3).
6. **Average structure:** the context's ladder and noise floor. Empty means the group of the lowest
   rung above the noise floor is the card's group, or `measure_symmetry_lowering` finds no site
   shifted above the noise floor (local-symmetry).
7. **Conclude.** "No order beyond thermal motion" only when steps 2–6 all came back empty; name
   each analysis and its result. When one found order, report that order: the analysis, the
   statistic, its effect size and σ.

## Judgment calls

- Whether every symmetry-inequivalent site needs its directions tested, or only those whose
  kurtosis or anisotropy stands out (to confirm with the user).
- Whether a polyhedron the data barely weight (a partial carrying a few % of the signal) must be
  tested before "no order" (to confirm with the user).
- How to word "no order" when a test is undefined: too few points per sphere cell, or a
  zone-boundary point an odd supercell axis cannot hold (to confirm with the user).

## Acceptance checklist

- [ ] Every analysis of steps 2–6 ran, or the answer names what was not tested.
- [ ] Directions per site, not pooled; anisotropy not called order.
- [ ] Split sites read against `analyze_modulations` before they are called disorder.
- [ ] Clusters tested for every element with a finite same-element shell.
- [ ] Order between cells at the zone boundary and at general wavevectors; no long wave called order.
- [ ] The average structure's group read above the noise floor.
- [ ] No literature expectation used as evidence either way.

## Known failure modes

- **The literature for the data** (validation, the thermal-null dev case): a perovskite whose
  composition is known for off-centred B cations, thermal motion only in the configuration. The
  answer is "no order"; its O ellipsoids are anisotropic (map 7.2σ, peak 2.3σ), which is not order.
- **"No order" from the pages alone** (validation arm A, the breathing dev case): no page shows a
  cluster's breathing; `analyze_clusters` finds it at z = 22.
- **Order called disorder** (validation arm A, the zone-boundary dev case): ordered antiphase Ti
  displacements give the same split cloud as random off-centring; only `analyze_modulations` tells
  them apart.
- *Firm rule in the app:* an answer that claims a superstructure, a cluster distortion, breathing,
  rigidity or a preferred direction that the question's analyses did not find gets a warning under
  it; and an answer that concludes "no order" before the superstructure test, a per-site direction
  map and the cluster analyses ran is sent back to run them (the evidence gate, THY-18).

## Report

Each analysis with its result (statistic, effect size, σ); the order found, or "none at ≥ 3σ";
what was not tested and why; the average structure's group above the noise floor.

## References

- [order-signatures](order-signatures.md) — how each kind of order shows, and hides, in each analysis, with the dev cases' numbers.
