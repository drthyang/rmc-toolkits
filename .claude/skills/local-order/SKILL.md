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

**Status: Interviewed 2026-10-24.** The steps are the user's own answers of that day
(docs/SKILL_INTERVIEWS.md §5).

Each kind of local order hides from every analysis but one or two. The average structure folds the
supercell into one cell; each site's displacement cloud is centred on its own mean; an RMC
configuration that the data did not shape looks like independent atoms. So one empty analysis
says nothing about the others: "no order" is the conclusion that every test the Workbench has came
back empty. What was not run is "not tested", never "absent".

## Steps

1. **Sites:** `list_sites`. Note the split sites (two positions, by the app's rule: a resolved
   PC1 whose excess kurtosis is ≤ −0.5 and at least 3 sampling errors below zero; `list_sites`
   gives them a `delta_A`), the other sites with negative kurtosis (flat-topped), large positive
   kurtosis (heavy tails, an anharmonic well) or large anisotropy.
2. **Directions, per site:** `analyze_displacement_directions` with `site_ref` for every
   symmetry-inequivalent site (one per orbit), not only those whose kurtosis or anisotropy stands
   out; never pooled over sites whose local frames differ. Empty means no peak at ≥ 3σ. A
   significant map or tensor test alone is anisotropy, which thermal ellipsoids have: it is not
   order.
3. **Order between cells:** `analyze_modulations` (every element). Empty means no zone-boundary
   point and no general wavevector called ordered, and no diffuse short-range order: the diffuse
   verdict, when it passes its calibrated bar, counts as order between cells. A long wave (every
   |q_i| ≤ 0.25) is a box-scale correlation, not order. A split site whose element has an ordered
   point is order with that q (a superstructure), not random disorder.
4. **Clusters:** for every element whose own first shell is closed and holds 0.5–4.5 neighbours of
   its own kind (`find_gr_peaks` on the A–A pair: a dimer 1, a trimer 2, a tetrahedron 3, an
   octahedron 4 per member), `analyze_clusters` in its `first_shell.window_A`.
   Empty means no shape distortion or pattern beyond independent atoms at the clusters' own size,
   no breathing, no rigidity, and a regular average cluster. An extended network (the tool says
   so) has no clusters to test: test its polyhedra with `polyhedral_distortion` (Baur's DI and
   Robinson's indices, against the average structure's polyhedra).
5. **Polyhedra:** `find_gr_peaks` for each bonded pair: the whole first shell's CN and peaks against
   the average structure. A split shell the average structure does not explain is a distorted
   polyhedron (coordination-shells); a split explained by a split centre site is that site's
   disorder (step 1–3). A partial under 5 % of the signal is "not constrained by the data": say
   so, and leave it out of the conclusion.
6. **Average structure:** the context's ladder and noise floor. Empty means the group of the lowest
   rung clearly above the noise floor (from 1.25 × `noise_floor_A`) is the card's group, or
   `measure_symmetry_lowering` finds no site shifted above the noise floor (local-symmetry).
7. **Conclude.** "No order beyond thermal motion" only when steps 2–6 all came back empty; name
   each analysis and its result, and the partials left out as not constrained by the data. When a
   test could not run (too few points per sphere cell, a zone-boundary point an odd supercell axis
   cannot hold), say "no order found in the tests that could run; ⟨test⟩ not tested because
   ⟨reason⟩": never a bare "no order". When one found order, report that order: the analysis, the
   statistic, its effect size and σ.

## Acceptance checklist

- [ ] Every analysis of steps 2–6 ran, or the answer names what was not tested and why.
- [ ] Directions for every symmetry-inequivalent site, not pooled; anisotropy not called order.
- [ ] Split sites read against `analyze_modulations` before they are called disorder.
- [ ] Clusters tested for every element with a finite same-element shell; an extended network's
      polyhedra with `polyhedral_distortion`.
- [ ] Order between cells at the zone boundary, at general wavevectors and in the diffuse
      verdict; no long wave called order.
- [ ] The average structure's group read above the noise floor.
- [ ] A partial the data barely weight (step 5) reported as "not constrained by the data" and
      left out of the conclusion.
- [ ] Never a bare "no order": a test that could not run is named with its reason.
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
- **"No order" before the whole list ran.** *Firm rule in the app:* an answer that concludes "no
  order" is refused until every applicable analysis of steps 2–6 has run, the diffuse verdict and
  `polyhedral_distortion` included (the evidence gate, extended).

## Report

Each analysis with its result (statistic, effect size, σ); the order found, or "none at ≥ 3σ";
what was not tested and why ("⟨test⟩ not tested because ⟨reason⟩"); the partials not constrained
by the data; the average structure's group above the noise floor.

## References

- [order-signatures](order-signatures.md) — how each kind of order shows, and hides, in each analysis, with the dev cases' numbers.
