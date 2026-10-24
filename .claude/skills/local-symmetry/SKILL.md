---
name: local-symmetry
description: >-
  Read the average structure's symmetry correctly: the space group at a given tolerance
  (detect_symmetry and the context's tolerance ladder), what changes from one rung to the next,
  which changes lie below the sampling noise, and the facts that follow from a point group
  (inversion centre, chirality, polarity). Use for "what is the space group?", "what breaks at a
  tighter tolerance?", "is it still cubic?", "is it chiral / centrosymmetric / polar?". NOT for
  bond windows (coordination-shells) or per-site displacements (displacement-anisotropy).
---

# Local symmetry: the tolerance ladder

**Status: Interviewed 2026-10-24.** The steps are the user's own answers of that day
(docs/SKILL_INTERVIEWS.md §3).

The finder analyses the **average** structure: the configuration folded into one cell, each
reference site at its mean position. It says nothing directly about local distortions, which
average out: those are the displacement and cluster analyses' job.

## Steps

1. **Read the ladder** in the context (`symmetry.ladder`): each rung is the group that holds
   from one tolerance (Å, the largest position mismatch allowed) to the next. Loose tolerances
   give high symmetry; tighter ones lower it. Note `noise_floor_A`: the site means carry
   sampling noise (rms displacement / √copies), so a change of group below that tolerance is
   noise, not a distortion.
2. **At a tolerance the ladder does not show**, call `detect_symmetry` with `tolerance_A`. It
   reports the group, its number, its point group with the point group's facts (centrosymmetric,
   Sohncke, polar, crystal system), the operations and the Wyckoff orbits.
3. **What breaks between two groups.** Describe a lowering as the operations lost (a three-fold
   axis means cubic is lost; a lattice centring halves or quarters the cell's translations), plus
   the group–subgroup index and its type, t (translationengleiche) or k (klassengleiche); give
   irreps only on request. No Workbench tool states the index or the type yet: work them out from
   the two groups' operations, which `detect_symmetry` lists. Take inversion, chirality and
   polarity from the point-group facts, never from memory: a change that lowers the symmetry
   cannot "lose" an inversion centre the higher group never had.
4. **Is the change real?** A rung below `noise_floor_A` is within the noise of the site means. A
   rung at or above 1.25 × `noise_floor_A` (the app's margin) is a distortion of the average
   structure; say how large (the tolerance at which it appears). A rung between the floor and
   1.25 × it is "at the limit of the sampling: neither a distortion nor ruled out", never a
   lowering. Against a published low-temperature phase whose doubled cell the box can't hold,
   compare only the part the box can represent (the t-subgroup or the point group), and say that
   the box can neither confirm nor rule out the cell-doubled phase.
5. **The group to report.** "The" group is the lowest rung at or above 1.25 × `noise_floor_A`,
   with the Detected SG card's 0.2 Å group beside it when they differ: the card's default can
   hide a real lowering (a 0.07 Å polar shift of every Ta reads Pm-3m at 0.2 Å, P4mm above the
   0.026 Å floor). When they differ, size the lowering with `measure_symmetry_lowering`: which
   sites move from the parent positions, by how much (Å), along which direction, against the
   noise floor. That shift is shared by every copy of the site, so the displacement analyses
   cannot see it. Name a group as "consistent with ⟨group⟩ at ⟨tolerance⟩ Å", never "the
   structure is ⟨group⟩"; with no space-group number, give the point group and the crystal class.
6. **Local versus average.** A cubic average with locally distorted clusters is possible and
   common: for local distortion, use `analyze_clusters` and the displacement analyses, not the
   ladder.
7. **A lower-symmetry CIF** (the Detected SG card's export) only for a rung at or above 1.25 ×
   `noise_floor_A`. First check that the model reproduces the site means within the rung's
   tolerance, and that its Wyckoff sites and occupancies make sense. Label it "average-structure
   model of an RMC configuration", never a refined structure.

## Point-group facts (what the tools report)

- **Centrosymmetric** (an inversion centre): -1, 2/m, mmm, 4/m, 4/mmm, -3, -3m, 6/m, 6/mmm, m-3,
  m-3m. Every other point group has none.
- **Sohncke** (proper rotations only; only these allow a chiral structure): 1, 2, 222, 4, 422, 3,
  32, 6, 622, 23, 432. Any group with a mirror, a glide, -1 or a -4/-3/-6 axis is achiral.
- **Polar** (a direction no operation reverses): 1, 2, m, mm2, 4, 4mm, 3, 3m, 6, 6mm.
- So F-43m (-43m) and P-42₁m (-42m) are both non-centrosymmetric, achiral and non-polar; P2₁3
  (23) and P2₁2₁2₁ (222) are Sohncke groups, non-centrosymmetric and non-polar.

## Acceptance checklist

- [ ] Each group is named "consistent with ⟨group⟩ at ⟨tolerance⟩ Å", with its number when the
      card gives one, else its point group and crystal class.
- [ ] What changes between rungs is stated as lost operations, with the group–subgroup index and
      its type (t or k).
- [ ] Inversion, chirality and polarity come from the point-group facts the tools report.
- [ ] Rungs below the noise floor are called noise; a rung between the floor and the app's margin
      is "at the limit of the sampling", never a lowering.
- [ ] The group reported is the lowest rung at or above the margin, with the card's group beside
      it when it differs; a lowering is sized with `measure_symmetry_lowering`.
- [ ] Against a published cell-doubled phase the box can't hold, only the part the box can
      represent is compared, and the phase is neither confirmed nor ruled out.
- [ ] The average structure is not used as evidence for or against a local distortion.

## Known failure modes

- **"P-42₁m is chiral and centrosymmetry is lost"** (spot check, qwen3.8-27b on LM Studio):
  F-43m has no inversion centre to lose, and P-42₁m (-42m) contains improper operations, so it
  is achiral. *Firm rule in the app:* `detect_symmetry` and the context's ladder state the point
  groups' facts, and an answer that calls a group chiral, centrosymmetric or polar against them,
  or says an inversion centre is lost from a group without one, gets a warning under it.
- **Lower rungs read as symmetry breaking** (review): the Demo's P-42₁m … P1 rungs lie below the
  0.031 Å noise floor; the average structure is consistent with F-43m.
- **The 0.2 Å card read as the group** (validation arm A, the symmetry-lowering dev case): every Ta
  shifted 0.07 Å along [001] lies inside the card's tolerance, so the card says Pm-3m — the wrong
  group, called centrosymmetric and non-polar. The ladder and the noise floor show P4mm, and
  `measure_symmetry_lowering` gives Ta 0.065 Å along [0 0 1].
- **A group–subgroup relation that doesn't hold.** *Firm rule in the app:* an answer that claims a
  group–subgroup relation that doesn't hold gets a warning under it.
- **A cell doubling the box can't support.** *Firm rule in the app:* an answer that claims a cell
  doubling the box's average structure can't support gets a warning under it.

## Report

The group at the tolerance asked, as "consistent with ⟨group⟩ at ⟨tolerance⟩ Å" (with the ladder
around it); what changes between the rungs asked about, as operations lost with the index and its
type; which changes clear the app's margin and which are only at the limit of the sampling; the
point-group facts that matter for the question; and, for local questions, the analyses that can
answer them.

## References

- [lacunar-spinel-symmetry](lacunar-spinel-symmetry.md) — the Demo's ladder and the published phases of GaTa₄Se₈ and GaNb₄Se₈ with their point groups.

Further reading for Claude Code (not loaded in the app): `docs/algorithms/run-dashboard.md`
(Detected SG, the tolerance ladder, the CIF).
