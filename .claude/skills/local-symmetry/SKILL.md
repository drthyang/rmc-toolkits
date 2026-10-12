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

**Status: Draft: awaiting the user's interview.** Written from the 2026-10-10 spot check (the
space-group question with the chirality error), the review that followed, the Detected SG card's
rules and the symmetry finder's notes. A step marked "(to confirm with the user)" is not yet the
user's own practice: say in the answer which choice you made.

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
3. **What breaks between two groups.** Compare their point groups and centring: which
   operations are lost (a three-fold axis means cubic is lost; a lattice centring halves or
   quarters the cell's translations), which survive. Take inversion, chirality and polarity from
   the point-group facts, never from memory: a change that lowers the symmetry cannot "lose" an
   inversion centre the higher group never had.
4. **Is the change real?** A rung below `noise_floor_A` is within the noise of the site means. A
   rung above it is a distortion of the average structure; say how large (the tolerance at which
   it appears) and compare with the published low-temperature structure when there is one.
5. **Local versus average.** A cubic average with locally distorted clusters is possible and
   common: for local distortion, use `analyze_clusters` and the displacement analyses, not the
   ladder.

## Point-group facts (what the tools report)

- **Centrosymmetric** (an inversion centre): -1, 2/m, mmm, 4/m, 4/mmm, -3, -3m, 6/m, 6/mmm, m-3,
  m-3m. Every other point group has none.
- **Sohncke** (proper rotations only; only these allow a chiral structure): 1, 2, 222, 4, 422, 3,
  32, 6, 622, 23, 432. Any group with a mirror, a glide, -1 or a -4/-3/-6 axis is achiral.
- **Polar** (a direction no operation reverses): 1, 2, m, mm2, 4, 4mm, 3, 3m, 6, 6mm.
- So F-43m (-43m) and P-42₁m (-42m) are both non-centrosymmetric, achiral and non-polar; P2₁3
  (23) and P2₁2₁2₁ (222) are Sohncke groups, non-centrosymmetric and non-polar.

## Judgment calls

- The tolerance to quote as "the" space group of an RMC configuration: the card's 0.2 Å default,
  the lowest rung above the noise floor, or both (to confirm with the user).
- Whether a rung just above the noise floor is reported as a distortion or as "at the limit of
  the sampling" (to confirm with the user).
- How to name a non-standard result ("≥" lower bounds, a crystal class without a number) in a
  report (to confirm with the user).

## Acceptance checklist

- [ ] Each group is named with its tolerance (and number when the card gives one).
- [ ] What changes between rungs is stated as lost operations (axes, centring, mirrors).
- [ ] Inversion, chirality and polarity come from the point-group facts the tools report.
- [ ] Rungs below the noise floor are called noise, not symmetry breaking.
- [ ] The average structure is not used as evidence for or against a local distortion.

## Known failure modes

- **"P-42₁m is chiral and centrosymmetry is lost"** (spot check, qwen3.8-27b on LM Studio):
  F-43m has no inversion centre to lose, and P-42₁m (-42m) contains improper operations, so it
  is achiral. *Firm rule in the app:* `detect_symmetry` and the context's ladder state the point
  groups' facts, and an answer that calls a group chiral, centrosymmetric or polar against them,
  or says an inversion centre is lost from a group without one, gets a warning under it.
- **Lower rungs read as symmetry breaking** (review): the Demo's P-42₁m … P1 rungs lie below the
  0.031 Å noise floor; the average structure is F-43m.

## Report

The group at the tolerance asked (with the ladder around it); what changes between the rungs
asked about, as operations; which changes are above the noise floor; the point-group facts that
matter for the question; and, for local questions, the analyses that can answer them.

## References

- [lacunar-spinel-symmetry](lacunar-spinel-symmetry.md) — the Demo's ladder and the published phases of GaTa₄Se₈ and GaNb₄Se₈ with their point groups.

Further reading for Claude Code (not loaded in the app): `docs/algorithms/run-dashboard.md`
(Detected SG, the tolerance ladder, the CIF).
