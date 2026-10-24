---
name: coordination-shells
description: >-
  Choose bond-length windows from the partial g(r) before any bond-angle, cluster or
  coordination analysis: find the whole first coordination shell (both peaks of a split shell),
  compare coordination numbers with the ideal polyhedron, never compute in a window that cuts a
  shell, and read angle distributions correctly. Use for bond angles, polyhedra, coordination
  numbers, and "is this polyhedron regular or distorted?". NOT for convergence
  (rmc-fit-assessment) or displacement directions (displacement-anisotropy).
---

# Coordination shells and bond windows

**Status: Interviewed 2026-10-24.** The steps are the user's own answers of that day
(docs/SKILL_INTERVIEWS.md §2).

A "bond" here is a distance window and nothing else: the engine reports geometry, not chemistry.
Every angle, coordination number and cluster shape depends on the window, so the window is the
analysis.

## Steps

1. **Name the polyhedron.** The central atom B and its neighbours A (Se–Ta–Se for the TaSe₆
   octahedron of GaTa₄Se₈, Se–Ga–Se for GaSe₄, Ta–Ta–Ta for the Ta₄ cluster) and the ideal
   coordination number (6, 4, 3).
2. **Find the shell:** `find_gr_peaks` for the A–B pair. Read `first_shell`, not the first peak:
   `window_A` (onset to the first deep minimum), `peaks_in_shell_A`, `mean_d_A` (from r²g, not
   the g peak), `cn` both ways, `cn_if_cut_0p2A_earlier_later`, `trough_g_fraction`. A shell
   with two peaks is ONE split shell, e.g. two bond lengths of one distorted polyhedron: the
   TaSe₆ shell has Ta–Se peaks at 2.52 and 2.74 Å; the first-peak window 2.21–2.65 Å holds
   3.25 Se per Ta, the whole shell 6. When the two bond lengths are well resolved, the trough
   between them is deep; `find_gr_peaks` then joins a later peak within 0.5 Å of the first whose
   summed CN is an ideal polyhedral count (4, 6, 8, 12), and reports the first peak alone as
   `first_shell.first_peak_only` with its own CN. A Jahn–Teller CuF₆ (4 × 1.93 + 2 × 2.25 Å) holds
   4 F in the first peak and 6 in `first_shell.window_A`: one octahedron with short and long
   bonds, never a four-coordinate Cu. `first_shell.window_A` is always the whole polyhedron. A
   split shell whose summed CN is ideal is one polyhedron with two bond lengths: quote
   `first_shell.split`'s long − short ± its sampling error (`long_minus_short_A` ±
   `long_minus_short_se_A`).
3. **Check the count.** Compare the shell's CN with the ideal. Short means the window cuts the
   shell. A shell not separated from the next (the "not separated" warning) has a CN that
   depends on the cut: cut at the minimum by default (`first_shell.window_A`), and quote beside
   it the CN range for a cut 0.2 Å earlier and later (`cn_if_cut_0p2A_earlier_later`).
4. **Choose the window.** Use `first_shell.window_A` unless the user asked for something else,
   and cross-check it against the average structure's bond lengths: flag when the two disagree
   by more than the shell's width (no Workbench tool reports an A–B bond length of the average
   structure yet). The Workbench asks by itself before an analysis runs in a window that cuts
   the shell (the first peak alone of a joined split shell included).
5. **Compute** `compute_bond_angles` (or `analyze_clusters`) in that window. The app refuses a
   window that cuts the shell: the angles of a cut shell are not the polyhedron's. Pass
   `window_confirmed: true` only when the user chose or typed that window. A disordered or
   glassy system has no exemption: a spread coordination goes through the user's own window
   (`window_confirmed`), and the report says the user chose it.
6. **Read the angles.** The angle modes and the 10° histogram are raw counts, which geometry
   depletes near 0° and 180° (as sin θ): an octahedron's trans angle shows up as a broad shoulder
   below 180°, not as "no 180° angles". Read and report the shape from
   `histogram_10deg_sin_corrected` (1 = random directions); give raw counts only when asked.
   Compare cis ≈ 90° and trans ≈ 180° for an octahedron; 109.5° for a tetrahedron; 60° inside an
   M₄ tetrahedron.
7. **Distortion.** A distortion needs `analyze_clusters`' effect size at z ≥ 3, the same sign in
   at least 3 independent runs, and supporting angle spreads; a g(r) peak shift alone never proves
   one (a peak's position or width cannot tell a distortion from thermal spread). For a cluster
   or polyhedron distortion (Jahn–Teller, dimers, trimers), run `analyze_clusters` (shape against
   a shuffle null at each cluster's own size, size against independent atoms) where the
   element's own first shell is a finite cluster's (0.5–4.5 neighbours of its own kind). For a
   corner-sharing polyhedron or another extended network, use angle spreads plus a standard
   polyhedral distortion index such as Baur's (`polyhedral_distortion`: Baur's DI and, for
   octahedra and tetrahedra, Robinson's quadratic elongation and bond-angle variance, over the
   configuration's polyhedra and the average structure's), never the cluster null.

## Acceptance checklist

- [ ] The window comes from `find_gr_peaks`' first shell (or the user's choice), stated in Å and
      cross-checked against the average structure's bond lengths.
- [ ] The CN in that window is quoted and compared with the ideal polyhedron; for a shell not
      separated from the next, beside the CN range for a cut 0.2 Å earlier and later.
- [ ] A split shell whose summed CN is a polyhedron is one polyhedron with two bond lengths, not a
      lower coordination (the first peak alone only when the user chose it), with its long −
      short ± its sampling error.
- [ ] No angle, CN or cluster number from a window that cuts the shell.
- [ ] Angles read and reported from `histogram_10deg_sin_corrected`, raw counts only when asked;
      trans angles not called missing.
- [ ] Bond lengths and spreads quoted from the tool (`bond_ab.mean_A`, `sd_A`), never invented.
- [ ] A distortion claim rests on the cluster effect size at z ≥ 3 with the same sign in at least 3
      independent runs and supporting angle spreads (an extended network: angle spreads and a
      distortion index), never on a g(r) peak shift.
- [ ] Every CN and angle carries its window in Å.

## Known failure modes

- **A window that cut the TaSe₆ shell** (spot check, 2026-10-10): asked for "Se–Ta–Se angles
  against an ideal octahedron", gemma4 used 2.3–2.7 Å, got CN 3.75 and called the octahedra
  "incomplete" (✗). qwen3.8 started at 2.33–2.71 Å (CN 3.86) and widened itself. *Firm rule in
  the app:* `compute_bond_angles` refuses, unless `window_confirmed` is set, a window that holds
  fewer partner atoms than the whole first shell by at least max(0.5, 15 % of the shell's CN);
  with no shell to compare (none with a CN, or a window apart from it), one in which fewer than
  60 % of the central atoms share the most common bond count, or most have none. It returns the
  numbers that show the cut, not the angles.
- **"There are no near-180° angles"** (spot check, after the fix): 20 % of the angles lay at
  150–170° in raw counts. The tool now returns the sin θ-corrected histogram too.
- **"± 0.02 Å" invented** for GaSe₄ (review): the tool now gives the measured spread (sd 0.095 Å).
- **"Peak shifted to 2.93 Å: distorted"** (review): the Ta–Ta g(r) peak sits below the shell
  mean (3.00 Å from r²g); the cluster analysis found the Ta₄ clusters 20 % *less* distorted than
  independent atoms at their own size.
- **A Jahn–Teller octahedron read as four-coordinate** (tool review, 2026-10-14, L1): the deep
  trough between 4 × 1.93 and 2 × 2.25 Å closed the shell at 2.10 Å (CN 4), and the long bonds
  looked like "the next shell". *Firm rule in the app:* `find_gr_peaks` joins the split shell
  (CN 6) and offers the first peak alone as `first_peak_only`, the Workbench asks before an
  analysis runs in it, and an answer that gives an element a coordination number its whole first
  shell contradicts gets a warning under it.
- **A CN or an angle without its window.** *Firm rule in the app:* every CN or angle in an answer
  carries its window in Å; an answer that gives one without it gets a warning under it.

## Report

The window in Å and why; the CN against the ideal (for a shell not separated from the next,
beside the CN range for a cut 0.2 Å earlier and later); for a split shell, `first_shell.split`'s
long − short ± its sampling error; the angle shape from `histogram_10deg_sin_corrected` with what
it means for the polyhedron; bond lengths with their measured spread; for a distortion, the
effect size and its z; a window the user chose, said to be the user's; and any constraint that
shapes the shell.

## References

- [lacunar-spinel-shells](lacunar-spinel-shells.md) — the shells, windows, coordination numbers and published distortions of GaTa₄Se₈ and GaNb₄Se₈.

Further reading for Claude Code (not loaded in the app): `docs/algorithms/bond-geometry.md`.
