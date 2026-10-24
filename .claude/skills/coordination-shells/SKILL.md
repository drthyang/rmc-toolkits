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

**Status: Draft: awaiting the user's interview.** Written from the 2026-10-10 spot check (the
Se–Ta–Se question that went wrong), the review that followed, the Bond Geometry algorithm notes
and the tools' checks; the split-shell rule of step 2 was added in validation batch 1 (the
2026-10-14 tool review's L1, approved by the user that day). A step marked "(to confirm with the
user)" is not yet the user's own practice: say in the answer which choice you made.

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
   bonds, never a four-coordinate Cu. `first_shell.window_A` is always the whole polyhedron.
3. **Check the count.** Compare the shell's CN with the ideal. Short means the window cuts the
   shell. A shallow trough (the "not separated" warning) means the CN depends on the cut: quote
   the range for a cut 0.2 Å earlier or later, and ask whether the crystal has A–B distances in
   the gap.
4. **Choose the window.** Use `first_shell.window_A` unless the user asked for something else.
   When the data offer two defensible windows (a split shell, a shell not separated from the
   next) and the question does not say which, ask with `ask_user`, offering each window with its
   CN and your default (the whole shell). The Workbench asks by itself before an analysis runs in
   a window that cuts the shell (the first peak alone of a joined split shell included).
5. **Compute** `compute_bond_angles` (or `analyze_clusters`) in that window. The app refuses a
   window that cuts the shell: the angles of a cut shell are not the polyhedron's. Pass
   `window_confirmed: true` only when the user chose or typed that window.
6. **Read the angles.** The angle modes and the 10° histogram are raw counts, which geometry
   depletes near 0° and 180° (as sin θ): an octahedron's trans angle shows up as a broad shoulder
   below 180°, not as "no 180° angles". Read `histogram_10deg_sin_corrected` (1 = random
   directions) for the shape. Compare cis ≈ 90° and trans ≈ 180° for an octahedron; 109.5° for
   a tetrahedron; 60° inside an M₄ tetrahedron.
7. **Distortion.** A g(r) peak's position or width cannot tell a distortion from thermal spread.
   For a cluster or polyhedron distortion (Jahn–Teller, dimers, trimers), run `analyze_clusters`
   (shape against a shuffle null at each cluster's own size, size against independent atoms) where
   the element's own first shell is a finite cluster's (0.5–4.5 neighbours of its own kind; an
   extended network has none), or compare angle spreads; quote the effect size with its z.

## Judgment calls

- The window when the shell is not separated from the next (cut at the minimum, or report a CN
  range) (to confirm with the user).
- Whether a split first shell is reported as one polyhedron or as short and long bonds (to
  confirm with the user).
- The split-shell rule's reach (a second peak within 0.5 Å of the first) and its ideal counts
  (4, 6, 8, 12, within half an atom) (batch 1; to confirm with the user).
- The CN margin the app calls a window cut by, fewer partner atoms than the whole shell by at
  least max(0.5, 15 % of its CN) (the user's choice of rule, 2026-10-19), and with no shell to
  compare the 60 % modal-coordination share (to confirm with the user).

## Acceptance checklist

- [ ] The window comes from `find_gr_peaks`' first shell (or the user's choice), stated in Å.
- [ ] The CN in that window is quoted and compared with the ideal polyhedron.
- [ ] A split shell whose summed CN is a polyhedron is one polyhedron with two bond lengths, not a
      lower coordination (the first peak alone only when the user chose it).
- [ ] No angle, CN or cluster number from a window that cuts the shell.
- [ ] Angles read with the sin θ effect in mind; trans angles not called missing.
- [ ] Bond lengths and spreads quoted from the tool (`bond_ab.mean_A`, `sd_A`), never invented.
- [ ] A distortion claim rests on `analyze_clusters` or angle spreads, not on a g(r) peak shift.

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

## Report

The window and why; the CN against the ideal; the angle modes or spreads with what they mean for
the polyhedron; bond lengths with their measured spread; for a distortion, the effect size and
its z; and any constraint that shapes the shell.

## References

- [lacunar-spinel-shells](lacunar-spinel-shells.md) — the shells, windows, coordination numbers and published distortions of GaTa₄Se₈ and GaNb₄Se₈.

Further reading for Claude Code (not loaded in the app): `docs/algorithms/bond-geometry.md`.
