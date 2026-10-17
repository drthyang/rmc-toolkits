# Skill interviews

A skill is the user's method in their own words (`.claude/skills/<name>/SKILL.md`), read by the AI
Copilot (`read_skill`) and by Claude Code alike. The five skills here are **drafts**: the first four written on
2026-10-11 from the Workbench's own record (the llm README's spot check and referee-style review,
`docs/algorithms/*`, the analyses' code and the papers the README cites), not yet from the user.
This page is the guide for turning each into the user's method. Each answer lands in one of three
places:

- **the skill's text:** the order of steps, the judgment calls (each marked "(to confirm with the
  user)" in the draft), the acceptance checklist, what to report;
- **code:** a firm rule the app enforces. Either a tool refuses or flags a result (as
  `compute_bond_angles` now refuses a window that cuts the coordination shell, in
  `src/copilotTools.js`), or a check on the answer itself flags a contradiction
  (`src/llm/copilot/answerAudit.js`: "converged" against an improving run, a group called chiral
  or centrosymmetric against its point group, a preferred direction no analysis found);
- **an eval scenario:** a scripted-model replay on the Demo run
  (`src/__tests__/copilotSkillReplays.test.js`): the right behaviour passes and the old wrong one
  is caught, so a later change cannot quietly break it.

How to run one: ask a few questions at a time, starting from the draft ("the draft says X — is that
how you do it?"). Ask for a recent run the method was used on, and for one time it went wrong.
Write the answers into the skill in the user's terms, mark what is a firm rule, add a replay for
each failure they describe, remove the draft status line, and show them the diff.

## 1. RMC fit assessment (`rmc-fit-assessment`, draft)

1. The draft judges convergence only from the χ² history's recent slope (the watchdog), and the
   history is the last `.log` column. Is that how you decide a run has converged, or do you look
   at every fit term, the total, or the configurations themselves?
2. The draft treats "improving, but slowly" as not converged. When is such a run good enough to
   analyse, and when must it continue?
3. The draft runs `compare_fit` on every dataset and reads where the misfit sits in r or Q. Which
   misfit level, or which concentration in one range, makes you distrust a run, for x-ray and for
   neutron data?
4. The draft names peaks the closest-approach constraints shape (g(r) still high at the cut). How
   do you choose those constraints, and when do you rerun with them loosened?
5. The draft calls ~3.5 accepted moves per atom "light sampling". How many moves per atom, or how
   many independent configurations, do you want before you quote a distribution?
6. Overfitting: the draft lists an Rw below the noise, claims about a weakly weighted partial, and
   box-scale waves. What else do you look for, and which of these would you refuse to report?
7. What must a fit report always contain for you, and what must it never claim?
8. Firm rules: anything the Copilot must never say about a fit without asking you? (Today it
   flags "converged" against an improving history.)
9. A recent run you assessed this way; one time a fit you trusted turned out wrong.

## 2. Coordination shells and bond windows (`coordination-shells`, draft)

1. The draft takes the window from `find_gr_peaks`' first shell: onset to the first deep minimum,
   both peaks of a split shell. Is that how you choose windows, or do you set them from the
   average structure, from bond-valence lengths, or by eye on the partial g(r)?
2. The draft treats the TaSe₆ shell (peaks at 2.52 and 2.74 Å) as one octahedron. When do you
   report the short and the long bonds separately instead?
3. When the shell is not separated from the next (g(r) still above zero at the minimum), the draft
   quotes the CN range for a cut 0.2 Å earlier or later. Is that what you do, or do you cut at the
   minimum and say so?
4. The app now refuses a window in which fewer than 60 % of the central atoms share the most
   common bond count. Is 60 % the right line, and are there systems (disordered sites, glasses)
   where a spread coordination is real and the rule should not apply?
5. The draft reads angles with the sin θ effect in mind (raw counts deplete 180°) and compares the
   cis/trans split with an ideal octahedron (80 % at 90°, 20 % at 180°). Which form do you read
   or publish: raw counts, the per-degree density, or sin θ-corrected?
6. For a distortion the draft goes to `analyze_clusters` (against independent atoms) rather than a
   g(r) peak shift. Which statistics convince you that a polyhedron or a cluster is distorted?
7. Firm rules: besides the cut-shell refusal, anything the Copilot must never do with bond windows?
8. A recent run you analysed this way; one time a window misled you.

## 3. Local symmetry (`local-symmetry`, draft)

1. The draft reads the Detected SG card's tolerance ladder and calls every rung below the sampling
   noise floor (rms displacement / √copies) noise. Is that your criterion, or do you require a
   margin above the floor?
2. Which tolerance do you quote as "the" space group of an RMC configuration: the card's 0.2 Å
   default, the lowest rung above the noise floor, or both?
3. The draft names what breaks between two groups as lost operations (three-fold axes, a lattice
   centring). Is that the description you want, or the group–subgroup relation (index, irreps)?
4. The app now states each point group's facts (inversion centre, Sohncke, polar) and flags an
   answer that contradicts them. Are there other symmetry statements you want checked the same
   way (polarity, a group–subgroup claim, the cell doubling)?
5. How do you weigh the average structure against the published low-temperature phases
   (GaTa₄Se₈ P-42₁m below 50 K; GaNb₄Se₈ P2₁3 below 50 K, P2₁2₁2₁ below 29 K) when the RMC box
   cannot hold the doubled cell?
6. When would you export a lower-symmetry CIF from a rung, and what do you check in it first?
7. A recent run you read the ladder on; one time the ladder suggested a distortion that was not
   there (or missed one that was).

## 4. Displacement anisotropy (`displacement-anisotropy`, draft)

1. The draft ranks sites by Mardia kurtosis, then U_iso, and reads anisotropy as rms1/rms3. Which
   numbers do you look at first for "the most disordered site", and which do you report?
2. The draft analyses directions one site at a time and treats a pooled element map as unable to
   rule out a distortion. Do you ever pool, and how do you report it?
3. The draft calls a direction "preferred" only at ≥ 3σ (calibrated, corrected over all cells) and
   asks for the enhancement as well. Is 3σ your line? Which enhancement is too small to matter?
4. The sphere's resolution: the tool lets the engine choose, the page uses ν = 10. Must a result
   hold at two resolutions before you trust it?
5. For a cluster-level question (a trigonal versus a tetragonal distortion), the draft goes to
   `analyze_clusters`' axis families against the shuffle null. Is that the test you would use?
6. Per-axis kurtosis is shown only on axes the PCA resolves. Do you use it, and how do you read a
   negative kurtosis (a split site)?
7. Firm rules: anything about displacements the Copilot must never claim? (Today it flags a
   preferred direction that no analysis of the question found.)
8. A recent site you analysed this way; one time a direction or an anisotropy turned out to be
   noise.

## 5. Local order (`local-order`, draft, added in validation batch 1)

Written on 2026-10-17 from the 2026-10-14 tool review's proposal (approved by the user that day)
and the validation's dev cases: how to conclude "no order".

1. The draft says "no order" only when every analysis that could show one came back empty:
   directions per site, split sites against `analyze_modulations`, clusters, order between cells,
   a lowering of the average structure, split coordination shells. Is that your list? What would
   you add or drop?
2. Must every symmetry-inequivalent site have its directions tested, or only those whose
   kurtosis or anisotropy stands out?
3. A partial the data barely weight (a few % of the signal): must its polyhedra be tested before
   "no order", or is "not constrained by the data" the answer?
4. How do you word "no order" when a test is undefined (too few points per sphere cell, a
   zone-boundary point an odd supercell axis cannot hold)?
5. Firm rules: should the Copilot refuse a "no order" answer until each analysis ran? (Part B of
   the validation cycle, the evidence gate, does this for the analyses the draft lists.)
6. A recent run you called ordered or disordered; one time an order was there that you first
   missed.

## Batch 1 changes to the four drafts (2026-10-17)

Each was approved by the user on 2026-10-14 as a proposal of the tool review
([copilot-review.md](algorithms/copilot-review.md)); each is still a draft item to confirm.

- `coordination-shells`: a split shell whose two bond lengths are well resolved (a deep trough
  between them) is one polyhedron when a second peak within 0.5 Å brings the CN to 4, 6, 8 or 12
  (a Jahn–Teller octahedron is never four-coordinate). Is 0.5 Å the right reach, and are those the
  right counts?
- `displacement-anisotropy`: a split site is called disorder only after `analyze_modulations`; a
  low-index family is named only within 15°; a preferred direction needs the peak test (the map
  test is anisotropy). Is 15° your cone?
- `local-symmetry`: report the lowest rung above the noise floor, with the card's 0.2 Å group
  beside it, and size a lowering with `measure_symmetry_lowering`. Is that the group you quote?
- `rmc-fit-assessment`: a box-scale wave (every |q_i| ≤ 0.25) voids a claim of order at its own q,
  never a zone-boundary point or general mode that passes its own test. Is that your rule?
