# Skill interviews

A skill is the user's method in their own words (`.claude/skills/<name>/SKILL.md`), read by the AI
Copilot (`read_skill`) and by Claude Code alike. The five skills here began as **drafts**: the first four written on
2026-10-11 from the Workbench's own record (the llm README's spot check and referee-style review,
`docs/algorithms/*`, the analyses' code and the papers the README cites), not yet from the user.
All five were interviewed on 2026-10-24 (§1–§5).
This page is the guide for turning each into the user's method. Each answer lands in one of three
places:

- **the skill's text:** the order of steps, the judgment calls (each marked "(to confirm with the
  user)" in the draft), the acceptance checklist, what to report;
- **code:** a firm rule the app enforces. Either a tool refuses or flags a result (as
  `compute_bond_angles` now refuses a window that cuts the coordination shell, in
  `src/copilotTools.js`), or a check on the answer itself flags a contradiction
  (`src/llm/copilot/verification/facts.js`: "converged" against an improving run, a group called chiral
  or centrosymmetric against its point group, a preferred direction no analysis found);
- **an eval scenario:** a scripted-model replay on the Demo run
  (`src/__tests__/copilotSkillReplays.test.js`): the right behaviour passes and the old wrong one
  is caught, so a later change cannot quietly break it.

How to run one: ask a few questions at a time, starting from the draft ("the draft says X — is that
how you do it?"). Ask for a recent run the method was used on, and for one time it went wrong.
Write the answers into the skill in the user's terms, mark what is a firm rule, add a replay for
each failure they describe, remove the draft status line, and show them the diff.

## 1. RMC fit assessment (`rmc-fit-assessment`, interviewed 2026-10-24)

Answered by the user on 2026-10-24. Each answer, and where it landed: the skill's text, code (the
fit-assessment code PR, its own PR off main) or an eval scenario (its replays in
`src/__tests__/copilotSkillReplays.test.js`, in the fit-assessment code PR). Until that PR lands,
the skill lint holds each line today's code does not apply as an entry marked "awaits the
fit-assessment code PR".

1. *Convergence: the last `.log` column, or every fit term, the total, the configurations?*
   Judge every fit term and the total; say "converged" only when none is still falling beyond its
   noise over the last tenth of the run. **Landed:** the skill (step 2, the checklist) and code
   (every-term convergence, in the fit-assessment code PR; lint entry: the convergence rule).
2. *When is "improving, but slowly" good enough?* First-shell quantities (distances, CN) may be
   answered with a "still improving" caveat; for subtle effects (distortions, directions, order)
   the run must continue. **Landed:** the skill (step 2) and an eval scenario (in the
   fit-assessment code PR).
3. *Which misfit makes you distrust a run, for X-ray and for neutron data?* No fixed Rw cutoff.
   Distrust the run when most of the misfit sits in one r or Q range that the box size or the
   scale doesn't explain, or when residual features rise above the data's noise; the same rule for
   X-ray and neutron. **Landed:** the skill (step 3, the checklist).
4. *Constraints: when is a peak constraint-shaped, and when do you rerun?* g(r) above 5 % of the
   peak at the cut means the peak is constraint-limited: interpret nothing on its low-r side, and
   suggest rerunning with the constraint loosened. **Landed:** the skill (step 4, the checklist)
   and code (the 5 % constant, in the fit-assessment code PR; today's code flags it at 2 %: lint
   entry).
5. *How much sampling before a distribution?* At least 3 independent runs or configurations, and
   at least 10 accepted moves per atom past convergence, before quoting a distribution or a small
   effect; below that the answer says so (3.5 moves per atom is too light). **Landed:** the skill
   (step 5, the checklist, the report) and an eval scenario (in the fit-assessment code PR).
6. *Overfitting: what do you refuse to report?* As a finding, never: an Rw below the noise
   presented as a better fit; a detailed claim about a partial under 5 % of the signal; order at a
   box-scale q; any feature that differs between independent runs. **Landed:** the skill (step 6)
   and an eval scenario (in the fit-assessment code PR).
7. *What must a report always contain, and never claim?* Always the datasets and their weights,
   the convergence status, each dataset's Rw and where its misfit sits, the constraint-limited
   pairs, and the sampling with the number of runs; never "refinement", never esds as if
   parameters were refined, never order at box-scale q. **Landed:** the skill (the report).
8. *Firm rules?* Besides the existing converged-versus-improving flag: (a) a claim of order at a
   box-scale wavevector is flagged; (b) a distribution quoted from a single configuration without
   saying so is flagged. **Landed:** the skill (two firm-rule sentences, its known failure modes)
   and code (the two answer checks, in the fit-assessment code PR; lint entries until then, then
   mapped in FIRM_RULES).
9. *A recent run; a fit you trusted that turned out wrong.* Not answered on 2026-10-24: open.

The long waves, in the user's words (answering the batch-1 item below and the fact sync's item 5):
"a long wave is a correlation the finite box can't represent, so it says nothing about order at
its own q; a point or mode that passes its own test at another q is judged on that test alone."
**Landed:** the skill (step 7), verbatim. The draft's last sentence there (a long wave that
dominates every element says the box is small or lightly sampled), which the answers do not
cover, is gone with it.

## 2. Coordination shells and bond windows (`coordination-shells`, interviewed 2026-10-24)

Answered by the user on 2026-10-24. Each answer, and where it landed: the skill's text or code (a
later coordination code PR). Until that PR lands, the skill lint holds the line today's code does
not apply as an entry marked "awaits the coordination code PR".

1. *How do you choose windows?* The window defaults to `find_gr_peaks`' first shell,
   cross-checked against the average structure's bond lengths; flag when the two disagree by more
   than the shell's width. **Landed:** the skill (step 4, the checklist). No Workbench tool reports
   an A–B bond length of the average structure yet, and the skill says so.
2. *When do you report the short and long bonds separately?* A split shell whose summed CN is
   ideal is one polyhedron with two bond lengths; the report always quotes `first_shell.split`'s
   long − short ± its sampling error. **Landed:** the skill (step 2, the checklist, the report).
   This keeps `find_gr_peaks`' join as the app applies it (reach 0.5 Å, ideal counts 4, 6, 8,
   12); the reach itself was not discussed.
3. *A shell not separated from the next: a CN range, or a cut at the minimum?* Cut at the minimum
   by default, and quote the CN range for a cut 0.2 Å earlier and later beside it. **Landed:** the
   skill (step 3, the checklist, the report).
4. *Is the CN margin the right line; any exemption?* The cut rule stays as it is: max(0.5, 15 % of
   the shell's CN), with the 60 % modal-count fallback where there is no shell. No exemption in
   code for disordered or glassy systems: those go through the user's `window_confirmed`, and the
   report says the user chose the window. **Landed:** the skill (step 5, the firm rule unchanged,
   the report).
5. *Which angle form do you read or publish?* Read and report the shape from
   `histogram_10deg_sin_corrected`; give raw counts only when asked. **Landed:** the skill (step
   6, the checklist, the report).
6. *Which statistics convince you of a distortion?* `analyze_clusters`' effect size at z ≥ 3,
   the same sign in ≥ 3 independent runs, and supporting angle spreads; a g(r) peak shift alone
   never proves one. **Landed:** the skill (step 7, the checklist).
7. *Firm rules besides the cut-shell refusal?* (the user's answer 8) Every CN or angle in an
   answer carries its window in Å; an answer that gives one without it is flagged. **Landed:** the
   skill (a firm-rule sentence) and code (its answer check, a later coordination code PR; lint
   entry until then).
8. *A recent run; a window that misled you.* Not answered on 2026-10-24: open.

Also answered (the user's answer 7): corner-sharing or extended networks (outside the 0.5–4.5
own-kind rule) use angle spreads plus a standard polyhedral distortion index such as Baur's, never
the cluster null. No Workbench tool computes the index yet, and the skill says so. **Landed:** the
skill (step 7).

## 3. Local symmetry (`local-symmetry`, interviewed 2026-10-24)

Answered by the user on 2026-10-24. Each answer, and where it landed: the skill's text or code (a
later local-symmetry code PR). Until that PR lands, the skill lint holds the lines today's code
does not apply as entries marked "awaits the local-symmetry code PR".

1. *A margin above the noise floor?* The margin is 1.25 × the noise floor. A rung between the
   floor and 1.25 × it is reported as "at the limit of the sampling: neither a distortion nor
   ruled out", never as a lowering. **Landed:** the skill (step 4, the checklist, the report); the
   lint binds the margin to the app's `NOISE_FLOOR_MARGIN` (1.25) wherever the skill states it,
   and its last JUDGMENT entry (local-symmetry:70) is gone.
2. *Which group do you quote?* The lowest rung at or above 1.25 × the floor, with the card's 0.2 Å
   group beside it when they differ, sized by `measure_symmetry_lowering`. **Landed:** the skill
   (step 5, the checklist).
3. *Lost operations, or the group–subgroup relation?* The operations lost, plus the group–subgroup
   index and its type, t (translationengleiche) or k (klassengleiche); irreps only on request.
   **Landed:** the skill (step 3, the checklist, the report). No Workbench tool states the index or
   the type yet, and the skill says so.
4. *Other symmetry statements to check?* (a) A claimed group–subgroup relation that doesn't hold is
   flagged; (b) a claim about a cell doubling that the box's average structure can't support is
   flagged. **Landed:** the skill (two firm-rule sentences) and code (their answer checks, a later
   local-symmetry code PR; lint entries until then).
5. *Against published low-temperature phases the box can't hold?* Compare only the part the box can
   represent, the t-subgroup or point group, and state that the box can neither confirm nor rule
   out the cell-doubled phase. **Landed:** the skill (step 4, the checklist).
6. *When do you export a lower-symmetry CIF, and what do you check first?* Only for a rung at or
   above 1.25 × the floor. First check that the model reproduces the site means within the rung's
   tolerance, and that the Wyckoff sites and occupancies make sense. Label it "average-structure
   model of an RMC configuration", never a refined structure. **Landed:** the skill (step 7). The
   app's CIF header words it differently: "Symmetry-averaged RMC average structure, written by
   rmc-toolkits." and `_audit_creation_method` 'rmc-toolkits symmetry-averaged RMC structure'
   (`src/cifWriter.js`). It never says "refined"; aligning it with the user's label is a code
   question.
7. *A recent run; a ladder that misled you.* Not answered on 2026-10-24: open.

Also answered: naming. "Consistent with ⟨group⟩ at ⟨tolerance⟩ Å", never "the structure is
⟨group⟩"; with no space-group number, give the point group and the crystal class. **Landed:** the
skill (step 5, the checklist, the report, and the reference's wording of the Demo's group).

## 4. Displacement anisotropy (`displacement-anisotropy`, interviewed 2026-10-24)

Answered by the user on 2026-10-24. Each answer, and where it landed: the skill's text or code (a
later displacement-anisotropy code PR). Until that PR lands, the skill lint holds the lines today's
code does not apply as entries marked "awaits the displacement-anisotropy code PR".

1. *Which numbers first for "the most disordered site", and which do you report?* Split sites
   first (the app's rule, by `delta_A`), then U_iso and anisotropy (rms1/rms3); report each with
   U_iso, anisotropy, kurtosis and, when split, `delta_A`. **Landed:** the skill (step 1, the
   checklist, the report). `list_sites` and the run context list the other sites by |kurtosis|:
   the skill ranks them by their U_iso and anisotropy (`largest_U_iso_refs`,
   `most_anisotropic_refs`), and keeps |kurtosis| for "the most anharmonic site", the draft's
   reading, which the answer did not address.
2. *Do you ever pool, and how do you report it?* Only for an element-wide summary, always with the
   frames caveat; an isotropic pooled map is never evidence against a distortion. **Landed:** the
   skill (step 2, the checklist).
3. *Is 3σ your line; which enhancement is too small?* Three tiers. A finding: ≥ 3σ (unrounded)
   and ≥ 1.2× isotropic. A hint: ≥ 3σ under 1.2×, or 2.5–3σ. Below 2.5σ: nothing to report.
   **Landed:** the skill (step 4, the checklist, the report) and code (the tier labels with their
   1.2× and 2.5σ constants, a later displacement-anisotropy code PR: today
   `analyze_displacement_directions` calls a direction preferred at ≥ 3σ whatever its
   enhancement, and labels no hint; lint entries until then). The 1.2× judgment call is a rule now.
4. *Must a result hold at two resolutions?* A direction is a finding only if it holds at the
   engine's Auto resolution and at the neighbouring ν; hints don't need it. **Landed:** the skill
   (step 5) and code (a later code PR: today one analysis at one resolution backs a preferred
   direction; a lint entry until then). The answer names no side, so the skill reads ν − 1 or
   ν + 1. The skill also reads a direction that holds only at the Auto ν as a hint, a case the
   answer leaves open.
5. *Is `analyze_clusters`' axis-family test the one you would use?* Yes: its axis families
   against the shuffle null, at z ≥ 3 with the same sign in at least 3 independent runs (the bar
   of §2, answer 6). **Landed:** the skill (step 7, the checklist); the lint binds z ≥ 3 to
   `CLUSTER_Z` and lists the runs as no gate (every tool reads one run).
6. *Per-axis kurtosis; a negative kurtosis?* A site is split only by the app's rule; a
   negative-kurtosis site that misses it is "flat-topped, not split". Quote a per-axis kurtosis
   only on resolved axes. **Landed:** the skill (step 1, the checklist) and its reference
   (`significance-reading.md`).
7. *Firm rules?* (the user's answer 8) An answer that calls a split or flat-topped site
   "disorder" when `analyze_modulations` has not run for its element is flagged. **Landed:** the
   skill (a firm-rule sentence) and code (its answer check, a later code PR; a lint entry until
   then). Today the evidence gate sends back a split site called random disorder until
   `analyze_modulations` ran in the question, for any element, and reads no flat-topped site.
8. *A recent site; a direction or an anisotropy that turned out to be noise.* Not answered on
   2026-10-24: open.

Also answered (the user's answer 7): keep the 15° cone, 15° included, and name a family only when
the peak cell's radius is also ≤ 15°; every reported direction carries "± R°", its cell's radius.
**Landed:** the skill (step 6, the checklist, the report); the tool already applies both
(`LOW_INDEX_CONE_DEG`, `directionResolved`) and gives the radius (`peak.radius_deg`).

## 5. Local order (`local-order`, interviewed 2026-10-24)

Written on 2026-10-17 from the 2026-10-14 tool review's proposal (approved by the user that day)
and the validation's dev cases: how to conclude "no order". Answered by the user on 2026-10-24.
Each answer, and where it landed: the skill's text or code (a later local-order answer-check code
PR; the extended networks' tool is its own PR). Until they land, the skill lint holds the lines
today's code does not apply as entries marked "awaits the local-order answer-check code PR" or
"awaits the polyhedral_distortion PR".

1. *Is the draft's list yours; what would you add or drop?* Keep it, and add (a)
   `analyze_modulations`' diffuse short-range-order verdict (#112) and (b), for an extended
   network, `polyhedral_distortion`. **Landed:** the skill (steps 3 and 4, the checklist, the
   firm rule) and its reference (`order-signatures.md`: a row for diffuse order).
   `polyhedral_distortion` is being built in its own PR: until it is on main the skill names it
   "where the Workbench offers it", and the lint holds it as a tool that awaits that PR.
2. *Every symmetry-inequivalent site, or only the outliers?* Every symmetry-inequivalent site.
   **Landed:** the skill (step 2, the checklist).
3. *A partial the data barely weight?* A partial under 5 % of the signal is reported as "not
   constrained by the data", and the "no order" conclusion explicitly excludes it. **Landed:** the
   skill (steps 5 and 7, the checklist, the report); the lint lists the 5 % as no gate (the context
   gives the weights, no module decides at 5 %), as in rmc-fit-assessment.
4. *How do you word "no order" when a test is undefined?* "No order found in the tests that could
   run; ⟨test⟩ not tested because ⟨reason⟩", never a bare "no order". **Landed:** the skill (step
   7, the checklist, the report).
5. *Should the Copilot refuse a "no order" answer until each analysis ran?* Yes: a "no order"
   answer is refused until every applicable analysis on the list has run, the two additions
   included. This extends the evidence gate. **Landed:** the skill (a firm-rule sentence beside the
   one that states today's gate) and code (the extended gate, a later answer-check code PR; a lint
   entry until then). Today the gate asks for `analyze_modulations`, one per-site directions map
   and `analyze_clusters` for each cluster element.
6. *A recent run; an order you first missed.* Not answered on 2026-10-24: open.

Also answered (the user's answer 6, closing the fact sync's item 11): the diffuse verdict, when it
passes its calibrated bar, counts as order between cells, and arm A+ and the answer checks read
it. **Landed:** the skill (step 3) and code (a later answer-check code PR; today neither reads it,
`validation/armA.js` leaving it as "the user's call"; a lint entry until then). The tool says an
in-phase correlation as neighbouring cells moving together, as RMC's long waves and acoustic-like
thermal motion give, "not a superstructure": whether an in-phase verdict counts as order is open
for that PR.

## Batch 1 changes to the four drafts (2026-10-17)

Each was approved by the user on 2026-10-14 as a proposal of the tool review
([copilot-review.md](algorithms/copilot-review.md)); each is still a draft item to confirm.

- `coordination-shells`: a split shell whose two bond lengths are well resolved (a deep trough
  between them) is one polyhedron when a second peak within 0.5 Å brings the CN to 4, 6, 8 or 12
  (a Jahn–Teller octahedron is never four-coordinate). Is 0.5 Å the right reach, and are those the
  right counts? **Answered 2026-10-24** (§2, answer 2: the join kept; the reach not discussed).
- `displacement-anisotropy`: a split site is called disorder only after `analyze_modulations`; a
  low-index family is named only within 15°; a preferred direction needs the peak test (the map
  test is anisotropy). Is 15° your cone? **Answered 2026-10-24** (§4: the cone kept, 15° included;
  the modulation test a firm rule; the peak test the first of the three tiers).
- `local-symmetry`: report the lowest rung above the noise floor, with the card's 0.2 Å group
  beside it, and size a lowering with `measure_symmetry_lowering`. Is that the group you quote?
  **Answered 2026-10-24** (§3, answer 2: at or above 1.25 × the floor).
- `rmc-fit-assessment`: a box-scale wave (every |q_i| ≤ 0.25) voids a claim of order at its own q,
  never a zone-boundary point or general mode that passes its own test. Is that your rule?
  **Answered 2026-10-24** (§1, the long waves).

## Fact sync (2026-10-19): what it left for the interviews

The user decided on 2026-10-19 that the skills' facts are fixed now and the interviews come after
v1. The fact sync changed only numbers and rule statements, so that they match the tools. The
skill lint (`src/__tests__/skillLint/tables.js`) binds each of them to its constant or probe. The
sync changed no method, no step order and no judgment call. It found these choices and leaves them
to you. Line numbers are in the synced text.

1. `local-symmetry/SKILL.md:70`, a judgment call: how to word a rung just above the noise floor.
   The app reads a group only from 1.25 × the floor (`lowerGroupTolerance`, and the Workbench's
   tolerance question), and step 4 now says so. Is 1.25 × your margin? How do you report a rung
   between the floor and 1.25 × it? This is the one entry left in the lint's allow-list, a
   JUDGMENT. **Answered 2026-10-24** (§3, answer 1: "at the limit of the sampling"; the entry is
   gone).
2. `displacement-anisotropy/SKILL.md:29–32`, `significance-reading.md:19–21` and
   `local-order/SKILL.md:26–29`: a split site is now stated as the app's rule. The rule is a
   resolved PC1 whose excess kurtosis is ≤ −0.5 and at least 3 sampling errors below zero;
   `list_sites` gives such a site a `delta_A`. Any negative Mardia kurtosis no longer counts. Is
   that your line? How do you report a flat-topped site that is not split? **Answered
   2026-10-24** (§4, answer 6: the rule kept; such a site is "flat-topped, not split").
3. `coordination-shells/SKILL.md:71–73` and `:91–95`: the cut rule is now the CN margin you chose.
   Are 0.5 atoms and 15 % the right values? Should the counts fallback (60 %, or most central
   atoms unbonded) stay where there is no shell to compare? **Answered 2026-10-24** (§2, answer
   4: as it is).
4. `coordination-shells/SKILL.md:57–61` (step 7) and `local-order/SKILL.md:38–40` (step 4):
   `analyze_clusters` runs only for an element whose own first shell holds 0.5–4.5 neighbours of
   its own kind. The evidence gate and arm A+ read the same rule. For a corner-sharing
   polyhedron, an extended network, the skill offers only "compare angle spreads". Is that your
   method there? **Answered 2026-10-24** (§2: angle spreads plus a polyhedral distortion index).
5. `rmc-fit-assessment/SKILL.md:52–55` (step 7): the draft's reason why a long wave cannot void a
   point that passes its own test is gone. The draft said "the test already leaves long waves out
   of the background it compares with", which does not hold:
   - the ratio's background is the median over every q, long waves included;
   - the percentile does leave long waves out, but it is no part of the order rule.
   Which reason do you give, if any? **Answered 2026-10-24** (§1, the long waves).
6. `displacement-anisotropy/SKILL.md:47–48` (step 3): the tool decides a preferred direction on
   the unrounded peak σ, and the skill now says that a σ just below 3 is never rounded to 3σ.
   Should a near miss (2.9σ, say) be reported as a hint, as the 1.2× judgment call does for a small
   enhancement? **Answered 2026-10-24** (§4, answer 3: 2.5–3σ is a hint).
7. `displacement-anisotropy/SKILL.md:72–73`: is a 3–5σ direction under 1.2× isotropic a finding
   or a hint? No module decides this; the lint lists the number as no gate. **Answered
   2026-10-24** (§4, answer 3: a hint; the 1.2× line now awaits the code PR in the lint).
8. `coordination-shells/SKILL.md:29–39` and its report (`:111`): `find_gr_peaks` now splits a
   joined shell by number, in `first_shell.split`: each part's CN and share, and long − short
   with its sampling error. The skill does not mention it. Should a Jahn–Teller report quote it?
   **Answered 2026-10-24** (§2, answer 2: always).
9. `displacement-anisotropy/SKILL.md:53–57` (step 4): the cone now reads "within 15°, 15°
   included", as the tools apply it, and since plan item 2.5 also needs the peak's cell to be no
   wider than 15° (`peak.radius_deg`; the tool names no family above it). The draft's reason said
   15° was half of 35.3°; it is under half (35.3° / 2 = 17.6°). The question stays the one from
   batch 1: is 15° your cone? **Answered 2026-10-24** (§4, answer 7: 15°, 15° included).
10. `displacement-anisotropy/SKILL.md:112–113` (the report): the tool now gives the peak direction
   ± its cell's radius (plan item 2.5). Should a report quote the radius beside the direction?
   **Answered 2026-10-24** (§4, answer 7: every reported direction carries it).
11. `local-order/SKILL.md:34–37` (step 3): `analyze_modulations` now also reports short-range order
   between cells spread over many q (its diffuse verdict, plan item 4.4). Neither arm A+ nor the
   answer checks read it ("the user's call", validation/armA.js). Does it count as order between
   cells before "no order"? **Answered 2026-10-24** (§5, answer 6: it counts, and arm A+ and the
   answer checks read it, in a later answer-check code PR).
