### Skills: the facts synced with the tools (the user's decision of 2026-10-19)

- The batch-1 skill drafts now state the numbers and rules that the tools apply. Nothing else in
  them changed: no method, step order or judgment call.
  - The Demo's Ta–Se first peak is 2.21–2.65 Å, holding 3.25 Se per Ta. The drafts said
    2.33–2.71 Å and 3.9, from before shells.js moved the first peak to the inter-peak minimum.
  - A split site is the app's rule: a resolved PC1 whose excess kurtosis is ≤ −0.5 and at least 3
    sampling errors below zero. A negative Mardia kurtosis alone is not one.
  - A preferred direction needs the peak test at 3σ or more on the unrounded σ.
  - The low-index cone is "within 15°", 15° included. That is under half the 35.3° between <110>
    and <111>, not half of it.
  - A family is named only when the peak cell's radius is 15° or less (plan item 2.5).
  - compute_bond_angles' cut rule is the CN margin (max(0.5, 15 % of the shell's CN)), as #113
    made it. The counts are its fallback when there is no shell to compare.
  - analyze_clusters runs for elements whose own first shell holds 0.5–4.5 neighbours of their
    own kind. It reads the shape against a shuffle null at each cluster's own size, and the size
    against independent atoms.
  - The average structure's group is read from 1.25 × the noise floor.
  - The order rule is ratio ≥ 3 and ≥ 3σ, against independent cells (F1).
  - The cluster tool's z is no tail probability.
  - The Displacement Directions page opens at Auto, the engine's resolution, with its numbers from
    the raw histogram (the skill said ν = 10).
  - A wrong reason in rmc-fit-assessment step 7 (long waves "left out of the background") is
    removed.
- The skill lint's allow-list holds no FACT entry any more, #113's three "awaits the skill
  fact-sync PR" entries included. One JUDGMENT stays: how to word a rung just above the noise
  floor. The bindings follow the new text:
  - the CN margin is bound to `CN_CHANGE_ABS` and `CN_CHANGE_REL`, with #113's
    compute_bond_angles asserts;
  - the counts fallback to `MODAL_COORDINATION_SHARE` and `UNBONDED_APEX_FRACTION`;
  - the split rule to `SPLIT_SITE_KURTOSIS` and `SPLIT_SITE_SE`;
  - the cell radius to `directionResolved`.
- `docs/SKILL_INTERVIEWS.md` lists the eleven choices the sync found and left to the user. The
  read_skill and over-claim record snapshots are re-recorded. The 32k budget pins move by 853 tokens
  (dock) and 854 (arm B). At 48k, 6,886 and 9,615 tokens are left to spare.
- **rmc-fit-assessment interviewed (2026-10-24).** The skill is now the user's own answers. Its
  status reads "Interviewed 2026-10-24", and `docs/SKILL_INTERVIEWS.md` §1 records each answer and
  where it landed. The rules it now states:
  - convergence judged on every fit term and the total;
  - a "still improving" caveat only for first-shell quantities;
  - no fixed Rw cutoff;
  - a peak constraint-limited above 5 %;
  - at least 3 independent runs and 10 accepted moves per atom past convergence before a
    distribution;
  - the user's list of what is never a finding, and the report's always and nevers;
  - two new firm rules;
  - the user's wording for long waves.

  The code half (every-term convergence, the 5 % constant, the two answer checks, the replays) is
  its own PR. Until it lands, four skill-lint entries are marked "awaits the fit-assessment code
  PR". `copilotSkills.test.js` tells an interviewed skill from a draft.
- **coordination-shells interviewed (2026-10-24).** The skill is now the user's own answers, and
  §2 records where each landed. The rules it now states:
  - the window cross-checked against the average structure's bond lengths;
  - a split shell's long − short ± its sampling error always quoted;
  - a shell not separated from the next cut at the minimum, with the 0.2 Å range beside it;
  - the cut rule kept, with no exemption for disordered or glassy systems (the user's window
    instead);
  - angles read from the sin θ-corrected histogram;
  - a distortion needs the cluster effect at z ≥ 3, the same sign in ≥ 3 independent runs, and
    angle spreads;
  - extended networks get angle spreads plus a polyhedral distortion index;
  - a new firm rule: every CN or angle carries its window. It waits for a later coordination code
    PR, as a skill-lint entry.
- **local-symmetry interviewed (2026-10-24).** The skill is now the user's own answers, and §3
  records where each landed. The rules it now states:
  - the 1.25 × margin, with a rung below it "at the limit of the sampling", never a lowering;
  - "the" group is the lowest rung at or above the margin;
  - a lowering is described by its lost operations plus the group–subgroup index and its t or k
    type;
  - against a published cell-doubled phase, only what the box can represent is compared;
  - CIF export rules;
  - groups named as "consistent with ⟨group⟩ at ⟨tolerance⟩ Å";
  - two new firm rules, which wait for a later local-symmetry code PR as skill-lint entries.

  The lint's last JUDGMENT entry is gone: the margin is bound to `NOISE_FLOOR_MARGIN`.
- **displacement-anisotropy interviewed (2026-10-24).** The skill is now the user's own answers,
  and §4 records where each landed. The rules it now states:
  - "the most disordered site": split sites first, then U_iso and anisotropy;
  - a negative-kurtosis site that misses the split rule is "flat-topped, not split";
  - an element pooled only for an element-wide summary, with the frames caveat;
  - three tiers for a direction: a finding at ≥ 3σ and ≥ 1.2× isotropic, a hint at ≥ 3σ under
    1.2× or at 2.5–3σ, nothing below 2.5σ;
  - a finding holds at the Auto resolution and at a neighbouring ν;
  - every reported direction carries ± its cell's radius, and the 15° cone stays;
  - a cluster distortion needs z ≥ 3 with the same sign in at least 3 independent runs;
  - a new firm rule: a split or flat-topped site is called disorder only after
    `analyze_modulations` ran for its element.

  The tiers, the two resolutions and the firm rule wait for a later displacement-anisotropy code
  PR, as five skill-lint entries. The skill is now the largest the Copilot reads (8,896 of
  read_skill's 9,000 characters). The 32k budget pins move by 868 tokens, and at 48k 5,685 (dock)
  and 8,415 (arm B) tokens are left to spare.
- **local-order interviewed (2026-10-24): all five skills are now the user's own.** §5 records
  where each answer landed. The rules it now states:
  - the list adds `analyze_modulations`' diffuse verdict and, for an extended network,
    `polyhedral_distortion`;
  - the diffuse verdict, past its calibrated bar, counts as order between cells;
  - directions for every symmetry-inequivalent site, not only the outliers;
  - a partial under 5 % of the signal is "not constrained by the data", and "no order" leaves it
    out;
  - never a bare "no order": a test that could not run is named with its reason;
  - a new firm rule: "no order" is refused until every applicable analysis has run.

  Three skill-lint entries wait for their code: the extended evidence gate and the diffuse
  verdict's reading (a later local-order answer-check code PR), and `polyhedral_distortion` (its
  own PR). No snapshot or budget pin moves.
