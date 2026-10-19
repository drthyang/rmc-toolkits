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
  - compute_bond_angles' cut rule is the CN margin (max(0.5, 15 % of the shell's CN)). The counts
    are its fallback when there is no shell to compare.
  - analyze_clusters runs for elements whose own first shell holds 0.5–4.5 neighbours of their
    own kind. It reads the shape against a shuffle null at each cluster's own size, and the size
    against independent atoms.
  - The average structure's group is read from 1.25 × the noise floor.
  - The order rule is ratio ≥ 3 and ≥ 3σ.
  - The cluster tool's z is no tail probability.
  - The page's resolution is "its default", not ν = 10.
  - A wrong reason in rmc-fit-assessment step 7 (long waves "left out of the background") is
    removed.
- The skill lint's allow-list holds no FACT entry any more. One JUDGMENT stays: how to word a rung
  just above the noise floor. The bindings follow the new text: the CN margin is bound to
  `CN_CHANGE_ABS` and `CN_CHANGE_REL`, the counts fallback to `MODAL_COORDINATION_SHARE` and
  `UNBONDED_APEX_FRACTION`, and the split rule to `SPLIT_SITE_KURTOSIS` and `SPLIT_SITE_SE`.
- `docs/SKILL_INTERVIEWS.md` lists the nine choices the sync found and left to the user. The read_skill
  and over-claim record snapshots are re-recorded. The 32k budget pins move by 660 tokens (dock)
  and 659 (arm B). At 48k, 7,088 and 9,818 tokens are left to spare.
