### Validation: arm A+ reads the app's rules (architecture candidate 8; not a protocol change)

- Arm A+ (the engine sweep, `validation/armA.js`) applied three rules of its own where the app has
  one. It now calls the app's functions, held in `SWEEP_RULES`, so the control arm and the Copilot
  read a run the same way:
  - **The split site.** A+ used arm A's rule (PC1/PC2 ≥ 1.3 and κ₁ ≤ −0.5). It now uses
    `isSplitSite` (PC1 resolved, κ₁ ≤ −0.5 and at least 3 sampling errors below zero) with
    `list_sites`' static displacement δ = √(rms1² − rms⊥²) ± its error (`splitDisplacement`). Its
    readings say `list_sites` with κ₁, its z and δ. On synthetic prolate Gaussian sites, arm A's
    rule called 3.6 % (216 copies) to 14 % (100 copies) split; `isSplitSite` called none.
  - **The cluster elements.** A+ tried every same-element first shell peaking below 3.6 Å. It now
    tries the evidence gate's elements: those whose own first shell holds 0.5–4.5 neighbours. The
    rule is one function, `clusterElementsOf` (`llm/context/shells.js`, which now holds
    `CLUSTER_CN_MAX`), and `evidence.js` imports it.
  - **The lower group's tolerance.** A+ read the space group at the sampling noise floor. It now
    reads it at 1.25 × the floor: `lowerGroupTolerance` (`llm/context/runContext.js`), which is
    `measure_symmetry_lowering`'s default. The margin is `NOISE_FLOOR_MARGIN`, one constant read by
    the tool, the Workbench's tolerance question (`clarify.js`, until now its own `NOISE_MARGIN`)
    and A+.
- Arm A's rules are the frozen protocol's and are unchanged. On the dev set every arm A and D
  record is identical, and so is every score and claim. Every A+ record moves in its readings; on
  the ten paired cases:
  - The group is now read at 1.25 × the floor, with the same group in every case.
  - Eight cluster readings of extended networks (CN 6–12) are gone: the O of four perovskites,
    FeSi's Fe and Si, rock-salt Mg and O.
  - In dev-jahn-teller, two F sites are no longer split (2.8 sampling errors). The F claim stays.
  - dev-offcentre's δ is now 0.202 Å, where it was 0.201 (planted 0.2).
  - The judgment subset's five A+ records move in the same ways, with unchanged scores (their
    readings stay in the records and pins).
  - The moved records are re-recorded in `arm-a-dev.json`.
- What the model sees is unchanged: the model-visible snapshots pass as they were. Live runs stay
  comparable, and `PROTOCOL_VERSION` stays 2026-10-19.9.
- The skill lint gained a `rule` check (`RULES`). It flags a skill sentence that states one of the
  three rules otherwise. Its allow-list now holds five sentences, awaiting the user's skill review:
  - displacement-anisotropy SKILL.md:27 and significance-reading.md:17: a split site is a negative
    Mardia kurtosis;
  - coordination-shells SKILL.md:51: any polyhedron distortion goes to `analyze_clusters`;
  - local-symmetry SKILL.md:38: the floor itself is the line;
  - local-symmetry SKILL.md:59: the margin is left to the user.
- Known-answer tests: `validation/__tests__/aplusRules.test.js`. Decision record:
  `validation/findings/decisions/2026-10-19-aplus-app-rules.json`.
