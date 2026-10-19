### Copilot: `list_sites` shows the split (plan 4.1)

- Each `list_sites` row now carries, where it applies, the excess kurtosis along each resolved
  principal axis (`kurtosis_axes`, null where an axis is not resolved), PC1 as a direct-lattice
  direction (`pc1_uvw`, when PC1 is resolved), and for a split site its static displacement with its
  sampling error (`delta_A`, `delta_se_A`), so an answer quotes the engine's number by reference
  instead of working out √(rms1² − rms⊥²) by hand. A site is split by the rule the cluster summary
  already used (PC1 resolved, PC1 excess kurtosis ≤ −0.5 and three sampling errors below zero), now
  one function, `isSplitSite`, for both tools.
- δ² = λ₁ − (λ₂ + λ₃)/2: rms⊥ is the RMS of the two other axes, unbiased where PC2 and PC3 are
  equal (rms2 alone, the larger of two equal sample variances, reads δ low). Its error comes from
  the eigenvalues' sampling errors, (κ + 2)·m₂²/n each, by the delta method
  ([ai-assistant.md](algorithms/ai-assistant.md), Step 20). Known answer over 1,000 seeds for each of
  five planted splits (equal and unequal occupation, 216–512 copies, anisotropic thermal motion;
  cubic, hexagonal and monoclinic cells): (δ̂ − δ)/SE has a mean of −0.01 to 0.07, an sd of
  0.99–1.04 and 94–96 % within ±1.96, and PC1 lies within 7.6° of the planted [u v w]. No Gaussian
  ellipsoid is called split. On the dev cases the planted 0.2 and 0.1 Å splits read
  0.202 ± 0.0028 and 0.0989 ± 0.0029 Å.
- Split sites come first in the table, so trimming to the result budget never drops one before
  another site, and a note check names them; the tool's description lists the new fields. What the
  model sees changes: the `list_sites` result and description snapshots. Arms A and A+ do not read
  `list_sites`, and their dev records are unchanged. Approved by the manager on the user's
  delegation, 2026-10-19.

### Copilot validation: one file per decision record

- Decisions on the validation itself are kept one JSON file per record,
  `validation/findings/decisions/<id>.json` (the id: the date and a slug), read in id order by
  `parseFindingFiles` and validated by the findings schema, so parallel changes never edit the same
  file. `decisions.jsonl` is gone: its one record (the recomposed held-out breathing case) moved
  there unchanged, beside plan 4.1's, `2026-10-19-list-sites-split.json`. `run.js findings <folder>`
  validates such a folder.
