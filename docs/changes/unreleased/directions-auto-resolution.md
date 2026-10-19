### Displacement Directions: the page and the tool at the same settings, and the direction's resolution (plan item 2.5; protocol 2026-10-19.17)

- **The page opens at Auto**, the engine's recommended resolution (the finest ν whose cells still
  average at least 12 displacements), which is the tiling `analyze_displacement_directions` reads.
  Until now it opened at ν = 10. The Auto choice is labelled with the ν the engine chose for the site
  on screen (`Auto · ν=2 (42)`), also while another resolution is picked. A picked resolution
  overrides Auto and stays, like every option on the page: across a site change, a tab switch and
  another run, until the user picks again or a Copilot step's Show in sets it.
- **Every readout comes from the raw histogram.** Both engines pick the peak cell on the smoothed
  map, so with smoothing on, the peak's cell, direction and σ moved with the Smoothing slider (on a
  dev case at ν = 10, two passes turned a −0.28σ peak into −10.4σ). The page now asks for the
  readouts at smoothing 0 (the tool's request) and draws the same request smoothed: the smoothing
  only redraws the map, its colorbar and the hover's enhancement. A smoothing change asks for the
  map only. The map test, the ± asymmetry and the Bingham anisotropy never depended on it.
- **The map opens unsmoothed** (smoothing 0, was 2). Auto's cells already average at least 12
  displacements, and on its coarse tilings two passes left about 30 % of a cell's excess in place.
  This is UI only: no readout, and no score of the validation's arms, depends on it.
- **The direction's resolution is reported.** The peak direction is the peak cell's centre, so it
  names a direction only to within the cell. Both engines (`orientation.py`, `workers/orientation.js`)
  now return `peakCellRadiusDeg`, the peak cell's circumscribed radius: the largest angle from its
  centre to its corners, so every direction counted in the cell lies within it — 37.4° at ν = 1,
  19.3–20.9° at ν = 2, under 15° from ν = 3. The page prints the peak "± R°, the cell radius at
  ν = n"; `analyze_displacement_directions` gives `peak.radius_deg` and the same wording, and above
  15° (`LOW_INDEX_CONE_DEG`) says the direction is not resolved finer than R° and names no <100>,
  <110> or <111> family. The orientation golden pins the field in both runtimes.
- One module holds the page's defaults and both requests (`src/directionSettings.js`), and
  `src/useDirectionHistogram.js` asks for both; `OrientationView` draws `result` and prints its
  readouts from `readout` (`DirectionReadouts`).
- **Validation:** arm A, and A+ through the same reading, reads a split site's direction map at
  the page's own defaults (`ARM_A_RULES.DIRECTION_SETTINGS` is `DIRECTION_DEFAULTS`); the reading
  records the ν it read and the peak cell's radius. Arm A takes a significant peak whatever its
  radius, as a user reading the page; arm A+ applies the tool's bound and takes PC1 where the peak
  cell is wider than 15° (`SWEEP_RULES.directionResolved`). `PROTOCOL_VERSION` 2026-10-19.17. On the
  dev set nine A and A+ records move, every one through that reading. One judgment case's arm A
  score falls: a significant peak on Auto's 12-cell tiling names a direction 31.7° off its cube axis
  (beyond the 15° tolerance) where PC1 named it before. Its A+ score is unchanged, since A+ keeps
  PC1 there. The moved records are re-recorded in `validation/baseline/arm-a-dev.json`, the others
  as recorded. The test set was not run.
- Known answers: on every site of the dev thermal case the page's default readout equals the
  tool's result field for field, and arm A reads exactly that; its flattened O sites clear 3σ on
  the map test there and stayed below it on the old default. A planted [0 0 1] split site of 216
  copies reads ν = 1, ± 37.4° and the unresolved wording; at 4000 copies, ν = 5 and resolved.
- What the model sees changes in `analyze_displacement_directions` only: its description names the
  cell radius (+21 characters a request) and its result gives `peak.radius_deg`, a longer note and
  the new wording. The model-visible snapshots, the record snapshots and the context-budget pins are
  updated for it. The displacement-anisotropy skill still says the page defaults to ν = 10 (two
  sentences), for the skill fact-sync PR.
- Decision record: `validation/findings/decisions/2026-10-19-directions-auto-resolution.json`.
  Approved by the manager on the user's delegation, 2026-10-19.
