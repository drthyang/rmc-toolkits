### Displacement Directions: the page and the tool at the same settings (plan item 2.5; protocol 2026-10-19.16)

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
- One module holds the page's defaults and both requests (`src/directionSettings.js`), and
  `src/useDirectionHistogram.js` asks for both; `OrientationView` draws `result` and prints its
  readouts from `readout` (`DirectionReadouts`).
- **Validation:** arm A, and A+ through the same reading, reads a split site's direction map at
  the page's own defaults (`ARM_A_RULES.DIRECTION_SETTINGS` is `DIRECTION_DEFAULTS`); the reading
  records the ν it read. `PROTOCOL_VERSION` 2026-10-19.16. On the dev set nine A and A+ records
  move, every one through that reading, and one judgment case's scores fall in both arms: a
  significant peak on Auto's 12-cell tiling names a direction 31.7° off its cube axis (beyond the
  15° tolerance) where PC1 named it before. The moved records are re-recorded in
  `validation/baseline/arm-a-dev.json`, the others as recorded. The test set was not run.
- Known answer: on every site of the dev thermal case the page's default readout equals the tool's
  result field for field, and arm A reads exactly that; its flattened O sites clear 3σ on the map
  test there and stayed below it on the old default.
- What the model sees is unchanged: the tool, its result and the model-visible snapshots are as
  they were. The displacement-anisotropy skill still says the page defaults to ν = 10 (two
  sentences), awaiting the user's skill review.
- Decision record: `validation/findings/decisions/2026-10-19-directions-auto-resolution.json`.
  Approved by the manager on the user's delegation, 2026-10-19.
