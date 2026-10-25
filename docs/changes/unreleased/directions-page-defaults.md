### Displacement Directions: the page opens at ν = 10, smoothed 4×

- The Displacement Directions page now opens at the resolution ν = 10 (1002 cells) with the map
  smoothed 4× (four neighbour passes), the user's choice. Auto, the engine's resolution and the
  tiling the AI Copilot reads, is still offered and still labelled with the ν it chose. Every
  readout under the sphere (the peak, its direction ± the cell radius, σ, the map, asymmetry and
  anisotropy tests) still comes from the raw histogram at smoothing 0: the smoothing redraws only
  the map. The page and the Copilot now agree only when Auto is picked.
- `analyze_displacement_directions` is unchanged (Auto, unsmoothed), and nothing the model is sent
  changes. The validation's arm A keeps the settings frozen at protocol 2026-10-24.6 (Auto,
  unsmoothed), now pinned in `validation/armA.js`; its dev records are identical (no protocol bump,
  no re-record). Decision record `2026-10-24-directions-page-defaults`.
