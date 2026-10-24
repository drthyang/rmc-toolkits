### Modulations: the independent-cell null (F1) and diffuse order between cells (plan item 4.4; protocol 2026-10-19.16)

- **F1.** The modulation engines (`modulations.py`, `workers/modulations.js`) test an element's
  power against its own cells placed at random: every correlation inside a cell is kept (the
  ligands of a Jahn–Teller octahedron, the vertices of a cluster), every correlation between cells
  destroyed. The null's cross-spectral matrix is white, N/(N − 1) Σ_c y(c) y(c)ᵀ at every q ≠ 0, so
  the power is a weighted sum of χ²₁ (zone boundary) or χ²₂/2 (general q) terms with its
  eigenvalues as weights (`independent_cell_null` / `independentCellNull`); the tail is a
  Lugannani–Rice saddlepoint, 0 to +8 % above the exact tail from 10⁻² to 10⁻⁹ (`null_tail` /
  `nullTail`). Each element reports `nullDof`. The null of independent sites it replaces found
  order in 8.8–31.8 % of fresh configurations of the dev Jahn–Teller, two-feature and octahedral generator
  kinds (1,000 each); the calibrated null in 0.1–0.6 %, each of its three tests at its 0.135 %. It
  also stops holding an independent site to one component: a weak order on one isotropic site is
  found in 99 % of configurations where it was found in 44 %, and both dev weak orders are now far above the
  3σ bar.
- **4.4.** Each element gains `diffuse` and `verdict.diffuseOrder`: the ranks of the power over the
  compared modes (no long wave, no zone-boundary point, none found ordered) regressed on the 13
  neighbour-cell lags cos(2π q·R), χ²_r under independent cells, Šidák over the elements; diffuse
  order needs 3σ and a neighbour-cell correlation |Γ| ≥ 0.1. Its Finding carries one Evidence
  record in the cluster engine's form, `clears` recomputable. Parseval fixes the total power, so
  the test reads how the power is spread. Planted short-range order (neighbour correlation ∓0.27)
  is found in every configuration, in antiphase or in phase; no configuration without order between cells
  reads as it beyond its bar (`tests/test_modulation_calibration.py`: a fast CI subset, and 1,000
  configurations per family with `RMC_TOOLKITS_FULL_CALIBRATION=1`). One strong mode never makes it. The
  Demo's three elements are correlated in phase beyond the long waves.
- `analyze_modulations` shows a `diffuse` row [lag, correlation, σ] per element (F2's rule on the
  correlation's 0.1 bar), says antiphase correlation as short-range order and in-phase correlation
  in one note as neighbouring cells moving together, and reads its σ "against independent cells".
  A nearest miss at σ = −∞ reads "< −8σ" (it printed "nullσ"). Its result budget is 4,500
  characters. What the model sees changes in this tool's result only: both tool-result snapshots
  (`analyze_modulations.json`, `analyze_modulations.element.json`).
- Arm A+ reads cellDoubling, generalOrder and longWave only (no feature kind for diffuse order).
  On the dev set one A+ record moves: a weak order now above the bar is claimed (re-recorded in
  `validation/baseline/arm-a-dev.json`); every arm A and D record and every other score is
  unchanged. The marginal judgment situation needs its amplitude re-set from the calibrated null,
  the user's call. Approved by the manager on the user's delegation, 2026-10-19
  (`validation/findings/decisions/2026-10-19-modulation-null.json`, `2026-10-19-diffuse-order.json`).
