### Modulations: the polarisation of an ordered mode (plan item 4.3), F2 and F3 (protocol 2026-10-19.12)

- `displacement_modulations` (`modulations.py`, `workers/modulations.js`) gives every zone-boundary
  point and strongest mode its `polarisation`: the leading eigenvector of the real part of the
  element's cross-spectral tensor at q, (1/m) Σ_s Re[U_s U_sᴴ] / (n_s tr C_s), whose trace is the
  power. It comes as a Cartesian unit `direction` and a direct-lattice `uvw` (largest component 1,
  the form of `list_sites`' `pc1_uvw`), both in the sense that makes the closest crystal axis acute,
  with its `share` λ₁/tr (1: one direction; ⅓: none preferred). The share is not a significance:
  one site at a zone-boundary point gives share 1 by construction. `element_spectrum` /
  `elementSpectrum` return the power and the tensor field. The parity fixture pins every
  polarisation, and its generator refuses an eigenvalue gap or a sense decided near its tolerance.
- Known answers: both sites of an element displaced 0.1 Å along the direct-lattice [1 0 1] of a
  monoclinic cell (14° from the Cartesian (1, 0, 1)) at q = (½, 0, 0) are found within 5°, with a
  share above 0.95; over 60 seeds the scatter is the sampling error σ/(δ√(mn)) rad per component.
  The existing planted modes are read along their axes.
- F3: the Šidák family is every element of the configuration with displacements, whichever element
  the call names. Asking about one element after seeing the data is selective inference, and
  naming it shrank the family to its own points: a weak doubling at 2.98σ over the 21 points of
  three elements read 3.30σ, "ordered", over its own 7. The payload's `family` names its
  `elements`. A call without an element is unchanged.
- `analyze_modulations` shows each ordered point's polarisation (a `polarisation` row: [u v w] to
  two decimals, share, the nearest <100>/<110>/<111> within 15°) and says it in the point's check.
  F2: a ratio or σ below its bar (3, 3σ) is never printed rounded onto it. Where the usual rounding
  reads at or above the bar, it is shown with the fewest significant digits that stay below (2.98,
  not 3), and nearer than six digits as `'<3'`; the no-order sentence adds "below the 3σ bar". The
  note says the correction is over every element's points. The tool gets its own 4,000-character
  result budget (`MODULATION_RESULT_BUDGET`). What the model sees changes in this tool's result
  only: the tool-result snapshot and a new sample of a call naming one element.
- Arm A+ reads a modulation claim's direction from the engine's polarisation, no longer from the
  PC1 of the element's most prolate site. On the dev set every arm A and D record, and every
  score, is unchanged; two A+ records move in their claim's direction and reading (re-recorded in
  `validation/baseline/arm-a-dev.json`). Approved by the manager on the user's delegation,
  2026-10-19 (`validation/findings/decisions/2026-10-19-mode-polarisation.json`).
