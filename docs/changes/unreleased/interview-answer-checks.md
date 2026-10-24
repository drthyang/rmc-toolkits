### Copilot: the interview answers of four skills as checks — coordination shells, local symmetry, displacement anisotropy, local order (protocol 2026-10-24.2)

- **coordination-shells.** Every coordination number or bond angle in an answer carries its window
  in Å: an answer that gives one without it — in its sentence, or stated once for the pair or the
  whole answer — gets a `window` finding under it (the user's rule, 2026-10-24). The coordination
  check reads "4 F per Cu" as Cu's count of F, and a clause that states the whole shell's window as
  a claim about the shell (before, any window in Å made the clause partial). `find_gr_peaks` gives
  the average structure's first shell for the pair, `first_shell.average_structure_d_A` (its
  shortest and longest bond between the site means), and warns when a bond lies outside the window
  by more than the first peak's FWHM: no Demo or dev pair does (the Demo's Ta–Se 2.52–2.76 Å in
  2.21–3.23 Å); a 4 + 1 pyramid whose shell closes before its fifth bond does. The evidence gate
  backs distorted octahedra or polyhedra with `polyhedral_distortion` where the run offers it, else
  `analyze_clusters` (EVIDENCE_RULES says so), an `anyOf` row now asks for the first analysis the
  turn offers, and `analyze_clusters`' method keeps it to finite clusters. Approved by the user,
  2026-10-24 (`validation/findings/decisions/2026-10-24-coordination-shells-interview.json`).
- **local-symmetry.** A group–subgroup relation the answer claims that does not hold gets a
  `subgroup` finding: read against the run's own groups where it knows both (the ladder's rungs are
  nested; the index is |G|/|H| in the cell; t when the translations per cell stay, k when the point
  group stays, else a general subgroup, never maximal — the Demo's P-42₁m in F-43m is index 12,
  neither t nor k), and by the point groups alone for a published phase the run does not hold
  ("F-43m → P4/mmm" cannot be: 4/mmm has an inversion centre -43m lacks). A cell doubling the box
  cannot hold — odd along an axis it doubles — gets a `cell-doubling` finding, as does a doubling
  along axes no ordered zone-boundary point of `analyze_modulations` doubles. The symmetry-averaged
  CIF calls itself an "average-structure model of an RMC configuration", never a refined structure.
  Approved by the user, 2026-10-24
  (`validation/findings/decisions/2026-10-24-local-symmetry-interview.json`).
- **displacement-anisotropy.** `analyze_displacement_directions` labels its peak in the user's
  three tiers: a **finding** at ≥ 3σ (unrounded) and ≥ 1.2× isotropic that holds at the
  neighbouring resolution (ν − 1; ν + 1 at ν = 1), which the tool reads for a candidate — the peak
  test passes there and the two peaks lie on one axis within the two cells' radii; a **hint** at
  ≥ 3σ under 1.2×, at 2.5–3σ or not held; **nothing** below 2.5σ. The second reading costs 0.7 ms
  on the Demo and about 240 tokens per finding. On the Demo, Ta site 17 and Se site 34 are
  findings and Ta site 25 (4.47σ at ν = 2, 2.33σ at ν = 1) is now a hint; every dev split site is
  a finding. The direction check and arm A+ take a finding only (no dev verdict moves). An answer
  that calls a split or flat-topped site disorder before `analyze_modulations` tested its element
  gets a `disorder` finding (ordered antiphase displacements give the same split cloud), and the
  evidence gate reads a flat-topped site called disordered as it reads a split site. Approved by
  the user, 2026-10-24
  (`validation/findings/decisions/2026-10-24-displacement-anisotropy-interview.json`).
- **local-order.** "No order" now also needs `polyhedral_distortion` on each centre of octahedra or
  tetrahedra no finite cluster holds (the Demo's Ga, a perovskite's B site), where the run offers
  it; its diffuse verdict is `analyze_modulations`' own. The diffuse verdict is read as the user
  reads it: only antiphase (or general-q) order past its bar is order; in-phase correlation is
  "neighbouring cells moving together (correlated motion), not order", as `analyze_modulations` now
  says it, and never supports "no correlation". Three checks on the answer: `no-order` ("no order"
  when order was found), `correlated-motion` (order between cells when only in-phase correlation, or
  none, was found) and `no-correlation` (uncorrelated displacements when a correlation passed, in
  phase included). Arm A+ claims an antiphase verdict as a modulation with no q, the verdict's one
  kind for order between cells; no dev case passes the diffuse test, so no dev record moves.
  Approved by the user, 2026-10-24
  (`validation/findings/decisions/2026-10-24-local-order-interview.json`).
