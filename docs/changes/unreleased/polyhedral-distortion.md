### Copilot: `polyhedral_distortion`, a polyhedron's distortion indices (protocol 2026-10-24.3)

- A new engine, `rmc_toolkits/distortion.py` (`polyhedral_distortion`, exported from
  `rmc_toolkits`; source of truth) with its straight port `workers/distortion.js`: a centre
  element's coordination polyhedra — its ligands within an inclusive window, from the bond-angle
  engine's own neighbour search, the polyhedra of the most common count measured and the rest
  counted — by Baur's bond-length distortion index DI = (1/n) Σ |lᵢ − l̄| / l̄ (*Acta Cryst.* B30,
  1974) and, for octahedra and tetrahedra, Robinson's quadratic elongation ⟨λ⟩ = (1/n) Σ (lᵢ/l₀)²
  (l₀ of the regular polyhedron of the same volume) and bond-angle variance σ²θ over the cis
  angles (12 of an octahedron, the trans pairs being the matching with the largest angle sum;
  6 of a tetrahedron; θ₀ = 90° or arccos(−1/3)) (*Science* 172, 1971). Each is summarised (mean,
  sd, 10th, 50th and 90th percentile) over the configuration's polyhedra and over the average
  structure's (the site means of the average-structure basis), so a static distortion separates
  from thermal spread.
- The PCA worker's `distortion` request and `GET /api/distortion` (`rMax` ≤ 8 Å; a window holding
  more than 10⁶ bonds is refused before any is stored); the request contract gains the route, its
  limits and 17 cases. Known answers: regular polyhedra 0, 1 and 0; the 4 × 1.93 + 2 × 2.25 Å
  Jahn–Teller octahedron DI = 0.069831 exactly; trigonal and flattened polyhedra and an off-centred
  atom in closed form; a thermal-only network whose mean DI is √(2/π)·√(5/6)·s/l̄ for a bond spread
  s, and whose average octahedron is regular. A parity golden (`distortion_fixture.json`) holds the
  port to 1e-9, refusals included.
- The AI Copilot's `polyhedral_distortion(centre, ligand, r_min, r_max, window_confirmed?)` renders
  the configuration's and the average structure's indices. It refuses a window that cuts the
  coordination shell by `compute_bond_angles`' rule (the Workbench asks first), unless the user chose
  it. Its description says only what it computes, which is what the validation's arm C reads; the
  method advice (where the window comes from, octahedra or tetrahedra isolated or in an extended
  network, the static/local reading) goes to the dock and arm B. On the Demo the TaSe6 octahedra
  carry a static 3 + 3 distortion (average structure DI 0.0415, configuration 0.0469), while the GaSe4
  tetrahedra are regular on average (0.0010 against 0.0278).
- What the model sees: the tool's definition (+1,745 characters in the dock, +1,741 in arm B, +1,246
  in arm C), the first request's tool lists, and a new tool-results sample. The 32k budget pins move
  by +698 (dock) and +693 (arm B) tokens; at the live runs' 48k the dock keeps 7,041 tokens and arm B
  9,776, with the batch-1 skill drafts. Arms A and A+ do not read the tool, and their dev records are
  unchanged. Approved by the
  manager on the user's delegation, 2026-10-24 (the user's interview answer 7 for coordination-shells;
  `validation/findings/decisions/2026-10-24-polyhedral-distortion.json`).
