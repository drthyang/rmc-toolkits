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
