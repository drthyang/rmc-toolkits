### Shells: a joined shell's split by number (plan item 4.2, protocol 2026-10-19.15)

- The shell catalogue (`src/llm/context/shells.js`) gives a shell joined across a deep trough its
  `split`: the parts cut at the minimum of the smoothed g between each two successive resolved
  peaks (the first is the first peak alone), the last ending at the shell's boundary; each part's
  CN over exactly its window (the parts sum to the shell's), its share and its mean distance and
  spread from r²g; `longMinusShort`, the last part's mean minus the first's, with
  `longMinusShortSe`, the configuration's sampling error (each part's spread over √ its bonds in
  the box, the bonds taken as independent; not the fit's uncertainty); and `direction`, the count
  the join decided by. A 2 + 2 + 2 octahedron keeps three parts and long − short is the longest
  minus the shortest: in the E⊗e Jahn–Teller problem that measures the amplitude whatever the mix
  of the Q₂ and Q₃ modes, and the shares show the mix. A shell whose peaks overlap without a deep
  trough has no split (its minimum cuts through both bond groups).
- Known answers (`shells.test.js`): 4 × 1.93 + 2 × 2.25 Å at σ 0.05 Å on 0.005–0.02 Å grids reads
  4 + 2, shares 2/3 and 1/3, long − short 0.32 Å within 0.002 Å; 2 + 2 + 2 three thirds and 0.37 Å;
  the reverse direction; the cut's bias (+0.001–0.0017 Å at a 0.07 Å spread, +0.006 Å at the edge
  of the join); and the sampling error's 95 % coverage, 0.953 over 300 synthetic configurations.
- `find_gr_peaks` returns it as `first_shell.split` (`parts` with window, mean, sd, CN in the
  join's direction and share; `long_minus_short_A` ± `long_minus_short_se_A`), declares it in its
  dataShape, names it in its description, and its joined-shell note quotes the parts, long − short
  and the shares by path instead of the peaks' separation, so the model quotes the split by
  reference. dev-jahn-teller's Cu–F reads 3.98 F per Cu at 1.929 Å + 2.01 at 2.251 Å, long − short
  0.321 ± 0.0043 Å (planted 0.32). What the model sees: the three tool-list snapshots (+123
  characters each; the budget pins +48 and +49 tokens); the Demo joins no shell.
- Arm A+ takes a Jahn–Teller claim's magnitude from the catalogue's split where it joins the
  claimed shell (else the g(r) maxima's separation arm A reads). On the dev set one record moves:
  dev-jahn-teller's A+ claim reads 0.321 Å where it read 0.31, with the same score (re-recorded in
  `validation/baseline/arm-a-dev.json`). Approved by the manager on the user's delegation,
  2026-10-19 (`validation/findings/decisions/2026-10-19-jt-split.json`).
