### Copilot: one shell reading per pair, per run (architecture candidate 6; protocol 2026-10-19.8)

- `src/llm/context/shells.js` reads the first coordination shell of every pair in the partials file
  once (`shellCatalogue(run)`, memoised per partials data, structure and run settings; frozen
  Shells), and every reader takes it from there: the run context's `pair_correlations`,
  `find_gr_peaks`, the Workbench's questions (`buildShellCatalog` is gone), the evidence gate, the
  answer checks (`cnMatches`), `compute_bond_angles`' window refusal (`countsCut`) and arm A+'s
  cluster windows. Before, five readers searched their own r ranges, rounded their own way and
  integrated the coordination number their own way.
- One rule per quantity: peaks below 8 Å and the closing minimum to max(8, r₁ + 3) Å, whatever range a
  caller passes; windows rounded to 0.01 Å and every CN integrated over exactly the window quoted
  with it (a count can move in its third digit: the Demo's Ta–Ta 2.76 → 2.75); separation and the
  closest-approach pile-up on one 2 % fraction (`find_gr_peaks`' FWHM edge clause is dropped; it
  moved nothing on the Demo or the dev cases); half an atom for every CN match.
- **The first peak alone ends at the minimum between the first two resolved peaks**, and exists only
  with two or more; r ± FWHM reached into the second peak. The Demo's Ta–Se first peak: 2.21–2.65 Å,
  3.25 Se per Ta (it read 2.33–2.71 Å, 3.86–3.9). Known answer: a 3 + 3 split with overlapping widths
  gives 3 ± 0.1 on three grids, where r ± FWHM gives 3.4–4.0. Approved by the manager on the user's
  delegation, 2026-10-19.
- `find_gr_peaks`: `first_peak_window_A` leaves; `first_shell.first_peak_only` (window, CN, trough)
  for every shell of two or more resolved peaks, its CN in the split-shell note; the shell no longer
  depends on `r_max` (a first peak above it gets a note and the shell, no longer "no peak"). The
  context adds `first_peak_to_A` / `first_peak_cn` for the same shells and lists peaks to 8 Å.
- The window-cut rule is unchanged: the Workbench asks on the CN margin, `compute_bond_angles` refuses
  on the coordination counts (60 % modal share, more than half unbonded; the Demo's 47 % text
  verbatim). Both readings live in the module; `validation/run.js shell-cuts` tabulates them on the
  dev windows (149: agree on 134, the counts alone refuse 11, the margin alone asks on 4). The
  CN-margin rule for the tool is a recommendation awaiting the user.
- The bond-angle engines (`triplets.py`, `workers/triplets.js`) return `coordinationSummary`
  `{mean, modal, modalShare, unbondedShare}`, the lowest count winning a tie; the Bond Geometry page
  and `compute_bond_angles` read it instead of recomputing it. Both golden fixtures pin it.
- The skill lint gains a value check (Demo values a skill quotes, against the tool): the four places
  the coordination-shells skill quotes the old first peak are allow-listed with both values until the
  user's skill review, a precondition for the freeze. What the model sees moves (`find_gr_peaks`, the
  context); the `compute_bond_angles` snapshots do not; the dock's worst turn grows by 74 tokens.
  Arms A and A+ on dev: scores and claims unchanged, three A+ cluster windows 0.01 Å lower at their
  start (re-recorded).
