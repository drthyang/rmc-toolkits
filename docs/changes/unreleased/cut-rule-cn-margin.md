### Copilot: one cut rule — compute_bond_angles refuses on the Workbench's CN margin (protocol 2026-10-19.18)

- Whether a bond window cuts its coordination shell is decided once, by the shell module's
  `windowCut` (`src/llm/context/shells.js`), for the Workbench's question before an analysis runs and
  for `compute_bond_angles`' refusal: a window that meets its first shell is cut when it holds fewer
  partner atoms than the whole shell by at least max(0.5, 15 % of its CN), whatever its coordination
  counts. Until now the tool refused on the counts (more than half the central atoms unbonded, or
  fewer than 60 % with the most common count) while the Workbench asked on the CN margin; the counts
  now decide only where the run reads no shell to compare the window with: no partials file, no
  closed shell for the pair, no density for a CN, or a window apart from the first shell (the module
  reads first shells only). Approved by the user, 2026-10-19
  (`validation/findings/decisions/2026-10-19-shell-cut-rule.json`).
- What changes: a whole shell is never refused (dev-offcentre's Nb–Nb, which runs into the next
  shell, 52 % of its Nb at the modal count), nor a window within the margin (the top 0.2 Å early);
  a sub-shell first peak is refused unless the user chose it (the Demo's Ta–Se, 3.25 of 6;
  dev-jahn-teller's four short Cu–F). On the 149 natural windows of the decision record, 15 change
  verdict (11 no longer refused, 4 refused now), and the Workbench asks on every window the tool
  refuses (`validation/run.js shell-cuts` shows both beside the old counts rule).
- What the model reads: the refusal names the CN shortfall by path ("it holds 3.75 Se per Ta
  (shell_ab.cn_in_window) of the shell's 6 (shell_ab.cn, in shell_ab.window_A)…"); every result adds
  the shell reading (`shell_ab`, and `shell_bc` for a B–C window of its own) and the counts' `modal`,
  `modal_share` and `unbonded_share` as information, with no note for a count spread; the
  coordination answer check counts a window as the whole shell when its CN is the shell's within
  the margin. The tool's description and schema are unchanged; five `compute_bond_angles` tool
  results and four Workbench tool messages grow by about 200–250 characters.
- `buildCopilotTools` takes the run's `plotFiles` and `runSettings`, so the tool reads the shell
  catalogue the Workbench reads; `window()` loses its `counts` option and the llm module exports
  `windowCut` instead of `countsCut` and its two thresholds, which stay in the shell module as the
  fallback.
- The coordination-shells skill still states the counts rule: the skill lint binds its 60 % and
  "most" to the CN margin and holds the three lines on its allow-list, with both sides, until the
  skill fact-sync PR. Arms A, A+ and D on dev are unchanged (A+ reads no bond-angle window).
