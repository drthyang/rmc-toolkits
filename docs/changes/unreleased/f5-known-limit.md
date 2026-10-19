### Shells: a split M₃/M₄ cluster shell does not join (F5, a known limit)

- A split first shell joins only when its summed CN is an ideal polyhedral count (4, 6, 8, 12), so a
  metal triangle or tetrahedron whose edges a distortion splits stays cut at its first trough. Accepting
  CN 2 and 3 was tried and left out: on the Demo and the dev cases it helps no shell and breaks one
  (dev-cell-wrap's V dimers read as an extended network; arm A+ there 1.00 → 0.83). The limit is
  documented with its evidence in `docs/algorithms/ai-assistant.md` and pinned by a known-answer test.
  Decided by the manager on the user's delegation, 2026-10-19
  (`validation/findings/decisions/2026-10-19-f5-shell-join.json`); the user may reopen it.
