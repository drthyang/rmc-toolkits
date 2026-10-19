### Copilot validation: the judgment dev subset (plan 1.6, part ii)

- Five dev cases where arm A+'s fixed joining rules are ambiguous or wrong, so a model has room to
  add. Each is marked `judgment: true` with a `situation`, and `devOnly` until its held-out match
  exists:
  - two features at once: Jahn–Teller octahedra whose centre is also off-centred;
  - a signal at 2.5–4σ: one weak cell-doubling order, one seed on each side of the 3σ bar;
  - two tools read differently: eight-site <111> off-centring;
  - a choice only the user can make: elongated octahedral clusters whose split shell the Workbench
    joins and asks about. Since protocol 2026-10-19.8 A+ reads that joined shell too, so on the
    score this case is the subset's control.
- Their parameters sit in `validation/judgmentCases.js`, apart from `validation/cases.js`, which
  builds them:
  - new options on existing kinds: `jahn_teller` `offCentre`, `zone_boundary` `weak` (no decoy),
    `offcentre` `wells: 8`;
  - one new kind, `elongated_clusters`.

  Without an option, a kind draws exactly as before: every earlier case, dev and test, builds byte
  for byte as it did.
- **A spec for the held-out author:** `validation/judgment-spec.md`. It gives the situations, the
  kinds and options, what makes each a judgment case and the knobs, with no dev value, result or
  arm score. `validation/__tests__/judgment.test.js` scans it, the docs section and this note for
  the dev cases' ids and values.
- **The report:** each split with judgment cases gets a table of each column's mean score per
  situation, over the subset and over the rest (`report.js` `renderJudgment`). The truth of a
  judgment case says so (`truth.judgment`), and a transcript's planted truth names the situation.
- **Two situations through the harness, scripted** (`judgment.test.js`): a window that cuts the
  octahedral clusters' split shell is asked about and the skip keeps the whole shell; a
  superstructure claimed below the bar scores as a claim while the answer check flags it.
- **The split test:** a judgment case may use a paired kind's options, and no test case is one yet.
  `TEST_SET_SHA256` is unchanged. The holdout scan never lets a judgment case remove a test token.
- Arms A, A+ and D on the five cases are appended to `validation/baseline/arm-a-dev.json`. The
  earlier records are unchanged.
- `PROTOCOL_VERSION` 2026-10-19.9 (2026-10-19.8 is the shell catalogue's): a dev batch now holds the
  subset. Approved by the manager on the user's delegation
  (`validation/findings/decisions/2026-10-19-judgment-dev-cases.json`).
