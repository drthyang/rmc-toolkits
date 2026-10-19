### Validation: the judgment subset's held-out matches join the test set (protocol 2026-10-19.11)

- Five judgment cases join the test set (sealed since the holdout gate:
  `validation/holdoutLoader.js`), one per dev judgment case, each with its situation, kind and
  option: two features at once (`jahn_teller` with `offCentre`), a
  signal at the edge of significance (`zone_boundary`, `weak: true`, `decoy: null`; two seeds of
  one amplitude, one on each side of the 3σ bar), two tools read differently (`offcentre`,
  `wells: 8`) and a choice only the user can make (`elongated_clusters`).
- Written blind from `validation/judgment-spec.md` alone: no dev judgment parameter, dev record,
  arm pin or result was read, and no arm or model ran on them. Each seed is the first of its own
  sequence that a rule fixed before any configuration was built accepts, from the case's engine
  readings; the property each rule checks is in docs/algorithms/copilot-validation.md (protocol
  2026-10-19.11), and the seeds and readings are sealed with the cases.
- `PROTOCOL_VERSION` is 2026-10-19.11 and `TEST_SET_SHA256` changes: results on the earlier test
  set are not comparable with this one.
