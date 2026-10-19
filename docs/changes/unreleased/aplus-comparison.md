### Copilot validation: the report beside the case score (plan 1.6, part i)

- `run.js report` takes each split's case score apart (`validation/comparison.js`). It adds three
  sections under the split's tables:
  - **the score's parts:** each column's feature-only score (the planted features alone), its
    spurious-claim test and its symmetry fields, with each live column minus A+ on the cases both
    have;
  - **every component against A+:** arm A, A+, each live column and its difference;
  - **every disagreement with the planted truth:** a missed feature and what took its place (wrong
    kind or wrong element), a direction, magnitude or q off, a spurious claim, a symmetry field
    wrong or unanswered, no verdict. Each gives the answer, the truth and the columns and runs that
    made it.
- The case score, every existing export of `validation/scoring.js` and the arm A/A+ pins are
  unchanged. New in `scoring.js`: `scoreComponents`, the components read back from the score,
  whose mean is the case score; and `restrictedScore`. `validation/__tests__/comparison.test.js`
  runs the dev baselines and synthetic B runs through the report: B matching A+, over-claiming,
  missing a feature, and better than A+ on the JT long axis. It also checks that the dev baselines
  re-score exactly as recorded. A split with no live run or no A+ record says so.
- On the dev baselines the feature-only score is A 0.25 and A+ 0.92, against case scores of 0.825
  and 0.99. A+'s one disagreement is the JT long axis.
- Not a protocol change: the report reads records and changes none, so `PROTOCOL_VERSION` stays
  2026-10-19.7.
