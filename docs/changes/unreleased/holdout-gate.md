### Copilot validation: the held-out test set sealed before the freeze

- The test cases' definitions move out of `validation/cases.js` into one sealed module,
  `validation/holdout/testCases.js`, with the test set's arm A / A+ / D records
  (`holdout/arm-a-test.json`, formerly `baseline/`). Builders and agents never open, print, grep
  or import that folder. `cases.js` keeps only the public ids (`TEST_CASE_IDS`), so `splitOf`, the
  findings schema and the records still name a test case. A generator helper only a test case used
  moved with it.
- One opt-in loader, `validation/holdoutLoader.js` (`loadTestCases`, `casesFor`): a dynamic import
  that runs only on `run.js`'s `--final` path (`arm-a`, `live` and `write` with `--split test`; the
  command line needs `--final` whatever the environment) or when `RMC_VALIDATION_HOLDOUT=1` is set;
  anything else gets `HoldoutSealedError` before importing. CI's frontend test step sets the flag.
  Every consumer goes through it (`run.js`, `split.test.js`, `scoring.test.js`, `caseBlind.test.js`);
  `answerCorpus.test.js` and the findings need only the public ids. A test-set file written inside
  the repository goes under the sealed folder or a gitignored output folder.
- Locally the tests that need the definitions are skipped (`it.skipIf`), and the holdout scan looks
  for the ids only; CI runs everything. Nothing a run sends changes and no case changes:
  `buildCase` gives byte-identical output for every case, dev and test, before and after, and
  `TEST_SET_SHA256` holds (checked with the flag). No protocol bump.
- New `validation/__tests__/holdoutGate.test.js`: no module of `src/**` or `validation/**` but the
  loader names the sealed folder in code (a small lexer reads every string and template literal:
  static and dynamic imports, `require`, globs, paths); no test file or model-visible snapshot
  writes a held-out id and no test title is built from the held-out list; the loader and the CLI
  refuse without `--final` or the flag.
- The holdout scan widens: it now reads the validation itself outside the sealed folder (its js
  and md files and its findings, the decision records among them) and
  `docs/algorithms/copilot-validation.md` (ids allowed there, nothing else), besides every md under
  `docs/` (the change notes and the CHANGELOG). New token kinds: the seed (next to the word
  "seed"), the supercell (next to its element), a thermal width a case sets itself and a breathing
  sd; a lattice parameter and a thermal width count next to a cation, not the shared anion. Its two
  lists move to `validation/__tests__/holdoutLists.js`, with four id-only entries (the protocol
  document, the id manifest, the one decision record on a test case, the lists themselves).
- Every held-out value is scrubbed from the docs and records: copilot-validation.md (the Cases
  table, protocol notes 2026-10-19.2 and .11, the Arm A tables' test readings),
  the recomposed-case decision record, the 2026-10-19.11 change note, `judgment-spec.md`, the
  validation's comments, AGENTS.md and, once and on the manager's instruction, two passages of
  docs/CHANGELOG.md. Each now points at the sealed module and `TEST_SET_SHA256`; the scrubbed
  facts are kept behind the seal. Git history keeps the old text: anyone who read it is closed out
  of science, threshold and scoring work on the validation.
