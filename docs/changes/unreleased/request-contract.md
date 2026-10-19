### Requests: one contract for the worker and Flask adapters (architecture candidate 9)

- `requestAnalysis` (useSiteCloud.js) is the analysis request seam, and its doc comment is the
  interface. `kind` must be a key of `ANALYSIS_ROUTES`, which moved into `workers/requestGuards.js`;
  the PCA worker's kinds are its keys. `params` are flat scalars that cannot override kind, text,
  dir or id. A refusal names the parameter. The comment also gives the cancel semantics.
  `pcaKdeWorker.js` exports `inProcessEngine`, the worker's handler in-thread behind the same
  interface. `openRunFolder` and the eight test helpers that rebuilt `{ kind, ...params }` use it.
- One hand-written fixture, `src/__tests__/fixtures/analysis_requests.json`, says what each kind
  accepts and refuses. It holds the 128-atom Nb/Se configuration, the six routes, 19 limits (each
  with the Copilot tool property that sends it) and 130 cases. `analysisRequests.test.js` replays it
  on the worker and in-process adapters and checks the Flask routing with axios mocked. It also
  runs the cancel table on the worker, Flask and in-process adapters and on useCopilotRun's engine,
  and checks the tool schemas against the limits (the schemas are unchanged). Its coverage checks
  give every bound an at-bound case and a past-bound case. `tests/test_analysis_requests.py`
  replays the cases on Flask. `test_backend_contract.py` checks that every limit is in its
  REFERENCE.md row. That check found that the `/api/triplets` and `/api/clusters` rows left out
  the lower bound 0 of `r12Min`, `r23Min` and `rMin`; the rows now state it.
- Three behaviour changes, each pinned by a fixture case:
  1. An unknown or missing kind is refused at `requestAnalysis` on every adapter. Flask used to
     answer it with a silent KDE.
  2. The worker's modulations request reads `'ALL'` and `' all '` as every element, as Flask does.
     It used to refuse them ("no site holds All").
  3. An in-process request aborted while it runs now rejects with AbortError when it finishes. It
     used to resolve.

  Every accepted request returns what it returned before.
- 19 cases where the adapters still differ are kept, each with its reason, and the JS test fixes
  their number. Eight differ in outcome: a fractional `grid` (Flask 400, worker rounds), a blank
  `grid` (48 or 8), a negative `minAmplitude`, numeric text for `bwScale`, `extent` and
  `referenceNumber` (accepted by Flask only), a null `probability` (accepted by Flask only) and a
  blank `binWidth` (accepted by Flask only). Ten are refused by both but named differently by the
  worker's engines (`numeric bandwidth`, `reference number`, `bond12`). One is read by the worker
  only (`clusterThreshold`).
- The tests the fixture replaces were removed: the worker and route request tables in
  workerGuards, pcaKdeWorker, tripletsWorker, clusters, orientationFixes, modulations and
  test_backend_validation. The NaN-result, work-budget, cache-key and kde/slice tests stay.
