### Flask mode: the AI Copilot gets the plot data a picked folder gives it

- In Flask mode each Dashboard chart fetched `/api/plot/data` for itself, so the files the Dashboard
  publishes for the Copilot carried no `plotData`. On the Demo run the dock offered 8 analyses
  instead of 11, without `assess_convergence`, `find_gr_peaks` and `compare_fit`. Its run context had
  no dataset titles or Rwp, no `pair_correlations` and no `convergence`.
- Once the Copilot is open, the Dashboard now fetches every plot file's payload, hidden ones
  included, and attaches `plotData` (or `parseError`) to the files it publishes. The χ² history is
  fetched once, for the charted log: the server already concatenates every log of its run, so no copy
  is spliced in per log. For the same files the run module now gives the browser path's tools, cards
  and run context.
- A new module, `src/plotDataCache.js`, sends one request per file and version, and InteractivePlot
  and the Dashboard both read through it. A chart already drawn costs no second request. Masked
  points (JSON `null`) are `NaN` again, as the browser parser gives them. A failed request is not
  kept, so the next reader tries again.
- A Live Data poll asks again only for the files whose version moved, and a file keeps its previous
  data until the new answer lands. A load of another folder drops the answers still in flight for the
  folder before. A failure is listed in the Problems section, hidden or collapsed charts included.
- The χ² chart's `refreshKey` is now the version of every log of its run. A Flask chart therefore
  follows a later restart log as it grows; before, it followed only the first log.
- `/api/structure` keeps the `.rmc6f`'s species order on the wire. The app's JSON provider sorts
  object keys, so `elementCounts` arrived as Ga, Se, Ta and the Copilot labelled the partial weights
  `Se-Ta`; for the same files picked in the browser they were `Ta-Se`. The route is now serialized in
  build order, so `elementCounts`, `atomIndices` and each basis row's `elementCounts` follow the
  `.rmc6f` in both runtimes. `tests/test_parsers_demo_run.py` pins the wire order.
- `flaskCopilotRun.test.jsx` builds the Flask payloads from `plot_parity_fixture.json` and serves the
  Demo's model and logs beside them. On those files the Flask Dashboard's published run gives the
  browser path's 11 tools with identical definitions, the same cards and the same context text, and
  the dataset titles and Rwp the fixture pins. It also checks one request per file and version, Live
  Data, a folder switch and a failure. The browser path's model-visible snapshots pass unchanged. In
  Flask mode the model now gets the three tools and the context blocks above.
