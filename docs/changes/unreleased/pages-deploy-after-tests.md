### CI: the public dashboard deploys only after the tests pass

- `.github/workflows/pages.yml` no longer deploys on every push to `main`. It starts when the Tests
  workflow finishes for a push to `main`, and builds and deploys that tested commit only if the run
  passed. Before, a push with failing tests still went live.
- A manual run is still possible and has the same condition: it is refused when the commit has no
  passed Tests run on `main`.
- A commit that is no longer the head of `main` is skipped, so an older build cannot land after a
  newer one. The tests stay defined once, in `tests.yml`, and run once per push.
