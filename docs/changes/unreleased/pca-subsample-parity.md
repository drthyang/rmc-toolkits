### PCA Ellipsoid: one fit subsample in both runtimes (plan 2.6)

- Above `MAX_PCA_FIT_POINTS` (20 000 points: a pooled element cloud, or one site in a box of 28 or
  more cells per edge) the KDE volume is fitted on a subsample. Python drew it with numpy's
  generator and the browser with mulberry32, so Flask and the static app fitted different points.
  Python now draws the browser's subsample, index for index: `pca_kde.fit_subsample_indices` and
  `mulberry32_uint32`, a bit-exact port of `fitSubsampleIndices` / `mulberry32` in
  `workers/pcaKde.js`. The draw seeds mulberry32 with the seed modulo 2³², runs a partial
  Fisher–Yates with `floor(u·(n−i))` in float64, and keeps the first 20 000 positions sorted
  ascending. Python computes it in uint64 arithmetic, so it runs on numpy 1.22.
- A pooled cloud is now site-major in Python too (`displacement_cloud`: sites by reference number,
  file order within a site), as the browser builds it. RMCProfile writes each element cell by cell,
  so the two orders differed, and the subsample picks rows by position. Below the cap only
  summation round-off changes (≤ 3e-13 relative on the demo's pooled Ga; single sites are
  bit-identical).
- Python's numbers change above the cap, because it fits a different subset. On a 27 000-copy
  Gaussian site the eigenvalues move up to 0.5% and the density up to 6% of its peak. On the demo's
  pooled Se cloud (32 000 atoms) they move 1.8% and 11%. They now equal the browser's. The page asks
  for single sites, and the AI Copilot and the validation's arms never request a KDE volume. None of
  the validation cases exceeds the cap: the largest cloud is 2 560 atoms.
- New golden `pca_kde_fixture.json` (`tests/generate_pca_kde_fixture.py` → `pcaKdeParity.test.js`;
  `tests/test_pca_kde_fixture.py` checks it is current). It pins the mulberry32 stream and the drawn
  indices exactly, then two volumes above the cap to the golden rule: a 21 952-copy synthetic site
  that both runtimes build from the same stream, and the demo's pooled Se. The existing
  site-statistics, orientation, cluster and modulation goldens are unchanged.
- The browser's mulberry32 now wraps its state modulo 2³² explicitly. This makes no difference
  below 4.9 million draws, and the app draws at most 20 000.
