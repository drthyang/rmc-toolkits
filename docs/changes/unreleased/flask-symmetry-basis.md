### Detected SG and the symmetry tools for server folders: the basis from Flask (plan item 4.5)

- `/api/structure` returns `basis`, the average-structure basis that `structureFromRmc6f` adds to a
  browser-loaded run's structure: one row per reference number, with `el`, `referenceNumber`,
  `frac`, `dispA`, `count`, `elementCounts`, `mean` and `covFrac`. A new module,
  `rmc_toolkits/structure_basis.py` (`structure_basis`, `read_structure_basis`, both exported),
  computes it with the browser's steps in the browser's order. Each copy's unit-cell fraction is
  folded with C's `fmod`, and its per-axis circular mean is `frac`. Offsets are wrapped by
  `Math.round`'s rule (a tie goes towards +∞). `mean` is the arithmetic mean unwrapped about `frac`,
  `covFrac` the population covariance (÷ n), `dispA` the rms displacement through the full metric,
  and `el` the first atom's element.
- A file of coords-only atoms gets `basis: []`. A basis holding NaN or Infinity (a coordinate too
  large to fold into the cell) is a computed result, so the route answers 400 rather than writing
  null.
- In Flask mode a server folder now gets the Detected SG card, its CIF button, the Copilot's
  `symmetry` context and the `detect_symmetry` / `measure_symmetry_lowering` tools. Until now only a
  browser-loaded run had a basis. The browser path is unchanged, and so is everything the model is
  sent for one: the model-visible snapshots pass unchanged. For a server folder the model now also
  gets the symmetry block and the two tools.
- In Flask mode the Dashboard's structure follows App's `configEpoch`, as the analysis pages do. A
  Load of the same folder, or a Live Data poll, that sees a new `.rmc6f` makes it re-read the
  structure in place; a read that lands mid-write keeps the previous one, with its *previous read*
  chip. So the card and the Copilot's symmetry tools read the configuration the Copilot's Flask
  analyses read. The plots still follow Live Data. Of overlapping structure requests, only the
  latest one's answer is applied.
- New golden `structure_basis_fixture.json` (`tests/generate_structure_basis_fixture.py` →
  `structureBasisParity.test.js`; `tests/test_structure_basis_fixture.py` checks it is current). Its
  cases are the demo run, a synthetic oblique-cell configuration and a coords-only file. The
  synthetic one has a site straddling the cell edge, a mixed site whose first atom is the minority
  species, a split, a correlated, a still and a one-copy site, and coords-only atoms. The browser
  matches the golden to 1e-12 (measured: 2e-16). On the demo run the server's basis gives the
  browser's ladder, orbits and CIF byte for byte.
