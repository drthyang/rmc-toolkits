### Engines: the cluster and modulation engines return their own verdicts (architecture candidate 5)

- `cluster_distortion` (`rmc_toolkits/clusters.py`, `workers/clusters.js`) always returns a
  `verdict`. Each question is a Finding `{answer, because, evidence}`: distorted, less distorted,
  breathing, rigid, a static distortion of the average cluster, split member sites, a trigonal
  lean, favoured gyration axes, and what the window did. `answer` is null where a question does
  not apply or its statistic is undefined, never for "not significant". Each Evidence record
  names its null, observed and null values, z, change, bar and side; its `clears` can be
  recomputed from the record. The rules are the Copilot tool's, ported rule for rule and echoed
  in the payload (`CLUSTER_RULES`). They are read from the unrounded statistics. The split rule
  (`isSplitSite`) moved into the cluster engine: `list_sites` imports it, and the engine reads the
  split member sites from its own clouds.
- `displacement_modulations` (`modulations.py`, `workers/modulations.js`) gives every zone-boundary
  point and strongest mode its calibrated `sigma`, its `kind`, whether it is `ordered`, `because`
  and its `bar`. σ is the tool's χ²/Šidák arithmetic, unchanged; +∞ is reported as `SIGMA_MAX` = 40
  and −∞ as null, as the tool showed them. Each element gets a `verdict` {cellDoubling,
  generalOrder, longWave}, with the q that decided each and the nearest miss. The payload also
  carries the `family` the corrections counted and the `rules`. One long-wave predicate replaces
  two (they agree on every grid point).
- `analyze_clusters` and `analyze_modulations` now only shape the data and say each finding through
  a sentence table keyed by (finding, because); their `facts` are the answers. What the model sees
  is unchanged: the model-visible snapshots pass as they were. Today's summaries and the new
  renderers give the same data and checks on 1,010 tool results. `analyze_clusters` makes one
  engine request instead of two. Before, when its second (`sites`) request failed, the split-site
  note was dropped silently, and split sites could read as rigid.
- The two engines agree: both fixtures now pin the verdicts and their vocabularies. Their
  generators refuse a case where a decided statistic sits within 10⁻⁶ of its threshold. New known
  answers: a pinched edge, rigid tetrahedra, split F₆ octahedra, rigid dimers, and order at
  q = (3/8, 0, 0), alone and under a long wave. Staleness tests (`test_clusters_fixture.py`,
  `test_modulations_fixture.py`) catch stale fixtures, and the Flask routes return the engine
  verdicts.
- Arm A+ reads the engines' verdicts, so its cluster claims follow the tool's rules where it had
  its own: trigonality counts for tetrahedra; rigidity defers to the engine's split member sites
  (not arm A's split site) and to a pattern beside a narrower edge spread, and collective
  breathing no longer blocks it. On the dev set every score and claim is unchanged; only its
  cluster and modulation readings change. Only a held-out case would exercise the rule changes.
  Approved by the manager on the user's delegation, 2026-10-19.
